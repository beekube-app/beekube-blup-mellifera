"""Préparation des notes d'un critère : échelle, moyenne par visite, effets fixes (spec §2.1)."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Optional

import numpy as np
import scipy.linalg
import scipy.sparse as sp

from .schemas import BlupV2Input, CriterionIn, QueenIn, RecordIn


@dataclass
class CriterionData:
    criterion: CriterionIn
    y: np.ndarray                         # (n,) moyenne par visite, échelle [0, 1]
    queen_ids: list[int]                  # reine de chaque observation
    X: sp.csr_matrix                      # (n, p) effets fixes, colonnes linéairement indépendantes
    x_names: list[str]
    fixed_factors: list[str]              # facteurs retenus parmi cg, evaluator, hive_class, age, f_workers
    cg_labels: list[tuple[int, int]]      # (rucher, année) de chaque observation
    factor_labels: dict[str, list]        # niveau de chaque facteur retenu (hors covariable) par observation
    recorded_queens: list[int]            # reines notées (ordre des colonnes de Z_pe)
    n_records: dict[int, int] = field(default_factory=dict)      # notes brutes par reine
    n_evaluators: dict[int, int] = field(default_factory=dict)
    n_years: dict[int, int] = field(default_factory=dict)
    n_raw_records: int = 0
    n_contemporary_groups: int = 0
    scale_mixed: bool = False
    n_farm_obs: Optional[int] = None      # observations de l'exploitation (les premières) ; None : toutes
    farm_X: Optional[sp.csr_matrix] = None    # effets fixes de ces seules observations ; None : X

    @property
    def n_obs(self) -> int:
        return int(self.y.shape[0])

    @property
    def farm_queens(self) -> list[int]:
        """Reines de l'exploitation notées sur le critère (sans les ancêtres Index Mellifera)."""
        return [q for q in self.recorded_queens if q in self.n_records]

    @property
    def prior_X(self) -> sp.csr_matrix:
        """Effets fixes servant aux a priori des variances : exploitation seule (spec volet 3 §5)."""
        return self.X if self.farm_X is None else self.farm_X

    @property
    def prior_y(self) -> np.ndarray:
        return self.y if self.n_farm_obs is None else self.y[: self.n_farm_obs]


def normalized_value(r: RecordIn, kind: str) -> float:
    if kind == "note":
        return (r.value - r.scale_min) / (r.scale_max - r.scale_min)
    return float(r.value)


def age_class(eval_year: int, born: Optional[int]) -> str:
    if born is None:
        return "?"
    age = eval_year - born
    return "0" if age <= 0 else ("1" if age == 1 else "2+")


def _independent_columns(X: np.ndarray, tol: float = 1e-8) -> list[int]:
    if X.shape[1] == 0:
        return []
    _, R, piv = scipy.linalg.qr(X, mode="economic", pivoting=True)
    d = np.abs(np.diag(R))
    rank = int(np.sum(d > tol * max(d[0], 1.0)))
    return sorted(piv[:rank].tolist())


def _design(n: int, factors_in: list[tuple[str, list]], covariates: list[tuple[str, np.ndarray]]):
    """Matrice des effets fixes : constante, facteurs à au moins deux niveaux, covariables ; colonnes
    linéairement indépendantes. Retourne (X dense, noms, facteurs retenus, niveaux par facteur retenu)."""
    cols, names, factors = [np.ones(n)], ["intercept"], []
    labels: dict[str, list] = {}
    for fname, values in factors_in:
        levels = sorted(set(values), key=repr)
        if len(levels) < 2:
            continue
        factors.append(fname)
        labels[fname] = values
        pos = {lv: k for k, lv in enumerate(levels)}
        codes = np.array([pos[v] for v in values])
        for k, lv in enumerate(levels[1:], start=1):
            cols.append((codes == k).astype(float))
            names.append(f"{fname}={lv}")
    for cname, values in covariates:
        factors.append(cname)
        cols.append(values)
        names.append(cname)
    Xd = np.column_stack(cols) if n > 0 else np.zeros((0, 1))
    keep = _independent_columns(Xd) if n > 0 else [0]
    Xd = Xd[:, keep]
    names = [names[k] for k in keep]
    factors = [f for f in factors if any(nm == f or nm.startswith(f + "=") for nm in names)]
    labels = {f: v for f, v in labels.items() if f in factors}
    return Xd, names, factors, labels


def _covariate(name: str, values: np.ndarray, rows: np.ndarray) -> list[tuple[str, np.ndarray]]:
    """Covariable retenue si elle varie sur ses lignes (`rows`) ; nulle ailleurs."""
    if rows.any() and np.var(values[rows]) > 1e-12:
        return [(name, np.where(rows, values, 0.0))]
    return []


def prepare_criterion(
    inp: BlupV2Input, criterion_id: int, f_workers: Optional[dict[int, float]] = None,
    extra_queens: Iterable[QueenIn] = (), extra_records: Iterable[RecordIn] = (),
) -> CriterionData:
    """extra_queens / extra_records : ancêtres Index Mellifera et leurs notes (spec volet 3 §5), ajoutés au
    modèle seulement si l'exploitation a des notes sur ce critère. Leurs observations suivent celles de
    l'exploitation ; elles ont leurs propres groupes (ruchers négatifs), leurs propres classes d'âge et leur
    propre covariable de consanguinité. Statistiques et a priori : exploitation seule."""
    crit = next(c for c in inp.criteria if c.id == criterion_id)
    born = {q.id: q.born for q in inp.queens}
    recs = [r for r in inp.records if r.criterion == criterion_id]
    ext = [r for r in extra_records if r.criterion == criterion_id] if recs else []
    born.update({q.id: q.born for q in extra_queens})
    scales = {(r.scale_min, r.scale_max) for r in recs}
    groups: dict[tuple, list[RecordIn]] = {}
    for r in recs + ext:
        # une visite = un rucher et une année : des notes de même visit mais de rucher ou d'année
        # différents restent des observations distinctes (groupe de contemporaines sans ambiguïté)
        if r.visit is not None:
            key = (r.queen, "v", r.visit, r.apiary, r.date.year)
        else:
            key = (r.queen, "d", r.date, r.apiary)
        groups.setdefault(key, []).append(r)
    farm_ids = {r.queen for r in recs}
    y, qids, cg, ev, hc, ag = [], [], [], [], [], []
    for key, rs in groups.items():
        first = rs[0]
        y.append(float(np.mean([normalized_value(r, crit.kind) for r in rs])))
        qids.append(key[0])
        cg.append((first.apiary, first.date.year))
        ev.append(first.evaluator if first.evaluator is not None else -1)
        external = key[0] not in farm_ids         # niveaux propres aux notes externes : aucun effet partagé
        hc.append("x" + first.hive_class if external else first.hive_class)
        age = age_class(first.date.year, born.get(key[0]))
        ag.append("x" + age if external else age)
    n = len(y)
    n_farm = sum(1 for q in qids if q in farm_ids)          # observations de l'exploitation : les n_farm premières
    farm = np.arange(n) < n_farm
    fw = np.array([f_workers.get(q, 0.0) for q in qids]) if f_workers is not None else None
    covs = []
    if fw is not None and n > 0:
        covs = _covariate("f_workers", fw, farm) + _covariate("f_workers_ext", fw, ~farm)
    factors_in = [("cg", cg), ("evaluator", ev), ("hive_class", hc), ("age", ag)]
    Xd, names, factors, labels = _design(n, factors_in, covs)
    if n_farm < n:      # a priori des variances : plan de l'exploitation seule, comme sans notes externes
        farm_cov = _covariate("f_workers", fw[:n_farm], farm[:n_farm]) if fw is not None else []
        Xf = _design(n_farm, [(f, v[:n_farm]) for f, v in factors_in], farm_cov)[0]
        farm_X = sp.csr_matrix(Xf)
    else:
        farm_X = None
    recorded = list(dict.fromkeys(qids))
    n_rec: dict[int, int] = {}
    evs: dict[int, set] = {}
    yrs: dict[int, set] = {}
    for r in recs:
        n_rec[r.queen] = n_rec.get(r.queen, 0) + 1
        evs.setdefault(r.queen, set())
        if r.evaluator is not None:
            evs[r.queen].add(r.evaluator)
        yrs.setdefault(r.queen, set()).add(r.date.year)
    X = sp.csr_matrix(Xd)
    return CriterionData(
        criterion=crit,
        y=np.array(y, dtype=float),
        queen_ids=qids,
        X=X,
        x_names=names,
        fixed_factors=factors,
        cg_labels=cg,
        factor_labels=labels,
        recorded_queens=recorded,
        n_records=n_rec,
        n_evaluators={q: len(s) for q, s in evs.items()},
        n_years={q: len(s) for q, s in yrs.items()},
        n_raw_records=len(recs),
        n_contemporary_groups=len(set(cg[:n_farm])),
        scale_mixed=crit.kind == "note" and len(scales) > 1,
        n_farm_obs=n_farm,
        farm_X=farm_X if farm_X is not None else X,
    )
