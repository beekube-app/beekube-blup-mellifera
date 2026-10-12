"""Préparation des notes d'un critère : échelle, moyenne par visite, effets fixes (spec §2.1)."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import scipy.linalg
import scipy.sparse as sp

from .schemas import BlupV2Input, CriterionIn, RecordIn


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

    @property
    def n_obs(self) -> int:
        return int(self.y.shape[0])


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


def prepare_criterion(
    inp: BlupV2Input, criterion_id: int, f_workers: Optional[dict[int, float]] = None
) -> CriterionData:
    crit = next(c for c in inp.criteria if c.id == criterion_id)
    born = {q.id: q.born for q in inp.queens}
    recs = [r for r in inp.records if r.criterion == criterion_id]
    scales = {(r.scale_min, r.scale_max) for r in recs}
    groups: dict[tuple, list[RecordIn]] = {}
    for r in recs:
        # une visite = un rucher et une année : des notes de même visit mais de rucher ou d'année
        # différents restent des observations distinctes (groupe de contemporaines sans ambiguïté)
        if r.visit is not None:
            key = (r.queen, "v", r.visit, r.apiary, r.date.year)
        else:
            key = (r.queen, "d", r.date, r.apiary)
        groups.setdefault(key, []).append(r)
    y, qids, cg, ev, hc, ag = [], [], [], [], [], []
    for key, rs in groups.items():
        first = rs[0]
        y.append(float(np.mean([normalized_value(r, crit.kind) for r in rs])))
        qids.append(key[0])
        cg.append((first.apiary, first.date.year))
        ev.append(first.evaluator if first.evaluator is not None else -1)
        hc.append(first.hive_class)
        ag.append(age_class(first.date.year, born.get(key[0])))
    n = len(y)
    cols, names, factors = [np.ones(n)], ["intercept"], []
    labels: dict[str, list] = {}
    for fname, values in (("cg", cg), ("evaluator", ev), ("hive_class", hc), ("age", ag)):
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
    if f_workers is not None and n > 0:
        fw = np.array([f_workers.get(q, 0.0) for q in qids])
        if np.var(fw) > 1e-12:
            factors.append("f_workers")
            cols.append(fw)
            names.append("f_workers")
    Xd = np.column_stack(cols) if n > 0 else np.zeros((0, 1))
    keep = _independent_columns(Xd) if n > 0 else [0]
    Xd = Xd[:, keep]
    names = [names[k] for k in keep]
    factors = [f for f in factors if any(nm == f or nm.startswith(f + "=") for nm in names)]
    labels = {f: v for f, v in labels.items() if f in factors}
    recorded = list(dict.fromkeys(qids))
    n_rec: dict[int, int] = {}
    evs: dict[int, set] = {}
    yrs: dict[int, set] = {}
    for r in recs:
        n_rec[r.queen] = n_rec.get(r.queen, 0) + 1
        evs.setdefault(r.queen, set()).add(r.evaluator)
        yrs.setdefault(r.queen, set()).add(r.date.year)
    return CriterionData(
        criterion=crit,
        y=np.array(y, dtype=float),
        queen_ids=qids,
        X=sp.csr_matrix(Xd),
        x_names=names,
        fixed_factors=factors,
        cg_labels=cg,
        factor_labels=labels,
        recorded_queens=recorded,
        n_records=n_rec,
        n_evaluators={q: len(s) for q, s in evs.items()},
        n_years={q: len(s) for q, s in yrs.items()},
        n_raw_records=len(recs),
        n_contemporary_groups=len(set(cg)),
        scale_mixed=crit.kind == "note" and len(scales) > 1,
    )
