"""Profil standard (spec §2.1) : y = Xb + Z_W a + Z_pe pe + e, a ~ N(0, A σ²_A)."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import scipy.sparse as sp

from .kinship import relationship_diagonal
from .mme import EXACT_LIMIT, solve_mme
from .pedigree import BeePedigree
from .records import CriterionData
from .reliability import approx_pev, effective_records
from .schemas import CriterionParameters, WarningOut
from .variance import (
    DEFAULT_H2,
    DEFAULT_REPEATABILITY,
    NO_VARIATION_SIGMA2_P,
    VarianceComponents,
    phenotypic_variance,
    prior_components,
    reml,
    should_estimate,
)

FEW_RECORDS_QUEENS = 10      # avertissement few_records en dessous de ce nombre de reines notées
PEDIGREE_MIN_RELIABILITY = 1e-3


@dataclass
class FitParams:
    estimate_variances: str = "auto"      # "auto" | "always" | "never"
    exact_limit: int = EXACT_LIMIT


@dataclass
class CriterionFit:
    criterion_id: int
    parameters: CriterionParameters
    ebv: dict[int, float]                 # valeur publiée (nœud W ; breeder : Â_W + Â_Q du nœud W)
    reliability: dict[int, float]
    ebv_queen: dict[int, float]           # valeur reproductrice de la reine (nœud Q) -> index_queen
    status: dict[int, str]                # evaluated | pedigree | none
    sigma_index: float                    # écart-type génétique du nœud publié (√σ²_a)
    sigma_queen: float                    # √σ²_A (nœud Q)
    n_records: dict[int, int] = field(default_factory=dict)
    n_evaluators: dict[int, int] = field(default_factory=dict)
    n_years: dict[int, int] = field(default_factory=dict)
    warnings: list[WarningOut] = field(default_factory=list)
    ebv_worker: Optional[dict[int, float]] = None   # breeder uniquement
    sigma_worker: Optional[float] = None            # breeder uniquement


def incidence(rows: list[int], n_cols: int) -> sp.csr_matrix:
    n = len(rows)
    return sp.csr_matrix((np.ones(n), (np.arange(n), np.asarray(rows, dtype=int))), shape=(n, n_cols))


def _variance_source(cdata: CriterionData) -> tuple[str, float, float]:
    c = cdata.criterion
    if c.h2 is not None or c.repeatability is not None:
        h2 = c.h2 if c.h2 is not None else min(DEFAULT_H2, c.repeatability)
        rep = c.repeatability if c.repeatability is not None else max(DEFAULT_REPEATABILITY, h2)
        return "input", h2, rep
    return "prior", DEFAULT_H2, DEFAULT_REPEATABILITY


def empty_fit(cdata: CriterionData, ped: BeePedigree) -> CriterionFit:
    cid = cdata.criterion.id
    params = CriterionParameters(
        source="prior", h2=DEFAULT_H2, repeatability=DEFAULT_REPEATABILITY, sigma2_a=0.0, sigma2_pe=0.0,
        sigma2_e=0.0, n_records=0, n_queens=0, n_contemporary_groups=0, reliability_method="none",
    )
    return CriterionFit(
        criterion_id=cid, parameters=params, ebv={}, reliability={}, ebv_queen={},
        status={q: "none" for q in ped.queen_ids}, sigma_index=0.0, sigma_queen=0.0,
        warnings=[WarningOut(code="few_records", message="aucune note pour ce critère", criterion=cid)],
    )


def fit_standard(cdata: CriterionData, ped: BeePedigree, Ainv: sp.csr_matrix, params: FitParams) -> CriterionFit:
    cid = cdata.criterion.id
    if cdata.n_obs == 0:
        return empty_fit(cdata, ped)
    warnings: list[WarningOut] = []
    diag = relationship_diagonal(ped)
    n_nodes = ped.n_nodes
    w_obs = [ped.w_node(q) for q in cdata.queen_ids]
    rec_pos = {q: k for k, q in enumerate(cdata.recorded_queens)}
    Zw = incidence(w_obs, n_nodes)
    Zpe = incidence([rec_pos[q] for q in cdata.queen_ids], len(rec_pos))
    farm_queens = cdata.farm_queens       # a priori et compteurs : exploitation seule (spec volet 3 §5)
    mean_aww = float(np.mean([diag[ped.w_node(q)] for q in farm_queens]))
    raw_p = phenotypic_variance(cdata.prior_X, cdata.prior_y, floor=0.0)
    no_variation = raw_p < NO_VARIATION_SIGMA2_P
    if no_variation:   # pas de REML : il partirait vers des variances de 1e-300 (index aberrants)
        warnings.append(WarningOut(
            code="no_variation", criterion=cid,
            message="notes sans variation une fois retirés les effets fixes : le critère ne classe pas les reines"))
    sigma2_p = max(raw_p, 1e-8)       # plancher de phenotypic_variance (σ²_P sert de diviseur)
    source, h2, rep = _variance_source(cdata)
    comp = prior_components(sigma2_p, h2, rep, mean_aww)
    n_q = len(cdata.recorded_queens)      # dimension de Z_pe (ancêtres notés compris)
    n_farm = len(farm_queens)
    if source == "prior" and not no_variation and should_estimate(
            params.estimate_variances, n_farm, cdata.n_raw_records, cdata.n_contemporary_groups):
        I_pe = sp.identity(n_q, format="csr")
        res = reml(cdata.X, [Zw, Zpe], [Ainv, I_pe], cdata.y, [comp.sigma2_A, comp.sigma2_pe], comp.sigma2_e)
        if res.converged:
            comp = VarianceComponents(sigma2_A=res.variances[0], sigma2_pe=res.variances[1],
                                      sigma2_e=res.sigma2_e, mean_aww=mean_aww)
            source = "reml"
        elif res.genetic_at_bound:
            warnings.append(WarningOut(
                code="low_heritability", criterion=cid,
                message="Héritabilité estimée nulle sur vos données : valeur par défaut appliquée"))
        else:
            warnings.append(WarningOut(
                code="reml_not_converged", criterion=cid,
                message=f"REML non convergé ({res.n_iter} évaluations) : paramètres a priori utilisés"))
    Ginv = sp.block_diag([Ainv / comp.sigma2_A, sp.identity(n_q, format="csr") / comp.sigma2_pe], format="csr")
    sol = solve_mme(cdata.X, [Zw, Zpe], Ginv, comp.sigma2_e, cdata.y, need_inverse=True,
                    exact_limit=params.exact_limit)
    a_hat = sol.u[0]
    if sol.pev is not None:
        pev, method = sol.pev[0], "exact"
    else:
        n_eff = effective_records({"mean": [0] * cdata.n_obs, **cdata.factor_labels}, w_obs, n_nodes)
        pev, method = approx_pev(ped, n_eff, comp.sigma2_A, comp.sigma2_pe, comp.sigma2_e), "approx"
        warnings.append(WarningOut(
            code="reliability_approx", criterion=cid,
            message=f"{sol.n_equations} équations (> {params.exact_limit}) : fiabilités approchées"))
    rel_nodes = np.clip(1.0 - pev / (diag * comp.sigma2_A), 0.0, 1.0)
    ebv, ebv_q, rel, status = {}, {}, {}, {}
    for q in ped.queen_ids:
        w, qn = ped.w_node(q), ped.q_node(q)
        ebv[q], ebv_q[q], rel[q] = float(a_hat[w]), float(a_hat[qn]), float(rel_nodes[w])
        if q in cdata.n_records:
            status[q] = "evaluated"
        elif rel[q] >= PEDIGREE_MIN_RELIABILITY:
            status[q] = "pedigree"
        else:
            status[q] = "none"
    if n_farm < FEW_RECORDS_QUEENS:
        warnings.append(WarningOut(code="few_records", criterion=cid,
                                   message=f"seulement {n_farm} reine(s) notée(s) : index peu fiables"))
    if cdata.scale_mixed:
        warnings.append(WarningOut(code="scale_mixed", criterion=cid,
                                   message="plusieurs échelles de notes : notes ramenées sur [0, 1]"))
    parameters = CriterionParameters(
        source=source, h2=comp.h2, repeatability=comp.repeatability, sigma2_a=comp.sigma2_a,
        sigma2_pe=comp.sigma2_pe, sigma2_e=comp.sigma2_e, profile_used="standard",
        n_records=cdata.n_raw_records, n_queens=n_farm, n_contemporary_groups=cdata.n_contemporary_groups,
        reliability_method=method,
    )
    return CriterionFit(
        criterion_id=cid, parameters=parameters, ebv=ebv, reliability=rel, ebv_queen=ebv_q, status=status,
        sigma_index=float(np.sqrt(comp.sigma2_a)), sigma_queen=float(np.sqrt(comp.sigma2_A)),
        n_records=dict(cdata.n_records), n_evaluators=dict(cdata.n_evaluators), n_years=dict(cdata.n_years),
        warnings=warnings,
    )
