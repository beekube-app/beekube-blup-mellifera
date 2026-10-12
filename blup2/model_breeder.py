"""Profil breeder (spec §2.2) : effets ouvrières et reine corrélés [BB éq. 1–5 ; Bienefeld et al. 2007].

y = Xb + Z_W a_W + Z_Q a_Q + Z_pe pe + e
[a_W ; a_Q] ~ N(0, G0 ⊗ A), G0 = [[σ²_AW, rG σ_AW σ_AQ], [rG σ_AW σ_AQ, σ²_AQ]]
a_W est porté par le nœud W de la colonie, a_Q par le nœud Q de sa reine (même A).

Le solveur du plan standard est réutilisé tel quel : effets aléatoires [a_W (tous les nœuds) ;
a_Q (tous les nœuds) ; pe] et G⁻¹ = blocdiag(G0⁻¹ ⊗ A⁻¹, I / σ²_pe).

Valeurs publiées (T = a_W + a_Q, valeur totale d'un nœud) :
- ebv        = T̂ du nœud W = Â_W + Â_Q du groupe d'ouvrières (fille moyenne de l'accouplement, BB éq. 2) ;
- ebv_worker = Â_W du nœud W (composante ouvrières) ;
- ebv_queen  = T̂ du nœud Q (valeur reproductrice totale de la reine, utilisée pour les accouplements).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp

from .kinship import relationship, relationship_diagonal
from .mme import solve_mme
from .model_standard import (
    FEW_RECORDS_QUEENS,
    PEDIGREE_MIN_RELIABILITY,
    CriterionFit,
    FitParams,
    empty_fit,
    fit_standard,
    incidence,
)
from .pedigree import BeePedigree
from .records import CriterionData
from .schemas import CriterionParameters, WarningOut
from .variance import (
    DEFAULT_H2,
    DEFAULT_REPEATABILITY,
    phenotypic_variance,
)

DEFAULT_QUEEN_RATIO = 0.5     # σ²_AQ / σ²_AW [BB p. 10]
DEFAULT_RG = -0.5             # corrélation génétique reine / ouvrières [BB p. 10]
RG_MAX = 0.95                 # |rG| ≤ 0,95 : G0 nettement définie positive
MIN_KNOWN_MATING_SHARE = 0.5  # spec §2.2 : repli standard en dessous de 50 %


@dataclass
class BreederComponents:
    sigma2_AW: float     # variance de base de l'effet ouvrières
    sigma2_AQ: float     # variance de base de l'effet reine
    rg: float
    sigma2_pe: float
    sigma2_e: float
    mean_aww: float      # moyennes sur les colonies notées : a(W, W), a(Q, Q), a(W, Q)
    mean_aqq: float
    mean_awq: float

    @property
    def cov(self) -> float:
        return self.rg * float(np.sqrt(self.sigma2_AW * self.sigma2_AQ))

    @property
    def g0(self) -> np.ndarray:
        return np.array([[self.sigma2_AW, self.cov], [self.cov, self.sigma2_AQ]])

    @property
    def sigma2_T(self) -> float:   # variance de base de la valeur totale T = a_W + a_Q
        return self.sigma2_AW + self.sigma2_AQ + 2.0 * self.cov

    @property
    def sigma2_colony(self) -> float:   # variance génétique d'une note : var(a_W(W) + a_Q(Q)) moyenne
        return (self.sigma2_AW * self.mean_aww + self.sigma2_AQ * self.mean_aqq
                + 2.0 * self.cov * self.mean_awq)

    @property
    def sigma2_p(self) -> float:
        return self.sigma2_colony + self.sigma2_pe + self.sigma2_e

    @property
    def h2(self) -> float:
        return self.sigma2_colony / self.sigma2_p

    @property
    def repeatability(self) -> float:
        return (self.sigma2_colony + self.sigma2_pe) / self.sigma2_p


def clamp_rg(rg: float) -> float:
    return float(np.clip(rg, -RG_MAX, RG_MAX))


def prior_breeder(sigma2_p: float, h2: float, repeatability: float, queen_ratio: float, rg: float,
                  mean_aww: float, mean_aqq: float, mean_awq: float) -> BreederComponents:
    """σ²_colonie = h² σ²_P = σ²_AW (ā_WW + ρ ā_QQ + 2 rG √ρ ā_WQ), ρ = σ²_AQ / σ²_AW ;
    σ²_pe = (r − h²) σ²_P ; σ²_e = (1 − r) σ²_P (spec §2.3)."""
    repeatability = max(repeatability, h2)
    rg = clamp_rg(rg)
    k = mean_aww + queen_ratio * mean_aqq + 2.0 * rg * np.sqrt(queen_ratio) * mean_awq
    s2aw = h2 * sigma2_p / k
    return BreederComponents(
        sigma2_AW=s2aw, sigma2_AQ=queen_ratio * s2aw, rg=rg,
        sigma2_pe=max(repeatability - h2, 1e-6) * sigma2_p, sigma2_e=(1.0 - repeatability) * sigma2_p,
        mean_aww=mean_aww, mean_aqq=mean_aqq, mean_awq=mean_awq,
    )


def genetic_inverse(g0: np.ndarray, Ainv: sp.spmatrix) -> sp.csr_matrix:
    """(G0 ⊗ A)⁻¹ = G0⁻¹ ⊗ A⁻¹, ordre [a_W (tous les nœuds) ; a_Q (tous les nœuds)]."""
    return sp.kron(sp.csr_matrix(np.linalg.inv(g0)), Ainv, format="csr")


def known_mating_share(cdata: CriterionData, ped: BeePedigree) -> float:
    """Part des reines notées dont l'accouplement (nœud sire) est connu."""
    if not cdata.recorded_queens:
        return 0.0
    known = sum(1 for q in cdata.recorded_queens if ped.sire_of_mating(q) is not None)
    return known / len(cdata.recorded_queens)


def _breeder_source(cdata: CriterionData) -> tuple[str, float, float, float, float]:
    c = cdata.criterion
    given = (c.h2, c.repeatability, c.queen_ratio, c.rg)
    # même règle qu'en profil standard : répétabilité seule fournie -> h² = min(0,25 ; répétabilité)
    h2 = c.h2 if c.h2 is not None else (
        DEFAULT_H2 if c.repeatability is None else min(DEFAULT_H2, c.repeatability))
    rep = c.repeatability if c.repeatability is not None else max(DEFAULT_REPEATABILITY, h2)
    ratio = c.queen_ratio if c.queen_ratio is not None else DEFAULT_QUEEN_RATIO
    rg = c.rg if c.rg is not None else DEFAULT_RG
    source = "input" if any(v is not None for v in given) else "prior"
    return source, h2, rep, ratio, rg


def fit_breeder(cdata: CriterionData, ped: BeePedigree, Ainv: sp.csr_matrix, params: FitParams) -> CriterionFit:
    cid = cdata.criterion.id
    share = known_mating_share(cdata, ped)
    if cdata.n_obs == 0:
        fit = empty_fit(cdata, ped)
        fit.parameters = fit.parameters.model_copy(
            update={"profile_used": "breeder", "known_mating_share": round(share, 3)})
        return fit
    if share < MIN_KNOWN_MATING_SHARE:
        fit = fit_standard(cdata, ped, Ainv, params)
        fit.parameters = fit.parameters.model_copy(update={"known_mating_share": round(share, 3)})
        fit.warnings.insert(0, WarningOut(
            code="breeder_fallback", criterion=cid,
            message=f"{round(100 * share)} % des reines notées ont un accouplement connu (< 50 %) : "
                    f"profil standard utilisé pour ce critère"))
        return fit
    warnings: list[WarningOut] = []
    diag = relationship_diagonal(ped)
    n = ped.n_nodes
    w_obs = [ped.w_node(q) for q in cdata.queen_ids]
    q_obs = [ped.q_node(q) for q in cdata.queen_ids]
    rec_pos = {q: k for k, q in enumerate(cdata.recorded_queens)}
    n_q = len(rec_pos)
    Zw, Zq = incidence(w_obs, n), incidence(q_obs, n)
    Zpe = incidence([rec_pos[q] for q in cdata.queen_ids], n_q)
    rec_w = [ped.w_node(q) for q in cdata.recorded_queens]
    rec_q = [ped.q_node(q) for q in cdata.recorded_queens]
    mean_aww = float(np.mean(diag[rec_w]))
    mean_aqq = float(np.mean(diag[rec_q]))
    mean_awq = float(np.mean([relationship(ped, w, q) for w, q in zip(rec_w, rec_q)]))
    sigma2_p = phenotypic_variance(cdata.X, cdata.y)
    source, h2, rep, ratio, rg = _breeder_source(cdata)
    comp = prior_breeder(sigma2_p, h2, rep, ratio, rg, mean_aww, mean_aqq, mean_awq)
    Ginv = sp.block_diag([genetic_inverse(comp.g0, Ainv), sp.identity(n_q, format="csr") / comp.sigma2_pe],
                         format="csr")
    sol = solve_mme(cdata.X, [Zw, Zq, Zpe], Ginv, comp.sigma2_e, cdata.y, need_inverse=True,
                    exact_limit=params.exact_limit)
    a_w, a_q = sol.u[0], sol.u[1]
    total = a_w + a_q
    prior_var = diag * comp.sigma2_T
    if sol.cinv is not None:   # PEV(T) = PEV(a_W) + PEV(a_Q) + 2 cov, bloc 2 × 2 de σ²_e C⁻¹
        p, idx, cinv = sol.n_fixed, np.arange(n), sol.cinv
        pev = cinv[p + idx, p + idx] + cinv[p + n + idx, p + n + idx] + 2.0 * cinv[p + idx, p + n + idx]
        method = "exact"
    else:   # fiabilité approchée : tâche 2
        raise NotImplementedError(f"{sol.n_equations} équations (> {params.exact_limit}) : "
                                  "fiabilité approchée du profil breeder non encore écrite")
    rel_nodes = np.clip(1.0 - pev / prior_var, 0.0, 1.0)
    ebv, ebv_w, ebv_q, rel, status = {}, {}, {}, {}, {}
    for q in ped.queen_ids:
        w, qn = ped.w_node(q), ped.q_node(q)
        ebv[q], ebv_w[q], ebv_q[q] = float(total[w]), float(a_w[w]), float(total[qn])
        rel[q] = float(rel_nodes[w])
        if q in cdata.n_records:
            status[q] = "evaluated"
        elif rel[q] >= PEDIGREE_MIN_RELIABILITY:
            status[q] = "pedigree"
        else:
            status[q] = "none"
    if n_q < FEW_RECORDS_QUEENS:
        warnings.append(WarningOut(code="few_records", criterion=cid,
                                   message=f"seulement {n_q} reine(s) notée(s) : index peu fiables"))
    if cdata.scale_mixed:
        warnings.append(WarningOut(code="scale_mixed", criterion=cid,
                                   message="plusieurs échelles de notes : notes ramenées sur [0, 1]"))
    parameters = CriterionParameters(
        source=source, h2=comp.h2, repeatability=comp.repeatability, sigma2_a=comp.sigma2_AW,
        sigma2_pe=comp.sigma2_pe, sigma2_e=comp.sigma2_e, sigma2_aq=comp.sigma2_AQ, rg=comp.rg,
        profile_used="breeder", known_mating_share=round(share, 3), n_records=cdata.n_raw_records,
        n_queens=n_q, n_contemporary_groups=cdata.n_contemporary_groups, reliability_method=method,
    )
    return CriterionFit(
        criterion_id=cid, parameters=parameters, ebv=ebv, reliability=rel, ebv_queen=ebv_q, status=status,
        sigma_index=float(np.sqrt(comp.sigma2_T * mean_aww)), sigma_queen=float(np.sqrt(comp.sigma2_T)),
        n_records=dict(cdata.n_records), n_evaluators=dict(cdata.n_evaluators), n_years=dict(cdata.n_years),
        warnings=warnings, ebv_worker=ebv_w, sigma_worker=float(np.sqrt(comp.sigma2_AW * mean_aww)),
    )
