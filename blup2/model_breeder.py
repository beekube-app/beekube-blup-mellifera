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
import scipy.optimize
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .kinship import relationship, relationship_diagonal
from .mme import inverse_entries, solve_mme
from .model_standard import (
    FEW_RECORDS_QUEENS,
    PEDIGREE_MIN_RELIABILITY,
    CriterionFit,
    FitParams,
    empty_fit,
    external_equations,
    fit_standard,
    incidence,
)
from .pedigree import BeePedigree
from .records import CriterionData
from .reliability import approx_pev, effective_records
from .schemas import CriterionParameters, WarningOut
from .variance import (
    DEFAULT_H2,
    DEFAULT_REPEATABILITY,
    NO_VARIATION_SIGMA2_P,
    phenotypic_variance,
    should_estimate,
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


def neg2_reml_spd(X, Z_list: list, Gstar_inv, logdet_gstar: float, y: np.ndarray) -> tuple[float, float]:
    """Même valeur que variance.neg2_reml_loglik (σ²_e profilée), mais factorisation sans pivotage
    (C est symétrique définie positive) : remplissage et temps divisés par deux environ sur le
    système à deux effets."""
    X = sp.csr_matrix(X)
    W = sp.hstack([X] + [sp.csr_matrix(Z) for Z in Z_list], format="csr")
    p, N = X.shape[1], len(y)
    pad = sp.block_diag([sp.csr_matrix((p, p)), sp.csr_matrix(Gstar_inv)], format="csc")
    C = (W.T @ W).tocsc() + pad
    rhs = W.T @ y
    lu = spla.splu(C, permc_spec="MMD_AT_PLUS_A", diag_pivot_thresh=0.0, options={"SymmetricMode": True})
    sol = lu.solve(rhs)
    ssr = max(float(y @ y - sol @ rhs), 1e-300)
    logdet_c = float(np.sum(np.log(np.abs(lu.U.diagonal()))))
    sigma2_e = ssr / (N - p)
    return (N - p) * (1.0 + np.log(sigma2_e)) + logdet_gstar + logdet_c, sigma2_e


@dataclass
class BreederReml:
    sigma2_AW: float
    sigma2_AQ: float
    rg: float
    sigma2_pe: float
    sigma2_e: float
    converged: bool
    n_iter: int
    genetic_at_bound: bool = False   # σ²_AW et σ²_AQ au bord inférieur (variance génétique nulle)


def reml_breeder(X, Zw, Zq, Zpe, Ainv, y: np.ndarray, start: BreederComponents, estimate_partition: bool,
                 max_evals: int = 400, bound: float = 1e4) -> BreederReml:
    """REML sans dérivées (σ²_e profilée, Nelder-Mead) du modèle à deux effets corrélés.

    θ = (log σ²_AW/σ²_e, log σ²_pe/σ²_e) et, si estimate_partition, (log ρ, atanh(rG / 0,95)) ;
    sinon ρ = σ²_AQ/σ²_AW et rG restent à leur valeur de départ. −2 log L est celle du plan
    standard (variance.neg2_reml_loglik) avec G* = blocdiag(G0* ⊗ A, I γ_pe), log|G*| =
    n_nœuds log|G0*| + n_pe log γ_pe (+ 2 log|A|, constant). Rapports bornés à [1/bound, bound].
    Non convergé (estimation inutilisable) : échec de l'optimiseur, σ²_AW ou σ²_AQ à un bord
    (inférieur ou supérieur ; σ²_AQ = ρ σ²_AW : bord de σ²_AW, ou de ρ si estimate_partition),
    σ²_pe au bord supérieur ou |rG| ≥ 0,94. σ²_pe au bord inférieur (nulle) est valide.
    """
    n = Ainv.shape[0]
    n_pe = Zpe.shape[1]
    lo, hi = -np.log(bound), np.log(bound)
    rho0 = start.sigma2_AQ / start.sigma2_AW
    z_max = float(np.arctanh(0.999))

    def system(theta: np.ndarray):
        """(γ_W, γ_pe, ρ, rG, G*⁻¹, log|G*|) pour θ ; rapports γ = variance / σ²_e."""
        g_w, g_pe = np.exp(np.clip(theta[0], lo, hi)), np.exp(np.clip(theta[1], lo, hi))
        if estimate_partition:
            rho = np.exp(np.clip(theta[2], lo, hi))
            rg = RG_MAX * np.tanh(np.clip(theta[3], -z_max, z_max))
        else:
            rho, rg = rho0, start.rg
        c = rg * g_w * np.sqrt(rho)
        g0 = np.array([[g_w, c], [c, rho * g_w]])
        Gi = sp.block_diag([genetic_inverse(g0, Ainv), sp.identity(n_pe, format="csr") / g_pe], format="csr")
        logdet = n * float(np.linalg.slogdet(g0)[1]) + n_pe * float(np.log(g_pe))
        return g_w, g_pe, rho, rg, Gi, logdet

    def f(theta: np.ndarray) -> float:
        *_, Gi, logdet = system(theta)
        return neg2_reml_spd(X, [Zw, Zq, Zpe], Gi, logdet, y)[0]

    theta0 = [np.log(start.sigma2_AW / start.sigma2_e), np.log(start.sigma2_pe / start.sigma2_e)]
    if estimate_partition:
        theta0 += [np.log(rho0), np.arctanh(start.rg / RG_MAX)]
    theta0 = np.clip(np.array(theta0, dtype=float), lo + 1.0, hi - 1.0)
    res = scipy.optimize.minimize(f, theta0, method="Nelder-Mead",
                                  options={"xatol": 1e-2, "fatol": 1e-2, "maxfev": max_evals})
    theta = np.clip(res.x, lo, hi)      # un seul θ borné (comme variance.reml), testé partout ci-dessous
    g_w, g_pe, rho, rg, Gi, logdet = system(theta)
    _, s2e = neg2_reml_spd(X, [Zw, Zq, Zpe], Gi, logdet, y)
    # Même règle que le profil standard : σ²_pe nulle est une estimation valide (répétabilité = h²).
    # Non validé sur données réelles (aucun jeu n'atteint le seuil du profil breeder) : appliqué par
    # cohérence avec le profil standard.
    genetic = [theta[0]] + ([theta[2]] if estimate_partition else [])     # log σ²_AW/σ²_e, log ρ
    unusable = (any(t <= lo + 0.05 or t >= hi - 0.05 for t in genetic) or theta[1] >= hi - 0.05
                or abs(rg) >= 0.94)
    genetic_at_bound = bool(theta[0] <= lo + 0.05 and theta[0] + float(np.log(rho)) <= lo + 0.05)
    return BreederReml(sigma2_AW=float(g_w * s2e), sigma2_AQ=float(rho * g_w * s2e), rg=float(rg),
                       sigma2_pe=float(g_pe * s2e), sigma2_e=float(s2e),
                       converged=bool(res.success) and not unusable, n_iter=int(res.nfev),
                       genetic_at_bound=genetic_at_bound)


def known_mating_share(cdata: CriterionData, ped: BeePedigree) -> float:
    """Part des reines notées de l'exploitation dont l'accouplement (nœud sire) est connu."""
    queens = cdata.farm_queens
    if not queens:
        return 0.0
    known = sum(1 for q in queens if ped.sire_of_mating(q) is not None)
    return known / len(queens)


def _published_share(share: float) -> float:
    """Part publiée (3 décimales). En repli (< 50 %), jamais 0,5 : l'arrondi d'une part de 0,4995 à
    0,5 contredirait le seuil."""
    rounded = round(share, 3)
    return min(rounded, 0.499) if share < MIN_KNOWN_MATING_SHARE else rounded


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


def _approx_pev_total(cdata: CriterionData, ped: BeePedigree, comp: BreederComponents, w_obs: list[int],
                      diag: np.ndarray) -> np.ndarray:
    """PEV approchée de T = a_W + a_Q pour chaque nœud (> exact_limit équations).

    Modèle à un effet équivalent (même variance génétique d'une note, σ²_A = σ²_colonie / ā_WW) traité
    par l'approximation du plan standard (reliability.approx_pev), puis, au nœud W de chaque colonie,
    r²(T_W) ≈ r²(g) × corr²(T_W, g), g = a_W(W) + a_Q(Q) étant la valeur génétique que mesure la note.
    Les autres nœuds (Q, S) gardent r²(g)."""
    n = ped.n_nodes
    s2a_eq = comp.sigma2_colony / comp.mean_aww
    n_eff = effective_records({"mean": [0] * cdata.n_obs, **cdata.factor_labels}, w_obs, n)
    rel_eq = 1.0 - approx_pev(ped, n_eff, s2a_eq, comp.sigma2_pe, comp.sigma2_e) / (diag * s2a_eq)
    rho = comp.sigma2_AQ / comp.sigma2_AW
    c = comp.rg * np.sqrt(rho)
    shrink = np.ones(n)
    for q in ped.queen_ids:
        w, qn = ped.w_node(q), ped.q_node(q)
        aww, aqq, awq = diag[w], diag[qn], relationship(ped, w, qn)
        cov_tg = aww + c * (awq + aww) + rho * awq        # cov(T_W, g) / σ²_AW
        var_g = aww + rho * aqq + 2.0 * c * awq           # var(g) / σ²_AW
        shrink[w] = cov_tg ** 2 / (aww * (1.0 + rho + 2.0 * c) * var_g)
    return (1.0 - np.clip(rel_eq, 0.0, 1.0) * shrink) * diag * comp.sigma2_T


def fit_breeder(cdata: CriterionData, ped: BeePedigree, Ainv: sp.csr_matrix, params: FitParams) -> CriterionFit:
    """Évalue un critère en profil breeder ; point d'entrée du module (même signature que fit_standard).
    Critère sans note : sortie vide (empty_fit), profile_used « breeder », known_mating_share 0. Moins de
    50 % des reines notées à accouplement connu : fit_standard, breeder_fallback en tête des avertissements
    et known_mating_share publiée. Sinon : ebv = T̂ du nœud W (Â_W + Â_Q), ebv_worker = Â_W du nœud W,
    ebv_queen = T̂ du nœud Q (valeur reproductrice totale de la reine, utilisée aux accouplements)."""
    cid = cdata.criterion.id
    share = known_mating_share(cdata, ped)
    if cdata.n_obs == 0:
        fit = empty_fit(cdata, ped)
        fit.parameters = fit.parameters.model_copy(
            update={"profile_used": "breeder", "known_mating_share": _published_share(share)})
        return fit
    if share < MIN_KNOWN_MATING_SHARE:
        fit = fit_standard(cdata, ped, Ainv, params)
        fit.parameters = fit.parameters.model_copy(update={"known_mating_share": _published_share(share)})
        pct = min(round(100 * share), 49)      # jamais « 50 % » à côté de « (< 50 %) »
        fit.warnings.insert(0, WarningOut(
            code="breeder_fallback", criterion=cid,
            message=f"{pct} % des reines notées ont un accouplement connu (< 50 %) : "
                    f"profil standard utilisé pour ce critère"))
        return fit
    warnings: list[WarningOut] = []
    diag = relationship_diagonal(ped)
    n = ped.n_nodes
    w_obs = [ped.w_node(q) for q in cdata.queen_ids]
    q_obs = [ped.q_node(q) for q in cdata.queen_ids]
    rec_pos = {q: k for k, q in enumerate(cdata.recorded_queens)}
    n_q = len(rec_pos)                    # dimension de Z_pe (ancêtres notés compris)
    farm_queens = cdata.farm_queens       # a priori et compteurs : exploitation seule (spec volet 3 §5)
    n_farm = len(farm_queens)
    Zw, Zq = incidence(w_obs, n), incidence(q_obs, n)
    Zpe = incidence([rec_pos[q] for q in cdata.queen_ids], n_q)
    rec_w = [ped.w_node(q) for q in farm_queens]
    rec_q = [ped.q_node(q) for q in farm_queens]
    mean_aww = float(np.mean(diag[rec_w]))
    mean_aqq = float(np.mean(diag[rec_q]))
    mean_awq = float(np.mean([relationship(ped, w, q) for w, q in zip(rec_w, rec_q)]))
    raw_p = phenotypic_variance(cdata.prior_X, cdata.prior_y, floor=0.0)
    no_variation = raw_p < NO_VARIATION_SIGMA2_P
    if no_variation:   # pas de REML : il partirait vers des variances de 1e-300 (index aberrants)
        warnings.append(WarningOut(
            code="no_variation", criterion=cid,
            message="notes sans variation une fois retirés les effets fixes : le critère ne classe pas les reines"))
    sigma2_p = max(raw_p, 1e-8)       # plancher de phenotypic_variance (σ²_P sert de diviseur)
    source, h2, rep, ratio, rg = _breeder_source(cdata)
    comp = prior_breeder(sigma2_p, h2, rep, ratio, rg, mean_aww, mean_aqq, mean_awq)
    mode = params.estimate_variances
    if source == "prior" and not no_variation and should_estimate(
            mode, n_farm, cdata.n_raw_records, cdata.n_contemporary_groups):
        res = reml_breeder(cdata.X, Zw, Zq, Zpe, Ainv, cdata.y, comp, estimate_partition=(mode == "always"))
        if res.converged:
            comp = BreederComponents(sigma2_AW=res.sigma2_AW, sigma2_AQ=res.sigma2_AQ, rg=res.rg,
                                     sigma2_pe=res.sigma2_pe, sigma2_e=res.sigma2_e, mean_aww=mean_aww,
                                     mean_aqq=mean_aqq, mean_awq=mean_awq)
            source = "reml"
        elif res.genetic_at_bound:
            warnings.append(WarningOut(
                code="low_heritability", criterion=cid,
                message="Héritabilité estimée nulle sur vos données : valeur par défaut appliquée"))
        else:
            warnings.append(WarningOut(
                code="reml_not_converged", criterion=cid,
                message=f"REML non convergé ({res.n_iter} évaluations) : paramètres a priori utilisés"))
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
    elif sol.n_equations - external_equations(cdata, params, 2) <= params.exact_limit:
        w = np.array([ped.w_node(q) for q in ped.queen_ids], dtype=int)
        a, b, k = sol.n_fixed + w, sol.n_fixed + n + w, len(w)
        v = inverse_entries(sol, np.r_[a, b, a], np.r_[a, b, b], comp.sigma2_e)
        pev, method = np.full(n, np.nan), "exact"          # seules les reines de l'exploitation sont lues
        pev[w] = v[:k] + v[k:2 * k] + 2.0 * v[2 * k:]
    else:
        pev = _approx_pev_total(cdata, ped, comp, w_obs, diag)
        method = "approx"
        warnings.append(WarningOut(
            code="reliability_approx", criterion=cid,
            message=f"{sol.n_equations} équations (> {params.exact_limit}) : fiabilités approchées"))
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
    if n_farm < FEW_RECORDS_QUEENS:
        warnings.append(WarningOut(code="few_records", criterion=cid,
                                   message=f"seulement {n_farm} reine(s) notée(s) : index peu fiables"))
    if cdata.scale_mixed:
        warnings.append(WarningOut(code="scale_mixed", criterion=cid,
                                   message="plusieurs échelles de notes : notes ramenées sur [0, 1]"))
    parameters = CriterionParameters(
        source=source, h2=comp.h2, repeatability=comp.repeatability, sigma2_a=comp.sigma2_AW,
        sigma2_pe=comp.sigma2_pe, sigma2_e=comp.sigma2_e, sigma2_aq=comp.sigma2_AQ, rg=comp.rg,
        profile_used="breeder", known_mating_share=_published_share(share), n_records=cdata.n_raw_records,
        n_queens=n_farm, n_contemporary_groups=cdata.n_contemporary_groups, reliability_method=method,
    )
    return CriterionFit(
        criterion_id=cid, parameters=parameters, ebv=ebv, reliability=rel, ebv_queen=ebv_q, status=status,
        sigma_index=float(np.sqrt(comp.sigma2_T * mean_aww)), sigma_queen=float(np.sqrt(comp.sigma2_T)),
        n_records=dict(cdata.n_records), n_evaluators=dict(cdata.n_evaluators), n_years=dict(cdata.n_years),
        warnings=warnings, ebv_worker=ebv_w, sigma_worker=float(np.sqrt(comp.sigma2_AW * mean_aww)),
    )
