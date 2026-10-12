"""Composantes de variance : a priori (spec §2.3) et REML.

Le REML est maximisé sans dérivées (DF-REML, Graser, Smith & Tier 1987 ; Meyer 1989) :
−2 log L = (N − p)·(1 + log σ̂²_e) + log|G*| + log|C|, σ̂²_e = (yᵀy − solᵀ RHS) / (N − p),
G = σ²_e G*, C = WᵀW + diag(0, G*⁻¹). log|C| est lu sur la factorisation LU creuse.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.optimize
import scipy.sparse as sp
import scipy.sparse.linalg as spla

DEFAULT_H2 = 0.25
DEFAULT_REPEATABILITY = 0.30   # REML de l'exploitation la plus riche (spec §2.3, validation §4.1)
REML_MIN_QUEENS, REML_MIN_RECORDS, REML_MIN_CG = 100, 300, 2
NO_VARIATION_SIGMA2_P = 1e-12  # σ²_P brute en dessous : notes sans variation (avertissement no_variation)


@dataclass
class VarianceComponents:
    sigma2_A: float      # variance génétique de base de A
    sigma2_pe: float
    sigma2_e: float
    mean_aww: float      # a_WW moyen des colonies notées

    @property
    def sigma2_a(self) -> float:   # variance génétique de la colonie (nœud W)
        return self.sigma2_A * self.mean_aww

    @property
    def sigma2_p(self) -> float:
        return self.sigma2_a + self.sigma2_pe + self.sigma2_e

    @property
    def h2(self) -> float:
        return self.sigma2_a / self.sigma2_p

    @property
    def repeatability(self) -> float:
        return (self.sigma2_a + self.sigma2_pe) / self.sigma2_p


@dataclass
class RemlResult:
    variances: list[float]
    sigma2_e: float
    converged: bool
    n_iter: int
    genetic_at_bound: bool = False   # variance génétique (premier effet) au bord inférieur


def phenotypic_variance(X: sp.spmatrix, y: np.ndarray, floor: float = 1e-8) -> float:
    """Variance résiduelle des notes après effets fixes (moindres carrés), au moins `floor`
    (floor = 0 : valeur brute, pour détecter un critère sans variation)."""
    Xd = X.toarray()
    n, p = Xd.shape
    if n == 0:
        return floor
    b, *_ = np.linalg.lstsq(Xd, y, rcond=None)
    r = y - Xd @ b
    dof = n - p if n > p else n
    return max(float(r @ r) / dof, floor)


def prior_components(sigma2_p: float, h2: float, repeatability: float, mean_aww: float) -> VarianceComponents:
    """σ²_a = h² σ²_P ; σ²_pe = (r − h²) σ²_P ; σ²_e = (1 − r) σ²_P ; σ²_A = σ²_a / a_WW moyen."""
    repeatability = max(repeatability, h2)
    s2a = h2 * sigma2_p
    s2pe = max(repeatability - h2, 1e-6) * sigma2_p
    s2e = (1.0 - repeatability) * sigma2_p
    return VarianceComponents(sigma2_A=s2a / mean_aww, sigma2_pe=s2pe, sigma2_e=s2e, mean_aww=mean_aww)


def should_estimate(mode: str, n_queens: int, n_records: int, n_cg: int) -> bool:
    if mode == "always":
        return True
    if mode == "never":
        return False
    return n_queens >= REML_MIN_QUEENS and n_records >= REML_MIN_RECORDS and n_cg >= REML_MIN_CG


def neg2_reml_loglik(
    X: sp.spmatrix, Z_list: list[sp.spmatrix], Gstar_inv: sp.spmatrix, logdet_gstar: float, y: np.ndarray,
    sigma2_e: float | None = None,
) -> tuple[float, float]:
    """(−2 log L REML à une constante près, σ²_e utilisé). σ²_e = None : valeur profilée σ̂²_e."""
    X = sp.csr_matrix(X)
    W = sp.hstack([X] + [sp.csr_matrix(Z) for Z in Z_list], format="csr")
    p = X.shape[1]
    N = len(y)
    pad = sp.block_diag([sp.csr_matrix((p, p)), sp.csr_matrix(Gstar_inv)], format="csc")
    C = (W.T @ W).tocsc() + pad
    rhs = W.T @ y
    lu = spla.splu(C, permc_spec="MMD_AT_PLUS_A", options={"SymmetricMode": True})
    sol = lu.solve(rhs)
    ssr = max(float(y @ y - sol @ rhs), 1e-300)
    logdet_c = float(np.sum(np.log(np.abs(lu.U.diagonal()))))
    if sigma2_e is None:
        sigma2_e = ssr / (N - p)
        return (N - p) * (1.0 + np.log(sigma2_e)) + logdet_gstar + logdet_c, sigma2_e
    return (N - p) * np.log(sigma2_e) + logdet_gstar + logdet_c + ssr / sigma2_e, sigma2_e


def reml(
    X: sp.spmatrix,
    Z_list: list[sp.spmatrix],
    K_inv_list: list[sp.spmatrix],
    y: np.ndarray,
    start: list[float],
    sigma2_e_start: float,
    max_evals: int = 400,
    bound: float = 1e4,
) -> RemlResult:
    """REML pour k effets aléatoires indépendants u_k ~ N(0, K_k σ²_k), e ~ N(0, I σ²_e).

    Contrat : le premier effet de Z_list (et de K_inv_list) est l'effet génétique ; les suivants
    (environnement permanent…) sont secondaires.
    Paramètres optimisés : θ_k = log(σ²_k / σ²_e) (Nelder-Mead), σ²_e profilée. Rapports bornés à
    [1/bound, bound]. Non convergé (estimation inutilisable) si l'optimiseur échoue, si la variance
    génétique atteint la borne inférieure (genetic_at_bound) ou si une variance quelconque atteint
    la borne supérieure (variance infinie). Une variance secondaire à la borne inférieure (nulle)
    est une estimation valide.
    """
    sizes = [K.shape[0] for K in K_inv_list]
    lo, hi = -np.log(bound), np.log(bound)

    def f(theta: np.ndarray) -> float:
        theta = np.clip(theta, lo, hi)
        g = np.exp(theta)
        Gi = sp.block_diag([K / gk for K, gk in zip(K_inv_list, g)], format="csr")
        logdet = float(sum(q * t for q, t in zip(sizes, theta)))
        return neg2_reml_loglik(X, Z_list, Gi, logdet, y)[0]

    theta0 = np.log(np.maximum(np.asarray(start, dtype=float), 1e-12) / max(sigma2_e_start, 1e-12))
    theta0 = np.clip(theta0, lo + 1.0, hi - 1.0)
    res = scipy.optimize.minimize(f, theta0, method="Nelder-Mead",
                                  options={"xatol": 1e-2, "fatol": 1e-2, "maxfev": max_evals})
    theta = np.clip(res.x, lo, hi)
    g = np.exp(theta)
    Gi = sp.block_diag([K / gk for K, gk in zip(K_inv_list, g)], format="csr")
    _, s2e = neg2_reml_loglik(X, Z_list, Gi, float(sum(q * t for q, t in zip(sizes, theta))), y)
    # Une variance secondaire nulle (ex. environnement permanent) est une estimation REML valide
    # (bord de l'espace des paramètres, fréquent sur données réelles : répétabilité ≈ h²) ; seules
    # une variance génétique (premier effet) nulle ou une variance infinie rendent l'estimation
    # inutilisable.
    genetic_at_bound = bool(theta[0] <= lo + 0.05)
    unusable = genetic_at_bound or bool(np.any(theta >= hi - 0.05))
    return RemlResult(variances=[float(gk * s2e) for gk in g], sigma2_e=float(s2e),
                      converged=bool(res.success) and not unusable, n_iter=int(res.nfev),
                      genetic_at_bound=genetic_at_bound)


em_reml = reml   # nom de l'interface prévue (contexte des plans) ; l'algorithme est un DF-REML
