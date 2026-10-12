"""Équations du modèle mixte (MME) creuses.

LHS = WᵀW + diag(0, σ²_e G⁻¹), RHS = Wᵀy, W = [X | Z_1 | … | Z_k].
Les solutions sont les BLUE / BLUP ; σ²_e · C⁻¹ donne les variances d'erreur de prédiction (PEV).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

EXACT_LIMIT = 6000


@dataclass
class MMESolution:
    beta: np.ndarray
    u: list[np.ndarray]                    # solutions des effets aléatoires, un tableau par bloc
    pev: Optional[list[np.ndarray]]        # σ²_e · diag(C⁻¹) par bloc (None si non calculée)
    cinv: Optional[np.ndarray]             # σ²_e · C⁻¹ dense (exact uniquement), ordre [b, u_1, …, u_k]
    solve: Callable[[np.ndarray], np.ndarray]  # v -> C⁻¹ v (factorisation creuse, sans σ²_e)
    n_equations: int
    n_fixed: int
    sizes: list[int]
    exact: bool
    ssr: float                             # yᵀy − solᵀ·RHS (mise à jour EM de σ²_e)


def solve_mme(
    X: sp.spmatrix,
    Z_list: list[sp.spmatrix],
    Ginv: sp.spmatrix,
    sigma2_e: float,
    y: np.ndarray,
    need_inverse: bool,
    exact_limit: int = EXACT_LIMIT,
) -> MMESolution:
    X = sp.csr_matrix(X)
    W = sp.hstack([X] + [sp.csr_matrix(Z) for Z in Z_list], format="csr")
    p = X.shape[1]
    sizes = [Z.shape[1] for Z in Z_list]
    q = sum(sizes)
    if Ginv.shape != (q, q):
        raise ValueError(f"Ginv de dimension {Ginv.shape}, attendu {(q, q)}")
    pad = sp.block_diag([sp.csr_matrix((p, p)), sigma2_e * sp.csr_matrix(Ginv)], format="csc")
    LHS = (W.T @ W).tocsc() + pad
    rhs = W.T @ y
    lu = spla.splu(LHS, permc_spec="MMD_AT_PLUS_A", options={"SymmetricMode": True})
    sol = lu.solve(rhs)
    n_eq = p + q
    exact = n_eq <= exact_limit
    cinv = None
    pev = None
    if need_inverse and exact:
        cinv = sigma2_e * np.linalg.inv(LHS.toarray())
        d = np.diag(cinv)
        pev, off = [], p
        for k in sizes:
            pev.append(d[off:off + k].copy())
            off += k
    u, off = [], p
    for k in sizes:
        u.append(sol[off:off + k])
        off += k
    return MMESolution(
        beta=sol[:p], u=u, pev=pev, cinv=cinv, solve=lu.solve, n_equations=n_eq, n_fixed=p,
        sizes=sizes, exact=exact, ssr=float(y @ y - sol @ rhs),
    )
