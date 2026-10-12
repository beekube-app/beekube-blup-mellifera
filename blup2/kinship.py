"""Parenté abeille (Brascamp & Bijma 2014) en représentation creuse.

Port de docs/analyses/bee_kinship.py. La D bloc-diagonale de BB (covariance cov(δ_FS) entre
pleines sœurs d'un même couple mère × sire) est rendue diagonale en prenant le groupe
d'ouvrières W(mère) comme « parent de bloc » : W(d) = ½ a_d + ½ a_s + η, var(η) = cov(δ_FS),
et chaque fille (ou nœud S) du couple = W(d) + ε, var(ε) = var(δ) − cov(δ_FS). La loi jointe est
identique à celle de BB ; A = T Φ Tᵀ avec T = (I − C)⁻¹ et Φ diagonale, d'où
A⁻¹ = (I − C)ᵀ Φ⁻¹ (I − C), creuse.

Sire inconnu (fécondation naturelle / inconnue, BB annexe 1 cas 3) : la reine reçoit
var(δ) = ¾ − ¼ F_mère ; le groupe d'ouvrières reçoit var(δ) = 1/(2D) (moyenne de D gamètes
paternels issus de colonies non apparentées : limite S → ∞ de l'éq. 21a), au lieu de 0 dans
bee_kinship.py (matrice singulière).
"""
from __future__ import annotations

from typing import Optional

import numpy as np
import scipy.sparse as sp

from .pedigree import BeePedigree

_PHI_MIN = 1e-10


def _dot(r1: dict, r2: dict, phi: np.ndarray) -> float:
    if len(r1) > len(r2):
        r1, r2 = r2, r1
    return float(sum(v * r2[k] * phi[k] for k, v in r1.items() if k in r2))


def _p12(S: int, D: int) -> tuple[float, float]:
    p1 = 1.0 / D
    p2 = min(1.0, 1.0 / D + 1.0 / S)  # hypothèse de Poisson (éq. 22a, 23a)
    return p1, p2


def _build(ped: BeePedigree) -> dict:
    if "core" in ped.cache:
        return ped.cache["core"]
    n = ped.n_nodes
    T: list[dict] = [None] * n  # type: ignore[list-item]
    phi = np.zeros(n)
    diag = np.zeros(n)
    F = np.zeros(n)
    parents: list[list[tuple[int, float]]] = [None] * n  # type: ignore[list-item]
    sire_attr: dict[int, tuple[float, float, int]] = {}  # nœud sire -> (Fs, ass, S)
    block: dict[int, Optional[tuple[float, float]]] = {}  # nœud Q mère -> (var_q, covFS)

    def a(i: int, j: int) -> float:
        return _dot(T[i], T[j], phi)

    for i in range(n):
        kind, d, s = ped.kind[i], ped.dam[i], ped.sire[i]
        if kind == "W":
            Fd = diag[d] - 1.0
            if s is not None:
                Fs, ass, S = sire_attr[s]
                var_q = 0.25 * (1 - Fd) + 0.25 * (1 - Fs) + 0.25 * (S - 1) / S * (1 + Fs - ass)
                p1, p2 = _p12(S, ped.D[i])
                cov = p1 * (1 - Fs) / 4.0 + (p2 - 1.0 / S) * (1 + Fs - ass) / 4.0  # éq. 21
                block[d] = (var_q, cov)
                par, ph = [(d, 0.5), (s, 0.5)], cov
            else:
                block[d] = None
                par, ph = [(d, 0.5)], 1.0 / (2.0 * ped.D[i])
        else:
            nS = ped.S[i] if kind == "S" else 1
            if d is None:
                par, ph = [], 1.0 / nS  # cas 1 (reine) et 2 (sire de base)
            elif s is not None:
                var_q, cov = block[d]  # type: ignore[misc]
                par, ph = [(ped.w_node(ped.queen_of_node(d)), 1.0)], (var_q - cov) / nS
            else:
                Fd = diag[d] - 1.0
                par, ph = [(d, 0.5)], (0.75 - 0.25 * Fd) / nS  # cas 3 et 6
        if ph < -1e-12:
            raise ValueError(f"variance de ségrégation négative au nœud {i} ({ph:.3g}) : pedigree incohérent")
        ph = max(ph, _PHI_MIN)
        row = {i: 1.0}
        for p, c in par:
            for k, v in T[p].items():
                row[k] = row.get(k, 0.0) + c * v
        T[i], phi[i], parents[i] = row, ph, par
        diag[i] = sum(v * v * phi[k] for k, v in row.items())
        if d is not None and s is not None:
            F[i] = 0.5 * a(d, s)  # éq. 25
        if kind == "S":
            if d is not None and s is not None:
                cov = block[d][1]  # type: ignore[index]
                ass = 0.25 * diag[d] + 0.25 * diag[s] + 0.5 * a(d, s) + cov
            elif d is not None:
                ass = 0.25 * diag[d]
            else:
                ass = 0.0
            sire_attr[i] = (F[i], ass, ped.S[i])
        elif kind == "Q":
            sire_attr[i] = (F[i], diag[i], 1)
    core = {"T": T, "phi": phi, "diag": diag, "F": F, "parents": parents}
    ped.cache["core"] = core
    return core


def generalized_parents(ped: BeePedigree) -> tuple[list[list[tuple[int, float]]], np.ndarray]:
    """(parents[i] = [(nœud, coefficient)], variances de ségrégation φ) : a_i = Σ c a_p + δ_i."""
    core = _build(ped)
    return core["parents"], core["phi"]


def relationship_diagonal(ped: BeePedigree) -> np.ndarray:
    return _build(ped)["diag"]


def inbreeding(ped: BeePedigree) -> np.ndarray:
    """F de chaque nœud = ½ a(mère, sire) (éq. 25) : F de la reine (Q), des ouvrières (W),
    des RPM (S). 0 si un parent est inconnu."""
    return _build(ped)["F"]


def relationship(ped: BeePedigree, i: int, j: int) -> float:
    core = _build(ped)
    return _dot(core["T"][i], core["T"][j], core["phi"])


def _factor_rows(ped: BeePedigree, nodes) -> sp.csr_matrix:
    core = _build(ped)
    T, sq = core["T"], np.sqrt(core["phi"])
    rows, cols, vals = [], [], []
    for r, i in enumerate(nodes):
        for k, v in T[i].items():
            rows.append(r)
            cols.append(k)
            vals.append(v * sq[k])
    return sp.csr_matrix((vals, (rows, cols)), shape=(len(nodes), ped.n_nodes))


def relationship_matrix(ped: BeePedigree, nodes=None) -> np.ndarray:
    """A dense, sur tous les nœuds (nodes=None) ou sur la liste de nœuds donnée."""
    if nodes is None:
        nodes = range(ped.n_nodes)
    L = _factor_rows(ped, list(nodes))
    return np.asarray((L @ L.T).todense())


def inverse_relationship(ped: BeePedigree) -> sp.csr_matrix:
    if "ainv" in ped.cache:
        return ped.cache["ainv"]
    parents, phi = generalized_parents(ped)
    n = ped.n_nodes
    rows, cols, vals = list(range(n)), list(range(n)), [1.0] * n
    for i, par in enumerate(parents):
        for p, c in par:
            rows.append(i)
            cols.append(p)
            vals.append(-c)
    IC = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
    Ainv = (IC.T @ sp.diags(1.0 / phi) @ IC).tocsr()
    ped.cache["ainv"] = Ainv
    return Ainv


def hypothetical_sire_relationship(
    ped: BeePedigree, A: Optional[np.ndarray], mother_queen_id: int, drone_queen_id: int, n_drones: int
) -> float:
    """a(Q_M, sire_P) pour un accouplement planifié M × P avec le sire hypothétique « P, S = 1, D ».

    Avec S = 1 et P comme RPM, le nœud sire est Q(P) lui-même (même règle que build_pedigree
    pour drone_queens = [P]) : a(Q_M, sire_P) = a(Q_M, Q_P). D n'intervient que dans la
    covariance entre pleines sœurs, pas dans cette parenté (paramètre gardé pour l'interface).
    A : matrice complète des nœuds (n_nodes × n_nodes) ou None (calcul creux à la demande).
    """
    i, j = ped.q_node(mother_queen_id), ped.q_node(drone_queen_id)
    if A is not None:
        return float(A[i, j])
    return relationship(ped, i, j)
