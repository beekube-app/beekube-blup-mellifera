"""Approximation de la PEV au-delà de EXACT_LIMIT équations (esprit Misztal & Wiggans 1988).

Information propre de chaque groupe d'ouvrières (notes après absorption des groupes de
contemporaines et de l'environnement permanent), puis propagation dans le pedigree généralisé
(a_i = Σ c a_p + δ_i) : passe montante (descendants -> parents) et passe descendante
(parents -> descendants), en ignorant les boucles. Exact sur un pedigree en arbre.
"""
from __future__ import annotations

import numpy as np

from .kinship import generalized_parents, relationship_diagonal
from .pedigree import BeePedigree


def effective_records(factor_labels: dict[str, list], w_nodes: list[int], n_nodes: int) -> np.ndarray:
    """Nombre effectif de notes par nœud W après absorption des effets fixes (Misztal & Wiggans) :
    pour chaque facteur, n_i ← n_i · Σ_g (n_ig − n_ig² / n_g) / n_i, n_ig = notes de la colonie i dans
    le niveau g, n_g = notes du niveau g ; les facteurs sont absorbés successivement (produit)."""
    n_i = np.zeros(n_nodes)
    for w in w_nodes:
        n_i[w] += 1.0
    ratio = np.ones(n_nodes)
    for labels in factor_labels.values():
        n_g: dict = {}
        n_ig: dict = {}
        for g, w in zip(labels, w_nodes):
            n_g[g] = n_g.get(g, 0) + 1
            n_ig[(w, g)] = n_ig.get((w, g), 0) + 1
        kept = np.zeros(n_nodes)
        for (w, g), c in n_ig.items():
            kept[w] += c - c * c / n_g[g]
        with np.errstate(invalid="ignore", divide="ignore"):
            ratio *= np.where(n_i > 0, kept / np.maximum(n_i, 1e-12), 0.0)
    return n_i * ratio


def approx_pev(ped: BeePedigree, n_eff: np.ndarray, sigma2_A: float, sigma2_pe: float,
               sigma2_e: float) -> np.ndarray:
    parents, phi = generalized_parents(ped)
    diag = relationship_diagonal(ped)
    n = ped.n_nodes
    r_pe, r_e = sigma2_pe / sigma2_A, sigma2_e / sigma2_A
    U = np.where(n_eff > 0, 1.0 / (r_pe + r_e / np.maximum(n_eff, 1e-12)), 0.0)
    up: list[list[float]] = [[0.0] * len(par) for par in parents]
    prior_part = [[c * c * diag[p] for p, c in par] for par in parents]
    for i in range(n - 1, -1, -1):
        if U[i] <= 0:
            continue
        tot = sum(prior_part[i])
        for k, (p, c) in enumerate(parents[i]):
            other = phi[i] + tot - prior_part[i][k]
            m = c * c / (1.0 / U[i] + other)
            up[i][k] = m
            U[p] += m
    P = np.zeros(n)
    for i in range(n):
        pv = phi[i]
        for k, (p, c) in enumerate(parents[i]):
            pv += c * c / max(P[p] - up[i][k], 1e-12)
        P[i] = U[i] + 1.0 / pv
    return sigma2_A / P
