"""Connexion (spec §3) : composantes du graphe groupe de contemporaines — reine — apparentées.

minor_share = part des reines notées (au moins une note) hors de la composante qui en compte le
plus, arrondie à 3 décimales (0 s'il n'y a qu'une composante).
"""
from __future__ import annotations

import numpy as np

from .kinship import relationship_matrix
from .pedigree import BeePedigree
from .schemas import BlupV2Input, Connectivity

RELATED_MIN = 0.125


class _UnionFind:
    def __init__(self) -> None:
        self.parent: dict = {}

    def find(self, x):
        self.parent.setdefault(x, x)
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra


def connectivity(inp: BlupV2Input, ped: BeePedigree, threshold: float = RELATED_MIN) -> Connectivity:
    groups_of: dict[int, list[tuple[int, int]]] = {}
    for r in inp.records:
        g = (r.apiary, r.date.year)
        lst = groups_of.setdefault(r.queen, [])
        if g not in lst:
            lst.append(g)
    if not groups_of:
        return Connectivity(n_components=0, components=[])
    uf = _UnionFind()
    queens = list(groups_of)
    for q in queens:
        for g in groups_of[q]:
            uf.union(groups_of[q][0], g)
    A = relationship_matrix(ped, [ped.q_node(q) for q in queens])
    np.fill_diagonal(A, 0.0)
    rows, cols = np.nonzero(np.triu(A) >= threshold)
    for i, j in zip(rows, cols):
        uf.union(groups_of[queens[i]][0], groups_of[queens[j]][0])
    comps: dict = {}
    n_queens: dict = {}
    for q in queens:
        root = uf.find(groups_of[q][0])
        n_queens[root] = n_queens.get(root, 0) + 1
        for g in groups_of[q]:
            comps.setdefault(uf.find(g), set()).add(g[0])
    components = sorted(sorted(c) for c in comps.values())
    minor_share = round((len(queens) - max(n_queens.values())) / len(queens), 3)
    return Connectivity(n_components=len(components), minor_share=minor_share, components=components)
