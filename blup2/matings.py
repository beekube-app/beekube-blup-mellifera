"""Accouplements (spec §4) : index attendu de la descendance, F attendu, score."""
from __future__ import annotations

from typing import Optional

import numpy as np

from .indexes import IndexTable, is_weighted
from .kinship import relationship_matrix
from .pedigree import BeePedigree
from .schemas import BlupV2Input, MatingOut, MatingPairIn


def _queen_vectors(table: IndexTable, ids: list[int]) -> tuple[dict[int, np.ndarray], np.ndarray, bool]:
    """index_queen par critère (null -> 100) et index total reine attendu t (moyenne pondérée)."""
    vec = {c.id: np.array([v if (v := table.index_queen[c.id][q]) is not None else 100.0 for q in ids])
           for c in table.criteria}
    weighted = [c for c in table.criteria if is_weighted(c)]
    den = sum(c.weight for c in weighted)
    if den <= 0:
        return vec, np.full(len(ids), 100.0), False
    t = sum(c.weight * vec[c.id] for c in weighted) / den
    return vec, t, True


def _expected_index(table: IndexTable, cid: int, m: int, p: int) -> Optional[float]:
    im, ip = table.index_queen[cid][m], table.index_queen[cid][p]
    if im is None and ip is None:
        return None
    return 100.0 + 0.5 * ((im if im is not None else 100.0) - 100.0) + 0.5 * ((ip if ip is not None else 100.0) - 100.0)


def _mating(inp: BlupV2Input, table: IndexTable, m: int, p: int, total: Optional[float], f_exp: float) -> MatingOut:
    opts = inp.options
    score = (total if total is not None else 100.0) - opts.inbreeding_penalty * f_exp
    return MatingOut(
        mother=m,
        drone_source=p,
        expected_total_index=None if total is None else round(float(total), 1),
        expected_indexes={str(c.id): (None if (v := _expected_index(table, c.id, m, p)) is None else round(v, 1))
                          for c in table.criteria},
        expected_inbreeding=round(float(f_exp), 4),
        score=round(float(score), 1),
        warnings=["high_inbreeding"] if f_exp >= opts.max_inbreeding else [],
    )


def suggest_matings(inp: BlupV2Input, ped: BeePedigree, table: IndexTable) -> list[MatingOut]:
    n_sug = inp.options.n_suggestions
    mothers = [q.id for q in inp.queens if not q.died and not q.eliminated and q.available_as_mother]
    sources = [q.id for q in inp.queens if not q.died and not q.eliminated and q.available_as_drone_source]
    if n_sug == 0 or not mothers or not sources:
        return []
    ids = list(dict.fromkeys(mothers + sources))
    pos = {q: k for k, q in enumerate(ids)}
    A = relationship_matrix(ped, [ped.q_node(q) for q in ids])
    _, t, has_total = _queen_vectors(table, ids)
    im = np.array([pos[m] for m in mothers])
    ip = np.array([pos[p] for p in sources])
    total = 0.5 * t[im][:, None] + 0.5 * t[ip][None, :]
    F = 0.5 * A[np.ix_(im, ip)]
    score = total - inp.options.inbreeding_penalty * F
    same = np.array(mothers)[:, None] == np.array(sources)[None, :]
    score[same] = -np.inf
    out: list[MatingOut] = []
    k = min(n_sug, len(sources))
    for r, m in enumerate(mothers):
        best = np.argsort(-score[r], kind="stable")[:k]
        for c in best:
            if not np.isfinite(score[r, c]):
                continue
            out.append(_mating(inp, table, m, sources[c], float(total[r, c]) if has_total else None, float(F[r, c])))
    return out


def evaluate_pairs(inp: BlupV2Input, ped: BeePedigree, table: IndexTable, pairs: list[MatingPairIn]) -> list[MatingOut]:
    if not pairs:
        return []
    ids = list(dict.fromkeys([x for pr in pairs for x in (pr.mother, pr.drone_source)]))
    pos = {q: k for k, q in enumerate(ids)}
    A = relationship_matrix(ped, [ped.q_node(q) for q in ids])
    _, t, has_total = _queen_vectors(table, ids)
    out = []
    for pr in pairs:
        i, j = pos[pr.mother], pos[pr.drone_source]
        total = 0.5 * t[i] + 0.5 * t[j] if has_total else None
        out.append(_mating(inp, table, pr.mother, pr.drone_source, total, 0.5 * float(A[i, j])))
    return out
