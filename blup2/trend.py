"""Tendance génétique (spec §3) : moyenne des index par année de naissance."""
from __future__ import annotations

import numpy as np

from .indexes import IndexTable
from .schemas import BlupV2Input, TrendPoint


def genetic_trend(inp: BlupV2Input, table: IndexTable) -> dict[str, list[TrendPoint]]:
    out: dict[str, list[TrendPoint]] = {}
    for c in inp.criteria:
        by_year: dict[int, list[float]] = {}
        for q in inp.queens:
            v = table.index[c.id][q.id]
            if v is not None:
                by_year.setdefault(q.born, []).append(v)
        out[str(c.id)] = [TrendPoint(birth_year=y, mean_index=round(float(np.mean(v)), 1), n=len(v))
                          for y, v in sorted(by_year.items())]
    return out
