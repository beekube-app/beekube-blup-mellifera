"""Index base 100 (spec §3) : base mobile, direction, statuts, index total pondéré."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from .model_standard import CriterionFit
from .pedigree import BeePedigree
from .schemas import BaseInfo, BlupV2Input, CriterionIndex, CriterionIn, QueenOut

INDEX_POINTS_PER_SD = 10.0


@dataclass
class IndexTable:
    criteria: list[CriterionIn]
    status: dict[int, dict[int, str]]                       # critère -> reine -> statut
    index: dict[int, dict[int, Optional[float]]]            # index publié (nœud W)
    index_queen: dict[int, dict[int, Optional[float]]]      # index du nœud Q
    index_worker: dict[int, dict[int, Optional[float]]]     # breeder uniquement
    total: dict[int, Optional[float]] = field(default_factory=dict)
    total_reliability: dict[int, Optional[float]] = field(default_factory=dict)
    base: BaseInfo = field(default_factory=lambda: BaseInfo(birth_years=[], n_queens=0))


def base_birth_years(inp: BlupV2Input) -> list[int]:
    years = sorted({q.born for q in inp.queens})
    return years[-inp.options.base_years:]


def _scale(values: dict[int, float], ref: list[int], sigma: float, direction: int,
           keep: list[int]) -> dict[int, Optional[float]]:
    m = float(np.mean([values[q] for q in ref])) if ref else 0.0
    out: dict[int, Optional[float]] = {}
    for q in keep:
        out[q] = 100.0 if sigma <= 0 else 100.0 + INDEX_POINTS_PER_SD * direction * (values[q] - m) / sigma
    return out


def is_weighted(c: CriterionIn) -> bool:
    return (not c.elimination) and c.weight > 0


def compute_indexes(inp: BlupV2Input, fits: dict[int, Optional[CriterionFit]]) -> IndexTable:
    years = base_birth_years(inp)
    born = {q.id: q.born for q in inp.queens}
    qids = [q.id for q in inp.queens]
    table = IndexTable(criteria=list(inp.criteria), status={}, index={}, index_queen={}, index_worker={},
                       base=BaseInfo(birth_years=years, n_queens=sum(1 for q in qids if born[q] in years)))
    for c in inp.criteria:
        fit = fits.get(c.id)
        if fit is None:
            table.status[c.id] = {q: "failed" for q in qids}
            table.index[c.id] = {q: None for q in qids}
            table.index_queen[c.id] = {q: None for q in qids}
            table.index_worker[c.id] = {q: None for q in qids}
            continue
        st = {q: fit.status.get(q, "none") for q in qids}
        informed = [q for q in qids if st[q] != "none"]
        ref = [q for q in informed if born[q] in years] or informed
        table.status[c.id] = st
        idx = _scale(fit.ebv, ref, fit.sigma_index, c.direction, informed)
        idq = _scale(fit.ebv_queen, ref, fit.sigma_queen, c.direction, informed)
        idw = (_scale(fit.ebv_worker, ref, fit.sigma_worker or 0.0, c.direction, informed)
               if fit.ebv_worker is not None else {})
        table.index[c.id] = {q: idx.get(q) for q in qids}
        table.index_queen[c.id] = {q: idq.get(q) for q in qids}
        table.index_worker[c.id] = {q: idw.get(q) for q in qids}
    for q in qids:
        num = den = rnum = 0.0
        for c in inp.criteria:
            fit = fits.get(c.id)
            v = table.index[c.id][q]
            if fit is None or v is None or not is_weighted(c):
                continue
            num += c.weight * v
            rnum += c.weight * fit.reliability.get(q, 0.0)
            den += c.weight
        table.total[q] = num / den if den > 0 else None
        table.total_reliability[q] = rnum / den if den > 0 else None
    return table


def _r(v: Optional[float], nd: int) -> Optional[float]:
    return None if v is None else round(float(v), nd)


def queen_outputs(inp: BlupV2Input, ped: BeePedigree, fits: dict[int, Optional[CriterionFit]],
                  table: IndexTable, F) -> list[QueenOut]:
    out = []
    for q in inp.queens:
        indexes = {}
        for c in inp.criteria:
            fit = fits.get(c.id)
            st = table.status[c.id][q.id]
            informed = fit is not None and st != "none"
            indexes[str(c.id)] = CriterionIndex(
                status=st,
                ebv=_r(fit.ebv[q.id], 6) if informed else None,
                index=_r(table.index[c.id][q.id], 1),
                reliability=_r(fit.reliability[q.id], 3) if informed else None,
                index_queen=_r(table.index_queen[c.id][q.id], 1),
                index_worker=_r(table.index_worker[c.id][q.id], 1),
                n_records=fit.n_records.get(q.id, 0) if fit else 0,
                n_evaluators=fit.n_evaluators.get(q.id, 0) if fit else 0,
                n_years=fit.n_years.get(q.id, 0) if fit else 0,
            )
        out.append(QueenOut(
            id=q.id,
            inbreeding=round(float(F[ped.q_node(q.id)]), 4),
            inbreeding_workers=round(float(F[ped.w_node(q.id)]), 4),
            total_index=_r(table.total[q.id], 1),
            total_reliability=_r(table.total_reliability[q.id], 3),
            indexes=indexes,
        ))
    return out
