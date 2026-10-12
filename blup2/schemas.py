"""Contrat d'API v2 (spec §5) : modèles pydantic d'entrée et de sortie."""
from __future__ import annotations

import datetime as dt
import re
from collections import Counter
from typing import Literal, Optional

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

BORN_MIN = 1900   # born : entre 1900 et l'année courante + 1 (spec §5.1)


class _In(BaseModel):
    model_config = ConfigDict(extra="forbid", allow_inf_nan=False)


class OptionsIn(_In):
    base_years: int = Field(5, ge=1, le=50)
    estimate_variances: Literal["auto", "always", "never"] = "auto"
    default_n_drones: int = Field(12, ge=1, le=100)
    inbreeding_penalty: float = Field(100.0, ge=0, le=1e6)
    max_inbreeding: float = Field(0.0625, ge=0, le=1)
    n_suggestions: int = Field(10, ge=0, le=1000)


class CriterionIn(_In):
    id: int
    name: str = ""
    direction: Literal[1, -1] = 1
    kind: Literal["note", "percentage", "harvest"] = "note"
    weight: float = Field(1.0, ge=0, le=1e6)
    elimination: bool = False
    h2: Optional[float] = Field(None, gt=0, lt=1)
    repeatability: Optional[float] = Field(None, gt=0, lt=1)
    queen_ratio: Optional[float] = Field(None, gt=0)
    rg: Optional[float] = Field(None, ge=-1, le=1)

    @model_validator(mode="after")
    def _repeatability_ge_h2(self) -> "CriterionIn":
        if self.h2 is not None and self.repeatability is not None and self.repeatability < self.h2:
            raise ValueError(f"critère {self.id} : repeatability ({self.repeatability}) < h2 ({self.h2})")
        return self


class MatingIn(_In):
    type: Literal["natural", "station", "insemination", "unknown"] = "unknown"
    drone_queens: list[int] = Field(default_factory=list)
    drone_queens_mother: Optional[int] = None
    n_drone_queens: Optional[int] = Field(None, ge=1)
    n_drones: Optional[int] = Field(None, ge=1)
    external_ref: Optional[str] = None
    station: Optional[int] = None
    year: Optional[int] = Field(None, ge=1000, le=9999)


class QueenIn(_In):
    id: int
    born: int
    mother: Optional[int] = None
    mother_ref_ext: Optional[str] = None
    mating: Optional[MatingIn] = None
    variety: Optional[str] = None
    died: bool = False
    eliminated: bool = False
    available_as_mother: bool = True
    available_as_drone_source: bool = True

    @field_validator("born", mode="before")
    @classmethod
    def _born_year(cls, v):
        if isinstance(v, bool) or not isinstance(v, int):
            raise ValueError("born doit être une année (entier)")
        last = dt.date.today().year + 1
        if not BORN_MIN <= v <= last:
            raise ValueError(f"born ({v}) doit être compris entre {BORN_MIN} et {last}")
        return v


class RecordIn(_In):
    queen: int
    criterion: int
    value: float
    scale_min: float = 0.0
    scale_max: float = 1.0
    date: dt.date
    evaluator: Optional[int] = None
    apiary: int
    hive_type: Optional[str] = None
    hive_class: Literal["hive", "nuc", "split"] = "hive"
    visit: Optional[int] = None

    @field_validator("date", mode="before")
    @classmethod
    def _iso_date(cls, v):
        if isinstance(v, dt.date) and not isinstance(v, dt.datetime):
            return v                      # relecture d'un modèle déjà validé (model_dump)
        if not isinstance(v, str) or not re.fullmatch(r"\d{4}-\d{2}-\d{2}", v):
            raise ValueError("date au format ISO YYYY-MM-DD attendue")
        return dt.date.fromisoformat(v)

    @model_validator(mode="after")
    def _value_in_scale(self) -> "RecordIn":
        if not self.scale_max > self.scale_min:
            raise ValueError(f"scale_max ({self.scale_max}) doit être > scale_min ({self.scale_min})")
        if not (self.scale_min <= self.value <= self.scale_max):
            raise ValueError(
                f"note {self.value} hors de l'échelle [{self.scale_min}, {self.scale_max}] "
                f"(reine {self.queen}, critère {self.criterion})"
            )
        return self


class BlupV2Input(_In):
    exploitation: int
    profile: Literal["standard", "breeder"] = "standard"
    options: OptionsIn = Field(default_factory=OptionsIn)
    criteria: list[CriterionIn]
    queens: list[QueenIn]
    records: list[RecordIn] = Field(default_factory=list)

    @model_validator(mode="after")
    def _cross_checks(self):
        errors: list[str] = []
        qids = [q.id for q in self.queens]
        cids = [c.id for c in self.criteria]
        dq = sorted(i for i, n in Counter(qids).items() if n > 1)
        if dq:
            errors.append(f"identifiants de reines dupliqués : {dq}")
        dc = sorted(i for i, n in Counter(cids).items() if n > 1)
        if dc:
            errors.append(f"identifiants de critères dupliqués : {dc}")
        qset, cmap = set(qids), {c.id: c for c in self.criteria}
        unknown_q = sorted({r.queen for r in self.records if r.queen not in qset})
        if unknown_q:
            errors.append(f"records[].queen absent de queens : {unknown_q}")
        unknown_c = sorted({r.criterion for r in self.records if r.criterion not in cmap})
        if unknown_c:
            errors.append(f"records[].criterion absent de criteria : {unknown_c}")
        for r in self.records:
            c = cmap.get(r.criterion)
            if c is not None and c.kind != "note" and not (r.scale_min == 0 and r.scale_max == 1):
                errors.append(
                    f"critère {c.id} de type {c.kind} : scale_min = 0 et scale_max = 1 attendus "
                    f"(reine {r.queen})"
                )
                break
        if errors:
            raise ValueError(" ; ".join(errors))
        from .pedigree import find_cycle, listed_ids

        cycle = find_cycle(self.queens)
        if cycle:
            raise ValueError(f"pedigree cyclique : reines concernées {listed_ids(cycle)}")
        return self


class MatingPairIn(_In):
    mother: int
    drone_source: int


class MatingV2Input(BlupV2Input):
    pairs: list[MatingPairIn] = Field(default_factory=list)

    @model_validator(mode="after")
    def _pairs_known(self):
        qset = {q.id for q in self.queens}
        bad = sorted({i for p in self.pairs for i in (p.mother, p.drone_source) if i not in qset})
        if bad:
            raise ValueError(f"pairs : reines absentes de queens : {bad}")
        same = sorted({p.mother for p in self.pairs if p.mother == p.drone_source})
        if same:
            raise ValueError(f"pairs : une reine ne peut pas être accouplée avec ses propres mâles : {same}")
        return self


# ---------------------------------------------------------------- sorties
class CriterionParameters(BaseModel):
    source: Literal["input", "reml", "prior"]
    h2: float
    repeatability: float
    sigma2_a: float
    sigma2_pe: float
    sigma2_e: float
    sigma2_aq: Optional[float] = None
    rg: Optional[float] = None
    profile_used: Literal["standard", "breeder"] = "standard"
    known_mating_share: Optional[float] = None   # profil breeder : part des reines notées à accouplement connu
    n_records: int
    n_queens: int
    n_contemporary_groups: int
    reliability_method: Literal["exact", "approx", "none"]


class CriterionIndex(BaseModel):
    status: Literal["evaluated", "pedigree", "none", "failed"]
    ebv: Optional[float] = None
    index: Optional[float] = None
    reliability: Optional[float] = None
    index_queen: Optional[float] = None
    index_worker: Optional[float] = None
    n_records: int = 0
    n_evaluators: int = 0
    n_years: int = 0


class QueenOut(BaseModel):
    id: int
    inbreeding: float
    inbreeding_workers: float
    total_index: Optional[float] = None
    total_reliability: Optional[float] = None
    indexes: dict[str, CriterionIndex]


class MatingOut(BaseModel):
    mother: int
    drone_source: int
    expected_total_index: Optional[float] = None
    expected_indexes: dict[str, Optional[float]]
    expected_inbreeding: float
    score: float
    warnings: list[str] = Field(default_factory=list)


class TrendPoint(BaseModel):
    birth_year: int
    mean_index: float
    n: int


class Connectivity(BaseModel):
    n_components: int
    minor_share: float = 0.0     # part des reines notées hors de la composante qui en compte le plus
    components: list[list[int]]


class WarningOut(BaseModel):
    code: Literal[
        "disconnected_groups",
        "breeder_fallback",
        "few_records",
        "reml_not_converged",
        "low_heritability",
        "scale_mixed",
        "unknown_parent_added",
        "self_reference_removed",
        "no_variation",
        "reliability_approx",
        "criterion_failed",
    ]
    message: str
    criterion: Optional[int] = None


class BaseInfo(BaseModel):
    birth_years: list[int]
    n_queens: int


class BlupV2Output(BaseModel):
    status: Literal["success"] = "success"
    engine_version: str
    profile: Literal["standard", "breeder"]
    parameters: dict[str, CriterionParameters]
    base: BaseInfo
    queens: list[QueenOut]
    matings: list[MatingOut]
    trend: dict[str, list[TrendPoint]]
    connectivity: Connectivity
    warnings: list[WarningOut]


class MatingV2Output(BaseModel):
    status: Literal["success"] = "success"
    engine_version: str
    profile: Literal["standard", "breeder"]
    matings: list[MatingOut]
    warnings: list[WarningOut]


class ErrorV2Out(BaseModel):
    status: Literal["error"] = "error"
    errors: Optional[list[dict]] = None
    message: Optional[str] = None
