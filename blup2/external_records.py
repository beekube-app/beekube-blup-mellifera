"""Notes Index Mellifera des ancêtres résolus et de leurs contemporaines (spec volet 3 §4), lues dans la base
du crawler.

Un BLUP n'apprend d'une note que par contraste avec les autres notes de son groupe (testeur × année) : la note
d'un ancêtre seul dans son groupe n'apporte rien. On lit donc aussi les notes des autres reines du même groupe
(même testeur, ou même éleveur pour le radar, même année, même caractère) : ces contemporaines entrent au
pedigree comme les ancêtres, sans jamais sortir. Inclure les reines parmi lesquelles les reproductrices ont été
choisies rend aussi l'estimation peu sensible à cette sélection.

Seuls les critères qui déclarent un caractère standard (`CriterionIn.trait`) ayant un nom Index Mellifera, et
qui ont au moins une note de l'exploitation, reçoivent des notes. Nettoyage comme `mellifera/export.py::_records`
(sans en importer le code : le moteur ne dépend pas du paquet mellifera), chaque motif d'écart compté. Les notes
retenues sont ramenées sur [0, 1] dans le sens du critère ; rucher et évaluateur = un entier négatif par
testeur (ou éleveur, pour le radar).
"""
from __future__ import annotations

import datetime as dt
import json
import math
import re
import sqlite3
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Optional, Union

from .schemas import BORN_MIN, BlupV2Input, RecordIn
from .traits import STANDARD_TRAITS

MAX_CONTEMPORARIES = 1000   # borne le coût : environ quatre équations par contemporaine notée
DISCARD_LABELS = {          # motif -> libellé de l'avertissement external_records_added, dans cet ordre
    "harvest": "en kg",
    "no_date": "sans date",
    "out_of_scale": "hors échelle",
    "non_numeric": "non numériques",
    "no_apiary": "sans rucher",
    "before_birth": "antérieures à la naissance",
    "no_queen": "sans fiche de reine",
    "alone": "sans contemporaine",
}
_SELECT = (
    "SELECT e.queen_code, e.source, e.tester_breeder_id, e.eval_date, e.value, t.name AS trait, t.kind, "
    "t.scale_min, t.scale_max, t.direction, q.breeder_id, q.birth_year "
    "FROM evaluation e JOIN criterion c ON c.mode = e.mode AND c.code = e.criterion_code "
    "JOIN trait t ON t.id = c.trait_id LEFT JOIN queen q ON q.code = e.queen_code "
    "WHERE e.source IN ('season_table', 'radar') AND t.name IN (SELECT value FROM json_each(?)) AND ")
_ANCESTORS = _SELECT + "e.queen_code IN (SELECT value FROM json_each(?)) ORDER BY e.id"
_GROUPS = _SELECT + (
    "((e.source = 'season_table' AND e.tester_breeder_id IN (SELECT value FROM json_each(?))) "
    "OR (e.source = 'radar' AND q.breeder_id IN (SELECT value FROM json_each(?)))) ORDER BY e.id")
_ISO = re.compile(r"\d{4}-\d{2}-\d{2}")

Cell = tuple            # (code, critère, date, source, rucher Index Mellifera)


@dataclass
class Notes:
    cells: dict[Cell, list[float]] = field(default_factory=dict)    # valeurs sur [0, 1], sens du critère
    contemporaries: list[str] = field(default_factory=list)          # codes, ordre croissant
    discarded: Counter = field(default_factory=Counter)
    capped: bool = False


def wanted_traits(inp: BlupV2Input) -> dict[str, int]:
    """Nom Index Mellifera -> critère de l'exploitation qui le porte (critères notés dans l'exploitation)."""
    noted = {r.criterion for r in inp.records}
    out = {}
    for c in inp.criteria:
        name = STANDARD_TRAITS[c.trait].mellifera if c.trait is not None else None
        if name is not None and c.id in noted:
            out[name] = c.id
    return out


def _number(v) -> Optional[float]:
    if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
        return None
    return float(v)


def _date(v) -> Optional[dt.date]:
    if not isinstance(v, str) or not _ISO.fullmatch(v):
        return None
    try:
        return dt.date.fromisoformat(v)
    except ValueError:
        return None


def _valid_year(y) -> bool:
    return isinstance(y, int) and not isinstance(y, bool) and BORN_MIN <= y <= dt.date.today().year + 1


def _clean(e, wanted: dict[str, int], direction: dict[int, int]) -> Union[str, tuple]:
    """Motif d'écart, ou (critère, date, source, rucher, valeur sur [0, 1] dans le sens du critère)."""
    if e["kind"] == "harvest":
        return "harvest"
    value = _number(e["value"])
    if value is None:
        return "non_numeric"
    date = _date(e["eval_date"])
    if date is None:
        return "no_date"
    if e["kind"] == "percentage":
        value /= 100
    if e["kind"] == "note":
        lo, hi = _number(e["scale_min"]), _number(e["scale_max"])
        if lo is None or hi is None or not hi > lo:
            return "out_of_scale"
    else:
        lo, hi = 0.0, 1.0
    if not lo <= value <= hi:
        return "out_of_scale"
    apiary = e["tester_breeder_id"] if e["source"] == "season_table" else e["breeder_id"]
    if isinstance(apiary, bool) or not isinstance(apiary, int):
        return "no_apiary"
    v = (value - lo) / (hi - lo)
    cid = wanted[e["trait"]]
    if (-1 if e["direction"] == -1 else 1) != direction[cid]:
        v = 1.0 - v                     # note du site dans le sens contraire du critère de l'exploitation
    return cid, date, e["source"], apiary, v


def read_notes(conn: sqlite3.Connection, inp: BlupV2Input, ancestors: dict[str, int]) -> Notes:
    """ancestors : code Index Mellifera -> année de naissance (ancêtres résolus)."""
    out = Notes()
    wanted = wanted_traits(inp)
    if not wanted or not ancestors:
        return out
    direction = {c.id: c.direction for c in inp.criteria}
    names = json.dumps(sorted(wanted))
    cells: dict[Cell, list[float]] = defaultdict(list)
    for e in conn.execute(_ANCESTORS, (names, json.dumps(sorted(ancestors)))):
        res = _clean(e, wanted, direction)
        if isinstance(res, tuple) and res[1].year < ancestors[e["queen_code"]]:
            res = "before_birth"
        if isinstance(res, str):
            out.discarded[res] += 1
            continue
        cells[(e["queen_code"], *res[:4])].append(res[4])
    groups = {(k[1], k[3], k[4], k[2].year) for k in cells}         # (critère, source, rucher, année)
    if groups:
        testers = sorted({g[2] for g in groups if g[1] == "season_table"})
        breeders = sorted({g[2] for g in groups if g[1] == "radar"})
        for e in conn.execute(_GROUPS, (names, json.dumps(testers), json.dumps(breeders))):
            code = e["queen_code"]
            if code in ancestors:
                continue
            res = _clean(e, wanted, direction)
            if isinstance(res, str) or (res[0], res[2], res[3], res[1].year) not in groups:
                continue                # hors des groupes des ancêtres : ignorée sans être comptée
            if code is None or not _valid_year(e["birth_year"]):
                out.discarded["no_queen"] += 1
            elif res[1].year < e["birth_year"]:
                out.discarded["before_birth"] += 1
            else:
                cells[(code, *res[:4])].append(res[4])
    season = {(k[0], k[1], k[2].year) for k in cells if k[3] == "season_table"}
    cells = {k: v for k, v in cells.items()       # le tableau de saison couvre la même reine, le même critère, l'année
             if not (k[3] == "radar" and (k[0], k[1], k[2].year) in season)}
    others = sorted({k[0] for k in cells if k[0] not in ancestors})
    if len(others) > MAX_CONTEMPORARIES:
        out.capped = True
        kept = set(others[:MAX_CONTEMPORARIES])
        cells = {k: v for k, v in cells.items() if k[0] in ancestors or k[0] in kept}
    members: dict[tuple, set] = defaultdict(set)
    for k in cells:
        members[(k[1], k[3], k[4], k[2].year)].add(k[0])
    live = {g for g, qs in members.items() if len(qs) >= 2 and qs & ancestors.keys()}
    for k in list(cells):
        if (k[1], k[3], k[4], k[2].year) not in live:
            if k[0] in ancestors:
                out.discarded["alone"] += 1
            del cells[k]
    out.cells = cells
    out.contemporaries = sorted({k[0] for k in cells if k[0] not in ancestors})
    return out


def build_records(notes: Notes, ids: dict[str, int], inp: BlupV2Input) -> list[RecordIn]:
    """Un RecordIn par cellule (valeurs moyennées). ids : code -> id négatif (ancêtres et contemporaines).
    Rucher = évaluateur = un entier négatif par (source, testeur ou éleveur), sous ceux de l'exploitation et
    sous -1 (évaluateur inconnu dans records.py)."""
    used = [r.apiary for r in inp.records] + [r.evaluator for r in inp.records if r.evaluator is not None]
    start = min(-1, min(used, default=-1))
    apiaries = {a: start - k for k, a in enumerate(sorted({(k[3], k[4]) for k in notes.cells}), start=1)}
    records = []
    for (code, cid, date, source, apiary), values in sorted(notes.cells.items()):
        x = apiaries[(source, apiary)]
        mean = min(1.0, max(0.0, sum(values) / len(values)))
        records.append(RecordIn(queen=ids[code], criterion=cid, value=mean, date=date, apiary=x, evaluator=x))
    return records


def records_message(records: list[RecordIn], ancestor_ids: set[int], notes: Notes) -> Optional[str]:
    """Message de external_records_added ; None si rien n'a été ajouté ni écarté."""
    if not records and not notes.discarded:
        return None
    queens = {r.queen for r in records}
    msg = (f"{len(records)} note(s) Index Mellifera ajoutée(s) sur {len({r.criterion for r in records})} critère(s) : "
           f"{len(queens & ancestor_ids)} ancêtre(s) et {len(queens - ancestor_ids)} contemporaine(s)")
    parts = [f"{notes.discarded[k]} {label}" for k, label in DISCARD_LABELS.items() if notes.discarded[k]]
    if parts:
        msg += f" ; écartées : {', '.join(parts)}"
    if notes.capped:
        msg += f" ; plafond de {MAX_CONTEMPORARIES} contemporaines atteint"
    return msg
