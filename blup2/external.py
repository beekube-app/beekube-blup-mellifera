"""Références externes Index Mellifera (spec Index Mellifera §8) : mother_ref_ext et mating.external_ref
résolus dans la base du crawler (lecture seule), ascendance ajoutée au pedigree comme reines de base.

L'entrée n'est jamais modifiée : expand_external renvoie des copies des reines dont les liens sont
complétés et, à part, les ancêtres (id négatifs) que seul build_pedigree reçoit. Aucune donnée de la base
ne fait échouer le calcul : liens erronés coupés, valeurs inattendues nettoyées, base illisible ou
erreur imprévue = calcul sans résolution ; le nombre d'ancêtres est plafonné (MAX_ANCESTORS).
"""
from __future__ import annotations

import datetime as dt
import logging
import os
import sqlite3
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from .pedigree import find_cycle, listed_ids
from .schemas import BORN_MIN, BlupV2Input, MatingIn, QueenIn, WarningOut

logger = logging.getLogger(__name__)

ENV_DB = "MELLIFERA_DB"
MAX_ANCESTORS = 1000    # plafond du nombre d'ancêtres ajoutés : borne le coût de la parenté (relecture I2)
_MATING_TYPES = ("natural", "station", "insemination", "unknown")
_COLUMNS = ("code, name, birth_year, mother_code, drone_mother_code, n_drone_queens, single_drone, "
            "race_description, mating_type")


@dataclass
class Expansion:
    queens: list[QueenIn]                                        # reines de l'entrée, liens complétés
    ancestors: list[QueenIn] = field(default_factory=list)      # ancêtres Index Mellifera (id négatifs)
    warnings: list[WarningOut] = field(default_factory=list)


def normalize_code(ref: str) -> str:
    """Code international saisi : espaces retirés, majuscules."""
    return "".join(ref.split()).upper()


def _valid_year(y) -> bool:
    return isinstance(y, int) and BORN_MIN <= y <= dt.date.today().year + 1


def _mother_ref(q: QueenIn) -> Optional[str]:
    if q.mother is None and q.mother_ref_ext and q.mother_ref_ext.strip():
        return q.mother_ref_ext
    return None


def _drone_ref(q: QueenIn) -> Optional[str]:
    """Côté mâle externe : accouplement non naturel, sans RPM ni mère des RPM connues dans l'exploitation."""
    m = q.mating
    if (m is not None and m.type != "natural" and not m.drone_queens and m.drone_queens_mother is None
            and m.external_ref and m.external_ref.strip()):
        return m.external_ref
    return None


def _connect(path: str) -> sqlite3.Connection:
    if not Path(path).is_file():
        raise sqlite3.OperationalError(f"base introuvable : {path}")
    conn = sqlite3.connect(f"file:{Path(path).resolve()}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    return conn


class _Base:
    def __init__(self, conn: sqlite3.Connection):
        self.conn = conn
        self.by_code: dict[str, Optional[sqlite3.Row]] = {}

    def row(self, code: str) -> Optional[sqlite3.Row]:
        if code not in self.by_code:
            self.by_code[code] = self.conn.execute(
                f"SELECT {_COLUMNS} FROM queen_mating WHERE code = ?", (code,)).fetchone()
        return self.by_code[code]

    def find(self, ref: str) -> Optional[sqlite3.Row]:
        """Code exact après normalisation (ou tel que saisi), à défaut nom affiché s'il est unique en base."""
        for code in dict.fromkeys((normalize_code(ref), ref.strip())):
            r = self.row(code)
            if r is not None:
                return r
        rows = self.conn.execute(f"SELECT {_COLUMNS} FROM queen_mating WHERE name = ? LIMIT 2",
                                 (ref.strip(),)).fetchall()
        return rows[0] if len(rows) == 1 else None


def _ancestry(base: _Base, roots: list[str], depth: int) -> tuple[dict[str, sqlite3.Row], int, bool]:
    """Ancêtres atteints en au plus `depth` générations (racines = génération 1), en largeur : une reine
    atteinte par plusieurs chemins n'est ajoutée qu'une fois ; la remontée s'arrête à MAX_ANCESTORS.
    Retourne (code → ligne, nombre sans année, plafond atteint)."""
    rows: dict[str, sqlite3.Row] = {}
    no_year: set[str] = set()
    frontier = list(roots)
    for _ in range(depth):
        nxt = []
        for code in frontier:
            if code in rows or code in no_year:
                continue
            if len(rows) >= MAX_ANCESTORS:
                return rows, len(no_year), True
            r = base.row(code)
            if r is None:
                continue
            if _valid_year(r["birth_year"]):
                rows[code] = r
                nxt += [p for p in (r["mother_code"], r["drone_mother_code"]) if p]
            else:
                no_year.add(code)
        frontier = nxt
    return rows, len(no_year), False


def _ancestor(r: sqlite3.Row, ids: dict[str, int]) -> tuple[QueenIn, int]:
    """QueenIn d'un ancêtre (accouplement construit comme à l'export) ; liens hors de l'ensemble coupés,
    auto-références retirées. Retourne (reine, nombre d'auto-références coupées)."""
    code, born = r["code"], r["birth_year"]
    mother, drone = r["mother_code"], r["drone_mother_code"]
    cut = 0
    if mother == code:
        mother, cut = None, cut + 1
    if drone == code:
        drone, cut = None, cut + 1
    mtype = r["mating_type"] if r["mating_type"] in _MATING_TYPES else "unknown"
    if drone and mtype == "unknown":
        mtype = "insemination"          # côté mâle connu : sinon le pedigree ignorerait l'accouplement
    mating: dict = {"type": mtype, "year": born}
    if drone:
        mating["external_ref"] = drone  # même nœud sire pour toutes les reines accouplées avec cette lignée
        if drone in ids:
            mating["drone_queens_mother"] = ids[drone]
    n_rpm = r["n_drone_queens"]
    if isinstance(n_rpm, int) and not isinstance(n_rpm, bool) and n_rpm >= 1:
        mating["n_drone_queens"] = n_rpm
    if r["single_drone"]:
        mating["n_drones"] = 1
    variety = r["race_description"] if isinstance(r["race_description"], str) else None
    q = QueenIn(id=ids[code], born=born, mother=ids.get(mother) if mother else None, mating=MatingIn(**mating),
                variety=variety, available_as_mother=False, available_as_drone_source=False)
    return q, cut


def _cut_cycles(ancestors: list[QueenIn]) -> int:
    """Coupe, tant qu'il en reste, un lien d'un ancêtre pris dans un cycle : la mère d'abord, puis la mère des
    RPM, sinon tout l'accouplement. Les ancêtres ne citent que des ancêtres : un cycle ne passe jamais par une
    reine de l'exploitation. Chaque tour retire un lien : la boucle se termine."""
    pos = {q.id: k for k, q in enumerate(ancestors)}
    cut = 0
    while True:
        cycle = find_cycle(ancestors)
        if not cycle:
            return cut
        members = set(cycle)
        cuts = ([(qid, {"mother": None}) for qid in cycle if ancestors[pos[qid]].mother in members]
                + [(qid, {"mating": ancestors[pos[qid]].mating.model_copy(
                    update={"drone_queens_mother": None, "external_ref": None})})
                   for qid in cycle if ancestors[pos[qid]].mating is not None
                   and ancestors[pos[qid]].mating.drone_queens_mother in members]
                + [(qid, {"mating": None}) for qid in cycle if ancestors[pos[qid]].mating is not None])
        if not cuts:                    # impossible : un membre d'un cycle a toujours un lien dans ce cycle
            return cut
        qid, update = cuts[0]
        ancestors[pos[qid]] = ancestors[pos[qid]].model_copy(update=update)
        cut += 1


def expand_external(inp: BlupV2Input, db_path: Optional[str] = None) -> Expansion:
    """Résout les références externes de l'entrée. db_path : chemin de la base (défaut : MELLIFERA_DB)."""
    refs = [(q, _mother_ref(q), _drone_ref(q)) for q in inp.queens]
    n_refs = sum((m is not None) + (d is not None) for _, m, d in refs)
    depth = inp.options.external_pedigree_depth
    out = Expansion(queens=list(inp.queens))
    if n_refs == 0 or depth == 0:
        return out
    path = db_path if db_path is not None else os.environ.get(ENV_DB, "")
    if not path:        # résolution non activée : aucun avertissement (les SDK antérieurs à 2.1.0 l'ignorent)
        return out
    try:
        conn = _connect(path)
        try:
            return _expand(inp, refs, depth, _Base(conn))
        finally:
            conn.close()
    except Exception as exc:     # sqlite3.Error attendu ; toute autre erreur est un bogue, mais jamais un 500
        if isinstance(exc, sqlite3.Error):
            logger.warning("base Index Mellifera illisible (%s) : %s", path, exc)
        else:
            logger.exception("résolution Index Mellifera en échec (%s)", path)
        out.warnings.append(WarningOut(code="external_pedigree_unavailable", message=(
            f"base Index Mellifera illisible : {n_refs} référence(s) externe(s) non résolue(s)")))
        return out


def _expand(inp: BlupV2Input, refs: list, depth: int, base: _Base) -> Expansion:
    found: dict[str, Optional[str]] = {}            # référence saisie → code en base
    for _, m, d in refs:
        for ref in (m, d):
            if ref is not None and ref not in found:
                r = base.find(ref)
                found[ref] = r["code"] if r is not None else None
    roots = sorted({c for c in found.values() if c is not None})
    rows, no_year, capped = _ancestry(base, roots, depth)

    used = {i for q in inp.queens for i in (q.id, q.mother) if i is not None}
    used |= {i for q in inp.queens if q.mating is not None
             for i in (*q.mating.drone_queens, q.mating.drone_queens_mother) if i is not None}
    start = min(0, min(used, default=0))            # id négatifs sans collision avec ceux de l'entrée
    ids = {code: start - k for k, code in enumerate(rows, start=1)}

    ancestors, cut = [], 0
    for r in rows.values():
        q, c = _ancestor(r, ids)
        ancestors.append(q)
        cut += c
    cut += _cut_cycles(ancestors)

    queens, linked = [], 0
    for q, m, d in refs:
        update = {}
        if m is not None and found[m] in ids:
            update["mother"] = ids[found[m]]
            linked += 1
        if d is not None and found[d] in ids:
            mtype = "insemination" if q.mating.type == "unknown" else q.mating.type
            update["mating"] = q.mating.model_copy(update={       # code en base : même nœud sire que les ancêtres
                "drone_queens_mother": ids[found[d]], "type": mtype, "external_ref": found[d]})
            linked += 1
        queens.append(q.model_copy(update=update) if update else q)

    out = Expansion(queens=queens, ancestors=ancestors)
    if any(c is not None for c in found.values()):
        msg = f"{linked} référence(s) Index Mellifera rattachée(s), {len(ancestors)} ancêtre(s) ajouté(s) au pedigree"
        if cut:
            msg += f", {cut} lien(s) erroné(s) de la base coupé(s)"
        if no_year:
            msg += f", {no_year} reine(s) sans année de naissance valide écartée(s)"
        if capped:
            msg += f", remontée arrêtée au plafond de {MAX_ANCESTORS} ancêtres"
        out.warnings.append(WarningOut(code="external_ref_resolved", message=msg))
    unknown = sorted(ref for ref, code in found.items() if code is None)
    if unknown:
        out.warnings.append(WarningOut(code="external_ref_unknown", message=(
            f"{len(unknown)} code(s) Index Mellifera introuvable(s), ignoré(s) jusqu'au prochain crawl : "
            f"{listed_ids(unknown)}")))
    return out
