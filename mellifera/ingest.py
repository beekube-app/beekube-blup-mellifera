"""Écriture des pages parsées dans la base (idempotente) et identification des reines évaluées (spec §6)."""
from __future__ import annotations

import logging
import re
import sqlite3
from typing import Optional

from .db import now_iso
from .parse.common import PARSER_VERSION, is_code, parse_number
from .parse.countries import CountryEntry, parse_countries
from .parse.country import BreederEntry, parse_country
from .parse.queen import QueenPage, parse_queen
from .parse.queens_list import QueensListPage, parse_queens_list
from .parse.season import SeasonPage, parse_season
from .urls import kind_of, path_ids

log = logging.getLogger("mellifera")

_TRIGRAM = re.compile(r"\(([A-Z0-9]{2,6})\)")


def mating_type_for(category_raw: Optional[str], description: Optional[str]) -> str:
    """Type d'accouplement d'un code : catégorie de la légende d'abord, description en repli."""
    for text in (category_raw, description):
        t = (text or "").lower()
        if "insemination" in t:
            return "insemination"
        if "mating station" in t or "belegstelle" in t:
            return "station"
    return "unknown"


def _upsert_mating_code(conn: sqlite3.Connection, code: Optional[str], description: Optional[str],
                        category_raw: Optional[str]) -> None:
    if not code:
        return
    conn.execute(
        "INSERT INTO mating_code (code, description, category_raw, mating_type) VALUES (?, ?, ?, ?) "
        "ON CONFLICT (code) DO UPDATE SET description = COALESCE(excluded.description, description), "
        "category_raw = COALESCE(excluded.category_raw, category_raw), mating_type = excluded.mating_type",
        (code, description, category_raw, mating_type_for(category_raw, description)))


def _ensure_criterion(conn: sqlite3.Connection, mode: str, code: str) -> None:
    """Colonne ou mode inconnu : critère créé sans trait (ignoré par l'export tant qu'il n'est pas mappé)."""
    conn.execute("INSERT OR IGNORE INTO criterion (mode, code) VALUES (?, ?)", (mode, code))


def store_countries(conn: sqlite3.Connection, entries: tuple[CountryEntry, ...]) -> None:
    for c in entries:
        conn.execute("INSERT INTO country (id, name, iso) VALUES (?, ?, ?) "
                     "ON CONFLICT (id) DO UPDATE SET name = excluded.name, iso = excluded.iso", (c.id, c.name, c.iso))


def _upsert_breeder(conn: sqlite3.Connection, breeder_id: int, trigram: str, country_id: Optional[int],
                    seen_at: str) -> None:
    conn.execute(
        "INSERT INTO breeder (id, trigram, country_id, first_seen, last_seen) VALUES (?, ?, ?, ?, ?) "
        "ON CONFLICT (id) DO UPDATE SET trigram = excluded.trigram, "
        "country_id = COALESCE(excluded.country_id, country_id), last_seen = MAX(last_seen, excluded.last_seen)",
        (breeder_id, trigram, country_id, seen_at, seen_at))


def store_country(conn: sqlite3.Connection, country_id: int, breeders: tuple[BreederEntry, ...], seen_at: str) -> None:
    for b in breeders:
        _upsert_breeder(conn, b.id, b.trigram, country_id, seen_at)


def store_queens_list(conn: sqlite3.Connection, page: QueensListPage, seen_at: str) -> None:
    """Éleveur + une ébauche par reine inconnue (la fiche, lue ensuite, fait foi et complète la ligne)."""
    _upsert_breeder(conn, page.breeder_id, page.trigram, page.country_id, seen_at)
    for r in page.rows:
        conn.execute(
            "INSERT OR IGNORE INTO queen (code, code_valid, name, birth_year, breeder_id, mother_code, "
            "drone_mother_code, serial_raw) VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            (r.code, int(is_code(r.code)), r.name, r.birth_year, page.breeder_id, r.mother_code,
             r.drone_mother_code, r.serial_raw))


def store_queen(conn: sqlite3.Connection, q: QueenPage, path: str, seen_at: str) -> None:
    _upsert_mating_code(conn, q.mating_code, q.mating_description, q.mating_category_raw)
    code = path_ids(path)[0]    # le code de l'URL fait foi : l'affiché peut différer (casse, saisie)
    if code != q.code:
        log.warning("%s affiche le code %r : fiche enregistrée sous %r", path, q.code, code)
    values = dict(
        code=code, code_valid=int(is_code(code)), name=q.name, birth_year=q.birth_year, breeder_id=q.breeder_id,
        mother_code=q.mother_code, drone_mother_code=q.drone_mother_code, n_drone_queens=q.n_drone_queens,
        mating_code=q.mating_code, single_drone=int(q.single_drone),
        insemination_date=q.insemination_date.isoformat() if q.insemination_date else None,
        race_code=q.race_code, race_description=q.race_description, lineage=q.lineage,
        generation_raw=q.generation_raw, coi=q.coi, coi_reliability=q.coi_reliability,
        not_found_raw=q.not_found_raw, cotation_count=q.cotation_count, serial_raw=q.serial_raw,
        origin_raw=q.origin_raw, data_origin_raw=q.data_origin_raw, mother_comment_raw=q.mother_comment,
        drone_comment_raw=q.drone_comment, parse_version=PARSER_VERSION, fetched_at=seen_at, parsed_at=now_iso(),
    )
    cols = ", ".join(values)
    marks = ", ".join("?" for _ in values)
    updates = ", ".join(f"{c} = excluded.{c}" for c in values if c != "code")
    conn.execute(f"INSERT INTO queen ({cols}) VALUES ({marks}) ON CONFLICT (code) DO UPDATE SET {updates}",
                 tuple(values.values()))
    for series in q.radar:
        for label, value in series.values.items():
            _ensure_criterion(conn, "radar", label)
            _upsert_evaluation(conn, dict(
                queen_name_raw=q.name, queen_birth_year_raw=q.birth_year, queen_code=code, queen_key=code,
                page_breeder_id=q.breeder_id, tester="", tester_breeder_id=None,
                eval_date=f"{series.year}-01-01", eval_year=series.year, mode="radar", criterion_code=label,
                value=value, value_raw=f"{value:g}", source="radar", url=path))


def _upsert_evaluation(conn: sqlite3.Connection, row: dict) -> None:
    cols = ", ".join(row)
    marks = ", ".join("?" for _ in row)
    conn.execute(f"INSERT INTO evaluation ({cols}) VALUES ({marks}) "
                 "ON CONFLICT (queen_key, source, mode, tester, eval_date, criterion_code) DO UPDATE SET "
                 "value = excluded.value, value_raw = excluded.value_raw, url = excluded.url, "
                 "tester_breeder_id = excluded.tester_breeder_id, eval_year = excluded.eval_year",
                 tuple(row.values()))


def lookup_queen(conn: sqlite3.Connection, name: str, birth_year: Optional[int],
                 page_breeder_id: Optional[int]) -> Optional[str]:
    """Étapes 3 et 4 de la spec §6 : (nom, année, éleveur), puis (nom, année) si unique en base."""
    if birth_year is None:
        return None
    m = _TRIGRAM.search(name)
    if m:
        rows = conn.execute("SELECT q.code FROM queen q JOIN breeder b ON b.id = q.breeder_id "
                            "WHERE q.name = ? AND q.birth_year = ? AND b.trigram = ?",
                            (name, birth_year, m.group(1))).fetchall()
    else:
        rows = conn.execute("SELECT code FROM queen WHERE name = ? AND birth_year = ? AND breeder_id = ?",
                            (name, birth_year, page_breeder_id)).fetchall()
    if len(rows) == 1:
        return rows[0]["code"]
    rows = conn.execute("SELECT code FROM queen WHERE name = ? AND birth_year = ?", (name, birth_year)).fetchall()
    return rows[0]["code"] if len(rows) == 1 else None


def unresolved_key(name: str, birth_year: Optional[int], page_breeder_id: Optional[int]) -> str:
    """Clé d'une ligne non résolue ; l'éleveur de la page évite de fusionner deux « A087 » d'élevages différents."""
    return f"name:{name}|{'' if birth_year is None else birth_year}|{'' if page_breeder_id is None else page_breeder_id}"


def store_season(conn: sqlite3.Connection, page: SeasonPage, path: str, page_breeder_id: int) -> None:
    for mc in page.mating_codes:
        _upsert_mating_code(conn, mc.code, mc.description, mc.category_raw)
    for e in page.evaluations:
        code = e.queen_code or lookup_queen(conn, e.queen_name_raw, e.queen_birth_year, page_breeder_id)
        for criterion, raw in e.values.items():
            _ensure_criterion(conn, e.mode, criterion)
            _upsert_evaluation(conn, dict(
                queen_name_raw=e.queen_name_raw, queen_birth_year_raw=e.queen_birth_year, queen_code=code,
                queen_key=code or unresolved_key(e.queen_name_raw, e.queen_birth_year, page_breeder_id),
                page_breeder_id=page_breeder_id, tester=e.tester, tester_breeder_id=e.tester_breeder_id,
                eval_date=e.eval_date.isoformat() if e.eval_date else "", eval_year=e.eval_year, mode=e.mode,
                criterion_code=criterion, value=parse_number(raw), value_raw=raw, source="season_table", url=path))


def resolve_pending(conn: sqlite3.Connection) -> int:
    """Relance les étapes 3 et 4 sur les évaluations non résolues ; une ligne résolue prend son code pour clé
    et fusionne avec la ligne identique déjà enregistrée sous ce code. Retourne le nombre de reines résolues."""
    pending = conn.execute("SELECT DISTINCT queen_key, queen_name_raw, queen_birth_year_raw, page_breeder_id "
                           "FROM evaluation WHERE queen_code IS NULL").fetchall()
    resolved = 0
    for p in pending:
        code = lookup_queen(conn, p["queen_name_raw"], p["queen_birth_year_raw"], p["page_breeder_id"])
        if code is None:
            continue
        conn.execute("UPDATE OR IGNORE evaluation SET queen_code = ?, queen_key = ? WHERE queen_key = ?",
                     (code, code, p["queen_key"]))
        # Ce qui reste sous l'ancienne clé doublonne une ligne déjà enregistrée sous le code.
        conn.execute("DELETE FROM evaluation WHERE queen_key = ?", (p["queen_key"],))
        resolved += 1
    return resolved


def ingest(conn: sqlite3.Connection, path: str, body: bytes, seen_at: str):
    """Parse la page selon son type et l'enregistre. Lève ParseError si la page ne se parse pas."""
    kind = kind_of(path)
    ids = path_ids(path)
    if kind == "countries":
        parsed = parse_countries(body)
        store_countries(conn, parsed)
    elif kind == "country":
        parsed = parse_country(body)
        store_country(conn, int(ids[0]), parsed, seen_at)
    elif kind == "queens_list":
        parsed = parse_queens_list(body)
        store_queens_list(conn, parsed, seen_at)
    elif kind == "queen":
        parsed = parse_queen(body)
        store_queen(conn, parsed, path, seen_at)
    elif kind == "season":
        parsed = parse_season(body)
        store_season(conn, parsed, path, int(ids[0]))
    else:
        raise ValueError(f"chemin hors liste blanche : {path}")
    return parsed
