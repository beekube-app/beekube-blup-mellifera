"""Parcours du site : crawl complet, incrémental, fermeture sur les parents, reprise, relecture du cache."""
from __future__ import annotations

import datetime as dt
import json
import logging
import sqlite3
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Callable, Optional

from .cache import RawCache
from .client import PoliteClient, TooManyErrors
from .db import now_iso
from .ingest import ingest, resolve_pending
from .parse.common import PARSER_VERSION, ParseError
from .parse.country import BreederEntry
from .urls import KINDS, countries_path, country_path, kind_of, queen_path, queens_list_path, season_path

log = logging.getLogger("mellifera")

MAX_PARSE_FAILURE = 0.05        # au-delà de 5 % d'échecs de parsing sur un type de page : sortie en erreur
RECENT_YEARS = 2                # incrémental : années ≥ Y − 2


@dataclass
class RunReport:
    run_id: int
    mode: str
    stats: dict
    stopped: bool = False                                        # arrêt sur erreurs réseau consécutives
    failing_kinds: dict[str, float] = field(default_factory=dict)  # type → taux d'échec de parsing (> 5 %)


class Crawler:
    def __init__(self, conn: sqlite3.Connection, client: Optional[PoliteClient], cache: RawCache, *,
                 today: Callable[[], dt.date] = dt.date.today):
        self.conn, self.client, self.cache = conn, client, cache
        self.today = today
        self.run_id = 0
        self.started_at = ""
        self.pages: dict[str, Counter] = defaultdict(Counter)    # type → statut → nombre
        self.parse_errors: Counter = Counter()
        self.parsed_ok: Counter = Counter()
        self.counters: Counter = Counter()

    # ------------------------------------------------------------------ passage
    def _start(self, mode: str) -> None:
        row = self.conn.execute("SELECT id, started_at FROM crawl_run WHERE mode = ? AND finished_at IS NULL "
                                "ORDER BY id DESC LIMIT 1", (mode,)).fetchone()
        if row is not None and mode != "reparse":
            self.run_id, self.started_at = row["id"], row["started_at"]
            log.info("reprise du passage %s n° %s", mode, self.run_id)
        else:
            self.started_at = now_iso()
            cur = self.conn.execute("INSERT INTO crawl_run (mode, started_at) VALUES (?, ?)", (mode, self.started_at))
            self.run_id = cur.lastrowid
        self.conn.commit()

    def _stats(self) -> dict:
        return {
            "pages": {k: dict(v) for k, v in sorted(self.pages.items())},
            "parse_errors": dict(self.parse_errors),
            "queens_added": self.counters["queens_added"],
            "queens_updated": self.counters["queens_updated"],
            "resumed_pages": self.counters["resumed_pages"],
            "invalid_codes": self.counters["invalid_codes"],          # codes impossibles à mettre dans une URL
            "nonconforming_codes": self.conn.execute(                 # codes hors format international (gardés)
                "SELECT COUNT(*) FROM queen WHERE code_valid = 0").fetchone()[0],
            "evaluations_resolved": self.counters["evaluations_resolved"],
            "evaluations_unresolved": self.conn.execute(
                "SELECT COUNT(DISTINCT queen_key) FROM evaluation WHERE queen_code IS NULL").fetchone()[0],
        }

    def _finish(self, mode: str, stopped: bool) -> RunReport:
        if not stopped:
            self.counters["evaluations_resolved"] += resolve_pending(self.conn)
        stats = self._stats()
        failing = {k: round(self.parse_errors[k] / (self.parse_errors[k] + self.parsed_ok[k]), 3)
                   for k in self.parse_errors
                   if self.parse_errors[k] / (self.parse_errors[k] + self.parsed_ok[k]) > MAX_PARSE_FAILURE}
        self.conn.execute("UPDATE crawl_run SET finished_at = ?, stats = ? WHERE id = ?",
                          (None if stopped else now_iso(), json.dumps(stats, sort_keys=True), self.run_id))
        self.conn.commit()
        return RunReport(self.run_id, mode, stats, stopped, failing)

    # ------------------------------------------------------------------ pages
    def _fetch(self, path: str) -> tuple[Optional[bytes], str]:
        """Corps de la page (None si absente ou refusée) et date de lecture. Une page déjà lue dans ce passage
        (reprise) est relue depuis le cache, sans requête."""
        kind = kind_of(path)
        row = self.conn.execute("SELECT status, cache_path, fetched_at, run_id FROM page WHERE url = ?",
                                (path,)).fetchone()
        if (row is not None and row["run_id"] == self.run_id and row["status"] != "error"
                and (row["status"] == "ok" or kind != "countries")):    # /all/ illisible : redemandée
            self.counters["resumed_pages"] += 1
            self.pages[kind][row["status"]] += 1
            body = self.cache.get(row["cache_path"]) if row["status"] == "ok" else None
            return body, row["fetched_at"]
        resp = self.client.get(path)
        fetched_at = now_iso()
        rel = sha = None
        if resp.status == "ok":
            rel, sha = self.cache.put(path, resp.body)
        self.conn.execute(
            "INSERT INTO page (url, kind, status, http_status, fetched_at, sha256, cache_path, run_id) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?) ON CONFLICT (url) DO UPDATE SET status = excluded.status, "
            "http_status = excluded.http_status, fetched_at = excluded.fetched_at, "
            "sha256 = COALESCE(excluded.sha256, sha256), cache_path = COALESCE(excluded.cache_path, cache_path), "
            "parse_error = NULL, run_id = excluded.run_id",
            (path, kind, resp.status, resp.http_status, fetched_at, sha, rel, self.run_id))
        self.pages[kind][resp.status] += 1
        if resp.status == "error":
            log.warning("%s illisible (%s) : page ignorée pour ce passage", path,
                        f"HTTP {resp.http_status}" if resp.http_status else "délai de lecture dépassé")
        return (resp.body if resp.status == "ok" else None), fetched_at

    def _parse(self, path: str, body: bytes, seen_at: str):
        """Parse et enregistre la page. Un échec (ParseError, ou autre exception d'un parseur devant une structure
        inattendue) annule les écritures de la page, est noté dans page.parse_error et n'arrête pas le crawl :
        sinon la reprise relirait la même page en cache et s'arrêterait au même endroit, indéfiniment."""
        kind = kind_of(path)
        self.conn.execute("SAVEPOINT ingest")
        try:
            parsed = ingest(self.conn, path, body, seen_at)
        except sqlite3.Error:
            raise                                   # base indisponible : ce n'est pas un défaut de la page
        except Exception as exc:
            self.conn.execute("ROLLBACK TO ingest")
            error = str(exc) if isinstance(exc, ParseError) else f"{type(exc).__name__}: {exc}"
            self.conn.execute("UPDATE page SET parse_error = ?, parse_version = ? WHERE url = ?",
                              (error, PARSER_VERSION, path))
            self.parse_errors[kind] += 1
            log.warning("échec de parsing %s : %s", path, error, exc_info=not isinstance(exc, ParseError))
            parsed = None
        else:
            self.conn.execute("UPDATE page SET parse_error = NULL, parse_version = ? WHERE url = ?",
                              (PARSER_VERSION, path))
            self.parsed_ok[kind] += 1
        self.conn.execute("RELEASE ingest")
        self.conn.commit()
        return parsed

    def _visit(self, path: str):
        body, seen_at = self._fetch(path)
        if body is None:
            self.conn.commit()
            return None
        return self._parse(path, body, seen_at)

    def _visit_queen(self, code: str, list_count: Optional[int] = None) -> None:
        path = queen_path(code)
        if path is None:
            self.counters["invalid_codes"] += 1
            return
        before = self.conn.execute("SELECT parsed_at FROM queen WHERE code = ?", (code,)).fetchone()
        if self._visit(path) is None:
            return
        self.counters["queens_updated" if before is not None and before["parsed_at"] else "queens_added"] += 1
        if list_count is not None:
            self.conn.execute("UPDATE queen SET list_cotation_count = ? WHERE code = ?", (list_count, code))
            self.conn.commit()

    # ------------------------------------------------------------------ parcours
    def _breeder(self, b: BreederEntry, min_year: Optional[int]) -> None:
        """Listes de reines (toutes les années, ou années ≥ min_year), fiches à lire, saisons."""
        if min_year is None:        # parcours complet : dernière saison connue (page qui existe)
            first = b.season_years[0] if b.season_years else self.today().year
        else:                       # incrémental : année courante ; sa navigation donne les autres années
            first = self.today().year
        first_page = self._visit(queens_list_path(b.id, first))
        years = set(first_page.years) if first_page is not None else set(b.season_years)
        lists = [first_page] if first_page is not None else []
        for y in sorted(years - {first}, reverse=True):
            if min_year is None or y >= min_year:
                page = self._visit(queens_list_path(b.id, y))
                if page is not None:
                    lists.append(page)
        for page in lists:
            for row in page.rows:
                known = self.conn.execute("SELECT parsed_at, list_cotation_count FROM queen WHERE code = ?",
                                          (row.code,)).fetchone()
                stale = (known is None or known["parsed_at"] is None
                         or known["list_cotation_count"] != row.cotation_count)
                if min_year is None or stale:
                    self._visit_queen(row.code, row.cotation_count)
        for y in b.season_years:
            if min_year is None or y >= min_year:
                self._visit(season_path(b.id, y))

    def _closure(self) -> None:
        """Lit toute reine citée comme mère ou « Drones » et absente de la base, jusqu'à épuisement."""
        while True:
            refs = [r[0] for r in self.conn.execute(
                "SELECT ref FROM (SELECT mother_code AS ref FROM queen UNION SELECT drone_mother_code FROM queen) "
                "WHERE ref IS NOT NULL AND ref NOT IN (SELECT code FROM queen WHERE parsed_at IS NOT NULL) "
                "ORDER BY ref")]
            todo = []
            for code in refs:
                path = queen_path(code)
                if path is None:
                    continue
                seen = self.conn.execute("SELECT run_id FROM page WHERE url = ?", (path,)).fetchone()
                if seen is None or seen["run_id"] != self.run_id:
                    todo.append(code)
            if not todo:
                return
            for code in todo:
                self._visit_queen(code)

    def _countries(self, new_breeder: Callable[[BreederEntry], bool], min_year: Optional[int]) -> None:
        countries = self._visit(countries_path())
        if countries is None:       # liste des pays refusée, absente ou illisible : rien ne peut être parcouru
            raise TooManyErrors(f"{countries_path()} illisible : aucun pays à parcourir")
        for c in countries:
            for b in self._visit(country_path(c.id)) or ():
                self._breeder(b, None if new_breeder(b) else min_year)

    def _run(self, mode: str, body: Callable[[], None]) -> RunReport:
        self._start(mode)
        try:
            body()
        except TooManyErrors as exc:
            log.error("arrêt : %s ; le prochain passage %s reprendra ici", exc, mode)
            self.conn.commit()
            return self._finish(mode, stopped=True)
        return self._finish(mode, stopped=False)

    def full(self) -> RunReport:
        def body():
            self._countries(lambda b: True, None)
            self._closure()
        return self._run("full", body)

    def incremental(self) -> RunReport:
        """Sans crawl complet abouti, ou avec un complet interrompu, reprend le complet : sinon l'incrémental
        (années récentes seulement) masquerait les pages jamais lues."""
        pending = self.conn.execute("SELECT 1 FROM crawl_run WHERE mode = 'full' AND finished_at IS NULL").fetchone()
        done = self.conn.execute("SELECT 1 FROM crawl_run WHERE mode = 'full' AND finished_at IS NOT NULL").fetchone()
        if pending is not None or done is None:
            log.warning("aucun crawl complet abouti : le passage incrémental est remplacé par le crawl complet")
            return self.full()

        def is_new(b: BreederEntry) -> bool:
            row = self.conn.execute("SELECT first_seen FROM breeder WHERE id = ?", (b.id,)).fetchone()
            return row is None or row["first_seen"] >= self.started_at

        def body():
            self._countries(is_new, self.today().year - RECENT_YEARS)
            self._closure()
        return self._run("incremental", body)

    def reparse(self, kind: Optional[str] = None) -> RunReport:
        """Reconstruit la base depuis le cache (aucune requête). Les évaluations des types relus sont
        effacées puis recréées ; les autres tables sont mises à jour en place (id des reines conservés)."""
        def body():
            kinds = [k for k in KINDS if kind is None or k == kind]
            if "season" in kinds:
                self.conn.execute("DELETE FROM evaluation WHERE source = 'season_table'")
            if "queen" in kinds:
                self.conn.execute("DELETE FROM evaluation WHERE source = 'radar'")
            for k in kinds:
                for row in self.conn.execute("SELECT url, cache_path, fetched_at FROM page "
                                             "WHERE kind = ? AND status = 'ok' ORDER BY url", (k,)).fetchall():
                    self.pages[k]["ok"] += 1
                    self._parse(row["url"], self.cache.get(row["cache_path"]), row["fetched_at"])
        return self._run("reparse", body)
