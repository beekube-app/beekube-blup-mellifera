"""Commandes : python -m mellifera crawl | reparse | status | export-blup | schedule."""
from __future__ import annotations

import argparse
import contextlib
import fcntl
import json
import logging
import time
from pathlib import Path
from typing import Iterator, Optional

from .cache import RawCache
from .client import MissingContact, PoliteClient
from .config import Settings
from .crawl import Crawler, RunReport
from .db import connect_readonly, open_db
from .urls import KINDS

EXIT_OK, EXIT_NETWORK, EXIT_PARSE, EXIT_LOCKED, EXIT_CONFIG = 0, 1, 2, 3, 4
log = logging.getLogger("mellifera")


class LockBusy(RuntimeError):
    pass


@contextlib.contextmanager
def crawl_lock(path: Path) -> Iterator[None]:
    """Verrou exclusif non bloquant : deux crawls (ou crawl et reparse) ne tournent jamais ensemble."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a") as f:
        try:
            fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise LockBusy(f"un crawl est déjà en cours ({path})") from exc
        try:
            yield
        finally:
            fcntl.flock(f, fcntl.LOCK_UN)


def _exit_code(report: RunReport) -> int:
    for kind, rate in sorted(report.failing_kinds.items()):
        log.error("échecs de parsing « %s » : %.1f %% (seuil 5 %%) : structure du site modifiée ?", kind, 100 * rate)
    log.info("passage %s n° %s : %s", report.mode, report.run_id, json.dumps(report.stats, sort_keys=True))
    if report.stopped:
        return EXIT_NETWORK
    return EXIT_PARSE if report.failing_kinds else EXIT_OK


def run_crawl(settings: Settings, mode: str, *, transport=None, sleep=None, clock=None) -> int:
    """Crawl complet ou incrémental sous verrou. transport/sleep/clock : injection pour les tests."""
    try:
        with crawl_lock(settings.lock_path):
            with PoliteClient(settings.base_url, settings.contact, settings.delay, settings.max_errors,
                              transport=transport, sleep=sleep or time.sleep, clock=clock or time.monotonic) as client:
                conn = open_db(settings.db_path)
                try:
                    crawler = Crawler(conn, client, RawCache(settings.cache_dir))
                    report = crawler.full() if mode == "full" else crawler.incremental()
                finally:
                    conn.close()
    except LockBusy as exc:
        log.error("%s", exc)
        return EXIT_LOCKED
    return _exit_code(report)


def run_reparse(settings: Settings, kind: Optional[str]) -> int:
    try:
        with crawl_lock(settings.lock_path):
            conn = open_db(settings.db_path)
            try:
                report = Crawler(conn, None, RawCache(settings.cache_dir)).reparse(kind)
            finally:
                conn.close()
    except LockBusy as exc:
        log.error("%s", exc)
        return EXIT_LOCKED
    return _exit_code(report)


def status_lines(settings: Settings) -> list[str]:
    conn = connect_readonly(settings.db_path)
    try:
        out = []
        run = conn.execute("SELECT * FROM crawl_run ORDER BY id DESC LIMIT 1").fetchone()
        if run is None:
            out.append("aucun passage")
        else:
            state = f"terminé le {run['finished_at']}" if run["finished_at"] else "inachevé (reprise au prochain passage)"
            out.append(f"dernier passage : n° {run['id']} {run['mode']}, commencé le {run['started_at']}, {state}")
        unfinished = [f"n° {r['id']} {r['mode']}" for r in conn.execute(
            "SELECT id, mode FROM crawl_run WHERE finished_at IS NULL ORDER BY id")]
        if unfinished:
            out.append("passages inachevés : " + ", ".join(unfinished))
        for table in ("country", "breeder", "queen", "evaluation", "mating_code", "criterion", "page"):
            out.append(f"{table} : {conn.execute(f'SELECT COUNT(*) FROM {table}').fetchone()[0]}")
        for r in conn.execute("SELECT kind, status, COUNT(*) AS n, SUM(parse_error IS NOT NULL) AS errors "
                              "FROM page GROUP BY kind, status ORDER BY kind, status"):
            out.append(f"pages {r['kind']} {r['status']} : {r['n']} (échecs de parsing : {r['errors']})")
        unmapped = conn.execute("SELECT COUNT(*) FROM criterion WHERE trait_id IS NULL").fetchone()[0]
        out.append(f"critères non mappés : {unmapped}")
        return out
    finally:
        conn.close()


def run_export(settings: Settings, args: argparse.Namespace) -> int:
    from .export import export_blup
    conn = connect_readonly(settings.db_path)
    try:
        inp, report = export_blup(conn, args.breeder, from_year=args.from_year,
                                  sources=[s.strip() for s in args.sources.split(",") if s.strip()],
                                  profile=args.profile)
    finally:
        conn.close()
    Path(args.output).write_text(inp.model_dump_json(indent=1, exclude_none=True), encoding="utf-8")
    for line in report.lines():
        print(line)
    return EXIT_OK


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="python -m mellifera", description="Crawler d'Index Mellifera")
    sub = p.add_subparsers(dest="command", required=True)
    crawl = sub.add_parser("crawl", help="lit le site (GET poli, liste blanche)")
    mode = crawl.add_mutually_exclusive_group(required=True)
    mode.add_argument("--full", action="store_const", const="full", dest="mode")
    mode.add_argument("--incremental", action="store_const", const="incremental", dest="mode")
    reparse = sub.add_parser("reparse", help="reconstruit la base depuis le cache, sans réseau")
    reparse.add_argument("--kind", choices=KINDS)
    sub.add_parser("status", help="dernier passage et volumes")
    exp = sub.add_parser("export-blup", help="exporte une sélection au format BlupV2Input")
    exp.add_argument("--breeder", action="append", required=True,
                     help="trigramme ou id Index Mellifera de l'éleveur (répétable)")
    exp.add_argument("--from", dest="from_year", type=int, help="année de naissance minimale des reines choisies")
    exp.add_argument("--sources", default="season_table", help="season_table, radar (séparés par des virgules)")
    exp.add_argument("--profile", choices=("standard", "breeder"), default="standard")
    exp.add_argument("-o", "--output", required=True)
    sub.add_parser("schedule", help="crawl incrémental hebdomadaire (MELLIFERA_SCHEDULE)")
    return p


def main(argv: Optional[list[str]] = None, settings: Optional[Settings] = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    args = build_parser().parse_args(argv)
    try:
        settings = settings or Settings.from_env()
        if args.command == "crawl":
            return run_crawl(settings, args.mode)
        if args.command == "reparse":
            return run_reparse(settings, args.kind)
        if args.command == "status":
            print("\n".join(status_lines(settings)))
            return EXIT_OK
        if args.command == "export-blup":
            return run_export(settings, args)
        from .schedule import run_schedule
        if not settings.contact:
            raise MissingContact("MELLIFERA_CONTACT est vide : renseignez un contact avant de planifier le crawl")
        run_schedule(settings.schedule, settings.timezone, lambda: run_crawl(settings, "incremental"))
        return EXIT_OK
    except (MissingContact, FileNotFoundError, ValueError) as exc:
        log.error("%s", exc)
        return EXIT_CONFIG
