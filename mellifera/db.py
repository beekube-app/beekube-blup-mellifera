"""Base SQLite (mode WAL) : schéma, référentiel des critères, connexions."""
from __future__ import annotations

import csv
import datetime as dt
import sqlite3
from pathlib import Path

HERE = Path(__file__).parent
TRAITS_CSV = HERE / "traits.csv"
CRITERIA_CSV = HERE / "criteria.csv"

SCHEMA = """
CREATE TABLE IF NOT EXISTS country (
    id INTEGER PRIMARY KEY,
    name TEXT NOT NULL,
    iso TEXT
);
CREATE TABLE IF NOT EXISTS breeder (
    id INTEGER PRIMARY KEY,
    trigram TEXT NOT NULL,
    country_id INTEGER,
    first_seen TEXT NOT NULL,
    last_seen TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS breeder_trigram ON breeder (trigram);
CREATE TABLE IF NOT EXISTS queen (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    code TEXT NOT NULL UNIQUE,
    code_valid INTEGER NOT NULL,
    name TEXT,
    birth_year INTEGER,
    breeder_id INTEGER,
    mother_code TEXT,
    drone_mother_code TEXT,
    n_drone_queens INTEGER,
    mating_code TEXT,
    single_drone INTEGER NOT NULL DEFAULT 0,
    insemination_date TEXT,
    race_code TEXT,
    race_description TEXT,
    lineage TEXT,
    generation_raw TEXT,
    coi REAL,
    coi_reliability REAL,
    not_found_raw TEXT,
    cotation_count INTEGER,
    list_cotation_count INTEGER,
    serial_raw TEXT,
    origin_raw TEXT,
    data_origin_raw TEXT,
    mother_comment_raw TEXT,
    drone_comment_raw TEXT,
    parse_version INTEGER,
    fetched_at TEXT,
    parsed_at TEXT
);
CREATE INDEX IF NOT EXISTS queen_name ON queen (name, birth_year);
CREATE TABLE IF NOT EXISTS mating_code (
    code TEXT PRIMARY KEY,
    description TEXT,
    category_raw TEXT,
    mating_type TEXT NOT NULL CHECK (mating_type IN ('natural', 'station', 'insemination', 'unknown')),
    mating_type_override TEXT CHECK (mating_type_override IN ('natural', 'station', 'insemination', 'unknown'))
);
CREATE VIEW IF NOT EXISTS queen_mating AS
    SELECT q.*, COALESCE(m.mating_type_override, m.mating_type, 'unknown') AS mating_type
    FROM queen q LEFT JOIN mating_code m ON m.code = q.mating_code;
CREATE TABLE IF NOT EXISTS trait (
    id INTEGER PRIMARY KEY,
    name TEXT NOT NULL UNIQUE,
    label TEXT NOT NULL,
    kind TEXT NOT NULL CHECK (kind IN ('note', 'percentage', 'harvest')),
    scale_min REAL,
    scale_max REAL,
    direction INTEGER NOT NULL CHECK (direction IN (1, -1)),
    CHECK (kind != 'note' OR scale_max > scale_min)
);
CREATE TABLE IF NOT EXISTS criterion (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    mode TEXT NOT NULL,
    code TEXT NOT NULL,
    trait_id INTEGER REFERENCES trait (id),
    UNIQUE (mode, code)
);
CREATE TABLE IF NOT EXISTS evaluation (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    queen_name_raw TEXT NOT NULL,
    queen_birth_year_raw INTEGER,
    queen_code TEXT,
    queen_key TEXT NOT NULL,
    page_breeder_id INTEGER,
    tester TEXT NOT NULL,
    tester_breeder_id INTEGER,
    eval_date TEXT NOT NULL,
    eval_year INTEGER,
    mode TEXT NOT NULL,
    criterion_code TEXT NOT NULL,
    value REAL,
    value_raw TEXT NOT NULL,
    source TEXT NOT NULL CHECK (source IN ('season_table', 'radar')),
    url TEXT NOT NULL,
    UNIQUE (queen_key, source, mode, tester, eval_date, criterion_code)
);
CREATE INDEX IF NOT EXISTS evaluation_queen ON evaluation (queen_code);
CREATE TABLE IF NOT EXISTS page (
    url TEXT PRIMARY KEY,
    kind TEXT NOT NULL,
    status TEXT NOT NULL CHECK (status IN ('ok', 'missing', 'forbidden', 'error')),
    http_status INTEGER,
    fetched_at TEXT NOT NULL,
    sha256 TEXT,
    cache_path TEXT,
    parse_error TEXT,
    parse_version INTEGER,
    run_id INTEGER
);
CREATE TABLE IF NOT EXISTS crawl_run (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    mode TEXT NOT NULL CHECK (mode IN ('full', 'incremental', 'reparse')),
    started_at TEXT NOT NULL,
    finished_at TEXT,
    stats TEXT
);
"""


def now_iso() -> str:
    """Horodatage UTC à la microseconde : « vu pour la première fois pendant ce passage » se décide par
    comparaison avec crawl_run.started_at, sans égalité possible entre deux passages."""
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="microseconds")


def connect(path: Path) -> sqlite3.Connection:
    conn = sqlite3.connect(path)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA busy_timeout=5000")
    return conn


def connect_readonly(path: Path) -> sqlite3.Connection:
    """Lecture seule (export, moteur) : n'écrit rien, même pendant un crawl (WAL)."""
    if not Path(path).is_file():
        raise FileNotFoundError(f"base introuvable : {path} (lancez d'abord « python -m mellifera crawl --full »)")
    conn = sqlite3.connect(f"file:{Path(path).resolve()}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    return conn


def _opt_float(s: str):
    return float(s) if s.strip() else None


def load_reference(conn: sqlite3.Connection) -> None:
    """Recharge traits.csv et criteria.csv (sources de vérité) ; les critères inconnus restent sans trait."""
    with TRAITS_CSV.open(encoding="utf-8") as f:
        for r in csv.DictReader(f):
            conn.execute(
                "INSERT INTO trait (id, name, label, kind, scale_min, scale_max, direction) VALUES (?, ?, ?, ?, ?, ?, ?) "
                "ON CONFLICT (id) DO UPDATE SET name = excluded.name, label = excluded.label, kind = excluded.kind, "
                "scale_min = excluded.scale_min, scale_max = excluded.scale_max, direction = excluded.direction",
                (int(r["id"]), r["name"], r["label"], r["kind"], _opt_float(r["scale_min"]),
                 _opt_float(r["scale_max"]), int(r["direction"])))
    conn.execute("UPDATE criterion SET trait_id = NULL")
    with CRITERIA_CSV.open(encoding="utf-8") as f:
        for r in csv.DictReader(f):
            trait = conn.execute("SELECT id FROM trait WHERE name = ?", (r["trait"],)).fetchone()
            if trait is None:
                raise ValueError(f"criteria.csv : trait inconnu « {r['trait']} » ({r['mode']}, {r['code']})")
            conn.execute("INSERT INTO criterion (mode, code, trait_id) VALUES (?, ?, ?) "
                         "ON CONFLICT (mode, code) DO UPDATE SET trait_id = excluded.trait_id",
                         (r["mode"], r["code"], trait["id"]))
    conn.commit()


def open_db(path: Path) -> sqlite3.Connection:
    """Base en écriture (crawler) : crée le dossier, le schéma et recharge le référentiel."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    conn = connect(path)
    conn.executescript(SCHEMA)
    load_reference(conn)
    return conn
