import re
import sqlite3

import pytest

from mellifera.db import connect_readonly, now_iso, open_db


@pytest.fixture
def conn(tmp_path):
    c = open_db(tmp_path / "m.sqlite")
    yield c
    c.close()


def count(conn, sql: str, *args) -> int:
    return conn.execute(f"SELECT COUNT(*) FROM {sql}", args).fetchone()[0]


# ------------------------------------------------------------------ schéma et référentiel
def test_schema_wal_and_reference_tables(conn):
    assert conn.execute("PRAGMA journal_mode").fetchone()[0] == "wal"
    assert count(conn, "trait") == 22
    assert count(conn, "criterion WHERE trait_id IS NULL") == 0
    mapped = dict(conn.execute("SELECT c.mode || '.' || c.code, t.name FROM criterion c JOIN trait t "
                               "ON t.id = c.trait_id").fetchall())
    assert mapped["radar.frameHolding"] == mapped["evjp.framHold"] == mapped["ceta.framHold"] == "framHold"


def test_reopening_keeps_unknown_criteria_and_reloads_mapping(tmp_path):
    c = open_db(tmp_path / "m.sqlite")
    c.execute("INSERT INTO criterion (mode, code) VALUES ('evjp', 'smr')")
    c.commit()
    c.close()
    c = open_db(tmp_path / "m.sqlite")
    assert count(c, "criterion WHERE mode = 'evjp' AND code = 'smr' AND trait_id IS NULL") == 1
    assert count(c, "criterion WHERE trait_id IS NOT NULL") == 48
    c.close()


def test_readonly_connection_for_status_and_export(tmp_path):
    with pytest.raises(FileNotFoundError, match="base introuvable"):
        connect_readonly(tmp_path / "absente.sqlite")
    open_db(tmp_path / "m.sqlite").close()
    ro = connect_readonly(tmp_path / "m.sqlite")
    assert ro.execute("SELECT COUNT(*) FROM trait").fetchone()[0] == 22
    with pytest.raises(sqlite3.OperationalError):
        ro.execute("DELETE FROM trait")
    ro.close()


def test_timestamps_are_utc_to_the_microsecond():
    # « éleveur vu pour la première fois pendant ce passage » se décide en comparant à crawl_run.started_at
    assert re.fullmatch(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d\.\d{6}\+00:00", now_iso())
