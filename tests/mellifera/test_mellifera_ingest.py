import pytest

from mellifera.db import open_db
from mellifera.ingest import ingest, mating_type_for, resolve_pending
from tests.mellifera.fake_site import E, Q, db_snapshot, fixture, queen_page, queens_list_page, season_page

SEEN = "2026-10-13T10:00:00+00:00"


@pytest.fixture
def conn(tmp_path):
    c = open_db(tmp_path / "m.sqlite")
    yield c
    c.close()


def load(conn, pages: dict[str, bytes]) -> None:
    for path, body in pages.items():
        ingest(conn, path, body, SEEN)
    conn.commit()


def count(conn, sql: str, *args) -> int:
    return conn.execute(f"SELECT COUNT(*) FROM {sql}", args).fetchone()[0]


FIXTURE_PAGES = {
    "/all/": fixture("all.html"),
    "/countrybyyear/8/": fixture("countrybyyear_8.html"),
    "/breederbyyear/145/2022/": fixture("breederbyyear_145_2022.html"),
    "/queen/FR-CBP-B-0276-2022/": fixture("queen_FR-CBP-B-0276-2022.html"),
    "/queen/FR-NG-B-4536-2024/": fixture("queen_FR-NG-B-4536-2024.html"),
    "/breeder_Year/145/2021/": fixture("breeder_Year_145_2021.html"),
    "/breeder_Year/145/2022/": fixture("breeder_Year_145_2022.html"),
    "/breeder_Year/984/2026/": fixture("breeder_Year_984_2026.html"),
}


# ------------------------------------------------------------------ type d'accouplement
@pytest.mark.parametrize("category, description, expected", [
    ("Instrumental insemination (FR)", "Thorigné (FR) instrumental insemination", "insemination"),
    ("Mating station (FR)", "Madame island (FR) mating station", "station"),
    (None, "Belegstelle Norderney", "station"),
    (None, "Instrumental insemination", "insemination"),
    ("Free mating", "Open mating", "unknown"),
    (None, None, "unknown"),
])
def test_mating_type_for(category, description, expected):
    assert mating_type_for(category, description) == expected


# ------------------------------------------------------------------ ingestion des fixtures
def test_ingest_fixtures(conn):
    load(conn, FIXTURE_PAGES)
    assert count(conn, "country") == 32
    assert count(conn, "breeder") == 85
    q = conn.execute("SELECT * FROM queen_mating WHERE code = 'FR-CBP-B-0276-2022'").fetchone()
    assert (q["mother_code"], q["drone_mother_code"], q["n_drone_queens"], q["mating_type"]) == (
        "FR-CBP-B-0006-2020", "FR-CBP-B-0943-2019", 9, "insemination")
    stub = conn.execute("SELECT * FROM queen WHERE code = 'FR-CBP-A-0242-2022'").fetchone()
    assert (stub["mother_code"], stub["parsed_at"]) == ("FR-CBP-A-0203-2020", None)   # ébauche issue de la liste
    radar = conn.execute("SELECT eval_date, value FROM evaluation WHERE source = 'radar' "
                         "AND criterion_code = 'wintering'").fetchall()
    assert [tuple(r) for r in radar] == [("2023-01-01", 5.0)]
    assert count(conn, "evaluation WHERE queen_code IS NULL") == 0
    assert conn.execute("SELECT mating_type FROM mating_code WHERE code = 'mdm'").fetchone()[0] == "station"


def test_ingest_is_idempotent(conn):
    load(conn, FIXTURE_PAGES)
    first = db_snapshot(conn)
    load(conn, FIXTURE_PAGES)
    assert db_snapshot(conn) == first


def test_same_evaluation_on_two_season_pages_is_stored_once(conn):
    rows = [E("FR-CBP-B-0104-2021", 2021, {"mildness": "5", "swarming": "6"})]
    load(conn, {"/breeder_Year/145/2022/": season_page(rows), "/breeder_Year/984/2026/": season_page(rows)})
    assert count(conn, "evaluation") == 2


def test_unknown_mode_and_column_create_unmapped_criteria(conn):
    load(conn, {"/breeder_Year/145/2022/": season_page([E("FR-CBP-B-0104-2021", 2021, {"smr": "40", "mildness": "5"},
                                                         mode="xyz")])})
    assert count(conn, "criterion WHERE mode = 'xyz' AND trait_id IS NULL") == 2
    assert count(conn, "evaluation WHERE mode = 'xyz'") == 2


def test_non_numeric_value_kept_raw(conn):
    load(conn, {"/breeder_Year/145/2022/": season_page([E("FR-CBP-B-0104-2021", 2021, {"mildness": "bof"})])})
    r = conn.execute("SELECT value, value_raw FROM evaluation").fetchone()
    assert (r["value"], r["value_raw"]) == (None, "bof")


# ------------------------------------------------------------------ identification des reines évaluées (spec §6)
def _list(*queens: Q) -> bytes:
    return queens_list_page(145, "CBP", 8, [2020], list(queens))


def test_row_by_code(conn):
    load(conn, {"/breeder_Year/145/2022/": season_page([E("FR-CBP-A-0041-2020", 2020, {"mildness": "4"})])})
    assert conn.execute("SELECT queen_code, queen_key FROM evaluation").fetchone()[:] == (
        "FR-CBP-A-0041-2020", "FR-CBP-A-0041-2020")


def test_row_by_name_with_trigram_resolved_in_db(conn):
    load(conn, {"/queen/FR-SPD-B-0131-2019/": queen_page(Q("FR-SPD-B-0131-2019", 2019, 12, name="B131(SPD)")),
                "/breederbyyear/12/2019/": queens_list_page(12, "SPD", 8, [2019], [])})
    load(conn, {"/breeder_Year/145/2022/": season_page([E("B131(SPD)", 2019, {"mildness": "4"})])})
    assert conn.execute("SELECT queen_code FROM evaluation").fetchone()[0] == "FR-SPD-B-0131-2019"


def test_row_by_bare_name_uses_page_breeder(conn):
    load(conn, {"/breederbyyear/145/2020/": _list(Q("FR-CBP-A-0087-2020", 2020, 145, name="A087")),
                "/breederbyyear/12/2020/": queens_list_page(12, "SPD", 8, [2020],
                                                            [Q("FR-SPD-A-0087-2020", 2020, 12, name="A087")])})
    load(conn, {"/breeder_Year/145/2022/": season_page([E("A087", 2020, {"mildness": "4"})])})
    assert conn.execute("SELECT queen_code FROM evaluation").fetchone()[0] == "FR-CBP-A-0087-2020"


def test_ambiguous_bare_name_left_unresolved(conn):
    load(conn, {"/breederbyyear/11/2020/": queens_list_page(11, "AAA", 8, [2020], [Q("FR-AAA-A-0087-2020", 2020, 11,
                                                                                       name="A087")]),
                "/breederbyyear/12/2020/": queens_list_page(12, "SPD", 8, [2020],
                                                            [Q("FR-SPD-A-0087-2020", 2020, 12, name="A087")])})
    load(conn, {"/breeder_Year/145/2022/": season_page([E("A087", 2020, {"mildness": "4"})])})
    r = conn.execute("SELECT queen_code, queen_key FROM evaluation").fetchone()
    assert (r["queen_code"], r["queen_key"]) == (None, "name:A087|2020|145")


def test_unresolved_row_resolved_after_queen_arrives_and_merged(conn):
    row = E("A087", 2020, {"mildness": "4"})
    load(conn, {"/breeder_Year/145/2022/": season_page([row])})
    assert count(conn, "evaluation WHERE queen_code IS NULL") == 1
    # la même évaluation, vue ailleurs sous son code (page d'un autre éleveur)
    load(conn, {"/breeder_Year/984/2026/": season_page([E("FR-CBP-A-0087-2020", 2020, {"mildness": "4"})])})
    load(conn, {"/breederbyyear/145/2020/": _list(Q("FR-CBP-A-0087-2020", 2020, 145, name="A087"))})
    assert resolve_pending(conn) == 1
    rows = conn.execute("SELECT queen_code, queen_key FROM evaluation").fetchall()
    assert [tuple(r) for r in rows] == [("FR-CBP-A-0087-2020", "FR-CBP-A-0087-2020")]
