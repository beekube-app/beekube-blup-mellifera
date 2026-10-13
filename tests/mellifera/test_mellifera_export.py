import subprocess
import sys
from pathlib import Path

import pytest

from blup2.engine import run_blup
from blup2.schemas import BlupV2Input
from mellifera.db import open_db
from mellifera.export import export_blup
from mellifera.ingest import ingest
from mellifera.urls import queen_path
from tests.mellifera.fake_site import E, Q, queen_page, queens_list_page, season_page

ROOT = Path(__file__).resolve().parents[2]
SEEN = "2026-10-13T10:00:00.000000+00:00"
OTHER = ("xx", "Some place", "Other")         # code d'accouplement sans type reconnu
H = Q("FR-NG-B-0009-2016", 2016, 200, name="H009")                                     # autre éleveur (ascendance)
G = Q("FR-CBP-B-0050-2018", 2018, 145, mother="FR-ZZZ-B-0001-2010", drones=H.code, number=6)
M = Q("FR-CBP-B-0100-2020", 2020, 145, mother=G.code)
A = Q("FR-CBP-B-0001-2022", 2022, 145, mother=M.code, drones="FR-OLD-B-0001-2019", number=8,
      radar={2023: [6, 5, 6, 5, 6]})
B = Q("FR-CBP-B-0002-2022", 2022, 145, name="B002(CBP)1dr", mother=M.code, drones=G.code, mating=OTHER)
X = Q("FR-CBP-B-0010-2019", 2019, 145, mother="FR-CBP-B-0011-2019")                   # cycle X ↔ Y
Y = Q("FR-CBP-B-0011-2019", 2019, 145, mother="FR-CBP-B-0010-2019")
Z = Q("FR-CBP-B-0012-2019", 2019, 145, mother="FR-CBP-B-0012-2019")                   # sa propre mère
P = Q("FR-CBP-B-0013-2018", 2018, 145, mother=A.code)                                 # mère plus jeune
W = Q("FR-CBP-B-0014-2019", 2019, 145)                    # sa mère n'a pas d'année valide (voir la fixture)
D = Q("FR-NG-B-0001-2024", 2024, 200, mother=A.code)                                  # hors sélection
QUEENS = (H, G, M, A, B, X, Y, Z, P, W, D)


@pytest.fixture
def conn(tmp_path):
    c = open_db(tmp_path / "m.sqlite")
    pages = {queen_path(q.code): queen_page(q) for q in QUEENS}
    pages["/breederbyyear/145/2022/"] = queens_list_page(145, "CBP", 8, [2022], [A, B])
    pages["/breederbyyear/200/2024/"] = queens_list_page(200, "NG", 8, [2024], [D])
    pages["/breeder_Year/145/2023/"] = season_page([
        E(A.code, 2022, {"mildness": "5", "hyg24h": "90", "smr": "3"}),
        E(A.code, 2022, {"mildness": "4", "honeyAmt": "48"}, mode="ceta"),
        E(B.code, 2022, {"mildness": "9", "swarming": "-", "framHold": "bof"}),
        E(M.code, 2020, {"mildness": "6"}, date="-"),
    ])
    for path, body in pages.items():
        ingest(c, path, body, SEEN)
    # reine sans année valide (ligne de liste illisible) citée comme mère
    c.execute("INSERT INTO queen (code, code_valid, birth_year, breeder_id) VALUES ('FR-CBP-B-0099-2017', 1, NULL, 145)")
    c.execute("UPDATE queen SET mother_code = 'FR-CBP-B-0099-2017' WHERE code = ?", (W.code,))
    c.commit()
    yield c
    c.close()


def by_code(conn, inp: BlupV2Input) -> dict:
    ids = dict(conn.execute("SELECT id, code FROM queen").fetchall())
    return {ids[q.id]: q for q in inp.queens}


def test_export_is_valid_and_runs(conn):
    inp, report = export_blup(conn, ["CBP"])
    assert isinstance(inp, BlupV2Input) and inp.exploitation == 0
    out = run_blup(inp)
    assert {q.id for q in out.queens} == {q.id for q in inp.queens}
    assert all(q.id > 0 for q in inp.queens)


def test_selection_and_ancestry(conn):
    inp, report = export_blup(conn, ["CBP"])
    queens = by_code(conn, inp)
    assert D.code not in queens                                   # autre éleveur, pas un ancêtre
    assert H.code in queens                                       # ancêtre d'un autre éleveur
    assert (queens[H.code].available_as_mother, queens[H.code].available_as_drone_source) == (False, False)
    assert queens[A.code].available_as_mother is True
    assert report.n_selected == 9 and report.n_queens == 10


def test_pedigree_links(conn):
    inp, report = export_blup(conn, ["CBP"])
    queens = by_code(conn, inp)
    a, b, g = queens[A.code], queens[B.code], queens[G.code]
    assert a.mother == queens[M.code].id
    assert (a.mating.type, a.mating.external_ref, a.mating.n_drone_queens, a.mating.year) == (
        "insemination", "FR-OLD-B-0001-2019", 8, 2022)                # côté mâle absent de la base
    assert (g.mother, g.mother_ref_ext) == (None, "FR-ZZZ-B-0001-2010")
    assert g.mating.drone_queens_mother == queens[H.code].id
    assert (b.mating.type, b.mating.drone_queens_mother, b.mating.n_drones) == ("insemination", g.id, 1)
    assert report.counts["mating_type_inferred"] >= 1
    assert queens[W.code].mother_ref_ext == "FR-CBP-B-0099-2017"  # mère sans année valide
    assert report.counts["invalid_birth_year"] == 1


def test_errors_in_source_data_are_cut(conn):
    inp, report = export_blup(conn, ["CBP"])
    queens = by_code(conn, inp)
    assert queens[Z.code].mother is None and report.counts["self_reference_cut"] == 1
    assert queens[P.code].mother is None and report.counts["younger_parent_cut"] == 1
    assert report.counts["cycle_link_cut"] == 1
    assert (queens[X.code].mother is None) != (queens[Y.code].mother is None)   # un seul des deux liens coupé


def test_records_scales_and_exclusions(conn):
    inp, report = export_blup(conn, ["CBP"])
    names = {c.id: c.name for c in inp.criteria}
    assert set(names.values()) == {"mildness", "hyg24h"}
    kinds = {c.name: c.kind for c in inp.criteria}
    assert kinds == {"mildness": "note", "hyg24h": "percentage"}
    a_id = by_code(conn, inp)[A.code].id
    recs = {(names[r.criterion], r.queen): r for r in inp.records}
    mild = recs[("mildness", a_id)]
    assert (mild.value, mild.scale_min, mild.scale_max) == (4.5, 0.0, 6.0)   # evjp 5 et ceta 4 : moyenne
    assert (mild.apiary, mild.evaluator, mild.hive_class, str(mild.date)) == (145, 145, "hive", "2023-07-06")
    assert recs[("hyg24h", a_id)].value == 0.9
    c = report.counts
    assert (c["merged_duplicates"], c["harvest"], c["unmapped_criterion"], c["out_of_scale"], c["non_numeric"],
            c["no_date"]) == (1, 1, 1, 1, 1, 1)


def test_radar_source(conn):
    inp, report = export_blup(conn, ["CBP"], sources=["season_table", "radar"])
    names = {c.id: c.name for c in inp.criteria}
    a_id = by_code(conn, inp)[A.code].id
    radar = [r for r in inp.records if r.queen == a_id and str(r.date) == "2023-01-01"]
    assert {names[r.criterion] for r in radar} == {"swarming", "framHold", "winter", "broVital"}
    assert all(r.apiary == 145 and r.evaluator is None for r in radar)   # éleveur de la reine
    assert report.counts["radar_covered"] == 1                           # douceur déjà dans le tableau 2023


def test_from_year_filter_and_unknown_breeder(conn):
    inp, _ = export_blup(conn, ["CBP"], from_year=2022)
    queens = by_code(conn, inp)
    assert queens[A.code].available_as_mother and not queens[M.code].available_as_mother
    assert X.code not in queens
    with pytest.raises(ValueError, match="éleveur inconnu"):
        export_blup(conn, ["XYZ"])
    with pytest.raises(ValueError, match="sources inconnues"):
        export_blup(conn, ["CBP"], sources=["tables"])
    with pytest.raises(ValueError, match="aucune reine"):
        export_blup(conn, ["CBP"], from_year=2030)


def test_breeder_by_trigram_any_case_or_by_id(conn):
    reference, report = export_blup(conn, ["CBP"])
    assert report.lines()[0] == "éleveurs : CBP (145, pays inconnu)"      # pays absent de cette mini-base
    for ref in (["cbp"], ["145"], [" CBP ", "145"]):
        assert export_blup(conn, ref)[0] == reference


def test_trigram_shared_by_two_breeders_is_refused(conn):
    conn.execute("INSERT INTO country (id, name, iso) VALUES (3, 'Germany', 'DE')")
    conn.execute("INSERT INTO breeder (id, trigram, country_id, first_seen, last_seen) VALUES (999, 'CBP', 3, ?, ?)",
                 (SEEN, SEEN))
    with pytest.raises(ValueError, match=r"plusieurs éleveurs : 145 \(pays inconnu\), 999 \(Germany\)"):
        export_blup(conn, ["CBP"])
    inp, report = export_blup(conn, ["145"])
    assert report.n_selected == 9


def test_export_does_not_import_numpy():
    """L'image du service mellifera n'installe ni numpy ni scipy : l'export n'importe que blup2.schemas et
    blup2.pedigree."""
    code = ("import sys, mellifera.export; "
            "heavy = sorted(m for m in sys.modules if m.split('.')[0] in ('numpy', 'scipy')); "
            "assert not heavy, heavy")
    subprocess.run([sys.executable, "-c", code], cwd=ROOT, check=True)


class _CountingConn:
    """Connexion qui note les lignes lues de la table evaluation."""

    def __init__(self, conn):
        self._conn, self.evaluation_rows = conn, []

    def __getattr__(self, name):
        return getattr(self._conn, name)

    def execute(self, sql, *args):
        rows = self._conn.execute(sql, *args).fetchall()
        if "FROM evaluation" in sql:
            self.evaluation_rows += rows
        return iter_cursor(rows)


class iter_cursor(list):
    def fetchall(self):
        return list(self)

    def fetchone(self):
        return self[0] if self else None


def test_export_reads_only_the_evaluations_of_exported_queens(conn):
    """Sur la base complète (des centaines de milliers d'évaluations), l'export ne lit que celles des reines
    exportées."""
    for i in range(50):
        conn.execute("INSERT INTO evaluation (queen_name_raw, queen_code, queen_key, tester, eval_date, eval_year, "
                     "mode, criterion_code, value, value_raw, source, url) VALUES (?, ?, ?, '', ?, 2024, 'breeder', "
                     "'mildness', 3, '3', 'season_table', '/breeder_Year/200/2024/')",
                     (D.code, D.code, D.code, f"2024-01-{i % 28 + 1:02d}T{i:02d}"))
    conn.commit()
    counting = _CountingConn(conn)
    inp, _ = export_blup(counting, ["CBP"])
    assert inp.records
    assert counting.evaluation_rows
    assert D.code not in {r["queen_code"] for r in counting.evaluation_rows}
