"""Volet 2 Index Mellifera : résolution dans le moteur, isolation des ancêtres (spec IM §8.1, §8.3)."""
import sqlite3

import pytest

from blup2.engine import run_blup, run_matings
from blup2.schemas import MatingPairIn, MatingV2Input
from tests.builders import criterion, mating, queen, record
from tests.external_base import CYCLES, inp, make_db, warning


@pytest.fixture
def db(tmp_path):
    return make_db(tmp_path / "mellifera.sqlite")


# ------------------------------------------------------------------ base absente ou illisible
@pytest.mark.parametrize("make", ["unset", "missing", "garbage", "no_table", "old_schema"])
def test_unavailable_base_gives_the_same_computation(tmp_path, monkeypatch, make):
    """Base illisible : avertissement, calcul identique. Variable absente : résolution non activée, aucun
    avertissement (la 2.1.0 peut tourner avant la publication du SDK qui connaît les nouveaux codes)."""
    path = tmp_path / "m.sqlite"
    if make == "garbage":
        path.write_bytes(b"pas une base sqlite" * 100)
    elif make == "no_table":
        sqlite3.connect(path).execute("CREATE TABLE other (x)").connection.close()
    elif make == "old_schema":           # queen_mating sans une colonne lue par le moteur
        make_db(path)
        conn = sqlite3.connect(path)
        conn.executescript("DROP VIEW queen_mating; CREATE VIEW queen_mating AS SELECT id, code, name, birth_year, "
                           "mother_code, drone_mother_code, 'unknown' AS mating_type FROM queen;")
        conn.close()
    if make != "unset":
        monkeypatch.setenv("MELLIFERA_DB", str(path))
    qs, recs = family(refs=False)
    out = run_blup(inp(qs + [queen(9, 2024, mother_ref_ext="FR-CBP-B-0276-2022")], records=recs))
    expected = [] if make == "unset" else [
        "base Index Mellifera illisible : 1 référence(s) externe(s) non résolue(s)"]
    assert [w.message for w in out.warnings if w.code.startswith("external")] == expected
    monkeypatch.delenv("MELLIFERA_DB", raising=False)
    ref = run_blup(inp(qs + [queen(9, 2024)], records=recs))
    assert [q.model_dump() for q in out.queens] == [q.model_dump() for q in ref.queens]


def test_no_reference_no_warning_even_without_base(monkeypatch):
    monkeypatch.delenv("MELLIFERA_DB", raising=False)
    qs, recs = family(refs=False)
    out = run_blup(inp(qs, records=recs))
    assert not [w for w in out.warnings if w.code.startswith("external")]


# ------------------------------------------------------------------ isolation dans le moteur (§8.3)
def family(refs=True):
    """Exploitation de 5 reines ; avec refs, 4 et 5 citent des lignées Index Mellifera (3 références)."""
    ext = dict(m4=dict(mother_ref_ext="FR-CBP-B-0276-2022"),
               m5=dict(mother_ref_ext="FR-CBP-B-0277-2022", mating=mating("unknown", ext="FR-CBP-B-0100-2020")))
    qs = [queen(1, 2021, mating=mating("station", station=1, n_rpm=8)), queen(2, 2022, mother=1),
          queen(3, 2022, mother=1), queen(4, 2023, **(ext["m4"] if refs else {})),
          queen(5, 2023, **(ext["m5"] if refs else {}))]
    vals = {1: 0.4, 2: 0.6, 3: 0.7, 4: 0.8, 5: 0.5}
    return qs, [record(k, v, date=f"{2022 + (k > 3)}-05-01") for k, v in vals.items()]


def test_ancestors_never_reach_outputs_base_or_suggestions(db, monkeypatch):
    qs, recs = family()
    monkeypatch.setenv("MELLIFERA_DB", db)
    with_db = run_blup(inp(qs, records=recs))
    monkeypatch.delenv("MELLIFERA_DB")
    without = run_blup(inp(qs, records=recs))
    assert [q.id for q in with_db.queens] == [q.id for q in without.queens] == [1, 2, 3, 4, 5]
    assert with_db.base == without.base
    assert {(m.mother, m.drone_source) for m in with_db.matings} == {(m.mother, m.drone_source) for m in without.matings}
    assert all(i > 0 for m in with_db.matings for i in (m.mother, m.drone_source))
    assert all(p.birth_year in (2021, 2022, 2023) for pts in with_db.trend.values() for p in pts)
    assert [w.code for w in with_db.warnings if w.code == "unknown_parent_added"] == []
    assert warning(with_db, "external_ref_resolved") == [
        "3 référence(s) Index Mellifera rattachée(s), 6 ancêtre(s) ajouté(s) au pedigree"]


def test_queens_of_sister_lineages_become_related(db, monkeypatch):
    qs = [queen(1, 2024, mother_ref_ext="FR-CBP-B-0276-2022", mating=mating("station", station=1)),
          queen(2, 2024, mother_ref_ext="FR-CBP-B-0277-2022", mating=mating("station", station=2))]
    pair = MatingV2Input(exploitation=1, criteria=[criterion(1)], queens=qs,
                         pairs=[MatingPairIn(mother=1, drone_source=2)])
    monkeypatch.delenv("MELLIFERA_DB", raising=False)
    assert run_matings(pair).matings[0].expected_inbreeding == 0.0
    monkeypatch.setenv("MELLIFERA_DB", db)
    assert run_matings(pair).matings[0].expected_inbreeding > 0.0


def test_errors_in_the_base_never_fail_the_computation(tmp_path, monkeypatch):
    monkeypatch.setenv("MELLIFERA_DB", make_db(tmp_path / "m.sqlite", CYCLES))
    qs = [queen(1, 2024, mother_ref_ext="FR-C-B-0001-2020", mating=mating("unknown", ext="FR-C-B-0004-2017")),
          queen(2, 2024, mother=1)]
    out = run_blup(inp(qs))
    assert {q.id for q in out.queens} == {1, 2}
    assert warning(out, "external_ref_resolved")[0].endswith("3 lien(s) erroné(s) de la base coupé(s)")


# ------------------------------------------------------------------ relecture finale
def test_cycle_between_farm_and_ancestors_falls_back_to_a_computation_without_them(db, monkeypatch):
    """C1 (filet de sécurité) : si le pedigree avec ancêtres est cyclique, le calcul se fait sans eux."""
    import blup2.engine as engine
    real = engine.build_pedigree

    def flaky(queens, n, extra_queens=()):
        if list(extra_queens):
            raise ValueError("pedigree cyclique : reines concernées -1, 5")
        return real(queens, n, extra_queens)
    monkeypatch.setattr(engine, "build_pedigree", flaky)
    monkeypatch.setenv("MELLIFERA_DB", db)
    qs, recs = family()
    out = run_blup(inp(qs, records=recs))
    assert [q.id for q in out.queens] == [1, 2, 3, 4, 5]
    assert warning(out, "external_pedigree_unavailable") == [
        "ascendance Index Mellifera incompatible avec le pedigree de l'exploitation : calcul sans les ancêtres"]
    assert warning(out, "external_ref_resolved") == []
