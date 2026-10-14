"""Contrat entre la base du crawler (mellifera/db.py) et le moteur (blup2/external.py)."""
from blup2.external import expand_external
from blup2.schemas import BlupV2Input, CriterionIn, MatingIn, QueenIn
from mellifera.db import open_db
from mellifera.ingest import ingest
from mellifera.urls import queen_path
from tests.mellifera.fake_site import Q, queen_page

SEEN = "2026-10-14T10:00:00.000000+00:00"
GRAND = Q("FR-CBP-B-0100-2020", 2020, 145)
DAUGHTER = Q("FR-CBP-B-0276-2022", 2022, 145, mother=GRAND.code, drones="FR-NG-B-0009-2018", number=8)


def crawler_db(path):
    conn = open_db(path)
    for q in (GRAND, DAUGHTER):
        ingest(conn, queen_path(q.code), queen_page(q), SEEN)
    conn.commit()
    return conn


def beekube_input():
    return BlupV2Input(exploitation=1, criteria=[CriterionIn(id=1)], queens=[
        QueenIn(id=1, born=2024, mother_ref_ext=DAUGHTER.code,
                mating=MatingIn(type="unknown", external_ref=GRAND.code))])


def test_engine_reads_the_crawler_database(tmp_path):
    crawler_db(tmp_path / "m.sqlite").close()
    exp = expand_external(beekube_input(), str(tmp_path / "m.sqlite"))
    assert [w.code for w in exp.warnings] == ["external_ref_resolved"]
    daughter = next(a for a in exp.ancestors if a.id == exp.queens[0].mother)
    assert (daughter.born, daughter.mother) == (2022, exp.queens[0].mating.drone_queens_mother)
    assert daughter.mating.type == "insemination" and daughter.mating.n_drone_queens == 8
    assert daughter.mating.external_ref == "FR-NG-B-0009-2018"


def test_engine_reads_while_the_crawler_writes(tmp_path):
    """Calcul pendant un crawl (WAL) : la lecture seule n'attend pas l'écrivain et voit le dernier état validé."""
    writer = crawler_db(tmp_path / "m.sqlite")
    writer.execute("BEGIN IMMEDIATE")
    writer.execute("UPDATE queen SET birth_year = 1800 WHERE code = ?", (DAUGHTER.code,))
    exp = expand_external(beekube_input(), str(tmp_path / "m.sqlite"))
    assert len(exp.ancestors) == 2
    writer.rollback()
    writer.close()
