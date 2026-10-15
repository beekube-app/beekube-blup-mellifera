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


# ------------------------------------------------------------------ volet 3 : notes et contemporaines
def test_engine_reads_the_notes_of_the_crawler_database(tmp_path):
    """Requêtes de blup2/external_records.py sur une base créée par mellifera/db.py et peuplée par l'ingestion."""
    from blup2.schemas import RecordIn
    from mellifera.urls import season_path
    from tests.mellifera.fake_site import E, season_page

    conn = crawler_db(tmp_path / "m.sqlite")
    radar = Q(DAUGHTER.code, 2022, 145, mother=GRAND.code, drones="FR-NG-B-0009-2018", number=8,
              radar={2024: [4, 3, 5, 2, 1]})                       # douceur 3 en 2024
    peers = [Q("FR-T-B-0001-2022", 2022, 300), Q("FR-CBP-B-0300-2023", 2023, 145, radar={2024: [1, 6, 1, 1, 1]})]
    for q in (radar, *peers):
        ingest(conn, queen_path(q.code), queen_page(q), SEEN)
    ingest(conn, season_path(300, 2023), season_page(
        [E(DAUGHTER.code, 2022, {"mildness": "5"}, tester_id=300, date="July 6, 2023"),
         E(peers[0].code, 2022, {"mildness": "3"}, tester_id=300, date="July 6, 2023")]), SEEN)
    conn.commit()
    conn.close()
    src = beekube_input().model_copy(update={
        "criteria": [CriterionIn(id=1, trait="gentleness")],
        "records": [RecordIn(queen=1, criterion=1, value=0.5, date="2025-06-01", apiary=1)]})
    exp = expand_external(BlupV2Input.model_validate(src.model_dump()), str(tmp_path / "m.sqlite"))
    assert [w.message for w in exp.warnings if w.code == "external_records_added"] == [
        "4 note(s) Index Mellifera ajoutée(s) sur 1 critère(s) : 1 ancêtre(s) et 2 contemporaine(s)"]
    assert [(a.id, a.born) for a in exp.ancestors] == [(-1, 2020), (-2, 2022)]      # GRAND, DAUGHTER
    assert [(q.id, q.born) for q in exp.contemporaries] == [(-3, 2023), (-4, 2022)]
    assert sorted((r.queen, round(r.value, 4), r.date.isoformat()) for r in exp.records) == [
        (-4, 0.5, "2023-07-06"), (-3, 1.0, "2024-01-01"), (-2, 0.5, "2024-01-01"), (-2, 0.8333, "2023-07-06")]


def test_notes_are_read_while_the_crawler_writes(tmp_path):
    """Lecture des notes pendant une transaction d'écriture du crawler (WAL) : dernier état validé."""
    from blup2.schemas import RecordIn
    from mellifera.urls import season_path
    from tests.mellifera.fake_site import E, season_page

    writer = crawler_db(tmp_path / "m.sqlite")
    ingest(writer, queen_path("FR-T-B-0001-2022"), queen_page(Q("FR-T-B-0001-2022", 2022, 300)), SEEN)
    ingest(writer, season_path(300, 2023), season_page(
        [E(DAUGHTER.code, 2022, {"mildness": "5"}, tester_id=300, date="July 6, 2023"),
         E("FR-T-B-0001-2022", 2022, {"mildness": "3"}, tester_id=300, date="July 6, 2023")]), SEEN)
    writer.commit()
    writer.execute("BEGIN IMMEDIATE")
    writer.execute("UPDATE evaluation SET value = 0")
    src = beekube_input().model_copy(update={
        "criteria": [CriterionIn(id=1, trait="gentleness")],
        "records": [RecordIn(queen=1, criterion=1, value=0.5, date="2025-06-01", apiary=1)]})
    exp = expand_external(BlupV2Input.model_validate(src.model_dump()), str(tmp_path / "m.sqlite"))
    assert sorted(round(r.value, 4) for r in exp.records) == [0.5, 0.8333]
    writer.rollback()
    writer.close()
