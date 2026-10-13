import datetime as dt

import pytest

from mellifera.parse.common import ParseError
from mellifera.parse.queen import parse_queen
from mellifera.parse.queens_list import parse_queens_list
from tests.mellifera.fake_site import Q, fixture, queen_page


# ------------------------------------------------------------------ /breederbyyear/
def test_queens_list():
    page = parse_queens_list(fixture("breederbyyear_145_2022.html"))
    assert (page.breeder_id, page.trigram, page.country_id) == (145, "CBP", 8)
    assert page.years == tuple(range(2024, 2012, -1))
    assert len(page.rows) == 16
    first = page.rows[0]
    assert (first.code, first.birth_year, first.mother_code, first.drone_mother_code, first.cotation_count) == (
        "FR-CBP-A-0242-2022", 2022, "FR-CBP-A-0203-2020", "FR-SPD-B-0109-2018", 1)


def test_queens_list_code_comes_from_the_link_not_the_displayed_name():
    page = parse_queens_list(fixture("breederbyyear_145_2022.html"))
    row = next(r for r in page.rows if r.name == "FR-CBP-B-0293-2022")
    assert row.code == "FR-CBP-B-0193-2022"


# ------------------------------------------------------------------ /queen/
def test_queen_with_cotation_and_radar():
    q = parse_queen(fixture("queen_FR-CBP-B-0276-2022.html"))
    assert (q.code, q.name, q.birth_year, q.breeder_id) == ("FR-CBP-B-0276-2022", "FR-CBP-B-0276-2022", 2022, 145)
    assert (q.mother_code, q.drone_mother_code, q.n_drone_queens) == ("FR-CBP-B-0006-2020", "FR-CBP-B-0943-2019", 9)
    assert q.coi == pytest.approx(0.0078125) and q.coi_reliability == pytest.approx(0.754)
    assert q.not_found_raw == "31/126"
    assert q.insemination_date == dt.date(2022, 6, 15)
    assert (q.mating_code, q.mating_description, q.mating_category_raw) == (
        "ins", "Instrumental insemination", "Instrumental insemination (ins)")
    assert (q.race_code, q.race_description, q.lineage) == ("B", "buckfast", None)
    assert q.cotation_count == 1
    assert [(s.year, s.values) for s in q.radar] == [(2023, {
        "swarming": 6.0, "mildness": 6.0, "frameHolding": 6.0, "wintering": 5.0, "broodVitality": 6.0})]


def test_queen_without_cotation():
    q = parse_queen(fixture("queen_FR-NG-B-4536-2024.html"))
    assert (q.code, q.name, q.breeder_id, q.single_drone) == ("FR-NG-B-4536-2024", "B4536(NG)", 984, False)
    # données du site conservées telles quelles : lien mère 0098, commentaire 0298 (spec §2.2)
    assert (q.mother_code, q.mother_comment) == ("FR-NG-B-0098-2022", "FR-NG-B-0298-2022")
    assert (q.drone_mother_code, q.n_drone_queens) == ("FR-CBP-B-0104-2021", 3)
    assert (q.mating_code, q.mating_category_raw) == ("inng", "Instrumental insemination (FR)")
    assert q.lineage == "PrimGreen00"
    assert q.cotation_count == 0 and q.radar == ()


def test_queen_unknown_mother_single_drone_and_multi_year_radar():
    q = parse_queen(queen_page(Q("FR-PJ-B-0562-2020", 2020, 7, name="B562(PJ)1dr", drones="FR-PJ-B-0001-2018",
                                 number=1, radar={2021: [5, 4, 5, 5, 6], 2022: [6, 5, 5, 4, 6]})))
    assert q.mother_code is None and q.single_drone is True and q.n_drone_queens == 1
    assert [s.year for s in q.radar] == [2021, 2022]
    assert q.radar[1].values["wintering"] == 4.0


def test_queen_page_without_card_is_a_parse_error():
    with pytest.raises(ParseError):
        parse_queen(b"<html><body><h1>Page not found</h1></body></html>")
