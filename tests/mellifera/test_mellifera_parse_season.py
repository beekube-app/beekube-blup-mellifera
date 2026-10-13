import datetime as dt

import pytest

from mellifera.parse.common import ParseError
from mellifera.parse.countries import parse_countries
from mellifera.parse.country import parse_country
from mellifera.parse.queen import parse_queen
from mellifera.parse.queens_list import parse_queens_list
from mellifera.parse.season import parse_season
from tests.mellifera.fake_site import E, fixture, season_page

PARSERS = [parse_countries, parse_country, parse_queens_list, parse_queen, parse_season]


# ------------------------------------------------------------------ /breeder_Year/
def test_season_two_tables_with_different_columns():
    page = parse_season(fixture("breeder_Year_145_2022.html"))
    assert {e.mode for e in page.evaluations} == {"evjp", "ceta"}
    assert all(e.queen_code is not None for e in page.evaluations)
    ev = next(e for e in page.evaluations if e.mode == "evjp" and e.queen_code == "FR-CBP-A-0203-2020")
    assert (ev.tester, ev.tester_breeder_id, ev.eval_year, ev.eval_date) == ("CBP", 145, 2023, dt.date(2023, 7, 6))
    assert ev.values["mildness"] == "6" and "hyg24h" not in ev.values      # « - » = absent
    ceta = next(e for e in page.evaluations if e.mode == "ceta" and e.queen_code == "FR-CBP-A-0087-2020")
    assert (ceta.values["honeyAmt"], ceta.values["hyg24h"]) == ("48", "100")
    assert {m.code: m.category_raw for m in page.mating_codes}["mdm"] == "Mating station (FR)"


def test_season_names_resolved_by_page_links():
    page = parse_season(fixture("breeder_Year_145_2021.html"))
    ev = page.evaluations[0]
    assert (ev.queen_name_raw, ev.queen_code, ev.mode) == ("D992(CBP)", "FR-CBP-D-0992-2019", "evjp")
    assert {e.mode for e in page.evaluations} == {"evjp", "evdk"}


def test_season_link_wins_over_code_looking_name():
    body = season_page([E("FR-CBP-B-0293-2022", 2022, {"mildness": "5"})],
                       links=[("FR-CBP-B-0193-2022", "FR-CBP-B-0293-2022")])
    assert parse_season(body).evaluations[0].queen_code == "FR-CBP-B-0193-2022"


def test_season_ambiguous_or_unlinked_name_left_unresolved():
    body = season_page([E("A087", 2020, {"mildness": "5"}), E("FR-CBP-A-0041-2020", 2020, {"mildness": "4"})],
                       links=[("FR-CBP-A-0087-2020", "A087"), ("FR-XX-A-0087-2020", "A087")])
    evs = parse_season(body).evaluations
    assert evs[0].queen_code is None                    # deux liens différents pour « A087 »
    assert evs[1].queen_code == "FR-CBP-A-0041-2020"    # pas de lien : texte au format code


def test_season_row_with_missing_cells_is_a_parse_error():
    body = season_page([E("A087", 2020, {"mildness": "5"})]).replace(b'<td class="td01">5</td>', b"")
    with pytest.raises(ParseError):
        parse_season(body)


# ------------------------------------------------------------------ robustesse de tous les parseurs
@pytest.mark.parametrize("parse", PARSERS)
@pytest.mark.parametrize("body", [b"", b"   ", b"<html><body><div class='cadree'><h1>FR-"])
def test_empty_or_truncated_page_is_a_parse_error(parse, body):
    with pytest.raises(ParseError):
        parse(body)


@pytest.mark.parametrize("parse, name", [
    (parse_countries, "all.html"), (parse_country, "countrybyyear_8.html"),
    (parse_queens_list, "breederbyyear_145_2022.html"), (parse_queen, "queen_FR-CBP-B-0276-2022.html"),
    (parse_season, "breeder_Year_145_2022.html"),
])
def test_truncated_fixture_is_a_parse_error(parse, name):
    body = fixture(name)
    with pytest.raises(ParseError):
        parse(body[: len(body) // 2])
