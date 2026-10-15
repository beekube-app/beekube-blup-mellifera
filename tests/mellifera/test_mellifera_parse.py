import datetime as dt

import pytest

from mellifera.parse.common import is_code, labelled_fields, make_soup, parse_django_date
from mellifera.parse.countries import parse_countries
from mellifera.parse.country import parse_country
from tests.mellifera.fake_site import countries_page, country_page, fixture


# ------------------------------------------------------------------ outils communs
@pytest.mark.parametrize("text, expected", [
    ("Sept. 7, 2020", dt.date(2020, 9, 7)), ("June 6, 2024", dt.date(2024, 6, 6)),
    ("Nov. 5, 2021", dt.date(2021, 11, 5)), ("May 3, 2021", dt.date(2021, 5, 3)),
    ("-", None), ("", None), ("Feb. 30, 2021", None),
])
def test_parse_django_date(text, expected):
    assert parse_django_date(text) == expected


@pytest.mark.parametrize("code, ok", [
    ("FR-CBP-B-0276-2022", True), ("LU-PJ-V-005H-2015", True), ("B4536(NG)", False), ("A087", False),
    # formes réelles du site : éleveur avec minuscules ou vide, race à plusieurs lettres, numéro avec minuscules
    ("DE-ABg-B-0003-2017", True), ("DE-AS-Msw-0017-2013", True), ("DE-AH-B-082b-2007", True),
    ("DE-AS-B-0mix-2013", True), ("DE---B-0072-1993", True), ("DE-ABg-B-0003-17", False), ("FR-CBP-B-2022", False),
])
def test_is_code(code, ok):
    assert is_code(code) is ok


def test_labelled_fields_ignore_unclosed_paragraphs_and_stop_at_headings():
    soup = make_soup(b"<html><body><div><p><strong>A : </strong>1<p><strong> -- B : </strong> deux "
                     b"<a href='/x/'>lien</a><h1>Details</h1>parasite</div></body></html>")
    fields = labelled_fields(soup.find("div"))
    assert [(f.label, f.text, f.hrefs) for f in fields] == [("A", "1", ()), ("B", "deux lien", ("/x/",))]


# ------------------------------------------------------------------ /all/
def test_countries():
    countries = parse_countries(fixture("all.html"))
    assert len(countries) == 32
    assert sum(c.n_breeders for c in countries) == 754
    by_id = {c.id: c for c in countries}
    assert (by_id[8].name, by_id[8].iso, by_id[8].n_breeders) == ("France", "FR", 85)
    assert by_id[35].iso == "RE"                    # drapeau « RE_HNhsPa6.png »
    assert (by_id[0].name, by_id[0].iso) == ("Unknown", None)


# ------------------------------------------------------------------ /countrybyyear/
def test_country():
    breeders = parse_country(fixture("countrybyyear_8.html"))
    assert len(breeders) == 85
    first = breeders[0]
    assert (first.id, first.trigram, first.season_years) == (1140, "AAD", (2024, 2022, 2021))
    cbp = next(b for b in breeders if b.trigram == "CBP")
    assert cbp.id == 145 and 2022 in cbp.season_years


def test_synthetic_country_pages_parse_like_the_site():
    countries = parse_countries(countries_page([(8, "France", "FR", 2), (0, "Unknown", None, 0)]))
    assert [(c.id, c.name, c.n_breeders) for c in countries] == [(8, "France", 2), (0, "Unknown", 0)]
    breeders = parse_country(country_page([(145, "CBP", [2022, 2021]), (200, "NG", [])]))
    assert [(b.id, b.trigram, b.season_years) for b in breeders] == [(145, "CBP", (2022, 2021)), (200, "NG", ())]
    assert parse_country(country_page([])) == ()           # pays sans éleveur (« Unknown ») : vide, pas une erreur
