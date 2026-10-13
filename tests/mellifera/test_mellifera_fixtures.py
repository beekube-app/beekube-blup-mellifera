"""Les fixtures commitées sont anonymisées : vérifications de forme (aucune liste de vrais noms ici)."""
import re

import pytest
from bs4 import BeautifulSoup

from tests.mellifera.fake_site import FIXTURES

ALL = sorted(p.name for p in FIXTURES.glob("*.html"))


def test_fixtures_present():
    assert ALL == ["all.html", "breeder_Year_145_2021.html", "breeder_Year_145_2022.html",
                   "breeder_Year_984_2026.html", "breederbyyear_145_2022.html", "countrybyyear_8.html",
                   "queen_FR-CBP-B-0276-2022.html", "queen_FR-NG-B-4536-2024.html"]


@pytest.mark.parametrize("name", ALL)
def test_rating_authors_are_anonymous(name):
    soup = BeautifulSoup((FIXTURES / name).read_bytes(), "html.parser")
    for a in soup.find_all("a", href=re.compile(r"^/rating/\d+$")):
        assert a.get_text(" ", strip=True).endswith(" - Testeur anonyme")


@pytest.mark.parametrize("name", [n for n in ALL if n.startswith("queen_")])
def test_import_comment_is_neutral(name):
    html = (FIXTURES / name).read_text(encoding="utf-8")
    m = re.search(r"<strong>Import : </strong>(?:(?!<strong>).)*<strong> -- Comment : </strong>\s*([^<]*?)\s*<",
                  html, re.S)
    assert m is not None and m.group(1) == "No"


@pytest.mark.parametrize("name", [n for n in ALL if n.startswith("breeder_Year_")])
def test_evaluation_comments_are_empty(name):
    soup = BeautifulSoup((FIXTURES / name).read_bytes(), "html.parser")
    tables = soup.find_all("table")
    assert tables
    for table in tables:
        rows = table.find_all("tr")
        header = [c.get_text(strip=True) for c in rows[0].find_all("td")]
        if header[-1:] != ["Comment"]:
            continue
        for tr in rows[1:]:
            assert tr.find_all("td")[-1].get_text(strip=True) == ""
