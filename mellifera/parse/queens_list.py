"""Page `/breederbyyear/<éleveur>/<année>/` (« Queens list ») : reines nées cette année chez l'éleveur."""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Optional

from bs4 import Tag

from .common import ParseError, make_soup, node_text, parse_int, queen_code_from_href

_BREEDER = re.compile(r"^/breeder/(\d+)/$")
_COUNTRY = re.compile(r"^/countrybyyear/(\d+)/$")
_YEAR = re.compile(r"^/breederbyyear/(\d+)/(\d{4})/$")


@dataclass(frozen=True)
class QueenListRow:
    code: str                           # code de l'URL /queen/<code>/ (le nom affiché peut différer)
    name: str
    birth_year: Optional[int]
    mother_code: Optional[str]
    drone_mother_code: Optional[str]    # colonne « Drones (RAM Mother) »
    cotation_count: Optional[int]
    serial_raw: Optional[str]


@dataclass(frozen=True)
class QueensListPage:
    breeder_id: int
    trigram: str
    country_id: Optional[int]
    years: tuple[int, ...]              # années proposées par la navigation de la page
    rows: tuple[QueenListRow, ...]


def _code(col: Tag) -> Optional[str]:
    a = col.find("a", href=True)
    return queen_code_from_href(a["href"]) if a is not None else None


def parse_queens_list(body: bytes) -> QueensListPage:
    soup = make_soup(body)
    breeder = soup.find("a", href=_BREEDER)
    if breeder is None:
        raise ParseError("éleveur absent du fil d'Ariane")
    country = soup.find("a", href=_COUNTRY)
    years = sorted({int(_YEAR.match(a["href"]).group(2)) for a in soup.find_all("a", href=_YEAR)}, reverse=True)
    rows: list[QueenListRow] = []
    for li in soup.find_all("li", class_="list-group-item"):
        if "bg-light" in (li.get("class") or []):
            continue                    # ligne d'en-tête
        row = li.find("div", class_="row")
        cols = row.find_all("div", recursive=False) if row is not None else []
        if len(cols) < 7:
            continue
        code = _code(cols[0])
        if code is None:
            continue
        serial = node_text(cols[5])
        rows.append(QueenListRow(code=code, name=node_text(cols[0]) or code, birth_year=parse_int(node_text(cols[1])),
                                 mother_code=_code(cols[2]), drone_mother_code=_code(cols[3]),
                                 cotation_count=parse_int(node_text(cols[4])),
                                 serial_raw=None if serial in ("", "None") else serial))
    return QueensListPage(breeder_id=int(_BREEDER.match(breeder["href"]).group(1)), trigram=node_text(breeder),
                          country_id=int(_COUNTRY.match(country["href"]).group(1)) if country is not None else None,
                          years=tuple(years), rows=tuple(rows))
