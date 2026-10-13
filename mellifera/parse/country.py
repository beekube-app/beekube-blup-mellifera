"""Page `/countrybyyear/<pays>/` : éleveurs du pays et leurs saisons de pedigree."""
from __future__ import annotations

import re
from dataclasses import dataclass

from .common import ParseError, make_soup, node_text

_BREEDER = re.compile(r"^/breeder/(\d+)/$")
_SEASON = re.compile(r"^/breeder_Year/(\d+)/(\d{4})/$")


@dataclass(frozen=True)
class BreederEntry:
    id: int
    trigram: str
    season_years: tuple[int, ...]       # années des pages /breeder_Year/<id>/<année>/, décroissantes


def parse_country(body: bytes) -> tuple[BreederEntry, ...]:
    soup = make_soup(body)
    if soup.find("ul") is None and soup.find("li") is None:
        raise ParseError("liste des éleveurs absente")
    out: list[BreederEntry] = []
    for li in soup.find_all("li", class_="list-group-item"):
        a = li.find("a", href=_BREEDER)
        if a is None:
            continue
        years = {int(_SEASON.match(s["href"]).group(2)) for s in li.find_all("a", href=_SEASON)}
        out.append(BreederEntry(int(_BREEDER.match(a["href"]).group(1)), node_text(a), tuple(sorted(years, reverse=True))))
    return tuple(out)
