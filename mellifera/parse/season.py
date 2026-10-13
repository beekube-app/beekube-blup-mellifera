"""Page `/breeder_Year/<éleveur>/<année>/` : tableaux d'évaluations et légende des codes de fécondation."""
from __future__ import annotations

import datetime as dt
from collections import defaultdict
from dataclasses import dataclass
from typing import Optional

from .common import ParseError, is_code, make_soup, node_text, norm, parse_django_date, parse_int, queen_code_from_href

# Colonnes descriptives ; toutes les autres sont des critères. « Comment » n'est pas conservée.
FIXED_COLUMNS = frozenset({"Queen", "BirthYear", "Tester", "YearEva", "tester", "year", "Comment"})


@dataclass(frozen=True)
class SeasonEvaluation:
    queen_name_raw: str
    queen_birth_year: Optional[int]
    queen_code: Optional[str]           # résolu sur la page (lien de même texte, sinon texte au format code)
    mode: str                           # protocole du tableau : « evjp », « ceta », « evdk »…
    tester: str                         # trigramme du testeur
    tester_breeder_id: Optional[int]
    eval_year: Optional[int]
    eval_date: Optional[dt.date]
    values: dict[str, str]              # critère → valeur brute ; cellules « - » et vides exclues


@dataclass(frozen=True)
class MatingCodeEntry:
    code: str
    description: str
    category_raw: Optional[str]         # title de la légende : « Mating station (FR) »


@dataclass(frozen=True)
class SeasonPage:
    evaluations: tuple[SeasonEvaluation, ...]
    mating_codes: tuple[MatingCodeEntry, ...]


def _link_codes(soup) -> dict[str, set[str]]:
    found: dict[str, set[str]] = defaultdict(set)
    for a in soup.find_all("a", href=True):
        code = queen_code_from_href(a["href"])
        if code:
            found[node_text(a)].add(code)
    return found


def _resolve(name: str, links: dict[str, set[str]]) -> Optional[str]:
    codes = links.get(name, set())
    if len(codes) == 1:
        return next(iter(codes))        # 1. lien de la page au texte identique (prime sur le format)
    return name if is_code(name) and not codes else None   # 2. texte au format code


def parse_season(body: bytes) -> SeasonPage:
    soup = make_soup(body)
    links = _link_codes(soup)
    evaluations: list[SeasonEvaluation] = []
    for h5 in soup.find_all("h5"):
        title = node_text(h5)
        if not title.startswith("mode") or ":" not in title:
            continue
        mode = norm(title.split(":", 1)[1])
        table = h5.find_next("table")
        rows = table.find_all("tr") if table is not None else []
        if not rows:
            continue
        header = [node_text(c) for c in rows[0].find_all(["td", "th"])]
        if "Queen" not in header:
            raise ParseError(f"tableau « {mode} » sans colonne Queen")
        for tr in rows[1:]:
            cells = [node_text(c) for c in tr.find_all(["td", "th"])]
            if len(cells) != len(header):
                raise ParseError(f"tableau « {mode} » : {len(cells)} cellules pour {len(header)} colonnes")
            row = dict(zip(header, cells))
            name = row["Queen"]
            if not name:
                continue
            evaluations.append(SeasonEvaluation(
                queen_name_raw=name,
                queen_birth_year=parse_int(row.get("BirthYear")),
                queen_code=_resolve(name, links),
                mode=mode,
                tester=row.get("Tester", ""),
                tester_breeder_id=parse_int(row.get("tester")),
                eval_year=parse_int(row.get("YearEva")),
                eval_date=parse_django_date(row.get("year")),
                values={k: v for k, v in row.items() if k not in FIXED_COLUMNS and v not in ("", "-")},
            ))
    legend = soup.find(id="insemCodeLegend")
    codes: list[MatingCodeEntry] = []
    for row in legend.find_all("div", class_="row") if legend is not None else []:
        cols = row.find_all("div", recursive=False)
        if len(cols) >= 2 and node_text(cols[0]):
            codes.append(MatingCodeEntry(node_text(cols[0]), node_text(cols[1]), norm(cols[0].get("title")) or None))
    return SeasonPage(tuple(evaluations), tuple(codes))
