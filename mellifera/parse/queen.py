"""Page `/queen/<code>/` : fiche d'une reine (pedigree, accouplement, consanguinité, graphique des cotations)."""
from __future__ import annotations

import datetime as dt
import re
from dataclasses import dataclass
from typing import Optional

from .common import (LabelledField, ParseError, field_after, labelled_fields, make_soup, node_text, norm,
                     parse_django_date, parse_int, parse_number, queen_code_from_href)

_BREEDER = re.compile(r"^/breederbyyear/(\d+)/\d{4}/$")
_RATING = re.compile(r"^/rating/\d+$")
_RADAR_LABELS = re.compile(r"labels:\s*\[([^\]]*)\]")
_RADAR_SET = re.compile(r"label:\s*'(\d{4})'.*?data:\s*\[([^\]]*)\]", re.S)
_EMPTY = ("", "No", "None", "Unknown")


@dataclass(frozen=True)
class RadarSeries:
    year: int
    values: dict[str, float]            # libellé du graphique (« frameHolding ») → valeur


@dataclass(frozen=True)
class QueenPage:
    code: str
    name: str
    birth_year: Optional[int]
    breeder_id: Optional[int]
    mother_code: Optional[str]
    mother_comment: Optional[str]
    drone_mother_code: Optional[str]    # « Drones : X » : mère des reines productrices de mâles
    n_drone_queens: Optional[int]       # « Number : n »
    drone_comment: Optional[str]
    single_drone: bool                  # suffixe « 1dr » du nom
    generation_raw: Optional[str]
    coi: Optional[float]
    coi_reliability: Optional[float]    # fraction (« 70.6% » → 0.706)
    not_found_raw: Optional[str]
    serial_raw: Optional[str]
    insemination_date: Optional[dt.date]
    mating_code: Optional[str]
    mating_description: Optional[str]
    mating_category_raw: Optional[str]  # title du paragraphe : « Instrumental insemination (FR) »
    race_code: Optional[str]
    race_description: Optional[str]
    lineage: Optional[str]
    origin_raw: Optional[str]
    data_origin_raw: Optional[str]
    cotation_count: int
    radar: tuple[RadarSeries, ...]


def _value(f: Optional[LabelledField]) -> Optional[str]:
    return None if f is None or f.text in _EMPTY else f.text


def _queen_code(f: Optional[LabelledField]) -> Optional[str]:
    for h in f.hrefs if f is not None else ():
        code = queen_code_from_href(h)
        if code:
            return code
    return None


def _radar(soup) -> tuple[RadarSeries, ...]:
    for script in soup.find_all("script"):
        src = script.string or ""
        m = _RADAR_LABELS.search(src)
        if m is None or "datasets" not in src:
            continue
        labels = re.findall(r"'([^']+)'", m.group(1))
        out = []
        for year, data in _RADAR_SET.findall(src):
            values = {lab: v for lab, v in zip(labels, (parse_number(x) for x in data.split(","))) if v is not None}
            if values:
                out.append(RadarSeries(int(year), values))
        return tuple(out)
    return ()


def _mating_category(card) -> Optional[str]:
    for strong in card.find_all("strong"):
        if norm(strong.get_text(" ")).rstrip(" :") == "Insemination Code":
            p = strong.find_parent("p")
            return norm(p.get("title")) or None if p is not None else None
    return None


def parse_queen(body: bytes) -> QueenPage:
    soup = make_soup(body)
    card = soup.find("div", class_="cadree")
    if card is None:
        raise ParseError("fiche absente")
    fields = labelled_fields(card)
    code_f = field_after(fields, "International Code", "International Code")
    if code_f is None or not code_f.text:
        raise ParseError("code international absent")
    h1 = card.find("h1")
    name = node_text(h1) if h1 is not None else ""
    name = name or code_f.text
    breeder_f = field_after(fields, "Breeder", "Breeder")
    breeder_id = next((int(m.group(1)) for h in (breeder_f.hrefs if breeder_f else ()) if (m := _BREEDER.match(h))),
                      None)
    reliability = parse_number(_value(field_after(fields, "COI", "Fiability")))
    n_dq = parse_int(_value(field_after(fields, "Drones", "Number", stop="Breeder")))
    return QueenPage(
        code=code_f.text,
        name=name,
        birth_year=parse_int(_value(field_after(fields, "Birth Year", "Birth Year"))),
        breeder_id=breeder_id,
        mother_code=_queen_code(field_after(fields, "Mother", "Mother")),
        mother_comment=_value(field_after(fields, "Mother", "Comment", stop="Drones")),
        drone_mother_code=_queen_code(field_after(fields, "Drones", "Drones")),
        n_drone_queens=n_dq if n_dq is not None and n_dq >= 1 else None,
        drone_comment=_value(field_after(fields, "Drones", "Comment", stop="Breeder")),
        single_drone=name.endswith("1dr"),
        generation_raw=_value(field_after(fields, "Generation", "Generation")),
        coi=parse_number(_value(field_after(fields, "COI", "COI"))),
        coi_reliability=round(reliability / 100, 6) if reliability is not None else None,
        not_found_raw=_value(field_after(fields, "COI", "Not found")),
        serial_raw=_value(field_after(fields, "General comment", "In a serial")),
        insemination_date=parse_django_date(_value(field_after(fields, "Insemination Date", "Insemination Date"))),
        mating_code=_value(field_after(fields, "Insemination Code", "Insemination Code")),
        mating_description=_value(field_after(fields, "Insemination Code", "Description", stop="Declared Race")),
        mating_category_raw=_mating_category(card),
        race_code=_value(field_after(fields, "Declared Race", "Race Code")),
        race_description=_value(field_after(fields, "Lineage", "Description", stop="Import")),
        lineage=_value(field_after(fields, "Declared Race", "Lineage")),
        origin_raw=_value(field_after(fields, "Import", "Geographic")),
        data_origin_raw=_value(field_after(fields, "Data origin", "Data origin")),
        cotation_count=len(soup.find_all("a", href=_RATING)),
        radar=_radar(soup),
    )
