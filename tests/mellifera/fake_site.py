"""Pages synthétiques au balisage d'Index Mellifera et faux site (httpx.MockTransport) pour les tests."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import httpx

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(name: str) -> bytes:
    return (FIXTURES / name).read_bytes()


# ------------------------------------------------------------------ faux site et horloge
class FakeSite:
    """Sert des pages par chemin. `failures[chemin]` : réponses d'erreur renvoyées avant la vraie page, code HTTP
    ou (code, Retry-After)."""

    def __init__(self, pages: Optional[dict[str, bytes]] = None):
        self.pages: dict[str, bytes] = dict(pages or {})
        self.failures: dict[str, list] = {}
        self.redirects: dict[str, str] = {}
        self.requests: list[httpx.Request] = []

    @property
    def paths(self) -> list[str]:
        return [r.url.path for r in self.requests]

    def handler(self, request: httpx.Request) -> httpx.Response:
        self.requests.append(request)
        path = request.url.path
        if request.method != "GET":
            return httpx.Response(405)
        pending = self.failures.get(path)
        if pending:
            status, retry_after = (pending.pop(0), None) if isinstance(pending[0], int) else pending.pop(0)
            return httpx.Response(status, headers={"Retry-After": retry_after} if retry_after else {})
        if path in self.redirects:
            return httpx.Response(302, headers={"Location": self.redirects[path]})
        if path in self.pages:
            return httpx.Response(200, content=self.pages[path])
        return httpx.Response(404)

    def transport(self) -> httpx.MockTransport:
        return httpx.MockTransport(self.handler)


class FakeClock:
    """Horloge simulée : sleep() avance le temps et note chaque attente."""

    def __init__(self):
        self.now = 1000.0
        self.sleeps: list[float] = []

    def clock(self) -> float:
        return self.now

    def sleep(self, seconds: float) -> None:
        self.sleeps.append(round(seconds, 6))
        self.now += seconds


# ------------------------------------------------------------------ pages synthétiques (balisage du site)
def _page(body: str) -> bytes:
    return f"<!DOCTYPE html><html><head><title>Index Mellifera</title></head><body>{body}</body></html>".encode()


def countries_page(countries: list[tuple[int, str, Optional[str], int]]) -> bytes:
    """[(id, nom, ISO ou None, nombre d'éleveurs)]"""
    blocks = "".join(
        f'<div class="col-lg-2"><a href="/countrybyyear/{cid}/">'
        f'<img src="/media/flags/{iso or "XX_abc"}.png" alt="flag"></a><h3>{name}</h3><h3>{n} Breeders</h3></div>'
        for cid, name, iso, n in countries)
    return _page(f'<div class="row">{blocks}</div>')


def _season_links(breeder_id: int, years: list[int]) -> str:
    return "".join(f'<a href="/breeder_Year/{breeder_id}/{y}/">{y},</a>' for y in years)


def country_page(breeders: list[tuple[int, str, list[int]]]) -> bytes:
    """[(id éleveur, trigramme, années de saison)]"""
    rows = "".join(
        f'<li class="list-group-item"><div class="row"><div class="col-3"><a href="/breeder/{bid}/">{tri}</a></div>'
        f'<div class="col-7">{_season_links(bid, years)}</div>'
        f'<div class="col-1"><a href="/profile/{bid}/">Profil</a></div></div></li>'
        for bid, tri, years in breeders)
    head = '<li class="list-group-item bg-light"><div class="row"><div class="col-3"><span>Breeder</span></div></div></li>'
    return _page(f'<ul class="list-group">{head}{rows}</ul>')


@dataclass
class Q:
    """Reine synthétique : sert à la fois aux listes et aux fiches."""
    code: str
    year: int
    breeder_id: int
    name: Optional[str] = None
    mother: Optional[str] = None
    drones: Optional[str] = None
    number: Optional[int] = None
    cotations: int = 0
    mating: tuple[str, str, str] = ("ins", "Instrumental insemination", "Instrumental insemination (ins)")
    radar: dict[int, list[float]] = field(default_factory=dict)   # année → 5 valeurs

    @property
    def label(self) -> str:
        return self.name or self.code


def _link(code: Optional[str], empty: str) -> str:
    return f'<a href="/queen/{code}/">{code}</a>' if code else empty


def queens_list_page(breeder_id: int, trigram: str, country_id: int, years: list[int], queens: list[Q]) -> bytes:
    nav = "".join(f'<a href="/breederbyyear/{breeder_id}/{y}/">{y}</a> ' for y in years)
    head = ('<li class="list-group-item bg-light"><div class="row">' + "".join(
        f'<div class="col">{h}</div>' for h in ("Name", "Birth", "Mother", "Drones", "Cotation", "Serial", "Queen"))
        + "</div></li>")
    rows = "".join(
        '<li class="list-group-item"><div class="row">'
        f'<div class="col"><a href="/queen/{q.code}/">{q.label}</a></div><div class="col">{q.year}</div>'
        f'<div class="col">{_link(q.mother, "Unknowne")}</div><div class="col">{_link(q.drones, "Unknown")}</div>'
        f'<div class="col">{q.cotations}</div><div class="col">None</div>'
        f'<div class="col"><a href="/queen/{q.code}/">{q.code}</a></div></div></li>'
        for q in queens)
    crumbs = f'<a href="/countrybyyear/{country_id}/">FR</a> / <a href="/breeder/{breeder_id}/">{trigram}</a>'
    return _page(f"{crumbs}<div>{nav}</div><ul>{head}{rows}</ul>")


def queen_page(q: Q) -> bytes:
    radar = ""
    if q.radar:
        sets = "".join(f"{{ label: '{y}', data: [{', '.join(f'{v:g}' for v in vals)}] }}," for y, vals in q.radar.items())
        radar = ("<script>var marksData = { labels: ['swarming', 'mildness', 'frameHolding', 'wintering', "
                 f"'broodVitality'], datasets: [{sets}] }};</script>")
    code, desc, title = q.mating
    ratings = "".join(f'<a href="/rating/{i + 1}" class="list-group-item">{q.code} - 01/07/2023 - Testeur anonyme</a>'
                      for i in range(q.cotations))
    return _page(
        f'<div class="cadree"><h1>{q.label}</h1>'
        f"<p><strong>International Code : </strong>{q.code}</p>"
        f"<p><strong>Birth Year : </strong>{q.year}</p>"
        f"<p><strong>Mother : </strong>{_link(q.mother, 'Unknown')}</p>"
        f"<p><strong>Drones : </strong>{_link(q.drones, 'Unknown')}"
        f"<strong> -- Number : </strong>{q.number if q.number is not None else 'None'}</p>"
        f'<p><strong>Breeder : </strong><a href="/breederbyyear/{q.breeder_id}/{q.year}/">Queens list {q.year}.</a></p>'
        f"<h1>Details</h1>"
        f'<p title="{title}"><strong>Insemination Code : </strong>{code}<strong> -- Description </strong>{desc}</p>'
        f"<p><strong>Declared Race : </strong>No<strong> -- Race Code : </strong>B"
        f"<strong> -- Lineage : </strong>Unknown<strong> -- Description </strong>buckfast</p>"
        f"<p><strong>Import : </strong>no<strong> -- Comment : </strong>No</p>{radar}</div>"
        f'<div class="list-group">{ratings}</div>')


@dataclass
class E:
    """Ligne d'un tableau d'évaluations."""
    queen: str                          # texte de la cellule Queen
    year: int
    values: dict[str, str]
    mode: str = "evjp"
    tester: str = "CBP"
    tester_id: int = 145
    date: str = "July 6, 2023"


def season_page(evaluations: list[E], links: list[tuple[str, str]] = (),
                legend: list[tuple[str, str, str]] = ()) -> bytes:
    """links : [(code, texte du lien)] présents sur la page ; legend : [(code, description, title)]."""
    anchors = "".join(f'<a href="/queen/{c}/">{t}</a> ' for c, t in links)
    tables = ""
    for mode in dict.fromkeys(e.mode for e in evaluations):
        rows = [e for e in evaluations if e.mode == mode]
        crit = list(dict.fromkeys(k for e in rows for k in e.values))
        head = "".join(f'<td class="tdh01">{h}</td>' for h in
                       ["Queen", "BirthYear", "Tester", "YearEva", "tester", "year", *crit, "Comment"])
        body = "".join(
            "<tr>" + "".join(f'<td class="td01">{v}</td>' for v in
                             [e.queen, e.year, e.tester, e.date[-4:], e.tester_id, e.date,
                              *(e.values.get(c, "-") for c in crit), ""]) + "</tr>"
            for e in rows)
        tables += f"<h5>mode : {mode}</h5><table><tr>{head}</tr>{body}</table>"
    leg = "".join(f'<div class="row"><div class="col-2" title="{t}">{c}</div><div class="col-10">{d}</div></div>'
                  for c, d, t in legend)
    return _page(f'{anchors}<h3>Evaluations</h3>{tables}<div id="insemCodeLegend">{leg}</div>')
