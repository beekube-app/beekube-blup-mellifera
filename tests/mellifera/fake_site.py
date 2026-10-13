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
