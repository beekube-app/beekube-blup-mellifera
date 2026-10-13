"""Chemins des pages publiques lues par le crawler : construction, type de page, liste blanche."""
from __future__ import annotations

import re
from typing import Optional

# Liste blanche : seuls ces chemins sont demandés au site (spec §4.3). Tout le reste est refusé avant l'envoi.
_KINDS = (
    ("countries", re.compile(r"^/all/$")),
    ("country", re.compile(r"^/countrybyyear/(\d+)/$")),
    ("queens_list", re.compile(r"^/breederbyyear/(\d+)/(\d{4})/$")),
    ("season", re.compile(r"^/breeder_Year/(\d+)/(\d{4})/$")),
    ("queen", re.compile(r"^/queen/([A-Za-z0-9()_-][A-Za-z0-9()._-]*)/$")),
)
KINDS = tuple(k for k, _ in _KINDS)        # ordre de relecture (reparse) : parents avant enfants
# Caractères admis dans un code mis en URL ; pas de point en tête (« .. » serait normalisé en « / » par httpx).
_SAFE_CODE = re.compile(r"^[A-Za-z0-9()_-][A-Za-z0-9()._-]*$")


def kind_of(path: str) -> Optional[str]:
    """Type de page d'un chemin autorisé, None s'il est hors liste blanche."""
    for kind, rx in _KINDS:
        if rx.match(path):
            return kind
    return None


def path_ids(path: str) -> tuple[str, ...]:
    """Identifiants contenus dans le chemin (« /breeder_Year/145/2022/ » → ('145', '2022'))."""
    for _, rx in _KINDS:
        m = rx.match(path)
        if m:
            return m.groups()
    raise ValueError(f"chemin hors liste blanche : {path}")


def countries_path() -> str:
    return "/all/"


def country_path(country_id: int) -> str:
    return f"/countrybyyear/{country_id}/"


def queens_list_path(breeder_id: int, year: int) -> str:
    return f"/breederbyyear/{breeder_id}/{year}/"


def season_path(breeder_id: int, year: int) -> str:
    return f"/breeder_Year/{breeder_id}/{year}/"


def queen_path(code: str) -> Optional[str]:
    """None si le code contient des caractères qu'on ne met pas dans une URL (compté par le crawler)."""
    return f"/queen/{code}/" if _SAFE_CODE.match(code or "") else None
