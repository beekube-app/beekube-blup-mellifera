"""Outils communs aux parseurs : erreur, champs étiquetés, codes de reine, dates et nombres."""
from __future__ import annotations

import datetime as dt
import re
from dataclasses import dataclass
from typing import Optional

from bs4 import BeautifulSoup, Comment, NavigableString, Tag

PARSER_VERSION = 1

# Code international : PAYS-ÉLEVEUR-RACE-NUMÉRO-ANNÉE (FR-CBP-B-0276-2022, LU-PJ-V-005H-2015). Formes réelles
# du site : éleveur avec minuscules (ABg) ou réduit à « - » (vieilles reines allemandes : DE---B-0072-1993), race à
# plusieurs lettres (Msw), numéro avec minuscules (082b, 0mix). La casse compte : elle est conservée telle quelle.
CODE_RE = re.compile(r"^[A-Z]{2}-(?:[A-Za-z0-9]{1,6}|-)-[A-Za-z0-9]{1,3}-[0-9A-Za-z]{1,6}-\d{4}$")
_QUEEN_HREF = re.compile(r"^/queen/([^/]+)/$")
_HEADINGS = {"h1", "h2", "h3", "h4", "h5", "h6"}
_LABEL_TRIM = re.compile(r"^[\s:\-–]+|[\s:\-–]+$")
_MONTHS = {"jan": 1, "feb": 2, "mar": 3, "march": 3, "apr": 4, "april": 4, "may": 5, "jun": 6, "june": 6,
           "jul": 7, "july": 7, "aug": 8, "sep": 9, "sept": 9, "oct": 10, "nov": 11, "dec": 12}


class ParseError(ValueError):
    """La page ne contient pas la structure attendue."""


def make_soup(body: bytes) -> BeautifulSoup:
    """Page complète uniquement : vide ou sans « </html> » final (réponse tronquée) → ParseError."""
    if not body or not body.strip():
        raise ParseError("page vide")
    if b"</html>" not in body[-200:].lower():
        raise ParseError("page tronquée (pas de </html> final)")
    return BeautifulSoup(body, "html.parser")


def norm(s: Optional[str]) -> str:
    """Espaces (y compris insécables) réduits à un seul, bords retirés."""
    return re.sub(r"\s+", " ", s or "").strip()


def node_text(node: Tag) -> str:
    return norm(node.get_text(" "))


def queen_code_from_href(href: Optional[str]) -> Optional[str]:
    m = _QUEEN_HREF.match(href or "")
    return m.group(1) if m else None


def is_code(s: Optional[str]) -> bool:
    return bool(CODE_RE.match(s or ""))


def parse_django_date(s: Optional[str]) -> Optional[dt.date]:
    """« Sept. 7, 2020 », « June 6, 2024 », « May 3, 2021 » → date ; sinon None."""
    m = re.fullmatch(r"([A-Za-z]+)\.?\s+(\d{1,2}),\s+(\d{4})", norm(s))
    if not m:
        return None
    month = _MONTHS.get(m.group(1).lower())
    if month is None:
        return None
    try:
        return dt.date(int(m.group(3)), month, int(m.group(2)))
    except ValueError:
        return None


def parse_number(s: Optional[str]) -> Optional[float]:
    t = norm(s).replace(",", ".").rstrip("%")
    if t in ("", "-", "None"):
        return None
    try:
        return float(t)
    except ValueError:
        return None


def parse_int(s: Optional[str]) -> Optional[int]:
    v = parse_number(s)
    return int(v) if v is not None and v.is_integer() else None


@dataclass(frozen=True)
class LabelledField:
    """Un <strong>libellé</strong> suivi de sa valeur (texte et liens jusqu'au <strong> suivant)."""

    label: str
    text: str
    hrefs: tuple[str, ...]


def labelled_fields(container: Tag) -> list[LabelledField]:
    """Champs « <strong>Libellé : </strong> valeur » d'un bloc, dans l'ordre du document.

    Indépendant de l'imbrication : les <p> mal fermés du site ne changent pas le découpage.
    Le texte des <script> et <style> est ignoré.
    """
    elems = list(container.descendants)
    starts = [i for i, e in enumerate(elems) if isinstance(e, Tag) and e.name == "strong"]
    out: list[LabelledField] = []
    for k, i in enumerate(starts):
        strong = elems[i]
        own = {id(x) for x in strong.descendants}
        end = starts[k + 1] if k + 1 < len(starts) else len(elems)
        parts: list[str] = []
        hrefs: list[str] = []
        for e in elems[i + 1:end]:
            if isinstance(e, Tag) and e.name in _HEADINGS:
                break                      # un titre de section clôt le champ (« Details », « Cotations »)
            if id(e) in own or isinstance(e, Comment):
                continue
            if isinstance(e, NavigableString):
                if e.parent is not None and e.parent.name in ("script", "style"):
                    continue
                parts.append(str(e))
            elif e.name == "a" and e.get("href"):
                hrefs.append(e["href"])
        label = _LABEL_TRIM.sub("", norm(strong.get_text(" ")))
        out.append(LabelledField(label, _LABEL_TRIM.sub("", norm(" ".join(parts))), tuple(hrefs)))
    return out


def field_after(fields: list[LabelledField], anchor: str, label: str,
                stop: Optional[str] = None) -> Optional[LabelledField]:
    """Premier champ `label` situé après le champ `anchor` (et avant `stop` s'il est donné).

    `label == anchor` renvoie le champ `anchor` lui-même.
    """
    seen = False
    for f in fields:
        if not seen:
            if f.label == anchor:
                seen = True
                if label == anchor:
                    return f
            continue
        if stop is not None and f.label == stop:
            return None
        if f.label == label:
            return f
    return None
