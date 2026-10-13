"""Page `/all/` : liste des pays."""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Optional

from .common import ParseError, make_soup, node_text, parse_int

_HREF = re.compile(r"^/countrybyyear/(\d+)/$")
_FLAG = re.compile(r"/flags/([A-Za-z]{2})(?:_\w+)?\.\w+$")   # « DE.png », « RE_HNhsPa6.png »


@dataclass(frozen=True)
class CountryEntry:
    id: int
    name: str
    iso: Optional[str]                  # None pour « Unknown » (id 0), sans drapeau
    n_breeders: Optional[int]


def parse_countries(body: bytes) -> tuple[CountryEntry, ...]:
    soup = make_soup(body)
    out: dict[int, CountryEntry] = {}
    for a in soup.find_all("a", href=_HREF):
        cid = int(_HREF.match(a["href"]).group(1))
        box = a.parent
        img = a.find("img")
        m = _FLAG.search(img.get("src", "")) if img is not None else None
        h3 = box.find_all("h3") if box is not None else []
        name = node_text(h3[0]) if h3 else ""
        count = parse_int(node_text(h3[1]).split(" ")[0]) if len(h3) > 1 else None
        out.setdefault(cid, CountryEntry(cid, name, m.group(1).upper() if m else None, count))
    if not out:
        raise ParseError("aucun pays")
    return tuple(out.values())
