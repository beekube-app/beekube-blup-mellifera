"""Produit les fixtures de test anonymisées à partir des pages brutes d'Index Mellifera.

    .venv/bin/python scripts/make_mellifera_fixtures.py tests/mellifera/fixtures/raw/index-mellifera tests/mellifera/fixtures

Les pages brutes (ignorées par git) contiennent des noms de personnes : auteurs des cotations, commentaire
d'import, commentaires des tableaux d'évaluations. Ces trois zones sont neutralisées ; le reste de la page est
conservé octet pour octet, pour que les tests portent sur le vrai balisage du site.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

PAGES = {   # page brute → fixture commitée
    "all_.html": "all.html",
    "countrybyyear_8_.html": "countrybyyear_8.html",
    "breederbyyear_145_2022.html": "breederbyyear_145_2022.html",
    "queen_FR-CBP-B-0276-2022_.html": "queen_FR-CBP-B-0276-2022.html",
    "queen_FR-NG-B-4536-2024_.html": "queen_FR-NG-B-4536-2024.html",
    "breeder_Year_145_2021.html": "breeder_Year_145_2021.html",
    "breeder_Year_145_2022.html": "breeder_Year_145_2022.html",
    "breeder_Year_984_2026_.html": "breeder_Year_984_2026.html",
}
ANONYMOUS = "Testeur anonyme"
_RATING = re.compile(r'(<a href="/rating/\d+"[^>]*>)((?:(?!</a>).)*)(</a>)', re.S)
_IMPORT_COMMENT = re.compile(r"(<strong>Import : </strong>(?:(?!<strong>).)*<strong> -- Comment : </strong>)"
                             r"((?:(?!<strong>|</p>).)*)", re.S)
_LAST_CELL = re.compile(r'(<td class="td\d+">)((?:(?!</td>).)*)(</td></tr>)', re.S)


def _rating(m: re.Match) -> str:
    parts = [p.strip() for p in re.sub(r"\s+", " ", m.group(2)).split(" - ")]
    return f"{m.group(1)}{' - '.join(parts[:2] + [ANONYMOUS])}{m.group(3)}"


def anonymize(html: str) -> str:
    html = _RATING.sub(_rating, html)
    html = _IMPORT_COMMENT.sub(lambda m: f"{m.group(1)}\n    No\n    ", html)
    return _LAST_CELL.sub(lambda m: f"{m.group(1)}{m.group(3)}", html)   # colonne Comment des tableaux


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print(__doc__)
        return 2
    raw, out = Path(argv[1]), Path(argv[2])
    out.mkdir(parents=True, exist_ok=True)
    for src, dst in PAGES.items():
        (out / dst).write_text(anonymize((raw / src).read_text(encoding="utf-8")), encoding="utf-8")
        print(f"{src} → {out / dst}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
