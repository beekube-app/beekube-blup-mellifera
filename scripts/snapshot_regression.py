"""Écrit l'instantané de non-régression tests/fixtures/expected-example-v2.json.

À relancer uniquement lors d'un changement de méthode volontaire (et à justifier dans le commit) :
    .venv/bin/python scripts/snapshot_regression.py
"""
from __future__ import annotations

import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from blup2.engine import run_blup  # noqa: E402
from blup2.schemas import BlupV2Input  # noqa: E402

SRC = os.path.join(ROOT, "tests", "fixtures", "input-example-v2.json")
DST = os.path.join(ROOT, "tests", "fixtures", "expected-example-v2.json")


def snapshot(out) -> dict:
    return {
        "engine_version": out.engine_version,
        "parameters": {k: {"source": p.source, "h2": round(p.h2, 4), "n_queens": p.n_queens,
                           "n_records": p.n_records} for k, p in out.parameters.items()},
        "queens": {str(q.id): {"total_index": q.total_index, "inbreeding": q.inbreeding,
                               "indexes": {c: [ci.status, ci.index, ci.reliability]
                                           for c, ci in q.indexes.items()}} for q in out.queens},
        "warnings": sorted({w.code for w in out.warnings}),
        "n_components": out.connectivity.n_components,
    }


def main() -> int:
    with open(SRC, encoding="utf-8") as f:
        inp = BlupV2Input.model_validate(json.load(f))
    with open(DST, "w", encoding="utf-8") as f:
        json.dump(snapshot(run_blup(inp)), f, ensure_ascii=False, indent=1, sort_keys=True)
        f.write("\n")
    print(f"{DST} écrit")
    return 0


if __name__ == "__main__":
    sys.exit(main())
