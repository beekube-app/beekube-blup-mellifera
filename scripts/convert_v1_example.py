"""Convertit data/input-example.json (entrée v1 de POST /blup) en fixture v2 : tests/fixtures/input-example-v2.json.

    .venv/bin/python scripts/convert_v1_example.py

Règles :
- critères : `evaluate` -> note, poids 1, sens +1 ; `evaluate_elimination` -> elimination = true ;
- reines : born = int, mother = queenbee_parent ; drone_parent (reine dont les mâles ont fécondé
  cette reine, spec §1) -> insémination avec cette seule RPM, sinon fécondation inconnue ;
- notes v1 = note / notationMax : grille en quarts (mars-mai) -> value = note × 4 sur [0, 4],
  grille en cinquièmes -> value = note × 5 sur [0, 5] (conversion sans perte), sinon [0, 1] ;
- date = 20YY-MM-01, evaluator = user, hive_class : dadant-6 -> nuc, dadant-div-10 -> split, sinon hive.
"""
from __future__ import annotations

import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "data", "input-example.json")
DST = os.path.join(ROOT, "tests", "fixtures", "input-example-v2.json")
HIVE_CLASS = {"dadant-6": "nuc", "dadant-div-10": "split"}


def _is_multiple(x: float, k: int) -> bool:
    return abs(x * k - round(x * k)) < 1e-9


def scale_of(note: float, month: int) -> tuple[float, float, float]:
    """(value, scale_min, scale_max) restituant exactement la note v1 après normalisation."""
    quarters, fifths = _is_multiple(note, 4), _is_multiple(note, 5)
    if quarters and (month <= 5 or not fifths):
        return round(note * 4), 0.0, 4.0
    if fifths:
        return round(note * 5), 0.0, 5.0
    return note, 0.0, 1.0


def convert(v1: dict) -> dict:
    crits = [{"id": c, "name": f"critère {c}", "direction": 1, "kind": "note", "weight": 1,
              "elimination": False} for c in v1["evaluate"]]
    crits += [{"id": c, "name": f"critère {c}", "direction": 1, "kind": "note", "weight": 1,
               "elimination": True} for c in v1["evaluate_elimination"]]
    queens, records = [], []
    for q in v1["data"]:
        mating = {"type": "unknown"}
        if q.get("drone_parent") is not None:
            mating = {"type": "insemination", "drone_queens": [q["drone_parent"]]}
        queens.append({"id": q["queenbee"], "born": int(q["born"]), "mother": q.get("queenbee_parent"),
                       "mating": mating, "variety": q.get("variety"), "died": bool(q.get("died", False))})
        for c in crits:
            for e in q.get("evaluate", {}).get(str(c["id"]), []):
                month = int(e["month"])
                value, smin, smax = scale_of(float(e["note"]), month)
                records.append({"queen": q["queenbee"], "criterion": c["id"], "value": value,
                                "scale_min": smin, "scale_max": smax,
                                "date": f"20{e['year']}-{month:02d}-01", "evaluator": e["user"],
                                "apiary": e["apiary"], "hive_type": e["beehiveType"],
                                "hive_class": HIVE_CLASS.get(e["beehiveType"], "hive")})
    return {"exploitation": v1["exploitation"], "profile": "standard", "options": {},
            "criteria": crits, "queens": queens, "records": records}


def main() -> int:
    with open(SRC, encoding="utf-8") as f:
        v2 = convert(json.load(f))
    os.makedirs(os.path.dirname(DST), exist_ok=True)
    with open(DST, "w", encoding="utf-8") as f:
        json.dump(v2, f, ensure_ascii=False, indent=1)
        f.write("\n")
    print(f"{DST} : {len(v2['queens'])} reines, {len(v2['records'])} notes")
    return 0


if __name__ == "__main__":
    sys.exit(main())
