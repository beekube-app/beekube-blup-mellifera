"""Convertit l'ancien data/input-example.json au format de POST /blup et l'anonymise.

Usage : python3 scripts/anonymize_example.py data/input-example.json data/input-example.json

Les identifiants (exploitation, reines et parents, évaluateurs, ruchers) sont renumérotés ;
la table de correspondance n'est ni écrite ni affichée. Notes, mois, années, types de ruche
et liens de parenté sont conservés.
"""
import json
import sys


def anonymize(src):
    criteria = [str(c) for c in src["evaluate"]]
    elimination = [str(c) for c in src["evaluate_elimination"]]
    keys = criteria + elimination
    queens = src["data"]
    if queens and "evaluate" in queens[0]:
        raise SystemExit("Le fichier est déjà au format de l'API.")

    def records(queen):
        for key in keys:
            for rec in queen.get(key, []):
                yield rec

    queen_ids = sorted(
        {q["queenbee"] for q in queens}
        | {q[k] for q in queens for k in ("queenbee_parent", "drone_parent") if q.get(k) is not None}
    )
    user_ids = sorted({r["user"] for q in queens for r in records(q)})
    apiary_ids = sorted(
        {r["apiary"] for q in queens for r in records(q)}
        | {q["apiary"] for q in queens if q.get("apiary") is not None}
    )
    queen_map = {old: i + 1 for i, old in enumerate(queen_ids)}
    user_map = {old: i + 1 for i, old in enumerate(user_ids)}
    apiary_map = {old: 101 + i for i, old in enumerate(apiary_ids)}

    def convert_record(rec):
        return {
            "note": rec["note"],
            "user": user_map[rec["user"]],
            "apiary": apiary_map[rec["apiary"]],
            "beehiveType": rec["beehiveType"],
            "month": rec["month"],
            "year": rec["year"],
        }

    data = []
    for q in queens:
        data.append({
            "queenbee": queen_map[q["queenbee"]],
            "queenbee_parent": queen_map.get(q.get("queenbee_parent")),
            "drone_parent": queen_map.get(q.get("drone_parent")),
            "apiary_default": apiary_map.get(q.get("apiary")),
            "born": str(q["born"]),
            "variety": "Buckfast",
            "died": False,
            "evaluate": {key: [convert_record(r) for r in q.get(key, [])] for key in keys},
        })

    return {
        "exploitation": 1,
        "evaluate": src["evaluate"],
        "evaluate_elimination": src["evaluate_elimination"],
        "data": data,
    }


def main(argv):
    if len(argv) != 3:
        print(__doc__)
        return 2
    with open(argv[1], encoding="utf-8") as f:
        src = json.load(f)
    converted = anonymize(src)
    with open(argv[2], "w", encoding="utf-8") as f:
        json.dump(converted, f, indent=2, ensure_ascii=False)
        f.write("\n")
    print(f"{len(converted['data'])} reines écrites dans {argv[2]}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
