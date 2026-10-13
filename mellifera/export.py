"""Export d'une sélection de la base au format BlupV2Input (spec §7)."""
from __future__ import annotations

import datetime as dt
import sqlite3
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Iterable, Optional

from blup2.pedigree import find_cycle
from blup2.schemas import BORN_MIN, BlupV2Input, QueenIn

SOURCES = ("season_table", "radar")
LABELS = {
    "invalid_birth_year": "reines sans année de naissance valide (exclues)",
    "self_reference_cut": "liens auto-référents coupés",
    "younger_parent_cut": "liens vers un parent plus jeune que la reine coupés",
    "cycle_link_cut": "liens coupés pour rompre un cycle",
    "mating_type_inferred": "accouplements « unknown » à côté mâle connu passés en insemination",
    "unmapped_criterion": "valeurs de critères non mappés (écartées)",
    "harvest": "récoltes en kg (écartées)",
    "non_numeric": "valeurs non numériques (écartées)",
    "no_date": "évaluations sans date (écartées)",
    "out_of_scale": "valeurs hors échelle (écartées)",
    "no_apiary": "évaluations sans testeur ni éleveur (écartées)",
    "radar_covered": "valeurs du graphique déjà couvertes par un tableau (écartées)",
    "merged_duplicates": "doublons fusionnés par moyenne",
}


@dataclass
class ExportReport:
    breeders: list[str] = field(default_factory=list)   # « CBP (145, France) »
    n_queens: int = 0
    n_selected: int = 0
    n_records: int = 0
    counts: Counter = field(default_factory=Counter)   # motifs d'exclusion et corrections, par nom

    def lines(self) -> list[str]:
        out = [f"éleveurs : {', '.join(self.breeders)}",
               f"reines : {self.n_queens} (dont {self.n_selected} sélectionnées)", f"enregistrements : {self.n_records}"]
        out += [f"  {LABELS.get(k, k)} : {v}" for k, v in sorted(self.counts.items())]
        return out


def _valid_year(y: Optional[int]) -> bool:
    return y is not None and BORN_MIN <= y <= dt.date.today().year + 1


def _breeders(conn: sqlite3.Connection, refs: Iterable[str]) -> list[sqlite3.Row]:
    """Éleveurs désignés par trigramme (sans tenir compte de la casse) ou par id Index Mellifera. Un trigramme
    porté par plusieurs éleveurs (pays différents) est refusé : sinon deux programmes seraient fusionnés."""
    sql = ("SELECT b.id, b.trigram, c.name AS country FROM breeder b LEFT JOIN country c ON c.id = b.country_id "
           "WHERE {} ORDER BY b.id")
    out: dict[int, sqlite3.Row] = {}
    for ref in (r.strip() for r in refs):
        if ref.isdigit():
            rows = conn.execute(sql.format("b.id = ?"), (int(ref),)).fetchall()
        else:
            rows = conn.execute(sql.format("UPPER(b.trigram) = UPPER(?)"), (ref,)).fetchall()
        if not rows:
            raise ValueError(f"éleveur inconnu : « {ref} »")
        if len(rows) > 1:
            listed = ", ".join(f"{r['id']} ({r['country'] or 'pays inconnu'})" for r in rows)
            raise ValueError(f"trigramme « {ref} » porté par plusieurs éleveurs : {listed} ; "
                             f"désignez l'éleveur par son id (--breeder {rows[0]['id']})")
        out[rows[0]["id"]] = rows[0]
    return list(out.values())


def _select(conn: sqlite3.Connection, breeder_ids: list[int], from_year: Optional[int]) -> dict[str, sqlite3.Row]:
    """Reines des éleveurs choisis, nées à partir de from_year."""
    marks = ", ".join("?" for _ in breeder_ids)
    sql = f"SELECT * FROM queen_mating WHERE breeder_id IN ({marks})"
    args: list = list(breeder_ids)
    if from_year is not None:
        sql += " AND birth_year >= ?"
        args.append(from_year)
    return {r["code"]: r for r in conn.execute(sql, args)}


def _ancestry(conn: sqlite3.Connection, selected: dict[str, sqlite3.Row]) -> dict[str, sqlite3.Row]:
    """Sélection + toute son ascendance présente en base (mères et « Drones »), quel que soit l'éleveur."""
    rows = dict(selected)
    frontier = list(selected.values())
    while frontier:
        refs = {c for r in frontier for c in (r["mother_code"], r["drone_mother_code"]) if c and c not in rows}
        frontier = []
        for code in sorted(refs):
            r = conn.execute("SELECT * FROM queen_mating WHERE code = ?", (code,)).fetchone()
            if r is not None:
                rows[code] = r
                frontier.append(r)
    return rows


def _queens(rows: dict[str, sqlite3.Row], selected: set[str], report: ExportReport) -> list[dict]:
    kept = {}
    for code, r in rows.items():
        if _valid_year(r["birth_year"]):
            kept[code] = r
        else:
            report.counts["invalid_birth_year"] += 1
    ids = {code: r["id"] for code, r in kept.items()}
    queens = []
    for code, r in sorted(kept.items(), key=lambda kv: kv[1]["id"]):
        mother = r["mother_code"]
        drone = r["drone_mother_code"]
        if mother == code:
            mother = None
            report.counts["self_reference_cut"] += 1
        if drone == code:
            drone = None
            report.counts["self_reference_cut"] += 1
        q: dict = {"id": r["id"], "born": r["birth_year"], "variety": r["race_description"]}
        if mother in ids and kept[mother]["birth_year"] > r["birth_year"]:
            mother = None                                   # mère plus jeune que sa fille : saisie erronée
            report.counts["younger_parent_cut"] += 1
        if mother in ids:
            q["mother"] = ids[mother]
        elif mother:
            q["mother_ref_ext"] = mother
        mtype = r["mating_type"]
        if drone and mtype == "unknown":
            mtype = "insemination"
            report.counts["mating_type_inferred"] += 1
        if drone in ids and kept[drone]["birth_year"] > r["birth_year"]:
            drone = None
            report.counts["younger_parent_cut"] += 1
        mating: dict = {"type": mtype, "year": r["birth_year"]}
        if drone in ids:
            mating["drone_queens_mother"] = ids[drone]
        elif drone:
            mating["external_ref"] = drone
        if r["n_drone_queens"]:
            mating["n_drone_queens"] = r["n_drone_queens"]
        if r["single_drone"]:
            mating["n_drones"] = 1
        q["mating"] = mating
        if code not in selected:
            q["available_as_mother"] = q["available_as_drone_source"] = False
        queens.append(q)
    _cut_cycles(queens, report)
    return queens


def _cut_cycles(queens: list[dict], report: ExportReport) -> None:
    """Coupe, tant qu'il en reste, un lien de la reine d'id maximal prise dans un cycle (mère d'abord)."""
    by_id = {q["id"]: q for q in queens}
    while True:
        cycle = find_cycle([QueenIn.model_validate(q) for q in queens])
        if not cycle:
            return
        members = set(cycle)
        for qid in sorted(cycle, reverse=True):
            q = by_id[qid]
            if q.get("mother") in members:
                q.pop("mother")
            elif q["mating"].get("drone_queens_mother") in members:
                q["mating"].pop("drone_queens_mother")
            else:
                continue
            report.counts["cycle_link_cut"] += 1
            break
        else:
            raise RuntimeError(f"cycle sans lien à couper : {sorted(cycle)}")


def _records(conn: sqlite3.Connection, codes: dict[int, str], breeders: dict[int, Optional[int]],
             sources: Iterable[str], report: ExportReport) -> tuple[list[dict], list[dict]]:
    marks = ", ".join("?" for _ in sources)
    rows = conn.execute(
        "SELECT e.*, t.id AS trait_id, t.name AS trait, t.kind, t.scale_min, t.scale_max, t.direction "
        "FROM evaluation e LEFT JOIN criterion c ON c.mode = e.mode AND c.code = e.criterion_code "
        f"LEFT JOIN trait t ON t.id = c.trait_id WHERE e.queen_code IS NOT NULL AND e.source IN ({marks}) ORDER BY e.id",
        list(sources)).fetchall()
    id_of = {code: qid for qid, code in codes.items()}
    cells: dict[tuple, list[float]] = defaultdict(list)
    traits: dict[int, sqlite3.Row] = {}
    for e in rows:
        qid = id_of.get(e["queen_code"])
        if qid is None:
            continue
        if e["trait_id"] is None:
            report.counts["unmapped_criterion"] += 1
            continue
        if e["kind"] == "harvest":
            report.counts["harvest"] += 1
            continue
        if e["value"] is None:
            report.counts["non_numeric"] += 1
            continue
        if not e["eval_date"]:
            report.counts["no_date"] += 1
            continue
        lo, hi = (e["scale_min"], e["scale_max"]) if e["kind"] == "note" else (0.0, 1.0)
        value = e["value"] / 100 if e["kind"] == "percentage" else e["value"]
        if not lo <= value <= hi:
            report.counts["out_of_scale"] += 1
            continue
        apiary = e["tester_breeder_id"] if e["source"] == "season_table" else breeders.get(qid)
        if apiary is None:
            report.counts["no_apiary"] += 1
            continue
        traits[e["trait_id"]] = e
        cells[(qid, e["trait_id"], apiary, e["eval_date"], e["source"], e["tester_breeder_id"])].append(value)
    season_years = {(k[0], k[1], k[3][:4]) for k in cells if k[4] == "season_table"}
    records = []
    for (qid, tid, apiary, date, source, evaluator), values in sorted(cells.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][3], kv[0][4], kv[0][2])):
        if source == "radar" and (qid, tid, date[:4]) in season_years:
            report.counts["radar_covered"] += 1
            continue
        if len(values) > 1:
            report.counts["merged_duplicates"] += len(values) - 1
        t = traits[tid]
        lo, hi = (t["scale_min"], t["scale_max"]) if t["kind"] == "note" else (0.0, 1.0)
        records.append({"queen": qid, "criterion": tid, "value": round(sum(values) / len(values), 4),
                        "scale_min": lo, "scale_max": hi, "date": date, "evaluator": evaluator, "apiary": apiary,
                        "hive_class": "hive"})
    criteria = [{"id": tid, "name": t["trait"], "kind": t["kind"], "direction": t["direction"]}
                for tid, t in sorted(traits.items())]
    return criteria, records


def export_blup(conn: sqlite3.Connection, breeders: Iterable[str], *, from_year: Optional[int] = None,
                sources: Iterable[str] = ("season_table",), profile: str = "standard") -> tuple[BlupV2Input, ExportReport]:
    """`breeders` : trigrammes ou id Index Mellifera. ValueError : éleveur inconnu ou ambigu, source inconnue,
    aucune reine sélectionnée."""
    sources = tuple(sources)
    bad = [s for s in sources if s not in SOURCES]
    if bad or not sources:
        raise ValueError(f"sources inconnues : {bad} (attendu : {', '.join(SOURCES)})")
    chosen = _breeders(conn, breeders)
    report = ExportReport(breeders=[f"{b['trigram']} ({b['id']}, {b['country'] or 'pays inconnu'})" for b in chosen])
    selected = _select(conn, [b["id"] for b in chosen], from_year)
    if not selected:
        raise ValueError(f"aucune reine en base pour {', '.join(report.breeders)}"
                         + (f" née à partir de {from_year}" if from_year is not None else ""))
    rows = _ancestry(conn, selected)
    queens = _queens(rows, set(selected), report)
    codes = {r["id"]: code for code, r in rows.items()}
    breeders = {r["id"]: r["breeder_id"] for r in rows.values()}
    kept = {q["id"] for q in queens}
    criteria, records = _records(conn, {i: c for i, c in codes.items() if i in kept}, breeders, sources, report)
    report.n_queens = len(queens)
    report.n_selected = sum(1 for q in queens if q.get("available_as_mother", True))
    report.n_records = len(records)
    inp = BlupV2Input.model_validate({"exploitation": 0, "profile": profile, "criteria": criteria,
                                      "queens": queens, "records": records})
    return inp, report
