"""Données simulées pour les tests de calibration et de performance."""
from __future__ import annotations

import os
import sys

import numpy as np

_ANALYSES = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "docs", "analyses")
if _ANALYSES not in sys.path:
    sys.path.insert(0, _ANALYSES)

SCALE_MIN, SCALE_MAX = -20.0, 30.0  # y simulé ~ 5 ± 2 ; normalisé = (y + 20) / 50


def scenario_pedigree(rng: np.random.Generator, S: int = 8, D: int = 12):
    """Pedigree de simulate.gen_pedigree : (pedigree de référence BeePedigree, individus,
    reines au format v2 (dicts), correspondance id simulé -> id de reine v2)."""
    from simulate import build_bb, gen_pedigree

    inds, mate = gen_pedigree(rng)
    P = build_bb(inds, S, D)
    queens = [x for x in inds if x["kind"] == "Q"]
    sires = {x["id"]: x for x in inds if x["kind"] == "S"}
    qid = {x["id"]: 1000 + k for k, x in enumerate(queens)}
    station = {sid: 100 + k for k, sid in enumerate(sires)}

    def mating_of(q):
        s = sires[mate[q]]
        return {"type": "station", "station": station[s["id"]], "year": 2000 + s["gen"],
                "drone_queens_mother": qid[s["dam"]] if s["dam"] is not None else None,
                "n_drone_queens": S, "n_drones": D}

    q_in = [{"id": qid[x["id"]], "born": 2015 + x["gen"],
             "mother": qid[x["dam"]] if x["dam"] is not None else None,
             "mating": mating_of(x["id"])} for x in queens]
    return P, inds, q_in, qid


def simulate_scenario(seed: int, s2e: float = 1.5, s2pe: float = 0.5, s2ap: float = 0.5, S: int = 8,
                      D: int = 12, p_notes: float = 0.7, n_ap: int = 4, true_parameters: bool = False):
    """Scénario de docs/analyses/simulate.py (4 générations, sires = groupes de 8 RPM, D = 12).

    Vraies valeurs tirées avec la matrice A de référence (bee_kinship.BeePedigree, σ²_A = 1).
    true_parameters : transmet h² et répétabilité vrais des colonies notées
    (h² = ā_WW / (ā_WW + σ²_pe + σ²_e), répétabilité = (ā_WW + σ²_pe) / (ā_WW + σ²_pe + σ²_e)) ;
    sinon h2 / repeatability = None (paramètres a priori ou REML).
    Retourne (entrée v2 en dict, vraie valeur du nœud W de chaque reine en unités normalisées,
    ensemble des reines notées).
    """
    rng = np.random.default_rng(seed)
    P, inds, q_in, qid = scenario_pedigree(rng, S, D)
    A = P.build()
    n = len(inds)
    L = np.linalg.cholesky(A + 1e-12 * np.eye(n))
    a_true = L @ rng.standard_normal(n)
    idx = {x["id"]: i for i, x in enumerate(inds)}
    queens = [x for x in inds if x["kind"] == "Q"]
    apiary_of = {x["id"]: int(rng.integers(n_ap)) for x in queens}
    ap_eff = rng.standard_normal(n_ap) * np.sqrt(s2ap)
    records, truth, recorded = [], {}, set()
    for x in queens:
        q = x["id"]
        truth[qid[q]] = a_true[idx["W_" + q]] / (SCALE_MAX - SCALE_MIN)
        if rng.random() >= p_notes:
            continue
        recorded.add(qid[q])
        pe = rng.standard_normal() * np.sqrt(s2pe)
        for k in range(int(rng.integers(2, 6))):
            y = 5.0 + a_true[idx["W_" + q]] + pe + ap_eff[apiary_of[q]] + rng.standard_normal() * np.sqrt(s2e)
            records.append({"queen": qid[q], "criterion": 1, "value": float(y), "scale_min": SCALE_MIN,
                            "scale_max": SCALE_MAX, "date": f"{2016 + x['gen']}-0{4 + k}-10",
                            "evaluator": 1, "apiary": apiary_of[q] + 1})
    h2 = repeatability = None
    if true_parameters:
        aww = float(np.mean([A[idx["W_" + x["id"]], idx["W_" + x["id"]]] for x in queens if qid[x["id"]] in recorded]))
        s2p = aww + s2pe + s2e
        h2, repeatability = aww / s2p, (aww + s2pe) / s2p
    crit = {"id": 1, "name": "sim", "direction": 1, "kind": "note", "weight": 1,
            "h2": h2, "repeatability": repeatability}
    inp = {"exploitation": 1, "profile": "standard", "options": {"n_suggestions": 0},
           "criteria": [crit], "queens": q_in, "records": records}
    return inp, truth, recorded


def big_population(seed: int = 7, n_queens: int = 3000, records_per_queen: int = 10, n_gen: int = 6,
                   p_controlled: float = 0.5, n_apiaries: int = 60, n_evaluators: int = 20) -> dict:
    """Population de taille réelle pour le test de performance (valeurs aléatoires)."""
    rng = np.random.default_rng(seed)
    per_gen = n_queens // n_gen
    queens, prev = [], []
    nid = 1
    for g in range(n_gen):
        cur = []
        for _ in range(per_gen):
            q = {"id": nid, "born": 2015 + g, "mother": int(rng.choice(prev)) if prev else None}
            if prev and rng.random() < p_controlled:
                q["mating"] = {"type": "insemination", "drone_queens": [int(rng.choice(prev))],
                               "n_drones": 8}
            else:
                q["mating"] = {"type": "natural"}
            queens.append(q)
            cur.append(nid)
            nid += 1
        prev = cur
    records = []
    for q in queens:
        ap = int(rng.integers(n_apiaries)) + 1
        ev = int(rng.integers(n_evaluators)) + 1
        for k in range(records_per_queen):
            records.append({"queen": q["id"], "criterion": 1, "value": int(rng.integers(1, 6)),
                            "scale_min": 1, "scale_max": 5, "date": f"{q['born'] + 1}-{4 + k % 6:02d}-{1 + k:02d}",
                            "evaluator": ev, "apiary": ap, "hive_class": "hive", "visit": None})
    return {"exploitation": 1, "profile": "standard", "options": {},
            "criteria": [{"id": 1, "name": "perf", "direction": 1, "kind": "note", "weight": 1}],
            "queens": queens, "records": records}
