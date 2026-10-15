"""Validation statistique du volet 3 Index Mellifera (spec volet 3 §7) et coût. Tests lents :
    .venv/bin/python -m pytest -m slow
"""
import time

import numpy as np
import pytest
from scipy.stats import spearmanr

from blup2.engine import run_blup
from blup2.schemas import BlupV2Input
from tests.simdata import external_scenario, write_crawler_db

N_REPLICATES = 20


def _indexes(monkeypatch, inp, db, trait):
    monkeypatch.setenv("MELLIFERA_DB", db)
    crit = {**inp["criteria"][0], **({"trait": "gentleness"} if trait else {})}
    out = run_blup(BlupV2Input.model_validate({**inp, "criteria": [crit]}))
    assert out.parameters["1"].reliability_method == "exact"   # le gain mesuré n'est pas un effet de l'approximation
    return {q.id: q.indexes["1"] for q in out.queens}


def _replicates(tmp_path, monkeypatch, top_half=False):
    """Par répétition : gain de fiabilité moyen des reines liées, corrélations de Pearson et de Spearman avec
    la vraie valeur (avec, sans notes externes), écart maximal des EBV des reines sans lien."""
    res = []
    for seed in range(1, N_REPLICATES + 1):
        inp, rows, notes, truth, free = external_scenario(seed, top_half=top_half)
        db = write_crawler_db(tmp_path / f"m{seed}.sqlite", rows, notes)
        a = _indexes(monkeypatch, inp, db, trait=True)
        b = _indexes(monkeypatch, inp, db, trait=False)            # ascendance du volet 2, sans les notes
        ids = list(truth)
        t = np.array([truth[q] for q in ids])
        ea, eb = (np.array([x[q].ebv for q in ids]) for x in (a, b))
        gain = float(np.mean([a[q].reliability - b[q].reliability for q in ids]))
        res.append((gain, np.corrcoef(t, ea)[0, 1], np.corrcoef(t, eb)[0, 1], spearmanr(t, ea)[0],
                    spearmanr(t, eb)[0], max(abs(a[q].ebv - b[q].ebv) for q in free)))
    return np.array(res)


@pytest.mark.slow
def test_notes_of_external_lineages_improve_their_daughters(tmp_path, monkeypatch):
    """(a) et (b) de la spec volet 3 §7, variances fixées."""
    r = _replicates(tmp_path, monkeypatch)
    # prototype (20 répétitions) : gain de fiabilité 0,068 (min 0,047) ; corrélation 0,472 -> 0,517
    assert (r[:, 0] > 0).all() and r[:, 0].mean() >= 0.05
    assert r[:, 1].mean() - r[:, 2].mean() >= 0.02
    assert r[:, 5].max() <= 1e-9                    # reines sans lien ni groupe commun : EBV inchangées


@pytest.mark.slow
def test_selected_external_notes_do_not_degrade_the_ranking(tmp_path, monkeypatch):
    """(c) de la spec volet 3 §7 : la base ne garde que la meilleure moitié des reines externes notées."""
    r = _replicates(tmp_path, monkeypatch, top_half=True)
    # prototype : gain de fiabilité 0,033 ; Spearman 0,452 sans notes, 0,472 avec
    assert (r[:, 0] > 0).all()
    assert r[:, 3].mean() >= r[:, 4].mean()


def _big_crawler_db(path, rng, n_queens=9000, years=range(2010, 2025), n_testers=40):
    """Base synthétique de taille réelle : une génération par an, mère et mère des sœurs de l'année d'avant,
    80 % des reines notées l'année suivant leur naissance par un testeur tiré au hasard."""
    per_year = n_queens // len(years)
    codes = {y: [f"FR-B{k % 50:02d}-B-{k + 1:04d}-{y}" for k in range(per_year)] for y in years}
    rows, notes = [], []
    for y in years:
        prev, prev2 = codes.get(y - 1), codes.get(y - 2)
        for c in codes[y]:
            rows.append((c, y, str(rng.choice(prev)) if prev else None, str(rng.choice(prev2)) if prev2 else None, 8))
            if rng.random() < 0.8 and y + 1 <= 2025:
                notes.append((c, "season_table", "evjp", "mildness", float(rng.integers(0, 7)), f"{y + 1}-07-01",
                              300 + int(rng.integers(n_testers))))
    return write_crawler_db(path, rows, notes), codes


@pytest.mark.slow
def test_cost_on_a_crawler_base_of_9000_queens(tmp_path, monkeypatch):
    """Spec volet 3 §9 : 300 reines de l'exploitation citant 100 lignées d'une base de 9 000 reines."""
    rng = np.random.default_rng(3)
    db, codes = _big_crawler_db(tmp_path / "big.sqlite", rng)
    mothers = rng.choice(codes[2023], 100, replace=False)
    queens, records = [], []
    for k in range(300):
        queens.append({"id": k + 1, "born": 2024, "mother_ref_ext": str(mothers[k % 100]),
                       "mating": {"type": "unknown", "external_ref": str(rng.choice(codes[2022]))}})
        ap = int(rng.integers(5)) + 1
        records += [{"queen": k + 1, "criterion": 1, "value": float(rng.integers(1, 6)), "scale_min": 1,
                     "scale_max": 5, "date": f"2025-0{4 + j}-10", "evaluator": 1, "apiary": ap} for j in range(3)]
    crit = {"id": 1, "name": "perf", "direction": 1, "kind": "note", "weight": 1, "h2": 0.3, "repeatability": 0.4}
    inp = {"exploitation": 1, "profile": "standard", "options": {"n_suggestions": 0}, "queens": queens,
           "records": records}
    monkeypatch.setenv("MELLIFERA_DB", db)
    times, outs = [], []
    for c in (crit, crit, {**crit, "trait": "gentleness"}):          # premier calcul : préchauffage
        t0 = time.perf_counter()
        outs.append(run_blup(BlupV2Input.model_validate({**inp, "criteria": [c]})))
        times.append(time.perf_counter() - t0)
    added = [w.message for w in outs[2].warnings if w.code == "external_records_added"]
    print(f"\nsans notes {times[1]:.2f} s, avec notes {times[2]:.2f} s ; {added}")
    assert added and "plafond de 1000 contemporaines atteint" in added[0]
    assert [o.parameters["1"].reliability_method for o in outs] == ["exact"] * 3    # spec volet 3 §5
    assert times[2] < 30.0
