"""Calibration par simulation (scénario docs/analyses/simulate.py) et performance. Tests lents :
    .venv/bin/python -m pytest -m slow
"""
import time

import numpy as np
import pytest

from blup2.engine import run_blup
from blup2.schemas import BlupV2Input
from tests.simdata import big_population, simulate_scenario

N_REPLICATES = 20


def _replicate(seed):
    inp, truth, recorded = simulate_scenario(seed, true_parameters=True)
    out = run_blup(BlupV2Input.model_validate(inp))
    ids = list(truth)
    ebv = {q.id: q.indexes["1"].ebv for q in out.queens}
    rel = {q.id: q.indexes["1"].reliability for q in out.queens}
    t = np.array([truth[q] for q in ids])
    e = np.array([ebv[q] if ebv[q] is not None else 0.0 for q in ids])
    r = np.array([rel[q] if rel[q] is not None else 0.0 for q in ids])
    corr = float(np.corrcoef(t, e)[0, 1])
    slope = float(np.cov(t, e)[0, 1] / np.var(e, ddof=1))        # régression vraie ~ estimée
    unrec = [k for k, q in enumerate(ids) if q not in recorded]
    corr_unrec = float(np.corrcoef(t[unrec], e[unrec])[0, 1])
    return corr, slope, float(np.mean(r)), corr_unrec


@pytest.mark.slow
def test_calibration_on_simulated_populations():
    res = np.array([_replicate(seed) for seed in range(1, N_REPLICATES + 1)])
    corr, slope, mean_rel, corr_unrec = res.mean(axis=0)
    # moteur v2, CG fixe rucher × année : 0,571 ; prototype rucher aléatoire : 0,60–0,61 (pente 1,018, r² moyen 0,352, corr² 0,381)
    assert corr >= 0.55
    assert 0.8 <= slope <= 1.2                                  # EBV sans biais d'échelle
    assert abs(mean_rel - float(np.mean(res[:, 0] ** 2))) <= 0.08  # fiabilité annoncée ≈ corr²
    assert corr_unrec > 0.3                                     # reines non notées : information pedigree


@pytest.mark.slow
def test_reml_estimates_heritability_on_simulated_populations():
    h2 = []
    for seed in range(1, 11):
        inp, _, _ = simulate_scenario(seed)
        inp["options"]["estimate_variances"] = "always"
        p = run_blup(BlupV2Input.model_validate(inp)).parameters["1"]
        if p.source == "reml":
            h2.append(p.h2)
    # vrai h² des colonies ≈ 0,17 (ā_WW ≈ 0,42 ; σ²_pe = 0,5 ; σ²_e = 1,5)
    # prototype : 0,144, 10 convergences sur 10
    assert len(h2) >= 7 and 0.08 <= float(np.mean(h2)) <= 0.30


@pytest.mark.slow
def test_performance_3000_queens_30000_records():
    inp = BlupV2Input.model_validate(big_population())
    assert len(inp.queens) == 3000 and len(inp.records) == 30000
    t0 = time.perf_counter()
    out = run_blup(inp)
    elapsed = time.perf_counter() - t0
    assert elapsed < 60.0                                       # prototype : 17,7 s
    assert out.parameters["1"].reliability_method == "approx"   # 6000 nœuds + 3000 pe > 6000 éq.
    assert len(out.queens) == 3000
