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


def _slope(t, e):
    return float(np.cov(t, e)[0, 1] / np.var(e, ddof=1))        # régression vraie ~ estimée


def _breeder_replicate(seed):
    from blup2.kinship import inbreeding, inverse_relationship
    from blup2.model_breeder import fit_breeder
    from blup2.model_standard import FitParams
    from blup2.pedigree import build_pedigree
    from blup2.records import prepare_criterion
    from tests.simdata import simulate_breeder_scenario

    inp, truth, _, _ = simulate_breeder_scenario(seed)
    v = BlupV2Input.model_validate(inp)
    ped, _ = build_pedigree(v.queens, v.options.default_n_drones)
    F = inbreeding(ped)
    cdata = prepare_criterion(v, 1, f_workers={q: float(F[ped.w_node(q)]) for q in ped.queen_ids})
    fit = fit_breeder(cdata, ped, inverse_relationship(ped), FitParams())
    assert fit.parameters.profile_used == "breeder" and fit.parameters.reliability_method == "exact"
    assert fit.parameters.source == "input"                      # paramètres vrais fournis, pas de REML
    ids = list(truth["total"])

    def arr(d):
        return np.array([d[q] for q in ids])

    t, e = arr(truth["total"]), arr(fit.ebv)
    return (float(np.corrcoef(t, e)[0, 1]), _slope(t, e), float(np.mean(arr(fit.reliability))),
            _slope(arr(truth["worker"]), arr(fit.ebv_worker)),
            _slope(arr(truth["queen_total"]), arr(fit.ebv_queen)))


@pytest.mark.slow
def test_breeder_calibration_on_simulated_populations():
    res = np.array([_breeder_replicate(seed) for seed in range(1, N_REPLICATES + 1)])
    corr, slope, mean_rel, slope_worker, slope_queen = res.mean(axis=0)
    # prototype (20 répétitions, 11 s) : corrélation 0,546, pente 1,082, r² moyen 0,252, corr² 0,307,
    # pente reine (nœud Q) 1,050, pente ouvrières 1,120 (60 répétitions : 1,025 / 1,016 / 1,072)
    assert corr >= 0.45
    assert 0.8 <= slope <= 1.2                                   # valeur totale publiée sans biais d'échelle
    assert 0.8 <= slope_queen <= 1.2                             # valeur reine utilisée aux accouplements
    assert abs(mean_rel - float(np.mean(res[:, 0] ** 2))) <= 0.08
    assert 0.7 <= slope_worker <= 1.4                            # composante ouvrières : plus bruitée


@pytest.mark.slow
def test_breeder_full_reml_when_always():
    # REML bivarié (ρ et rG libres) : peu identifiable sur ~1000 colonies ; on vérifie seulement
    # des estimations admissibles et un repli propre (prototype : 4/4 convergent, ρ 0,10–1,29)
    from tests.simdata import simulate_breeder_scenario

    kw = dict(n_gen=7, n_base=150, n_sires=10, n_dams=50)
    sources = []
    for seed in (1, 2):
        inp, _, _, _ = simulate_breeder_scenario(seed, true_parameters=False, **kw)
        inp["options"]["estimate_variances"] = "always"
        out = run_blup(BlupV2Input.model_validate(inp))
        p = out.parameters["1"]
        sources.append(p.source)
        if p.source == "reml":
            assert abs(p.rg) <= 0.95 and p.sigma2_aq > 0 and 0.0 < p.h2 < 1.0
        else:
            assert any(w.code in ("reml_not_converged", "low_heritability") for w in out.warnings)
    assert "reml" in sources


@pytest.mark.slow
def test_breeder_performance_3000_queens_30000_records():
    d = big_population(p_controlled=0.7)                         # ≥ 50 % d'accouplements connus
    d["profile"] = "breeder"
    inp = BlupV2Input.model_validate(d)
    t0 = time.perf_counter()
    out = run_blup(inp)
    elapsed = time.perf_counter() - t0
    assert elapsed < 60.0                                        # prototype : 30 s (REML auto, 15 380 éq.)
    p = out.parameters["1"]
    assert p.profile_used == "breeder" and p.reliability_method == "approx"
    assert p.source == "reml"                                    # sinon un REML sauté passerait le seuil de 60 s
    assert len(out.queens) == 3000
