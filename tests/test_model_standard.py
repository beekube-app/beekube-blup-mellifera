import numpy as np
import pytest

from blup2.kinship import inverse_relationship, relationship_matrix
from blup2.model_standard import FitParams, fit_standard
from blup2.pedigree import build_pedigree
from blup2.records import prepare_criterion
from blup2.schemas import BlupV2Input
from tests.builders import criterion, mating, queen, record
from tests.simdata import simulate_scenario


def _fit(inp, **kw):
    v = BlupV2Input.model_validate(inp) if isinstance(inp, dict) else inp
    ped, _ = build_pedigree(v.queens, v.options.default_n_drones)
    cdata = prepare_criterion(v, v.criteria[0].id)
    return fit_standard(cdata, ped, inverse_relationship(ped), FitParams(**kw)), ped, cdata


def _family_input(h2=0.3, rep=0.5):
    # mère 1 (non notée) et ses 4 filles notées ; 6 et 7 : reines de base notées (sans elles, la
    # moyenne de la famille serait absorbée par l'intercept et la mère resterait sans information) ;
    # reine 9 sans lien ni note
    qs = [queen(1, 2020, mating=mating("station", station=1, n_rpm=8)), queen(9, 2020)]
    qs += [queen(k, 2022, mother=1) for k in range(2, 6)] + [queen(6, 2022), queen(7, 2022)]
    vals = {2: 0.8, 3: 0.6, 4: 0.7, 5: 0.9, 6: 0.3, 7: 0.5}
    recs = [record(k, v) for k, v in vals.items()] + [record(k, v - 0.1, date="2024-06-01") for k, v in vals.items()]
    return BlupV2Input(exploitation=1, criteria=[criterion(h2=h2, repeatability=rep)], queens=qs, records=recs)


def test_ebv_matches_dense_henderson_solution():
    inp = _family_input()
    fit, ped, cdata = _fit(inp)
    c = fit.parameters
    s2A = c.sigma2_a / np.mean([relationship_matrix(ped, [ped.w_node(q)])[0, 0] for q in cdata.recorded_queens])
    # modèle dense : y = 1 μ + Z_W a + Z_pe pe + e ; V = Z_W A Z_Wᵀ σ²_A + Z_pe Z_peᵀ σ²_pe + I σ²_e
    A = relationship_matrix(ped)
    Zw = np.zeros((cdata.n_obs, ped.n_nodes))
    Zw[np.arange(cdata.n_obs), [ped.w_node(q) for q in cdata.queen_ids]] = 1
    rec = {q: k for k, q in enumerate(cdata.recorded_queens)}
    Zp = np.zeros((cdata.n_obs, len(rec)))
    Zp[np.arange(cdata.n_obs), [rec[q] for q in cdata.queen_ids]] = 1
    V = Zw @ A @ Zw.T * s2A + Zp @ Zp.T * c.sigma2_pe + np.eye(cdata.n_obs) * c.sigma2_e
    X = cdata.X.toarray()
    Vi = np.linalg.inv(V)
    b = np.linalg.solve(X.T @ Vi @ X, X.T @ Vi @ cdata.y)
    a = A @ Zw.T @ Vi @ (cdata.y - X @ b) * s2A
    for q in ped.queen_ids:
        assert fit.ebv[q] == pytest.approx(a[ped.w_node(q)], abs=1e-9)
        assert fit.ebv_queen[q] == pytest.approx(a[ped.q_node(q)], abs=1e-9)


def test_input_parameters_are_used():
    fit, _, _ = _fit(_family_input(0.3, 0.5))
    p = fit.parameters
    assert p.source == "input" and p.h2 == pytest.approx(0.3) and p.repeatability == pytest.approx(0.5)
    assert p.reliability_method == "exact" and p.profile_used == "standard"
    assert p.n_records == 12 and p.n_queens == 6 and p.n_contemporary_groups == 1


def test_statuses():
    fit, _, _ = _fit(_family_input())
    assert [fit.status[q] for q in (2, 3, 4, 5, 6, 7)] == ["evaluated"] * 6
    assert fit.status[1] == "pedigree"          # mère de reines notées
    assert fit.status[9] == "none"              # sans lien ni note
    assert fit.reliability[9] == pytest.approx(0.0, abs=1e-9) and fit.ebv[9] == pytest.approx(0.0, abs=1e-12)


def test_better_daughters_rank_higher():
    fit, _, _ = _fit(_family_input())
    assert fit.ebv[5] > fit.ebv[2] > fit.ebv[4] > fit.ebv[3]     # moyennes 0,85 > 0,75 > 0,65 > 0,55
    assert fit.ebv[1] > 0.0 > fit.ebv[6]          # famille de 1 au-dessus des reines de base 6 et 7
    assert 0.0 < fit.reliability[3] < 1.0


def test_few_records_warning_and_prior_source():
    inp = _family_input()
    inp.criteria[0].h2 = None
    inp.criteria[0].repeatability = None
    fit, _, _ = _fit(inp)                          # auto : 6 reines < 100 -> a priori
    assert fit.parameters.source == "prior"
    assert fit.parameters.h2 == pytest.approx(0.25) and fit.parameters.repeatability == pytest.approx(0.30)
    assert [w.code for w in fit.warnings] == ["few_records"]


def test_criterion_without_records_gives_status_none():
    v = BlupV2Input(exploitation=1, criteria=[criterion()], queens=[queen(1), queen(2)], records=[])
    fit, _, _ = _fit(v)
    assert set(fit.status.values()) == {"none"} and fit.ebv == {}
    assert fit.parameters.reliability_method == "none" and fit.warnings[0].code == "few_records"


def test_reml_is_used_when_always():
    inp, _, _ = simulate_scenario(1)
    fit, _, _ = _fit(inp, estimate_variances="always")
    assert fit.parameters.source in ("reml", "prior")
    if fit.parameters.source == "prior":
        assert any(w.code in ("reml_not_converged", "low_heritability") for w in fit.warnings)
    else:
        assert 0.0 < fit.parameters.h2 < 1.0


def test_null_heritability_gives_low_heritability_and_prior():
    # notes sans composante génétique (bruit pur), 193 reines notées : REML en auto, la variance
    # génétique finit au bord inférieur -> a priori et low_heritability, pas reml_not_converged
    inp, _, _ = simulate_scenario(1)
    rng = np.random.default_rng(1)
    for r in inp["records"]:
        r["value"] = float(5.0 + rng.standard_normal())
    fit, _, _ = _fit(inp)
    assert fit.parameters.source == "prior" and fit.parameters.repeatability == pytest.approx(0.30)
    assert [w.code for w in fit.warnings] == ["low_heritability"]
    assert fit.warnings[0].criterion == 1 and "Héritabilité estimée nulle" in fit.warnings[0].message


def test_constant_notes_give_no_variation_without_reml():
    # notes toutes identiques (σ²_P = 0 après effets fixes) : sans garde, le REML part vers des
    # variances de 1e-300 et des index de ±1e139 ; attendu no_variation, a priori, EBV nulles
    inp, _, _ = simulate_scenario(1)
    for r in inp["records"]:
        r["value"] = 5.0
    fit, _, _ = _fit(inp)
    assert fit.parameters.source == "prior"
    assert [(w.code, w.criterion) for w in fit.warnings] == [("no_variation", 1)]
    assert fit.sigma_index > 0.0
    assert max(abs(e) for e in fit.ebv.values()) < 1e-6 * fit.sigma_index    # index 100 pour toutes


def test_reliability_is_one_minus_pev_over_prior_variance_on_simulation():
    inp, truth, recorded = simulate_scenario(4, true_parameters=True)
    fit, _, _ = _fit(inp)
    ebv = np.array([fit.ebv[q] for q in truth])
    tv = np.array([truth[q] for q in truth])
    assert np.corrcoef(ebv, tv)[0, 1] > 0.3
    rec = [fit.reliability[q] for q in recorded]
    unrec = [fit.reliability[q] for q in truth if q not in recorded]
    assert np.mean(rec) > np.mean(unrec) > 0.0


@pytest.mark.parametrize("seed", [1, 2])
def test_approximate_reliabilities_track_exact_ones(seed):
    inp, _, _ = simulate_scenario(seed, true_parameters=True)
    exact, _, _ = _fit(inp, estimate_variances="never")
    approx, _, _ = _fit(inp, estimate_variances="never", exact_limit=0)
    assert exact.parameters.reliability_method == "exact" and approx.parameters.reliability_method == "approx"
    q = list(exact.reliability)
    r_ex = np.array([exact.reliability[k] for k in q])
    r_ap = np.array([approx.reliability[k] for k in q])
    # seuils issus du prototype : corrélation 0,92-0,99, surestimation moyenne 0,06-0,12
    assert np.corrcoef(r_ex, r_ap)[0, 1] >= 0.85
    assert abs(float(np.mean(r_ap - r_ex))) <= 0.15
    assert np.allclose([exact.ebv[k] for k in q], [approx.ebv[k] for k in q])   # mêmes solutions
    assert any(w.code == "reliability_approx" for w in approx.warnings)


def test_repeatability_alone_caps_h2_at_the_repeatability():
    inp = _family_input(None, 0.1)
    fit, _, _ = _fit(inp)
    p = fit.parameters
    assert p.source == "input" and p.h2 <= 0.1 and p.repeatability == pytest.approx(0.1, abs=1e-6)
    assert [w.code for w in fit.warnings] == ["few_records"]      # seul avertissement : peu de reines notées


def test_scale_mixed_warning_when_two_scales_for_one_criterion():
    inp = _family_input()
    inp.records[0].scale_min, inp.records[0].scale_max = 0.0, 4.0
    inp.records[0].value = 3.0
    fit, _, _ = _fit(inp)
    assert any(w.code == "scale_mixed" and w.criterion == 1 for w in fit.warnings)
