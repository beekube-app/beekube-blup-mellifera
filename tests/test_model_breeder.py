"""Profil breeder : modèle à deux effets corrélés (ouvrières sur W, reine sur Q) [BB éq. 1–5]."""
import numpy as np
import pytest
import scipy.optimize

from blup2 import model_breeder
from blup2.kinship import inverse_relationship, relationship_matrix
from blup2.model_breeder import fit_breeder, genetic_inverse
from blup2.model_standard import FitParams, fit_standard
from blup2.pedigree import build_pedigree
from blup2.records import prepare_criterion
from blup2.schemas import BlupV2Input
from tests.builders import criterion, mating, queen, record
from tests.simdata import simulate_breeder_scenario


def _fit(inp, fn=fit_breeder, **kw):
    v = BlupV2Input.model_validate(inp) if isinstance(inp, dict) else inp
    ped, _ = build_pedigree(v.queens, v.options.default_n_drones)
    cdata = prepare_criterion(v, v.criteria[0].id)
    return fn(cdata, ped, inverse_relationship(ped), FitParams(**kw)), ped, cdata


def _station(s):
    return mating("station", station=s, n_rpm=8)


def _two_colonies(**crit):
    # deux reines de base fécondées sur deux stations différentes (sires de base distincts, S = 8,
    # D = 12) : colonies non apparentées ; une note chacune, 0,2 et 0,6 ; aucun effet fixe hors moyenne
    c = dict(h2=0.3, repeatability=0.5, queen_ratio=0.5, rg=-0.5)
    c.update(crit)
    return BlupV2Input(exploitation=1, criteria=[criterion(**c)],
                       queens=[queen(1, 2022, mating=_station(1)), queen(2, 2022, mating=_station(2))],
                       records=[record(1, 0.2), record(2, 0.6)])


def _family(**crit):
    # mère 1 et ses 4 filles (sœurs) ; 6 et 7 reines de base ; 9 sans lien ni note ;
    # toutes les reines notées ont un accouplement connu (station)
    c = dict(h2=0.3, repeatability=0.5, queen_ratio=0.5, rg=-0.5)
    c.update(crit)
    qs = [queen(1, 2020, mating=_station(1)), queen(9, 2020)]
    qs += [queen(k, 2022, mother=1, mating=_station(2)) for k in range(2, 6)]
    qs += [queen(6, 2022, mating=_station(3)), queen(7, 2022, mating=_station(3))]
    vals = {2: 0.8, 3: 0.6, 4: 0.7, 5: 0.9, 6: 0.3, 7: 0.5}
    recs = [record(k, v) for k, v in vals.items()] + [record(k, v - 0.1, date="2024-06-01") for k, v in vals.items()]
    return BlupV2Input(exploitation=1, criteria=[criterion(**c)], queens=qs, records=recs)


def test_two_unrelated_colonies_match_the_hand_calculation():
    # a_WW = ¼ + 1/(4S) + 1/(2D) = 0,322917 ; a_QQ = 1 ; a_WQ = ½ ; ρ = 0,5 ; rG = −0,5 ; √ρ = 0,707107
    # σ²_P = (0,2² + 0,2²) / 1 = 0,08 ; var. génétique d'une note / σ²_AW :
    #   k = a_WW + ρ a_QQ + 2 rG √ρ a_WQ = 0,469363  ->  σ²_AW = h² σ²_P / k = 0,051133
    # Deux colonies indépendantes, moyenne inconnue : T̂ = cov(T, y) / σ²_P × (y − ȳ), et
    # cov(T, y) / σ²_AW = a_WW + rG √ρ (a_WQ + a_WW) + ρ a_WQ = 0,281972 (T = a_W + a_Q du nœud W)
    fit, _, _ = _fit(_two_colonies())
    p = fit.parameters
    assert p.profile_used == "breeder" and p.source == "input"
    assert p.sigma2_a == pytest.approx(0.3 * 0.08 / 0.469363, rel=1e-5)
    assert p.sigma2_aq == pytest.approx(0.5 * p.sigma2_a) and p.rg == pytest.approx(-0.5)
    assert p.h2 == pytest.approx(0.3) and p.repeatability == pytest.approx(0.5)
    assert fit.ebv[1] == pytest.approx(0.3 * 0.281972 / 0.469363 * -0.2, rel=1e-5)        # −0,036045
    # composante ouvrières : cov(a_W(W), y) / σ²_AW = a_WW + rG √ρ a_WQ = 0,146140
    assert fit.ebv_worker[1] == pytest.approx(0.3 * 0.146140 / 0.469363 * -0.2, rel=1e-5)  # −0,018682
    # reine (nœud Q) : cov(T_Q, y) / σ²_AW = a_QW + rG √ρ (a_QQ + a_QW) + ρ a_QQ = 0,469670
    assert fit.ebv_queen[1] == pytest.approx(0.3 * 0.469670 / 0.469363 * -0.2, rel=1e-5)   # −0,060039
    # symétrie : notes opposées autour de la moyenne -> valeurs opposées
    assert fit.ebv[2] == pytest.approx(-fit.ebv[1]) and fit.ebv_queen[2] == pytest.approx(-fit.ebv_queen[1])
    # r² = var(T̂) / var(T) = cov² / (2 σ²_P a_WW σ²_T), σ²_T / σ²_AW = 1 + ρ + 2 rG √ρ = 0,792893
    #    = 0,281972² × 0,3 / (2 × 0,469363 × 0,322917 × 0,792893) = 0,099240
    assert fit.reliability[1] == pytest.approx(0.099240, rel=1e-4)
    assert fit.sigma_index == pytest.approx(np.sqrt(0.792893 * p.sigma2_a * 0.322917), rel=1e-5)
    assert fit.sigma_worker == pytest.approx(np.sqrt(p.sigma2_a * 0.322917), rel=1e-5)
    assert fit.sigma_queen == pytest.approx(np.sqrt(0.792893 * p.sigma2_a), rel=1e-5)


def test_matches_dense_henderson_solution():
    inp = _family()
    fit, ped, cdata = _fit(inp)
    p = fit.parameters
    cov = p.rg * np.sqrt(p.sigma2_a * p.sigma2_aq)
    G0 = np.array([[p.sigma2_a, cov], [cov, p.sigma2_aq]])
    A = relationship_matrix(ped)
    n, N = ped.n_nodes, cdata.n_obs
    Zw, Zq = np.zeros((N, n)), np.zeros((N, n))
    Zw[np.arange(N), [ped.w_node(q) for q in cdata.queen_ids]] = 1
    Zq[np.arange(N), [ped.q_node(q) for q in cdata.queen_ids]] = 1
    rec = {q: k for k, q in enumerate(cdata.recorded_queens)}
    Zp = np.zeros((N, len(rec)))
    Zp[np.arange(N), [rec[q] for q in cdata.queen_ids]] = 1
    Z, G = np.hstack([Zw, Zq]), np.kron(G0, A)
    # V = Z (G0 ⊗ A) Zᵀ + Z_pe Z_peᵀ σ²_pe + I σ²_e ; u = G Zᵀ V⁻¹ (y − X b)
    V = Z @ G @ Z.T + Zp @ Zp.T * p.sigma2_pe + np.eye(N) * p.sigma2_e
    X, Vi = cdata.X.toarray(), np.linalg.inv(V)
    b = np.linalg.solve(X.T @ Vi @ X, X.T @ Vi @ cdata.y)
    u = G @ Z.T @ Vi @ (cdata.y - X @ b)
    P = Vi - Vi @ X @ np.linalg.solve(X.T @ Vi @ X, X.T @ Vi)
    pev_all = G - G @ Z.T @ P @ Z @ G
    sigma2_T = G0.sum()
    for q in ped.queen_ids:
        w, qn = ped.w_node(q), ped.q_node(q)
        assert fit.ebv[q] == pytest.approx(u[w] + u[n + w], abs=1e-9)
        assert fit.ebv_worker[q] == pytest.approx(u[w], abs=1e-9)
        assert fit.ebv_queen[q] == pytest.approx(u[qn] + u[n + qn], abs=1e-9)
        # PEV de T = a_W + a_Q du nœud W : bloc 2 × 2 de l'inverse
        pev = pev_all[w, w] + pev_all[n + w, n + w] + 2 * pev_all[w, n + w]
        expected = min(max(1 - pev / (A[w, w] * sigma2_T), 0.0), 1.0)
        assert fit.reliability[q] == pytest.approx(expected, abs=1e-9)


def test_genetic_inverse_is_the_inverse_of_g0_kron_a():
    ped, _ = build_pedigree(_family().queens, 12)
    G0 = np.array([[1.0, -0.35], [-0.35, 0.5]])
    Gi = genetic_inverse(G0, inverse_relationship(ped)).toarray()
    G = np.kron(G0, relationship_matrix(ped))
    assert np.abs(Gi - Gi.T).max() < 1e-12
    assert np.abs(G @ Gi - np.eye(2 * ped.n_nodes)).max() < 1e-8


def test_tends_to_the_standard_model_when_the_queen_variance_vanishes():
    inp = _family(queen_ratio=1e-8, rg=0.0)
    b, _, _ = _fit(inp)
    s, _, _ = _fit(inp, fit_standard)
    for q in s.ebv:
        assert b.ebv[q] == pytest.approx(s.ebv[q], abs=1e-6)
        assert b.ebv_worker[q] == pytest.approx(s.ebv[q], abs=1e-6)
        assert b.ebv_queen[q] == pytest.approx(s.ebv_queen[q], abs=1e-6)
        assert b.reliability[q] == pytest.approx(s.reliability[q], abs=1e-6)
    assert b.status == s.status
    assert b.sigma_index == pytest.approx(s.sigma_index, rel=1e-4)
    assert b.sigma_queen == pytest.approx(s.sigma_queen, rel=1e-4)


def test_result_does_not_depend_on_input_order():
    inp = _family()
    rev = inp.model_copy(update={"queens": list(reversed(inp.queens)), "records": list(reversed(inp.records))})
    a, _, _ = _fit(inp)
    b, _, _ = _fit(rev)
    for q in a.ebv:
        assert b.ebv[q] == pytest.approx(a.ebv[q], abs=1e-12)
        assert b.ebv_queen[q] == pytest.approx(a.ebv_queen[q], abs=1e-12)
        assert b.reliability[q] == pytest.approx(a.reliability[q], abs=1e-12)


def test_statuses_and_ranking():
    fit, _, _ = _fit(_family())
    assert [fit.status[q] for q in (2, 3, 4, 5, 6, 7)] == ["evaluated"] * 6
    assert fit.status[1] == "pedigree" and fit.status[9] == "none"
    assert fit.ebv[5] > fit.ebv[2] > fit.ebv[4] > fit.ebv[3]
    assert fit.ebv_queen[1] > 0.0 > fit.ebv_queen[6]


def test_defaults_and_clamped_correlation():
    inp = _family(h2=None, repeatability=None, queen_ratio=None, rg=None)
    fit, _, _ = _fit(inp)                                       # auto : 6 reines < 100 -> a priori
    p = fit.parameters
    assert p.source == "prior" and p.h2 == pytest.approx(0.25) and p.repeatability == pytest.approx(0.30)
    assert p.sigma2_aq / p.sigma2_a == pytest.approx(0.5) and p.rg == pytest.approx(-0.5)
    fit, _, _ = _fit(_family(rg=-1.0))
    assert fit.parameters.rg == pytest.approx(-0.95)          # |rG| ≤ 0,95 : G0 définie positive
    fit, _, _ = _fit(_family(rg=None))
    assert fit.parameters.source == "input" and fit.parameters.rg == pytest.approx(-0.5)
    # répétabilité seule fournie, inférieure au h² par défaut : h² = min(0,25 ; répétabilité), comme
    # le profil standard (sinon h² > répétabilité) ; σ²_pe tombe au plancher (≈ 1e-6 σ²_P)
    fit, _, _ = _fit(_family(h2=None, repeatability=0.2, queen_ratio=None, rg=None))
    p = fit.parameters
    assert p.source == "input" and p.h2 == pytest.approx(0.2, abs=1e-4)
    assert p.repeatability == pytest.approx(0.2, abs=1e-4) and p.h2 <= p.repeatability


def _known_matings(n_queens, n_known):
    # n_queens reines notées (une note chacune), les n_known premières à accouplement connu (station),
    # les autres en fécondation naturelle
    qs = [queen(k, 2022, mating=_station(1 + k % 5) if k <= n_known else mating("natural"))
          for k in range(1, n_queens + 1)]
    return BlupV2Input(exploitation=1, criteria=[criterion(h2=0.3, repeatability=0.5)], queens=qs,
                       records=[record(k, 0.5 * k / n_queens) for k in range(1, n_queens + 1)])


def test_fallback_to_standard_below_half_known_matings():
    fit, _, _ = _fit(_known_matings(5, 2))                      # 2 / 5 = 40 % < 50 %
    assert fit.parameters.profile_used == "standard" and fit.ebv_worker is None
    assert fit.warnings[0].code == "breeder_fallback" and fit.warnings[0].criterion == 1
    assert "40 %" in fit.warnings[0].message
    fit, _, _ = _fit(_known_matings(5, 3))                      # 3 / 5 = 60 %
    assert fit.parameters.profile_used == "breeder" and fit.ebv_worker is not None
    assert all(w.code != "breeder_fallback" for w in fit.warnings)
    # seuil exact : « moins de 50 % » est strict, 2 reines sur 4 (50 %) gardent le profil breeder
    fit, _, _ = _fit(_known_matings(4, 2))
    assert fit.parameters.profile_used == "breeder" and fit.ebv_worker is not None
    assert fit.parameters.known_mating_share == 0.5
    assert all(w.code != "breeder_fallback" for w in fit.warnings)


def test_fallback_never_displays_half():
    # une part juste sous 50 % s'arrondit à 50 % (message) ou à 0,5 (valeur publiée) : le chiffre
    # affiché contredirait le repli (« 50 % ... (< 50 %) »)
    fit, _, _ = _fit(_known_matings(101, 50))                   # 50 / 101 = 0,49505 -> arrondi à 50 %
    assert fit.parameters.profile_used == "standard" and fit.warnings[0].code == "breeder_fallback"
    assert fit.warnings[0].message.startswith("49 % ") and fit.parameters.known_mating_share == 0.495
    fit, _, _ = _fit(_known_matings(1001, 500))                 # 500 / 1001 = 0,49950 -> arrondi à 0,5
    assert fit.parameters.profile_used == "standard" and fit.warnings[0].code == "breeder_fallback"
    assert fit.warnings[0].message.startswith("49 % ") and fit.parameters.known_mating_share == 0.499


def test_known_mating_share_is_published():
    # part des reines notées à accouplement connu, arrondie à 3 décimales, y compris en repli
    def inp(n_known):
        qs = [queen(k, 2022, mating=_station(k) if k <= n_known else mating("natural")) for k in range(1, 7)]
        return BlupV2Input(exploitation=1, criteria=[criterion(h2=0.3, repeatability=0.5)], queens=qs,
                           records=[record(k, 0.1 * k) for k in range(1, 7)])
    fit, _, _ = _fit(inp(2))                                    # 2 / 6 : repli breeder_fallback
    assert fit.parameters.profile_used == "standard" and fit.parameters.known_mating_share == 0.333
    fit, _, _ = _fit(inp(4))                                    # 4 / 6 : profil breeder
    assert fit.parameters.profile_used == "breeder" and fit.parameters.known_mating_share == 0.667
    empty = BlupV2Input(exploitation=1, criteria=[criterion()], queens=[queen(1), queen(2)], records=[])
    assert _fit(empty)[0].parameters.known_mating_share == 0.0             # aucune reine notée
    assert _fit(inp(4), fit_standard)[0].parameters.known_mating_share is None   # profil standard


def test_criterion_without_records():
    v = BlupV2Input(exploitation=1, criteria=[criterion()], queens=[queen(1), queen(2)], records=[])
    fit, _, _ = _fit(v)
    assert set(fit.status.values()) == {"none"} and fit.parameters.profile_used == "breeder"
    assert [w.code for w in fit.warnings] == ["few_records"]


@pytest.mark.parametrize("seed", [1, 2])
def test_approximate_reliabilities_track_exact_ones(seed):
    inp, _, _, _ = simulate_breeder_scenario(seed)
    exact, _, _ = _fit(inp, estimate_variances="never")
    approx, _, _ = _fit(inp, estimate_variances="never", exact_limit=0)
    assert exact.parameters.reliability_method == "exact" and approx.parameters.reliability_method == "approx"
    q = list(exact.reliability)
    r_ex = np.array([exact.reliability[k] for k in q])
    r_ap = np.array([approx.reliability[k] for k in q])
    # prototype (graines 1 à 6) : corrélation 0,81-0,87, surestimation moyenne 0,04-0,05
    assert np.corrcoef(r_ex, r_ap)[0, 1] >= 0.75
    assert abs(float(np.mean(r_ap - r_ex))) <= 0.10
    # seule la PEV change : les trois valeurs publiées sont identiques
    for field in ("ebv", "ebv_worker", "ebv_queen"):
        assert np.allclose([getattr(exact, field)[k] for k in q], [getattr(approx, field)[k] for k in q])
    assert any(w.code == "reliability_approx" for w in approx.warnings)
    assert all(w.code != "reliability_approx" for w in exact.warnings)


def test_auto_reml_keeps_the_queen_worker_partition():
    # ≥ 100 reines notées : REML en auto, mais ρ et rG restent a priori (seuls σ²_AW et σ²_pe estimés)
    inp, _, recorded, _ = simulate_breeder_scenario(1, true_parameters=False)
    assert len(recorded) >= 100
    fit, _, _ = _fit(inp)
    p = fit.parameters
    assert p.source == "reml"
    assert p.sigma2_aq / p.sigma2_a == pytest.approx(0.5) and p.rg == pytest.approx(-0.5)


def test_full_reml_when_always():
    inp, _, _, _ = simulate_breeder_scenario(4, true_parameters=False)
    fit, _, _ = _fit(inp, estimate_variances="always")
    p = fit.parameters
    assert p.source in ("reml", "prior")
    if p.source == "prior":
        assert any(w.code in ("reml_not_converged", "low_heritability") for w in fit.warnings)
    else:
        assert abs(p.rg) <= 0.95 and p.sigma2_aq > 0 and 0.0 < p.h2 < 1.0


@pytest.mark.parametrize("mode, seed, n_params", [("auto", 1, 2), ("always", 4, 4)])
def test_reml_estimates_the_partition_only_when_always(monkeypatch, mode, seed, n_params):
    # espion sur l'optimiseur appelé par blup2.model_breeder : auto n'optimise que
    # θ = (log σ²_AW/σ²_e, log σ²_pe/σ²_e), always y ajoute (log ρ, atanh rG/0,95)
    real, sizes = scipy.optimize.minimize, []

    def spy(fun, x0, *args, **kwargs):
        sizes.append(len(x0))
        return real(fun, x0, *args, **kwargs)

    monkeypatch.setattr(model_breeder.scipy.optimize, "minimize", spy)
    inp, _, recorded, _ = simulate_breeder_scenario(seed, true_parameters=False)
    assert len(recorded) >= 100
    fit, _, _ = _fit(inp, estimate_variances=mode)
    assert sizes == [n_params]
    p = fit.parameters
    partition_is_prior = p.sigma2_aq / p.sigma2_a == pytest.approx(0.5) and p.rg == pytest.approx(-0.5)
    if mode == "auto":                              # seuls σ²_AW et σ²_pe estimés : ρ et rG a priori
        assert p.source == "reml" and partition_is_prior
    elif p.source == "reml":                        # graine 4 : ρ ≈ 0,40, rG ≈ −0,90 estimés
        assert not partition_is_prior
    else:                                           # estimation inutilisable (|rG| ≥ 0,94...) : repli propre
        assert any(w.code in ("reml_not_converged", "low_heritability") for w in fit.warnings)


@pytest.mark.parametrize("mode, rg_theta, success", [("auto", None, False), ("always", 3.0, True)])
def test_unusable_reml_falls_back_to_the_prior(monkeypatch, mode, rg_theta, success):
    # optimiseur simulé : échec (success=False) ou |rG| = 0,95 tanh(3) = 0,945 ≥ 0,94 -> a priori
    # et un seul reml_not_converged (ni low_heritability, ni paramètres estimés)
    def fake(fun, x0, *args, **kwargs):
        x = np.array(x0, dtype=float)
        if rg_theta is not None:
            x[3] = rg_theta
        return scipy.optimize.OptimizeResult(x=x, success=success, nfev=7, fun=0.0)

    inp, _, _, _ = simulate_breeder_scenario(1, true_parameters=False)
    prior, _, _ = _fit(inp, estimate_variances="never")
    monkeypatch.setattr(model_breeder.scipy.optimize, "minimize", fake)
    fit, _, _ = _fit(inp, estimate_variances=mode)
    p = fit.parameters
    assert p.source == "prior" and p.profile_used == "breeder"
    assert p.h2 == pytest.approx(0.25) and p.repeatability == pytest.approx(0.30)
    assert p.sigma2_aq / p.sigma2_a == pytest.approx(0.5) and p.rg == pytest.approx(-0.5)
    assert p.sigma2_a == pytest.approx(prior.parameters.sigma2_a)
    assert p.sigma2_pe == pytest.approx(prior.parameters.sigma2_pe)
    assert [(w.code, w.criterion) for w in fit.warnings] == [("reml_not_converged", 1)]
    assert "7 évaluations" in fit.warnings[0].message
    assert all(fit.ebv[k] == pytest.approx(prior.ebv[k]) for k in prior.ebv)


def test_spd_factorisation_gives_the_reference_likelihood():
    # neg2_reml_spd (sans pivotage) == variance.neg2_reml_loglik (plan standard) sur le système breeder
    import scipy.sparse as sp

    from blup2.model_breeder import neg2_reml_spd
    from blup2.model_standard import incidence
    from blup2.variance import neg2_reml_loglik

    inp, _, _, _ = simulate_breeder_scenario(3)
    v = BlupV2Input.model_validate(inp)
    ped, _ = build_pedigree(v.queens, 12)
    cdata = prepare_criterion(v, 1)
    n = ped.n_nodes
    rec = {q: k for k, q in enumerate(cdata.recorded_queens)}
    Z = [incidence([ped.w_node(q) for q in cdata.queen_ids], n),
         incidence([ped.q_node(q) for q in cdata.queen_ids], n),
         incidence([rec[q] for q in cdata.queen_ids], len(rec))]
    g0 = np.array([[0.4, -0.15], [-0.15, 0.25]])
    Gi = sp.block_diag([genetic_inverse(g0, inverse_relationship(ped)), sp.identity(len(rec)) / 0.3],
                       format="csr")
    logdet = n * float(np.linalg.slogdet(g0)[1]) + len(rec) * float(np.log(0.3))
    ref, s2e_ref = neg2_reml_loglik(cdata.X, Z, Gi, logdet, cdata.y)
    val, s2e = neg2_reml_spd(cdata.X, Z, Gi, logdet, cdata.y)
    assert val == pytest.approx(ref, rel=1e-9) and s2e == pytest.approx(s2e_ref, rel=1e-9)


def test_null_permanent_environment_is_a_usable_estimate():
    # σ²_pe vraie nulle (répétabilité = h²) : σ²_pe finit au bord inférieur, estimation valide
    inp, _, recorded, _ = simulate_breeder_scenario(8, s2pe=0.0, true_parameters=False)
    assert len(recorded) >= 100
    fit, _, _ = _fit(inp)
    p = fit.parameters
    assert p.source == "reml" and p.sigma2_pe < 1e-3 * p.sigma2_e
    assert all(w.code not in ("reml_not_converged", "low_heritability") for w in fit.warnings)


def test_null_genetic_variance_gives_low_heritability():
    # σ²_AW et σ²_AQ vraies quasi nulles : σ²_AW au bord inférieur -> a priori et low_heritability
    inp, _, _, _ = simulate_breeder_scenario(1, s2aw=1e-8, s2aq=0.5e-8, true_parameters=False)
    fit, _, _ = _fit(inp)
    assert fit.parameters.source == "prior" and fit.parameters.profile_used == "breeder"
    assert [w.code for w in fit.warnings] == ["low_heritability"]


def test_constant_notes_give_no_variation_without_reml():
    inp, _, _, _ = simulate_breeder_scenario(1, true_parameters=False)
    for r in inp["records"]:
        r["value"] = r["scale_min"]
    fit, _, _ = _fit(inp)
    assert fit.parameters.source == "prior" and fit.parameters.profile_used == "breeder"
    assert [(w.code, w.criterion) for w in fit.warnings] == [("no_variation", 1)]
    assert max(abs(e) for e in fit.ebv.values()) < 1e-6 * fit.sigma_index
