import numpy as np
import pytest
import scipy.sparse as sp

from blup2.variance import (
    NO_VARIATION_SIGMA2_P,
    VarianceComponents,
    neg2_reml_loglik,
    phenotypic_variance,
    prior_components,
    reml,
    should_estimate,
)


def test_prior_components_follow_h2_and_repeatability():
    # σ²_P = 2, h² = 0,25, r = 0,40 : σ²_a = 0,5 ; σ²_pe = 0,3 ; σ²_e = 1,2 ; σ²_A = 0,5 / 0,3125 = 1,6
    c = prior_components(2.0, 0.25, 0.40, 0.3125)
    assert (c.sigma2_A, c.sigma2_pe, c.sigma2_e) == pytest.approx((1.6, 0.3, 1.2))
    assert (c.sigma2_a, c.h2, c.repeatability) == pytest.approx((0.5, 0.25, 0.40))


def test_repeatability_below_h2_is_raised_to_h2():
    c = prior_components(1.0, 0.3, 0.2, 1.0)
    assert c.repeatability == pytest.approx(0.3, abs=1e-5) and c.sigma2_pe == pytest.approx(1e-6)


def test_phenotypic_variance_is_residual_after_fixed_effects():
    # 2 groupes {0, 2} et {10, 12} : résidus ±1, ddl = 4 − 2 -> 4 / 2 = 2
    X = sp.csr_matrix(np.array([[1, 0], [1, 0], [1, 1], [1, 1]], dtype=float))
    assert phenotypic_variance(X, np.array([0.0, 2.0, 10.0, 12.0])) == pytest.approx(2.0)


def test_phenotypic_variance_of_constant_notes():
    # notes identiques : plancher 1e-8 par défaut (σ²_P sert de diviseur), 0 avec floor=0 (no_variation)
    X = sp.csr_matrix(np.array([[1, 0], [1, 0], [1, 1], [1, 1]], dtype=float))
    y = np.full(4, 0.6)
    assert phenotypic_variance(X, y) == 1e-8
    assert phenotypic_variance(X, y, floor=0.0) < NO_VARIATION_SIGMA2_P


@pytest.mark.parametrize("mode,args,expected", [
    ("auto", (100, 300, 2), True), ("auto", (99, 300, 2), False), ("auto", (100, 299, 2), False),
    ("auto", (100, 300, 1), False), ("always", (1, 1, 1), True), ("never", (999, 9999, 9), False)])
def test_should_estimate_thresholds(mode, args, expected):
    assert should_estimate(mode, *args) is expected


def _dense_neg2(X, Z_list, G_list, s2e, y):
    # −2 log L_REML = log|V| + log|Xᵀ V⁻¹ X| + yᵀ P y (à une constante près), V = Σ Z G Zᵀ + σ²_e I
    V = s2e * np.eye(len(y)) + sum(Z @ G @ Z.T for Z, G in zip(Z_list, G_list))
    Vi = np.linalg.inv(V)
    XtViX = X.T @ Vi @ X
    P = Vi - Vi @ X @ np.linalg.solve(XtViX, X.T @ Vi)
    return np.linalg.slogdet(V)[1] + np.linalg.slogdet(XtViX)[1] + y @ P @ y


def test_neg2_loglik_matches_dense_formula_up_to_a_constant():
    rng = np.random.default_rng(3)
    n, q = 30, 8
    X = np.column_stack([np.ones(n), rng.integers(0, 2, n)]).astype(float)
    Z = np.zeros((n, q))
    Z[np.arange(n), rng.integers(0, q, n)] = 1.0
    M = rng.standard_normal((q, q))
    K = M @ M.T / q + np.eye(q)
    y = rng.standard_normal(n)
    Kinv = np.linalg.inv(K)
    diffs = []
    for s2u, s2e in [(0.5, 1.0), (2.0, 0.7), (0.1, 3.0)]:
        Gstar_inv = sp.csr_matrix(Kinv * s2e / s2u)          # G* = G / σ²_e
        logdet_gstar = np.linalg.slogdet(K * s2u / s2e)[1]
        ours, _ = neg2_reml_loglik(sp.csr_matrix(X), [sp.csr_matrix(Z)], Gstar_inv, logdet_gstar, y, s2e)
        diffs.append(ours - _dense_neg2(X, [Z], [K * s2u], s2e, y))
    # la constante (N − p) log 1 + log|XᵀX| ne dépend pas des variances
    assert max(diffs) - min(diffs) == pytest.approx(0.0, abs=1e-8)


def test_reml_recovers_simulated_variances():
    # 400 colonies × 4 notes, effet colonie σ²_u = 1, résiduelle σ²_e = 2 ; K = I
    rng = np.random.default_rng(11)
    q, k = 400, 4
    u = rng.standard_normal(q)
    rows = np.repeat(np.arange(q), k)
    y = 3.0 + u[rows] + rng.standard_normal(q * k) * np.sqrt(2.0)
    X = sp.csr_matrix(np.ones((q * k, 1)))
    Z = sp.csr_matrix((np.ones(q * k), (np.arange(q * k), rows)), shape=(q * k, q))
    res = reml(X, [Z], [sp.identity(q, format="csr")], y, start=[0.5], sigma2_e_start=1.0)
    assert res.converged
    assert res.variances[0] == pytest.approx(1.0, abs=0.25)
    assert res.sigma2_e == pytest.approx(2.0, abs=0.25)


def test_reml_without_signal_is_not_converged():
    # aucune variance entre colonies : le rapport σ²_u / σ²_e tend vers la borne basse
    rng = np.random.default_rng(5)
    q, k = 50, 2
    rows = np.repeat(np.arange(q), k)
    y = rng.standard_normal(q * k)
    y = y - np.bincount(rows, y)[rows] / k + rng.standard_normal(q * k) * 1e-3   # moyennes par colonie nulles
    X = sp.csr_matrix(np.ones((q * k, 1)))
    Z = sp.csr_matrix((np.ones(q * k), (np.arange(q * k), rows)), shape=(q * k, q))
    res = reml(X, [Z], [sp.identity(q, format="csr")], y, start=[0.5], sigma2_e_start=1.0, max_evals=400)
    assert res.converged is False


def test_variance_components_properties():
    c = VarianceComponents(sigma2_A=2.0, sigma2_pe=0.5, sigma2_e=1.0, mean_aww=0.25)
    assert c.sigma2_a == 0.5 and c.sigma2_p == 2.0 and c.h2 == 0.25 and c.repeatability == 0.5


def test_em_reml_is_the_reml_entry_point():
    from blup2.variance import em_reml
    assert em_reml is reml


def test_reml_with_null_permanent_environment_is_usable():
    # 100 familles × 4 colonies × 3 notes : effet famille σ²_f = 1, aucun effet colonie propre
    # (environnement permanent nul, répétabilité = h² comme sur l'exploitation 798) : la variance
    # secondaire est au bord, l'estimation reste utilisable et la variance principale est retrouvée.
    rng = np.random.default_rng(21)
    nf, nc, k = 100, 4, 3
    fam = np.repeat(np.arange(nf), nc * k)
    col = np.repeat(np.arange(nf * nc), k)
    e = rng.standard_normal(nf * nc * k)
    e = e - np.bincount(col, e)[col] / k                    # aucune variance entre colonies d'une famille
    y = 2.0 + rng.standard_normal(nf)[fam] + e
    n = len(y)
    X = sp.csr_matrix(np.ones((n, 1)))
    Zf = sp.csr_matrix((np.ones(n), (np.arange(n), fam)), shape=(n, nf))
    Zc = sp.csr_matrix((np.ones(n), (np.arange(n), col)), shape=(n, nf * nc))
    res = reml(X, [Zf, Zc], [sp.identity(nf, format="csr"), sp.identity(nf * nc, format="csr")], y,
               start=[0.5, 0.3], sigma2_e_start=1.0)
    assert res.variances[1] < 0.01
    assert res.converged
    assert res.variances[0] == pytest.approx(1.0, abs=0.3)


def test_reml_reports_a_genetic_variance_at_the_lower_bound():
    # genetic_at_bound : le premier effet de Z_list (génétique) finit au bord inférieur (fit_standard
    # émet alors low_heritability) ; faux dès que l'effet existe
    rng = np.random.default_rng(5)
    q, k = 50, 2
    rows = np.repeat(np.arange(q), k)
    y = rng.standard_normal(q * k)
    y = y - np.bincount(rows, y)[rows] / k + rng.standard_normal(q * k) * 1e-3   # moyennes par colonie nulles
    X = sp.csr_matrix(np.ones((q * k, 1)))
    Z = sp.csr_matrix((np.ones(q * k), (np.arange(q * k), rows)), shape=(q * k, q))
    K = [sp.identity(q, format="csr")]
    null = reml(X, [Z], K, y, start=[0.5], sigma2_e_start=1.0)
    assert null.genetic_at_bound is True and null.converged is False
    res = reml(X, [Z], K, y + rng.standard_normal(q)[rows], start=[0.5], sigma2_e_start=1.0)
    assert res.genetic_at_bound is False and res.converged is True
