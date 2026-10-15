import numpy as np
import scipy.sparse as sp

from blup2.mme import solve_mme


def _problem(seed=0, n=40, p=3, q=12):
    rng = np.random.default_rng(seed)
    X = np.column_stack([np.ones(n)] + [rng.integers(0, 2, n) for _ in range(p - 1)]).astype(float)
    Z = np.zeros((n, q))
    Z[np.arange(n), rng.integers(0, q, n)] = 1.0
    M = rng.standard_normal((q, q))
    G = M @ M.T / q + np.eye(q)
    y = rng.standard_normal(n)
    return X, Z, G, y


def _dense_reference(X, Z, G, s2e, y):
    # Henderson : V = Z G Zᵀ + σ²_e I ; b = (Xᵀ V⁻¹ X)⁻¹ Xᵀ V⁻¹ y ; u = G Zᵀ V⁻¹ (y − X b)
    V = Z @ G @ Z.T + s2e * np.eye(len(y))
    Vi = np.linalg.inv(V)
    b = np.linalg.solve(X.T @ Vi @ X, X.T @ Vi @ y)
    u = G @ Z.T @ Vi @ (y - X @ b)
    # PEV = G − G Zᵀ P Z G, P = V⁻¹ − V⁻¹ X (Xᵀ V⁻¹ X)⁻¹ Xᵀ V⁻¹
    P = Vi - Vi @ X @ np.linalg.inv(X.T @ Vi @ X) @ X.T @ Vi
    pev = np.diag(G - G @ Z.T @ P @ Z @ G)
    return b, u, pev


def test_solutions_and_pev_match_henderson_dense_formulas():
    X, Z, G, y = _problem()
    s2e = 0.7
    sol = solve_mme(sp.csr_matrix(X), [sp.csr_matrix(Z)], sp.csr_matrix(np.linalg.inv(G)), s2e, y, True)
    b, u, pev = _dense_reference(X, Z, G, s2e, y)
    assert np.allclose(sol.beta, b, atol=1e-9)
    assert np.allclose(sol.u[0], u, atol=1e-9)
    assert sol.exact and np.allclose(sol.pev[0], pev, atol=1e-9)


def test_two_random_blocks_are_split():
    X, Z, G, y = _problem(seed=1)
    Z2 = np.eye(len(y))[:, :5]
    Ginv = sp.block_diag([sp.csr_matrix(np.linalg.inv(G)), sp.identity(5) / 0.3], format="csr")
    sol = solve_mme(X, [sp.csr_matrix(Z), sp.csr_matrix(Z2)], Ginv, 1.0, y, True)
    assert [len(u) for u in sol.u] == [12, 5] and [len(p) for p in sol.pev] == [12, 5]
    assert sol.n_equations == 3 + 12 + 5


def test_no_inverse_above_exact_limit():
    X, Z, G, y = _problem()
    sol = solve_mme(X, [sp.csr_matrix(Z)], sp.csr_matrix(np.linalg.inv(G)), 0.7, y, True, exact_limit=10)
    assert sol.exact is False and sol.pev is None and sol.cinv is None
    b, u, _ = _dense_reference(X, Z, G, 0.7, y)
    assert np.allclose(sol.u[0], u, atol=1e-9)        # les solutions restent exactes


def test_wrong_ginv_dimension_is_rejected():
    X, Z, G, y = _problem()
    import pytest
    with pytest.raises(ValueError):
        solve_mme(X, [sp.csr_matrix(Z)], sp.identity(3), 1.0, y, False)


def test_inverse_entries_match_the_dense_inverse_above_the_limit():
    """Spec volet 3 §5 : PEV exactes de quelques nœuds par résolutions du système factorisé."""
    from blup2.mme import inverse_entries
    X, Z, G, y = _problem(seed=2)
    Ginv = sp.csr_matrix(np.linalg.inv(G))
    dense = solve_mme(X, [sp.csr_matrix(Z)], Ginv, 0.7, y, True)
    big = solve_mme(X, [sp.csr_matrix(Z)], Ginv, 0.7, y, True, exact_limit=10)
    rows, cols = [3, 5, 14, 3], [3, 5, 14, 14]
    got = inverse_entries(big, rows, cols, 0.7, chunk=2)
    assert np.allclose(got, dense.cinv[rows, cols], atol=1e-12)
