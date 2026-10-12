import numpy as np
import pytest

from blup2.kinship import relationship_diagonal, relationship_matrix
from blup2.pedigree import build_pedigree
from blup2.reliability import approx_pev, effective_records
from tests.builders import queen


def test_effective_records_absorb_fixed_effects():
    # colonie 0 : 2 notes dans le groupe a (n_a = 3) et 1 dans b (n_b = 1)
    # facteur g : 2 − 4/3 + 1 − 1/1 = 2/3 ; colonie 1 : 1 − 1/3 = 2/3
    n = effective_records({"g": ["a", "a", "b", "a"]}, [0, 0, 0, 1], 2)
    assert n.tolist() == pytest.approx([2 / 3, 2 / 3])


def test_effective_records_without_factor_is_the_count():
    assert effective_records({}, [0, 0, 1], 3).tolist() == [2.0, 1.0, 0.0]


def test_single_unrelated_colony_is_exact():
    # reine de base, 4 notes, aucun effet fixe : PEV(W) = 1 / (1/a_WW σ²_A + 1 / (σ²_pe + σ²_e / n))
    ped, _ = build_pedigree([queen(1)], 12)
    w = ped.w_node(1)
    n_eff = np.zeros(ped.n_nodes)
    n_eff[w] = 4.0
    s2A, s2pe, s2e = 1.0, 0.5, 2.0
    aww = 0.25 + 1 / 24                 # W de reine de base, sire inconnu
    expected = 1.0 / (1.0 / (aww * s2A) + 1.0 / (s2pe + s2e / 4.0))
    assert approx_pev(ped, n_eff, s2A, s2pe, s2e)[w] == pytest.approx(expected)


def test_no_information_gives_prior_variance():
    ped, _ = build_pedigree([queen(1), queen(2, mother=1)], 12)
    pev = approx_pev(ped, np.zeros(ped.n_nodes), 1.0, 0.5, 2.0)
    assert pev == pytest.approx(relationship_diagonal(ped))


def test_approx_pev_matches_exact_on_a_tree_pedigree():
    # Mère 1 et deux filles notées (sire inconnu : le pedigree généralisé est un arbre).
    # Exact : PEV = σ²_A · (A⁻¹ + diag(U))⁻¹ avec U_i = 1 / (σ²_pe/σ²_A + σ²_e/(σ²_A n_i)), c'est-à-dire
    # l'information propre de chaque colonie après absorption des effets fixes (n_i = n_eff donné).
    ped, _ = build_pedigree([queen(1), queen(2, mother=1), queen(3, mother=1)], 12)
    s2A, s2pe, s2e = 1.3, 0.5, 2.0
    n_eff = np.zeros(ped.n_nodes)
    n_eff[ped.w_node(2)] = 4.0
    n_eff[ped.w_node(3)] = 2.5
    U = np.where(n_eff > 0, 1.0 / (s2pe / s2A + s2e / (s2A * np.maximum(n_eff, 1e-12))), 0.0)
    A = relationship_matrix(ped)
    exact = s2A * np.diag(np.linalg.inv(np.linalg.inv(A) + np.diag(U)))
    assert approx_pev(ped, n_eff, s2A, s2pe, s2e) == pytest.approx(exact, abs=1e-8)
