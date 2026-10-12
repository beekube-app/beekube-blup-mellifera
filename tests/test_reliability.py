import numpy as np
import pytest

from blup2.kinship import relationship_diagonal
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
