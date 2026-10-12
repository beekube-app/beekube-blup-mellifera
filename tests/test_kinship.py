"""Parenté abeille : comparaison au code de référence docs/analyses/bee_kinship.py et invariants
de Brascamp & Bijma 2014 (BB)."""
import numpy as np
import pytest

from blup2.kinship import (
    generalized_parents,
    hypothetical_sire_relationship,
    inbreeding,
    inverse_relationship,
    relationship,
    relationship_diagonal,
    relationship_matrix,
)
from blup2.pedigree import build_pedigree
from blup2.schemas import QueenIn
from tests.builders import beebreed_queens, mating, queen
from tests.simdata import scenario_pedigree  # ajoute docs/analyses au sys.path

from verify_paper import beebreed_pedigree  # noqa: E402  (référence vérifiée)


def _beebreed_map(ped):
    name = {}
    for k in range(1, 5):
        name[f"Q0_{k}"], name[f"W_Q0_{k}"] = ped.q_node(k), ped.w_node(k)
    for k in range(1, 8):
        name[f"Q1_{k}"], name[f"W_Q1_{k}"] = ped.q_node(10 + k), ped.w_node(10 + k)
    for k in range(1, 7):
        name[f"Q2_{k}"], name[f"W_Q2_{k}"] = ped.q_node(20 + k), ped.w_node(20 + k)
    for k in range(1, 3):
        name[f"Q3_{k}"] = ped.q_node(30 + k)
    name["S0a"], name["S0b"] = ped.sire_of_mating(1), ped.sire_of_mating(3)
    name["S1"], name["S2"] = ped.sire_of_mating(11), ped.sire_of_mating(21)
    return name


def test_matches_reference_on_beebreed_pedigree():
    ped, _ = build_pedigree(beebreed_queens(), 12)
    ref = beebreed_pedigree(8, 12)
    A_ref = ref.build()
    name = _beebreed_map(ped)
    ids = [i.id for i in ref.order]
    ours = relationship_matrix(ped, [name[x] for x in ids])
    theirs = A_ref[np.ix_([ref.inds[x].idx for x in ids], [ref.inds[x].idx for x in ids])]
    assert np.abs(ours - theirs).max() < 1e-12


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_matches_reference_on_simulated_pedigrees(seed):
    P, inds, q_in, qid = scenario_pedigree(np.random.default_rng(seed))
    A_ref = P.build()
    ped, _ = build_pedigree([QueenIn.model_validate(q) for q in q_in], 12)
    ref_idx, our_idx = [], []
    for x in inds:
        if x["kind"] == "Q":
            ref_idx.append(P.inds[x["id"]].idx)
            our_idx.append(ped.q_node(qid[x["id"]]))
        elif x["kind"] == "W":
            ref_idx.append(P.inds[x["id"]].idx)
            our_idx.append(ped.w_node(qid[x["dam"]]))
    ours = relationship_matrix(ped, our_idx)
    assert np.abs(ours - A_ref[np.ix_(ref_idx, ref_idx)]).max() < 1e-10


def test_inverse_times_matrix_is_identity():
    ped, _ = build_pedigree(beebreed_queens(), 12)
    A = relationship_matrix(ped)
    assert np.abs(A @ inverse_relationship(ped).toarray() - np.eye(ped.n_nodes)).max() < 1e-9


def test_queen_diagonal_is_one_plus_inbreeding():
    ped, _ = build_pedigree(beebreed_queens(), 12)
    d, F = relationship_diagonal(ped), inbreeding(ped)
    q = [ped.q_node(x) for x in ped.queen_ids]
    assert np.abs(d[q] - 1.0 - F[q]).max() < 1e-12
    assert F[ped.q_node(31)] > 0.0          # Q3_1 : fille de Q2_1 × S2, parents apparentés


@pytest.mark.parametrize("S,D", [(8, 12), (8, 16), (4, 12)])
def test_full_sisters_of_base_parents(S, D):
    # BB éq. 24a avec F = 0 et a_ss = 0 : a_FS = 1/4 + 1/(2D) + 1/(4S) ; S = 8, D = 12 -> 0,3229
    qs = [queen(1, mating=mating("station", station=1, n_rpm=S, n_drones=D)),
          queen(2, mother=1), queen(3, mother=1)]
    ped, _ = build_pedigree(qs, 12)
    assert relationship(ped, ped.q_node(2), ped.q_node(3)) == pytest.approx(0.25 + 1 / (2 * D) + 1 / (4 * S))
    # le groupe d'ouvrières de la mère a la même variance que la parenté entre pleines sœurs
    w1 = ped.w_node(1)
    assert relationship(ped, w1, w1) == pytest.approx(0.25 + 1 / (2 * D) + 1 / (4 * S))
    assert relationship(ped, ped.q_node(1), w1) == pytest.approx(0.5)


def test_single_drone_queen_sisters():
    # RPM unique (S = 1) de base : 1/4 + 1/4 a_PP + cov(δ_FS) = 1/4 + 1/4 + 1/(4D) = 0,5208 (D = 12)
    qs = [queen(1, mating=mating("insemination", [9])), queen(9), queen(2, mother=1), queen(3, mother=1)]
    ped, _ = build_pedigree(qs, 12)
    assert relationship(ped, ped.q_node(2), ped.q_node(3)) == pytest.approx(0.5 + 1 / 48)
    assert relationship(ped, ped.q_node(2), ped.q_node(9)) == pytest.approx(0.5)   # père = RPM


def _inbred_line():
    # 2 et 3 : filles de 1 (fécondation inconnue) -> a(2,3) = 1/4 ; 2 inséminée avec les mâles de 3 ;
    # 4 = fille de 2 -> F_4 = 1/2 a(2,3) = 0,125 ; 5 = fille de 4 (fécondation naturelle).
    return [queen(1), queen(2, mother=1, mating=mating("insemination", [3])), queen(3, mother=1),
            queen(4, mother=2, mating=mating("natural")), queen(5, mother=4)]


def test_unknown_sire_mendelian_variance():
    ped, _ = build_pedigree(_inbred_line(), 12)
    parents, phi = generalized_parents(ped)
    F = inbreeding(ped)
    q4, q5 = ped.q_node(4), ped.q_node(5)
    assert F[q4] == pytest.approx(0.125)
    # BB annexe 1 cas 3 : var(δ) = 3/4 − 1/4 F_mère = 0,75 − 0,03125
    assert parents[q5] == [(q4, 0.5)]
    assert phi[q5] == pytest.approx(0.75 - 0.25 * 0.125)
    assert relationship(ped, q5, q5) == pytest.approx(1.0)
    # groupe d'ouvrières, sire inconnu : var(δ) = 1/(2D) -> a_WW = 1/4 (1 + F_4) + 1/24
    w4 = ped.w_node(4)
    assert phi[w4] == pytest.approx(1 / 24)
    assert relationship(ped, w4, w4) == pytest.approx(0.25 * 1.125 + 1 / 24)


def test_worker_inbreeding_is_half_queen_sire_relationship():
    ped, _ = build_pedigree(_inbred_line(), 12)
    F = inbreeding(ped)
    assert F[ped.w_node(2)] == pytest.approx(0.5 * 0.25)     # BB éq. 25 : ½ a(Q2, Q3)
    assert F[ped.w_node(1)] == 0.0                           # sire inconnu


def test_planned_and_recorded_mating_give_the_same_inbreeding():
    planned, _ = build_pedigree([queen(1), queen(2, mother=1), queen(3, mother=1)], 12)
    f_planned = 0.5 * hypothetical_sire_relationship(planned, None, 2, 3, 12)
    recorded, _ = build_pedigree(
        [queen(1), queen(2, mother=1, mating=mating("insemination", [3])), queen(3, mother=1)], 12)
    assert f_planned == pytest.approx(0.125)
    assert inbreeding(recorded)[recorded.w_node(2)] == pytest.approx(f_planned)
    A = relationship_matrix(planned)
    assert hypothetical_sire_relationship(planned, A, 2, 3, 12) == pytest.approx(0.25)


def test_result_does_not_depend_on_input_order():
    # nièce avant tante (docs/analyses/test_order.py) : l'ordre d'entrée ne change pas A
    qs = beebreed_queens()
    p1, _ = build_pedigree(qs, 12)
    p2, _ = build_pedigree(list(reversed(qs)), 12)
    n1 = [p1.q_node(q) for q in p1.queen_ids] + [p1.w_node(q) for q in p1.queen_ids]
    n2 = [p2.q_node(q) for q in p1.queen_ids] + [p2.w_node(q) for q in p1.queen_ids]
    assert np.abs(relationship_matrix(p1, n1) - relationship_matrix(p2, n2)).max() < 1e-12
