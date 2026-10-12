import pytest

from blup2.indexes import base_birth_years, compute_indexes, queen_outputs
from blup2.kinship import inbreeding
from blup2.pedigree import build_pedigree
from blup2.schemas import BlupV2Input, OptionsIn
from blup2.trend import genetic_trend
from tests.builders import criterion, queen
from tests.fits import make_fit


def _inp(criteria=None, base_years=5):
    qs = [queen(1, 2018), queen(2, 2023), queen(3, 2024), queen(4, 2024)]
    return BlupV2Input(exploitation=1, options=OptionsIn(base_years=base_years),
                       criteria=criteria or [criterion(1)], queens=qs)


def test_base_is_the_last_birth_years_present():
    assert base_birth_years(_inp(base_years=2)) == [2023, 2024]
    assert base_birth_years(_inp(base_years=5)) == [2018, 2023, 2024]


def test_index_is_100_plus_10_per_genetic_sd_around_the_base_mean():
    # base 2 ans = reines 2, 3, 4 : m = (0 + 1 + 2) / 3 = 1 ; σ = 2 -> index = 100 + 10 (ebv − 1) / 2
    fit = make_fit(1, {1: -1.0, 2: 0.0, 3: 1.0, 4: 2.0}, sigma=2.0)
    t = compute_indexes(_inp(base_years=2), {1: fit})
    assert [t.index[1][q] for q in (1, 2, 3, 4)] == pytest.approx([90.0, 95.0, 100.0, 105.0])
    assert t.base.birth_years == [2023, 2024] and t.base.n_queens == 3


def test_direction_minus_one_reverses_the_index():
    fit = make_fit(1, {1: -1.0, 2: 0.0, 3: 1.0, 4: 2.0}, sigma=2.0)
    t = compute_indexes(_inp([criterion(1, direction=-1)], base_years=2), {1: fit})
    assert [t.index[1][q] for q in (1, 2, 3, 4)] == pytest.approx([110.0, 105.0, 100.0, 95.0])


def test_queens_without_information_have_no_index_and_are_excluded_from_the_base():
    # reine 4 au statut none : m = moyenne de 2 et 3 = 0,5
    fit = make_fit(1, {1: 0.0, 2: 0.0, 3: 1.0, 4: 9.0},
                   status={1: "pedigree", 2: "evaluated", 3: "evaluated", 4: "none"})
    t = compute_indexes(_inp(base_years=2), {1: fit})
    assert t.index[1][4] is None and t.status[1][4] == "none"
    assert t.index[1][2] == pytest.approx(95.0) and t.index[1][3] == pytest.approx(105.0)


def test_queen_index_uses_queen_sd():
    fit = make_fit(1, {q: 0.0 for q in (1, 2, 3, 4)}, ebv_queen={1: 0.0, 2: 0.0, 3: 0.0, 4: 4.0},
                   sigma=1.0, sigma_queen=4.0)
    t = compute_indexes(_inp(base_years=1), {1: fit})       # base 2024 : m_Q = (0 + 4) / 2 = 2
    assert t.index_queen[1][4] == pytest.approx(105.0) and t.index_queen[1][3] == pytest.approx(95.0)


def test_total_index_is_weighted_over_non_elimination_criteria():
    crits = [criterion(1, weight=3), criterion(2, weight=1), criterion(3, weight=5, elimination=True),
             criterion(4, weight=0)]
    e = {1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0}
    fits = {1: make_fit(1, {**e, 4: 1.0}, reliability={q: 0.8 for q in e}),
            2: make_fit(2, {**e, 4: -1.0}, reliability={q: 0.4 for q in e}),
            3: make_fit(3, {**e, 4: 50.0}), 4: make_fit(4, {**e, 4: 50.0})}
    t = compute_indexes(_inp(crits, base_years=1), fits)    # base : reines 3 et 4
    # critère 1 : m = 0,5 -> reine 4 : 105 ; critère 2 : m = −0,5 -> 95 ; total = (3·105 + 1·95) / 4 = 102,5
    assert t.total[4] == pytest.approx(102.5)
    assert t.total_reliability[4] == pytest.approx((3 * 0.8 + 0.4) / 4)


def test_failed_criterion_is_reported_and_left_out_of_the_total():
    fits = {1: make_fit(1, {1: 0.0, 2: 0.0, 3: 0.0, 4: 1.0}), 2: None}
    t = compute_indexes(_inp([criterion(1), criterion(2)], base_years=1), fits)
    assert t.status[2] == {1: "failed", 2: "failed", 3: "failed", 4: "failed"}
    assert t.total[4] == pytest.approx(t.index[1][4])


def test_queen_outputs_are_rounded_and_limited_to_input_queens():
    inp = BlupV2Input(exploitation=1, criteria=[criterion(1)], queens=[queen(1, 2024, mother=77)])
    ped, _ = build_pedigree(inp.queens, 12)
    fit = make_fit(1, {1: 0.123456789}, sigma=3.0, reliability={1: 0.98765})
    t = compute_indexes(inp, {1: fit})
    out = queen_outputs(inp, ped, {1: fit}, t, inbreeding(ped))
    assert [q.id for q in out] == [1]                    # la mère 77 ajoutée n'est pas restituée
    ci = out[0].indexes["1"]
    assert ci.ebv == 0.123457 and ci.reliability == 0.988 and ci.index == 100.0
    assert out[0].inbreeding == 0.0 and out[0].total_index == 100.0


def test_trend_is_mean_index_by_birth_year():
    fit = make_fit(1, {1: -1.0, 2: 0.0, 3: 1.0, 4: 2.0}, sigma=2.0)
    inp = _inp(base_years=2)
    tr = genetic_trend(inp, compute_indexes(inp, {1: fit}))["1"]
    assert [(p.birth_year, p.mean_index, p.n) for p in tr] == [(2018, 90.0, 1), (2023, 95.0, 1), (2024, 102.5, 2)]
