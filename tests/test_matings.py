import pytest

from blup2.connectivity import connectivity
from blup2.indexes import compute_indexes
from blup2.matings import evaluate_pairs, suggest_matings
from blup2.pedigree import build_pedigree
from blup2.schemas import BlupV2Input, MatingPairIn, OptionsIn
from tests.builders import criterion, queen, record
from tests.fits import make_fit


def _setup(n_suggestions=10, penalty=100.0, max_f=0.0625, flags=None):
    flags = flags or {}
    # 1 : reine de base ; 2, 3 : filles de 1 (sœurs, a = 1/4 avec sire inconnu) ; 4 : sans lien
    qs = [queen(1, 2020), queen(2, 2022, mother=1, **flags.get(2, {})), queen(3, 2022, mother=1, **flags.get(3, {})),
          queen(4, 2022, **flags.get(4, {}))]
    inp = BlupV2Input(exploitation=1, criteria=[criterion(1)], queens=qs,
                      options=OptionsIn(n_suggestions=n_suggestions, inbreeding_penalty=penalty,
                                        max_inbreeding=max_f))
    ped, _ = build_pedigree(inp.queens, 12)
    # index_queen : base = années 2020 et 2022 (reines 1 à 4), m = (0 + 1 + 2 − 3) / 4 = 0, σ = 1
    # -> 1 : 100, 2 : 110, 3 : 120, 4 : 70
    fit = make_fit(1, {1: 0.0, 2: 1.0, 3: 2.0, 4: -3.0})
    return inp, ped, compute_indexes(inp, {1: fit})


def test_expected_index_and_inbreeding_of_a_pair():
    inp, ped, t = _setup()
    m = evaluate_pairs(inp, ped, t, [MatingPairIn(mother=2, drone_source=3)])[0]
    # index attendu = ½ 110 + ½ 120 = 115 ; F = ½ a(Q2, Q3) = ½ · ¼ = 0,125 ; score = 115 − 100 · 0,125
    assert m.expected_total_index == pytest.approx(115.0)
    assert m.expected_indexes == {"1": pytest.approx(115.0)}
    assert m.expected_inbreeding == pytest.approx(0.125)
    assert m.score == pytest.approx(102.5)
    assert m.warnings == ["high_inbreeding"]                # 0,125 ≥ 0,0625


def test_unrelated_pair_has_no_inbreeding():
    inp, ped, t = _setup()
    m = evaluate_pairs(inp, ped, t, [MatingPairIn(mother=2, drone_source=4)])[0]
    assert m.expected_inbreeding == 0.0 and m.score == pytest.approx(90.0) and m.warnings == []


def test_suggestions_are_ranked_per_mother_and_exclude_self():
    inp, ped, t = _setup(n_suggestions=2)
    out = suggest_matings(inp, ped, t)
    by_mother = {}
    for m in out:
        by_mother.setdefault(m.mother, []).append(m)
    assert set(by_mother) == {1, 2, 3, 4} and all(len(v) == 2 for v in by_mother.values())
    assert all(m.mother != m.drone_source for m in out)
    # mère 2 : 3 -> 115 − 12,5 = 102,5 ; 1 -> 105 − 100·½·a(Q2,Q1)=105 − 25 = 80 ; 4 -> 90 -> [3, 4]
    assert [m.drone_source for m in by_mother[2]] == [3, 4]
    assert by_mother[2][0].score >= by_mother[2][1].score


def test_penalty_changes_the_ranking():
    inp, ped, t = _setup(n_suggestions=1, penalty=300.0)
    best = {m.mother: m.drone_source for m in suggest_matings(inp, ped, t)}
    assert best[2] == 4                    # 115 − 37,5 = 77,5 < 90


def test_unavailable_dead_and_eliminated_queens_are_not_candidates():
    flags = {2: {"died": True}, 3: {"available_as_drone_source": False}, 4: {"eliminated": True}}
    inp, ped, t = _setup(flags=flags)
    out = suggest_matings(inp, ped, t)
    # mères possibles {1, 3}, source possible {1} ; 1 × 1 exclu -> seul 3 × 1
    assert [(m.mother, m.drone_source) for m in out] == [(3, 1)]


def test_no_suggestion_when_disabled():
    inp, ped, t = _setup(n_suggestions=0)
    assert suggest_matings(inp, ped, t) == []


def test_connectivity_components():
    # ruchers 1 et 2 reliés par la reine 1 notée dans les deux ; rucher 3 relié par la sœur 3 (a = 1/4) ;
    # rucher 4 : reine 4 sans lien -> 2 composantes
    qs = [queen(1, 2020), queen(2, 2022, mother=1), queen(3, 2022, mother=1), queen(4, 2022)]
    recs = [record(1, 0.5, apiary=1), record(1, 0.5, apiary=2, date="2024-06-01"), record(2, 0.5, apiary=2),
            record(3, 0.5, apiary=3), record(4, 0.5, apiary=4)]
    inp = BlupV2Input(exploitation=1, criteria=[criterion()], queens=qs, records=recs)
    ped, _ = build_pedigree(inp.queens, 12)
    c = connectivity(inp, ped)
    assert c.n_components == 2 and c.components == [[1, 2, 3], [4]]


def _two_apiaries(n_main, n_minor):
    # reines de base non apparentées : n_main notées au rucher 1, n_minor au rucher 2
    n = n_main + n_minor
    qs = [queen(k, 2022) for k in range(1, n + 1)]
    recs = [record(k, 0.5, apiary=1 if k <= n_main else 2) for k in range(1, n + 1)]
    return BlupV2Input(exploitation=1, criteria=[criterion()], queens=qs, records=recs)


def test_connectivity_minor_share():
    # part des reines notées hors de la composante qui en compte le plus, arrondie à 3 décimales
    for (n_main, n_minor), expected in {(99, 1): 0.01, (70, 30): 0.3, (2, 1): 0.333, (5, 0): 0.0}.items():
        inp = _two_apiaries(n_main, n_minor)
        c = connectivity(inp, build_pedigree(inp.queens, 12)[0])
        assert (c.n_components, c.minor_share) == ((2 if n_minor else 1), expected)


def test_connectivity_without_records():
    inp = BlupV2Input(exploitation=1, criteria=[criterion()], queens=[queen(1)])
    ped, _ = build_pedigree(inp.queens, 12)
    assert connectivity(inp, ped).n_components == 0


def test_suggestion_count_when_the_mother_is_also_a_drone_source():
    # 4 sources dont la mère elle-même : au plus 3 suggestions par mère, même si n_suggestions = 10
    inp, ped, t = _setup(n_suggestions=10)
    out = suggest_matings(inp, ped, t)
    assert len(out) == 4 * 3
    inp, ped, t = _setup(n_suggestions=3)
    assert len(suggest_matings(inp, ped, t)) == 4 * 3


def test_connectivity_same_apiary_two_years():
    # les nœuds sont des groupes (rucher, année) : deux années d'un même rucher sans reine commune
    # ni apparentée ne sont pas comparables -> deux composantes, le rucher figure dans les deux
    qs = [queen(1, 2020), queen(2, 2020)]
    recs = [record(1, 0.5, apiary=1, date="2023-06-01"), record(2, 0.5, apiary=1, date="2024-06-01")]
    inp = BlupV2Input(exploitation=1, criteria=[criterion()], queens=qs, records=recs)
    c = connectivity(inp, build_pedigree(inp.queens, 12)[0])
    assert c.n_components == 2 and c.minor_share == 0.5 and c.components == [[1], [1]]
    # une reine notée les deux années relie les groupes
    inp = BlupV2Input(exploitation=1, criteria=[criterion()], queens=qs,
                      records=recs + [record(1, 0.5, apiary=1, date="2024-07-01")])
    c = connectivity(inp, build_pedigree(inp.queens, 12)[0])
    assert c.n_components == 1 and c.components == [[1]]
