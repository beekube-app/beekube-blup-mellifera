import pytest
from pydantic import ValidationError

from blup2.pedigree import build_pedigree, find_cycle
from blup2.schemas import BlupV2Input
from tests.builders import beebreed_queens, mating, queen


def test_father_of_queen_is_mating_of_her_mother():
    # M inséminée avec les mâles de P : P est le père des ouvrières de M ET des filles de M,
    # pas le père de M (correction de blup.r:324).
    qs = [queen(1), queen(2), queen(3, mother=1, mating=mating("insemination", [2])),
          queen(4, mother=3)]
    ped, warnings = build_pedigree(qs, 12)
    assert warnings == []
    q3, q4, w3 = ped.q_node(3), ped.q_node(4), ped.w_node(3)
    assert ped.dam[q3] == ped.q_node(1) and ped.sire[q3] is None      # mère 1 en fécondation inconnue
    assert ped.sire_of_mating(3) == ped.q_node(2)                      # RPM unique : sire = Q(2)
    assert (ped.dam[w3], ped.sire[w3]) == (q3, ped.q_node(2))
    assert (ped.dam[q4], ped.sire[q4]) == (q3, ped.q_node(2))
    assert ped.D[w3] == 12 and ped.D[q4] == 12


def test_natural_or_unknown_mating_has_no_sire_even_with_drone_queens():
    qs = [queen(1), queen(2, mating=mating("natural", [1])), queen(3, mating=mating("unknown"))]
    ped, _ = build_pedigree(qs, 12)
    assert ped.sire_of_mating(2) is None and ped.sire_of_mating(3) is None
    assert ped.sire[ped.w_node(2)] is None


def test_node_kinds_and_counts():
    ped, _ = build_pedigree(beebreed_queens(), 12)
    kinds = [ped.kind.count(k) for k in "QSW"]
    assert kinds == [19, 4, 19]   # 4 + 7 + 6 + 2 = 19 reines, 4 groupes de RPM (S0a, S0b, S1, S2), 19 W
    s1 = ped.sire_of_mating(11)
    assert ped.kind[s1] == "S" and ped.S[s1] == 8
    assert ped.dam[s1] == ped.q_node(1)         # S1 : RPM filles de Q0_1 ...
    assert ped.sire[s1] == ped.sire_of_mating(1)  # ... et de l'accouplement de Q0_1 (S0a)


def test_same_station_year_and_mother_share_one_sire_node():
    m = dict(type="station", station=7, year=2024, mother=1, n_rpm=6)
    qs = [queen(1), queen(2, mating=mating(**m)), queen(3, mating=mating(**m)),
          queen(4, mating=mating(**{**m, "year": 2025}))]
    ped, _ = build_pedigree(qs, 12)
    assert ped.sire_of_mating(2) == ped.sire_of_mating(3) != ped.sire_of_mating(4)


def test_external_ref_identifies_the_sire_group():
    qs = [queen(1, mating=mating("station", ext="LIGNEE-A", n_rpm=8)),
          queen(2, mating=mating("station", ext="LIGNEE-A", n_rpm=8))]
    ped, _ = build_pedigree(qs, 12)
    s = ped.sire_of_mating(1)
    assert s == ped.sire_of_mating(2) and ped.dam[s] is None and ped.S[s] == 8


def test_several_listed_drone_queens_with_common_mother_form_a_group():
    qs = [queen(1), queen(2, mother=1), queen(3, mother=1), queen(4, mating=mating("insemination", [2, 3]))]
    ped, _ = build_pedigree(qs, 12)
    s = ped.sire_of_mating(4)
    assert ped.kind[s] == "S" and ped.S[s] == 2 and ped.dam[s] == ped.q_node(1)


def test_absent_parents_are_added_as_base_queens_with_warning():
    qs = [queen(10, mother=3257, mating=mating("insemination", [3301], n_drones=8))]
    ped, warnings = build_pedigree(qs, 12)
    assert ped.added_queen_ids == [3257, 3301]
    assert ped.dam[ped.q_node(3257)] is None and ped.sire[ped.q_node(3257)] is None
    assert [w.code for w in warnings] == ["unknown_parent_added"]
    assert warnings[0].message.endswith("(3257, 3301)")
    assert ped.queen_ids == [10]
    assert ped.D[ped.w_node(10)] == 8          # D vient de l'accouplement de la reine 10


def test_nodes_are_in_topological_order_whatever_the_input_order():
    qs = list(reversed(beebreed_queens()))
    ped, _ = build_pedigree(qs, 12)
    for i in range(ped.n_nodes):
        for p in (ped.dam[i], ped.sire[i]):
            assert p is None or p < i


def test_cycles_are_detected_through_mother_and_mating_links():
    assert find_cycle([queen(1, mother=2), queen(2, mother=1)]) == [1, 2]
    # 2 est fille de 1 et productrice de mâles pour 1 : W(1) dépend de Q(2) qui dépend de W(1)
    assert find_cycle([queen(1, mating=mating("insemination", [2])), queen(2, mother=1)]) == [1, 2]
    # fécondations croisées (1 avec les mâles de 2, 2 avec ceux de 1) : pas de cycle
    assert find_cycle([queen(1, mating=mating("insemination", [2])),
                       queen(2, mating=mating("insemination", [1]))]) == []


def test_cyclic_pedigree_is_a_validation_error():
    d = {"exploitation": 1, "criteria": [], "queens": [
        {"id": 1, "born": 2020, "mother": 3}, {"id": 2, "born": 2020, "mother": 1},
        {"id": 3, "born": 2020, "mother": 2}, {"id": 4, "born": 2020, "mother": 1}]}
    with pytest.raises(ValidationError) as exc:
        BlupV2Input.model_validate(d)
    assert "pedigree cyclique : reines concernées 1, 2, 3" in str(exc.value)


def test_shared_external_ref_uses_the_first_declared_mother_of_drone_queens():
    # même lot de mâles (external_ref) déclaré avec deux mères de RPM différentes : un seul nœud sire,
    # dont la mère est celle du premier accouplement rencontré (ordre d'entrée)
    qs = [queen(1), queen(2),
          queen(3, mating=mating("station", mother=1, n_rpm=6, ext="lot-A")),
          queen(4, mating=mating("station", mother=2, n_rpm=6, ext="lot-A"))]
    ped, _ = build_pedigree(qs, 12)
    s = ped.sire_of_mating(3)
    assert ped.sire_of_mating(4) == s
    assert ped.dam[s] == ped.q_node(1) and ped.S[s] == 6


def test_self_references_are_removed_with_one_warning():
    # 5 : sa propre mère ; 6 : parmi ses RPM (avec 2) ; 7 : mère de ses RPM ; 8 : sa seule RPM
    qs = [queen(1), queen(2, mother=1), queen(5, mother=5),
          queen(6, mating=mating("insemination", [6, 2])),
          queen(7, mating=mating("station", mother=7, n_rpm=8, station=3)),
          queen(8, mother=1, mating=mating("insemination", [8]))]
    assert find_cycle(qs) == []
    # les cycles de longueur ≥ 2 restent détectés
    assert find_cycle(qs + [queen(30, mother=31), queen(31, mother=30)]) == [30, 31]
    ped, warnings = build_pedigree(qs, 12)
    assert [(w.code, w.criterion) for w in warnings] == [("self_reference_removed", None)]
    assert warnings[0].message.endswith("(5, 6, 7, 8)")
    assert ped.dam[ped.q_node(5)] is None
    assert ped.sire_of_mating(6) == ped.q_node(2)               # RPM restante : 2
    s7 = ped.sire_of_mating(7)
    assert ped.kind[s7] == "S" and ped.dam[s7] is None and ped.S[s7] == 8
    assert ped.sire_of_mating(8) is None                        # plus aucune RPM connue
    assert ped.dam[ped.q_node(8)] == ped.q_node(1)              # la mère, elle, est gardée
    assert ped.queen_ids == [1, 2, 5, 6, 7, 8] and ped.added_queen_ids == []


def test_self_reference_warning_lists_at_most_20_queens():
    _, warnings = build_pedigree([queen(k, mother=k) for k in range(1, 26)], 12)
    assert len(warnings) == 1 and warnings[0].code == "self_reference_removed"
    assert warnings[0].message.startswith("25 reine(s)")
    assert warnings[0].message.endswith("(" + ", ".join(str(k) for k in range(1, 21)) + " … et 5 autres)")


def test_unknown_parent_and_cycle_messages_are_capped():
    qs = [queen(k, mother=1000 + k) for k in range(1, 26)]
    _, warnings = build_pedigree(qs, 12)
    msg = [w for w in warnings if w.code == "unknown_parent_added"][0].message
    assert msg.startswith("25 reine(s)") and msg.endswith(" … et 5 autres)")
    cyc = [queen(k, mother=k % 25 + 1) for k in range(1, 26)]
    with pytest.raises(ValidationError) as exc:
        BlupV2Input.model_validate({"exploitation": 1, "criteria": [],
                                    "queens": [q.model_dump() for q in cyc]})
    assert "… et 5 autres" in str(exc.value)
    with pytest.raises(ValueError, match="… et 5 autres"):
        build_pedigree(cyc, 12)


def test_extra_queens_enter_the_graph_but_not_queen_ids():
    """Ancêtres externes (id négatifs) : nœuds du pedigree comme des reines de base, hors de queen_ids et de
    l'avertissement unknown_parent_added."""
    anc = [queen(-1, 2015, mother=-2, mating=mating("insemination", mother=-3, ext="FR-X-B-0003-2013")),
           queen(-2, 2013), queen(-3, 2013)]
    qs = [queen(10, 2020, mother=-1), queen(11, 2020, mother=3257)]
    ped, warnings = build_pedigree(qs, 12, extra_queens=anc)
    assert ped.queen_ids == [10, 11]
    assert sorted(ped.added_queen_ids) == [-3, -2, -1, 3257]
    assert [w.message for w in warnings if w.code == "unknown_parent_added"] == [
        "1 reine(s) citée(s) comme parent mais absente(s) de queens, ajoutée(s) comme reines de base : (3257)"]
    q10 = ped.q_node(10)
    assert ped.dam[q10] == ped.q_node(-1)
    s = ped.sire[q10]                                   # père de 10 = sire de l'accouplement de -1
    assert ped.kind[s] == "S" and ped.dam[s] == ped.q_node(-3)
