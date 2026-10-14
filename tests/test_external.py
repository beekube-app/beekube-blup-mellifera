"""Volet 2 Index Mellifera : résolution des références externes (spec IM §8.1, §8.2, §8.4)."""
import sqlite3

import pytest

from blup2.external import expand_external, normalize_code
from blup2.pedigree import build_pedigree
from tests.builders import mating, queen
from tests.external_base import CYCLES, inp, make_db, warning


@pytest.fixture
def db(tmp_path):
    return make_db(tmp_path / "mellifera.sqlite")


# ------------------------------------------------------------------ résolution d'une référence
@pytest.mark.parametrize("ref", ["FR-CBP-B-0276-2022", " fr-cbp-b-0276-2022 ", "fr-cbp-b-0276-2022\t", "B276(CBP)"])
def test_mother_ref_resolved_by_code_or_unique_name(db, ref):
    src = inp([queen(1, 2024, mother_ref_ext=ref)])
    exp = expand_external(src, db)
    mother = exp.queens[0].mother
    assert mother is not None and mother < 0
    assert {a.id: a.born for a in exp.ancestors}[mother] == 2022
    assert src.queens[0].mother is None                     # l'entrée n'est jamais modifiée
    assert warning(exp, "external_ref_resolved") == [
        "1 référence(s) Index Mellifera rattachée(s), 5 ancêtre(s) ajouté(s) au pedigree"]


def test_normalize_code():
    assert normalize_code(" fr-cbp-b 0276-2022\n") == "FR-CBP-B0276-2022"


@pytest.mark.parametrize("ref", ["TWIN", "FR-XXX-B-0001-2020", "B276"])
def test_ambiguous_name_or_unknown_code_is_reported_and_ignored(db, ref):
    exp = expand_external(inp([queen(1, 2024, mother_ref_ext=ref)]), db)
    assert exp.queens[0].mother is None and exp.ancestors == []
    assert warning(exp, "external_ref_unknown") == [
        f"1 code(s) Index Mellifera introuvable(s), ignoré(s) jusqu'au prochain crawl : {ref}"]
    assert warning(exp, "external_ref_resolved") == []


def test_unknown_codes_are_listed_at_most_twenty(db):
    qs = [queen(k, 2024, mother_ref_ext=f"FR-XXX-B-{k:04d}-2020") for k in range(1, 26)]
    msg = warning(expand_external(inp(qs), db), "external_ref_unknown")[0]
    assert msg.startswith("25 code(s)") and msg.endswith("… et 5 autres")


def test_known_mother_and_blank_refs_are_left_alone(db):
    qs = [queen(1, 2020), queen(2, 2024, mother=1, mother_ref_ext="FR-CBP-B-0276-2022"),
          queen(3, 2024, mother_ref_ext="  ")]
    exp = expand_external(inp(qs), db)
    assert [q.mother for q in exp.queens] == [None, 1, None]
    assert exp.ancestors == [] and exp.warnings == []


# ------------------------------------------------------------------ côté mâle
def test_drone_side_ref_becomes_an_insemination_with_known_drone_queens_mother(db):
    qs = [queen(1, 2024, mating=mating("unknown", ext="FR-CBP-B-0276-2022")),
          queen(2, 2024, mating=mating("unknown", ext="FR-CBP-B-0276-2022"))]
    exp = expand_external(inp(qs), db)
    m = exp.queens[0].mating
    assert m.type == "insemination" and m.drone_queens_mother < 0
    assert m.external_ref == "FR-CBP-B-0276-2022" and m.n_drone_queens is None
    ped, _ = build_pedigree(exp.queens, 12, extra_queens=exp.ancestors)
    s1, s2 = ped.sire_of_mating(1), ped.sire_of_mating(2)
    assert s1 == s2 and ped.kind[s1] == "S" and ped.S[s1] == 1          # nœud sire partagé, S = 1
    assert ped.dam[s1] == ped.q_node(m.drone_queens_mother)


def test_drone_side_ref_typed_differently_shares_the_sire_node_of_the_lineage(db):
    qs = [queen(1, 2024, mating=mating("unknown", ext=" fr-cbp-b-0276-2022")),
          queen(2, 2024, mating=mating("unknown", ext="B276(CBP)")),
          queen(3, 2024, mother_ref_ext="FR-CBP-B-0276-2022")]
    exp = expand_external(inp(qs), db)
    assert {q.mating.external_ref for q in exp.queens[:2]} == {"FR-CBP-B-0276-2022"}
    ped, _ = build_pedigree(exp.queens, 12, extra_queens=exp.ancestors)
    assert ped.sire_of_mating(1) == ped.sire_of_mating(2)


def test_station_type_is_kept_and_natural_or_known_drone_side_is_left_alone(db):
    qs = [queen(1, 2024, mating=mating("station", ext="FR-CBP-B-0276-2022")),
          queen(2, 2024, mating=mating("natural", ext="FR-CBP-B-0276-2022")),
          queen(3, 2024, mating=mating("insemination", mother=1, ext="FR-CBP-B-0276-2022")),
          queen(4, 2024, mating=mating("insemination", drone_queens=[1], ext="FR-CBP-B-0276-2022"))]
    exp = expand_external(inp(qs), db)
    assert exp.queens[0].mating.type == "station" and exp.queens[0].mating.drone_queens_mother < 0
    assert exp.queens[1:] == qs[1:]


# ------------------------------------------------------------------ ascendance
def test_ancestry_mating_built_like_the_export(db):
    exp = expand_external(inp([queen(1, 2024, mother_ref_ext="FR-CBP-B-0276-2022")]), db)
    by_born = {}
    for a in exp.ancestors:
        by_born.setdefault(a.born, []).append(a)
    a276 = next(a for a in exp.ancestors if a.id == exp.queens[0].mother)
    assert a276.mating.type == "insemination" and a276.mating.external_ref == "FR-NG-B-0009-2018"
    assert a276.mating.drone_queens_mother in {a.id for a in by_born[2018]}
    assert all(not a.available_as_mother and not a.available_as_drone_source for a in exp.ancestors)
    assert all(a.id < 0 for a in exp.ancestors) and len({a.id for a in exp.ancestors}) == len(exp.ancestors)


@pytest.mark.parametrize("depth, n", [(1, 1), (2, 3), (3, 4), (7, 5)])
def test_depth_limits_the_generations_and_cuts_links_beyond(db, depth, n):
    exp = expand_external(inp([queen(1, 2024, mother_ref_ext="FR-CBP-B-0276-2022")], depth=depth), db)
    assert len(exp.ancestors) == n
    ids = {a.id for a in exp.ancestors}
    for a in exp.ancestors:
        assert a.mother is None or a.mother in ids
        assert a.mating.drone_queens_mother is None or a.mating.drone_queens_mother in ids


def test_long_chain_is_cut_at_the_maximum_depth(tmp_path):
    rows = [(f"FR-L-B-{k:04d}-{1960 + k}", 1960 + k, f"FR-L-B-{k - 1:04d}-{1959 + k}" if k else None, None, None, None)
            for k in range(60)]
    db = make_db(tmp_path / "m.sqlite", rows)
    exp = expand_external(inp([queen(1, 2024, mother_ref_ext=rows[-1][0])], depth=30), db)
    assert len(exp.ancestors) == 30 and min(a.born for a in exp.ancestors) == 1990


def test_depth_zero_disables_the_resolution(db):
    src = inp([queen(1, 2024, mother_ref_ext="FR-CBP-B-0276-2022")], depth=0)
    exp = expand_external(src, db)
    assert exp.queens == src.queens and exp.ancestors == [] and exp.warnings == []


def test_shared_ancestor_is_added_once(db):
    qs = [queen(1, 2024, mother_ref_ext="FR-CBP-B-0276-2022"), queen(2, 2024, mother_ref_ext="FR-CBP-B-0277-2022"),
          queen(3, 2024, mother_ref_ext="B276(CBP)")]
    exp = expand_external(inp(qs), db)
    assert len(exp.ancestors) == 6                       # 0276, 0277, 0100, 0009, 0050, 0010
    assert exp.queens[0].mother == exp.queens[2].mother != exp.queens[1].mother
    assert warning(exp, "external_ref_resolved")[0].startswith("3 référence(s)")


def test_negative_ids_never_collide_with_input_ids(db):
    qs = [queen(-1, 2020), queen(-5, 2024, mother=-1, mating=mating("insemination", drone_queens=[-9])),
          queen(1, 2024, mother_ref_ext="FR-CBP-B-0276-2022")]
    exp = expand_external(inp(qs), db)
    assert max(a.id for a in exp.ancestors) == -10


def test_ancestor_without_valid_year_stops_the_climb(tmp_path):
    rows = [("FR-A-B-0003-2022", 2022, "FR-A-B-0002-2020", None, None, None),
            ("FR-A-B-0002-2020", None, "FR-A-B-0001-2018", None, None, None),
            ("FR-A-B-0001-2018", 2018, None, None, None, None),
            ("FR-A-B-0009-2022", 3022, None, None, None, None)]
    db = make_db(tmp_path / "m.sqlite", rows)
    exp = expand_external(inp([queen(1, 2024, mother_ref_ext="FR-A-B-0003-2022"),
                               queen(2, 2024, mother_ref_ext="FR-A-B-0009-2022")]), db)
    assert [a.born for a in exp.ancestors] == [2022] and exp.ancestors[0].mother is None
    assert exp.queens[1].mother is None
    assert warning(exp, "external_ref_resolved") == [
        "1 référence(s) Index Mellifera rattachée(s), 1 ancêtre(s) ajouté(s) au pedigree, "
        "2 reine(s) sans année de naissance valide écartée(s)"]


def test_cycles_and_self_references_in_the_base_are_cut(tmp_path):
    db = make_db(tmp_path / "m.sqlite", CYCLES)
    qs = [queen(1, 2024, mother_ref_ext="FR-C-B-0001-2020", mating=mating("unknown", ext="FR-C-B-0004-2017")),
          queen(2, 2024, mother=1)]
    exp = expand_external(inp(qs), db)
    assert len(exp.ancestors) == 4
    assert warning(exp, "external_ref_resolved") == [
        "2 référence(s) Index Mellifera rattachée(s), 4 ancêtre(s) ajouté(s) au pedigree, "
        "3 lien(s) erroné(s) de la base coupé(s)"]
    build_pedigree(exp.queens, 12, extra_queens=exp.ancestors)       # ne lève pas « pedigree cyclique »


# ------------------------------------------------------------------ relecture finale
def test_farm_queen_with_own_drone_queens_mother_never_shares_an_ancestor_sire_node(db):
    """C1 : la reine 7 n'est pas résolue (écart 5) et cite, tel quel, le code du père de 0276 : son nœud sire
    reste distinct de celui des ancêtres, et chacun garde sa propre mère des RPM."""
    qs = [queen(5, 2024, mother_ref_ext="FR-CBP-B-0276-2022"), queen(6, 2020),
          queen(7, 2024, mating=mating("insemination", mother=6, ext="FR-NG-B-0009-2018"))]
    exp = expand_external(inp(qs), db)
    ped, _ = build_pedigree(exp.queens, 12, extra_queens=exp.ancestors)
    a276 = next(a for a in exp.ancestors if a.id == exp.queens[0].mother)
    s_anc, s7 = ped.sire_of_mating(a276.id), ped.sire_of_mating(7)
    assert s_anc != s7
    assert ped.dam[s_anc] == ped.q_node(a276.mating.drone_queens_mother)
    assert ped.dam[s7] == ped.q_node(6)


def test_farm_queen_descending_from_the_lineage_as_drone_queens_mother_is_not_a_cycle(db):
    """C1 : la reine 7 a pour mère des RPM la reine 5, fille de 0276, et cite le code du père de 0276."""
    qs = [queen(5, 2024, mother_ref_ext="FR-CBP-B-0276-2022"),
          queen(7, 2025, mating=mating("insemination", mother=5, ext="FR-NG-B-0009-2018"))]
    exp = expand_external(inp(qs), db)
    build_pedigree(exp.queens, 12, extra_queens=exp.ancestors)      # ne lève pas « pedigree cyclique »


def test_unexpected_value_types_in_the_base_are_cleaned_not_fatal(db):
    """I1 : n_drone_queens hors domaine, mating_type hors enum, race_description non texte."""
    conn = sqlite3.connect(db)
    conn.executescript("""
        INSERT INTO mating_code VALUES ('bad', 'foo', NULL);
        UPDATE queen SET n_drone_queens = -1, race_description = X'00FF' WHERE code = 'FR-CBP-B-0276-2022';
        UPDATE queen SET mating_code = 'bad', n_drone_queens = 'abc' WHERE code = 'FR-CBP-B-0100-2020';
        """)
    conn.commit()
    conn.close()
    exp = expand_external(inp([queen(1, 2024, mother_ref_ext="FR-CBP-B-0276-2022")]), db)
    by_born = {a.born: a for a in exp.ancestors}
    assert by_born[2022].mating.n_drone_queens is None and by_born[2022].variety is None
    assert by_born[2020].mating.type == "unknown" and by_born[2020].mating.n_drone_queens is None
    assert warning(exp, "external_ref_resolved")[0].startswith("1 référence(s)")


def test_any_failure_of_the_resolution_leaves_the_computation_without_ancestors(db, monkeypatch):
    """I1 : la garantie « aucune donnée de la base ne fait échouer le calcul » tient sans condition."""
    def boom(*a, **k):
        raise RuntimeError("inattendu")
    monkeypatch.setattr("blup2.external._expand", boom)
    src = inp([queen(1, 2024, mother_ref_ext="FR-CBP-B-0276-2022")])
    exp = expand_external(src, db)
    assert exp.queens == src.queens and exp.ancestors == []
    assert warning(exp, "external_pedigree_unavailable") == [
        "base Index Mellifera illisible : 1 référence(s) externe(s) non résolue(s)"]


def test_climb_stops_at_the_ancestor_cap(tmp_path, monkeypatch):
    """I2 : le nombre total d'ancêtres est plafonné, la remontée s'arrête et le signale."""
    rows = [(f"FR-L-B-{k:04d}-{1960 + k}", 1960 + k, f"FR-L-B-{k - 1:04d}-{1959 + k}" if k else None, None, None, None)
            for k in range(60)]
    db = make_db(tmp_path / "m.sqlite", rows)
    monkeypatch.setattr("blup2.external.MAX_ANCESTORS", 10)
    exp = expand_external(inp([queen(1, 2024, mother_ref_ext=rows[-1][0])], depth=30), db)
    assert len(exp.ancestors) == 10 and min(a.born for a in exp.ancestors) == 2010
    assert all(a.mother is None or a.mother in {b.id for b in exp.ancestors} for a in exp.ancestors)
    assert warning(exp, "external_ref_resolved") == [
        "1 référence(s) Index Mellifera rattachée(s), 10 ancêtre(s) ajouté(s) au pedigree, "
        "remontée arrêtée au plafond de 10 ancêtres"]
