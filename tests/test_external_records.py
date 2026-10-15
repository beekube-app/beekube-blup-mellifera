"""Volet 3 Index Mellifera : notes des ancêtres résolus et de leurs contemporaines (spec volet 3 §4, §6)."""
import datetime as dt
import sqlite3

import pytest

import blup2.external_records as external_records
from blup2.external import expand_external
from tests.builders import criterion, queen, record
from tests.external_base import NOTES, inp, make_db, warning

D = dt.date


@pytest.fixture
def db(tmp_path):
    return make_db(tmp_path / "mellifera.sqlite", notes=NOTES)


def farm(criteria, apiary=1, evaluator=1, depth=7):
    """Reine 1 fille de FR-CBP-B-0276-2022, notée sur chaque critère. Ancêtres : 0276 = -1, 0100 = -2,
    0009 = -3, 0050 = -4, 0010 = -5 ; contemporaines ensuite, par code croissant."""
    recs = [record(1, 0.5, criterion=c.id, apiary=apiary, evaluator=evaluator) for c in criteria]
    return inp([queen(1, 2024, mother_ref_ext="FR-CBP-B-0276-2022")], records=recs, criteria=criteria, depth=depth)


def notes(exp):
    return sorted((r.queen, r.criterion, round(r.value, 6), r.date, r.apiary, r.evaluator) for r in exp.records)


# ------------------------------------------------------------------ rattachement et nettoyage
def test_notes_of_ancestors_and_contemporaries_are_attached_trait_by_trait(db):
    exp = expand_external(farm([criterion(1, trait="gentleness"), criterion(2, trait="swarming"), criterion(3)]), db)
    r6 = lambda x: round(x, 6)                                    # noqa: E731
    assert notes(exp) == [
        (-9, 1, r6(2 / 6), D(2021, 7, 1), -3, -3),               # FR-T-B-0003-2020, contemporaine de 0100
        (-8, 1, r6(4 / 6), D(2023, 7, 7), -3, -3),               # FR-T-B-0002-2022
        (-7, 1, r6(3 / 6), D(2023, 7, 6), -3, -3),               # FR-T-B-0001-2022
        (-7, 2, r6(4 / 6), D(2023, 7, 6), -3, -3),
        (-6, 1, r6(5 / 6), D(2024, 1, 1), -2, -2),               # FR-CBP-B-0300-2023, radar de l'éleveur 145
        (-2, 1, r6(5 / 6), D(2021, 7, 1), -3, -3),               # 0100 : (6 + 4) / 2, deux modes moyennés
        (-1, 1, r6(3 / 6), D(2024, 1, 1), -2, -2),               # 0276, radar 2024
        (-1, 1, r6(5 / 6), D(2023, 7, 6), -3, -3),               # 0276, saison 2023 (le radar 2023 est ignoré)
        (-1, 2, r6(3 / 6), D(2023, 7, 6), -3, -3),
    ]
    assert all((r.scale_min, r.scale_max, r.hive_class) == (0.0, 1.0, "hive") for r in exp.records)
    assert warning(exp, "external_records_added") == [
        "9 note(s) Index Mellifera ajoutée(s) sur 2 critère(s) : 2 ancêtre(s) et 4 contemporaine(s) ; "
        "écartées : 1 sans fiche de reine, 1 sans contemporaine"]
    assert warning(exp, "external_ref_resolved") == [   # l'ascendance du volet 2 est inchangée
        "1 référence(s) Index Mellifera rattachée(s), 5 ancêtre(s) ajouté(s) au pedigree"]


def test_contemporaries_are_built_like_the_ancestors(db):
    exp = expand_external(farm([criterion(1, trait="gentleness")]), db)
    assert [q.id for q in exp.contemporaries] == [-6, -7, -8, -9]
    assert {q.id: q.born for q in exp.contemporaries} == {-6: 2023, -7: 2022, -8: 2022, -9: 2020}
    assert {q.id: q.mother for q in exp.contemporaries}[-9] == -4      # fille de 0050, ancêtre résolu
    assert all(not q.available_as_mother and not q.available_as_drone_source for q in exp.contemporaries)


def test_each_discard_reason_is_counted(tmp_path):
    rows = [("FR-CBP-B-0276-2022", "season_table", "evjp", "mildness", 5, "2023-07-06", 300),   # gardée
            ("FR-CBP-B-0276-2022", "season_table", "evjp", "hyg24h", 80, "2023-07-06", 300),    # gardée : 0,8
            ("FR-T-B-0001-2022", "season_table", "evjp", "mildness", 3, "2023-07-06", 300),     # contemporaine
            ("FR-T-B-0001-2022", "season_table", "evjp", "hyg24h", 50, "2023-07-06", 300),
            ("FR-CBP-B-0276-2022", "season_table", "evjp", "honeyAmt", 30, "2023-07-06", 300),
            ("FR-CBP-B-0276-2022", "season_table", "evsp", "mildness", 4, "", 300),
            ("FR-CBP-B-0276-2022", "season_table", "evsp", "mildness", 7, "2023-07-07", 300),
            ("FR-CBP-B-0276-2022", "season_table", "evjp", "hyg24h", 150, "2023-07-07", 300),
            ("FR-CBP-B-0276-2022", "season_table", "evjp", "mildness", None, "2023-07-08", 300),
            ("FR-CBP-B-0276-2022", "season_table", "evjp", "mildness", 4, "2023-07-09", None),
            ("FR-CBP-B-0276-2022", "season_table", "evjp", "mildness", 4, "2021-07-09", 300)]
    db = make_db(tmp_path / "m.sqlite", notes=rows)
    cs = [criterion(1, trait="gentleness"), criterion(2, trait="hygienic_24h"), criterion(3, trait="honey_yield")]
    exp = expand_external(farm(cs), db)
    assert notes(exp) == [(-6, 1, 0.5, D(2023, 7, 6), -2, -2), (-6, 2, 0.5, D(2023, 7, 6), -2, -2),
                          (-1, 1, round(5 / 6, 6), D(2023, 7, 6), -2, -2), (-1, 2, 0.8, D(2023, 7, 6), -2, -2)]
    assert warning(exp, "external_records_added") == [
        "4 note(s) Index Mellifera ajoutée(s) sur 2 critère(s) : 1 ancêtre(s) et 1 contemporaine(s) ; écartées : "
        "1 en kg, 1 sans date, 2 hors échelle, 1 non numériques, 1 sans rucher, 1 antérieures à la naissance"]


def test_criterion_direction_opposite_to_the_site_inverts_the_note(db):
    same = {(q, c, d): v for q, c, v, d, *_ in notes(expand_external(farm([criterion(1, trait="gentleness")]), db))}
    flip = {(q, c, d): v for q, c, v, d, *_ in notes(
        expand_external(farm([criterion(1, trait="gentleness", direction=-1)]), db))}
    assert flip.keys() == same.keys() and all(flip[k] == pytest.approx(1 - same[k], abs=1e-6) for k in same)


def test_identifiers_stay_below_those_of_the_farm(db):
    exp = expand_external(farm([criterion(1, trait="gentleness")], apiary=-7, evaluator=-9), db)
    assert sorted({r.apiary for r in exp.records}) == [-11, -10]     # saison du testeur 300, radar de l'éleveur 145
    assert all(r.apiary == r.evaluator for r in exp.records)


def test_contemporaries_are_capped(db, monkeypatch):
    monkeypatch.setattr(external_records, "MAX_CONTEMPORARIES", 2)
    exp = expand_external(farm([criterion(1, trait="gentleness"), criterion(2, trait="swarming")]), db)
    assert [q.id for q in exp.contemporaries] == [-6, -7]             # 0300 et T-0001, par code croissant
    assert warning(exp, "external_records_added") == [
        "6 note(s) Index Mellifera ajoutée(s) sur 2 critère(s) : 1 ancêtre(s) et 2 contemporaine(s) ; écartées : "
        "1 sans fiche de reine, 2 sans contemporaine ; plafond de 2 contemporaines atteint"]


# ------------------------------------------------------------------ rien à ajouter
@pytest.mark.parametrize("crit", [
    criterion(1),                                   # critère sans caractère
    criterion(1, trait="colony_strength"),          # caractère sans nom Index Mellifera (BeeBreed seul)
    criterion(1, trait="fertility"),                # nom Index Mellifera absent de cette base
])
def test_nothing_is_added_without_a_known_trait(db, crit):
    exp = expand_external(farm([crit]), db)
    assert exp.records == [] and exp.contemporaries == [] and warning(exp, "external_records_added") == []
    assert len(exp.ancestors) == 5


def test_criterion_without_farm_records_receives_nothing(db):
    src = inp([queen(1, 2024, mother_ref_ext="FR-CBP-B-0276-2022")], records=[record(1, 0.5, criterion=1)],
              criteria=[criterion(1), criterion(2, trait="gentleness")])
    exp = expand_external(src, db)
    assert exp.records == [] and warning(exp, "external_records_added") == []


def test_depth_bounds_the_ancestors_whose_notes_are_read(db):
    exp = expand_external(farm([criterion(1, trait="gentleness")], depth=1), db)
    assert [a.id for a in exp.ancestors] == [-1]
    assert {r.queen for r in exp.records} == {-1, -2, -3, -4}         # 0276 et ses contemporaines 0300, T-0001, T-0002


# ------------------------------------------------------------------ erreurs
def test_base_without_note_tables_keeps_the_ancestry(tmp_path):
    db = make_db(tmp_path / "m.sqlite")                     # base du volet 2 : ni trait, ni evaluation
    exp = expand_external(farm([criterion(1, trait="gentleness")]), db)
    assert len(exp.ancestors) == 5 and exp.records == [] and exp.contemporaries == []
    assert warning(exp, "external_records_unavailable") == [
        "notes Index Mellifera illisibles : calcul avec l'ascendance mais sans leurs notes"]
    assert warning(exp, "external_pedigree_unavailable") == []


def test_base_without_note_tables_is_silent_when_no_trait_is_declared(tmp_path):
    exp = expand_external(farm([criterion(1)]), make_db(tmp_path / "m.sqlite"))
    assert [w.code for w in exp.warnings] == ["external_ref_resolved"]


def test_unexpected_values_in_the_base_are_discarded(tmp_path):
    db = make_db(tmp_path / "m.sqlite", notes=[
        ("FR-CBP-B-0276-2022", "season_table", "evjp", "mildness", 5, "2023-07-06", 300)])
    conn = sqlite3.connect(db)
    conn.executescript("UPDATE evaluation SET value = 'cinq'; UPDATE trait SET scale_max = 'x' WHERE name = 'swarming';"
                       "INSERT INTO evaluation (queen_name_raw, queen_code, queen_key, tester, tester_breeder_id, "
                       "eval_date, mode, criterion_code, value, value_raw, source, url) VALUES ('a', "
                       "'FR-CBP-B-0276-2022', 'k2', 'T', 300, '2023-13-45', 'evjp', 'swarming', 3, '3', "
                       "'season_table', '/x/'), ('a', 'FR-CBP-B-0276-2022', 'k3', 'T', 300, '2023-07-06', 'evjp', "
                       "'swarming', 3, '3', 'season_table', '/x/'), ('a', 'FR-CBP-B-0276-2022', 'k4', 'T', 'abc', "
                       "'2023-07-06', 'evjp', 'mildness', 3, '3', 'season_table', '/x/');")
    conn.close()
    exp = expand_external(farm([criterion(1, trait="gentleness"), criterion(2, trait="swarming")]), db)
    assert exp.records == []
    assert warning(exp, "external_records_added") == [
        "0 note(s) Index Mellifera ajoutée(s) sur 0 critère(s) : 0 ancêtre(s) et 0 contemporaine(s) ; écartées : "
        "1 sans date, 1 hors échelle, 1 non numériques, 1 sans rucher"]


# ------------------------------------------------------------------ points de relecture
def test_negative_farm_ids_keep_external_ids_below(db):
    src = inp([queen(-20, 2024, mother_ref_ext="FR-CBP-B-0276-2022")],
              records=[record(-20, 0.5, apiary=-30, evaluator=-40)], criteria=[criterion(1, trait="gentleness")])
    exp = expand_external(src, db)
    ext_ids = [q.id for q in exp.ancestors + exp.contemporaries]
    assert len(set(ext_ids)) == len(ext_ids) == 9 and max(ext_ids) < -20
    assert {r.queen for r in exp.records} <= set(ext_ids)
    assert max(r.apiary for r in exp.records) < -40


def test_cycles_and_self_references_between_contemporaries_are_cut(db):
    from blup2.pedigree import find_cycle
    conn = sqlite3.connect(db)
    conn.executescript("UPDATE queen SET mother_code = 'FR-T-B-0002-2022' WHERE code = 'FR-T-B-0001-2022';"
                       "UPDATE queen SET mother_code = 'FR-T-B-0001-2022' WHERE code = 'FR-T-B-0002-2022';"
                       "UPDATE queen SET mother_code = code WHERE code = 'FR-CBP-B-0300-2023';")
    conn.commit()
    conn.close()
    exp = expand_external(farm([criterion(1, trait="gentleness")]), db)
    queens = exp.ancestors + exp.contemporaries
    assert len(exp.contemporaries) == 4 and not find_cycle(queens)
    assert all(q.mother != q.id for q in queens)
    assert {q.id: q.mother for q in exp.contemporaries}[-9] == -4      # lien valide conservé


def test_depth_zero_disables_the_notes_too(db):
    exp = expand_external(farm([criterion(1, trait="gentleness")], depth=0), db)
    assert exp.ancestors == [] and exp.contemporaries == [] and exp.records == [] and exp.warnings == []
