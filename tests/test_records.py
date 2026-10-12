import numpy as np
import pytest

from blup2.records import age_class, prepare_criterion
from blup2.schemas import BlupV2Input
from tests.builders import criterion, queen, record


def _inp(records, crit=None, queens=None):
    return BlupV2Input(exploitation=1, criteria=[crit or criterion()],
                       queens=queens or [queen(1, 2022), queen(2, 2022), queen(3, 2023)], records=records)


def test_scale_is_normalized_to_unit_interval():
    # 3 sur [1, 5] -> (3 − 1) / 4 = 0,5 ; 4 sur [0, 4] -> 1,0
    d = prepare_criterion(_inp([record(1, 3, scale_min=1, scale_max=5), record(2, 4, scale_min=0, scale_max=4)]), 1)
    assert d.y.tolist() == pytest.approx([0.5, 1.0])
    assert d.scale_mixed is True


def test_percentage_values_are_kept_as_is():
    d = prepare_criterion(_inp([record(1, 0.3)], crit=criterion(kind="percentage")), 1)
    assert d.y.tolist() == pytest.approx([0.3])
    assert d.scale_mixed is False


def test_records_of_a_same_visit_are_averaged():
    # visite 7 : (0,2 + 0,6) / 2 = 0,4 ; la note sans visite est groupée par date
    d = prepare_criterion(_inp([record(1, 0.2, visit=7, evaluator=1), record(1, 0.6, visit=7, evaluator=2),
                                record(1, 0.9, date="2024-06-01")]), 1)
    assert sorted(d.y.tolist()) == pytest.approx([0.4, 0.9])
    assert d.n_records[1] == 3 and d.n_evaluators[1] == 2 and d.n_raw_records == 3


def test_same_visit_in_two_apiaries_stays_two_observations():
    # même identifiant de visite mais ruchers différents : 2 observations, 2 groupes de contemporaines
    recs = [record(1, 0.2, visit=7, apiary=1), record(1, 0.6, visit=7, apiary=2), record(2, 0.5, apiary=1)]
    d = prepare_criterion(_inp(recs), 1)
    assert d.n_obs == 3 and set(d.cg_labels) == {(1, 2024), (2, 2024)}


def test_same_date_without_visit_is_one_observation():
    d = prepare_criterion(_inp([record(1, 0.2), record(1, 0.4)]), 1)
    assert d.n_obs == 1 and d.y[0] == pytest.approx(0.3)


def test_contemporary_group_is_apiary_by_year():
    recs = [record(1, 0.1, apiary=1, date="2023-05-01"), record(2, 0.5, apiary=1, date="2024-05-01"),
            record(3, 0.9, apiary=2, date="2024-05-01")]
    d = prepare_criterion(_inp(recs), 1)
    assert set(d.cg_labels) == {(1, 2023), (1, 2024), (2, 2024)}
    assert d.n_contemporary_groups == 3
    assert "cg" in d.fixed_factors
    assert d.X.shape == (3, 3)          # intercept + 2 contrastes de groupe


def test_single_level_factors_are_dropped():
    d = prepare_criterion(_inp([record(1, 0.1), record(2, 0.5)]), 1)
    assert d.fixed_factors == [] and d.x_names == ["intercept"]


def test_evaluator_hive_class_and_age_factors():
    # 6 colonies, plan croisé : intercept + 2 évaluateurs + 1 classe + 1 âge = 5 colonnes indépendantes
    qs = [queen(k, 2022 if k <= 3 else 2023) for k in range(1, 7)]
    ev = [1, 2, None, 1, 2, None]
    hc = ["hive", "nuc", "hive", "nuc", "hive", "hive"]
    recs = [record(k, 0.1 * k, evaluator=ev[k - 1], hive_class=hc[k - 1]) for k in range(1, 7)]
    d = prepare_criterion(_inp(recs, queens=qs), 1)
    assert d.fixed_factors == ["evaluator", "hive_class", "age"]
    assert d.factor_labels["evaluator"] == [1, 2, -1, 1, 2, -1]   # évaluateur absent = niveau −1
    assert d.factor_labels["age"] == ["2+", "2+", "2+", "1", "1", "1"]
    assert d.X.shape == (6, 5)


def test_dependent_columns_are_removed():
    # évaluateur et rucher confondus : 1 seule colonne de contraste conservée, X de plein rang
    recs = [record(1, 0.1, evaluator=1, apiary=1), record(2, 0.5, evaluator=2, apiary=2)]
    d = prepare_criterion(_inp(recs), 1)
    X = d.X.toarray()
    assert X.shape[1] == 2 and np.linalg.matrix_rank(X) == 2


def test_age_classes():
    assert [age_class(2024, b) for b in (2024, 2023, 2022, 2019, None)] == ["0", "1", "2+", "2+", "?"]


def test_worker_inbreeding_covariate_only_when_it_varies():
    qs = [queen(k, 2022) for k in range(1, 6)]
    recs = [record(k, 0.1 * k) for k in range(1, 6)]
    d = prepare_criterion(_inp(recs, queens=qs), 1, f_workers={k: 0.0 for k in range(1, 6)})
    assert "f_workers" not in d.fixed_factors
    d = prepare_criterion(_inp(recs, queens=qs), 1, f_workers={1: 0.0, 2: 0.125, 3: 0.0625, 4: 0.0, 5: 0.25})
    assert "f_workers" in d.fixed_factors and d.x_names == ["intercept", "f_workers"]


def test_factor_whose_columns_are_all_dependent_is_not_reported():
    # 2 observations, rucher et évaluateur confondus : seul le premier facteur garde sa colonne
    recs = [record(1, 0.1, evaluator=1, apiary=1), record(2, 0.5, evaluator=2, apiary=2)]
    d = prepare_criterion(_inp(recs), 1)
    assert len(d.fixed_factors) == 1 and set(d.factor_labels) == set(d.fixed_factors)


def test_criterion_without_records():
    d = prepare_criterion(_inp([]), 1)
    assert d.n_obs == 0 and d.recorded_queens == [] and d.n_raw_records == 0


def test_only_records_of_the_requested_criterion():
    inp = BlupV2Input(exploitation=1, criteria=[criterion(1), criterion(2)], queens=[queen(1), queen(2)],
                      records=[record(1, 0.1, criterion=1), record(2, 0.5, criterion=2)])
    d = prepare_criterion(inp, 2)
    assert d.recorded_queens == [2] and d.y.tolist() == pytest.approx([0.5])


def test_none_evaluator_is_not_counted():
    d = prepare_criterion(_inp([record(1, 0.2, evaluator=None), record(1, 0.4, date="2024-06-01", evaluator=5)]), 1)
    assert d.n_evaluators[1] == 1
