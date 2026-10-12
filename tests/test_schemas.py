import copy

import pytest
from pydantic import ValidationError

from blup2.schemas import BlupV2Input, MatingV2Input

VALID = {
    "exploitation": 798,
    "profile": "standard",
    "options": {"base_years": 5, "estimate_variances": "auto", "default_n_drones": 12,
                "inbreeding_penalty": 100, "max_inbreeding": 0.0625, "n_suggestions": 10},
    "criteria": [{"id": 44, "name": "Douceur", "direction": 1, "kind": "note", "weight": 3,
                  "elimination": False, "h2": None, "repeatability": None, "queen_ratio": None, "rg": None}],
    "queens": [
        {"id": 4074, "born": 2021, "mother": 3257, "mother_ref_ext": None,
         "mating": {"type": "insemination", "drone_queens": [3301], "drone_queens_mother": None,
                    "n_drone_queens": None, "n_drones": 8, "external_ref": None, "station": None, "year": 2021},
         "variety": "Buckfast", "died": False, "eliminated": False,
         "available_as_mother": True, "available_as_drone_source": True}
    ],
    "records": [
        {"queen": 4074, "criterion": 44, "value": 3, "scale_min": 1, "scale_max": 5, "date": "2024-03-15",
         "evaluator": 848, "apiary": 1504, "hive_type": "dadant-10", "hive_class": "hive", "visit": 99812}
    ],
}


def mutate(path, value):
    d = copy.deepcopy(VALID)
    target = d
    for k in path[:-1]:
        target = target[k]
    target[path[-1]] = value
    return d


def assert_invalid(d, fragment):
    with pytest.raises(ValidationError) as exc:
        BlupV2Input.model_validate(d)
    assert fragment in str(exc.value)


def test_spec_example_is_valid_and_absent_parents_are_accepted():
    # 3257 (mère) et 3301 (RPM) sont absentes de queens : acceptées (ajoutées en base au pedigree)
    inp = BlupV2Input.model_validate(VALID)
    assert inp.queens[0].mating.n_drones == 8
    assert inp.records[0].date.year == 2024


def test_defaults_when_options_omitted():
    d = copy.deepcopy(VALID)
    del d["options"]
    inp = BlupV2Input.model_validate(d)
    assert (inp.options.base_years, inp.options.default_n_drones, inp.options.n_suggestions) == (5, 12, 10)
    assert inp.options.inbreeding_penalty == 100 and inp.options.max_inbreeding == 0.0625


@pytest.mark.parametrize("path,value,fragment", [
    (("profile",), "expert", "profile"),
    (("queens", 0, "mating", "type"), "artificial", "type"),
    (("criteria", 0, "kind"), "binary", "kind"),
    (("records", 0, "hive_class"), "ruchette", "hive_class"),
    (("criteria", 0, "direction"), 0, "direction"),
    (("queens", 0, "born"), 21, "born"),
    (("queens", 0, "born"), "2021", "born"),
    (("records", 0, "date"), "15/03/2024", "date"),
    (("records", 0, "date"), "2024-3-15", "date"),
    (("records", 0, "value"), 6, "hors de l'échelle"),
    (("records", 0, "scale_max"), 1, "scale_max"),
    (("options", "estimate_variances"), "sometimes", "estimate_variances"),
])
def test_enumerations_and_formats_rejected(path, value, fragment):
    assert_invalid(mutate(path, value), fragment)


def test_unknown_field_rejected():
    assert_invalid(mutate(("queens", 0, "colour"), "red"), "colour")


def test_record_queen_and_criterion_must_exist():
    assert_invalid(mutate(("records", 0, "queen"), 999), "records[].queen absent de queens : [999]")
    assert_invalid(mutate(("records", 0, "criterion"), 7), "records[].criterion absent de criteria : [7]")


def test_duplicate_ids_rejected():
    d = copy.deepcopy(VALID)
    d["queens"].append(copy.deepcopy(d["queens"][0]))
    assert_invalid(d, "identifiants de reines dupliqués : [4074]")


def test_percentage_requires_unit_scale():
    d = mutate(("criteria", 0, "kind"), "percentage")
    assert_invalid(d, "scale_min = 0 et scale_max = 1 attendus")
    d["records"][0].update(value=0.4, scale_min=0, scale_max=1)
    assert BlupV2Input.model_validate(d).records[0].value == 0.4


def test_repeatability_lower_than_h2_rejected():
    d = mutate(("criteria", 0, "h2"), 0.5)
    d["criteria"][0]["repeatability"] = 0.3
    assert_invalid(d, "repeatability (0.3) < h2 (0.5)")


def test_self_references_are_accepted():
    # reine qui est sa propre mère, sa propre RPM et la mère de ses RPM : pas de 422 (spec §5.1), le
    # pedigree retire ces liens et avertit (self_reference_removed, tâche 2)
    d = mutate(("queens", 0, "mother"), 4074)
    d["queens"][0]["mating"].update(drone_queens=[4074], drone_queens_mother=4074)
    q = BlupV2Input.model_validate(d).queens[0]
    assert q.mother == 4074 and q.mating.drone_queens == [4074]


def test_born_between_1900_and_next_year():
    import datetime as dt

    last = dt.date.today().year + 1
    for year in (1900, last):
        assert BlupV2Input.model_validate(mutate(("queens", 0, "born"), year)).queens[0].born == year
    for year in (1899, last + 1):
        assert_invalid(mutate(("queens", 0, "born"), year), f"born ({year}) doit être compris entre 1900 et {last}")


def test_mating_pairs_must_reference_known_queens():
    d = copy.deepcopy(VALID)
    d["pairs"] = [{"mother": 4074, "drone_source": 5}]
    with pytest.raises(ValidationError) as exc:
        MatingV2Input.model_validate(d)
    assert "pairs : reines absentes de queens : [5]" in str(exc.value)


def test_mating_pair_with_itself_rejected():
    d = copy.deepcopy(VALID)
    d["pairs"] = [{"mother": 4074, "drone_source": 4074}]
    with pytest.raises(ValidationError) as exc:
        MatingV2Input.model_validate(d)
    assert "ses propres mâles : [4074]" in str(exc.value)


def test_date_objects_are_accepted_but_json_dates_must_be_iso_strings():
    import datetime as dt

    inp = BlupV2Input.model_validate(VALID)
    again = BlupV2Input.model_validate(inp.model_dump())          # date déjà validée (dt.date)
    assert again.records[0].date == dt.date(2024, 3, 15)
    assert_invalid(mutate(("records", 0, "date"), 20240315), "date")
    assert_invalid(mutate(("records", 0, "date"), "2024-03-15T10:00"), "date")
