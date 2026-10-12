import json
import os

from app import BLUPInput

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLE = os.path.join(ROOT, "data", "input-example.json")


def load_example():
    with open(EXAMPLE, encoding="utf-8") as f:
        return json.load(f)


def count_rated(data, criterion):
    return sum(1 for q in data["data"] if q["evaluate"][criterion])


def test_example_is_valid_api_input():
    BLUPInput.model_validate(load_example())


def test_example_is_anonymized():
    data = load_example()
    assert data["exploitation"] == 1
    assert max(q["queenbee"] for q in data["data"]) <= 400
    users = {r["user"] for q in data["data"] for recs in q["evaluate"].values() for r in recs}
    assert users <= set(range(1, 11))


def test_example_keeps_reference_cases():
    data = load_example()
    assert len(data["data"]) == 328
    assert data["evaluate"] == [46, 139, 44, 47, 274, 228, 45]
    assert count_rated(data, "47") == 3
    assert count_rated(data, "228") == 0 and count_rated(data, "274") == 0
    assert count_rated(data, "46") == 65
    assert all(q["evaluate"]["196"] == [] for q in data["data"])
