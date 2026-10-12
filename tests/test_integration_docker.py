"""Non-régression de bout en bout : image Docker, gunicorn, vrai Rscript, jeu d'exemple.

Lancement : .venv/bin/pytest -m docker -v
Régénérer l'instantané (après un changement voulu de blup.r) : UPDATE_SNAPSHOT=1 .venv/bin/pytest -m docker -v
Attention : la fixture lance `docker compose up -d --build api` puis `docker compose down`.
"""
import json
import math
import os
import subprocess
import time
import urllib.error
import urllib.request

import pytest

pytestmark = pytest.mark.docker

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
API_URL = os.getenv("BLUP_API_URL", "http://localhost:8081")
EXAMPLE = os.path.join(ROOT, "data", "input-example.json")
SNAPSHOT = os.path.join(ROOT, "tests", "fixtures", "v1-example-output.json")
CRITERIA = ["46", "139", "44", "47", "274", "228", "45"]
METHODS = {"BLUP", "lm", "no_data", "failed"}


def post_json(url, payload, timeout=200):
    req = urllib.request.Request(url, data=json.dumps(payload).encode("utf-8"),
                                 headers={"Content-Type": "application/json"}, method="POST")
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.status, json.load(resp)
    except urllib.error.HTTPError as e:
        return e.code, json.load(e)


@pytest.fixture(scope="module")
def api_url():
    subprocess.run(["docker", "compose", "up", "-d", "--build", "api"], cwd=ROOT, check=True)
    try:
        deadline = time.time() + 180
        while True:
            try:
                with urllib.request.urlopen(API_URL + "/", timeout=5) as resp:
                    if resp.status == 200:
                        break
            except OSError:
                pass
            if time.time() > deadline:
                pytest.fail("API non joignable après 180 s")
            time.sleep(2)
        yield API_URL
    finally:
        subprocess.run(["docker", "compose", "down"], cwd=ROOT, check=False)


@pytest.fixture(scope="module")
def example_results(api_url):
    with open(EXAMPLE, encoding="utf-8") as f:
        payload = json.load(f)
    status, body = post_json(api_url + "/blup", payload)
    assert status == 200, body
    assert body["status"] == "success"
    return body["results"]


def has_index(queen):
    return any(v is not None for v in queen["blups"].values())


def test_every_queen_has_every_criterion(example_results):
    assert len(example_results["blup"]) == 328
    for queen in example_results["blup"]:
        assert set(queen["blups"]) == set(CRITERIA)
        assert set(queen["methods"]) == set(CRITERIA)
        assert set(queen["methods"].values()) <= METHODS
        for c in CRITERIA:
            assert (queen["blups"][c] is not None) == (queen["methods"][c] in {"BLUP", "lm"})


def test_criterion_with_one_note_per_queen_is_not_lost(example_results):
    rated = [q for q in example_results["blup"] if q["methods"]["47"] != "no_data"]
    assert len(rated) == 3
    assert all(q["methods"]["47"] in {"BLUP", "lm", "failed"} for q in rated)


def test_criteria_without_notes_are_no_data(example_results):
    for c in ("228", "274"):
        assert all(q["methods"][c] == "no_data" and q["blups"][c] is None
                   for q in example_results["blup"])


def test_heritabilities_cover_all_criteria_without_fake_se(example_results):
    assert set(example_results["heritabilities"]) == set(CRITERIA)
    assert all(h["se"] is None for h in example_results["heritabilities"].values())


def test_mating_candidates_use_any_criterion(example_results):
    queens = example_results["blup"]
    indexed = {q["queenbee"] for q in queens if has_index(q) and not q["died"]}
    first_criterion = {q["queenbee"] for q in queens if q["blups"]["46"] is not None}
    with_suggestions = {q["queenbee"] for q in queens if q["mating_suggestions"]}
    assert len(indexed) == 69 and len(first_criterion) == 65
    assert with_suggestions == indexed
    for q in queens:
        for s in q["mating_suggestions"]:
            assert s["queen"] in indexed and s["queen"] != q["queenbee"]


def assert_close(actual, expected, path="$"):
    if isinstance(expected, dict):
        assert isinstance(actual, dict) and set(actual) == set(expected), path
        for key in expected:
            assert_close(actual[key], expected[key], f"{path}.{key}")
    elif isinstance(expected, list):
        assert isinstance(actual, list) and len(actual) == len(expected), path
        for i, (a, e) in enumerate(zip(actual, expected)):
            assert_close(a, e, f"{path}[{i}]")
    elif isinstance(expected, float) or isinstance(actual, float):
        assert actual is not None and math.isclose(actual, expected, rel_tol=1e-6, abs_tol=2e-4), \
            f"{path}: {actual} != {expected}"
    else:
        assert actual == expected, f"{path}: {actual!r} != {expected!r}"


def test_matches_snapshot(example_results):
    if os.getenv("UPDATE_SNAPSHOT") == "1":
        os.makedirs(os.path.dirname(SNAPSHOT), exist_ok=True)
        with open(SNAPSHOT, "w", encoding="utf-8") as f:
            json.dump(example_results, f, indent=2, sort_keys=True)
            f.write("\n")
        pytest.skip("instantané régénéré")
    with open(SNAPSHOT, encoding="utf-8") as f:
        expected = json.load(f)
    assert_close(example_results, expected)
