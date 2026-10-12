"""Non-régression sur les données d'exemple (data/input-example.json converties au format v2)."""
import json
import os
import time

import pytest

from blup2.engine import run_blup
from blup2.schemas import BlupV2Input
from scripts.convert_v1_example import convert
from scripts.snapshot_regression import snapshot

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIX = os.path.join(ROOT, "tests", "fixtures")


def _load(name):
    with open(os.path.join(FIX, name), encoding="utf-8") as f:
        return json.load(f)


@pytest.fixture(scope="module")
def example():
    inp = BlupV2Input.model_validate(_load("input-example-v2.json"))
    t0 = time.perf_counter()
    out = run_blup(inp)
    return inp, out, time.perf_counter() - t0


def test_fixture_is_the_conversion_of_the_v1_example():
    with open(os.path.join(ROOT, "data", "input-example.json"), encoding="utf-8") as f:
        assert convert(json.load(f)) == _load("input-example-v2.json")


def test_example_runs_quickly_with_expected_structure(example):
    inp, out, elapsed = example
    assert elapsed < 10.0
    assert len(out.queens) == 328 and len(inp.records) == 1111
    assert out.base.birth_years == [2021, 2022, 2023, 2024]
    st = [q.indexes["44"].status for q in out.queens]
    # docs/analyses : 69 reines notées et 180 reines informées par le seul pedigree sur le critère 44
    assert (st.count("evaluated"), st.count("pedigree"), st.count("none")) == (69, 180, 79)
    p = out.parameters
    assert (p["44"].n_queens, p["44"].n_records, p["44"].n_contemporary_groups) == (69, 347, 6)
    assert (p["46"].n_queens, p["47"].n_queens, p["274"].n_queens) == (65, 3, 0)
    assert p["44"].source == "prior" and p["44"].reliability_method == "exact"   # 69 reines < 100
    assert (out.connectivity.n_components, out.connectivity.minor_share) == (3, 0.145)
    codes = {w.code for w in out.warnings}
    assert {"unknown_parent_added", "scale_mixed", "few_records", "disconnected_groups", "no_variation"} <= codes
    assert [w.criterion for w in out.warnings if w.code == "no_variation"] == [47]   # 3 notes, toutes à 1


def test_example_matches_snapshot(example):
    _, out, _ = example
    got, ref = snapshot(out), _load("expected-example-v2.json")
    assert got["parameters"] == ref["parameters"] and got["warnings"] == ref["warnings"]
    assert got["n_components"] == ref["n_components"]
    for qid, r in ref["queens"].items():
        g = got["queens"][qid]
        assert g["inbreeding"] == pytest.approx(r["inbreeding"], abs=1e-4)
        assert (g["total_index"] is None) == (r["total_index"] is None)
        if r["total_index"] is not None:
            assert g["total_index"] == pytest.approx(r["total_index"], abs=0.11)
        for c, (status, index, rel) in r["indexes"].items():
            gs, gi, gr = g["indexes"][c]
            assert gs == status
            assert (gi is None and index is None) or gi == pytest.approx(index, abs=0.11)
            assert (gr is None and rel is None) or gr == pytest.approx(rel, abs=0.002)
