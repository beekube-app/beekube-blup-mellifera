"""Profil breeder de bout en bout : moteur, API, contrat OpenAPI inchangé."""
import json

import numpy as np
import pytest

from app import app
from blup2.engine import run_blup
from blup2.schemas import BlupV2Input
from tests.simdata import simulate_breeder_scenario


def _scenario(**kw):
    inp, _, recorded, _ = simulate_breeder_scenario(1, **kw)
    inp["options"]["n_suggestions"] = 2
    return inp, recorded


def test_breeder_profile_end_to_end():
    inp, recorded = _scenario()
    out = run_blup(BlupV2Input.model_validate(inp))
    p = out.parameters["1"]
    assert out.profile == "breeder" and p.profile_used == "breeder" and p.source == "input"
    assert p.sigma2_aq == pytest.approx(0.5 * p.sigma2_a) and p.rg == pytest.approx(-0.5)
    assert all(w.code != "breeder_fallback" for w in out.warnings)
    for q in out.queens:
        ci = q.indexes["1"]
        assert (ci.status == "evaluated") == (q.id in recorded)
        assert ci.index is not None and ci.index_queen is not None and ci.index_worker is not None
    assert len(out.matings) > 0


def test_breeder_and_standard_profiles_differ_but_agree_on_ranking():
    inp, _ = _scenario()
    b = run_blup(BlupV2Input.model_validate(inp))
    inp["profile"] = "standard"
    s = run_blup(BlupV2Input.model_validate(inp))
    assert s.parameters["1"].profile_used == "standard" and s.parameters["1"].sigma2_aq is None
    ib = [q.indexes["1"].index for q in b.queens]
    is_ = [q.indexes["1"].index for q in s.queens]
    assert ib != is_
    assert np.corrcoef(ib, is_)[0, 1] > 0.8
    assert all(q.indexes["1"].index_worker is None for q in s.queens)


def test_breeder_fallback_is_reported_per_criterion():
    inp, _ = _scenario()
    for q in inp["queens"]:
        q["mating"] = None                       # aucun accouplement connu
    out = run_blup(BlupV2Input.model_validate(inp))
    assert out.parameters["1"].profile_used == "standard"
    w = [w for w in out.warnings if w.code == "breeder_fallback"]
    assert len(w) == 1 and w[0].criterion == 1 and w[0].message.startswith("0 %")


def test_breeder_request_on_current_beekube_data_falls_back_to_standard():
    # cas courant (validation sur données réelles : repli sur 169 critères sur 170) : fécondations
    # non saisies pour deux reines sur trois. Le profil breeder rend alors les index du profil
    # standard, avec breeder_fallback et la part d'accouplements connus du critère.
    inp, _ = _scenario()
    for q in inp["queens"]:
        if q["id"] % 3:
            q["mating"] = None
    inp["records"][0]["scale_max"] = 40.0        # deuxième échelle : le critère émet aussi scale_mixed
    b = run_blup(BlupV2Input.model_validate(inp))
    inp["profile"] = "standard"
    s = run_blup(BlupV2Input.model_validate(inp))
    p = b.parameters["1"]
    assert b.profile == "breeder" and p.profile_used == "standard" and p.sigma2_aq is None
    assert 0.0 < p.known_mating_share < 0.5 and s.parameters["1"].known_mating_share is None
    w = [w for w in b.warnings if w.code == "breeder_fallback"]
    assert len(w) == 1 and w[0].criterion == 1 and w[0].message.startswith(f"{round(100 * p.known_mating_share)} %")
    assert [w.code for w in b.warnings if w.criterion == 1] == ["breeder_fallback", "scale_mixed"]   # en tête
    assert [w.code for w in b.warnings if w.code != "breeder_fallback"] == [w.code for w in s.warnings]
    assert [w.code for w in s.warnings] == ["scale_mixed"]
    assert len(b.queens) == len(s.queens) > 0
    for qb, qs in zip(b.queens, s.queens, strict=True):
        assert qb.indexes["1"] == qs.indexes["1"] and qb.indexes["1"].index_worker is None
        assert qb.total_index == qs.total_index
    assert b.matings == s.matings


def test_v2_blup_breeder_over_http():
    inp, _ = _scenario()
    with app.test_client() as c:
        r = c.post("/v2/blup", json=inp)
    assert r.status_code == 200
    body = r.get_json()
    assert body["profile"] == "breeder" and body["parameters"]["1"]["profile_used"] == "breeder"
    assert body["queens"][0]["indexes"]["1"]["index_worker"] is not None
    json.dumps(body)
