import pytest

from blup2 import ENGINE_VERSION
from blup2.engine import run_blup, run_matings
from blup2.schemas import BlupV2Input, MatingV2Input, OptionsIn
from tests.builders import criterion, mating, queen, record


def _inp(**kw):
    qs = [queen(1, 2020, mating=mating("station", station=1, n_rpm=8)), queen(6, 2022), queen(7, 2022)]
    qs += [queen(k, 2022, mother=1) for k in range(2, 6)]
    vals = {2: 0.8, 3: 0.6, 4: 0.7, 5: 0.9, 6: 0.3, 7: 0.5}
    recs = [record(k, v) for k, v in vals.items()] + [record(k, v, criterion=2) for k, v in vals.items()]
    base = dict(exploitation=1, criteria=[criterion(1, weight=2), criterion(2, direction=-1)], queens=qs,
                records=recs, options=OptionsIn(n_suggestions=2))
    base.update(kw)
    return BlupV2Input(**base)


def test_run_blup_output_shape():
    out = run_blup(_inp())
    assert out.status == "success" and out.engine_version == ENGINE_VERSION == "2.1.0"
    assert set(out.parameters) == {"1", "2"} and [q.id for q in out.queens] == [1, 6, 7, 2, 3, 4, 5]
    assert set(out.queens[0].indexes) == {"1", "2"}
    assert len(out.matings) == 7 * 2 and set(out.trend) == {"1", "2"}
    assert out.base.birth_years == [2020, 2022]


def test_opposite_directions_give_mirrored_indexes():
    out = run_blup(_inp())
    for q in out.queens:
        i1, i2 = q.indexes["1"].index, q.indexes["2"].index
        if i1 is not None:
            assert i1 + i2 == pytest.approx(200.0, abs=0.11)


def test_failed_criterion_does_not_stop_the_others(monkeypatch):
    import blup2.engine as engine

    real = engine.fit_standard

    def flaky(cdata, ped, Ainv, params):
        if cdata.criterion.id == 2:
            raise RuntimeError("boum")
        return real(cdata, ped, Ainv, params)

    monkeypatch.setattr(engine, "fit_standard", flaky)
    out = run_blup(_inp())
    assert set(out.parameters) == {"1"}
    assert all(q.indexes["2"].status == "failed" and q.indexes["2"].index is None for q in out.queens)
    assert any(w.code == "criterion_failed" and w.criterion == 2 and "boum" in w.message for w in out.warnings)
    assert out.queens[3].total_index == out.queens[3].indexes["1"].index


def _break_import(monkeypatch, module_name):
    import builtins

    real = builtins.__import__

    def fake(name, globals=None, locals=None, fromlist=(), level=0):
        if name.endswith("model_breeder") or (level > 0 and name == "" and "model_breeder" in (fromlist or ())):
            raise ModuleNotFoundError(f"No module named {module_name!r}", name=module_name)
        return real(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", fake)


def test_breeder_profile_falls_back_to_standard_when_unavailable(monkeypatch):
    _break_import(monkeypatch, "blup2.model_breeder")
    out = run_blup(_inp(profile="breeder"))
    assert out.profile == "breeder"
    # avertissement global du moteur (sans critère, « non disponible »), pas le repli par critère
    # de fit_breeder (< 50 % d'accouplements connus) qui donnerait aussi profile_used = "standard"
    fallback = [w for w in out.warnings if w.code == "breeder_fallback"]
    assert len(fallback) == 1 and fallback[0].criterion is None and "non disponible" in fallback[0].message
    assert all(p.profile_used == "standard" for p in out.parameters.values())


def test_breeder_import_error_of_another_module_is_raised(monkeypatch):
    _break_import(monkeypatch, "scipy_missing")
    with pytest.raises(ModuleNotFoundError):
        run_blup(_inp(profile="breeder"))


def test_unknown_parent_and_disconnected_warnings():
    qs = [queen(1, 2022, mother=99), queen(2, 2022)]
    recs = [record(1, 0.5, apiary=1), record(2, 0.5, apiary=2)]
    out = run_blup(BlupV2Input(exploitation=1, criteria=[criterion()], queens=qs, records=recs))
    codes = [w.code for w in out.warnings]
    assert "unknown_parent_added" in codes and "disconnected_groups" in codes
    assert out.connectivity.n_components == 2


def test_disconnected_warning_needs_ten_percent_of_recorded_queens_outside_the_main_block():
    def run(n_minor):
        # 100 reines de base non apparentées notées : n_minor au rucher 2, les autres au rucher 1
        qs = [queen(k, 2022) for k in range(1, 101)]
        recs = [record(k, (k % 7) / 10, apiary=2 if k <= n_minor else 1) for k in range(1, 101)]
        return run_blup(BlupV2Input(exploitation=1, criteria=[criterion()], queens=qs, records=recs,
                                    options=OptionsIn(n_suggestions=0)))

    out = run(1)
    assert (out.connectivity.n_components, out.connectivity.minor_share) == (2, 0.01)
    assert all(w.code != "disconnected_groups" for w in out.warnings)
    w = [w for w in run(30).warnings if w.code == "disconnected_groups"]
    assert len(w) == 1 and w[0].criterion is None and "30 % des reines notées" in w[0].message


def test_run_matings_evaluates_requested_pairs():
    data = _inp().model_dump()
    data["pairs"] = [{"mother": 2, "drone_source": 3}, {"mother": 6, "drone_source": 7}]
    out = run_matings(MatingV2Input.model_validate(data))
    assert [(m.mother, m.drone_source) for m in out.matings] == [(2, 3), (6, 7)]
    # 2 et 3 : filles de 1 × même station (S = 8, D = 12) -> a = ¼ + 1/24 + 1/32 ; F = ½ a
    assert out.matings[0].expected_inbreeding == pytest.approx(round(0.5 * (0.25 + 1 / 24 + 1 / 32), 4))
    assert out.matings[1].expected_inbreeding == 0.0
