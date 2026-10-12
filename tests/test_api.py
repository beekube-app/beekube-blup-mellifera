import copy
import json

import pytest

from app import app
from tests.test_schemas import VALID


@pytest.fixture()
def client():
    app.config["TESTING"] = True
    with app.test_client() as c:
        yield c


def test_v2_blup_returns_engine_output(client):
    r = client.post("/v2/blup", json=VALID)
    assert r.status_code == 200
    body = r.get_json()
    assert body["status"] == "success" and body["engine_version"] == "2.0.0"
    assert [q["id"] for q in body["queens"]] == [4074]
    assert set(body) == {"status", "engine_version", "profile", "parameters", "base", "queens", "matings",
                         "trend", "connectivity", "warnings"}


def test_v2_blup_validation_error_is_422_v2_format(client):
    bad = copy.deepcopy(VALID)
    bad["records"][0]["value"] = 9          # hors de [1, 5]
    r = client.post("/v2/blup", json=bad)
    assert r.status_code == 422
    body = r.get_json()
    assert body["status"] == "error" and "hors de l'échelle" in json.dumps(body["errors"], ensure_ascii=False)


def test_v2_mating_evaluates_pairs(client):
    data = copy.deepcopy(VALID)
    data["queens"].append({"id": 5000, "born": 2022})
    data["pairs"] = [{"mother": 4074, "drone_source": 5000}]
    r = client.post("/v2/mating", json=data)
    assert r.status_code == 200
    m = r.get_json()["matings"]
    assert [(x["mother"], x["drone_source"], x["expected_inbreeding"]) for x in m] == [(4074, 5000, 0.0)]


def test_v2_mating_unknown_pair_queen_is_422(client):
    data = copy.deepcopy(VALID)
    data["pairs"] = [{"mother": 4074, "drone_source": 1}]
    r = client.post("/v2/mating", json=data)
    assert r.status_code == 422 and r.get_json()["status"] == "error"


def test_v1_validation_error_format_is_unchanged(client):
    r = client.post("/blup", json={"exploitation": "x"})
    assert r.status_code == 422
    assert isinstance(r.get_json(), list)        # format flask-openapi3 par défaut (liste d'erreurs)


def test_openapi_documents_v1_and_v2_routes():
    doc = json.loads(json.dumps(app.api_doc))
    paths = doc["paths"]
    assert {"/blup", "/v2/blup", "/v2/mating"} <= set(paths)
    for route in ("/v2/blup", "/v2/mating"):
        responses = paths[route]["post"]["responses"]
        assert responses["422"]["content"]["application/json"]["schema"] == {"$ref": "#/components/schemas/ErrorV2Out"}
        assert "500" in responses
    assert doc["info"]["version"] == "2.0.0"
    assert paths["/blup"]["post"]["operationId"] == "calulate_blup_post"     # SDK v1 inchangé
    # Contrat avec le SDK PHP (plan 04, tâche A1) : une seule classe BlupApi, méthodes nommées
    assert paths["/v2/blup"]["post"]["operationId"] == "blup_v2_post"
    assert paths["/v2/mating"]["post"]["operationId"] == "mating_v2_post"
    for route in ("/blup", "/v2/blup", "/v2/mating"):
        assert paths[route]["post"]["tags"] == ["blup"]


def test_v2_validation_error_does_not_echo_the_input(client):
    cyc = copy.deepcopy(VALID)
    cyc["queens"] = [{"id": 1, "born": 2020, "mother": 2}, {"id": 2, "born": 2020, "mother": 1}]
    cyc["records"] = []
    r = client.post("/v2/blup", json=cyc)
    assert r.status_code == 422
    errors = r.get_json()["errors"]
    assert errors and all("input" not in e for e in errors)


@pytest.mark.parametrize("path,runner", [("/v2/blup", "run_blup"), ("/v2/mating", "run_matings")])
def test_v2_unexpected_exception_is_500_in_error_format(client, monkeypatch, caplog, path, runner):
    import app as app_module

    def boom(_):
        raise RuntimeError("boum")

    monkeypatch.setattr(app_module, runner, boom)
    r = client.post(path, json=VALID)
    assert r.status_code == 500
    body = r.get_json()
    assert body["status"] == "error" and body["message"] == "RuntimeError: boum"
    assert "Traceback" not in json.dumps(body)
    assert any("boum" in (rec.exc_text or "") or rec.exc_info for rec in caplog.records)
