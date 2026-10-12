import pytest

import app as app_module
from v1_payloads import minimal_input


@pytest.fixture
def client():
    app_module.app.config["TESTING"] = True
    return app_module.app.test_client()


def test_welcome(client):
    resp = client.get("/")
    assert resp.status_code == 200
    assert resp.get_json() == {"status": "success", "results": "Beekube BLUP Mellifera"}


def test_missing_variety_is_rejected_with_422(client):
    payload = minimal_input()
    del payload["data"][0]["variety"]
    resp = client.post("/blup", json=payload)
    assert resp.status_code == 422
    assert resp.get_json()[0]["loc"] == ["data", 0, "variety"]


def test_minimal_input_is_valid():
    app_module.BLUPInput.model_validate(minimal_input())
