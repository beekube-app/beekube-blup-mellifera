import json
import os
import re
import subprocess
import threading
import time

import pytest

import app as app_module
from v1_payloads import minimal_input, valid_output

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


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


def fake_rscript(output=None, returncode=0, stdout="", stderr="", calls=None, delay=0.0):
    """Remplace subprocess.run : lit l'entrée cmd[2], écrit `output` dans cmd[3]."""
    def _run(cmd, **kwargs):
        with open(cmd[2]) as f:
            received = json.load(f)
        if calls is not None:
            calls.append({"cmd": list(cmd), "kwargs": kwargs, "input": received})
        if delay:
            time.sleep(delay)
        if output is not None:
            with open(cmd[3], "w") as f:
                json.dump(output(received) if callable(output) else output, f)
        return subprocess.CompletedProcess(cmd, returncode, stdout, stderr)
    return _run


def assert_error_response(resp, message_part):
    assert resp.status_code == 500
    body = resp.get_json()
    assert set(body) == {"status", "message", "output", "error", "results"}
    assert body["status"] == "error"
    assert message_part in body["message"]
    return body


def test_success_returns_validated_results(client, monkeypatch):
    monkeypatch.setattr(app_module.subprocess, "run", fake_rscript(valid_output()))
    resp = client.post("/blup", json=minimal_input())
    assert resp.status_code == 200
    body = resp.get_json()
    assert body["status"] == "success"
    assert body["results"] == valid_output()
    assert body["results"]["heritabilities"]["10"]["se"] is None


def test_rscript_gets_private_paths_and_timeout(client, monkeypatch):
    calls = []
    monkeypatch.setattr(app_module.subprocess, "run", fake_rscript(valid_output(), calls=calls))
    client.post("/blup", json=minimal_input())
    cmd = calls[0]["cmd"]
    assert cmd[0] == "Rscript"
    assert cmd[1] == app_module.BLUP_SCRIPT
    workdir = os.path.dirname(cmd[2])
    assert os.path.dirname(cmd[3]) == workdir
    assert not os.path.abspath(workdir).startswith(os.path.join(app_module.APP_DIR, "data"))
    assert calls[0]["kwargs"]["timeout"] == app_module.RSCRIPT_TIMEOUT
    assert calls[0]["input"] == minimal_input()
    assert not os.path.exists(workdir), "le répertoire temporaire doit être supprimé"


def test_each_request_uses_its_own_directory(client, monkeypatch):
    calls = []
    monkeypatch.setattr(app_module.subprocess, "run", fake_rscript(valid_output(), calls=calls))
    client.post("/blup", json=minimal_input())
    client.post("/blup", json=minimal_input())
    assert os.path.dirname(calls[0]["cmd"][2]) != os.path.dirname(calls[1]["cmd"][2])


def test_concurrent_requests_get_their_own_results(monkeypatch):
    def echo_exploitation(received):
        out = valid_output()
        out["blup"][0]["queenbee"] = str(received["exploitation"])
        return out

    monkeypatch.setattr(app_module.subprocess, "run", fake_rscript(echo_exploitation, delay=0.05))
    results = {}

    def worker(exploitation):
        payload = minimal_input()
        payload["exploitation"] = exploitation
        resp = app_module.app.test_client().post("/blup", json=payload)
        results[exploitation] = resp.get_json()["results"]["blup"][0]["queenbee"]

    threads = [threading.Thread(target=worker, args=(100 + i,)) for i in range(8)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert results == {100 + i: str(100 + i) for i in range(8)}


def test_nonzero_exit_returns_500_with_stderr(client, monkeypatch):
    monkeypatch.setattr(app_module.subprocess, "run",
                        fake_rscript(None, returncode=1, stdout="out", stderr="Error in lmer"))
    body = assert_error_response(client.post("/blup", json=minimal_input()), "code 1")
    assert body["error"] == "Error in lmer"
    assert body["output"] == "out"


def test_timeout_returns_500(client, monkeypatch):
    def too_slow(cmd, **kwargs):
        raise subprocess.TimeoutExpired(cmd, kwargs["timeout"], output=b"partial", stderr=b"slow")
    monkeypatch.setattr(app_module.subprocess, "run", too_slow)
    body = assert_error_response(client.post("/blup", json=minimal_input()), "timed out")
    assert body["error"] == "slow"


def test_missing_output_returns_500(client, monkeypatch):
    monkeypatch.setattr(app_module.subprocess, "run", fake_rscript(None))
    assert_error_response(client.post("/blup", json=minimal_input()), "Results file not found")


def test_non_conforming_output_returns_500(client, monkeypatch):
    bad = valid_output()
    bad["blup"][0]["blups"]["10"] = "abc"
    monkeypatch.setattr(app_module.subprocess, "run", fake_rscript(bad))
    assert_error_response(client.post("/blup", json=minimal_input()), "Invalid output")


def test_empty_data_returns_empty_results_without_r(client, monkeypatch):
    def must_not_run(cmd, **kwargs):
        raise AssertionError("Rscript ne doit pas être lancé")
    monkeypatch.setattr(app_module.subprocess, "run", must_not_run)
    payload = minimal_input()
    payload["data"] = []
    resp = client.post("/blup", json=payload)
    assert resp.status_code == 200
    assert resp.get_json() == {"status": "success", "results": {"blup": [], "heritabilities": {}}}


def test_gunicorn_timeout_exceeds_rscript_timeout():
    with open(os.path.join(ROOT, "Dockerfile-Python")) as f:
        dockerfile = f.read()
    cmd = json.loads(re.search(r"^CMD (\[.*\])$", dockerfile, re.M).group(1))
    assert int(cmd[cmd.index("--timeout") + 1]) == 180
    assert app_module.RSCRIPT_TIMEOUT < 180
    assert int(cmd[cmd.index("--workers") + 1]) >= 2
