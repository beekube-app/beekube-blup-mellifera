import pytest


@pytest.fixture(autouse=True)
def _no_mellifera_db(monkeypatch):
    """Un MELLIFERA_DB exporté dans le shell du développeur ne doit pas changer les calculs des tests."""
    monkeypatch.delenv("MELLIFERA_DB", raising=False)
