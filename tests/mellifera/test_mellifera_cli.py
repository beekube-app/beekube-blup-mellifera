import json
import subprocess
import sys
from pathlib import Path

import pytest

from blup2.engine import run_blup
from blup2.schemas import BlupV2Input
from mellifera.cli import (EXIT_CONFIG, EXIT_LOCKED, EXIT_NETWORK, EXIT_OK, EXIT_PARSE, crawl_lock, main, run_crawl,
                           status_lines)
from mellifera.config import Settings
from tests.mellifera.fake_site import DEMO, FakeClock, demo_site

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture
def settings(tmp_path):
    return Settings(data_dir=tmp_path / "mellifera", contact="c@example.org")


def crawl(settings, site, mode="full"):
    clock = FakeClock()
    return run_crawl(settings, mode, transport=site.transport(), sleep=clock.sleep, clock=clock.clock)


def test_crawl_status_export_end_to_end(settings, tmp_path, capsys):
    assert crawl(settings, demo_site()) == EXIT_OK
    lines = status_lines(settings)
    assert lines[0].startswith("dernier passage : n° 1 full") and "terminé" in lines[0]
    assert "queen : 6" in lines and "evaluation : 4" in lines
    out = tmp_path / "cbp.json"
    assert main(["export-blup", "--breeder", "CBP", "-o", str(out)], settings=settings) == EXIT_OK
    inp = BlupV2Input.model_validate(json.loads(out.read_text()))
    assert len(inp.queens) == 5 and len(inp.records) == 4
    run_blup(inp)
    assert "reines : 5 (dont 5 sélectionnées)" in capsys.readouterr().out


def test_parse_failure_threshold_exit_code(settings):
    site = demo_site()
    site.pages[f"/queen/{DEMO['C'].code}/"] = b"<html><body>oops</body></html>"
    assert crawl(settings, site) == EXIT_PARSE


def test_network_stop_exit_code(settings):
    site = demo_site()
    site.failures["/all/"] = [503] * 10
    assert crawl(settings, site) == EXIT_NETWORK


def test_lock_held_by_a_crawl(settings, tmp_path):
    assert crawl(settings, demo_site()) == EXIT_OK
    with crawl_lock(settings.lock_path):                    # un crawl en cours (des heures pour le premier)
        assert crawl(settings, demo_site(), "incremental") == EXIT_LOCKED
        assert main(["reparse"], settings=settings) == EXIT_LOCKED
        assert main(["status"], settings=settings) == EXIT_OK          # lecture seule : jamais bloquée
        out = tmp_path / "cbp.json"
        assert main(["export-blup", "--breeder", "CBP", "-o", str(out)], settings=settings) == EXIT_OK
    assert crawl(settings, demo_site(), "incremental") == EXIT_OK


def test_missing_contact_refuses_to_start(tmp_path):
    s = Settings(data_dir=tmp_path / "m", contact="")
    assert main(["crawl", "--full"], settings=s) == EXIT_CONFIG
    assert main(["schedule"], settings=s) == EXIT_CONFIG
    assert not s.db_path.exists()


def test_invalid_environment_is_a_configuration_error(settings, monkeypatch):
    assert crawl(settings, demo_site()) == EXIT_OK
    monkeypatch.setenv("MELLIFERA_DATA", str(settings.data_dir))
    assert main(["status"]) == EXIT_OK                      # environnement valide : la base répond
    monkeypatch.setenv("MELLIFERA_TZ", "Mars/Olympus")
    assert main(["status"]) == EXIT_CONFIG
    monkeypatch.setenv("MELLIFERA_TZ", "Europe/Paris")
    monkeypatch.setenv("MELLIFERA_DELAY", "2s")
    assert main(["status"]) == EXIT_CONFIG


def test_status_and_export_without_database(settings, tmp_path):
    assert main(["status"], settings=settings) == EXIT_CONFIG
    assert main(["export-blup", "--breeder", "CBP", "-o", str(tmp_path / "x.json")], settings=settings) == EXIT_CONFIG


def test_reparse_command(settings):
    assert crawl(settings, demo_site()) == EXIT_OK
    assert main(["reparse"], settings=settings) == EXIT_OK
    assert main(["reparse", "--kind", "season"], settings=settings) == EXIT_OK


def test_cli_does_not_import_numpy():
    """L'image du service mellifera n'installe ni numpy ni scipy : la CLI se charge sans eux (l'export est
    vérifié par test_mellifera_export)."""
    code = ("import sys, mellifera.cli; "
            "heavy = sorted(m for m in sys.modules if m.split('.')[0] in ('numpy', 'scipy')); "
            "assert not heavy, heavy")
    subprocess.run([sys.executable, "-c", code], cwd=ROOT, check=True)


@pytest.mark.parametrize("status", [403, 404])
def test_unreadable_country_list_is_not_a_finished_run(settings, status):
    """/all/ refusée ou absente : rien n'a été parcouru, le passage n'est pas « terminé » et sort en erreur."""
    site = demo_site()
    site.failures["/all/"] = [status]
    assert crawl(settings, site) == EXIT_NETWORK
    assert "inachevé" in status_lines(settings)[0]
    site.requests.clear()
    assert crawl(settings, site) == EXIT_OK                 # la reprise redemande /all/
    assert site.paths[0] == "/all/"


def test_status_lists_unfinished_runs(settings):
    site = demo_site()
    site.failures[f"/queen/{DEMO['B'].code}/"] = [503] * 10
    assert crawl(settings, site) == EXIT_NETWORK
    assert "passages inachevés : n° 1 full" in status_lines(settings)
