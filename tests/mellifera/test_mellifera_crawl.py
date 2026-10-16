import datetime as dt
import logging

import pytest

import mellifera.ingest
from mellifera.cache import RawCache
from mellifera.client import PoliteClient
from mellifera.crawl import Crawler
from mellifera.db import open_db
from mellifera.urls import kind_of, queen_path
from tests.mellifera.fake_site import (DEMO, E, FakeClock, FakeSite, Q, country_page, db_snapshot, demo_site,
                                       queen_page, queens_list_page, season_page)

A, B, C = DEMO["A"], DEMO["B"], DEMO["C"]
QUEENS = tuple(DEMO.values())

TODAY = dt.date(2024, 6, 1)


def make_crawler(tmp_path, site: FakeSite, max_errors: int = 5):
    conn = open_db(tmp_path / "m.sqlite")
    clock = FakeClock()
    client = PoliteClient("https://index-mellifera.org", "c@example.org", 2.0, max_errors,
                          transport=site.transport(), sleep=clock.sleep, clock=clock.clock)
    return Crawler(conn, client, RawCache(tmp_path / "cache"), today=lambda: TODAY), conn


def parsed_codes(conn) -> set[str]:
    return {r[0] for r in conn.execute("SELECT code FROM queen WHERE parsed_at IS NOT NULL")}


def test_full_crawl_with_closure(tmp_path):
    site = demo_site()
    crawler, conn = make_crawler(tmp_path, site)
    report = crawler.full()
    assert not report.stopped and report.failing_kinds == {}
    assert parsed_codes(conn) == {q.code for q in QUEENS}
    assert conn.execute("SELECT status FROM page WHERE url = '/queen/FR-OLD-B-0001-2019/'").fetchone()[0] == "missing"
    assert all(kind_of(p) for p in site.paths)
    assert len(site.paths) == len(set(site.paths))              # aucune page demandée deux fois
    assert report.stats["evaluations_unresolved"] == 0          # « M100 » résolue après la fermeture
    assert report.stats["evaluations_resolved"] == 1
    assert conn.execute("SELECT COUNT(*) FROM evaluation WHERE queen_code = ?", (A.code,)).fetchone()[0] == 1
    assert conn.execute("SELECT finished_at IS NOT NULL FROM crawl_run").fetchone()[0] == 1


def test_resume_after_network_stop(tmp_path):
    site = demo_site()
    site.failures[queen_path(B.code)] = [503] * 5
    crawler, conn = make_crawler(tmp_path, site, max_errors=5)
    first = crawler.full()
    assert first.stopped
    assert conn.execute("SELECT finished_at FROM crawl_run").fetchone()[0] is None
    before = set(site.paths)
    conn.close()

    crawler, conn = make_crawler(tmp_path, site)
    site.requests.clear()
    second = crawler.full()
    assert second.run_id == first.run_id and not second.stopped
    already = before - {queen_path(B.code)}
    assert not already & set(site.paths)                        # rien de ce qui était lu n'est redemandé
    assert parsed_codes(conn) == {q.code for q in QUEENS}

    fresh, fresh_conn = make_crawler(tmp_path / "fresh", demo_site())
    fresh.full()
    assert db_snapshot(conn)["evaluation"] == db_snapshot(fresh_conn)["evaluation"]


def test_incremental_reads_only_new_or_changed(tmp_path):
    site = demo_site()
    crawler, conn = make_crawler(tmp_path, site)
    crawler.full()
    conn.close()

    b2 = Q(B.code, 2022, 145, mother=B.mother, cotations=2)           # nouvelle cotation
    n = Q("FR-CBP-B-0004-2022", 2022, 145)                          # nouvelle reine
    o = Q("FR-ABC-B-0001-2015", 2015, 300)                          # reine d'un nouvel éleveur
    site.pages.update({
        "/countrybyyear/8/": country_page([(145, "CBP", [2022]), (200, "NG", [2024]), (300, "ABC", [2015])]),
        "/breederbyyear/145/2022/": queens_list_page(145, "CBP", 8, [2022, 2021], [A, b2, n]),
        queen_path(b2.code): queen_page(b2), queen_path(n.code): queen_page(n), queen_path(o.code): queen_page(o),
        "/breederbyyear/300/2015/": queens_list_page(300, "ABC", 8, [2015], [o]),
        "/breeder_Year/300/2015/": season_page([E(o.code, 2015, {"mildness": "5"}, tester="ABC", tester_id=300)]),
    })
    site.requests.clear()
    crawler, conn = make_crawler(tmp_path, site)
    report = crawler.incremental()
    assert not report.stopped
    queens = {p for p in site.paths if kind_of(p) == "queen"}
    assert queens == {queen_path(c) for c in (B.code, n.code, o.code, "FR-OLD-B-0001-2019")}
    assert "/breederbyyear/145/2021/" not in site.paths             # année < Y − 2
    assert {"/breeder_Year/145/2022/", "/breeder_Year/200/2024/", "/breederbyyear/300/2015/",
            "/breeder_Year/300/2015/"} <= set(site.paths)           # nouvel éleveur : parcouru en entier
    assert conn.execute("SELECT cotation_count FROM queen WHERE code = ?", (B.code,)).fetchone()[0] == 2
    assert n.code in parsed_codes(conn) and o.code in parsed_codes(conn)


def test_reparse_rebuilds_the_same_database_without_network(tmp_path):
    site = demo_site()
    crawler, conn = make_crawler(tmp_path, site)
    crawler.full()
    before = db_snapshot(conn)
    report = Crawler(conn, None, RawCache(tmp_path / "cache"), today=lambda: TODAY).reparse()
    assert not report.stopped and report.failing_kinds == {}
    assert db_snapshot(conn) == before


def test_parse_failures_are_recorded_and_reported(tmp_path):
    site = demo_site()
    site.pages[queen_path(C.code)] = b"<html><body><h1>Server error page</h1></body></html>"
    crawler, conn = make_crawler(tmp_path, site)
    report = crawler.full()
    assert report.failing_kinds == {"queen": pytest.approx(1 / 6, abs=1e-3)}
    err = conn.execute("SELECT parse_error FROM page WHERE url = ?", (queen_path(C.code),)).fetchone()[0]
    assert "fiche absente" in err
    assert len(parsed_codes(conn)) == 5                             # le crawl a continué


def test_unexpected_parser_exception_does_not_block_the_crawl(tmp_path, monkeypatch):
    real = mellifera.ingest.store_queen

    def buggy(conn, q, path, seen_at):
        real(conn, q, path, seen_at)
        if q.code == C.code:
            raise AttributeError("'NoneType' object has no attribute 'get'")   # bogue face à une structure inattendue
    monkeypatch.setattr(mellifera.ingest, "store_queen", buggy)
    crawler, conn = make_crawler(tmp_path, demo_site())
    report = crawler.full()
    assert not report.stopped                                       # passage terminé : la reprise ne bute pas dessus
    assert report.failing_kinds == {"queen": pytest.approx(1 / 6, abs=1e-3)}
    err = conn.execute("SELECT parse_error FROM page WHERE url = ?", (queen_path(C.code),)).fetchone()[0]
    assert err.startswith("AttributeError")
    assert C.code not in parsed_codes(conn)                         # écritures de la page annulées


def test_nonconforming_codes_are_kept_and_counted(tmp_path):
    odd = Q("B4536(NG)", 2024, 200)                                 # code mis en URL, hors format international
    site = demo_site()
    site.pages["/breederbyyear/200/2024/"] = queens_list_page(200, "NG", 8, [2024], [DEMO["D"], odd])
    site.pages[queen_path(odd.code)] = queen_page(odd)
    crawler, conn = make_crawler(tmp_path, site)
    report = crawler.full()
    assert report.stats["nonconforming_codes"] == 1
    assert conn.execute("SELECT code_valid FROM queen WHERE code = ?", (odd.code,)).fetchone()[0] == 0


def test_incremental_resumes_an_interrupted_full_crawl(tmp_path):
    """Un crawl complet interrompu n'est pas masqué par l'incrémental : celui-ci reprend le complet."""
    site = demo_site()
    site.failures[queen_path(B.code)] = [503] * 5
    crawler, conn = make_crawler(tmp_path, site)
    first = crawler.full()
    assert first.stopped
    conn.close()
    crawler, conn = make_crawler(tmp_path, site)
    report = crawler.incremental()
    assert (report.mode, report.run_id, report.stopped) == ("full", first.run_id, False)
    assert parsed_codes(conn) == {q.code for q in QUEENS}
    assert conn.execute("SELECT COUNT(*) FROM crawl_run").fetchone()[0] == 1


def test_incremental_on_an_empty_database_runs_a_full_crawl(tmp_path):
    crawler, conn = make_crawler(tmp_path, demo_site())
    report = crawler.incremental()
    assert report.mode == "full" and not report.stopped
    assert parsed_codes(conn) == {q.code for q in QUEENS}


def test_queen_page_that_never_answers_does_not_block_the_full_crawl(tmp_path, caplog):
    """Une fiche qui ne répond jamais est notée en erreur et le passage complet se termine (crawl du 15/10 bloqué
    sur /queen/DE-LB-C-010b-2008/)."""
    site = demo_site()
    site.failures[queen_path(B.code)] = ["timeout"] * 50
    crawler, conn = make_crawler(tmp_path, site, max_errors=10)
    with caplog.at_level(logging.WARNING, logger="mellifera"):
        report = crawler.full()
    assert not report.stopped
    assert report.stats["pages"]["queen"]["error"] == 1
    assert conn.execute("SELECT status FROM page WHERE url = ?", (queen_path(B.code),)).fetchone()[0] == "error"
    assert parsed_codes(conn) == {q.code for q in QUEENS} - {B.code}
    assert site.paths.count(queen_path(B.code)) == 3
    assert f"{queen_path(B.code)} illisible" in caplog.text
    assert conn.execute("SELECT finished_at IS NOT NULL FROM crawl_run").fetchone()[0] == 1
