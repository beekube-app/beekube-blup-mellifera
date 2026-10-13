import httpx
import pytest

from mellifera.client import ForbiddenPath, MissingContact, PoliteClient, TooManyErrors
from tests.mellifera.fake_site import FakeClock, FakeSite

BASE = "https://index-mellifera.org"


def make_client(site: FakeSite, clock: FakeClock, **kw) -> PoliteClient:
    return PoliteClient(BASE, "contact@example.org", kw.pop("delay", 2.0), kw.pop("max_errors", 5),
                        transport=site.transport(), sleep=clock.sleep, clock=clock.clock, **kw)


def test_refuses_to_start_without_contact():
    with pytest.raises(MissingContact):
        PoliteClient(BASE, "  ")


def test_user_agent_language_cookie_and_get_only():
    site = FakeSite({"/all/": b"<html></html>"})
    resp = make_client(site, FakeClock()).get("/all/")
    assert (resp.status, resp.http_status, resp.body) == ("ok", 200, b"<html></html>")
    req = site.requests[0]
    assert req.method == "GET"
    assert req.headers["user-agent"] == "BeekubeCrawler/1.0 (+contact@example.org)"
    assert "django_language=en" in req.headers["cookie"]


@pytest.mark.parametrize("path", ["/rating/2759", "/accounts/login/", "/exportYearQueens/145/2022/", "/api/queens/",
                                  "/queen/../accounts/", "/queen/../", "/queenTree/FR-CBP-B-0276-2022/3/", "/all",
                                  "/welcome/"])
def test_paths_outside_whitelist_are_refused_before_sending(path):
    site = FakeSite()
    with pytest.raises(ForbiddenPath):
        make_client(site, FakeClock()).get(path)
    assert site.requests == []


def test_delay_between_requests():
    site = FakeSite({"/all/": b"x", "/countrybyyear/8/": b"y"})
    clock = FakeClock()
    client = make_client(site, clock)
    client.get("/all/")
    clock.now += 0.5                        # 0,5 s de traitement entre les deux requêtes
    client.get("/countrybyyear/8/")
    assert clock.sleeps == [1.5]


def test_backoff_doubles_on_429_and_5xx_and_respects_retry_after():
    site = FakeSite({"/all/": b"x"})
    site.failures["/all/"] = [429, 503, 500]
    clock = FakeClock()
    client = make_client(site, clock)
    assert client.get("/all/").status == "ok"
    assert clock.sleeps == [4.0, 8.0, 16.0]     # délai doublé à chaque erreur
    clock.now += 100
    site.pages["/countrybyyear/8/"] = b"y"
    client.get("/countrybyyear/8/")
    assert clock.sleeps[-1:] == [16.0] and len(clock.sleeps) == 3   # succès : délai revenu à 2 s


def test_retry_after_longer_than_backoff_is_respected():
    site = FakeSite({"/all/": b"x"})
    site.failures["/all/"] = [(429, "120")]
    clock = FakeClock()
    make_client(site, clock).get("/all/")
    assert clock.sleeps == [120.0]


def test_retry_after_as_http_date_falls_back_to_doubling():
    site = FakeSite({"/all/": b"x"})
    site.failures["/all/"] = [(503, "Wed, 21 Oct 2026 07:28:00 GMT")]
    clock = FakeClock()
    make_client(site, clock).get("/all/")
    assert clock.sleeps == [4.0]


def test_backoff_capped_at_five_minutes():
    site = FakeSite({"/all/": b"x"})
    site.failures["/all/"] = [503] * 9
    clock = FakeClock()
    make_client(site, clock, max_errors=20).get("/all/")
    assert max(clock.sleeps) == 300.0


def test_stops_after_max_consecutive_errors():
    site = FakeSite({"/all/": b"x"})
    site.failures["/all/"] = [503] * 10
    with pytest.raises(TooManyErrors):
        make_client(site, FakeClock(), max_errors=5).get("/all/")
    assert len(site.requests) == 5


def test_timeouts_count_as_errors():
    clock = FakeClock()

    def boom(request):
        raise httpx.ReadTimeout("lent", request=request)
    client = PoliteClient(BASE, "c", 2.0, 3, transport=httpx.MockTransport(boom), sleep=clock.sleep, clock=clock.clock)
    with pytest.raises(TooManyErrors):
        client.get("/all/")


def test_404_is_missing_not_an_error():
    site = FakeSite()
    client = make_client(site, FakeClock(), max_errors=2)
    for _ in range(3):
        assert client.get("/queen/FR-XX-B-0001-2020/").status == "missing"


def test_redirect_to_login_is_forbidden_and_not_followed():
    site = FakeSite()
    site.redirects["/breederbyyear/145/2022/"] = f"{BASE}/accounts/login/?next=/breederbyyear/145/2022/"
    resp = make_client(site, FakeClock()).get("/breederbyyear/145/2022/")
    assert (resp.status, resp.http_status) == ("forbidden", 302)
    assert site.paths == ["/breederbyyear/145/2022/"]


def test_redirect_inside_whitelist_is_followed():
    site = FakeSite({"/queen/FR-CBP-B-0276-2022/": b"fiche"})
    site.redirects["/queen/fr-cbp-b-0276-2022/"] = "/queen/FR-CBP-B-0276-2022/"
    resp = make_client(site, FakeClock()).get("/queen/fr-cbp-b-0276-2022/")
    assert (resp.path, resp.status, resp.body) == ("/queen/fr-cbp-b-0276-2022/", "ok", b"fiche")


def test_redirect_to_other_host_is_forbidden():
    site = FakeSite()
    site.redirects["/all/"] = "https://pedigreeapis.org/all/"
    assert make_client(site, FakeClock()).get("/all/").status == "forbidden"
    assert len(site.requests) == 1


def test_redirect_chain_is_cut_after_three_hops():
    site = FakeSite()
    for a, b in (("A", "B"), ("B", "C"), ("C", "D"), ("D", "E")):
        site.redirects[f"/queen/{a}/"] = f"/queen/{b}/"
    resp = make_client(site, FakeClock()).get("/queen/A/")
    assert (resp.status, resp.http_status) == ("forbidden", 302)
    assert site.paths == ["/queen/A/", "/queen/B/", "/queen/C/", "/queen/D/"]
