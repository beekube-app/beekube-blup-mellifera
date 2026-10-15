import hashlib

import pytest

from mellifera.cache import RawCache
from mellifera.config import Settings
from mellifera.urls import kind_of, path_ids, queen_path, queens_list_path, season_path


# ------------------------------------------------------------------ réglages
def test_settings_from_env():
    s = Settings.from_env({"MELLIFERA_DATA": "/x", "MELLIFERA_CONTACT": " me@example.org ", "MELLIFERA_DELAY": "3"})
    assert (str(s.db_path), s.contact, s.delay, s.max_errors, s.schedule) == (
        "/x/mellifera.sqlite", "me@example.org", 3.0, 10, "sun 02:00")      # 10 erreurs ≈ 20 min de panne
    assert (str(s.cache_dir), str(s.lock_path)) == ("/x/cache", "/x/crawl.lock")
    assert Settings.from_env({}).contact == ""


@pytest.mark.parametrize("env, variable", [
    ({"MELLIFERA_DELAY": "2s"}, "MELLIFERA_DELAY"),
    ({"MELLIFERA_DELAY": "0.2"}, "MELLIFERA_DELAY"),           # politesse : au moins 1 s entre deux requêtes
    ({"MELLIFERA_MAX_ERRORS": "0"}, "MELLIFERA_MAX_ERRORS"),
    ({"MELLIFERA_MAX_ERRORS": "cinq"}, "MELLIFERA_MAX_ERRORS"),
    ({"MELLIFERA_TZ": "Europe/Pariss"}, "MELLIFERA_TZ"),
])
def test_invalid_settings_name_the_variable(env, variable):
    with pytest.raises(ValueError, match=variable):
        Settings.from_env(env)


# ------------------------------------------------------------------ chemins autorisés
def test_kind_of_and_path_ids():
    assert kind_of("/all/") == "countries"
    assert kind_of("/countrybyyear/8/") == "country"
    assert kind_of(queens_list_path(145, 2022)) == "queens_list"
    assert (kind_of(season_path(145, 2022)), path_ids(season_path(145, 2022))) == ("season", ("145", "2022"))
    assert path_ids("/queen/B562(PJ)1dr/") == ("B562(PJ)1dr",)
    for path in ("/rating/2759", "/accounts/login/", "/all", "/queen/../", "/queen/./"):
        assert kind_of(path) is None
    with pytest.raises(ValueError):
        path_ids("/accounts/login/")


def test_queen_path_refuses_codes_unsafe_in_an_url():
    assert queen_path("FR-CBP-B-0276-2022") == "/queen/FR-CBP-B-0276-2022/"
    assert queen_path("B562(PJ)1dr") == "/queen/B562(PJ)1dr/"
    for bad in ("", None, "FR CBP", "FR/CBP", "..", ".", "été"):
        assert queen_path(bad) is None


# ------------------------------------------------------------------ cache brut
def test_raw_cache_keeps_one_compressed_version_per_url(tmp_path):
    cache = RawCache(tmp_path)
    rel, sha = cache.put("/all/", b"<html>1</html>")
    assert rel == cache.relative_path("/all/") and rel.endswith(".html.gz") and rel[:2] == rel[3:5]
    assert sha == hashlib.sha256(b"<html>1</html>").hexdigest()
    assert cache.get(rel) == b"<html>1</html>"
    assert cache.put("/all/", b"<html>2</html>")[0] == rel
    assert cache.get(rel) == b"<html>2</html>"
    assert [p.name for p in tmp_path.rglob("*") if p.is_file()] == [rel.split("/")[1]]   # ni .tmp ni ancienne version
