import datetime as dt
from zoneinfo import ZoneInfo

import pytest

from mellifera.schedule import next_run, parse_schedule, run_schedule

PARIS = ZoneInfo("Europe/Paris")


def test_parse_schedule():
    assert parse_schedule("sun 02:00") == (6, 2, 0)
    for bad in ("dimanche 2h", "sun 25:00", "xyz 02:00", ""):
        with pytest.raises(ValueError):
            parse_schedule(bad)


@pytest.mark.parametrize("now, expected", [
    (dt.datetime(2026, 10, 13, 9, 0, tzinfo=PARIS), dt.datetime(2026, 10, 18, 2, 0, tzinfo=PARIS)),   # mardi
    (dt.datetime(2026, 10, 18, 1, 59, tzinfo=PARIS), dt.datetime(2026, 10, 18, 2, 0, tzinfo=PARIS)),  # dimanche avant
    (dt.datetime(2026, 10, 18, 2, 0, tzinfo=PARIS), dt.datetime(2026, 10, 25, 2, 0, tzinfo=PARIS)),   # échéance passée
])
def test_next_run(now, expected):
    assert next_run(now, "sun 02:00") == expected


def test_run_schedule_waits_then_runs():
    t = {"now": dt.datetime(2026, 10, 17, 23, 0, tzinfo=PARIS)}
    sleeps, runs = [], []

    def sleep(s):
        sleeps.append(s)
        t["now"] += dt.timedelta(seconds=s)
    run_schedule("sun 02:00", "Europe/Paris", lambda: runs.append(t["now"]) or 0, now=lambda: t["now"], sleep=sleep,
                 iterations=1)
    assert runs == [dt.datetime(2026, 10, 18, 2, 0, tzinfo=PARIS)]
    assert max(sleeps) <= 3600


def test_run_schedule_survives_a_failing_job():
    t = {"now": dt.datetime(2026, 10, 17, 23, 0, tzinfo=PARIS)}
    runs = []

    def sleep(s):
        t["now"] += dt.timedelta(seconds=s)

    def job():
        runs.append(t["now"])
        t["now"] += dt.timedelta(minutes=30)
        if len(runs) == 1:
            raise RuntimeError("database is locked")
        return 0
    run_schedule("sun 02:00", "Europe/Paris", job, now=lambda: t["now"], sleep=sleep, iterations=2)
    assert runs == [dt.datetime(2026, 10, 18, 2, 0, tzinfo=PARIS), dt.datetime(2026, 10, 25, 2, 0, tzinfo=PARIS)]
