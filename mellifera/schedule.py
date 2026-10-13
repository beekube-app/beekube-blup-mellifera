"""Planification hebdomadaire du crawl incrémental (« sun 02:00 », fuseau Europe/Paris par défaut)."""
from __future__ import annotations

import datetime as dt
import logging
import re
import time
from typing import Callable, Optional
from zoneinfo import ZoneInfo

log = logging.getLogger("mellifera")
DAYS = ("mon", "tue", "wed", "thu", "fri", "sat", "sun")


def parse_schedule(spec: str) -> tuple[int, int, int]:
    """« sun 02:00 » → (jour de la semaine 0–6, heure, minute)."""
    m = re.fullmatch(r"\s*([a-z]{3})\s+(\d{1,2}):(\d{2})\s*", spec.lower())
    if not m or m.group(1) not in DAYS or int(m.group(2)) > 23 or int(m.group(3)) > 59:
        raise ValueError(f"MELLIFERA_SCHEDULE invalide : « {spec} » (attendu : « sun 02:00 »)")
    return DAYS.index(m.group(1)), int(m.group(2)), int(m.group(3))


def next_run(now: dt.datetime, spec: str) -> dt.datetime:
    """Prochaine échéance strictement après `now` (datetime avec fuseau), dans le fuseau de `now`."""
    weekday, hour, minute = parse_schedule(spec)
    days = (weekday - now.weekday()) % 7
    candidate = (now + dt.timedelta(days=days)).replace(hour=hour, minute=minute, second=0, microsecond=0)
    if candidate <= now:
        candidate += dt.timedelta(days=7)
    return candidate


def run_schedule(spec: str, timezone: str, job: Callable[[], int], *,
                 now: Optional[Callable[[], dt.datetime]] = None, sleep: Callable[[float], None] = time.sleep,
                 iterations: Optional[int] = None) -> None:
    """Boucle : attend l'échéance, lance `job`, recommence. `iterations` borne la boucle (tests)."""
    tz = ZoneInfo(timezone)
    now = now or (lambda: dt.datetime.now(tz))
    parse_schedule(spec)
    done = 0
    while iterations is None or done < iterations:
        target = next_run(now(), spec)
        log.info("prochain crawl incrémental : %s", target.isoformat())
        while (wait := (target - now()).total_seconds()) > 0:
            sleep(min(wait, 3600))
        try:
            code = job()
        except Exception:               # bogue ou base indisponible : on réessaie à l'échéance suivante
            log.exception("échec du crawl incrémental ; nouvel essai à l'échéance suivante")
        else:
            log.info("crawl incrémental terminé (code %s)", code)
        done += 1
