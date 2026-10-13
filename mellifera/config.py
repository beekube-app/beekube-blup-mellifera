"""Réglages lus dans l'environnement (variables MELLIFERA_*), validés dès le démarrage."""
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Mapping
from zoneinfo import ZoneInfo, ZoneInfoNotFoundError

BASE_URL = "https://index-mellifera.org"
MIN_DELAY = 1.0          # politesse : jamais moins d'une seconde entre deux requêtes


@dataclass(frozen=True)
class Settings:
    data_dir: Path
    contact: str
    delay: float = 2.0
    max_errors: int = 5
    schedule: str = "sun 02:00"
    timezone: str = "Europe/Paris"
    base_url: str = BASE_URL

    def __post_init__(self) -> None:
        if not self.delay >= MIN_DELAY:
            raise ValueError(f"MELLIFERA_DELAY doit valoir au moins {MIN_DELAY:g} s (reçu : {self.delay:g})")
        if self.max_errors < 1:
            raise ValueError(f"MELLIFERA_MAX_ERRORS doit valoir au moins 1 (reçu : {self.max_errors})")
        try:
            ZoneInfo(self.timezone)
        except (ZoneInfoNotFoundError, ValueError) as exc:    # ZoneInfoNotFoundError hérite de KeyError
            raise ValueError(f"MELLIFERA_TZ inconnu : « {self.timezone} »") from exc

    @property
    def db_path(self) -> Path:
        return self.data_dir / "mellifera.sqlite"

    @property
    def cache_dir(self) -> Path:
        return self.data_dir / "cache"

    @property
    def lock_path(self) -> Path:
        return self.data_dir / "crawl.lock"

    @classmethod
    def from_env(cls, env: Mapping[str, str] = os.environ) -> "Settings":
        """Lève ValueError en nommant la variable fautive."""
        def number(name: str, default: str, cast: Callable[[str], float]):
            raw = env.get(name, default)
            try:
                return cast(raw)
            except ValueError:
                raise ValueError(f"{name} invalide : « {raw} » (nombre attendu)") from None

        return cls(
            data_dir=Path(env.get("MELLIFERA_DATA", "data/mellifera")),
            contact=env.get("MELLIFERA_CONTACT", "").strip(),
            delay=number("MELLIFERA_DELAY", "2.0", float),
            max_errors=number("MELLIFERA_MAX_ERRORS", "5", int),
            schedule=env.get("MELLIFERA_SCHEDULE", "sun 02:00"),
            timezone=env.get("MELLIFERA_TZ", "Europe/Paris"),
        )
