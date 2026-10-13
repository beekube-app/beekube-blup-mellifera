"""Client HTTP poli : GET seulement, liste blanche, délai, ralentissement sur erreur, pas de connexion."""
from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Callable, Literal, Optional
from urllib.parse import urljoin, urlsplit

import httpx

from .urls import kind_of

MAX_BACKOFF = 300.0          # plafond du délai doublé (5 min)
MAX_RETRY_AFTER = 3600.0     # un Retry-After plus long est ramené à 1 h
MAX_REDIRECTS = 3
LOGIN_PATH = "/accounts/login/"


class MissingContact(ValueError):
    """MELLIFERA_CONTACT vide : le crawl refuse de démarrer."""


class ForbiddenPath(ValueError):
    """Chemin hors liste blanche : refusé avant l'envoi."""


class TooManyErrors(RuntimeError):
    """Trop d'erreurs réseau consécutives : arrêt propre du crawl."""


@dataclass(frozen=True)
class Response:
    path: str
    status: Literal["ok", "missing", "forbidden", "error"]
    http_status: Optional[int]
    body: bytes = b""


def _retry_after(resp: httpx.Response) -> float:
    try:
        return min(MAX_RETRY_AFTER, max(0.0, float(resp.headers.get("retry-after", "0"))))
    except ValueError:
        return 0.0               # date HTTP : ignorée, le doublement du délai s'applique


class PoliteClient:
    def __init__(self, base_url: str, contact: str, delay: float = 2.0, max_errors: int = 5, *,
                 transport: Optional[httpx.BaseTransport] = None,
                 sleep: Callable[[float], None] = time.sleep,
                 clock: Callable[[], float] = time.monotonic,
                 timeout: float = 30.0):
        if not contact.strip():
            raise MissingContact("MELLIFERA_CONTACT est vide : renseignez un contact (e-mail ou URL) avant de crawler")
        self.base_url = base_url.rstrip("/")
        self.user_agent = f"BeekubeCrawler/1.0 (+{contact.strip()})"
        self.delay = delay
        self.max_errors = max_errors
        self._sleep, self._clock = sleep, clock
        self._current_delay = delay
        self._errors = 0
        self._last: Optional[float] = None
        self._http = httpx.Client(
            base_url=self.base_url, transport=transport, timeout=timeout, follow_redirects=False,
            headers={"User-Agent": self.user_agent, "Accept-Language": "en", "Cookie": "django_language=en"},
        )

    def close(self) -> None:
        self._http.close()

    def __enter__(self) -> "PoliteClient":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def _wait(self) -> None:
        if self._last is not None:
            remaining = self._current_delay - (self._clock() - self._last)
            if remaining > 0:
                self._sleep(remaining)

    def _ok(self) -> None:
        self._errors = 0
        self._current_delay = self.delay

    def _failure(self, path: str, why: str, retry_after: float = 0.0) -> None:
        self._errors += 1
        if self._errors >= self.max_errors:
            raise TooManyErrors(f"{self._errors} erreurs consécutives (dernière : {why} sur {path})")
        self._current_delay = max(min(self._current_delay * 2, MAX_BACKOFF), retry_after)

    def get(self, path: str) -> Response:
        if kind_of(path) is None:
            raise ForbiddenPath(path)
        current, redirects = path, 0
        while True:
            self._wait()
            try:
                r = self._http.get(current)
            except httpx.TransportError as exc:          # délai dépassé, connexion refusée…
                self._failure(current, type(exc).__name__)
                continue
            finally:
                self._last = self._clock()               # le délai court depuis la fin de la réponse
                self._http.cookies.clear()               # jamais de cookie posé par le site (session, csrf)
            code = r.status_code
            if code == 200:
                self._ok()
                return Response(path, "ok", code, r.content)
            if code in (301, 302, 303, 307, 308):
                self._ok()
                target = urlsplit(urljoin(self.base_url + current, r.headers.get("location", "")))
                same_host = target.netloc in ("", urlsplit(self.base_url).netloc)
                if (not same_host or target.path.startswith(LOGIN_PATH) or kind_of(target.path) is None
                        or redirects >= MAX_REDIRECTS):
                    return Response(path, "forbidden", code)
                current, redirects = target.path, redirects + 1
                continue
            if code in (404, 410):
                self._ok()
                return Response(path, "missing", code)
            if code in (401, 403):
                self._ok()
                return Response(path, "forbidden", code)
            if code == 429 or code >= 500:
                self._failure(current, f"HTTP {code}", _retry_after(r))
                continue
            self._ok()
            return Response(path, "error", code)
