"""Cache brut des pages : une version par URL, compressée (gzip), adressée par sha1(url)."""
from __future__ import annotations

import gzip
import hashlib
import os
from pathlib import Path


class RawCache:
    def __init__(self, root: Path):
        self.root = Path(root)

    def relative_path(self, url: str) -> str:
        h = hashlib.sha1(url.encode("utf-8")).hexdigest()
        return f"{h[:2]}/{h}.html.gz"

    def put(self, url: str, body: bytes) -> tuple[str, str]:
        """Écrit la page (remplace la version précédente) ; retourne (chemin relatif, sha256 du contenu)."""
        rel = self.relative_path(url)
        path = self.root / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(".tmp")
        tmp.write_bytes(gzip.compress(body, mtime=0))
        os.replace(tmp, path)
        return rel, hashlib.sha256(body).hexdigest()

    def get(self, rel: str) -> bytes:
        return gzip.decompress((self.root / rel).read_bytes())
