"""Exporte la spécification OpenAPI de l'API dans openapi/openapi.json.

Usage : .venv/bin/python scripts/export_openapi.py
Le fichier est versionné ; tests/test_openapi.py vérifie qu'il est à jour.
"""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

OPENAPI_PATH = os.path.join(ROOT, "openapi", "openapi.json")
# Hôte public (celui codé par défaut dans le SDK PHP) : le fichier ne dépend pas de API_BASE_URL
PUBLIC_SERVER_URL = "https://mellifera.beekube.com"


def build_spec():
    from app import app

    spec = json.loads(json.dumps(app.api_doc))
    spec["servers"] = [{"url": PUBLIC_SERVER_URL}]
    return spec


def render_spec():
    return json.dumps(build_spec(), indent=2, ensure_ascii=False, sort_keys=True) + "\n"


def main():
    os.makedirs(os.path.dirname(OPENAPI_PATH), exist_ok=True)
    with open(OPENAPI_PATH, "w", encoding="utf-8") as f:
        f.write(render_spec())
    print(f"OpenAPI written to {OPENAPI_PATH}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
