"""Catalogue des caractères standard (spec volet 3 §2)."""
import csv
import json
import os
import typing

from blup2.schemas import CriterionIn
from blup2.traits import STANDARD_TRAITS, TraitName

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def site_traits():
    with open(os.path.join(ROOT, "mellifera", "traits.csv"), encoding="utf-8") as f:
        return {r["name"]: r for r in csv.DictReader(f)}


def test_every_site_trait_is_mapped_exactly_once_with_its_kind_and_direction():
    site = site_traits()
    cited = [t.mellifera for t in STANDARD_TRAITS.values() if t.mellifera is not None]
    assert sorted(cited) == sorted(site)                       # chaque nom existe et est cité une seule fois
    for t in STANDARD_TRAITS.values():
        if t.mellifera is not None:
            assert (t.kind, t.direction) == (site[t.mellifera]["kind"], int(site[t.mellifera]["direction"]))


def test_beebreed_and_kehrle_mappings_are_one_to_one():
    for source in ("beebreed", "kehrle"):
        names = [getattr(t, source) for t in STANDARD_TRAITS.values() if getattr(t, source) is not None]
        assert len(names) == len(set(names)), source


def test_catalogue_shape():
    assert len(STANDARD_TRAITS) == 30
    assert all(k == k.lower() and " " not in k for k in STANDARD_TRAITS)
    assert {t.kind for t in STANDARD_TRAITS.values()} == {"note", "percentage", "harvest"}
    assert [k for k, t in STANDARD_TRAITS.items() if t.direction == -1] == ["varroa_infestation"]


def test_contract_enumeration_is_the_catalogue():
    assert list(typing.get_args(TraitName)) == list(STANDARD_TRAITS)
    schema = CriterionIn.model_json_schema()["properties"]["trait"]
    enum = next(s["enum"] for s in schema["anyOf"] if "enum" in s)
    assert enum == list(STANDARD_TRAITS)
    assert "Index Mellifera mildness" in schema["description"] and "BeeBreed Defensive behavior" in schema["description"]
    with open(os.path.join(ROOT, "openapi", "openapi.json"), encoding="utf-8") as f:
        committed = json.load(f)["components"]["schemas"]["CriterionIn"]["properties"]["trait"]
    assert next(s["enum"] for s in committed["anyOf"] if "enum" in s) == list(STANDARD_TRAITS)
