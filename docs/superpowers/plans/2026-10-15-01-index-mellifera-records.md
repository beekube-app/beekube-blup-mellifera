# Index Mellifera, volet 3 : notes des lignées dans le BLUP et catalogue des caractères — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Quand un critère Beekube déclare un caractère standard (`CriterionIn.trait`) et que des reines de l'exploitation citent des lignées Index Mellifera, le moteur v2 lit dans la base du crawler les notes des ancêtres résolus et de leurs contemporaines, et les ajoute au modèle du critère, pour que la fiabilité des reines apparentées augmente, sans qu'aucune reine externe ni aucune note externe n'apparaisse dans les sorties.

**Architecture:** Nouveau catalogue figé `blup2/traits.py` (30 caractères, correspondances Index Mellifera, BeeBreed, Karl Kehrle), exposé par l'énumération `CriterionIn.trait` du contrat 2.2.0. Nouveau module `blup2/external_records.py` : lecture SQL des notes (ancêtres, puis groupes testeur × année), nettoyage par motif, contemporaines, construction des `RecordIn` à id négatifs. `blup2/external.py` l'appelle après l'ascendance du volet 2 et renvoie contemporaines et notes dans `Expansion`. `prepare_criterion` place ces notes après celles de l'exploitation, avec des effets fixes propres, et `CriterionData` garde le plan de l'exploitation seule pour les a priori et les statistiques. Enfin la méthode de fiabilité se décide sans les équations dues aux seules notes externes : au-delà de la limite à cause d'elles, les PEV des reines de l'exploitation restent exactes (`mme.inverse_entries`).

**Tech Stack:** Python ≥ 3.11 (3.14 en local, 3.12 en CI et dans l'image `api`) ; `sqlite3` de la bibliothèque standard (`json_each`, URI `mode=ro`) ; numpy, scipy (factorisation creuse déjà en place) ; pydantic 2.12.3 ; pytest 8.3.3. Aucune dépendance ajoutée.

**Spec:** `docs/superpowers/specs/2026-10-15-index-mellifera-records-design.md` (71b53ed, amendée en 6bc6823 : contemporaines, aucun effet fixe partagé, fiabilité exacte). Les volets 1 et 2 sont fusionnés dans `main` (632045e).

**Prérequis :**
- worktree `/Users/tups/Sites/beekube-blup-external`, branche `mellifera-records`, sans modification en cours ;
- `.venv` du worktree créé (Python 3.14, `requirements-dev.txt`) ;
- suite actuelle : `420 passed, 12 deselected`, vérifiée à l'étape 1 de la tâche 1.

**Prototype :** tout le code de ce plan a été écrit et exécuté avant rédaction : suite complète `462 passed, 15 deselected` sous Python 3.14 et 3.12, tests lents `9 passed, 468 deselected`, ruff (F, E9) sans erreur. Le test de la tâche 4 sur la méthode de fiabilité a été vu en échec sans le code de la tâche 4, profil `breeder` compris. Le plan a ensuite été rejoué de bout en bout dans une copie du dépôt : remplacements appliqués dans l'ordre, commandes lancées, sorties et commits conformes aux lignes « Expected ». Le code se transcrit tel quel. Un comportement différent de celui annoncé est un défaut du plan, à signaler.

## Écarts / précisions de spec

Points où ce plan précise la spec. Chacun est appliqué par la tâche indiquée.

1. **Message du doublon de caractère** : « critères 1 et 2 déclarent le même caractère standard gentleness » (sans accents graves autour du nom, comme les autres messages de validation) ; avec trois critères, « critères 1, 2 et 3 … ». Tâche 1.
2. **Motifs d'écart des contemporaines.** Une note de contemporaine qui échoue au nettoyage (en kg, sans date, hors échelle, non numérique, sans rucher) ne peut être rattachée à aucun groupe : elle est ignorée **sans être comptée**, comme les notes des groupes que ne partage aucun ancêtre. Seuls « sans fiche de reine » et « antérieures à la naissance » sont comptés pour les contemporaines. Les notes des ancêtres comptent tous les motifs. Tâche 3.
3. **Ordre plafond puis groupes.** Le plafond de 1000 contemporaines (codes croissants) s'applique avant le filtre des groupes : un groupe peut ensuite tomber sous deux reines et ses notes d'ancêtres être comptées « sans contemporaine ». Tâche 3.
4. **`external_records_added` absent** quand aucune note n'est ajoutée ni écartée (aucun groupe lu, ou caractères sans note dans la base). Il est retiré, comme `external_ref_resolved`, quand le repli du volet 2 écarte les ancêtres. Tâche 3.
5. **Consanguinité des ouvrières** : `f_workers` est calculé pour toutes les reines du pedigree (`ped.queen_ids` et `ped.added_queen_ids`), puis séparé en deux covariables, `f_workers` (exploitation) et `f_workers_ext` (notes externes), chacune retenue seulement si elle varie sur ses lignes. Tâches 2 et 3.
6. **Équations externes** (règle de §5) : `X.shape[1] − prior_X.shape[1]` (effets fixes propres) + reines externes notées (environnement permanent) + `n_genetic × external_nodes`. `external_nodes` est la différence de `n_nodes` entre le pedigree complet et un second pedigree construit avec les seuls ancêtres ; il n'est calculé que s'il y a des contemporaines. Dans la branche exacte par `inverse_entries`, seules les PEV des reines de l'exploitation sont lues ; les autres nœuds restent à NaN et ne sortent jamais. Tâche 4.
7. **Catalogue** : un `TraitName` (`Literal` des 30 identifiants) sert de type au champ ; la description OpenAPI liste chaque caractère avec ses correspondances (`trait_description()`). Tâche 1.
8. **Coût de la requête des groupes** sur la vraie base : la base du crawler n'a pas d'index sur `evaluation.tester_breeder_id`. Le test lent de la tâche 5 mesure le coût sur une base synthétique de 9 000 reines (calcul complet sous la seconde) ; la vraie base est à mesurer après le crawl complet, hors de ce plan.

## Global Constraints

- Objectif (spec, décisions) : « **de meilleures estimations des reines de l'exploitation** » ; « Les reines externes restent invisibles dans les sorties (isolement du volet 2 conservé). »
- Bit à bit (§1.2) : « Une exploitation sans `trait`, ou sans référence externe, ou sans `MELLIFERA_DB`, obtient des résultats **strictement identiques** à la version 2.1.0. »
- Robustesse (§1.3, §6) : « Aucune donnée de la base du crawler ne fait échouer un calcul » ; toute erreur de lecture des notes : `logger.exception`, avertissement `external_records_unavailable`, calcul avec l'ascendance sans les notes. « Jamais de 500. »
- Contrat (§3) : `CriterionIn.trait : Optional[Literal[<identifiants du catalogue>]] = None` ; deux critères sur un même caractère refusés ; aucune option nouvelle ; `MELLIFERA_DB` absente : rien ne change et aucun avertissement.
- Statistiques (§3, §5) : `n_records`, `n_queens`, `n_contemporary_groups`, compteurs par reine, a priori de variance, `should_estimate`, `few_records` « **ne comptent que l'exploitation** ».
- Aucun effet fixe partagé (§5) : ruchers et évaluateurs négatifs propres, niveaux `x…` de classe de ruche et d'âge, covariable `f_workers_ext`.
- Fiabilité (§5) : la méthode se décide sur les équations de l'exploitation et de son ascendance ; jamais de bascule en approché à cause des seules notes externes.
- Récoltes (§2) : les notes `harvest` d'Index Mellifera (kilos) ne sont **pas** chargées. Les valeurs d'élevage BeeBreed ne sont jamais injectées.
- `MAX_CONTEMPORARIES = 1000` (§4.4).
- Version `2.2.0` : `ENGINE_VERSION`, `app.py`, `METHODOLOGIE.md`, `README.fr.md`, `README.md`, snapshot OpenAPI régénéré. **`MELLIFERA_DB` n'est activé en production qu'après régénération et publication du SDK PHP 2.2.0.**
- `blup2/` n'importe rien de `mellifera/` : le contrat avec le crawler tient aux noms des tables `evaluation`, `criterion`, `trait`, `queen` et de leurs colonnes, couvert par le test de contrat.
- Tests : aucun test n'appelle le vrai site ; `MELLIFERA_CONTACT=` forcé dans les commandes. Noms de fichiers de test uniques ; jamais de `__init__.py` dans `tests/` ni dans `tests/mellifera/`. Helpers importés par `from tests.external_base import …` et `from tests.simdata import …`.
- Commandes lancées depuis la racine du worktree : `.venv/bin/pytest …`, `.venv/bin/python …`, `uvx ruff check --select F,E9 blup2 tests`.
- Git :
  - travailler dans le worktree `/Users/tups/Sites/beekube-blup-external`, branche `mellifera-records` ;
  - **toujours** `git add <chemins listés>`, jamais `git add -A` ni `git add .` ;
  - pas de push ;
  - messages en français, terminés par `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>` puis `Claude-Session: https://claude.ai/code/session_015J8yR1HR7mLSRtwnvEYap6`.
- Ne rien modifier dans le checkout `/Users/tups/Sites/beekube-blup` (branche `main`) : `ROADMAP.fr.md` et `docs/` y portent des modifications sans rapport. Ne pas toucher au Pi (crawl en cours).
- Données personnelles : les fixtures et bases de test n'ont ni nom ni adresse d'éleveur, seulement des id et des trigrammes.
- Langue : messages, commentaires, `METHODOLOGIE.md` et README FR en français, avec les accents ; README EN en anglais.

## Review Focus

1. **Exploitation qui envoie des id négatifs** (reines, ruchers, évaluateurs ; autorisés par `BlupV2Input`). Attendu : contemporaines, ruchers et évaluateurs externes toujours sous ceux de l'entrée, aucune collision. Test : tâche 3, `test_negative_farm_ids_keep_external_ids_below`.
2. **Liens erronés entre contemporaines dans la base** (reine mère d'elle-même, cycle entre deux contemporaines). Attendu : liens coupés, calcul réussi. Test : tâche 3, `test_cycles_and_self_references_between_contemporaries_are_cut`.
3. **REML forcé (`estimate_variances = "always"`) avec des notes externes.** Attendu : le REML voit toutes les notes, mais `n_records` et `n_queens` du critère restent ceux de l'exploitation. Test : tâche 3, `test_reml_with_external_notes_counts_only_the_farm`.
4. **Base en cours d'écriture par le crawl** pendant la lecture des notes. Attendu : lecture du dernier état validé (WAL), notes lues sans « database is locked ». Test : tâche 3, `test_notes_are_read_while_the_crawler_writes`.
5. **`external_pedigree_depth = 0`** avec des critères qui déclarent un caractère. Attendu : ni ascendance ni notes, sorties identiques au calcul sans base. Test : tâche 3, `test_depth_zero_disables_the_notes_too`.

---

## Carte des fichiers

| Fichier | Rôle | Tâche |
|---|---|---|
| `blup2/traits.py` | catalogue `STANDARD_TRAITS`, `TraitName`, description OpenAPI | 1 |
| `blup2/schemas.py` | `CriterionIn.trait`, refus des doublons, deux codes d'avertissement | 1 |
| `blup2/__init__.py`, `app.py`, `METHODOLOGIE.md` (version) | version 2.2.0 | 1 |
| `openapi/openapi.json` | snapshot régénéré | 1 |
| `tests/test_traits.py`, `tests/test_schemas.py`, `tests/test_api.py`, `tests/test_engine.py`, `tests/fixtures/expected-example-v2.json` | catalogue et contrat 2.2.0 | 1 |
| `blup2/records.py` | notes externes dans `prepare_criterion`, effets fixes propres, plan de l'exploitation seule | 2 |
| `tests/test_records.py` | niveaux propres, notes ignorées sans note de l'exploitation, covariable externe | 2 |
| `blup2/external_records.py` | lecture, nettoyage, contemporaines, `RecordIn` externes, message | 3 |
| `blup2/external.py` | appel après l'ascendance, `Expansion.contemporaries` et `.records`, repli | 3 |
| `blup2/engine.py` | câblage des reines et notes externes, repli | 3, 4 |
| `blup2/model_standard.py`, `blup2/model_breeder.py` | a priori et compteurs sur l'exploitation seule ; branche de fiabilité exacte | 3, 4 |
| `tests/external_base.py` | tables des notes dans la mini-base | 3 |
| `tests/test_external_records.py`, `tests/test_external_engine.py`, `tests/mellifera/test_mellifera_external_contract.py` | règles de §4, isolement, contrat avec la vraie base | 3, 4 |
| `blup2/mme.py`, `tests/test_mme.py` | `inverse_entries` | 4 |
| `tests/simdata.py`, `tests/test_external_calibration.py` | simulation (§7) et coût (§9), tests lents | 5 |
| `METHODOLOGIE.md` §3 bis, §7, §10, §11, `README.fr.md`, `README.md` | documentation | 5 |


---

### Task 1 : Catalogue des caractères et contrat 2.2.0

Le catalogue figé des caractères standard (spec §2) et le champ `trait` des critères (§3), avec le refus de deux critères sur un même caractère, les deux codes d'avertissement du volet et la version 2.2.0. Aucun calcul ne change : le moteur ne lit pas encore `trait`.

**Files:**
- Create: `blup2/traits.py`
- Modify: `blup2/schemas.py`
- Modify: `blup2/__init__.py`
- Modify: `app.py`
- Modify: `METHODOLOGIE.md`
- Modify: `openapi/openapi.json` (régénéré)
- Test (create): `tests/test_traits.py`
- Test (modify): `tests/test_api.py`
- Test (modify): `tests/test_engine.py`
- Test (modify): `tests/fixtures/expected-example-v2.json`
- Test (modify): `tests/test_schemas.py`

**Interfaces:**
- Consumes : `mellifera/traits.csv` (colonne `name`), lu par le seul test.
- Produces :
  - `blup2.traits.Trait` : dataclass figée `(label: str, kind: Literal["note", "percentage", "harvest"], mellifera: Optional[str] = None, beebreed: Optional[str] = None, kehrle: Optional[str] = None, direction: Literal[1, -1] = 1)` ;
  - `blup2.traits.STANDARD_TRAITS: dict[str, Trait]` (ordre du tableau de la spec) et `TraitName` (`Literal` de ses clés) ;
  - `blup2.traits.trait_description() -> str` ;
  - `CriterionIn.trait: Optional[TraitName] = None` ;
  - codes `external_records_added` et `external_records_unavailable` dans `WarningOut.code` ;
  - `blup2.ENGINE_VERSION == "2.2.0"`.

- [ ] **Step 1 : Vérifier le point de départ**

```bash
git branch --show-current
git status --short
MELLIFERA_CONTACT= .venv/bin/pytest -q -p no:cacheprovider | tail -1
```

Expected :

```text
mellifera-records
420 passed, 12 deselected in …
```

- [ ] **Step 2 : Écrire les tests du catalogue et du contrat**

Créer `tests/test_traits.py` :

```python
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
```

`tests/test_api.py` : remplacer :

```python
    assert r.status_code == 200
    body = r.get_json()
    assert body["status"] == "success" and body["engine_version"] == "2.1.0"
    assert [q["id"] for q in body["queens"]] == [4074]
    assert set(body) == {"status", "engine_version", "profile", "parameters", "base", "queens", "matings",
```

par :

```python
    assert r.status_code == 200
    body = r.get_json()
    assert body["status"] == "success" and body["engine_version"] == "2.2.0"
    assert [q["id"] for q in body["queens"]] == [4074]
    assert set(body) == {"status", "engine_version", "profile", "parameters", "base", "queens", "matings",
```

`tests/test_api.py` : remplacer :

```python
        assert responses["422"]["content"]["application/json"]["schema"] == {"$ref": "#/components/schemas/ErrorV2Out"}
        assert "500" in responses
    assert doc["info"]["version"] == "2.1.0"
    assert paths["/blup"]["post"]["operationId"] == "calulate_blup_post"     # SDK v1 inchangé
    # Contrat avec le SDK PHP (plan 04, tâche A1) : une seule classe BlupApi, méthodes nommées
```

par :

```python
        assert responses["422"]["content"]["application/json"]["schema"] == {"$ref": "#/components/schemas/ErrorV2Out"}
        assert "500" in responses
    assert doc["info"]["version"] == "2.2.0"
    assert paths["/blup"]["post"]["operationId"] == "calulate_blup_post"     # SDK v1 inchangé
    # Contrat avec le SDK PHP (plan 04, tâche A1) : une seule classe BlupApi, méthodes nommées
```

`tests/test_engine.py` : remplacer :

```python
def test_run_blup_output_shape():
    out = run_blup(_inp())
    assert out.status == "success" and out.engine_version == ENGINE_VERSION == "2.1.0"
    assert set(out.parameters) == {"1", "2"} and [q.id for q in out.queens] == [1, 6, 7, 2, 3, 4, 5]
    assert set(out.queens[0].indexes) == {"1", "2"}
```

par :

```python
def test_run_blup_output_shape():
    out = run_blup(_inp())
    assert out.status == "success" and out.engine_version == ENGINE_VERSION == "2.2.0"
    assert set(out.parameters) == {"1", "2"} and [q.id for q in out.queens] == [1, 6, 7, 2, 3, 4, 5]
    assert set(out.queens[0].indexes) == {"1", "2"}
```

`tests/fixtures/expected-example-v2.json` : remplacer :

```json
{
 "engine_version": "2.1.0",
 "n_components": 3,
 "parameters": {
```

par :

```json
{
 "engine_version": "2.2.0",
 "n_components": 3,
 "parameters": {
```

`tests/test_schemas.py` : remplacer :

```python
        assert WarningOut(code=code, message="m").code == code

```

par :

```python
        assert WarningOut(code=code, message="m").code == code


def test_criterion_trait_comes_from_the_catalogue():
    ok = BlupV2Input.model_validate(mutate(("criteria", 0, "trait"), "gentleness"))
    assert ok.criteria[0].trait == "gentleness"
    assert BlupV2Input.model_validate(VALID).criteria[0].trait is None
    assert_invalid(mutate(("criteria", 0, "trait"), "douceur"), "trait")


def test_two_criteria_on_the_same_trait_are_refused():
    d = copy.deepcopy(VALID)
    d["criteria"] = [dict(d["criteria"][0], id=k, trait="gentleness") for k in (44, 45, 46)]
    assert_invalid(d, "critères 44, 45 et 46 déclarent le même caractère standard gentleness")


def test_external_records_warning_codes_are_part_of_the_contract():
    from blup2.schemas import WarningOut

    for code in ("external_records_added", "external_records_unavailable"):
        assert WarningOut(code=code, message="m").code == code

```

- [ ] **Step 3 : Vérifier l'échec**

Run : `MELLIFERA_CONTACT= .venv/bin/pytest -q -p no:cacheprovider -rfE --tb=no --continue-on-collection-errors tests/test_traits.py tests/test_schemas.py tests/test_api.py tests/test_engine.py | tail -16`

`tests/test_traits.py` ne s'importe pas (`blup2.traits` n'existe pas) ; les tests du champ `trait` et des codes échouent ; la version attendue est 2.2.0 au lieu de 2.1.0.

Expected :

```text
................................FFFF....F...F.......                     [100%]
=========================== short test summary info ============================
FAILED tests/test_schemas.py::test_criterion_trait_comes_from_the_catalogue
FAILED tests/test_schemas.py::test_two_criteria_on_the_same_trait_are_refused
FAILED tests/test_schemas.py::test_external_records_warning_codes_are_part_of_the_contract
FAILED tests/test_api.py::test_v2_blup_returns_engine_output - AssertionError...
FAILED tests/test_api.py::test_openapi_documents_v1_and_v2_routes - Assertion...
FAILED tests/test_engine.py::test_run_blup_output_shape - AssertionError: ass...
ERROR tests/test_traits.py
6 failed, 46 passed, 1 error in …
```

- [ ] **Step 4 : Écrire le catalogue, le champ `trait` et la version 2.2.0**

Créer `blup2/traits.py` :

```python
"""Catalogue des caractères standard (spec volet 3 §2), figé dans le contrat et donc dans le SDK.

Un critère Beekube déclare son caractère (`CriterionIn.trait`) ; le moteur s'en sert pour rattacher les notes
Index Mellifera des ancêtres résolus (colonne `mellifera`, noms de `mellifera/traits.csv`). Les colonnes
BeeBreed et Karl Kehrle sont documentaires, pour les exports futurs de Beekube. Correspondances une pour une.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Optional


@dataclass(frozen=True)
class Trait:
    label: str
    kind: Literal["note", "percentage", "harvest"]
    mellifera: Optional[str] = None      # trait.name de la base du crawler (mellifera/traits.csv)
    beebreed: Optional[str] = None       # caractère publié par BeeBreed
    kehrle: Optional[str] = None         # champ du registre Karl Kehrle
    direction: Literal[1, -1] = 1        # 1 : valeur haute = mieux


STANDARD_TRAITS: dict[str, Trait] = {
    "honey_yield": Trait("Récolte de miel (kg)", "harvest", "honeyAmt", "Honey yield", "honeykgpromean"),
    "honey_yield_average": Trait("Récolte de miel moyenne (kg)", "harvest", "aveHonAmt"),
    "honey_score": Trait("Indice miel (note)", "note", "honeyIdx", None, "honeyindex"),
    "average_crop_score": Trait("Indice de récolte moyen (note)", "note", "aveCropIdx"),
    "spring_crop": Trait("Récolte de printemps (note)", "note", "sprCrop", None, "cropearly"),
    "mid_season_crop": Trait("Récolte de mi-saison (note)", "note", "midCrop"),
    "summer_crop": Trait("Récolte d'été (note)", "note", "sumCrop", None, "cropsumm"),
    "autumn_crop": Trait("Récolte d'automne (note)", "note", "falCrop", None, "croplate"),
    "pollen_crop": Trait("Récolte de pollen (note)", "note", "polCrop"),
    "gentleness": Trait("Douceur", "note", "mildness", "Defensive behavior", "temperofbees"),
    "calmness": Trait("Tenue au cadre", "note", "framHold", "Steadiness on comb", "calmbehaviour"),
    "swarming": Trait("Essaimage (note haute = essaime peu)", "note", "swarming", "Swarming drive", "swarmtendency"),
    "overwintering": Trait("Hivernage", "note", "winter", "Robustness in winter", "overwintering"),
    "spring_development": Trait("Développement de printemps", "note", None, "Development in Spring"),
    "colony_strength": Trait("Force de la colonie", "note", None, "Colony strength"),
    "brood_vitality": Trait("Vitalité du couvain", "note", "broVital", None, "vitalityofbrood"),
    "bee_vitality": Trait("Vitalité des abeilles", "note", "beeVital", None, "vitalityofbees"),
    "fertility": Trait("Fécondité", "note", "fertil", None, "fertility"),
    "comb_building": Trait("Bâtisse", "note", "bldCombs", None, "buildingcombs"),
    "propolis": Trait("Propolis", "note", "propolis", None, "buildingpropol"),
    "hygienic_24h": Trait("Test d'aiguille, couvain nettoyé à 24 h (%)", "percentage", "hyg24h",
                          "rate of opened cells", "hygprozafter24h"),
    "hygienic_48h": Trait("Test d'aiguille, couvain nettoyé à 48 h (%)", "percentage", "hyg48h", None,
                          "hygienicprozafter48h"),
    "hygienic_score": Trait("Indice de comportement hygiénique (note)", "note", "hygBehaIdx", None, "hygienicindex"),
    "varroa_score": Trait("Indice varroa (note)", "note", "varIdx", "Varroa index"),
    "varroa_infestation": Trait("Infestation varroa (mesure)", "percentage", None,
                                "Varroa infestation development", "varroanumber", direction=-1),
    "smr": Trait("SMR / VSH (%)", "percentage", None, "SMR", "varroavsh"),
    "recapping": Trait("Recapping (%)", "percentage", None, "Recapping"),
    "chalkbrood": Trait("Couvain plâtré (note haute = résistant)", "note", None, "Chalkbrood"),
    "nosema": Trait("Nosémose (note haute = résistant)", "note", None, "Nosemosis", "nosemaindex"),
    "cbpv": Trait("Paralysie chronique (note haute = résistant)", "note", None, "Chronic Bee Paralysis"),
}

TraitName = Literal[tuple(STANDARD_TRAITS)]


def trait_description() -> str:
    """Description OpenAPI de CriterionIn.trait : le catalogue et ses correspondances."""
    lines = ["Caractère standard du critère ; sert au moteur à rattacher les notes Index Mellifera des ancêtres "
             "résolus, et à Beekube pour ses exports. Catalogue (identifiant : libellé [type] ; Index Mellifera, "
             "BeeBreed, Karl Kehrle) :"]
    for key, t in STANDARD_TRAITS.items():
        refs = " ; ".join(f"{src} {name}" for src, name in (("Index Mellifera", t.mellifera), ("BeeBreed", t.beebreed),
                                                           ("Karl Kehrle", t.kehrle)) if name)
        sens = ", valeur haute = moins bien" if t.direction == -1 else ""
        lines.append(f"- {key} : {t.label} [{t.kind}{sens}]" + (f" ; {refs}" if refs else ""))
    return "\n".join(lines)
```

`blup2/schemas.py` : remplacer :

```python

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

BORN_MIN = 1900   # born : entre 1900 et l'année courante + 1 (spec §5.1)
```

par :

```python

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from .traits import TraitName, trait_description

BORN_MIN = 1900   # born : entre 1900 et l'année courante + 1 (spec §5.1)
```

`blup2/schemas.py` : remplacer :

```python
    queen_ratio: Optional[float] = Field(None, gt=0)
    rg: Optional[float] = Field(None, ge=-1, le=1)

    @model_validator(mode="after")
```

par :

```python
    queen_ratio: Optional[float] = Field(None, gt=0)
    rg: Optional[float] = Field(None, ge=-1, le=1)
    trait: Optional[TraitName] = Field(None, description=trait_description())

    @model_validator(mode="after")
```

`blup2/schemas.py` : remplacer :

```python
        if dc:
            errors.append(f"identifiants de critères dupliqués : {dc}")
        qset, cmap = set(qids), {c.id: c for c in self.criteria}
        unknown_q = sorted({r.queen for r in self.records if r.queen not in qset})
```

par :

```python
        if dc:
            errors.append(f"identifiants de critères dupliqués : {dc}")
        by_trait: dict[str, list[int]] = {}
        for c in self.criteria:
            if c.trait is not None:
                by_trait.setdefault(c.trait, []).append(c.id)
        for trait, ids in by_trait.items():
            if len(ids) > 1:
                errors.append(f"critères {', '.join(map(str, ids[:-1]))} et {ids[-1]} déclarent le même "
                              f"caractère standard {trait}")
        qset, cmap = set(qids), {c.id: c for c in self.criteria}
        unknown_q = sorted({r.queen for r in self.records if r.queen not in qset})
```

`blup2/schemas.py` : remplacer :

```python
        "external_ref_unknown",
        "external_pedigree_unavailable",
    ]
    message: str
```

par :

```python
        "external_ref_unknown",
        "external_pedigree_unavailable",
        "external_records_added",
        "external_records_unavailable",
    ]
    message: str
```

`blup2/__init__.py` : remplacer :

```python
"""Moteur BLUP v2 (abeille) : parenté Brascamp & Bijma, modèle animal, index base 100."""

ENGINE_VERSION = "2.1.0"

```

par :

```python
"""Moteur BLUP v2 (abeille) : parenté Brascamp & Bijma, modèle animal, index base 100."""

ENGINE_VERSION = "2.2.0"

```

`app.py` : remplacer :

```python
RSCRIPT_TIMEOUT = int(os.getenv('RSCRIPT_TIMEOUT', '170'))

info = Info(title="Beekube BLUP Mellifera API", version="2.1.0", description="SDK PHP pour l'API BLUP de Beekube - Calcul génétique des abeilles mellifères")
servers = [
    Server(url=API_BASE_URL),
```

par :

```python
RSCRIPT_TIMEOUT = int(os.getenv('RSCRIPT_TIMEOUT', '170'))

info = Info(title="Beekube BLUP Mellifera API", version="2.2.0", description="SDK PHP pour l'API BLUP de Beekube - Calcul génétique des abeilles mellifères")
servers = [
    Server(url=API_BASE_URL),
```

`METHODOLOGIE.md` : remplacer :

```markdown
lisent les index dans Beekube, et aux personnes qui relisent le code. Les références entre crochets
renvoient à la section « Références » ; [BB] = Brascamp & Bijma 2014, [BER08] = Bienefeld,
Ehrhardt & Reinhardt 2008. Version du moteur : `2.1.0`, profils `standard` et `breeder`.

## 1. Pourquoi un BLUP avec parenté
```

par :

```markdown
lisent les index dans Beekube, et aux personnes qui relisent le code. Les références entre crochets
renvoient à la section « Références » ; [BB] = Brascamp & Bijma 2014, [BER08] = Bienefeld,
Ehrhardt & Reinhardt 2008. Version du moteur : `2.2.0`, profils `standard` et `breeder`.

## 1. Pourquoi un BLUP avec parenté
```

- [ ] **Step 5 : Vérifier que le snapshot OpenAPI est périmé**

Run : `MELLIFERA_CONTACT= .venv/bin/pytest -q -p no:cacheprovider tests/test_openapi.py tests/test_traits.py -rfE --tb=no | tail -3`

Le snapshot commité n'a pas encore le champ `trait` : le test du snapshot et celui de l'énumération échouent.

Expected :

```text
FAILED tests/test_openapi.py::test_committed_openapi_is_up_to_date - Assertio...
FAILED tests/test_traits.py::test_contract_enumeration_is_the_catalogue - Key...
2 failed, 9 passed in …
```

- [ ] **Step 6 : Régénérer le snapshot**

```bash
.venv/bin/python scripts/export_openapi.py
git diff --stat openapi/openapi.json
```

Expected : une ligne « OpenAPI written to … », puis :

```text
 openapi/openapi.json | 50 ++++++++++++++++++++++++++++++++++++++++++++++++--
 1 file changed, 48 insertions(+), 2 deletions(-)
```

- [ ] **Step 7 : Vérifier que les tests passent**

Run : `MELLIFERA_CONTACT= .venv/bin/pytest -q -p no:cacheprovider tests/test_traits.py tests/test_schemas.py tests/test_api.py tests/test_engine.py | tail -1`

Expected :

```text
56 passed in …
```

- [ ] **Step 8 : Lancer la suite complète**

```bash
MELLIFERA_CONTACT= .venv/bin/pytest -q -p no:cacheprovider | tail -1
```

Expected :

```text
427 passed, 12 deselected in …
```

Puis ruff :

```bash
uvx ruff check --select F,E9 blup2 tests
```

Expected :

```text
All checks passed!
```

- [ ] **Step 9 : Commit**

```bash
git add blup2/traits.py blup2/schemas.py blup2/__init__.py app.py METHODOLOGIE.md openapi/openapi.json tests/test_traits.py tests/test_api.py tests/test_engine.py tests/fixtures/expected-example-v2.json tests/test_schemas.py
git status --short
git commit -q -F - <<'EOF'
Ajoute le catalogue des caractères standard et le champ trait des critères (contrat 2.2.0)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015J8yR1HR7mLSRtwnvEYap6
EOF
git log --oneline -1
```

Expected : `git status --short` ne liste que les 11 fichiers ci-dessus, puis :

```text
<hash> Ajoute le catalogue des caractères standard et le champ trait des critères (contrat 2.2.0)
```

---

### Task 2 : Notes externes dans la préparation d'un critère, effets fixes séparés

`prepare_criterion` accepte les reines et notes externes (spec §5) : leurs observations suivent celles de l'exploitation, avec leurs propres niveaux de classe de ruche et d'âge (préfixe `x`) et leur propre covariable de consanguinité ; elles ne sont ajoutées que si l'exploitation a des notes sur le critère. `CriterionData` garde le plan des effets fixes de l'exploitation seule pour les a priori. Le moteur ne passe encore aucune note : aucun calcul ne change.

**Files:**
- Modify: `blup2/records.py`
- Test (modify): `tests/test_records.py`

**Interfaces:**
- Consumes : rien des tâches précédentes (`QueenIn`, `RecordIn` existent).
- Produces :
  - `prepare_criterion(inp, criterion_id, f_workers: Optional[dict[int, float]] = None, extra_queens: Iterable[QueenIn] = (), extra_records: Iterable[RecordIn] = ()) -> CriterionData` ;
  - `CriterionData.n_farm_obs: Optional[int]` (observations de l'exploitation, les premières ; `None` : toutes) et `CriterionData.farm_X: Optional[sp.csr_matrix]` ;
  - propriétés `CriterionData.farm_queens -> list[int]` (reines de l'exploitation notées), `prior_X -> sp.csr_matrix`, `prior_y -> np.ndarray` (exploitation seule, ou `X` et `y` sans notes externes) ;
  - effets fixes nommés `hive_class=x…`, `age=x…`, covariables `f_workers` et `f_workers_ext`.

- [ ] **Step 1 : Écrire les tests de la préparation**

`tests/test_records.py` : remplacer :

```python
    assert d.n_evaluators[1] == 1

```

par :

```python
    assert d.n_evaluators[1] == 1


# ------------------------------------------------------------------ notes externes (spec volet 3 §5)
def _with_external():
    from blup2.schemas import RecordIn

    qs = [queen(k, 2022) for k in range(1, 4)]
    farm = [record(1, 0.2), record(2, 0.6, hive_class="nuc"), record(3, 0.9, evaluator=None)]
    ext = [RecordIn(queen=-1, criterion=1, value=0.5, date="2023-07-01", apiary=-2, evaluator=-2),
           RecordIn(queen=-2, criterion=1, value=0.8, date="2023-07-01", apiary=-2, evaluator=-2),
           RecordIn(queen=-1, criterion=2, value=0.1, date="2023-07-01", apiary=-2, evaluator=-2)]
    inp = BlupV2Input(exploitation=1, criteria=[criterion(1), criterion(2)], queens=qs, records=farm)
    return inp, [queen(-1, 2022), queen(-2, 2020)], ext


def test_external_records_follow_the_farm_with_their_own_levels():
    inp, anc, ext = _with_external()
    fw = {1: 0.0, 2: 0.1, 3: 0.2, -1: 0.0, -2: 0.3}
    d = prepare_criterion(inp, 1, f_workers=fw, extra_queens=anc, extra_records=ext)
    ref = prepare_criterion(inp, 1, f_workers=fw)
    assert d.queen_ids == [1, 2, 3, -1, -2] and d.n_farm_obs == 3
    assert d.factor_labels["age"] == ["2+", "2+", "2+", "x1", "x2+"]
    assert d.factor_labels["hive_class"] == ["hive", "nuc", "hive", "xhive", "xhive"]
    assert d.factor_labels["evaluator"] == [1, 1, -1, -2, -2]
    assert d.farm_queens == [1, 2, 3] and d.recorded_queens == [1, 2, 3, -1, -2]
    assert (d.n_records, d.n_evaluators, d.n_years) == (ref.n_records, ref.n_evaluators, ref.n_years)
    assert (d.n_raw_records, d.n_contemporary_groups) == (ref.n_raw_records, ref.n_contemporary_groups)
    assert (d.prior_X != ref.X).nnz == 0 and d.prior_y.tolist() == ref.y.tolist()


def test_external_records_are_ignored_when_the_farm_has_no_record_on_the_criterion():
    inp, anc, ext = _with_external()
    d = prepare_criterion(inp, 2, extra_queens=anc, extra_records=ext)
    assert d.n_obs == 0 and d.recorded_queens == []


def test_worker_inbreeding_of_external_records_has_its_own_covariate():
    from blup2.schemas import RecordIn

    qs = [queen(k, 2022) for k in range(1, 4)]
    inp = BlupV2Input(exploitation=1, criteria=[criterion(1)], queens=qs,
                      records=[record(k, 0.2 * k) for k in range(1, 4)])
    anc = [queen(-k, 2021) for k in range(1, 4)]
    ext = [RecordIn(queen=-k, criterion=1, value=0.1 * k, date="2023-07-01", apiary=-2, evaluator=-2)
           for k in range(1, 4)]
    fw = {1: 0.0, 2: 0.1, 3: 0.25, -1: 0.0, -2: 0.2, -3: 0.05}
    d = prepare_criterion(inp, 1, f_workers=fw, extra_queens=anc, extra_records=ext)
    X = d.X.toarray()
    assert X.shape[1] == 4 and d.x_names[-2:] == ["f_workers", "f_workers_ext"]   # + constante, groupe externe
    assert X[:, 2].tolist() == [0.0, 0.1, 0.25, 0.0, 0.0, 0.0]
    assert X[:, 3].tolist() == [0.0, 0.0, 0.0, 0.0, 0.2, 0.05]

```

- [ ] **Step 2 : Vérifier l'échec**

Run : `MELLIFERA_CONTACT= .venv/bin/pytest -q -p no:cacheprovider -rfE --tb=no --continue-on-collection-errors tests/test_records.py | tail -12`

Les trois nouveaux tests échouent : `prepare_criterion` n'accepte pas `extra_queens` ni `extra_records`.

Expected :

```text
...............FFF                                                       [100%]
=========================== short test summary info ============================
FAILED tests/test_records.py::test_external_records_follow_the_farm_with_their_own_levels
FAILED tests/test_records.py::test_external_records_are_ignored_when_the_farm_has_no_record_on_the_criterion
FAILED tests/test_records.py::test_worker_inbreeding_of_external_records_has_its_own_covariate
3 failed, 15 passed in …
```

- [ ] **Step 3 : Étendre `prepare_criterion` et `CriterionData`**

`blup2/records.py` : remplacer :

```python

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
```

par :

```python

from dataclasses import dataclass, field
from typing import Iterable, Optional

import numpy as np
```

`blup2/records.py` : remplacer :

```python
import scipy.sparse as sp

from .schemas import BlupV2Input, CriterionIn, RecordIn


```

par :

```python
import scipy.sparse as sp

from .schemas import BlupV2Input, CriterionIn, QueenIn, RecordIn


```

`blup2/records.py` : remplacer :

```python
    n_contemporary_groups: int = 0
    scale_mixed: bool = False

    @property
    def n_obs(self) -> int:
        return int(self.y.shape[0])


```

par :

```python
    n_contemporary_groups: int = 0
    scale_mixed: bool = False
    n_farm_obs: Optional[int] = None      # observations de l'exploitation (les premières) ; None : toutes
    farm_X: Optional[sp.csr_matrix] = None    # effets fixes de ces seules observations ; None : X

    @property
    def n_obs(self) -> int:
        return int(self.y.shape[0])

    @property
    def farm_queens(self) -> list[int]:
        """Reines de l'exploitation notées sur le critère (sans les ancêtres Index Mellifera)."""
        return [q for q in self.recorded_queens if q in self.n_records]

    @property
    def prior_X(self) -> sp.csr_matrix:
        """Effets fixes servant aux a priori des variances : exploitation seule (spec volet 3 §5)."""
        return self.X if self.farm_X is None else self.farm_X

    @property
    def prior_y(self) -> np.ndarray:
        return self.y if self.n_farm_obs is None else self.y[: self.n_farm_obs]


```

`blup2/records.py` : remplacer :

```python


def prepare_criterion(
    inp: BlupV2Input, criterion_id: int, f_workers: Optional[dict[int, float]] = None
) -> CriterionData:
    crit = next(c for c in inp.criteria if c.id == criterion_id)
    born = {q.id: q.born for q in inp.queens}
    recs = [r for r in inp.records if r.criterion == criterion_id]
    scales = {(r.scale_min, r.scale_max) for r in recs}
    groups: dict[tuple, list[RecordIn]] = {}
    for r in recs:
        # une visite = un rucher et une année : des notes de même visit mais de rucher ou d'année
        # différents restent des observations distinctes (groupe de contemporaines sans ambiguïté)
        if r.visit is not None:
            key = (r.queen, "v", r.visit, r.apiary, r.date.year)
        else:
            key = (r.queen, "d", r.date, r.apiary)
        groups.setdefault(key, []).append(r)
    y, qids, cg, ev, hc, ag = [], [], [], [], [], []
    for key, rs in groups.items():
        first = rs[0]
        y.append(float(np.mean([normalized_value(r, crit.kind) for r in rs])))
        qids.append(key[0])
        cg.append((first.apiary, first.date.year))
        ev.append(first.evaluator if first.evaluator is not None else -1)
        hc.append(first.hive_class)
        ag.append(age_class(first.date.year, born.get(key[0])))
    n = len(y)
    cols, names, factors = [np.ones(n)], ["intercept"], []
    labels: dict[str, list] = {}
    for fname, values in (("cg", cg), ("evaluator", ev), ("hive_class", hc), ("age", ag)):
        levels = sorted(set(values), key=repr)
        if len(levels) < 2:
```

par :

```python


def _design(n: int, factors_in: list[tuple[str, list]], covariates: list[tuple[str, np.ndarray]]):
    """Matrice des effets fixes : constante, facteurs à au moins deux niveaux, covariables ; colonnes
    linéairement indépendantes. Retourne (X dense, noms, facteurs retenus, niveaux par facteur retenu)."""
    cols, names, factors = [np.ones(n)], ["intercept"], []
    labels: dict[str, list] = {}
    for fname, values in factors_in:
        levels = sorted(set(values), key=repr)
        if len(levels) < 2:
```

`blup2/records.py` : remplacer :

```python
            cols.append((codes == k).astype(float))
            names.append(f"{fname}={lv}")
    if f_workers is not None and n > 0:
        fw = np.array([f_workers.get(q, 0.0) for q in qids])
        if np.var(fw) > 1e-12:
            factors.append("f_workers")
            cols.append(fw)
            names.append("f_workers")
    Xd = np.column_stack(cols) if n > 0 else np.zeros((0, 1))
    keep = _independent_columns(Xd) if n > 0 else [0]
```

par :

```python
            cols.append((codes == k).astype(float))
            names.append(f"{fname}={lv}")
    for cname, values in covariates:
        factors.append(cname)
        cols.append(values)
        names.append(cname)
    Xd = np.column_stack(cols) if n > 0 else np.zeros((0, 1))
    keep = _independent_columns(Xd) if n > 0 else [0]
```

`blup2/records.py` : remplacer :

```python
    factors = [f for f in factors if any(nm == f or nm.startswith(f + "=") for nm in names)]
    labels = {f: v for f, v in labels.items() if f in factors}
    recorded = list(dict.fromkeys(qids))
    n_rec: dict[int, int] = {}
```

par :

```python
    factors = [f for f in factors if any(nm == f or nm.startswith(f + "=") for nm in names)]
    labels = {f: v for f, v in labels.items() if f in factors}
    return Xd, names, factors, labels


def _covariate(name: str, values: np.ndarray, rows: np.ndarray) -> list[tuple[str, np.ndarray]]:
    """Covariable retenue si elle varie sur ses lignes (`rows`) ; nulle ailleurs."""
    if rows.any() and np.var(values[rows]) > 1e-12:
        return [(name, np.where(rows, values, 0.0))]
    return []


def prepare_criterion(
    inp: BlupV2Input, criterion_id: int, f_workers: Optional[dict[int, float]] = None,
    extra_queens: Iterable[QueenIn] = (), extra_records: Iterable[RecordIn] = (),
) -> CriterionData:
    """extra_queens / extra_records : ancêtres Index Mellifera et leurs notes (spec volet 3 §5), ajoutés au
    modèle seulement si l'exploitation a des notes sur ce critère. Leurs observations suivent celles de
    l'exploitation ; elles ont leurs propres groupes (ruchers négatifs), leurs propres classes d'âge et leur
    propre covariable de consanguinité. Statistiques et a priori : exploitation seule."""
    crit = next(c for c in inp.criteria if c.id == criterion_id)
    born = {q.id: q.born for q in inp.queens}
    recs = [r for r in inp.records if r.criterion == criterion_id]
    ext = [r for r in extra_records if r.criterion == criterion_id] if recs else []
    born.update({q.id: q.born for q in extra_queens})
    scales = {(r.scale_min, r.scale_max) for r in recs}
    groups: dict[tuple, list[RecordIn]] = {}
    for r in recs + ext:
        # une visite = un rucher et une année : des notes de même visit mais de rucher ou d'année
        # différents restent des observations distinctes (groupe de contemporaines sans ambiguïté)
        if r.visit is not None:
            key = (r.queen, "v", r.visit, r.apiary, r.date.year)
        else:
            key = (r.queen, "d", r.date, r.apiary)
        groups.setdefault(key, []).append(r)
    farm_ids = {r.queen for r in recs}
    y, qids, cg, ev, hc, ag = [], [], [], [], [], []
    for key, rs in groups.items():
        first = rs[0]
        y.append(float(np.mean([normalized_value(r, crit.kind) for r in rs])))
        qids.append(key[0])
        cg.append((first.apiary, first.date.year))
        ev.append(first.evaluator if first.evaluator is not None else -1)
        external = key[0] not in farm_ids         # niveaux propres aux notes externes : aucun effet partagé
        hc.append("x" + first.hive_class if external else first.hive_class)
        age = age_class(first.date.year, born.get(key[0]))
        ag.append("x" + age if external else age)
    n = len(y)
    n_farm = sum(1 for q in qids if q in farm_ids)          # observations de l'exploitation : les n_farm premières
    farm = np.arange(n) < n_farm
    fw = np.array([f_workers.get(q, 0.0) for q in qids]) if f_workers is not None else None
    covs = []
    if fw is not None and n > 0:
        covs = _covariate("f_workers", fw, farm) + _covariate("f_workers_ext", fw, ~farm)
    factors_in = [("cg", cg), ("evaluator", ev), ("hive_class", hc), ("age", ag)]
    Xd, names, factors, labels = _design(n, factors_in, covs)
    if n_farm < n:      # a priori des variances : plan de l'exploitation seule, comme sans notes externes
        farm_cov = _covariate("f_workers", fw[:n_farm], farm[:n_farm]) if fw is not None else []
        Xf = _design(n_farm, [(f, v[:n_farm]) for f, v in factors_in], farm_cov)[0]
        farm_X = sp.csr_matrix(Xf)
    else:
        farm_X = None
    recorded = list(dict.fromkeys(qids))
    n_rec: dict[int, int] = {}
```

`blup2/records.py` : remplacer :

```python
            evs[r.queen].add(r.evaluator)
        yrs.setdefault(r.queen, set()).add(r.date.year)
    return CriterionData(
        criterion=crit,
        y=np.array(y, dtype=float),
        queen_ids=qids,
        X=sp.csr_matrix(Xd),
        x_names=names,
        fixed_factors=factors,
```

par :

```python
            evs[r.queen].add(r.evaluator)
        yrs.setdefault(r.queen, set()).add(r.date.year)
    X = sp.csr_matrix(Xd)
    return CriterionData(
        criterion=crit,
        y=np.array(y, dtype=float),
        queen_ids=qids,
        X=X,
        x_names=names,
        fixed_factors=factors,
```

`blup2/records.py` : remplacer :

```python
        n_years={q: len(s) for q, s in yrs.items()},
        n_raw_records=len(recs),
        n_contemporary_groups=len(set(cg)),
        scale_mixed=crit.kind == "note" and len(scales) > 1,
    )

```

par :

```python
        n_years={q: len(s) for q, s in yrs.items()},
        n_raw_records=len(recs),
        n_contemporary_groups=len(set(cg[:n_farm])),
        scale_mixed=crit.kind == "note" and len(scales) > 1,
        n_farm_obs=n_farm,
        farm_X=farm_X if farm_X is not None else X,
    )

```

- [ ] **Step 4 : Vérifier que les tests passent**

Run : `MELLIFERA_CONTACT= .venv/bin/pytest -q -p no:cacheprovider tests/test_records.py | tail -1`

Expected :

```text
18 passed in …
```

- [ ] **Step 5 : Lancer la suite complète**

```bash
MELLIFERA_CONTACT= .venv/bin/pytest -q -p no:cacheprovider | tail -1
```

Expected :

```text
430 passed, 12 deselected in …
```

Puis ruff :

```bash
uvx ruff check --select F,E9 blup2 tests
```

Expected :

```text
All checks passed!
```

- [ ] **Step 6 : Commit**

```bash
git add blup2/records.py tests/test_records.py
git status --short
git commit -q -F - <<'EOF'
Sépare les effets fixes des notes externes dans la préparation d'un critère

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015J8yR1HR7mLSRtwnvEYap6
EOF
git log --oneline -1
```

Expected : `git status --short` ne liste que les 2 fichiers ci-dessus, puis :

```text
<hash> Sépare les effets fixes des notes externes dans la préparation d'un critère
```

---

### Task 3 : Lecture des notes des ancêtres et de leurs contemporaines

Le cœur du volet (spec §4, §5, §6) : `blup2/external_records.py` lit et nettoie les notes, choisit les contemporaines et construit les `RecordIn` externes ; `expand_external` l'appelle après l'ascendance, avec repli `external_records_unavailable` ; le moteur passe reines et notes externes à `build_pedigree` et à `prepare_criterion` ; les modèles calculent a priori, compteurs, REML et `few_records` sur l'exploitation seule. La mini-base des tests gagne les tables de notes.

**Files:**
- Create: `blup2/external_records.py`
- Modify: `blup2/external.py`
- Modify: `blup2/engine.py`
- Modify: `blup2/model_standard.py`
- Modify: `blup2/model_breeder.py`
- Test (modify): `tests/external_base.py`
- Test (create): `tests/test_external_records.py`
- Test (modify): `tests/mellifera/test_mellifera_external_contract.py`
- Test (modify): `tests/test_external_engine.py`

**Interfaces:**
- Consumes :
  - tâche 1 : `STANDARD_TRAITS[c.trait].mellifera`, `CriterionIn.trait`, codes `external_records_added` et `external_records_unavailable` ;
  - tâche 2 : `prepare_criterion(…, extra_queens=, extra_records=)`, `CriterionData.farm_queens`, `prior_X`, `prior_y` ;
  - volet 2 : `expand_external(inp) -> Expansion`, `_ancestor(row, ids) -> (QueenIn, …)`, `_cut_cycles(queens)`, `_Base.conn`, `_Base.row(code)`, `build_pedigree(queens, n_drones, extra_queens=)`, `Pedigree.added_queen_ids`.
- Produces :
  - `blup2.external_records` : `MAX_CONTEMPORARIES = 1000`, `DISCARD_LABELS: dict[str, str]`, dataclass `Notes(cells, contemporaries, discarded, capped)`, `wanted_traits(inp) -> dict[str, int]`, `read_notes(conn, inp, ancestors: dict[str, int]) -> Notes`, `build_records(notes, ids: dict[str, int], inp) -> list[RecordIn]`, `records_message(records, ancestor_ids: set[int], notes) -> Optional[str]` ;
  - `Expansion.contemporaries: list[QueenIn]` et `Expansion.records: list[RecordIn]` ;
  - helpers de test `tests.external_base` : `NOTES_SCHEMA`, `BREEDERS`, `PEERS`, `NOTES`, `make_db(path, rows=LINEAGE, notes=None)`, `inp(queens, depth=7, records=(), criteria=None, **kw)`.

- [ ] **Step 1 : Écrire la mini-base des notes et les tests de lecture, d'isolement et de contrat**

`tests/external_base.py` : remplacer :

```python
"""Mini-base Index Mellifera pour les tests du volet 2 (sous-ensemble du schéma de mellifera/db.py)."""
import sqlite3

```

par :

```python
"""Mini-base Index Mellifera pour les tests des volets 2 et 3 (sous-ensemble du schéma de mellifera/db.py)."""
import sqlite3

```

`tests/external_base.py` : remplacer :

```python
from tests.builders import criterion

# Colonnes lues par blup2/external.py ; le contrat avec la vraie base du crawler est vérifié par
# tests/mellifera/test_mellifera_external_contract.py.
SCHEMA = """
```

par :

```python
from tests.builders import criterion

# Colonnes lues par blup2/external.py et blup2/external_records.py ; le contrat avec la vraie base du crawler est vérifié par
# tests/mellifera/test_mellifera_external_contract.py.
SCHEMA = """
```

`tests/external_base.py` : remplacer :

```python
    id INTEGER PRIMARY KEY AUTOINCREMENT, code TEXT NOT NULL UNIQUE, name TEXT, birth_year INTEGER,
    mother_code TEXT, drone_mother_code TEXT, n_drone_queens INTEGER, mating_code TEXT,
    single_drone INTEGER NOT NULL DEFAULT 0, race_description TEXT);
CREATE INDEX queen_name ON queen (name, birth_year);
CREATE TABLE mating_code (code TEXT PRIMARY KEY, mating_type TEXT NOT NULL, mating_type_override TEXT);
```

par :

```python
    id INTEGER PRIMARY KEY AUTOINCREMENT, code TEXT NOT NULL UNIQUE, name TEXT, birth_year INTEGER,
    mother_code TEXT, drone_mother_code TEXT, n_drone_queens INTEGER, mating_code TEXT,
    single_drone INTEGER NOT NULL DEFAULT 0, race_description TEXT, breeder_id INTEGER);
CREATE INDEX queen_name ON queen (name, birth_year);
CREATE TABLE mating_code (code TEXT PRIMARY KEY, mating_type TEXT NOT NULL, mating_type_override TEXT);
```

`tests/external_base.py` : remplacer :

```python
INSERT INTO mating_code VALUES ('ins', 'insemination', NULL), ('nat', 'natural', NULL), ('xx', 'unknown', NULL);
"""

# Lignée CBP : 0276 ← mère 0100 (× drones 0009 NG) ← mère 0050 ← mère 0010 ; 0277 sœur de 0276.
```

par :

```python
INSERT INTO mating_code VALUES ('ins', 'insemination', NULL), ('nat', 'natural', NULL), ('xx', 'unknown', NULL);
"""

NOTES_SCHEMA = """
CREATE TABLE trait (id INTEGER PRIMARY KEY, name TEXT NOT NULL UNIQUE, label TEXT NOT NULL, kind TEXT NOT NULL,
    scale_min REAL, scale_max REAL, direction INTEGER NOT NULL);
CREATE TABLE criterion (id INTEGER PRIMARY KEY AUTOINCREMENT, mode TEXT NOT NULL, code TEXT NOT NULL,
    trait_id INTEGER REFERENCES trait (id), UNIQUE (mode, code));
CREATE TABLE evaluation (
    id INTEGER PRIMARY KEY AUTOINCREMENT, queen_name_raw TEXT NOT NULL, queen_birth_year_raw INTEGER,
    queen_code TEXT, queen_key TEXT NOT NULL, page_breeder_id INTEGER, tester TEXT NOT NULL,
    tester_breeder_id INTEGER, eval_date TEXT NOT NULL, eval_year INTEGER, mode TEXT NOT NULL,
    criterion_code TEXT NOT NULL, value REAL, value_raw TEXT NOT NULL, source TEXT NOT NULL, url TEXT NOT NULL,
    UNIQUE (queen_key, source, mode, tester, eval_date, criterion_code));
INSERT INTO trait VALUES (1, 'swarming', 'Essaimage', 'note', 0, 6, 1), (2, 'mildness', 'Douceur', 'note', 0, 6, 1),
    (19, 'hyg24h', 'Couvain nettoyé en 24 h (%)', 'percentage', NULL, NULL, 1),
    (21, 'honeyAmt', 'Quantité de miel (kg)', 'harvest', NULL, NULL, 1);
INSERT INTO criterion (mode, code, trait_id) VALUES ('evjp', 'mildness', 2), ('evsp', 'mildness', 2),
    ('radar', 'mildness', 2), ('evjp', 'swarming', 1), ('evjp', 'hyg24h', 19), ('evjp', 'honeyAmt', 21),
    ('evjp', 'comment', NULL);
"""
BREEDERS = {"CBP": 145, "NG": 200}         # éleveur Index Mellifera d'après le trigramme du code

# Lignée CBP : 0276 ← mère 0100 (× drones 0009 NG) ← mère 0050 ← mère 0010 ; 0277 sœur de 0276.
```

`tests/external_base.py` : remplacer :

```python


def make_db(path, rows=LINEAGE):
    conn = sqlite3.connect(path)
    conn.executescript(SCHEMA)
    conn.executemany("INSERT INTO queen (code, birth_year, mother_code, drone_mother_code, name, mating_code) "
                     "VALUES (?, ?, ?, ?, ?, ?)", rows)
    conn.commit()
    conn.close()
```

par :

```python


# Contemporaines : reines notées par les mêmes testeurs, les mêmes années, que des ancêtres de la lignée CBP.
PEERS = [
    ("FR-T-B-0001-2022", 2022, None, None, None, None),
    ("FR-T-B-0002-2022", 2022, None, None, None, None),
    ("FR-T-B-0003-2020", 2020, "FR-CBP-B-0050-2018", None, None, None),        # fille d'un ancêtre
    ("FR-CBP-B-0300-2023", 2023, None, None, None, None),                      # radar du même éleveur (145)
]

# Notes : code, source, mode, critère du site, valeur, date, testeur. Ancêtres 0276, 0100 et 0009 puis
# contemporaines ; 0009 (testeur 301 en 2019) n'a pas de contemporaine.
NOTES = [
    ("FR-CBP-B-0276-2022", "season_table", "evjp", "mildness", 5, "2023-07-06", 300),
    ("FR-CBP-B-0276-2022", "season_table", "evjp", "swarming", 3, "2023-07-06", 300),
    ("FR-CBP-B-0276-2022", "radar", "radar", "mildness", 4, "2023-01-01", None),     # couverte par la saison 2023
    ("FR-CBP-B-0276-2022", "radar", "radar", "mildness", 3, "2024-01-01", None),     # rucher : éleveur CBP (145)
    ("FR-CBP-B-0100-2020", "season_table", "evjp", "mildness", 6, "2021-07-01", 300),
    ("FR-CBP-B-0100-2020", "season_table", "evsp", "mildness", 4, "2021-07-01", 300),  # même clé : moyennée
    ("FR-NG-B-0009-2018", "season_table", "evjp", "mildness", 2, "2019-06-01", 301),
    ("FR-T-B-0001-2022", "season_table", "evjp", "mildness", 3, "2023-07-06", 300),
    ("FR-T-B-0001-2022", "season_table", "evjp", "swarming", 4, "2023-07-06", 300),
    ("FR-T-B-0002-2022", "season_table", "evjp", "mildness", 4, "2023-07-07", 300),
    ("FR-T-B-0003-2020", "season_table", "evjp", "mildness", 2, "2021-07-01", 300),
    ("FR-T-B-0003-2020", "season_table", "evjp", "mildness", 1, "2022-07-01", 300),  # autre année : ignorée
    ("FR-CBP-B-0300-2023", "radar", "radar", "mildness", 5, "2024-01-01", None),
    (None, "season_table", "evjp", "mildness", 4, "2023-07-08", 300),               # reine non résolue
]


def make_db(path, rows=LINEAGE, notes=None):
    """notes : None = base du volet 2 (sans tables de notes) ; sinon lignes au format de NOTES (PEERS ajoutées)."""
    conn = sqlite3.connect(path)
    conn.executescript(SCHEMA)
    conn.executemany("INSERT INTO queen (code, birth_year, mother_code, drone_mother_code, name, mating_code) "
                     "VALUES (?, ?, ?, ?, ?, ?)", list(rows) + (PEERS if notes is not None else []))
    for trigram, breeder in BREEDERS.items():
        conn.execute("UPDATE queen SET breeder_id = ? WHERE code LIKE ?", (breeder, f"FR-{trigram}-%"))
    if notes is not None:
        conn.executescript(NOTES_SCHEMA)
        conn.executemany(
            "INSERT INTO evaluation (queen_name_raw, queen_code, queen_key, tester, tester_breeder_id, eval_date, "
            "mode, criterion_code, value, value_raw, source, url) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, '/x/')",
            [(code or "X", code, code or f"name:X|{k}", "" if tester is None else f"T{tester}", tester, date, mode,
              crit, value, "" if value is None else str(value), source)
             for k, (code, source, mode, crit, value, date, tester) in enumerate(notes)])
    conn.commit()
    conn.close()
```

`tests/external_base.py` : remplacer :

```python


def inp(queens, depth=7, records=()):
    return BlupV2Input(exploitation=1, criteria=[criterion(1)], queens=queens, records=list(records),
                       options=OptionsIn(external_pedigree_depth=depth))


```

par :

```python


def inp(queens, depth=7, records=(), criteria=None, **kw):
    return BlupV2Input(exploitation=1, criteria=criteria or [criterion(1)], queens=queens, records=list(records),
                       options=OptionsIn(external_pedigree_depth=depth), **kw)


```

Créer `tests/test_external_records.py` :

```python
"""Volet 3 Index Mellifera : notes des ancêtres résolus et de leurs contemporaines (spec volet 3 §4, §6)."""
import datetime as dt
import sqlite3

import pytest

import blup2.external_records as external_records
from blup2.external import expand_external
from tests.builders import criterion, queen, record
from tests.external_base import NOTES, inp, make_db, warning

D = dt.date


@pytest.fixture
def db(tmp_path):
    return make_db(tmp_path / "mellifera.sqlite", notes=NOTES)


def farm(criteria, apiary=1, evaluator=1, depth=7):
    """Reine 1 fille de FR-CBP-B-0276-2022, notée sur chaque critère. Ancêtres : 0276 = -1, 0100 = -2,
    0009 = -3, 0050 = -4, 0010 = -5 ; contemporaines ensuite, par code croissant."""
    recs = [record(1, 0.5, criterion=c.id, apiary=apiary, evaluator=evaluator) for c in criteria]
    return inp([queen(1, 2024, mother_ref_ext="FR-CBP-B-0276-2022")], records=recs, criteria=criteria, depth=depth)


def notes(exp):
    return sorted((r.queen, r.criterion, round(r.value, 6), r.date, r.apiary, r.evaluator) for r in exp.records)


# ------------------------------------------------------------------ rattachement et nettoyage
def test_notes_of_ancestors_and_contemporaries_are_attached_trait_by_trait(db):
    exp = expand_external(farm([criterion(1, trait="gentleness"), criterion(2, trait="swarming"), criterion(3)]), db)
    r6 = lambda x: round(x, 6)                                    # noqa: E731
    assert notes(exp) == [
        (-9, 1, r6(2 / 6), D(2021, 7, 1), -3, -3),               # FR-T-B-0003-2020, contemporaine de 0100
        (-8, 1, r6(4 / 6), D(2023, 7, 7), -3, -3),               # FR-T-B-0002-2022
        (-7, 1, r6(3 / 6), D(2023, 7, 6), -3, -3),               # FR-T-B-0001-2022
        (-7, 2, r6(4 / 6), D(2023, 7, 6), -3, -3),
        (-6, 1, r6(5 / 6), D(2024, 1, 1), -2, -2),               # FR-CBP-B-0300-2023, radar de l'éleveur 145
        (-2, 1, r6(5 / 6), D(2021, 7, 1), -3, -3),               # 0100 : (6 + 4) / 2, deux modes moyennés
        (-1, 1, r6(3 / 6), D(2024, 1, 1), -2, -2),               # 0276, radar 2024
        (-1, 1, r6(5 / 6), D(2023, 7, 6), -3, -3),               # 0276, saison 2023 (le radar 2023 est ignoré)
        (-1, 2, r6(3 / 6), D(2023, 7, 6), -3, -3),
    ]
    assert all((r.scale_min, r.scale_max, r.hive_class) == (0.0, 1.0, "hive") for r in exp.records)
    assert warning(exp, "external_records_added") == [
        "9 note(s) Index Mellifera ajoutée(s) sur 2 critère(s) : 2 ancêtre(s) et 4 contemporaine(s) ; "
        "écartées : 1 sans fiche de reine, 1 sans contemporaine"]
    assert warning(exp, "external_ref_resolved") == [   # l'ascendance du volet 2 est inchangée
        "1 référence(s) Index Mellifera rattachée(s), 5 ancêtre(s) ajouté(s) au pedigree"]


def test_contemporaries_are_built_like_the_ancestors(db):
    exp = expand_external(farm([criterion(1, trait="gentleness")]), db)
    assert [q.id for q in exp.contemporaries] == [-6, -7, -8, -9]
    assert {q.id: q.born for q in exp.contemporaries} == {-6: 2023, -7: 2022, -8: 2022, -9: 2020}
    assert {q.id: q.mother for q in exp.contemporaries}[-9] == -4      # fille de 0050, ancêtre résolu
    assert all(not q.available_as_mother and not q.available_as_drone_source for q in exp.contemporaries)


def test_each_discard_reason_is_counted(tmp_path):
    rows = [("FR-CBP-B-0276-2022", "season_table", "evjp", "mildness", 5, "2023-07-06", 300),   # gardée
            ("FR-CBP-B-0276-2022", "season_table", "evjp", "hyg24h", 80, "2023-07-06", 300),    # gardée : 0,8
            ("FR-T-B-0001-2022", "season_table", "evjp", "mildness", 3, "2023-07-06", 300),     # contemporaine
            ("FR-T-B-0001-2022", "season_table", "evjp", "hyg24h", 50, "2023-07-06", 300),
            ("FR-CBP-B-0276-2022", "season_table", "evjp", "honeyAmt", 30, "2023-07-06", 300),
            ("FR-CBP-B-0276-2022", "season_table", "evsp", "mildness", 4, "", 300),
            ("FR-CBP-B-0276-2022", "season_table", "evsp", "mildness", 7, "2023-07-07", 300),
            ("FR-CBP-B-0276-2022", "season_table", "evjp", "hyg24h", 150, "2023-07-07", 300),
            ("FR-CBP-B-0276-2022", "season_table", "evjp", "mildness", None, "2023-07-08", 300),
            ("FR-CBP-B-0276-2022", "season_table", "evjp", "mildness", 4, "2023-07-09", None),
            ("FR-CBP-B-0276-2022", "season_table", "evjp", "mildness", 4, "2021-07-09", 300)]
    db = make_db(tmp_path / "m.sqlite", notes=rows)
    cs = [criterion(1, trait="gentleness"), criterion(2, trait="hygienic_24h"), criterion(3, trait="honey_yield")]
    exp = expand_external(farm(cs), db)
    assert notes(exp) == [(-6, 1, 0.5, D(2023, 7, 6), -2, -2), (-6, 2, 0.5, D(2023, 7, 6), -2, -2),
                          (-1, 1, round(5 / 6, 6), D(2023, 7, 6), -2, -2), (-1, 2, 0.8, D(2023, 7, 6), -2, -2)]
    assert warning(exp, "external_records_added") == [
        "4 note(s) Index Mellifera ajoutée(s) sur 2 critère(s) : 1 ancêtre(s) et 1 contemporaine(s) ; écartées : "
        "1 en kg, 1 sans date, 2 hors échelle, 1 non numériques, 1 sans rucher, 1 antérieures à la naissance"]


def test_criterion_direction_opposite_to_the_site_inverts_the_note(db):
    same = {(q, c, d): v for q, c, v, d, *_ in notes(expand_external(farm([criterion(1, trait="gentleness")]), db))}
    flip = {(q, c, d): v for q, c, v, d, *_ in notes(
        expand_external(farm([criterion(1, trait="gentleness", direction=-1)]), db))}
    assert flip.keys() == same.keys() and all(flip[k] == pytest.approx(1 - same[k], abs=1e-6) for k in same)


def test_identifiers_stay_below_those_of_the_farm(db):
    exp = expand_external(farm([criterion(1, trait="gentleness")], apiary=-7, evaluator=-9), db)
    assert sorted({r.apiary for r in exp.records}) == [-11, -10]     # saison du testeur 300, radar de l'éleveur 145
    assert all(r.apiary == r.evaluator for r in exp.records)


def test_contemporaries_are_capped(db, monkeypatch):
    monkeypatch.setattr(external_records, "MAX_CONTEMPORARIES", 2)
    exp = expand_external(farm([criterion(1, trait="gentleness"), criterion(2, trait="swarming")]), db)
    assert [q.id for q in exp.contemporaries] == [-6, -7]             # 0300 et T-0001, par code croissant
    assert warning(exp, "external_records_added") == [
        "6 note(s) Index Mellifera ajoutée(s) sur 2 critère(s) : 1 ancêtre(s) et 2 contemporaine(s) ; écartées : "
        "1 sans fiche de reine, 2 sans contemporaine ; plafond de 2 contemporaines atteint"]


# ------------------------------------------------------------------ rien à ajouter
@pytest.mark.parametrize("crit", [
    criterion(1),                                   # critère sans caractère
    criterion(1, trait="colony_strength"),          # caractère sans nom Index Mellifera (BeeBreed seul)
    criterion(1, trait="fertility"),                # nom Index Mellifera absent de cette base
])
def test_nothing_is_added_without_a_known_trait(db, crit):
    exp = expand_external(farm([crit]), db)
    assert exp.records == [] and exp.contemporaries == [] and warning(exp, "external_records_added") == []
    assert len(exp.ancestors) == 5


def test_criterion_without_farm_records_receives_nothing(db):
    src = inp([queen(1, 2024, mother_ref_ext="FR-CBP-B-0276-2022")], records=[record(1, 0.5, criterion=1)],
              criteria=[criterion(1), criterion(2, trait="gentleness")])
    exp = expand_external(src, db)
    assert exp.records == [] and warning(exp, "external_records_added") == []


def test_depth_bounds_the_ancestors_whose_notes_are_read(db):
    exp = expand_external(farm([criterion(1, trait="gentleness")], depth=1), db)
    assert [a.id for a in exp.ancestors] == [-1]
    assert {r.queen for r in exp.records} == {-1, -2, -3, -4}         # 0276 et ses contemporaines 0300, T-0001, T-0002


# ------------------------------------------------------------------ erreurs
def test_base_without_note_tables_keeps_the_ancestry(tmp_path):
    db = make_db(tmp_path / "m.sqlite")                     # base du volet 2 : ni trait, ni evaluation
    exp = expand_external(farm([criterion(1, trait="gentleness")]), db)
    assert len(exp.ancestors) == 5 and exp.records == [] and exp.contemporaries == []
    assert warning(exp, "external_records_unavailable") == [
        "notes Index Mellifera illisibles : calcul avec l'ascendance mais sans leurs notes"]
    assert warning(exp, "external_pedigree_unavailable") == []


def test_base_without_note_tables_is_silent_when_no_trait_is_declared(tmp_path):
    exp = expand_external(farm([criterion(1)]), make_db(tmp_path / "m.sqlite"))
    assert [w.code for w in exp.warnings] == ["external_ref_resolved"]


def test_unexpected_values_in_the_base_are_discarded(tmp_path):
    db = make_db(tmp_path / "m.sqlite", notes=[
        ("FR-CBP-B-0276-2022", "season_table", "evjp", "mildness", 5, "2023-07-06", 300)])
    conn = sqlite3.connect(db)
    conn.executescript("UPDATE evaluation SET value = 'cinq'; UPDATE trait SET scale_max = 'x' WHERE name = 'swarming';"
                       "INSERT INTO evaluation (queen_name_raw, queen_code, queen_key, tester, tester_breeder_id, "
                       "eval_date, mode, criterion_code, value, value_raw, source, url) VALUES ('a', "
                       "'FR-CBP-B-0276-2022', 'k2', 'T', 300, '2023-13-45', 'evjp', 'swarming', 3, '3', "
                       "'season_table', '/x/'), ('a', 'FR-CBP-B-0276-2022', 'k3', 'T', 300, '2023-07-06', 'evjp', "
                       "'swarming', 3, '3', 'season_table', '/x/'), ('a', 'FR-CBP-B-0276-2022', 'k4', 'T', 'abc', "
                       "'2023-07-06', 'evjp', 'mildness', 3, '3', 'season_table', '/x/');")
    conn.close()
    exp = expand_external(farm([criterion(1, trait="gentleness"), criterion(2, trait="swarming")]), db)
    assert exp.records == []
    assert warning(exp, "external_records_added") == [
        "0 note(s) Index Mellifera ajoutée(s) sur 0 critère(s) : 0 ancêtre(s) et 0 contemporaine(s) ; écartées : "
        "1 sans date, 1 hors échelle, 1 non numériques, 1 sans rucher"]


# ------------------------------------------------------------------ points de relecture
def test_negative_farm_ids_keep_external_ids_below(db):
    src = inp([queen(-20, 2024, mother_ref_ext="FR-CBP-B-0276-2022")],
              records=[record(-20, 0.5, apiary=-30, evaluator=-40)], criteria=[criterion(1, trait="gentleness")])
    exp = expand_external(src, db)
    ext_ids = [q.id for q in exp.ancestors + exp.contemporaries]
    assert len(set(ext_ids)) == len(ext_ids) == 9 and max(ext_ids) < -20
    assert {r.queen for r in exp.records} <= set(ext_ids)
    assert max(r.apiary for r in exp.records) < -40


def test_cycles_and_self_references_between_contemporaries_are_cut(db):
    from blup2.pedigree import find_cycle
    conn = sqlite3.connect(db)
    conn.executescript("UPDATE queen SET mother_code = 'FR-T-B-0002-2022' WHERE code = 'FR-T-B-0001-2022';"
                       "UPDATE queen SET mother_code = 'FR-T-B-0001-2022' WHERE code = 'FR-T-B-0002-2022';"
                       "UPDATE queen SET mother_code = code WHERE code = 'FR-CBP-B-0300-2023';")
    conn.commit()
    conn.close()
    exp = expand_external(farm([criterion(1, trait="gentleness")]), db)
    queens = exp.ancestors + exp.contemporaries
    assert len(exp.contemporaries) == 4 and not find_cycle(queens)
    assert all(q.mother != q.id for q in queens)
    assert {q.id: q.mother for q in exp.contemporaries}[-9] == -4      # lien valide conservé


def test_depth_zero_disables_the_notes_too(db):
    exp = expand_external(farm([criterion(1, trait="gentleness")], depth=0), db)
    assert exp.ancestors == [] and exp.contemporaries == [] and exp.records == [] and exp.warnings == []
```

`tests/mellifera/test_mellifera_external_contract.py` : remplacer :

```python
    writer.close()

```

par :

```python
    writer.close()


# ------------------------------------------------------------------ volet 3 : notes et contemporaines
def test_engine_reads_the_notes_of_the_crawler_database(tmp_path):
    """Requêtes de blup2/external_records.py sur une base créée par mellifera/db.py et peuplée par l'ingestion."""
    from blup2.schemas import RecordIn
    from mellifera.urls import season_path
    from tests.mellifera.fake_site import E, season_page

    conn = crawler_db(tmp_path / "m.sqlite")
    radar = Q(DAUGHTER.code, 2022, 145, mother=GRAND.code, drones="FR-NG-B-0009-2018", number=8,
              radar={2024: [4, 3, 5, 2, 1]})                       # douceur 3 en 2024
    peers = [Q("FR-T-B-0001-2022", 2022, 300), Q("FR-CBP-B-0300-2023", 2023, 145, radar={2024: [1, 6, 1, 1, 1]})]
    for q in (radar, *peers):
        ingest(conn, queen_path(q.code), queen_page(q), SEEN)
    ingest(conn, season_path(300, 2023), season_page(
        [E(DAUGHTER.code, 2022, {"mildness": "5"}, tester_id=300, date="July 6, 2023"),
         E(peers[0].code, 2022, {"mildness": "3"}, tester_id=300, date="July 6, 2023")]), SEEN)
    conn.commit()
    conn.close()
    src = beekube_input().model_copy(update={
        "criteria": [CriterionIn(id=1, trait="gentleness")],
        "records": [RecordIn(queen=1, criterion=1, value=0.5, date="2025-06-01", apiary=1)]})
    exp = expand_external(BlupV2Input.model_validate(src.model_dump()), str(tmp_path / "m.sqlite"))
    assert [w.message for w in exp.warnings if w.code == "external_records_added"] == [
        "4 note(s) Index Mellifera ajoutée(s) sur 1 critère(s) : 1 ancêtre(s) et 2 contemporaine(s)"]
    assert [(a.id, a.born) for a in exp.ancestors] == [(-1, 2020), (-2, 2022)]      # GRAND, DAUGHTER
    assert [(q.id, q.born) for q in exp.contemporaries] == [(-3, 2023), (-4, 2022)]
    assert sorted((r.queen, round(r.value, 4), r.date.isoformat()) for r in exp.records) == [
        (-4, 0.5, "2023-07-06"), (-3, 1.0, "2024-01-01"), (-2, 0.5, "2024-01-01"), (-2, 0.8333, "2023-07-06")]


def test_notes_are_read_while_the_crawler_writes(tmp_path):
    """Lecture des notes pendant une transaction d'écriture du crawler (WAL) : dernier état validé."""
    from blup2.schemas import RecordIn
    from mellifera.urls import season_path
    from tests.mellifera.fake_site import E, season_page

    writer = crawler_db(tmp_path / "m.sqlite")
    ingest(writer, queen_path("FR-T-B-0001-2022"), queen_page(Q("FR-T-B-0001-2022", 2022, 300)), SEEN)
    ingest(writer, season_path(300, 2023), season_page(
        [E(DAUGHTER.code, 2022, {"mildness": "5"}, tester_id=300, date="July 6, 2023"),
         E("FR-T-B-0001-2022", 2022, {"mildness": "3"}, tester_id=300, date="July 6, 2023")]), SEEN)
    writer.commit()
    writer.execute("BEGIN IMMEDIATE")
    writer.execute("UPDATE evaluation SET value = 0")
    src = beekube_input().model_copy(update={
        "criteria": [CriterionIn(id=1, trait="gentleness")],
        "records": [RecordIn(queen=1, criterion=1, value=0.5, date="2025-06-01", apiary=1)]})
    exp = expand_external(BlupV2Input.model_validate(src.model_dump()), str(tmp_path / "m.sqlite"))
    assert sorted(round(r.value, 4) for r in exp.records) == [0.5, 0.8333]
    writer.rollback()
    writer.close()

```

`tests/test_external_engine.py` : remplacer :

```python
from blup2.schemas import MatingPairIn, MatingV2Input
from tests.builders import criterion, mating, queen, record
from tests.external_base import CYCLES, inp, make_db, warning


```

par :

```python
from blup2.schemas import MatingPairIn, MatingV2Input
from tests.builders import criterion, mating, queen, record
from tests.external_base import CYCLES, NOTES, inp, make_db, warning


```

`tests/test_external_engine.py` : remplacer :

```python
    assert warning(out, "external_ref_resolved") == []

```

par :

```python
    assert warning(out, "external_ref_resolved") == []


# ------------------------------------------------------------------ volet 3 : notes des ancêtres
GENTLE = dict(trait="gentleness", h2=0.3, repeatability=0.4)          # variances fixées : comparaisons nettes


def run(monkeypatch, db, criteria, refs=True, profile="standard"):
    qs, recs = family(refs=refs)
    if db is None:
        monkeypatch.delenv("MELLIFERA_DB", raising=False)
    else:
        monkeypatch.setenv("MELLIFERA_DB", db)
    return run_blup(inp(qs, records=recs, criteria=criteria, profile=profile))


def dump(out):
    return [q.model_dump() for q in out.queens], out.parameters, out.base, out.matings, out.trend, out.connectivity


@pytest.mark.parametrize("case", ["no_trait", "no_reference", "no_base"])
def test_outputs_are_bit_identical_when_notes_cannot_apply(tmp_path, monkeypatch, case):
    db = make_db(tmp_path / "m.sqlite", notes=NOTES)
    plain = make_db(tmp_path / "plain.sqlite")                 # base du volet 2, sans notes
    crit = [criterion(1, h2=0.3, repeatability=0.4, **({} if case == "no_trait" else {"trait": "gentleness"}))]
    refs = case != "no_reference"
    out = run(monkeypatch, None if case == "no_base" else db, crit, refs=refs)
    ref = run(monkeypatch, None if case == "no_base" else plain, [criterion(1, h2=0.3, repeatability=0.4)], refs=refs)
    assert dump(out) == dump(ref)
    assert warning(out, "external_records_added") == []


def daughters(monkeypatch, db, criteria):
    """Six reines nées en 2023, notées en 2024 dans un même rucher : 1 à 3 filles de FR-CBP-B-0276-2022
    (notée avec ses contemporaines), 4 à 6 sans lien avec Index Mellifera."""
    qs = [queen(k, 2023, mother_ref_ext="FR-CBP-B-0276-2022" if k <= 3 else None) for k in range(1, 7)]
    recs = [record(k, v, date="2024-06-01") for k, v in zip(range(1, 7), (0.7, 0.6, 0.8, 0.4, 0.5, 0.3))]
    monkeypatch.setenv("MELLIFERA_DB", db)
    return run_blup(inp(qs, records=recs, criteria=criteria))


def test_notes_of_the_ancestors_raise_the_reliability_of_their_daughters(tmp_path, monkeypatch):
    db = make_db(tmp_path / "m.sqlite", notes=NOTES)
    with_notes = daughters(monkeypatch, db, [criterion(1, **GENTLE)])
    without = daughters(monkeypatch, db, [criterion(1, h2=0.3, repeatability=0.4)])
    a = {o.id: o.indexes["1"] for o in with_notes.queens}
    b = {o.id: o.indexes["1"] for o in without.queens}
    for q in (1, 2, 3):
        assert a[q].reliability > b[q].reliability
    assert warning(with_notes, "external_records_added") == [
        "7 note(s) Index Mellifera ajoutée(s) sur 1 critère(s) : 2 ancêtre(s) et 4 contemporaine(s) ; "
        "écartées : 1 sans fiche de reine, 1 sans contemporaine"]


def test_unrelated_queens_keep_their_ebv_and_statistics_stay_those_of_the_farm(tmp_path, monkeypatch):
    db = make_db(tmp_path / "m.sqlite", notes=NOTES)
    with_notes = run(monkeypatch, db, [criterion(1, **GENTLE)])
    without = run(monkeypatch, db, [criterion(1, h2=0.3, repeatability=0.4)])
    a = {o.id: o.indexes["1"] for o in with_notes.queens}
    b = {o.id: o.indexes["1"] for o in without.queens}
    for q in (1, 2, 3):                 # sans lien avec les lignées Index Mellifera ni groupe commun avec 4 et 5
        assert a[q].ebv == pytest.approx(b[q].ebv, abs=1e-9)
        assert a[q].reliability == pytest.approx(b[q].reliability, abs=1e-9)
    for q in range(1, 6):
        assert (a[q].n_records, a[q].n_evaluators, a[q].n_years) == (b[q].n_records, b[q].n_evaluators, b[q].n_years)
    p, p0 = with_notes.parameters["1"], without.parameters["1"]
    assert (p.n_records, p.n_queens, p.n_contemporary_groups, p.h2) == (p0.n_records, p0.n_queens,
                                                                         p0.n_contemporary_groups, p0.h2)
    assert [q.id for q in with_notes.queens] == [1, 2, 3, 4, 5] and with_notes.base == without.base
    assert with_notes.connectivity == without.connectivity
    assert all(i > 0 for m in with_notes.matings for i in (m.mother, m.drone_source))
    assert all(p.birth_year in (2021, 2022, 2023) for pts in with_notes.trend.values() for p in pts)


def test_breeder_profile_counts_only_the_farm(tmp_path, monkeypatch):
    db = make_db(tmp_path / "m.sqlite", notes=NOTES)
    out = run(monkeypatch, db, [criterion(1, **GENTLE)], profile="breeder")
    ref = run(monkeypatch, db, [criterion(1, h2=0.3, repeatability=0.4)], profile="breeder")
    p, p0 = out.parameters["1"], ref.parameters["1"]
    assert (p.n_queens, p.n_records, p.known_mating_share) == (p0.n_queens, p0.n_records, p0.known_mating_share)
    assert warning(out, "external_records_added")


def test_notes_unreadable_after_the_ancestry_keep_the_ancestry(tmp_path, monkeypatch):
    plain = make_db(tmp_path / "plain.sqlite")
    out = run(monkeypatch, plain, [criterion(1, **GENTLE)])
    ref = run(monkeypatch, plain, [criterion(1, h2=0.3, repeatability=0.4)])
    assert warning(out, "external_records_unavailable") == [
        "notes Index Mellifera illisibles : calcul avec l'ascendance mais sans leurs notes"]
    assert warning(out, "external_ref_resolved") == warning(ref, "external_ref_resolved") != []
    assert [q.model_dump() for q in out.queens] == [q.model_dump() for q in ref.queens]


def test_pedigree_fallback_also_drops_the_notes(tmp_path, monkeypatch):
    import blup2.engine as engine
    real = engine.build_pedigree

    def flaky(queens, n, extra_queens=()):
        if list(extra_queens):
            raise ValueError("pedigree cyclique : reines concernées -1, 5")
        return real(queens, n, extra_queens)
    monkeypatch.setattr(engine, "build_pedigree", flaky)
    db = make_db(tmp_path / "m.sqlite", notes=NOTES)
    out = run(monkeypatch, db, [criterion(1, **GENTLE)])
    assert [q.id for q in out.queens] == [1, 2, 3, 4, 5]
    assert all(q.indexes["1"].status != "failed" for q in out.queens)
    assert warning(out, "external_records_added") == [] and warning(out, "criterion_failed") == []
    assert warning(out, "external_pedigree_unavailable")


def test_matings_route_accepts_the_notes(tmp_path, monkeypatch):
    monkeypatch.setenv("MELLIFERA_DB", make_db(tmp_path / "m.sqlite", notes=NOTES))
    qs, recs = family()
    pair = MatingV2Input(exploitation=1, criteria=[criterion(1, **GENTLE)], queens=qs, records=recs,
                         pairs=[MatingPairIn(mother=4, drone_source=5)])
    out = run_matings(pair)
    assert out.matings[0].expected_indexes["1"] is not None
    assert warning(out, "external_records_added")



def test_reml_with_external_notes_counts_only_the_farm(tmp_path, monkeypatch):
    """Spec volet 3 §5 : le REML voit toutes les données, les compteurs publiés restent ceux de l'exploitation."""
    from blup2.schemas import BlupV2Input, OptionsIn
    monkeypatch.setenv("MELLIFERA_DB", make_db(tmp_path / "m.sqlite", notes=NOTES))
    qs, recs = family()
    out = run_blup(BlupV2Input(exploitation=1, criteria=[criterion(1, trait="gentleness")], queens=qs, records=recs,
                               options=OptionsIn(estimate_variances="always")))
    p = out.parameters["1"]
    assert warning(out, "criterion_failed") == [] and warning(out, "external_records_added")
    assert (p.n_records, p.n_queens) == (5, 5) and 0.0 <= p.h2 <= 1.0
    assert all(q.indexes["1"].status != "failed" for q in out.queens)

```

- [ ] **Step 2 : Vérifier l'échec**

Run : `MELLIFERA_CONTACT= .venv/bin/pytest -q -p no:cacheprovider -rfE --tb=no --continue-on-collection-errors tests/test_external_records.py tests/test_external_engine.py tests/mellifera/test_mellifera_external_contract.py | tail -30`

`tests/test_external_records.py` ne s'importe pas (`blup2.external_records` n'existe pas) ; les tests du moteur et du contrat qui attendent des notes échouent. Les tests de sorties identiques passent déjà : ils fixent le comportement à garder.

Expected :

```text
.............F.FF.FF..FF                                                 [100%]
=========================== short test summary info ============================
FAILED tests/test_external_engine.py::test_notes_of_the_ancestors_raise_the_reliability_of_their_daughters
FAILED tests/test_external_engine.py::test_breeder_profile_counts_only_the_farm
FAILED tests/test_external_engine.py::test_notes_unreadable_after_the_ancestry_keep_the_ancestry
FAILED tests/test_external_engine.py::test_matings_route_accepts_the_notes - ...
FAILED tests/test_external_engine.py::test_reml_with_external_notes_counts_only_the_farm
FAILED tests/mellifera/test_mellifera_external_contract.py::test_engine_reads_the_notes_of_the_crawler_database
FAILED tests/mellifera/test_mellifera_external_contract.py::test_notes_are_read_while_the_crawler_writes
ERROR tests/test_external_records.py
7 failed, 17 passed, 1 error in …
```

- [ ] **Step 3 : Écrire `external_records.py` et le câbler dans `external.py`, le moteur et les modèles**

Créer `blup2/external_records.py` :

```python
"""Notes Index Mellifera des ancêtres résolus et de leurs contemporaines (spec volet 3 §4), lues dans la base
du crawler.

Un BLUP n'apprend d'une note que par contraste avec les autres notes de son groupe (testeur × année) : la note
d'un ancêtre seul dans son groupe n'apporte rien. On lit donc aussi les notes des autres reines du même groupe
(même testeur, ou même éleveur pour le radar, même année, même caractère) : ces contemporaines entrent au
pedigree comme les ancêtres, sans jamais sortir. Inclure les reines parmi lesquelles les reproductrices ont été
choisies rend aussi l'estimation peu sensible à cette sélection.

Seuls les critères qui déclarent un caractère standard (`CriterionIn.trait`) ayant un nom Index Mellifera, et
qui ont au moins une note de l'exploitation, reçoivent des notes. Nettoyage comme `mellifera/export.py::_records`
(sans en importer le code : le moteur ne dépend pas du paquet mellifera), chaque motif d'écart compté. Les notes
retenues sont ramenées sur [0, 1] dans le sens du critère ; rucher et évaluateur = un entier négatif par
testeur (ou éleveur, pour le radar).
"""
from __future__ import annotations

import datetime as dt
import json
import math
import re
import sqlite3
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Optional, Union

from .schemas import BORN_MIN, BlupV2Input, RecordIn
from .traits import STANDARD_TRAITS

MAX_CONTEMPORARIES = 1000   # borne le coût : environ quatre équations par contemporaine notée
DISCARD_LABELS = {          # motif -> libellé de l'avertissement external_records_added, dans cet ordre
    "harvest": "en kg",
    "no_date": "sans date",
    "out_of_scale": "hors échelle",
    "non_numeric": "non numériques",
    "no_apiary": "sans rucher",
    "before_birth": "antérieures à la naissance",
    "no_queen": "sans fiche de reine",
    "alone": "sans contemporaine",
}
_SELECT = (
    "SELECT e.queen_code, e.source, e.tester_breeder_id, e.eval_date, e.value, t.name AS trait, t.kind, "
    "t.scale_min, t.scale_max, t.direction, q.breeder_id, q.birth_year "
    "FROM evaluation e JOIN criterion c ON c.mode = e.mode AND c.code = e.criterion_code "
    "JOIN trait t ON t.id = c.trait_id LEFT JOIN queen q ON q.code = e.queen_code "
    "WHERE e.source IN ('season_table', 'radar') AND t.name IN (SELECT value FROM json_each(?)) AND ")
_ANCESTORS = _SELECT + "e.queen_code IN (SELECT value FROM json_each(?)) ORDER BY e.id"
_GROUPS = _SELECT + (
    "((e.source = 'season_table' AND e.tester_breeder_id IN (SELECT value FROM json_each(?))) "
    "OR (e.source = 'radar' AND q.breeder_id IN (SELECT value FROM json_each(?)))) ORDER BY e.id")
_ISO = re.compile(r"\d{4}-\d{2}-\d{2}")

Cell = tuple            # (code, critère, date, source, rucher Index Mellifera)


@dataclass
class Notes:
    cells: dict[Cell, list[float]] = field(default_factory=dict)    # valeurs sur [0, 1], sens du critère
    contemporaries: list[str] = field(default_factory=list)          # codes, ordre croissant
    discarded: Counter = field(default_factory=Counter)
    capped: bool = False


def wanted_traits(inp: BlupV2Input) -> dict[str, int]:
    """Nom Index Mellifera -> critère de l'exploitation qui le porte (critères notés dans l'exploitation)."""
    noted = {r.criterion for r in inp.records}
    out = {}
    for c in inp.criteria:
        name = STANDARD_TRAITS[c.trait].mellifera if c.trait is not None else None
        if name is not None and c.id in noted:
            out[name] = c.id
    return out


def _number(v) -> Optional[float]:
    if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
        return None
    return float(v)


def _date(v) -> Optional[dt.date]:
    if not isinstance(v, str) or not _ISO.fullmatch(v):
        return None
    try:
        return dt.date.fromisoformat(v)
    except ValueError:
        return None


def _valid_year(y) -> bool:
    return isinstance(y, int) and not isinstance(y, bool) and BORN_MIN <= y <= dt.date.today().year + 1


def _clean(e, wanted: dict[str, int], direction: dict[int, int]) -> Union[str, tuple]:
    """Motif d'écart, ou (critère, date, source, rucher, valeur sur [0, 1] dans le sens du critère)."""
    if e["kind"] == "harvest":
        return "harvest"
    value = _number(e["value"])
    if value is None:
        return "non_numeric"
    date = _date(e["eval_date"])
    if date is None:
        return "no_date"
    if e["kind"] == "percentage":
        value /= 100
    if e["kind"] == "note":
        lo, hi = _number(e["scale_min"]), _number(e["scale_max"])
        if lo is None or hi is None or not hi > lo:
            return "out_of_scale"
    else:
        lo, hi = 0.0, 1.0
    if not lo <= value <= hi:
        return "out_of_scale"
    apiary = e["tester_breeder_id"] if e["source"] == "season_table" else e["breeder_id"]
    if isinstance(apiary, bool) or not isinstance(apiary, int):
        return "no_apiary"
    v = (value - lo) / (hi - lo)
    cid = wanted[e["trait"]]
    if (-1 if e["direction"] == -1 else 1) != direction[cid]:
        v = 1.0 - v                     # note du site dans le sens contraire du critère de l'exploitation
    return cid, date, e["source"], apiary, v


def read_notes(conn: sqlite3.Connection, inp: BlupV2Input, ancestors: dict[str, int]) -> Notes:
    """ancestors : code Index Mellifera -> année de naissance (ancêtres résolus)."""
    out = Notes()
    wanted = wanted_traits(inp)
    if not wanted or not ancestors:
        return out
    direction = {c.id: c.direction for c in inp.criteria}
    names = json.dumps(sorted(wanted))
    cells: dict[Cell, list[float]] = defaultdict(list)
    for e in conn.execute(_ANCESTORS, (names, json.dumps(sorted(ancestors)))):
        res = _clean(e, wanted, direction)
        if isinstance(res, tuple) and res[1].year < ancestors[e["queen_code"]]:
            res = "before_birth"
        if isinstance(res, str):
            out.discarded[res] += 1
            continue
        cells[(e["queen_code"], *res[:4])].append(res[4])
    groups = {(k[1], k[3], k[4], k[2].year) for k in cells}         # (critère, source, rucher, année)
    if groups:
        testers = sorted({g[2] for g in groups if g[1] == "season_table"})
        breeders = sorted({g[2] for g in groups if g[1] == "radar"})
        for e in conn.execute(_GROUPS, (names, json.dumps(testers), json.dumps(breeders))):
            code = e["queen_code"]
            if code in ancestors:
                continue
            res = _clean(e, wanted, direction)
            if isinstance(res, str) or (res[0], res[2], res[3], res[1].year) not in groups:
                continue                # hors des groupes des ancêtres : ignorée sans être comptée
            if code is None or not _valid_year(e["birth_year"]):
                out.discarded["no_queen"] += 1
            elif res[1].year < e["birth_year"]:
                out.discarded["before_birth"] += 1
            else:
                cells[(code, *res[:4])].append(res[4])
    season = {(k[0], k[1], k[2].year) for k in cells if k[3] == "season_table"}
    cells = {k: v for k, v in cells.items()       # le tableau de saison couvre la même reine, le même critère, l'année
             if not (k[3] == "radar" and (k[0], k[1], k[2].year) in season)}
    others = sorted({k[0] for k in cells if k[0] not in ancestors})
    if len(others) > MAX_CONTEMPORARIES:
        out.capped = True
        kept = set(others[:MAX_CONTEMPORARIES])
        cells = {k: v for k, v in cells.items() if k[0] in ancestors or k[0] in kept}
    members: dict[tuple, set] = defaultdict(set)
    for k in cells:
        members[(k[1], k[3], k[4], k[2].year)].add(k[0])
    live = {g for g, qs in members.items() if len(qs) >= 2 and qs & ancestors.keys()}
    for k in list(cells):
        if (k[1], k[3], k[4], k[2].year) not in live:
            if k[0] in ancestors:
                out.discarded["alone"] += 1
            del cells[k]
    out.cells = cells
    out.contemporaries = sorted({k[0] for k in cells if k[0] not in ancestors})
    return out


def build_records(notes: Notes, ids: dict[str, int], inp: BlupV2Input) -> list[RecordIn]:
    """Un RecordIn par cellule (valeurs moyennées). ids : code -> id négatif (ancêtres et contemporaines).
    Rucher = évaluateur = un entier négatif par (source, testeur ou éleveur), sous ceux de l'exploitation et
    sous -1 (évaluateur inconnu dans records.py)."""
    used = [r.apiary for r in inp.records] + [r.evaluator for r in inp.records if r.evaluator is not None]
    start = min(-1, min(used, default=-1))
    apiaries = {a: start - k for k, a in enumerate(sorted({(k[3], k[4]) for k in notes.cells}), start=1)}
    records = []
    for (code, cid, date, source, apiary), values in sorted(notes.cells.items()):
        x = apiaries[(source, apiary)]
        mean = min(1.0, max(0.0, sum(values) / len(values)))
        records.append(RecordIn(queen=ids[code], criterion=cid, value=mean, date=date, apiary=x, evaluator=x))
    return records


def records_message(records: list[RecordIn], ancestor_ids: set[int], notes: Notes) -> Optional[str]:
    """Message de external_records_added ; None si rien n'a été ajouté ni écarté."""
    if not records and not notes.discarded:
        return None
    queens = {r.queen for r in records}
    msg = (f"{len(records)} note(s) Index Mellifera ajoutée(s) sur {len({r.criterion for r in records})} critère(s) : "
           f"{len(queens & ancestor_ids)} ancêtre(s) et {len(queens - ancestor_ids)} contemporaine(s)")
    parts = [f"{notes.discarded[k]} {label}" for k, label in DISCARD_LABELS.items() if notes.discarded[k]]
    if parts:
        msg += f" ; écartées : {', '.join(parts)}"
    if notes.capped:
        msg += f" ; plafond de {MAX_CONTEMPORARIES} contemporaines atteint"
    return msg
```

`blup2/external.py` : remplacer :

```python
from typing import Optional

from .pedigree import find_cycle, listed_ids
from .schemas import BORN_MIN, BlupV2Input, MatingIn, QueenIn, WarningOut

logger = logging.getLogger(__name__)
```

par :

```python
from typing import Optional

from .external_records import build_records, read_notes, records_message
from .pedigree import find_cycle, listed_ids
from .schemas import BORN_MIN, BlupV2Input, MatingIn, QueenIn, RecordIn, WarningOut

logger = logging.getLogger(__name__)
```

`blup2/external.py` : remplacer :

```python
    queens: list[QueenIn]                                        # reines de l'entrée, liens complétés
    ancestors: list[QueenIn] = field(default_factory=list)      # ancêtres Index Mellifera (id négatifs)
    warnings: list[WarningOut] = field(default_factory=list)

```

par :

```python
    queens: list[QueenIn]                                        # reines de l'entrée, liens complétés
    ancestors: list[QueenIn] = field(default_factory=list)      # ancêtres Index Mellifera (id négatifs)
    contemporaries: list[QueenIn] = field(default_factory=list)  # reines notées dans les mêmes groupes (id négatifs)
    records: list[RecordIn] = field(default_factory=list)       # notes des ancêtres et des contemporaines
    warnings: list[WarningOut] = field(default_factory=list)

```

`blup2/external.py` : remplacer :

```python
            msg += f", remontée arrêtée au plafond de {MAX_ANCESTORS} ancêtres"
        out.warnings.append(WarningOut(code="external_ref_resolved", message=msg))
    unknown = sorted(ref for ref, code in found.items() if code is None)
    if unknown:
```

par :

```python
            msg += f", remontée arrêtée au plafond de {MAX_ANCESTORS} ancêtres"
        out.warnings.append(WarningOut(code="external_ref_resolved", message=msg))
    try:                # notes en échec : calcul avec l'ascendance, sans les notes (spec volet 3 §6)
        _add_notes(out, inp, base, rows, ids)
    except Exception:
        logger.exception("notes Index Mellifera illisibles")
        out.contemporaries, out.records = [], []
        out.warnings.append(WarningOut(code="external_records_unavailable", message=(
            "notes Index Mellifera illisibles : calcul avec l'ascendance mais sans leurs notes")))
    unknown = sorted(ref for ref, code in found.items() if code is None)
    if unknown:
```

`blup2/external.py` : remplacer :

```python
            f"{listed_ids(unknown)}")))
    return out

```

par :

```python
            f"{listed_ids(unknown)}")))
    return out


def _add_notes(out: Expansion, inp: BlupV2Input, base: _Base, rows: dict[str, sqlite3.Row], ids: dict[str, int]) -> None:
    """Notes des ancêtres et de leurs contemporaines (spec volet 3 §4). Les contemporaines prennent les id
    négatifs suivants et sont construites comme les ancêtres : liées à leur mère et à leur côté mâle quand ils
    sont dans l'ensemble."""
    if not ids:
        return
    notes = read_notes(base.conn, inp, {code: r["birth_year"] for code, r in rows.items()})
    start = min(ids.values())
    all_ids = dict(ids)
    for k, code in enumerate(notes.contemporaries, start=1):
        all_ids[code] = start - k
    queens = list(out.ancestors) + [_ancestor(base.row(code), all_ids)[0] for code in notes.contemporaries]
    _cut_cycles(queens)               # un cycle ne peut passer que par des contemporaines
    out.ancestors, out.contemporaries = queens[: len(out.ancestors)], queens[len(out.ancestors):]
    out.records = build_records(notes, all_ids, inp)
    msg = records_message(out.records, set(ids.values()), notes)
    if msg is not None:
        out.warnings.append(WarningOut(code="external_records_added", message=msg))

```

`blup2/engine.py` : remplacer :

```python
def _evaluate(inp: BlupV2Input):
    ext = expand_external(inp)      # ancêtres Index Mellifera : dans le pedigree seulement (spec IM §8.3)
    try:
        ped, warnings = build_pedigree(ext.queens, inp.options.default_n_drones, extra_queens=ext.ancestors)
        warnings.extend(ext.warnings)
    except ValueError:
        if not ext.ancestors:
            raise
        logger.exception("pedigree invalide avec les ancêtres Index Mellifera : calcul sans eux")
        ped, warnings = build_pedigree(inp.queens, inp.options.default_n_drones)
        warnings.extend(w for w in ext.warnings if w.code != "external_ref_resolved")
        warnings.append(WarningOut(code="external_pedigree_unavailable", message=(
            "ascendance Index Mellifera incompatible avec le pedigree de l'exploitation : calcul sans les ancêtres")))
    Ainv = inverse_relationship(ped)
    F = inbreeding(ped)
    f_workers = {q: float(F[ped.w_node(q)]) for q in ped.queen_ids}
    fit_fn = _fitter(inp.profile, warnings)
    params = FitParams(estimate_variances=inp.options.estimate_variances)
```

par :

```python
def _evaluate(inp: BlupV2Input):
    ext = expand_external(inp)      # ancêtres Index Mellifera : dans le pedigree seulement (spec IM §8.3)
    extra_queens = ext.ancestors + ext.contemporaries             # et leurs notes (spec volet 3 §5)
    extra_records = ext.records
    try:
        ped, warnings = build_pedigree(ext.queens, inp.options.default_n_drones, extra_queens=extra_queens)
        warnings.extend(ext.warnings)
    except ValueError:
        if not extra_queens:
            raise
        logger.exception("pedigree invalide avec les ancêtres Index Mellifera : calcul sans eux")
        ped, warnings = build_pedigree(inp.queens, inp.options.default_n_drones)
        extra_queens, extra_records = [], []
        warnings.extend(w for w in ext.warnings if w.code not in ("external_ref_resolved", "external_records_added"))
        warnings.append(WarningOut(code="external_pedigree_unavailable", message=(
            "ascendance Index Mellifera incompatible avec le pedigree de l'exploitation : calcul sans les ancêtres")))
    Ainv = inverse_relationship(ped)
    F = inbreeding(ped)
    f_workers = {q: float(F[ped.w_node(q)]) for q in (*ped.queen_ids, *ped.added_queen_ids)}
    fit_fn = _fitter(inp.profile, warnings)
    params = FitParams(estimate_variances=inp.options.estimate_variances)
```

`blup2/engine.py` : remplacer :

```python
    for c in inp.criteria:
        try:
            cdata = prepare_criterion(inp, c.id, f_workers=f_workers)
            fit = fit_fn(cdata, ped, Ainv, params)
            fits[c.id] = fit
```

par :

```python
    for c in inp.criteria:
        try:
            cdata = prepare_criterion(inp, c.id, f_workers=f_workers, extra_queens=extra_queens,
                                      extra_records=extra_records)
            fit = fit_fn(cdata, ped, Ainv, params)
            fits[c.id] = fit
```

`blup2/model_standard.py` : remplacer :

```python
    Zw = incidence(w_obs, n_nodes)
    Zpe = incidence([rec_pos[q] for q in cdata.queen_ids], len(rec_pos))
    mean_aww = float(np.mean([diag[ped.w_node(q)] for q in cdata.recorded_queens]))
    raw_p = phenotypic_variance(cdata.X, cdata.y, floor=0.0)
    no_variation = raw_p < NO_VARIATION_SIGMA2_P
    if no_variation:   # pas de REML : il partirait vers des variances de 1e-300 (index aberrants)
```

par :

```python
    Zw = incidence(w_obs, n_nodes)
    Zpe = incidence([rec_pos[q] for q in cdata.queen_ids], len(rec_pos))
    farm_queens = cdata.farm_queens       # a priori et compteurs : exploitation seule (spec volet 3 §5)
    mean_aww = float(np.mean([diag[ped.w_node(q)] for q in farm_queens]))
    raw_p = phenotypic_variance(cdata.prior_X, cdata.prior_y, floor=0.0)
    no_variation = raw_p < NO_VARIATION_SIGMA2_P
    if no_variation:   # pas de REML : il partirait vers des variances de 1e-300 (index aberrants)
```

`blup2/model_standard.py` : remplacer :

```python
    source, h2, rep = _variance_source(cdata)
    comp = prior_components(sigma2_p, h2, rep, mean_aww)
    n_q = len(cdata.recorded_queens)
    if source == "prior" and not no_variation and should_estimate(
            params.estimate_variances, n_q, cdata.n_raw_records, cdata.n_contemporary_groups):
        I_pe = sp.identity(n_q, format="csr")
        res = reml(cdata.X, [Zw, Zpe], [Ainv, I_pe], cdata.y, [comp.sigma2_A, comp.sigma2_pe], comp.sigma2_e)
```

par :

```python
    source, h2, rep = _variance_source(cdata)
    comp = prior_components(sigma2_p, h2, rep, mean_aww)
    n_q = len(cdata.recorded_queens)      # dimension de Z_pe (ancêtres notés compris)
    n_farm = len(farm_queens)
    if source == "prior" and not no_variation and should_estimate(
            params.estimate_variances, n_farm, cdata.n_raw_records, cdata.n_contemporary_groups):
        I_pe = sp.identity(n_q, format="csr")
        res = reml(cdata.X, [Zw, Zpe], [Ainv, I_pe], cdata.y, [comp.sigma2_A, comp.sigma2_pe], comp.sigma2_e)
```

`blup2/model_standard.py` : remplacer :

```python
        else:
            status[q] = "none"
    if n_q < FEW_RECORDS_QUEENS:
        warnings.append(WarningOut(code="few_records", criterion=cid,
                                   message=f"seulement {n_q} reine(s) notée(s) : index peu fiables"))
    if cdata.scale_mixed:
        warnings.append(WarningOut(code="scale_mixed", criterion=cid,
```

par :

```python
        else:
            status[q] = "none"
    if n_farm < FEW_RECORDS_QUEENS:
        warnings.append(WarningOut(code="few_records", criterion=cid,
                                   message=f"seulement {n_farm} reine(s) notée(s) : index peu fiables"))
    if cdata.scale_mixed:
        warnings.append(WarningOut(code="scale_mixed", criterion=cid,
```

`blup2/model_standard.py` : remplacer :

```python
        source=source, h2=comp.h2, repeatability=comp.repeatability, sigma2_a=comp.sigma2_a,
        sigma2_pe=comp.sigma2_pe, sigma2_e=comp.sigma2_e, profile_used="standard",
        n_records=cdata.n_raw_records, n_queens=n_q, n_contemporary_groups=cdata.n_contemporary_groups,
        reliability_method=method,
    )
```

par :

```python
        source=source, h2=comp.h2, repeatability=comp.repeatability, sigma2_a=comp.sigma2_a,
        sigma2_pe=comp.sigma2_pe, sigma2_e=comp.sigma2_e, profile_used="standard",
        n_records=cdata.n_raw_records, n_queens=n_farm, n_contemporary_groups=cdata.n_contemporary_groups,
        reliability_method=method,
    )
```

`blup2/model_breeder.py` : remplacer :

```python

def known_mating_share(cdata: CriterionData, ped: BeePedigree) -> float:
    """Part des reines notées dont l'accouplement (nœud sire) est connu."""
    if not cdata.recorded_queens:
        return 0.0
    known = sum(1 for q in cdata.recorded_queens if ped.sire_of_mating(q) is not None)
    return known / len(cdata.recorded_queens)


```

par :

```python

def known_mating_share(cdata: CriterionData, ped: BeePedigree) -> float:
    """Part des reines notées de l'exploitation dont l'accouplement (nœud sire) est connu."""
    queens = cdata.farm_queens
    if not queens:
        return 0.0
    known = sum(1 for q in queens if ped.sire_of_mating(q) is not None)
    return known / len(queens)


```

`blup2/model_breeder.py` : remplacer :

```python
    q_obs = [ped.q_node(q) for q in cdata.queen_ids]
    rec_pos = {q: k for k, q in enumerate(cdata.recorded_queens)}
    n_q = len(rec_pos)
    Zw, Zq = incidence(w_obs, n), incidence(q_obs, n)
    Zpe = incidence([rec_pos[q] for q in cdata.queen_ids], n_q)
    rec_w = [ped.w_node(q) for q in cdata.recorded_queens]
    rec_q = [ped.q_node(q) for q in cdata.recorded_queens]
    mean_aww = float(np.mean(diag[rec_w]))
    mean_aqq = float(np.mean(diag[rec_q]))
    mean_awq = float(np.mean([relationship(ped, w, q) for w, q in zip(rec_w, rec_q)]))
    raw_p = phenotypic_variance(cdata.X, cdata.y, floor=0.0)
    no_variation = raw_p < NO_VARIATION_SIGMA2_P
    if no_variation:   # pas de REML : il partirait vers des variances de 1e-300 (index aberrants)
```

par :

```python
    q_obs = [ped.q_node(q) for q in cdata.queen_ids]
    rec_pos = {q: k for k, q in enumerate(cdata.recorded_queens)}
    n_q = len(rec_pos)                    # dimension de Z_pe (ancêtres notés compris)
    farm_queens = cdata.farm_queens       # a priori et compteurs : exploitation seule (spec volet 3 §5)
    n_farm = len(farm_queens)
    Zw, Zq = incidence(w_obs, n), incidence(q_obs, n)
    Zpe = incidence([rec_pos[q] for q in cdata.queen_ids], n_q)
    rec_w = [ped.w_node(q) for q in farm_queens]
    rec_q = [ped.q_node(q) for q in farm_queens]
    mean_aww = float(np.mean(diag[rec_w]))
    mean_aqq = float(np.mean(diag[rec_q]))
    mean_awq = float(np.mean([relationship(ped, w, q) for w, q in zip(rec_w, rec_q)]))
    raw_p = phenotypic_variance(cdata.prior_X, cdata.prior_y, floor=0.0)
    no_variation = raw_p < NO_VARIATION_SIGMA2_P
    if no_variation:   # pas de REML : il partirait vers des variances de 1e-300 (index aberrants)
```

`blup2/model_breeder.py` : remplacer :

```python
    mode = params.estimate_variances
    if source == "prior" and not no_variation and should_estimate(
            mode, n_q, cdata.n_raw_records, cdata.n_contemporary_groups):
        res = reml_breeder(cdata.X, Zw, Zq, Zpe, Ainv, cdata.y, comp, estimate_partition=(mode == "always"))
        if res.converged:
```

par :

```python
    mode = params.estimate_variances
    if source == "prior" and not no_variation and should_estimate(
            mode, n_farm, cdata.n_raw_records, cdata.n_contemporary_groups):
        res = reml_breeder(cdata.X, Zw, Zq, Zpe, Ainv, cdata.y, comp, estimate_partition=(mode == "always"))
        if res.converged:
```

`blup2/model_breeder.py` : remplacer :

```python
        else:
            status[q] = "none"
    if n_q < FEW_RECORDS_QUEENS:
        warnings.append(WarningOut(code="few_records", criterion=cid,
                                   message=f"seulement {n_q} reine(s) notée(s) : index peu fiables"))
    if cdata.scale_mixed:
        warnings.append(WarningOut(code="scale_mixed", criterion=cid,
```

par :

```python
        else:
            status[q] = "none"
    if n_farm < FEW_RECORDS_QUEENS:
        warnings.append(WarningOut(code="few_records", criterion=cid,
                                   message=f"seulement {n_farm} reine(s) notée(s) : index peu fiables"))
    if cdata.scale_mixed:
        warnings.append(WarningOut(code="scale_mixed", criterion=cid,
```

`blup2/model_breeder.py` : remplacer :

```python
        sigma2_pe=comp.sigma2_pe, sigma2_e=comp.sigma2_e, sigma2_aq=comp.sigma2_AQ, rg=comp.rg,
        profile_used="breeder", known_mating_share=_published_share(share), n_records=cdata.n_raw_records,
        n_queens=n_q, n_contemporary_groups=cdata.n_contemporary_groups, reliability_method=method,
    )
    return CriterionFit(
```

par :

```python
        sigma2_pe=comp.sigma2_pe, sigma2_e=comp.sigma2_e, sigma2_aq=comp.sigma2_AQ, rg=comp.rg,
        profile_used="breeder", known_mating_share=_published_share(share), n_records=cdata.n_raw_records,
        n_queens=n_farm, n_contemporary_groups=cdata.n_contemporary_groups, reliability_method=method,
    )
    return CriterionFit(
```

- [ ] **Step 4 : Vérifier que les tests passent**

Run : `MELLIFERA_CONTACT= .venv/bin/pytest -q -p no:cacheprovider tests/test_external_records.py tests/test_external_engine.py tests/mellifera/test_mellifera_external_contract.py | tail -1`

Expected :

```text
41 passed in …
```

- [ ] **Step 5 : Lancer la suite complète**

```bash
MELLIFERA_CONTACT= .venv/bin/pytest -q -p no:cacheprovider | tail -1
```

Expected :

```text
459 passed, 12 deselected in …
```

Puis ruff :

```bash
uvx ruff check --select F,E9 blup2 tests
```

Expected :

```text
All checks passed!
```

- [ ] **Step 6 : Commit**

```bash
git add blup2/external_records.py blup2/external.py blup2/engine.py blup2/model_standard.py blup2/model_breeder.py tests/external_base.py tests/test_external_records.py tests/mellifera/test_mellifera_external_contract.py tests/test_external_engine.py
git status --short
git commit -q -F - <<'EOF'
Lit les notes Index Mellifera des ancêtres et de leurs contemporaines dans le moteur

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015J8yR1HR7mLSRtwnvEYap6
EOF
git log --oneline -1
```

Expected : `git status --short` ne liste que les 9 fichiers ci-dessus, puis :

```text
<hash> Lit les notes Index Mellifera des ancêtres et de leurs contemporaines dans le moteur
```

---

### Task 4 : Fiabilité exacte malgré les notes externes

Règle de spec §5 : la méthode de fiabilité se décide sans les équations dues aux seules notes externes. Quand elles font dépasser `EXACT_LIMIT`, les PEV des reines de l'exploitation sont lues par résolutions du système déjà factorisé (`inverse_entries`), sans inverser la matrice entière.

**Files:**
- Modify: `blup2/mme.py`
- Modify: `blup2/engine.py`
- Modify: `blup2/model_standard.py`
- Modify: `blup2/model_breeder.py`
- Test (modify): `tests/test_mme.py`
- Test (modify): `tests/test_external_engine.py`

**Interfaces:**
- Consumes :
  - tâche 2 : `CriterionData.prior_X`, `recorded_queens`, `farm_queens` ;
  - tâche 3 : `Expansion.contemporaries`, `Expansion.ancestors`, câblage `extra_queens` dans `engine._evaluate` ;
  - existant : `MMESolution.solve(B)`, `.n_equations`, `.n_fixed`, `Pedigree.w_node(q)`, `FitParams.exact_limit`.
- Produces :
  - `blup2.mme.inverse_entries(sol, rows, cols, sigma2_e: float, chunk: int = 256) -> np.ndarray` ;
  - `FitParams.external_nodes: int = 0` ;
  - `blup2.model_standard.external_equations(cdata, params, n_genetic: int) -> int`, importée par `model_breeder`.

- [ ] **Step 1 : Écrire les tests de `inverse_entries` et de la méthode de fiabilité**

`tests/test_mme.py` : remplacer :

```python
        solve_mme(X, [sp.csr_matrix(Z)], sp.identity(3), 1.0, y, False)

```

par :

```python
        solve_mme(X, [sp.csr_matrix(Z)], sp.identity(3), 1.0, y, False)


def test_inverse_entries_match_the_dense_inverse_above_the_limit():
    """Spec volet 3 §5 : PEV exactes de quelques nœuds par résolutions du système factorisé."""
    from blup2.mme import inverse_entries
    X, Z, G, y = _problem(seed=2)
    Ginv = sp.csr_matrix(np.linalg.inv(G))
    dense = solve_mme(X, [sp.csr_matrix(Z)], Ginv, 0.7, y, True)
    big = solve_mme(X, [sp.csr_matrix(Z)], Ginv, 0.7, y, True, exact_limit=10)
    rows, cols = [3, 5, 14, 3], [3, 5, 14, 14]
    got = inverse_entries(big, rows, cols, 0.7, chunk=2)
    assert np.allclose(got, dense.cinv[rows, cols], atol=1e-12)

```

`tests/test_external_engine.py` : remplacer :

```python



def test_reml_with_external_notes_counts_only_the_farm(tmp_path, monkeypatch):
```

par :

```python


@pytest.mark.parametrize("profile", ["standard", "breeder"])
def test_notes_never_switch_the_reliability_method(tmp_path, monkeypatch, profile):
    """Spec volet 3 §5 : la méthode de fiabilité se décide sur les équations de l'exploitation et de son
    ascendance ; au-delà de la limite à cause des seules notes externes, les PEV des reines de l'exploitation
    restent exactes, sans inverser la matrice entière."""
    import functools

    import blup2.engine as engine
    from blup2.model_standard import FitParams
    db = make_db(tmp_path / "m.sqlite", notes=NOTES)

    monkeypatch.setenv("MELLIFERA_DB", db)
    qs = [queen(k, 2023, mother_ref_ext="FR-CBP-B-0276-2022" if k <= 4 else None,
                mating=mating("unknown", ext="FR-CBP-B-0100-2020")) for k in range(1, 9)]   # accouplements connus
    recs = [record(k, 0.1 * k, date="2024-06-01", apiary=1 + k % 2) for k in range(1, 9)]

    def solve(limit, criteria):
        monkeypatch.setattr(engine, "FitParams", functools.partial(FitParams, exact_limit=limit))
        out = run_blup(inp(qs, records=recs, criteria=criteria, profile=profile))
        return out, {o.id: o.indexes["1"].reliability for o in out.queens}

    sizes = [int(warning(solve(0, c)[0], "reliability_approx")[0].split()[0])
             for c in ([criterion(1, h2=0.3, repeatability=0.4)], [criterion(1, **GENTLE)])]
    assert sizes[0] < sizes[1]                                    # les notes ajoutent des équations
    out, rel = solve(sizes[0], [criterion(1, **GENTLE)])
    assert out.parameters["1"].reliability_method == "exact" and warning(out, "reliability_approx") == []
    assert out.parameters["1"].profile_used == profile
    _, dense = solve(10 ** 6, [criterion(1, **GENTLE)])
    assert rel == pytest.approx(dense, abs=1e-9)
    out, _ = solve(sizes[0] - 1, [criterion(1, **GENTLE)])      # l'exploitation seule dépasse : approché
    assert out.parameters["1"].reliability_method == "approx"


def test_reml_with_external_notes_counts_only_the_farm(tmp_path, monkeypatch):
```

- [ ] **Step 2 : Vérifier l'échec**

Run : `MELLIFERA_CONTACT= .venv/bin/pytest -q -p no:cacheprovider -rfE --tb=no --continue-on-collection-errors tests/test_mme.py tests/test_external_engine.py | tail -8`

`inverse_entries` n'existe pas ; avec la limite fixée au nombre d'équations du calcul sans notes, le calcul avec notes passe en fiabilité approchée, dans les deux profils.

Expected :

```text
....F...................FF.                                              [100%]
=========================== short test summary info ============================
FAILED tests/test_mme.py::test_inverse_entries_match_the_dense_inverse_above_the_limit
FAILED tests/test_external_engine.py::test_notes_never_switch_the_reliability_method[standard]
FAILED tests/test_external_engine.py::test_notes_never_switch_the_reliability_method[breeder]
3 failed, 24 passed in …
```

- [ ] **Step 3 : Écrire `inverse_entries` et la branche exacte des deux modèles**

`blup2/mme.py` : remplacer :

```python
    )

```

par :

```python
    )


def inverse_entries(sol: MMESolution, rows, cols, sigma2_e: float, chunk: int = 256) -> np.ndarray:
    """σ²_e · C⁻¹[rows[k], cols[k]] par résolutions du système déjà factorisé, une par indice distinct de
    cols, par paquets de `chunk` : PEV exactes de quelques nœuds sans inverser C (spec volet 3 §5)."""
    rows, cols = np.asarray(rows, dtype=int), np.asarray(cols, dtype=int)
    uniq, pos = np.unique(cols, return_inverse=True)
    out = np.empty(len(rows))
    for start in range(0, len(uniq), chunk):
        block = uniq[start:start + chunk]
        E = np.zeros((sol.n_equations, len(block)))
        E[block, np.arange(len(block))] = 1.0
        C = sol.solve(E)
        sel = (pos >= start) & (pos < start + len(block))
        out[sel] = C[rows[sel], pos[sel] - start]
    return sigma2_e * out

```

`blup2/engine.py` : remplacer :

```python
    f_workers = {q: float(F[ped.w_node(q)]) for q in (*ped.queen_ids, *ped.added_queen_ids)}
    fit_fn = _fitter(inp.profile, warnings)
    params = FitParams(estimate_variances=inp.options.estimate_variances)
    fits: dict[int, Optional[CriterionFit]] = {}
    for c in inp.criteria:
```

par :

```python
    f_workers = {q: float(F[ped.w_node(q)]) for q in (*ped.queen_ids, *ped.added_queen_ids)}
    fit_fn = _fitter(inp.profile, warnings)
    external_nodes = 0
    if ext.contemporaries and extra_queens:    # nœuds des seules contemporaines (spec volet 3 §5)
        external_nodes = ped.n_nodes - build_pedigree(ext.queens, inp.options.default_n_drones,
                                                      extra_queens=ext.ancestors)[0].n_nodes
    params = FitParams(estimate_variances=inp.options.estimate_variances, external_nodes=external_nodes)
    fits: dict[int, Optional[CriterionFit]] = {}
    for c in inp.criteria:
```

`blup2/model_standard.py` : remplacer :

```python

from .kinship import relationship_diagonal
from .mme import EXACT_LIMIT, solve_mme
from .pedigree import BeePedigree
from .records import CriterionData
```

par :

```python

from .kinship import relationship_diagonal
from .mme import EXACT_LIMIT, inverse_entries, solve_mme
from .pedigree import BeePedigree
from .records import CriterionData
```

`blup2/model_standard.py` : remplacer :

```python
    estimate_variances: str = "auto"      # "auto" | "always" | "never"
    exact_limit: int = EXACT_LIMIT


```

par :

```python
    estimate_variances: str = "auto"      # "auto" | "always" | "never"
    exact_limit: int = EXACT_LIMIT
    external_nodes: int = 0               # nœuds du pedigree dus aux seules contemporaines (spec volet 3 §5)


def external_equations(cdata: CriterionData, params: FitParams, n_genetic: int) -> int:
    """Équations dues aux seules notes externes : effets fixes propres, environnement permanent des reines
    externes notées, nœuds des contemporaines (n_genetic blocs génétiques). La méthode de fiabilité se décide
    sur les autres (spec volet 3 §5)."""
    return (cdata.X.shape[1] - cdata.prior_X.shape[1] + len(cdata.recorded_queens) - len(cdata.farm_queens)
            + n_genetic * params.external_nodes)


```

`blup2/model_standard.py` : remplacer :

```python
    if sol.pev is not None:
        pev, method = sol.pev[0], "exact"
    else:
        n_eff = effective_records({"mean": [0] * cdata.n_obs, **cdata.factor_labels}, w_obs, n_nodes)
```

par :

```python
    if sol.pev is not None:
        pev, method = sol.pev[0], "exact"
    elif sol.n_equations - external_equations(cdata, params, 1) <= params.exact_limit:
        nodes = sol.n_fixed + np.array([ped.w_node(q) for q in ped.queen_ids], dtype=int)
        pev, method = np.full(n_nodes, np.nan), "exact"     # seules les reines de l'exploitation sont lues
        pev[nodes - sol.n_fixed] = inverse_entries(sol, nodes, nodes, comp.sigma2_e)
    else:
        n_eff = effective_records({"mean": [0] * cdata.n_obs, **cdata.factor_labels}, w_obs, n_nodes)
```

`blup2/model_breeder.py` : remplacer :

```python

from .kinship import relationship, relationship_diagonal
from .mme import solve_mme
from .model_standard import (
    FEW_RECORDS_QUEENS,
```

par :

```python

from .kinship import relationship, relationship_diagonal
from .mme import inverse_entries, solve_mme
from .model_standard import (
    FEW_RECORDS_QUEENS,
```

`blup2/model_breeder.py` : remplacer :

```python
    FitParams,
    empty_fit,
    fit_standard,
    incidence,
```

par :

```python
    FitParams,
    empty_fit,
    external_equations,
    fit_standard,
    incidence,
```

`blup2/model_breeder.py` : remplacer :

```python
        pev = cinv[p + idx, p + idx] + cinv[p + n + idx, p + n + idx] + 2.0 * cinv[p + idx, p + n + idx]
        method = "exact"
    else:
        pev = _approx_pev_total(cdata, ped, comp, w_obs, diag)
```

par :

```python
        pev = cinv[p + idx, p + idx] + cinv[p + n + idx, p + n + idx] + 2.0 * cinv[p + idx, p + n + idx]
        method = "exact"
    elif sol.n_equations - external_equations(cdata, params, 2) <= params.exact_limit:
        w = np.array([ped.w_node(q) for q in ped.queen_ids], dtype=int)
        a, b, k = sol.n_fixed + w, sol.n_fixed + n + w, len(w)
        v = inverse_entries(sol, np.r_[a, b, a], np.r_[a, b, b], comp.sigma2_e)
        pev, method = np.full(n, np.nan), "exact"          # seules les reines de l'exploitation sont lues
        pev[w] = v[:k] + v[k:2 * k] + 2.0 * v[2 * k:]
    else:
        pev = _approx_pev_total(cdata, ped, comp, w_obs, diag)
```

- [ ] **Step 4 : Vérifier que les tests passent**

Run : `MELLIFERA_CONTACT= .venv/bin/pytest -q -p no:cacheprovider tests/test_mme.py tests/test_external_engine.py | tail -1`

Expected :

```text
27 passed in …
```

- [ ] **Step 5 : Lancer la suite complète**

```bash
MELLIFERA_CONTACT= .venv/bin/pytest -q -p no:cacheprovider | tail -1
```

Expected :

```text
462 passed, 12 deselected in …
```

Puis ruff :

```bash
uvx ruff check --select F,E9 blup2 tests
```

Expected :

```text
All checks passed!
```

- [ ] **Step 6 : Commit**

```bash
git add blup2/mme.py blup2/engine.py blup2/model_standard.py blup2/model_breeder.py tests/test_mme.py tests/test_external_engine.py
git status --short
git commit -q -F - <<'EOF'
Garde des fiabilités exactes quand seules les notes externes dépassent la limite

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015J8yR1HR7mLSRtwnvEYap6
EOF
git log --oneline -1
```

Expected : `git status --short` ne liste que les 6 fichiers ci-dessus, puis :

```text
<hash> Garde des fiabilités exactes quand seules les notes externes dépassent la limite
```

---

### Task 5 : Validation statistique et documentation

La porte statistique de la spec §7 et la mesure de coût du §9, en tests lents, puis la documentation du volet. Ces tests valident le comportement livré par les tâches 3 et 4 : ils passent dès leur première exécution. Ce n'est pas un écart au TDD, c'est une validation ; les seuils viennent du prototype (20 populations).

**Files:**
- Modify: `README.fr.md`
- Modify: `README.md`
- Modify: `METHODOLOGIE.md`
- Test (modify): `tests/simdata.py`
- Test (create): `tests/test_external_calibration.py`

**Interfaces:**
- Consumes : `run_blup`, `BlupV2Input`, `tests.external_base.SCHEMA` et `NOTES_SCHEMA` (tâche 3), `simulate.gen_pedigree` (`docs/analyses`, déjà importé par `tests/simdata.py`), avertissements `external_records_added` et `reliability_approx`.
- Produces :
  - `tests.simdata.external_scenario(seed, s2e=1.5, s2pe=0.5, s2cg=0.5, S=8, D=12, p_ext=0.8, n_testers=6, p_farm=0.7, n_farm_ap=4, n_free=30, top_half=False)` ;
  - `tests.simdata.write_crawler_db(path, rows, notes) -> str` ;
  - `METHODOLOGIE.md` §3 bis « Notes des lignées Index Mellifera ».

- [ ] **Step 1 : Écrire le scénario simulé et les tests lents**

`tests/simdata.py` : remplacer :

```python
    return inp, truth, recorded, true

```

par :

```python
    return inp, truth, recorded, true


def external_scenario(seed: int, s2e: float = 1.5, s2pe: float = 0.5, s2cg: float = 0.5, S: int = 8, D: int = 12,
                      p_ext: float = 0.8, n_testers: int = 6, p_farm: float = 0.7, n_farm_ap: int = 4,
                      n_free: int = 30, top_half: bool = False):
    """Lignées Index Mellifera simulées (spec volet 3 §7) : pedigree de simulate.gen_pedigree sur 5
    générations ; les générations 0 à 3 forment la base du crawler (notes de saison sur 0 à 6, un groupe par
    testeur × année), la génération 4 est l'exploitation (filles de reines externes, accouplées sur des
    sœurs filles d'une reine externe). n_free reines sans lien, dans leurs propres ruchers, complètent
    l'exploitation. top_half : la base ne garde que la meilleure moitié des reines externes notées (sur la note).
    Vraies valeurs tirées avec la matrice A de référence (σ²_A = 1), variances transmises au critère.
    Retourne (entrée v2 en dict, lignes de la table queen, lignes de notes au format de
    tests/external_base.NOTES, vérité {reine liée: valeur du nœud W normalisée}, reines sans lien)."""
    from simulate import build_bb, gen_pedigree

    rng = np.random.default_rng(seed)
    inds, mate = gen_pedigree(rng, G=4)
    P = build_bb(inds, S, D)
    A = P.build()
    L = np.linalg.cholesky(A + 1e-12 * np.eye(len(inds)))
    a_true = L @ rng.standard_normal(len(inds))
    idx = {x["id"]: i for i, x in enumerate(inds)}
    sires = {x["id"]: x for x in inds if x["kind"] == "S"}
    queens = [x for x in inds if x["kind"] == "Q"]
    external = [x for x in queens if x["gen"] < 4]
    code = {x["id"]: f"FR-SIM-B-{k + 1:04d}-{2015 + x['gen']}" for k, x in enumerate(external)}
    width = SCALE_MAX - SCALE_MIN

    rows = [(code[x["id"]], 2015 + x["gen"], code.get(x["dam"]), code.get(sires[mate[x["id"]]]["dam"]), S)
            for x in external]
    cg = {}
    notes = []
    for x in external:
        if rng.random() >= p_ext:
            continue
        tester, year = int(rng.integers(n_testers)), 2016 + x["gen"]
        cg.setdefault((tester, year), rng.standard_normal() * np.sqrt(s2cg))
        y = (5.0 + a_true[idx["W_" + x["id"]]] + cg[(tester, year)] + rng.standard_normal() * np.sqrt(s2pe)
             + rng.standard_normal() * np.sqrt(s2e))
        value = float(np.clip((y - SCALE_MIN) / width * 6, 0, 6))
        notes.append((code[x["id"]], "season_table", "evjp", "mildness", value, f"{year}-07-01", 300 + tester))
    if top_half:
        cut = float(np.median([n[4] for n in notes]))
        notes = [n for n in notes if n[4] >= cut]

    farm = [x for x in queens if x["gen"] == 4]
    ap_eff = rng.standard_normal(n_farm_ap + 2) * np.sqrt(s2cg)
    q_in, records, truth, recorded = [], [], {}, []
    for k, x in enumerate(farm):
        q = 1000 + k
        truth[q] = a_true[idx["W_" + x["id"]]] / width
        q_in.append({"id": q, "born": 2019, "mother_ref_ext": code[x["dam"]],
                     "mating": {"type": "unknown", "external_ref": code[sires[mate[x["id"]]]["dam"]],
                                "n_drone_queens": S, "n_drones": D}})
        if rng.random() >= p_farm:
            continue
        recorded.append(x["id"])
        ap, pe = int(rng.integers(n_farm_ap)), rng.standard_normal() * np.sqrt(s2pe)
        for j in range(int(rng.integers(2, 6))):
            y = 5.0 + a_true[idx["W_" + x["id"]]] + pe + ap_eff[ap] + rng.standard_normal() * np.sqrt(s2e)
            records.append({"queen": q, "criterion": 1, "value": float(y), "scale_min": SCALE_MIN,
                            "scale_max": SCALE_MAX, "date": f"2020-0{4 + j}-10", "evaluator": 1, "apiary": ap + 1})
    free = list(range(5000, 5000 + n_free))
    for q in free:
        ap, g = n_farm_ap + int(rng.integers(2)), rng.standard_normal() * np.sqrt(0.5)
        q_in.append({"id": q, "born": 2019})
        for j in range(int(rng.integers(2, 6))):
            y = 5.0 + g + ap_eff[ap] + rng.standard_normal() * np.sqrt(s2pe + s2e)
            records.append({"queen": q, "criterion": 1, "value": float(y), "scale_min": SCALE_MIN,
                            "scale_max": SCALE_MAX, "date": f"2020-0{4 + j}-10", "evaluator": 1, "apiary": ap + 1})
    aww = float(np.mean([A[idx["W_" + q], idx["W_" + q]] for q in recorded]))
    s2p = aww + s2pe + s2e
    crit = {"id": 1, "name": "sim", "direction": 1, "kind": "note", "weight": 1, "h2": aww / s2p,
            "repeatability": (aww + s2pe) / s2p}
    inp = {"exploitation": 1, "profile": "standard", "options": {"n_suggestions": 0},
           "criteria": [crit], "queens": q_in, "records": records}
    return inp, rows, notes, truth, free


def write_crawler_db(path, rows, notes) -> str:
    """Base du crawler minimale (schéma de tests/external_base.py) pour external_scenario."""
    import sqlite3

    from tests.external_base import NOTES_SCHEMA, SCHEMA

    conn = sqlite3.connect(path)
    conn.executescript(SCHEMA + NOTES_SCHEMA)
    conn.executemany("INSERT INTO queen (code, birth_year, mother_code, drone_mother_code, n_drone_queens, "
                     "mating_code, breeder_id) VALUES (?, ?, ?, ?, ?, 'ins', 900)", rows)
    conn.executemany(
        "INSERT INTO evaluation (queen_name_raw, queen_code, queen_key, tester, tester_breeder_id, eval_date, mode, "
        "criterion_code, value, value_raw, source, url) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, '/x/')",
        [(c, c, c, f"T{t}", t, d, mode, crit, v, str(v), src) for c, src, mode, crit, v, d, t in notes])
    conn.commit()
    conn.close()
    return str(path)

```

Créer `tests/test_external_calibration.py` :

```python
"""Validation statistique du volet 3 Index Mellifera (spec volet 3 §7) et coût. Tests lents :
    .venv/bin/python -m pytest -m slow
"""
import time

import numpy as np
import pytest
from scipy.stats import spearmanr

from blup2.engine import run_blup
from blup2.schemas import BlupV2Input
from tests.simdata import external_scenario, write_crawler_db

N_REPLICATES = 20


def _indexes(monkeypatch, inp, db, trait):
    monkeypatch.setenv("MELLIFERA_DB", db)
    crit = {**inp["criteria"][0], **({"trait": "gentleness"} if trait else {})}
    out = run_blup(BlupV2Input.model_validate({**inp, "criteria": [crit]}))
    assert out.parameters["1"].reliability_method == "exact"   # le gain mesuré n'est pas un effet de l'approximation
    return {q.id: q.indexes["1"] for q in out.queens}


def _replicates(tmp_path, monkeypatch, top_half=False):
    """Par répétition : gain de fiabilité moyen des reines liées, corrélations de Pearson et de Spearman avec
    la vraie valeur (avec, sans notes externes), écart maximal des EBV des reines sans lien."""
    res = []
    for seed in range(1, N_REPLICATES + 1):
        inp, rows, notes, truth, free = external_scenario(seed, top_half=top_half)
        db = write_crawler_db(tmp_path / f"m{seed}.sqlite", rows, notes)
        a = _indexes(monkeypatch, inp, db, trait=True)
        b = _indexes(monkeypatch, inp, db, trait=False)            # ascendance du volet 2, sans les notes
        ids = list(truth)
        t = np.array([truth[q] for q in ids])
        ea, eb = (np.array([x[q].ebv for q in ids]) for x in (a, b))
        gain = float(np.mean([a[q].reliability - b[q].reliability for q in ids]))
        res.append((gain, np.corrcoef(t, ea)[0, 1], np.corrcoef(t, eb)[0, 1], spearmanr(t, ea)[0],
                    spearmanr(t, eb)[0], max(abs(a[q].ebv - b[q].ebv) for q in free)))
    return np.array(res)


@pytest.mark.slow
def test_notes_of_external_lineages_improve_their_daughters(tmp_path, monkeypatch):
    """(a) et (b) de la spec volet 3 §7, variances fixées."""
    r = _replicates(tmp_path, monkeypatch)
    # prototype (20 répétitions) : gain de fiabilité 0,068 (min 0,047) ; corrélation 0,472 -> 0,517
    assert (r[:, 0] > 0).all() and r[:, 0].mean() >= 0.05
    assert r[:, 1].mean() - r[:, 2].mean() >= 0.02
    assert r[:, 5].max() <= 1e-9                    # reines sans lien ni groupe commun : EBV inchangées


@pytest.mark.slow
def test_selected_external_notes_do_not_degrade_the_ranking(tmp_path, monkeypatch):
    """(c) de la spec volet 3 §7 : la base ne garde que la meilleure moitié des reines externes notées."""
    r = _replicates(tmp_path, monkeypatch, top_half=True)
    # prototype : gain de fiabilité 0,033 ; Spearman 0,452 sans notes, 0,472 avec
    assert (r[:, 0] > 0).all()
    assert r[:, 3].mean() >= r[:, 4].mean()


def _big_crawler_db(path, rng, n_queens=9000, years=range(2010, 2025), n_testers=40):
    """Base synthétique de taille réelle : une génération par an, mère et mère des sœurs de l'année d'avant,
    80 % des reines notées l'année suivant leur naissance par un testeur tiré au hasard."""
    per_year = n_queens // len(years)
    codes = {y: [f"FR-B{k % 50:02d}-B-{k + 1:04d}-{y}" for k in range(per_year)] for y in years}
    rows, notes = [], []
    for y in years:
        prev, prev2 = codes.get(y - 1), codes.get(y - 2)
        for c in codes[y]:
            rows.append((c, y, str(rng.choice(prev)) if prev else None, str(rng.choice(prev2)) if prev2 else None, 8))
            if rng.random() < 0.8 and y + 1 <= 2025:
                notes.append((c, "season_table", "evjp", "mildness", float(rng.integers(0, 7)), f"{y + 1}-07-01",
                              300 + int(rng.integers(n_testers))))
    return write_crawler_db(path, rows, notes), codes


@pytest.mark.slow
def test_cost_on_a_crawler_base_of_9000_queens(tmp_path, monkeypatch):
    """Spec volet 3 §9 : 300 reines de l'exploitation citant 100 lignées d'une base de 9 000 reines."""
    rng = np.random.default_rng(3)
    db, codes = _big_crawler_db(tmp_path / "big.sqlite", rng)
    mothers = rng.choice(codes[2023], 100, replace=False)
    queens, records = [], []
    for k in range(300):
        queens.append({"id": k + 1, "born": 2024, "mother_ref_ext": str(mothers[k % 100]),
                       "mating": {"type": "unknown", "external_ref": str(rng.choice(codes[2022]))}})
        ap = int(rng.integers(5)) + 1
        records += [{"queen": k + 1, "criterion": 1, "value": float(rng.integers(1, 6)), "scale_min": 1,
                     "scale_max": 5, "date": f"2025-0{4 + j}-10", "evaluator": 1, "apiary": ap} for j in range(3)]
    crit = {"id": 1, "name": "perf", "direction": 1, "kind": "note", "weight": 1, "h2": 0.3, "repeatability": 0.4}
    inp = {"exploitation": 1, "profile": "standard", "options": {"n_suggestions": 0}, "queens": queens,
           "records": records}
    monkeypatch.setenv("MELLIFERA_DB", db)
    times, outs = [], []
    for c in (crit, crit, {**crit, "trait": "gentleness"}):          # premier calcul : préchauffage
        t0 = time.perf_counter()
        outs.append(run_blup(BlupV2Input.model_validate({**inp, "criteria": [c]})))
        times.append(time.perf_counter() - t0)
    added = [w.message for w in outs[2].warnings if w.code == "external_records_added"]
    print(f"\nsans notes {times[1]:.2f} s, avec notes {times[2]:.2f} s ; {added}")
    assert added and "plafond de 1000 contemporaines atteint" in added[0]
    assert [o.parameters["1"].reliability_method for o in outs] == ["exact"] * 3    # spec volet 3 §5
    assert times[2] < 30.0
```

- [ ] **Step 2 : Lancer les tests lents du volet**

```bash
MELLIFERA_CONTACT= .venv/bin/pytest -q -p no:cacheprovider -m slow tests/test_external_calibration.py | tail -1
```

Trois tests : gain de fiabilité et de corrélation avec toutes les notes, rang non dégradé avec la meilleure moitié seulement, coût sur une base de 9 000 reines (plafond atteint, méthode exacte, moins de 30 s). Une dizaine de secondes en tout.

Expected (durées variables) :

```text
3 passed in …
```

- [ ] **Step 3 : Documenter le volet**

`README.fr.md` : remplacer :

```markdown
### Lignées Index Mellifera dans le moteur

Quand une reine cite une lignée Index Mellifera par son code (`mother_ref_ext` pour sa mère, `mating.external_ref` pour la mère des reines productrices de mâles), le moteur v2 lit la base du crawler en lecture seule et ajoute l'ascendance de cette lignée au pedigree (au plus `options.external_pedigree_depth` générations, 7 par défaut, 0 = désactivé, et au plus 1 000 ancêtres au total). Ces ancêtres servent au calcul de parenté et de consanguinité ; ils n'apparaissent dans aucune sortie et leurs évaluations ne sont jamais utilisées. Un code est cherché tel quel, puis sans espaces ni distinction de casse (la casse compte sur le site : `DE-ABg-B-0003-2017` ; deux codes qui ne diffèrent que par la casse rendent la saisie ambiguë), puis par nom affiché (`B4536(NG)`) s'il est unique ; un code introuvable est ignoré (avertissement `external_ref_unknown`).

Activation : `MELLIFERA_DB=/app/data/mellifera/mellifera.sqlite` dans `.env.production` du service `api`. Sans cette variable, la résolution est désactivée et aucun nouvel avertissement n'est émis ; si la base est absente ou illisible, une entrée qui cite des lignées reçoit l'avertissement `external_pedigree_unavailable`. La version 2.1.0 ajoute l'option `external_pedigree_depth` et trois codes d'avertissement (`external_ref_resolved`, `external_ref_unknown`, `external_pedigree_unavailable`) au contrat : **n'activez `MELLIFERA_DB` en production qu'après avoir régénéré et publié le SDK PHP** à partir de `openapi/openapi.json`.

## Structure du projet
```

par :

```markdown
### Lignées Index Mellifera dans le moteur

Quand une reine cite une lignée Index Mellifera par son code (`mother_ref_ext` pour sa mère, `mating.external_ref` pour la mère des reines productrices de mâles), le moteur v2 lit la base du crawler en lecture seule et ajoute l'ascendance de cette lignée au pedigree (au plus `options.external_pedigree_depth` générations, 7 par défaut, 0 = désactivé, et au plus 1 000 ancêtres au total). Ces ancêtres servent au calcul de parenté et de consanguinité ; ils n'apparaissent dans aucune sortie. Depuis la version 2.2.0, un critère qui déclare son caractère standard (`criteria[].trait`, liste figée dans le contrat) reçoit aussi les notes Index Mellifera de ces ancêtres et des reines notées dans les mêmes groupes testeur × année (au plus 1 000), sans effet fixe commun avec les notes de l'exploitation (voir `METHODOLOGIE.md`, §3 bis). Un code est cherché tel quel, puis sans espaces ni distinction de casse (la casse compte sur le site : `DE-ABg-B-0003-2017` ; deux codes qui ne diffèrent que par la casse rendent la saisie ambiguë), puis par nom affiché (`B4536(NG)`) s'il est unique ; un code introuvable est ignoré (avertissement `external_ref_unknown`).

Activation : `MELLIFERA_DB=/app/data/mellifera/mellifera.sqlite` dans `.env.production` du service `api`. Sans cette variable, la résolution est désactivée et aucun nouvel avertissement n'est émis ; si la base est absente ou illisible, une entrée qui cite des lignées reçoit l'avertissement `external_pedigree_unavailable`. La version 2.1.0 ajoute l'option `external_pedigree_depth` et trois codes d'avertissement (`external_ref_resolved`, `external_ref_unknown`, `external_pedigree_unavailable`) au contrat, la version 2.2.0 le champ `criteria[].trait` et deux codes (`external_records_added`, `external_records_unavailable`) : **n'activez `MELLIFERA_DB` en production qu'après avoir régénéré et publié le SDK PHP** à partir de `openapi/openapi.json`.

## Structure du projet
```

`README.md` : remplacer :

```markdown
### Index Mellifera lineages in the engine

When a queen cites an Index Mellifera lineage by its code (`mother_ref_ext` for her mother, `mating.external_ref` for the mother of the drone-producing queens), the v2 engine reads the crawler database read-only and adds that lineage's ancestry to the pedigree (at most `options.external_pedigree_depth` generations, 7 by default, 0 = disabled, and at most 1,000 ancestors in total). These ancestors only feed relationship and inbreeding; they never appear in any output and their evaluations are never used. A code is looked up as given, then without spaces and case-insensitively (case matters on the site: `DE-ABg-B-0003-2017`; two codes differing only by case make the input ambiguous), then by displayed name (`B4536(NG)`) when unique; an unknown code is ignored (`external_ref_unknown` warning).

Activation: `MELLIFERA_DB=/app/data/mellifera/mellifera.sqlite` in the `api` service's `.env.production`. Without this variable, resolution is disabled and no new warning is emitted; if the database is missing or unreadable, an input citing lineages gets the `external_pedigree_unavailable` warning. Version 2.1.0 adds the `external_pedigree_depth` option and three warning codes (`external_ref_resolved`, `external_ref_unknown`, `external_pedigree_unavailable`) to the contract: **enable `MELLIFERA_DB` in production only after regenerating and publishing the PHP SDK** from `openapi/openapi.json`.

## Project Structure
```

par :

```markdown
### Index Mellifera lineages in the engine

When a queen cites an Index Mellifera lineage by its code (`mother_ref_ext` for her mother, `mating.external_ref` for the mother of the drone-producing queens), the v2 engine reads the crawler database read-only and adds that lineage's ancestry to the pedigree (at most `options.external_pedigree_depth` generations, 7 by default, 0 = disabled, and at most 1,000 ancestors in total). These ancestors only feed relationship and inbreeding; they never appear in any output. Since version 2.2.0, a criterion that declares its standard trait (`criteria[].trait`, a list frozen in the contract) also receives the Index Mellifera records of these ancestors and of the queens recorded in the same tester × year groups (at most 1,000), with no fixed effect shared with the farm's records (see `METHODOLOGIE.md`, §3 bis). A code is looked up as given, then without spaces and case-insensitively (case matters on the site: `DE-ABg-B-0003-2017`; two codes differing only by case make the input ambiguous), then by displayed name (`B4536(NG)`) when unique; an unknown code is ignored (`external_ref_unknown` warning).

Activation: `MELLIFERA_DB=/app/data/mellifera/mellifera.sqlite` in the `api` service's `.env.production`. Without this variable, resolution is disabled and no new warning is emitted; if the database is missing or unreadable, an input citing lineages gets the `external_pedigree_unavailable` warning. Version 2.1.0 adds the `external_pedigree_depth` option and three warning codes (`external_ref_resolved`, `external_ref_unknown`, `external_pedigree_unavailable`) to the contract, version 2.2.0 the `criteria[].trait` field and two codes (`external_records_added`, `external_records_unavailable`): **enable `MELLIFERA_DB` in production only after regenerating and publishing the PHP SDK** from `openapi/openapi.json`.

## Project Structure
```

`METHODOLOGIE.md` : remplacer :

```markdown
  des deux est fournie, l'autre prend sa valeur par défaut ; répétabilité seule : h² = min(0,25,
  répétabilité)), et, pour le profil `breeder` seulement, `queen_ratio` (ρ) et `rg` (rG)
  facultatifs (§6).
- **Profil** : champ `profile` de la requête, `"standard"` (défaut) ou `"breeder"` (§5, §5 bis).

```

par :

```markdown
  des deux est fournie, l'autre prend sa valeur par défaut ; répétabilité seule : h² = min(0,25,
  répétabilité)), et, pour le profil `breeder` seulement, `queen_ratio` (ρ) et `rg` (rG)
  facultatifs (§6) ; caractère standard facultatif (`trait`, §3 bis).
- **Profil** : champ `profile` de la requête, `"standard"` (défaut) ou `"breeder"` (§5, §5 bis).

```

`METHODOLOGIE.md` : remplacer :

```markdown
  `MELLIFERA_DB` est configurée : son ascendance (au plus `options.external_pedigree_depth`
  générations, 7 par défaut, 0 = désactivé ; au plus 1 000 ancêtres au total) entre dans le pedigree
  comme reines de base, sans notes ;
  elle rapproche les reines de l'exploitation issues de lignées apparentées mais n'apparaît dans
  aucune sortie (reines, base 100, accouplements, tendance). Un accouplement `unknown` dont le côté mâle
```

par :

```markdown
  `MELLIFERA_DB` est configurée : son ascendance (au plus `options.external_pedigree_depth`
  générations, 7 par défaut, 0 = désactivé ; au plus 1 000 ancêtres au total) entre dans le pedigree
  comme reines de base (leurs notes : §3 bis) ;
  elle rapproche les reines de l'exploitation issues de lignées apparentées mais n'apparaît dans
  aucune sortie (reines, base 100, accouplements, tendance). Un accouplement `unknown` dont le côté mâle
```

`METHODOLOGIE.md` : remplacer :

```markdown
- un pedigree cyclique (une reine ancêtre d'elle-même par au moins une autre reine) est refusé
  (HTTP 422).

## 4. Parenté
```

par :

```markdown
- un pedigree cyclique (une reine ancêtre d'elle-même par au moins une autre reine) est refusé
  (HTTP 422).

## 3 bis. Notes des lignées Index Mellifera

Un critère peut déclarer son **caractère standard** (`trait`) : une liste de 30 caractères figée dans
le contrat (douceur, essaimage, tenue au cadre, hivernage, récoltes, test d'aiguille, varroa…), qui
indique pour chacun son nom sur Index Mellifera, BeeBreed et le registre Karl Kehrle. Le moteur ne lit
que les notes d'Index Mellifera ; les autres correspondances servent aux exports de Beekube. Deux
critères ne peuvent pas déclarer le même caractère.

Quand un critère déclare un caractère connu d'Index Mellifera et compte au moins une note de
l'exploitation, le moteur lit dans la base du crawler les notes des ancêtres résolus (§3) sur ce
caractère :
- **sources** : tableaux de saison (un groupe par testeur et par année) et radar annuel des fiches
  (un groupe par éleveur et par année) ; le radar d'une année est ignoré quand le tableau de saison
  note la même reine sur le même caractère la même année ;
- **nettoyage**, chaque écart compté dans l'avertissement `external_records_added` : récoltes en kg
  (pas d'échelle commune), notes sans date, hors échelle, non numériques, sans testeur ni éleveur,
  antérieures à la naissance de la reine, reines sans fiche (code inconnu ou année invalide), notes
  seules dans leur groupe ;
- **mise à l'échelle** : chaque note est ramenée sur [0, 1] avec l'échelle du site (0 à 6 pour les
  notes, pourcentages divisés par 100), puis inversée (1 − v) si le sens du critère de l'exploitation
  est contraire à celui du site ; les notes d'une même reine, même jour, même groupe sont moyennées.

**Contemporaines.** Un BLUP n'apprend d'une note que par comparaison avec les autres notes de son
groupe (testeur × année) : la note d'un ancêtre seul dans son groupe est absorbée par l'effet du
groupe et n'apporte rien. Le moteur lit donc aussi les notes des autres reines du même groupe (même
testeur, ou même éleveur pour le radar, même année, même caractère). Ces **contemporaines** entrent au
pedigree comme les ancêtres (liées à leur mère et à leur côté mâle quand ils y sont), au plus 1 000 par
calcul, retenues par code croissant (au-delà : « plafond de 1000 contemporaines atteint »). Inclure les
reines parmi lesquelles les reproductrices ont été choisies rend aussi l'estimation peu sensible à
cette sélection [Henderson 1975].

**Isolement.** Les notes externes ne partagent aucun effet fixe avec celles de l'exploitation : un
rucher et un évaluateur négatifs par testeur (ou éleveur), des classes d'âge et de ruche propres, une
covariable de consanguinité des ouvrières propre. Leur effet d'environnement permanent est estimé
comme celui des reines de l'exploitation. Les a priori de variance (§6), le choix d'estimer par REML
et les compteurs publiés (`n_records`, `n_queens`, `n_contemporary_groups`, nombre de notes,
d'évaluateurs et d'années par reine) portent sur l'exploitation seule. La méthode de fiabilité (§7) se
décide sans les équations dues aux notes externes. Ancêtres et contemporaines n'apparaissent dans
aucune sortie. Les reines de l'exploitation sans parenté avec les lignées externes ni groupe commun
gardent exactement leurs index ; les autres peuvent bouger, comme dans tout BLUP.

Notes illisibles (tables absentes, base d'une version antérieure) : le calcul garde l'ascendance sans
les notes (`external_records_unavailable`).

## 4. Parenté
```

`METHODOLOGIE.md` : remplacer :

```markdown
effets fixes, propagée dans le pedigree [Misztal & Wiggans 1988]. Sur simulation, l'approximation
suit le classement exact (corrélation 0,92 à 0,99) mais **surestime la fiabilité d'environ 0,1**.

Profil `breeder` : r² = 1 − PEV(T) / (a_ii σ²_T) pour la valeur totale T du nœud W publié (σ²_T :
```

par :

```markdown
effets fixes, propagée dans le pedigree [Misztal & Wiggans 1988]. Sur simulation, l'approximation
suit le classement exact (corrélation 0,92 à 0,99) mais **surestime la fiabilité d'environ 0,1**.
Les équations dues aux seules notes Index Mellifera (§3 bis) ne comptent pas dans cette limite : si
l'exploitation et son ascendance tiennent sous 6 000 équations, les PEV de ses reines restent exactes,
lues par résolutions du système déjà factorisé, une par reine, sans inverser la matrice entière.

Profil `breeder` : r² = 1 − PEV(T) / (a_ii σ²_T) pour la valeur totale T du nœud W publié (σ²_T :
```

`METHODOLOGIE.md` : remplacer :

```markdown
  le REML coûte environ 65 s par critère à 3 000 reines, jusqu'à environ 105 s au plafond de 400
  évaluations : à éviter au-delà d'un critère à cette taille.
- Non-régression sur le jeu d'exemple anonymisé (`data/input-example.json`).

```

par :

```markdown
  le REML coûte environ 65 s par critère à 3 000 reines, jusqu'à environ 105 s au plafond de 400
  évaluations : à éviter au-delà d'un critère à cette taille.
- Notes Index Mellifera (§3 bis, `pytest -m slow`, 20 populations simulées : générations 0 à 3 dans
  une base Index Mellifera simulée, 80 % des reines notées par 6 testeurs, génération 4 dans
  l'exploitation, 70 % notées, plus 30 reines sans lien dans leurs propres ruchers ; variances
  fournies). Avec les notes externes, la fiabilité des reines issues des lignées gagne 0,068 en moyenne
  (au moins 0,047 sur chaque population ; seuil : gain partout positif et ≥ 0,05 en moyenne) et la
  corrélation index / vraie valeur passe de 0,472 à 0,517 (seuil : + 0,02) ; les index des reines sans
  lien ne bougent pas (écart ≤ 1e-9). Si la base ne garde que la meilleure moitié des reines externes
  notées (sélection sur la note), le gain de fiabilité tombe à 0,033 et la corrélation de rang passe
  de 0,452 à 0,472 (seuil : pas de baisse). Coût, sur une base de 9 000 reines : 300 reines de
  l'exploitation citant 100 lignées reçoivent 1 789 notes (789 ancêtres, 1 000 contemporaines) ; le
  calcul reste exact et prend environ 0,45 à 0,5 s, contre 0,7 à 0,8 s sans les notes (l'inverse dense
  est alors remplacée par des résolutions).
- Non-régression sur le jeu d'exemple anonymisé (`data/input-example.json`).

```

`METHODOLOGIE.md` : remplacer :

```markdown
- La classe de ruche est la seule information de ruche utilisée (`hive_type` est reçu mais non
  modélisé).

## Références
```

par :

```markdown
- La classe de ruche est la seule information de ruche utilisée (`hive_type` est reçu mais non
  modélisé).
- **Notes Index Mellifera** (§3 bis) : le sens et l'échelle des notes du site ne sont pas vérifiés
  sur des cas réels ; les reines de l'exploitation ne sont jamais comparées aux lignées externes ; le
  gain réel reste à mesurer sur données réelles, par validation croisée, après le crawl complet.

## Références
```

- [ ] **Step 4 : Lancer la suite complète**

```bash
MELLIFERA_CONTACT= .venv/bin/pytest -q -p no:cacheprovider | tail -1
```

Expected :

```text
462 passed, 15 deselected in …
```

Puis ruff :

```bash
uvx ruff check --select F,E9 blup2 tests
```

Expected :

```text
All checks passed!
```

Puis tous les tests lents (environ deux minutes) :

```bash
MELLIFERA_CONTACT= .venv/bin/pytest -q -p no:cacheprovider -m "slow and not docker" | tail -1
```

Expected :

```text
9 passed, 468 deselected in …
```

- [ ] **Step 5 : Commit**

```bash
git add README.fr.md README.md METHODOLOGIE.md tests/simdata.py tests/test_external_calibration.py
git status --short
git commit -q -F - <<'EOF'
Valide les notes Index Mellifera par simulation et documente le volet 3

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015J8yR1HR7mLSRtwnvEYap6
EOF
git log --oneline -1
```

Expected : `git status --short` ne liste que les 5 fichiers ci-dessus, puis :

```text
<hash> Valide les notes Index Mellifera par simulation et documente le volet 3
```

---
