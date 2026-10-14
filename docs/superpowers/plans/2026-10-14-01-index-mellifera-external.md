# Index Mellifera, volet 2 : lignées externes dans le moteur v2 — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Quand une reine de Beekube cite une lignée Index Mellifera par son code (`mother_ref_ext` côté mère, `mating.external_ref` côté mâle), le moteur v2 lit la base du crawler en lecture seule et ajoute l'ascendance de cette lignée au pedigree, pour que la parenté et la consanguinité en tiennent compte, sans que ces ancêtres n'apparaissent dans aucune sortie.

**Architecture:** Nouveau module `blup2/external.py` : `expand_external(inp)` cherche les références dans la vue `queen_mating` de la base du crawler (`MELLIFERA_DB`, ouverte en `mode=ro`), remonte l'ascendance en largeur sur `external_pedigree_depth` générations, construit les ancêtres en `QueenIn` à id négatifs, coupe les liens erronés de la base et renvoie des copies des reines de l'entrée dont les liens sont complétés. `engine._evaluate` l'appelle avant `build_pedigree`, qui reçoit les ancêtres par un nouveau paramètre `extra_queens` et les traite comme des reines de base : ils n'entrent jamais dans `inp.queens`, donc dans aucune des fonctions qui produisent les sorties. Le contrat passe en 2.1.0 (option et trois codes d'avertissement) et le snapshot OpenAPI est régénéré.

**Tech Stack:** Python ≥ 3.11 (3.14 en local, 3.12 en CI et dans l'image `api`) ; `sqlite3` de la bibliothèque standard (URI `mode=ro`) ; pydantic 2.12.3 ; pytest 8.3.3. Aucune dépendance ajoutée.

**Spec:** `docs/superpowers/specs/2026-10-13-index-mellifera-design.md`, §8 (résolution dans le moteur) et §10 (tests du volet 2). Le volet 1 (§1 à §7, §9) est fusionné dans `main` (049decf).

**Prérequis :**
- worktree `/Users/tups/Sites/beekube-blup-external`, branche `mellifera-external`, créée depuis `main` 049decf, sans modification en cours ;
- `.venv` du worktree créé (Python 3.14, `requirements-dev.txt`) ;
- suite actuelle : `364 passed, 12 deselected`, vérifiée à la tâche 1.

**Prototype :** tout le code de ce plan a été écrit et exécuté avant rédaction : suite complète `402 passed, 12 deselected` sous Python 3.12 et 3.14, ruff (F, E9) sans erreur. Des mutations du code (cycles non coupés, type `unknown` conservé, id positifs, profondeur ignorée, erreurs de la base propagées, ancêtres ajoutés à `inp.queens`) font chacune échouer au moins un test. Le plan a ensuite été rejoué de bout en bout dans une copie du dépôt : fichiers extraits du plan, commandes lancées, sorties et commits conformes aux lignes « Expected ». Le code se transcrit tel quel. Un comportement différent de celui annoncé est un défaut du plan, à signaler.

## Écarts / précisions de spec

Points où ce plan précise ou s'écarte de la spec. Chacun est appliqué par la tâche indiquée.

1. **Variable `MELLIFERA_DB` absente : résolution désactivée, sans avertissement** (§8.1 dit : avertissement `external_pedigree_unavailable` si l'entrée contient des références). **Décision à confirmer par l'utilisateur avant la tâche 3.** Raison : la 2.1.0 doit pouvoir être déployée avant la publication du SDK PHP qui connaît les nouveaux codes. Coût si elle est mauvaise : en production, une variable oubliée passe inaperçue (aucun avertissement, résolution silencieusement absente). Sans variable, aucune sortie ne change ; l'avertissement n'est émis que si la variable est définie et que la base est absente ou illisible (fichier manquant, fichier qui n'est pas une base, table ou colonne manquante). Tâches 3 et 4.
2. **Version du moteur 2.1.0** (`ENGINE_VERSION`, `Info(version=)` de l'API, `METHODOLOGIE.md`) : le contrat s'enrichit (option, codes d'avertissement) sans rien retirer. Tâche 1.
3. **`external_pedigree_depth` borné à 30** (`Field(7, ge=0, le=30)`) : la spec ne donne pas de maximum ; 30 générations couvrent tout Index Mellifera et bornent le coût d'une entrée hostile. Tâche 1.
4. **`external_ref` normalisé au code en base.** Quand le côté mâle est résolu, `mating.external_ref` de la copie devient le code trouvé en base (`B276(CBP)`, ` fr-cbp-b-0276-2022` → `FR-CBP-B-0276-2022`). Le nœud sire est indexé par `external_ref` : sans cette normalisation, deux saisies différentes de la même lignée donneraient deux pères distincts, et le père d'une reine de l'exploitation différerait de celui des ancêtres. Tâche 3.
5. **Côté mâle : accouplement `natural` ignoré, `drone_queens_mother` saisi respecté.** La règle de §8.2 (« `drone_queens` vide et `external_ref` trouvé ») est restreinte aux accouplements non naturels sans `drone_queens_mother`. Une fécondation naturelle n'a pas de père identifiable (règle v2 du pedigree), et une mère des RPM déjà connue dans l'exploitation prime sur la référence externe. Un type `station` ou `insemination` est conservé ; seul `unknown` devient `insemination`. Tâche 3.
6. **Recherche d'un code** : code normalisé, puis code tel que saisi (sans les espaces de bord), puis nom affiché unique (`LIMIT 2`). Les codes hors format du crawler (`queen.code_valid = 0`) restent ainsi trouvables. Les lignes sont lues dans la vue `queen_mating` (colonnes explicites) et mises en cache par code. Tâche 3.
7. **Ancêtre sans année valide** (hors [`BORN_MIN`, année courante + 1]) : il n'est pas ajouté, la remontée s'arrête à lui et il est compté dans `external_ref_resolved` (« j reine(s) sans année de naissance valide écartée(s) »). Une référence qui pointe directement vers une telle reine reste non rattachée, sans être comptée comme introuvable. Tâche 3.
8. **Coupure des liens erronés.** Auto-références coupées à la construction de l'ancêtre ; cycles coupés par `find_cycle` (de `blup2/pedigree.py`) sur les seuls ancêtres, tant qu'il en reste : la mère d'abord, puis la mère des RPM (avec `external_ref`, qui identifierait sinon le même père), sinon tout l'accouplement. Un cycle « entre un ancêtre et une reine Beekube » (§8.2) ne se forme pas par les liens mère ou mère des RPM, puisque les ancêtres ne citent que des ancêtres ; la relecture finale a montré qu'il pouvait se former par le nœud sire d'un `external_ref` saisi tel quel par une reine de l'exploitation qui a sa propre mère des RPM (relecture C1) : `build_pedigree` donne depuis aux accouplements des ancêtres, et à ceux qui s'y rattachent, une clé sire à part (`xext`), et le moteur refait le calcul sans les ancêtres si le pedigree reste invalide. Les liens coupés sont comptés dans `external_ref_resolved`. Tâche 3.
9. **Id négatifs** : ils commencent sous le plus petit id déjà utilisé par l'entrée (reines, mères, RPM, mères des RPM), et non à −1, car `BlupV2Input` accepte des id négatifs. Tâche 3.
10. **Messages des avertissements** :
    - `external_ref_resolved` : « N référence(s) Index Mellifera rattachée(s), M ancêtre(s) ajouté(s) au pedigree[, k lien(s) erroné(s) de la base coupé(s)][, j reine(s) sans année de naissance valide écartée(s)] » ;
    - `external_ref_unknown` : « N code(s) Index Mellifera introuvable(s), ignoré(s) jusqu'au prochain crawl : … » (20 au plus, liste de `listed_ids`) ;
    - `external_pedigree_unavailable` : « base Index Mellifera illisible : N référence(s) externe(s) non résolue(s) ». Le chemin et l'erreur ne vont qu'au journal (`logger.warning`), pas dans la réponse.
    - Tâches 3 et 4.
11. **Isolation dans les tests** : `tests/conftest.py` retire `MELLIFERA_DB` de l'environnement de chaque test, pour qu'une variable exportée dans le shell du développeur ne change aucun calcul. Tâche 4.
12. **Contrat avec la vraie base du crawler** : un test construit la base par `mellifera.db.open_db` et `mellifera.ingest.ingest` (faux site du volet 1) et la lit avec `expand_external`, y compris pendant qu'une transaction d'écriture est ouverte (WAL). La mini-base des autres tests (`tests/external_base.py`) ne reprend que les colonnes lues. Tâche 3.

## Global Constraints

- Activation (§8.1) : « Variable `MELLIFERA_DB` (chemin de `mellifera.sqlite`) dans le service `api`, base ouverte en **lecture seule** (`mode=ro`). » Option « `options.external_pedigree_depth` (entier ≥ 0, défaut 7 ; 0 = désactivé) », bornée à 30 (écart 3).
- Entrée intacte (§8.2) : « Elle ne modifie pas `inp` : elle renvoie une **copie** des reines dont les liens sont complétés, la liste des ancêtres ajoutés et des avertissements. »
- Robustesse (§8.2) : « **Les erreurs de saisie d'Index Mellifera ne font jamais échouer le calcul**. » Une base illisible non plus : calcul sans résolution.
- Ancêtres (§8.2) : « `available_as_mother = available_as_drone_source = false` » ; « Les évaluations Index Mellifera ne sont **jamais** ajoutées. »
- Isolation (§8.3) : ancêtres passés « **uniquement** à `build_pedigree` (nouveau paramètre `extra_queens`) » ; « mêmes reines en sortie, avec ou sans base Index Mellifera ; même `base.birth_years` et même `base.n_queens` ; aucun id négatif dans `matings`, `queens` ou `trend` ; l'avertissement `unknown_parent_added` ne compte pas les ancêtres Index Mellifera ».
- Contrat (§8.4) : « Le snapshot OpenAPI est régénéré, et **`MELLIFERA_DB` n'est activé en production qu'après régénération et publication du SDK PHP**. »
- `blup2/` n'importe rien de `mellifera/` (l'image `api` n'a ni httpx ni BeautifulSoup) : le contrat avec le crawler tient au nom de la vue `queen_mating` et de ses colonnes.
- Tests : aucun test n'appelle le vrai site. Noms de fichiers de test uniques dans le dépôt ; jamais de `__init__.py` dans `tests/` ni dans `tests/mellifera/`. Les helpers s'importent par `from tests.external_base import …` et `from tests.mellifera.fake_site import …`.
- Commandes lancées depuis la racine du worktree : `.venv/bin/pytest …`, `.venv/bin/python …`.
- Git :
  - travailler dans le worktree `/Users/tups/Sites/beekube-blup-external`, branche `mellifera-external` ;
  - **toujours** `git add <chemins listés>`, jamais `git add -A` ni `git add .` ;
  - pas de push ;
  - messages en français (« Ajoute … »), terminés par `Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>` puis `Claude-Session: https://claude.ai/code/session_015J8yR1HR7mLSRtwnvEYap6`.
- Ne rien modifier dans le checkout `/Users/tups/Sites/beekube-blup` (branche `main`) : `ROADMAP.fr.md` et `docs/` y portent des modifications sans rapport. Ne pas toucher au Pi (crawl en cours).
- Langue : messages, commentaires, `METHODOLOGIE.md` et README FR en français, avec les accents ; README EN en anglais.

## Review Focus

1. **Même lignée saisie de deux façons** (`B276(CBP)` sur une reine, ` fr-cbp-b-0276-2022` sur une autre). Attendu : une seule reine en base, un seul nœud sire partagé, la même mère pour les deux. Tests : tâche 3, `test_drone_side_ref_typed_differently_shares_the_sire_node_of_the_lineage` et `test_shared_ancestor_is_added_once`.
2. **Base en cours d'écriture par le crawl hebdomadaire** pendant un calcul de l'API. Attendu : lecture immédiate du dernier état validé (WAL, `mode=ro`), jamais de « database is locked » ni d'échec du calcul. Test : tâche 3, `test_engine_reads_while_the_crawler_writes`.
3. **Base d'un ancien schéma ou corrompue** (volume monté vide, migration ratée, fichier tronqué). Attendu : calcul identique à un calcul sans référence, plus l'avertissement `external_pedigree_unavailable` ; jamais d'erreur 500. Test : tâche 4, `test_unavailable_base_gives_the_same_computation` (cas `garbage`, `no_table`, `old_schema`).
4. **Données erronées d'Index Mellifera** (reine mère d'elle-même, cycle entre trois reines, année de naissance absente ou 3022). Attendu : liens coupés, comptés dans `external_ref_resolved`, calcul réussi. Tests : tâche 3, `test_cycles_and_self_references_in_the_base_are_cut` et `test_ancestor_without_valid_year_stops_the_climb` ; tâche 4, `test_errors_in_the_base_never_fail_the_computation`.
5. **Utilisateur qui saisit des id négatifs** (autorisés par `BlupV2Input`). Attendu : aucun ancêtre ne reprend un id de l'entrée. Test : tâche 3, `test_negative_ids_never_collide_with_input_ids`.

---

## Carte des fichiers

| Fichier | Rôle | Tâche |
|---|---|---|
| `blup2/schemas.py` | option `external_pedigree_depth`, trois codes d'avertissement | 1 |
| `blup2/__init__.py`, `app.py`, `METHODOLOGIE.md` | version 2.1.0 | 1 |
| `openapi/openapi.json` | snapshot régénéré | 1 |
| `tests/test_schemas.py`, `tests/test_engine.py`, `tests/test_api.py`, `tests/fixtures/expected-example-v2.json` | contrat 2.1.0 | 1 |
| `blup2/pedigree.py` | `build_pedigree(…, extra_queens=())` | 2 |
| `tests/test_pedigree.py` | ancêtres hors de `queen_ids` | 2 |
| `blup2/external.py` | résolution, ascendance, coupures, avertissements | 3 |
| `tests/external_base.py` | mini-base Index Mellifera des tests | 3 |
| `tests/test_external.py` | règles de §8.2 | 3 |
| `tests/mellifera/test_mellifera_external_contract.py` | lecture de la vraie base du crawler | 3 |
| `blup2/engine.py` | appel de `expand_external` avant `build_pedigree` | 4 |
| `tests/conftest.py` | `MELLIFERA_DB` retiré de l'environnement des tests | 4 |
| `tests/test_external_engine.py` | base indisponible, isolation (§8.3), parenté | 4 |
| `METHODOLOGIE.md` §3, `README.fr.md`, `README.md` | documentation, activation | 4 |

---

### Task 1 : Contrat 2.1.0 (option, codes d'avertissement, version, OpenAPI)

Spec §8.1 et §8.4 ; écarts 2 et 3.

**Files:**
- Modify: `blup2/schemas.py` (`OptionsIn`, `WarningOut.code`)
- Modify: `blup2/__init__.py`, `app.py`, `METHODOLOGIE.md` (version)
- Modify: `openapi/openapi.json` (régénéré)
- Test: `tests/test_schemas.py`, `tests/test_engine.py`, `tests/test_api.py`, `tests/fixtures/expected-example-v2.json`

**Interfaces:**
- Consumes : rien.
- Produces : `OptionsIn.external_pedigree_depth: int` (défaut 7, 0 à 30) ; `WarningOut.code` accepte `"external_ref_resolved"`, `"external_ref_unknown"`, `"external_pedigree_unavailable"` ; `blup2.ENGINE_VERSION == "2.1.0"`.

- [ ] **Step 1 : Vérifier le point de départ**

```bash
git branch --show-current
git status --short
.venv/bin/pytest -q | tail -1
```

Expected : `mellifera-external`, aucune ligne de `git status`, puis `364 passed, 12 deselected in …`.

- [ ] **Step 2 : Écrire les tests (en échec)**

Ajouter à la fin de `tests/test_schemas.py`, précédé de deux lignes vides :

```python
def test_external_pedigree_depth_option():
    d = copy.deepcopy(VALID)
    del d["options"]
    assert BlupV2Input.model_validate(d).options.external_pedigree_depth == 7
    off = BlupV2Input.model_validate(mutate(("options", "external_pedigree_depth"), 0))
    assert off.options.external_pedigree_depth == 0
    for bad in (-1, 31, 2.5):
        assert_invalid(mutate(("options", "external_pedigree_depth"), bad), "external_pedigree_depth")


def test_external_warning_codes_are_part_of_the_contract():
    from blup2.schemas import WarningOut

    for code in ("external_ref_resolved", "external_ref_unknown", "external_pedigree_unavailable"):
        assert WarningOut(code=code, message="m").code == code
```

Passer les attentes de version à 2.1.0 (une occurrence par ligne, dans trois fichiers de test) :

```bash
perl -pi -e 's/"2\.0\.0"/"2.1.0"/' tests/test_engine.py tests/test_api.py tests/fixtures/expected-example-v2.json
git diff --stat
```

Expected : `4 files changed, 21 insertions(+), 4 deletions(-)` (dont `tests/test_schemas.py | 17 +`).

- [ ] **Step 3 : Vérifier l'échec**

Run : `.venv/bin/pytest tests/test_schemas.py tests/test_engine.py tests/test_api.py tests/test_regression.py -q | tail -1`
Expected : `5 failed, 47 passed` : les deux nouveaux tests de `test_schemas.py` (option inconnue, code d'avertissement refusé), `test_run_blup_output_shape` et deux tests de `test_api.py` (version `2.0.0` reçue). `test_regression.py` passe.

- [ ] **Step 4 : Ajouter l'option et les codes, passer le moteur en 2.1.0**

`blup2/schemas.py` : remplacer :

```python
    n_suggestions: int = Field(10, ge=0, le=1000)
```

par :

```python
    n_suggestions: int = Field(10, ge=0, le=1000)
    external_pedigree_depth: int = Field(7, ge=0, le=30)   # générations d'ancêtres Index Mellifera ; 0 = désactivé
```

`blup2/schemas.py` : remplacer :

```python
        "criterion_failed",
    ]
```

par :

```python
        "criterion_failed",
        "external_ref_resolved",
        "external_ref_unknown",
        "external_pedigree_unavailable",
    ]
```

```bash
perl -pi -e 's/2\.0\.0/2.1.0/' blup2/__init__.py app.py METHODOLOGIE.md
git diff --stat blup2/__init__.py app.py METHODOLOGIE.md
```

Expected : `3 files changed, 3 insertions(+), 3 deletions(-)`.

- [ ] **Step 5 : Vérifier que le snapshot OpenAPI est périmé**

Run : `.venv/bin/pytest tests/test_schemas.py tests/test_engine.py tests/test_api.py tests/test_regression.py tests/test_openapi.py -q | tail -3`
Expected : `1 failed, 58 passed`, l'échec étant `tests/test_openapi.py::test_committed_openapi_is_up_to_date` (snapshot périmé).

- [ ] **Step 6 : Régénérer le snapshot et lancer la suite**

```bash
.venv/bin/python scripts/export_openapi.py
git diff --stat openapi/openapi.json
.venv/bin/pytest -q | tail -1
```

Expected : `1 file changed, 12 insertions(+), 2 deletions(-)` (version, option, trois codes), puis `366 passed, 12 deselected in …`.

- [ ] **Step 7 : Commit**

```bash
git add blup2/schemas.py blup2/__init__.py app.py METHODOLOGIE.md openapi/openapi.json tests/test_schemas.py tests/test_engine.py tests/test_api.py tests/fixtures/expected-example-v2.json
git status --short
git commit -F - <<'EOF'
Passe le contrat en 2.1.0 : profondeur des lignées Index Mellifera et nouveaux avertissements

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015J8yR1HR7mLSRtwnvEYap6
EOF
```

Expected : `git status --short` ne montre, avant le commit, que les 9 fichiers modifiés (`M`).

---

### Task 2 : Ancêtres externes dans le pedigree (`extra_queens`)

Spec §8.3.

**Files:**
- Modify: `blup2/pedigree.py` (`build_pedigree`)
- Test: `tests/test_pedigree.py`

**Interfaces:**
- Consumes : rien de la tâche 1.
- Produces : `build_pedigree(queens: list[QueenIn], default_n_drones: int, extra_queens: Iterable[QueenIn] = ()) -> tuple[BeePedigree, list[WarningOut]]`. Les reines de `extra_queens` sont des nœuds du graphe, dans `ped.added_queen_ids`, jamais dans `ped.queen_ids` ni dans l'avertissement `unknown_parent_added`.

- [ ] **Step 1 : Écrire le test (en échec)**

Ajouter à la fin de `tests/test_pedigree.py`, précédé de deux lignes vides :

```python
def test_extra_queens_enter_the_graph_but_not_queen_ids():
    """Ancêtres externes (id négatifs) : nœuds du pedigree comme des reines de base, hors de queen_ids et de
    l'avertissement unknown_parent_added."""
    anc = [queen(-1, 2015, mother=-2, mating=mating("insemination", mother=-3, ext="FR-X-B-0003-2013")),
           queen(-2, 2013), queen(-3, 2013)]
    qs = [queen(10, 2020, mother=-1), queen(11, 2020, mother=3257)]
    ped, warnings = build_pedigree(qs, 12, extra_queens=anc)
    assert ped.queen_ids == [10, 11]
    assert sorted(ped.added_queen_ids) == [-3, -2, -1, 3257]
    assert [w.message for w in warnings if w.code == "unknown_parent_added"] == [
        "1 reine(s) citée(s) comme parent mais absente(s) de queens, ajoutée(s) comme reines de base : (3257)"]
    q10 = ped.q_node(10)
    assert ped.dam[q10] == ped.q_node(-1)
    s = ped.sire[q10]                                   # père de 10 = sire de l'accouplement de -1
    assert ped.kind[s] == "S" and ped.dam[s] == ped.q_node(-3)
```

- [ ] **Step 2 : Vérifier l'échec**

Run : `.venv/bin/pytest tests/test_pedigree.py -q | tail -2`
Expected : `1 failed, 14 passed`, l'échec étant `test_extra_queens_enter_the_graph_but_not_queen_ids` (`TypeError` : argument `extra_queens` inattendu).

- [ ] **Step 3 : Ajouter le paramètre**

`blup2/pedigree.py` : remplacer :

```python
from typing import TYPE_CHECKING, Any, Optional
```

par :

```python
from typing import TYPE_CHECKING, Any, Iterable, Optional
```

`blup2/pedigree.py` : remplacer :

```python
def build_pedigree(queens: list["QueenIn"], default_n_drones: int) -> tuple[BeePedigree, list["WarningOut"]]:
    from .schemas import WarningOut

    queens, self_refs = remove_self_references(queens)
    by_id, added, specs, deps, s_info = _graph(queens)
```

par :

```python
def build_pedigree(queens: list["QueenIn"], default_n_drones: int,
                   extra_queens: Iterable["QueenIn"] = ()) -> tuple[BeePedigree, list["WarningOut"]]:
    """extra_queens : ancêtres venus d'une base externe (id négatifs). Ils entrent dans le graphe comme les
    reines de base (added_queen_ids), jamais dans queen_ids ni dans l'avertissement unknown_parent_added."""
    from .schemas import WarningOut

    queens, self_refs = remove_self_references(queens)
    extra = list(extra_queens)
    by_id, added, specs, deps, s_info = _graph(queens + extra)
```

`blup2/pedigree.py` : remplacer :

```python
    ped.added_queen_ids = list(added)
```

par :

```python
    ped.added_queen_ids = list(added) + [q.id for q in extra]
```

- [ ] **Step 4 : Vérifier que les tests passent**

Run : `.venv/bin/pytest tests/test_pedigree.py -q | tail -1 && .venv/bin/pytest -q | tail -1`
Expected : `15 passed`, puis `367 passed, 12 deselected in …`.

- [ ] **Step 5 : Commit**

```bash
git add blup2/pedigree.py tests/test_pedigree.py
git status --short
git commit -F - <<'EOF'
Accepte des ancêtres externes dans le pedigree, hors des reines évaluées

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015J8yR1HR7mLSRtwnvEYap6
EOF
```

Expected : `git status --short` ne montre, avant le commit, que les 2 fichiers modifiés.

---

### Task 3 : Résolution des références dans la base du crawler (`blup2/external.py`)

Spec §8.2 et §8.4 ; écarts 1 et 4 à 10, 12.

**Files:**
- Create: `blup2/external.py`
- Create: `tests/external_base.py`
- Test: `tests/test_external.py`, `tests/mellifera/test_mellifera_external_contract.py`

**Interfaces:**
- Consumes : `OptionsIn.external_pedigree_depth` et les trois codes de `WarningOut` (tâche 1) ; `build_pedigree(…, extra_queens=…)` (tâche 2, dans les tests) ; `find_cycle(queens) -> list[int]` et `listed_ids(ids: list) -> str` (20 au plus, puis « … et N autres ») de `blup2/pedigree.py` (existants) ; vue `queen_mating` de `mellifera/db.py` (colonnes `code, name, birth_year, mother_code, drone_mother_code, n_drone_queens, single_drone, race_description, mating_type`).
- Produces :
  - `blup2.external.ENV_DB = "MELLIFERA_DB"` ;
  - `@dataclass Expansion(queens: list[QueenIn], ancestors: list[QueenIn] = [], warnings: list[WarningOut] = [])` ;
  - `normalize_code(ref: str) -> str` ;
  - `expand_external(inp: BlupV2Input, db_path: Optional[str] = None) -> Expansion` (défaut : `os.environ["MELLIFERA_DB"]`, vide ou absente = désactivé) ;
  - `tests/external_base.py` : `SCHEMA`, `LINEAGE`, `CYCLES`, `make_db(path, rows=LINEAGE) -> str` (chemin de la base créée), `inp(queens, depth=7, records=()) -> BlupV2Input`, `warning(obj, code) -> list[str]` (messages d'un code, sur une `Expansion` ou une sortie du moteur).

- [ ] **Step 1 : Écrire la mini-base et les tests (en échec)**

`tests/external_base.py` (nouveau) :

```python
"""Mini-base Index Mellifera pour les tests du volet 2 (sous-ensemble du schéma de mellifera/db.py)."""
import sqlite3

from blup2.schemas import BlupV2Input, OptionsIn
from tests.builders import criterion

# Colonnes lues par blup2/external.py ; le contrat avec la vraie base du crawler est vérifié par
# tests/mellifera/test_mellifera_external_contract.py.
SCHEMA = """
CREATE TABLE queen (
    id INTEGER PRIMARY KEY AUTOINCREMENT, code TEXT NOT NULL UNIQUE, name TEXT, birth_year INTEGER,
    mother_code TEXT, drone_mother_code TEXT, n_drone_queens INTEGER, mating_code TEXT,
    single_drone INTEGER NOT NULL DEFAULT 0, race_description TEXT);
CREATE INDEX queen_name ON queen (name, birth_year);
CREATE TABLE mating_code (code TEXT PRIMARY KEY, mating_type TEXT NOT NULL, mating_type_override TEXT);
CREATE VIEW queen_mating AS
    SELECT q.*, COALESCE(m.mating_type_override, m.mating_type, 'unknown') AS mating_type
    FROM queen q LEFT JOIN mating_code m ON m.code = q.mating_code;
INSERT INTO mating_code VALUES ('ins', 'insemination', NULL), ('nat', 'natural', NULL), ('xx', 'unknown', NULL);
"""

# Lignée CBP : 0276 ← mère 0100 (× drones 0009 NG) ← mère 0050 ← mère 0010 ; 0277 sœur de 0276.
LINEAGE = [
    # code, année, mère, « Drones », nom, code d'accouplement
    ("FR-CBP-B-0276-2022", 2022, "FR-CBP-B-0100-2020", "FR-NG-B-0009-2018", "B276(CBP)", "ins"),
    ("FR-CBP-B-0277-2022", 2022, "FR-CBP-B-0100-2020", "FR-NG-B-0009-2018", "B277(CBP)", "ins"),
    ("FR-CBP-B-0100-2020", 2020, "FR-CBP-B-0050-2018", None, "B100(CBP)", "nat"),
    ("FR-CBP-B-0050-2018", 2018, "FR-CBP-B-0010-2016", None, None, None),
    ("FR-CBP-B-0010-2016", 2016, None, None, None, None),
    ("FR-NG-B-0009-2018", 2018, None, None, "B9(NG)", "nat"),
    ("FR-ZZ-B-0001-2020", 2020, None, None, "TWIN", None),          # nom porté par deux reines
    ("FR-ZZ-B-0002-2020", 2020, None, None, "TWIN", None),
]


def make_db(path, rows=LINEAGE):
    conn = sqlite3.connect(path)
    conn.executescript(SCHEMA)
    conn.executemany("INSERT INTO queen (code, birth_year, mother_code, drone_mother_code, name, mating_code) "
                     "VALUES (?, ?, ?, ?, ?, ?)", rows)
    conn.commit()
    conn.close()
    return str(path)


def inp(queens, depth=7, records=()):
    return BlupV2Input(exploitation=1, criteria=[criterion(1)], queens=queens, records=list(records),
                       options=OptionsIn(external_pedigree_depth=depth))


def warning(exp_or_out, code):
    return [w.message for w in exp_or_out.warnings if w.code == code]


CYCLES = [
    ("FR-C-B-0001-2020", 2020, "FR-C-B-0002-2019", "FR-C-B-0001-2020", None, "ins"),   # sa propre RPM
    ("FR-C-B-0002-2019", 2019, "FR-C-B-0003-2018", None, None, None),
    ("FR-C-B-0003-2018", 2018, "FR-C-B-0002-2019", None, None, None),                  # cycle 2 ↔ 3
    ("FR-C-B-0004-2017", 2017, "FR-C-B-0004-2017", None, None, None),                  # sa propre mère
]
```

`tests/test_external.py` (nouveau) :

```python
"""Volet 2 Index Mellifera : résolution des références externes (spec IM §8.1, §8.2, §8.4)."""
import pytest

from blup2.external import expand_external, normalize_code
from blup2.pedigree import build_pedigree
from tests.builders import mating, queen
from tests.external_base import CYCLES, inp, make_db, warning


@pytest.fixture
def db(tmp_path):
    return make_db(tmp_path / "mellifera.sqlite")


# ------------------------------------------------------------------ résolution d'une référence
@pytest.mark.parametrize("ref", ["FR-CBP-B-0276-2022", " fr-cbp-b-0276-2022 ", "fr-cbp-b-0276-2022\t", "B276(CBP)"])
def test_mother_ref_resolved_by_code_or_unique_name(db, ref):
    src = inp([queen(1, 2024, mother_ref_ext=ref)])
    exp = expand_external(src, db)
    mother = exp.queens[0].mother
    assert mother is not None and mother < 0
    assert {a.id: a.born for a in exp.ancestors}[mother] == 2022
    assert src.queens[0].mother is None                     # l'entrée n'est jamais modifiée
    assert warning(exp, "external_ref_resolved") == [
        "1 référence(s) Index Mellifera rattachée(s), 5 ancêtre(s) ajouté(s) au pedigree"]


def test_normalize_code():
    assert normalize_code(" fr-cbp-b 0276-2022\n") == "FR-CBP-B0276-2022"


@pytest.mark.parametrize("ref", ["TWIN", "FR-XXX-B-0001-2020", "B276"])
def test_ambiguous_name_or_unknown_code_is_reported_and_ignored(db, ref):
    exp = expand_external(inp([queen(1, 2024, mother_ref_ext=ref)]), db)
    assert exp.queens[0].mother is None and exp.ancestors == []
    assert warning(exp, "external_ref_unknown") == [
        f"1 code(s) Index Mellifera introuvable(s), ignoré(s) jusqu'au prochain crawl : {ref}"]
    assert warning(exp, "external_ref_resolved") == []


def test_unknown_codes_are_listed_at_most_twenty(db):
    qs = [queen(k, 2024, mother_ref_ext=f"FR-XXX-B-{k:04d}-2020") for k in range(1, 26)]
    msg = warning(expand_external(inp(qs), db), "external_ref_unknown")[0]
    assert msg.startswith("25 code(s)") and msg.endswith("… et 5 autres")


def test_known_mother_and_blank_refs_are_left_alone(db):
    qs = [queen(1, 2020), queen(2, 2024, mother=1, mother_ref_ext="FR-CBP-B-0276-2022"),
          queen(3, 2024, mother_ref_ext="  ")]
    exp = expand_external(inp(qs), db)
    assert [q.mother for q in exp.queens] == [None, 1, None]
    assert exp.ancestors == [] and exp.warnings == []


# ------------------------------------------------------------------ côté mâle
def test_drone_side_ref_becomes_an_insemination_with_known_drone_queens_mother(db):
    qs = [queen(1, 2024, mating=mating("unknown", ext="FR-CBP-B-0276-2022")),
          queen(2, 2024, mating=mating("unknown", ext="FR-CBP-B-0276-2022"))]
    exp = expand_external(inp(qs), db)
    m = exp.queens[0].mating
    assert m.type == "insemination" and m.drone_queens_mother < 0
    assert m.external_ref == "FR-CBP-B-0276-2022" and m.n_drone_queens is None
    ped, _ = build_pedigree(exp.queens, 12, extra_queens=exp.ancestors)
    s1, s2 = ped.sire_of_mating(1), ped.sire_of_mating(2)
    assert s1 == s2 and ped.kind[s1] == "S" and ped.S[s1] == 1          # nœud sire partagé, S = 1
    assert ped.dam[s1] == ped.q_node(m.drone_queens_mother)


def test_drone_side_ref_typed_differently_shares_the_sire_node_of_the_lineage(db):
    qs = [queen(1, 2024, mating=mating("unknown", ext=" fr-cbp-b-0276-2022")),
          queen(2, 2024, mating=mating("unknown", ext="B276(CBP)")),
          queen(3, 2024, mother_ref_ext="FR-CBP-B-0276-2022")]
    exp = expand_external(inp(qs), db)
    assert {q.mating.external_ref for q in exp.queens[:2]} == {"FR-CBP-B-0276-2022"}
    ped, _ = build_pedigree(exp.queens, 12, extra_queens=exp.ancestors)
    assert ped.sire_of_mating(1) == ped.sire_of_mating(2)


def test_station_type_is_kept_and_natural_or_known_drone_side_is_left_alone(db):
    qs = [queen(1, 2024, mating=mating("station", ext="FR-CBP-B-0276-2022")),
          queen(2, 2024, mating=mating("natural", ext="FR-CBP-B-0276-2022")),
          queen(3, 2024, mating=mating("insemination", mother=1, ext="FR-CBP-B-0276-2022")),
          queen(4, 2024, mating=mating("insemination", drone_queens=[1], ext="FR-CBP-B-0276-2022"))]
    exp = expand_external(inp(qs), db)
    assert exp.queens[0].mating.type == "station" and exp.queens[0].mating.drone_queens_mother < 0
    assert exp.queens[1:] == qs[1:]


# ------------------------------------------------------------------ ascendance
def test_ancestry_mating_built_like_the_export(db):
    exp = expand_external(inp([queen(1, 2024, mother_ref_ext="FR-CBP-B-0276-2022")]), db)
    by_born = {}
    for a in exp.ancestors:
        by_born.setdefault(a.born, []).append(a)
    a276 = next(a for a in exp.ancestors if a.id == exp.queens[0].mother)
    assert a276.mating.type == "insemination" and a276.mating.external_ref == "FR-NG-B-0009-2018"
    assert a276.mating.drone_queens_mother in {a.id for a in by_born[2018]}
    assert all(not a.available_as_mother and not a.available_as_drone_source for a in exp.ancestors)
    assert all(a.id < 0 for a in exp.ancestors) and len({a.id for a in exp.ancestors}) == len(exp.ancestors)


@pytest.mark.parametrize("depth, n", [(1, 1), (2, 3), (3, 4), (7, 5)])
def test_depth_limits_the_generations_and_cuts_links_beyond(db, depth, n):
    exp = expand_external(inp([queen(1, 2024, mother_ref_ext="FR-CBP-B-0276-2022")], depth=depth), db)
    assert len(exp.ancestors) == n
    ids = {a.id for a in exp.ancestors}
    for a in exp.ancestors:
        assert a.mother is None or a.mother in ids
        assert a.mating.drone_queens_mother is None or a.mating.drone_queens_mother in ids


def test_long_chain_is_cut_at_the_maximum_depth(tmp_path):
    rows = [(f"FR-L-B-{k:04d}-{1960 + k}", 1960 + k, f"FR-L-B-{k - 1:04d}-{1959 + k}" if k else None, None, None, None)
            for k in range(60)]
    db = make_db(tmp_path / "m.sqlite", rows)
    exp = expand_external(inp([queen(1, 2024, mother_ref_ext=rows[-1][0])], depth=30), db)
    assert len(exp.ancestors) == 30 and min(a.born for a in exp.ancestors) == 1990


def test_depth_zero_disables_the_resolution(db):
    src = inp([queen(1, 2024, mother_ref_ext="FR-CBP-B-0276-2022")], depth=0)
    exp = expand_external(src, db)
    assert exp.queens == src.queens and exp.ancestors == [] and exp.warnings == []


def test_shared_ancestor_is_added_once(db):
    qs = [queen(1, 2024, mother_ref_ext="FR-CBP-B-0276-2022"), queen(2, 2024, mother_ref_ext="FR-CBP-B-0277-2022"),
          queen(3, 2024, mother_ref_ext="B276(CBP)")]
    exp = expand_external(inp(qs), db)
    assert len(exp.ancestors) == 6                       # 0276, 0277, 0100, 0009, 0050, 0010
    assert exp.queens[0].mother == exp.queens[2].mother != exp.queens[1].mother
    assert warning(exp, "external_ref_resolved")[0].startswith("3 référence(s)")


def test_negative_ids_never_collide_with_input_ids(db):
    qs = [queen(-1, 2020), queen(-5, 2024, mother=-1, mating=mating("insemination", drone_queens=[-9])),
          queen(1, 2024, mother_ref_ext="FR-CBP-B-0276-2022")]
    exp = expand_external(inp(qs), db)
    assert max(a.id for a in exp.ancestors) == -10


def test_ancestor_without_valid_year_stops_the_climb(tmp_path):
    rows = [("FR-A-B-0003-2022", 2022, "FR-A-B-0002-2020", None, None, None),
            ("FR-A-B-0002-2020", None, "FR-A-B-0001-2018", None, None, None),
            ("FR-A-B-0001-2018", 2018, None, None, None, None),
            ("FR-A-B-0009-2022", 3022, None, None, None, None)]
    db = make_db(tmp_path / "m.sqlite", rows)
    exp = expand_external(inp([queen(1, 2024, mother_ref_ext="FR-A-B-0003-2022"),
                               queen(2, 2024, mother_ref_ext="FR-A-B-0009-2022")]), db)
    assert [a.born for a in exp.ancestors] == [2022] and exp.ancestors[0].mother is None
    assert exp.queens[1].mother is None
    assert warning(exp, "external_ref_resolved") == [
        "1 référence(s) Index Mellifera rattachée(s), 1 ancêtre(s) ajouté(s) au pedigree, "
        "2 reine(s) sans année de naissance valide écartée(s)"]


def test_cycles_and_self_references_in_the_base_are_cut(tmp_path):
    db = make_db(tmp_path / "m.sqlite", CYCLES)
    qs = [queen(1, 2024, mother_ref_ext="FR-C-B-0001-2020", mating=mating("unknown", ext="FR-C-B-0004-2017")),
          queen(2, 2024, mother=1)]
    exp = expand_external(inp(qs), db)
    assert len(exp.ancestors) == 4
    assert warning(exp, "external_ref_resolved") == [
        "2 référence(s) Index Mellifera rattachée(s), 4 ancêtre(s) ajouté(s) au pedigree, "
        "3 lien(s) erroné(s) de la base coupé(s)"]
    build_pedigree(exp.queens, 12, extra_queens=exp.ancestors)       # ne lève pas « pedigree cyclique »
```

`tests/mellifera/test_mellifera_external_contract.py` (nouveau) :

```python
"""Contrat entre la base du crawler (mellifera/db.py) et le moteur (blup2/external.py)."""
from blup2.external import expand_external
from blup2.schemas import BlupV2Input, CriterionIn, MatingIn, QueenIn
from mellifera.db import open_db
from mellifera.ingest import ingest
from mellifera.urls import queen_path
from tests.mellifera.fake_site import Q, queen_page

SEEN = "2026-10-14T10:00:00.000000+00:00"
GRAND = Q("FR-CBP-B-0100-2020", 2020, 145)
DAUGHTER = Q("FR-CBP-B-0276-2022", 2022, 145, mother=GRAND.code, drones="FR-NG-B-0009-2018", number=8)


def crawler_db(path):
    conn = open_db(path)
    for q in (GRAND, DAUGHTER):
        ingest(conn, queen_path(q.code), queen_page(q), SEEN)
    conn.commit()
    return conn


def beekube_input():
    return BlupV2Input(exploitation=1, criteria=[CriterionIn(id=1)], queens=[
        QueenIn(id=1, born=2024, mother_ref_ext=DAUGHTER.code,
                mating=MatingIn(type="unknown", external_ref=GRAND.code))])


def test_engine_reads_the_crawler_database(tmp_path):
    crawler_db(tmp_path / "m.sqlite").close()
    exp = expand_external(beekube_input(), str(tmp_path / "m.sqlite"))
    assert [w.code for w in exp.warnings] == ["external_ref_resolved"]
    daughter = next(a for a in exp.ancestors if a.id == exp.queens[0].mother)
    assert (daughter.born, daughter.mother) == (2022, exp.queens[0].mating.drone_queens_mother)
    assert daughter.mating.type == "insemination" and daughter.mating.n_drone_queens == 8
    assert daughter.mating.external_ref == "FR-NG-B-0009-2018"


def test_engine_reads_while_the_crawler_writes(tmp_path):
    """Calcul pendant un crawl (WAL) : la lecture seule n'attend pas l'écrivain et voit le dernier état validé."""
    writer = crawler_db(tmp_path / "m.sqlite")
    writer.execute("BEGIN IMMEDIATE")
    writer.execute("UPDATE queen SET birth_year = 1800 WHERE code = ?", (DAUGHTER.code,))
    exp = expand_external(beekube_input(), str(tmp_path / "m.sqlite"))
    assert len(exp.ancestors) == 2
    writer.rollback()
    writer.close()
```

- [ ] **Step 2 : Vérifier l'échec**

Run : `.venv/bin/pytest tests/test_external.py tests/mellifera/test_mellifera_external_contract.py -q | tail -3`
Expected : `2 errors` de collecte, `ModuleNotFoundError: No module named 'blup2.external'` pour chacun des deux fichiers.

- [ ] **Step 3 : Écrire le module**

`blup2/external.py` (nouveau) :

```python
"""Références externes Index Mellifera (spec Index Mellifera §8) : mother_ref_ext et mating.external_ref
résolus dans la base du crawler (lecture seule), ascendance ajoutée au pedigree comme reines de base.

L'entrée n'est jamais modifiée : expand_external renvoie des copies des reines dont les liens sont
complétés et, à part, les ancêtres (id négatifs) que seul build_pedigree reçoit. Aucune donnée de la base
ne fait échouer le calcul : liens erronés coupés, base illisible = calcul sans résolution.
"""
from __future__ import annotations

import datetime as dt
import logging
import os
import sqlite3
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from .pedigree import find_cycle, listed_ids
from .schemas import BORN_MIN, BlupV2Input, MatingIn, QueenIn, WarningOut

logger = logging.getLogger(__name__)

ENV_DB = "MELLIFERA_DB"
_COLUMNS = ("code, name, birth_year, mother_code, drone_mother_code, n_drone_queens, single_drone, "
            "race_description, mating_type")


@dataclass
class Expansion:
    queens: list[QueenIn]                                        # reines de l'entrée, liens complétés
    ancestors: list[QueenIn] = field(default_factory=list)      # ancêtres Index Mellifera (id négatifs)
    warnings: list[WarningOut] = field(default_factory=list)


def normalize_code(ref: str) -> str:
    """Code international saisi : espaces retirés, majuscules."""
    return "".join(ref.split()).upper()


def _valid_year(y) -> bool:
    return isinstance(y, int) and BORN_MIN <= y <= dt.date.today().year + 1


def _mother_ref(q: QueenIn) -> Optional[str]:
    if q.mother is None and q.mother_ref_ext and q.mother_ref_ext.strip():
        return q.mother_ref_ext
    return None


def _drone_ref(q: QueenIn) -> Optional[str]:
    """Côté mâle externe : accouplement non naturel, sans RPM ni mère des RPM connues dans l'exploitation."""
    m = q.mating
    if (m is not None and m.type != "natural" and not m.drone_queens and m.drone_queens_mother is None
            and m.external_ref and m.external_ref.strip()):
        return m.external_ref
    return None


def _connect(path: str) -> sqlite3.Connection:
    if not Path(path).is_file():
        raise sqlite3.OperationalError(f"base introuvable : {path}")
    conn = sqlite3.connect(f"file:{Path(path).resolve()}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    return conn


class _Base:
    def __init__(self, conn: sqlite3.Connection):
        self.conn = conn
        self.by_code: dict[str, Optional[sqlite3.Row]] = {}

    def row(self, code: str) -> Optional[sqlite3.Row]:
        if code not in self.by_code:
            self.by_code[code] = self.conn.execute(
                f"SELECT {_COLUMNS} FROM queen_mating WHERE code = ?", (code,)).fetchone()
        return self.by_code[code]

    def find(self, ref: str) -> Optional[sqlite3.Row]:
        """Code exact après normalisation (ou tel que saisi), à défaut nom affiché s'il est unique en base."""
        for code in dict.fromkeys((normalize_code(ref), ref.strip())):
            r = self.row(code)
            if r is not None:
                return r
        rows = self.conn.execute(f"SELECT {_COLUMNS} FROM queen_mating WHERE name = ? LIMIT 2",
                                 (ref.strip(),)).fetchall()
        return rows[0] if len(rows) == 1 else None


def _ancestry(base: _Base, roots: list[str], depth: int) -> tuple[dict[str, sqlite3.Row], int]:
    """Ancêtres atteints en au plus `depth` générations (racines = génération 1), en largeur : une reine
    atteinte par plusieurs chemins n'est ajoutée qu'une fois. Retourne (code → ligne, nombre sans année)."""
    rows: dict[str, sqlite3.Row] = {}
    no_year: set[str] = set()
    frontier = []
    for code in roots:
        r = base.row(code)
        if code in rows or code in no_year or r is None:
            continue
        if _valid_year(r["birth_year"]):
            rows[code] = r
            frontier.append(code)
        else:
            no_year.add(code)
    for _ in range(depth - 1):
        nxt = []
        for code in frontier:
            for parent in (rows[code]["mother_code"], rows[code]["drone_mother_code"]):
                if not parent or parent in rows or parent in no_year:
                    continue
                r = base.row(parent)
                if r is None:
                    continue
                if _valid_year(r["birth_year"]):
                    rows[parent] = r
                    nxt.append(parent)
                else:
                    no_year.add(parent)
        frontier = nxt
    return rows, len(no_year)


def _ancestor(r: sqlite3.Row, ids: dict[str, int]) -> tuple[QueenIn, int]:
    """QueenIn d'un ancêtre (accouplement construit comme à l'export) ; liens hors de l'ensemble coupés,
    auto-références retirées. Retourne (reine, nombre d'auto-références coupées)."""
    code, born = r["code"], r["birth_year"]
    mother, drone = r["mother_code"], r["drone_mother_code"]
    cut = 0
    if mother == code:
        mother, cut = None, cut + 1
    if drone == code:
        drone, cut = None, cut + 1
    mtype = r["mating_type"] or "unknown"
    if drone and mtype == "unknown":
        mtype = "insemination"          # côté mâle connu : sinon le pedigree ignorerait l'accouplement
    mating: dict = {"type": mtype, "year": born}
    if drone:
        mating["external_ref"] = drone  # même nœud sire pour toutes les reines accouplées avec cette lignée
        if drone in ids:
            mating["drone_queens_mother"] = ids[drone]
    if r["n_drone_queens"]:
        mating["n_drone_queens"] = r["n_drone_queens"]
    if r["single_drone"]:
        mating["n_drones"] = 1
    q = QueenIn(id=ids[code], born=born, mother=ids.get(mother) if mother else None, mating=MatingIn(**mating),
                variety=r["race_description"], available_as_mother=False, available_as_drone_source=False)
    return q, cut


def _cut_cycles(ancestors: list[QueenIn]) -> int:
    """Coupe, tant qu'il en reste, un lien d'un ancêtre pris dans un cycle : la mère d'abord, puis la mère des
    RPM, sinon tout l'accouplement. Les ancêtres ne citent que des ancêtres : un cycle ne passe jamais par une
    reine de l'exploitation. Chaque tour retire un lien : la boucle se termine."""
    pos = {q.id: k for k, q in enumerate(ancestors)}
    cut = 0
    while True:
        cycle = find_cycle(ancestors)
        if not cycle:
            return cut
        members = set(cycle)
        cuts = ([(qid, {"mother": None}) for qid in cycle if ancestors[pos[qid]].mother in members]
                + [(qid, {"mating": ancestors[pos[qid]].mating.model_copy(
                    update={"drone_queens_mother": None, "external_ref": None})})
                   for qid in cycle if ancestors[pos[qid]].mating is not None
                   and ancestors[pos[qid]].mating.drone_queens_mother in members]
                + [(qid, {"mating": None}) for qid in cycle if ancestors[pos[qid]].mating is not None])
        if not cuts:                    # impossible : un membre d'un cycle a toujours un lien dans ce cycle
            return cut
        qid, update = cuts[0]
        ancestors[pos[qid]] = ancestors[pos[qid]].model_copy(update=update)
        cut += 1


def expand_external(inp: BlupV2Input, db_path: Optional[str] = None) -> Expansion:
    """Résout les références externes de l'entrée. db_path : chemin de la base (défaut : MELLIFERA_DB)."""
    refs = [(q, _mother_ref(q), _drone_ref(q)) for q in inp.queens]
    n_refs = sum((m is not None) + (d is not None) for _, m, d in refs)
    depth = inp.options.external_pedigree_depth
    out = Expansion(queens=list(inp.queens))
    if n_refs == 0 or depth == 0:
        return out
    path = db_path if db_path is not None else os.environ.get(ENV_DB, "")
    if not path:        # résolution non activée : aucun avertissement (les SDK antérieurs à 2.1.0 l'ignorent)
        return out
    try:
        conn = _connect(path)
        try:
            return _expand(inp, refs, depth, _Base(conn))
        finally:
            conn.close()
    except sqlite3.Error as exc:
        logger.warning("base Index Mellifera illisible (%s) : %s", path, exc)
        out.warnings.append(WarningOut(code="external_pedigree_unavailable", message=(
            f"base Index Mellifera illisible : {n_refs} référence(s) externe(s) non résolue(s)")))
        return out


def _expand(inp: BlupV2Input, refs: list, depth: int, base: _Base) -> Expansion:
    found: dict[str, Optional[str]] = {}            # référence saisie → code en base
    for _, m, d in refs:
        for ref in (m, d):
            if ref is not None and ref not in found:
                r = base.find(ref)
                found[ref] = r["code"] if r is not None else None
    roots = sorted({c for c in found.values() if c is not None})
    rows, no_year = _ancestry(base, roots, depth)

    used = {i for q in inp.queens for i in (q.id, q.mother) if i is not None}
    used |= {i for q in inp.queens if q.mating is not None
             for i in (*q.mating.drone_queens, q.mating.drone_queens_mother) if i is not None}
    start = min(0, min(used, default=0))            # id négatifs sans collision avec ceux de l'entrée
    ids = {code: start - k for k, code in enumerate(rows, start=1)}

    ancestors, cut = [], 0
    for r in rows.values():
        q, c = _ancestor(r, ids)
        ancestors.append(q)
        cut += c
    cut += _cut_cycles(ancestors)

    queens, linked = [], 0
    for q, m, d in refs:
        update = {}
        if m is not None and found[m] in ids:
            update["mother"] = ids[found[m]]
            linked += 1
        if d is not None and found[d] in ids:
            mtype = "insemination" if q.mating.type == "unknown" else q.mating.type
            update["mating"] = q.mating.model_copy(update={       # code en base : même nœud sire que les ancêtres
                "drone_queens_mother": ids[found[d]], "type": mtype, "external_ref": found[d]})
            linked += 1
        queens.append(q.model_copy(update=update) if update else q)

    out = Expansion(queens=queens, ancestors=ancestors)
    if any(c is not None for c in found.values()):
        msg = f"{linked} référence(s) Index Mellifera rattachée(s), {len(ancestors)} ancêtre(s) ajouté(s) au pedigree"
        if cut:
            msg += f", {cut} lien(s) erroné(s) de la base coupé(s)"
        if no_year:
            msg += f", {no_year} reine(s) sans année de naissance valide écartée(s)"
        out.warnings.append(WarningOut(code="external_ref_resolved", message=msg))
    unknown = sorted(ref for ref, code in found.items() if code is None)
    if unknown:
        out.warnings.append(WarningOut(code="external_ref_unknown", message=(
            f"{len(unknown)} code(s) Index Mellifera introuvable(s), ignoré(s) jusqu'au prochain crawl : "
            f"{listed_ids(unknown)}")))
    return out
```

- [ ] **Step 4 : Vérifier que les tests passent**

Run : `.venv/bin/pytest tests/test_external.py tests/mellifera/test_mellifera_external_contract.py -q | tail -1 && .venv/bin/pytest -q | tail -1`
Expected : `26 passed`, puis `393 passed, 12 deselected in …`.

- [ ] **Step 5 : Commit**

```bash
git add blup2/external.py tests/external_base.py tests/test_external.py tests/mellifera/test_mellifera_external_contract.py
git status --short
git commit -F - <<'EOF'
Résout les lignées Index Mellifera citées par code dans la base du crawler

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015J8yR1HR7mLSRtwnvEYap6
EOF
```

Expected : `git status --short` ne montre, avant le commit, que les 4 fichiers ajoutés (`A`).

---

### Task 4 : Branchement dans le moteur, isolation des ancêtres, documentation

Spec §8.1, §8.3 et §8.4 ; écarts 1 et 11.

**Files:**
- Modify: `blup2/engine.py` (`_evaluate`)
- Create: `tests/conftest.py`
- Test: `tests/test_external_engine.py`
- Modify: `METHODOLOGIE.md` (§3), `README.fr.md`, `README.md`

**Interfaces:**
- Consumes : `expand_external(inp) -> Expansion` (tâche 3) ; `build_pedigree(…, extra_queens=…)` (tâche 2) ; `tests/external_base.py` (`CYCLES`, `inp`, `make_db`, `warning`, tâche 3).
- Produces : `run_blup` et `run_matings` résolvent les références quand `MELLIFERA_DB` est définie ; sorties inchangées en dehors de la parenté.

- [ ] **Step 1 : Écrire les tests (en échec)**

`tests/conftest.py` (nouveau) :

```python
import pytest


@pytest.fixture(autouse=True)
def _no_mellifera_db(monkeypatch):
    """Un MELLIFERA_DB exporté dans le shell du développeur ne doit pas changer les calculs des tests."""
    monkeypatch.delenv("MELLIFERA_DB", raising=False)
```

`tests/test_external_engine.py` (nouveau) :

```python
"""Volet 2 Index Mellifera : résolution dans le moteur, isolation des ancêtres (spec IM §8.1, §8.3)."""
import sqlite3

import pytest

from blup2.engine import run_blup, run_matings
from blup2.schemas import MatingPairIn, MatingV2Input
from tests.builders import criterion, mating, queen, record
from tests.external_base import CYCLES, inp, make_db, warning


@pytest.fixture
def db(tmp_path):
    return make_db(tmp_path / "mellifera.sqlite")


# ------------------------------------------------------------------ base absente ou illisible
@pytest.mark.parametrize("make", ["unset", "missing", "garbage", "no_table", "old_schema"])
def test_unavailable_base_gives_the_same_computation(tmp_path, monkeypatch, make):
    """Base illisible : avertissement, calcul identique. Variable absente : résolution non activée, aucun
    avertissement (la 2.1.0 peut tourner avant la publication du SDK qui connaît les nouveaux codes)."""
    path = tmp_path / "m.sqlite"
    if make == "garbage":
        path.write_bytes(b"pas une base sqlite" * 100)
    elif make == "no_table":
        sqlite3.connect(path).execute("CREATE TABLE other (x)").connection.close()
    elif make == "old_schema":           # queen_mating sans une colonne lue par le moteur
        make_db(path)
        conn = sqlite3.connect(path)
        conn.executescript("DROP VIEW queen_mating; CREATE VIEW queen_mating AS SELECT id, code, name, birth_year, "
                           "mother_code, drone_mother_code, 'unknown' AS mating_type FROM queen;")
        conn.close()
    if make != "unset":
        monkeypatch.setenv("MELLIFERA_DB", str(path))
    qs, recs = family(refs=False)
    out = run_blup(inp(qs + [queen(9, 2024, mother_ref_ext="FR-CBP-B-0276-2022")], records=recs))
    expected = [] if make == "unset" else [
        "base Index Mellifera illisible : 1 référence(s) externe(s) non résolue(s)"]
    assert [w.message for w in out.warnings if w.code.startswith("external")] == expected
    monkeypatch.delenv("MELLIFERA_DB", raising=False)
    ref = run_blup(inp(qs + [queen(9, 2024)], records=recs))
    assert [q.model_dump() for q in out.queens] == [q.model_dump() for q in ref.queens]


def test_no_reference_no_warning_even_without_base(monkeypatch):
    monkeypatch.delenv("MELLIFERA_DB", raising=False)
    qs, recs = family(refs=False)
    out = run_blup(inp(qs, records=recs))
    assert not [w for w in out.warnings if w.code.startswith("external")]


# ------------------------------------------------------------------ isolation dans le moteur (§8.3)
def family(refs=True):
    """Exploitation de 5 reines ; avec refs, 4 et 5 citent des lignées Index Mellifera (3 références)."""
    ext = dict(m4=dict(mother_ref_ext="FR-CBP-B-0276-2022"),
               m5=dict(mother_ref_ext="FR-CBP-B-0277-2022", mating=mating("unknown", ext="FR-CBP-B-0100-2020")))
    qs = [queen(1, 2021, mating=mating("station", station=1, n_rpm=8)), queen(2, 2022, mother=1),
          queen(3, 2022, mother=1), queen(4, 2023, **(ext["m4"] if refs else {})),
          queen(5, 2023, **(ext["m5"] if refs else {}))]
    vals = {1: 0.4, 2: 0.6, 3: 0.7, 4: 0.8, 5: 0.5}
    return qs, [record(k, v, date=f"{2022 + (k > 3)}-05-01") for k, v in vals.items()]


def test_ancestors_never_reach_outputs_base_or_suggestions(db, monkeypatch):
    qs, recs = family()
    monkeypatch.setenv("MELLIFERA_DB", db)
    with_db = run_blup(inp(qs, records=recs))
    monkeypatch.delenv("MELLIFERA_DB")
    without = run_blup(inp(qs, records=recs))
    assert [q.id for q in with_db.queens] == [q.id for q in without.queens] == [1, 2, 3, 4, 5]
    assert with_db.base == without.base
    assert {(m.mother, m.drone_source) for m in with_db.matings} == {(m.mother, m.drone_source) for m in without.matings}
    assert all(i > 0 for m in with_db.matings for i in (m.mother, m.drone_source))
    assert all(p.birth_year in (2021, 2022, 2023) for pts in with_db.trend.values() for p in pts)
    assert [w.code for w in with_db.warnings if w.code == "unknown_parent_added"] == []
    assert warning(with_db, "external_ref_resolved") == [
        "3 référence(s) Index Mellifera rattachée(s), 6 ancêtre(s) ajouté(s) au pedigree"]


def test_queens_of_sister_lineages_become_related(db, monkeypatch):
    qs = [queen(1, 2024, mother_ref_ext="FR-CBP-B-0276-2022", mating=mating("station", station=1)),
          queen(2, 2024, mother_ref_ext="FR-CBP-B-0277-2022", mating=mating("station", station=2))]
    pair = MatingV2Input(exploitation=1, criteria=[criterion(1)], queens=qs,
                         pairs=[MatingPairIn(mother=1, drone_source=2)])
    monkeypatch.delenv("MELLIFERA_DB", raising=False)
    assert run_matings(pair).matings[0].expected_inbreeding == 0.0
    monkeypatch.setenv("MELLIFERA_DB", db)
    assert run_matings(pair).matings[0].expected_inbreeding > 0.0


def test_errors_in_the_base_never_fail_the_computation(tmp_path, monkeypatch):
    monkeypatch.setenv("MELLIFERA_DB", make_db(tmp_path / "m.sqlite", CYCLES))
    qs = [queen(1, 2024, mother_ref_ext="FR-C-B-0001-2020", mating=mating("unknown", ext="FR-C-B-0004-2017")),
          queen(2, 2024, mother=1)]
    out = run_blup(inp(qs))
    assert {q.id for q in out.queens} == {1, 2}
    assert warning(out, "external_ref_resolved")[0].endswith("3 lien(s) erroné(s) de la base coupé(s)")
```

- [ ] **Step 2 : Vérifier l'échec**

Run : `.venv/bin/pytest tests/test_external_engine.py -q | tail -3`
Expected : `7 failed, 2 passed` : passent seulement `test_unavailable_base_gives_the_same_computation[unset]` et `test_no_reference_no_warning_even_without_base`, car le moteur n'appelle pas encore la résolution.

- [ ] **Step 3 : Appeler la résolution avant le pedigree**

`blup2/engine.py` : remplacer :

```python
from .connectivity import connectivity
```

par :

```python
from .connectivity import connectivity
from .external import expand_external
```

`blup2/engine.py` : remplacer :

```python
    ped, warnings = build_pedigree(inp.queens, inp.options.default_n_drones)
```

par :

```python
    ext = expand_external(inp)      # ancêtres Index Mellifera : dans le pedigree seulement (spec IM §8.3)
    ped, warnings = build_pedigree(ext.queens, inp.options.default_n_drones, extra_queens=ext.ancestors)
    warnings.extend(ext.warnings)
```

- [ ] **Step 4 : Vérifier que les tests passent**

Run : `.venv/bin/pytest tests/test_external_engine.py -q | tail -1 && .venv/bin/pytest -q | tail -1`
Expected : `9 passed`, puis `402 passed, 12 deselected in …`.

- [ ] **Step 5 : Documenter**

`METHODOLOGIE.md` : insérer, après la ligne `` - un parent cité mais absent des données est ajouté comme reine de base (`unknown_parent_added`) ; `` :

```markdown
- une lignée Index Mellifera citée par son code (`mother_ref_ext` pour la mère, `mating.external_ref`
  pour la mère des RPM, convention d'Index Mellifera) est cherchée dans la base du crawler quand
  `MELLIFERA_DB` est configurée : son ascendance (au plus `options.external_pedigree_depth`
  générations, 7 par défaut, 0 = désactivé) entre dans le pedigree comme reines de base, sans notes ;
  elle rapproche les reines de l'exploitation issues de lignées apparentées mais n'apparaît dans
  aucune sortie (reines, base 100, accouplements, tendance). Un accouplement `unknown` dont le côté mâle
  est résolu devient une insémination. Les erreurs de la base (auto-références, cycles) sont coupées
  sans faire échouer le calcul (`external_ref_resolved`) ; les codes introuvables sont ignorés jusqu'au
  crawl suivant (`external_ref_unknown`) ; base configurée mais absente ou illisible :
  `external_pedigree_unavailable` ;
```

`README.fr.md` : insérer, juste avant `## Structure du projet` :

```markdown
### Lignées Index Mellifera dans le moteur

Quand une reine cite une lignée Index Mellifera par son code (`mother_ref_ext` pour sa mère, `mating.external_ref` pour la mère des reines productrices de mâles), le moteur v2 lit la base du crawler en lecture seule et ajoute l'ascendance de cette lignée au pedigree (au plus `options.external_pedigree_depth` générations, 7 par défaut, 0 = désactivé). Ces ancêtres servent au calcul de parenté et de consanguinité ; ils n'apparaissent dans aucune sortie et leurs évaluations ne sont jamais utilisées. Un code est cherché tel quel, en majuscules sans espaces, puis par nom affiché (`B4536(NG)`) s'il est unique ; un code introuvable est ignoré (avertissement `external_ref_unknown`).

Activation : `MELLIFERA_DB=/app/data/mellifera/mellifera.sqlite` dans `.env.production` du service `api`. Sans cette variable, la résolution est désactivée et aucun nouvel avertissement n'est émis ; si la base est absente ou illisible, une entrée qui cite des lignées reçoit l'avertissement `external_pedigree_unavailable`. La version 2.1.0 ajoute l'option `external_pedigree_depth` et trois codes d'avertissement (`external_ref_resolved`, `external_ref_unknown`, `external_pedigree_unavailable`) au contrat : **n'activez `MELLIFERA_DB` en production qu'après avoir régénéré et publié le SDK PHP** à partir de `openapi/openapi.json`.
```

`README.md` : insérer, juste avant `## Project Structure` :

```markdown
### Index Mellifera lineages in the engine

When a queen cites an Index Mellifera lineage by its code (`mother_ref_ext` for her mother, `mating.external_ref` for the mother of the drone-producing queens), the v2 engine reads the crawler database read-only and adds that lineage's ancestry to the pedigree (at most `options.external_pedigree_depth` generations, 7 by default, 0 = disabled). These ancestors only feed relationship and inbreeding; they never appear in any output and their evaluations are never used. A code is looked up as given, upper-cased without spaces, then by displayed name (`B4536(NG)`) when unique; an unknown code is ignored (`external_ref_unknown` warning).

Activation: `MELLIFERA_DB=/app/data/mellifera/mellifera.sqlite` in the `api` service's `.env.production`. Without this variable, resolution is disabled and no new warning is emitted; if the database is missing or unreadable, an input citing lineages gets the `external_pedigree_unavailable` warning. Version 2.1.0 adds the `external_pedigree_depth` option and three warning codes (`external_ref_resolved`, `external_ref_unknown`, `external_pedigree_unavailable`) to the contract: **enable `MELLIFERA_DB` in production only after regenerating and publishing the PHP SDK** from `openapi/openapi.json`.
```

Run : `.venv/bin/pytest -q | tail -1`
Expected : `402 passed, 12 deselected in …`.

- [ ] **Step 6 : Commit**

```bash
git add blup2/engine.py tests/conftest.py tests/test_external_engine.py METHODOLOGIE.md README.fr.md README.md
git status --short
git commit -F - <<'EOF'
Ajoute l'ascendance Index Mellifera au pedigree du moteur v2, sans l'exposer dans les sorties

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015J8yR1HR7mLSRtwnvEYap6
EOF
```

Expected : `git status --short` ne montre, avant le commit, que les 6 fichiers (`A`/`M`).

---

## Après le plan

- Vérification finale : `.venv/bin/pytest -q` sous 3.14 et sous 3.12 → `402 passed, 12 deselected`.
- Mise en production, dans cet ordre (hors de ce plan) :
  1. déployer l'API 2.1.0 **sans** `MELLIFERA_DB` : aucune sortie ne change ;
  2. régénérer le SDK PHP depuis `openapi/openapi.json`, le publier et l'intégrer dans Beekube ;
  3. ensuite seulement, ajouter `MELLIFERA_DB=/app/data/mellifera/mellifera.sqlite` à `.env.production` du service `api`, dont le volume `./data:/app/data` contient déjà la base du crawler.
- L'aide de Beekube devra préciser que le code saisi côté mâle est celui de la **mère des reines productrices de mâles**, comme sur Index Mellifera (§8.2).
