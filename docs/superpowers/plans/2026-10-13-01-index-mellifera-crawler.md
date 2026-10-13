# Index Mellifera, volet 1 : crawler, base SQLite, export BLUP — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Lire les pages publiques d'Index Mellifera (index-mellifera.org) par un crawl poli, complet une fois puis incrémental chaque semaine, les stocker dans une base SQLite (`data/mellifera/mellifera.sqlite`) et exporter les reines d'éleveurs choisis, avec toute leur ascendance, au format `BlupV2Input` du moteur v2.

**Architecture:** Paquet autonome `mellifera/` dans ce dépôt. Un client HTTP poli (`client.py`, liste blanche `urls.py`) écrit chaque page dans un cache brut compressé (`cache.py`) et la journalise dans la table `page`. Des parseurs purs (`parse/`, octets → dataclasses, BeautifulSoup) alimentent une base SQLite en mode WAL : `db.py` porte le schéma et le référentiel des critères, `ingest.py` les écritures idempotentes et l'identification des reines évaluées. `crawl.py` orchestre le parcours complet, l'incrémental, la fermeture sur les parents, la reprise et la relecture du cache. `export.py` produit un `BlupV2Input` validé en n'important que `blup2.schemas` et `blup2.pedigree` (ni numpy ni scipy). `cli.py` et `schedule.py` exposent `python -m mellifera` et le crawl hebdomadaire, servis par une image Python légère (`Dockerfile-Mellifera`, service Compose `mellifera`). Aucun fichier de `blup2/` n'est modifié : la résolution des références externes dans le moteur est le volet 2 (plan suivant).

**Tech Stack:** Python ≥ 3.11 (3.14 en local, 3.12 dans l'image, 3.12 et 3.14 en CI) ; httpx 0.28.1 (client, `MockTransport` dans les tests) ; beautifulsoup4 4.14.2 (`html.parser`) ; pydantic 2.12.3 (déjà présent) ; `sqlite3`, `zoneinfo`, `fcntl`, `gzip` de la bibliothèque standard ; pytest 8.3.3.

**Spec:** `docs/superpowers/specs/2026-10-13-index-mellifera-design.md`, §1 à §7, §9 et §10 (volet 1). La §8 (résolution dans `blup2`) relève du volet 2.

**Prérequis :**
- worktree `/Users/tups/Sites/beekube-blup-mellifera`, branche `index-mellifera` (créée depuis `main` 5f2170a ; commits de la spec et du `.gitignore` des pages brutes), sans modification en cours ;
- pages brutes téléchargées pendant l'étude dans `tests/mellifera/fixtures/raw/index-mellifera/` (ignorées par git, 14 fichiers, dont les 8 sources des fixtures) : nécessaires à la tâche 3 seulement ;
- Python ≥ 3.11 (`/opt/homebrew/bin/python3.14`) et un accès réseau pour pip ; Docker Compose ≥ 2.24 pour la tâche 13 ;
- suite actuelle : `210 passed, 12 deselected`, vérifiée à la tâche 1.

**Prototype :** tout le code de ce plan a été écrit et exécuté avant rédaction. Résultats : 144 tests du volet 1 ; suite complète `354 passed, 12 deselected` sous Python 3.12 et 3.14 ; image `mellifera` construite et lancée (imports sans numpy, codes de sortie). Le plan lui-même a ensuite été rejoué de bout en bout dans une copie du dépôt : fichiers extraits du plan, commandes lancées, sorties et commits conformes aux lignes « Expected ». Le code se transcrit tel quel. Un comportement différent de celui annoncé est un défaut du plan, à signaler.

## Écarts / précisions de spec

Points où ce plan précise ou s'écarte de la spec. Chacun est appliqué par la tâche indiquée. La spec n'est corrigée que sur deux faits (point 20).

1. **Ordre de résolution des lignes d'évaluation** (§6) : un lien `/queen/<code>/` de la page, de texte identique et à `href` unique, passe **avant** l'expression régulière du code. Sur `breeder_Year/145/2022`, la cellule « FR-CBP-B-0293-2022 » (au format code) pointe vers FR-CBP-B-0193-2022 : le lien fait foi, car c'est le code de l'URL. Ensuite seulement, un texte au format code et sans lien est pris tel quel. Les étapes 3 et 4 (recherche en base) suivent inchangées. Tâches 6 et 8.
2. **Colonne `mode`** (protocole du tableau : `evjp`, `ceta`, `evdk` ; `radar` pour le graphique) dans `evaluation`, et dans la clé d'unicité `(queen_key, source, mode, tester, eval_date, criterion_code)`. Une même cotation figure dans deux protocoles à la même date : D992, le 07/09/2020, a une douceur de 4 en evjp et de 5 en evdk. `eval_date` est non nul (`''` si la date manque), pour que l'unicité tienne. Tâche 7.
3. **Clé d'une ligne non résolue** : `name:<nom>|<année>|<éleveur de la page>`, avec la colonne `page_breeder_id`. Deux « A087 » d'élevages différents ne fusionnent pas, et l'éleveur de la page sert à l'étape 3. Tâches 7 et 8.
4. **Référentiel en deux tables**, au lieu d'une seule table `criterion`. `trait` contient les 22 caractères BLUP (nom, libellé, `kind` note | percentage | harvest, échelle, sens) et vient de `mellifera/traits.csv`. `criterion` associe la colonne d'un protocole à son trait et vient de `mellifera/criteria.csv` (48 lignes). Le critère BLUP exporté est le trait (`criteria[].id` = `trait.id`) : une caractéristique notée dans plusieurs protocoles (`evjp.mildness`, `ceta.mildness`, `radar.mildness`) forme un seul critère. Une colonne inconnue crée un `criterion` sans trait (au lieu de `kind = unknown`), ignoré par l'export. `open_db` recharge les deux CSV à chaque ouverture en écriture. Colonnes ajoutées au schéma de §5 : `country.iso`, `queen.code_valid`, `queen.list_cotation_count`, `page.parse_version` et `mating_code.category_raw`. Tâche 7.
5. **Échelle des notes [0, 6]** (§11) : c'est l'échelle du graphique. Le mapping se limite aux colonnes observées dans les fixtures. Restent non mappés : varNbr, vshPlus, nosemIdx, smr, recap, varVsh, organiza, thriftin, springDev, cropIdx, jelCrop, cellApp, jellyAmt, maladyResist, nosemMspBee et apiAmm. Tâche 7.
6. **Type d'accouplement.** La catégorie est lue d'abord, dans l'attribut `title` de la légende d'une page de saison ou du paragraphe de la fiche : « Instrumental insemination (FR) », « Mating station (FR) ». La description sert de repli. Contient « insemination » → insemination ; « mating station » ou « Belegstelle » → station ; sinon unknown. Ces valeurs sont stockées dans `mating_code (code, description, category_raw, mating_type, mating_type_override)`. `queen` n'a pas de colonne `mating_type` : la vue `queen_mating` calcule `COALESCE(override, type, 'unknown')`, si bien qu'une correction manuelle s'applique sans relecture. Tâches 7 et 8.
7. **Graphique des cotations (radar).** Lignes d'évaluation avec `mode = 'radar'`, `tester = ''` et `tester_breeder_id` nul. À l'export, `apiary` = éleveur de la reine (le site ne dit pas qui a noté) et `evaluator` est nul. Les libellés du graphique sont mappés aux traits des tableaux : frameHolding → framHold, wintering → winter, broodVitality → broVital. Tâches 8 et 10.
8. **Doublons et radar à l'export.** Les enregistrements de même (reine, trait, rucher, date, source) sont fusionnés par moyenne et comptés (cas evjp / evdk du point 2). Une valeur du graphique est écartée si un tableau couvre la même reine, le même trait et la même année. `--sources` vaut `season_table` par défaut. Tâche 10.
9. **Liens de pedigree erronés coupés à l'export** (§7, « auto-référents ou cycliques »). Est aussi coupé le lien vers un **parent plus jeune** que la reine (saisie erronée). Pour les cycles : tant que `find_cycle` en trouve, on coupe un lien de la reine d'id maximal du cycle, la mère d'abord, sinon la mère des RPM. Une reine sans année valide (hors [1900, Y+1]) est exclue, et ses filles la citent en `mother_ref_ext`. Tâche 10.
10. **Ébauches de reines.** Chaque ligne d'une « Queens list » crée une ligne `queen` minimale si elle n'existe pas encore : code, nom, année, éleveur, mère, Drones et série. La fiche la complète ensuite (`parsed_at` non nul = fiche lue). `list_cotation_count` garde le nombre de cotations affiché dans la liste : l'incrémental relit la fiche quand ce nombre change (§4.2). Tâches 8 et 9.
11. **Incrémental.** La première liste lue est celle de l'année courante : sa navigation donne les autres années, avec repli sur les années de saison si elle n'existe pas. Un éleveur est « nouveau » si `breeder.first_seen ≥ crawl_run.started_at` du passage. Les horodatages sont en UTC à la microseconde : à la seconde près, un éleveur vu dans la seconde du début du passage passait pour nouveau. Les parents absents (404) sont redemandés à chaque passage par la fermeture. C'est contraire à la lettre de §4.2 (« jamais relues »), mais voulu au rythme hebdomadaire : un ancêtre ajouté au site est récupéré. Tâche 9.
12. **Modules.** `client.py` remplace `http.py`, qui masquerait le module standard `http`. S'ajoutent `config.py` (réglages), `urls.py` (liste blanche et chemins), `cache.py`, `ingest.py` (écritures, séparées du schéma de `db.py`) et `schedule.py`. Tâches 1, 2, 8 et 11.
13. **Dépendances** dans `requirements-mellifera.txt` (pydantic, httpx, beautifulsoup4), inclus par `requirements-dev.txt`, au lieu de `requirements.txt`. L'image `api` reste inchangée ; l'image `mellifera` (python:3.12-slim) n'a ni R ni numpy. Tâches 1 et 13.
14. **Commandes** : `python -m mellifera <commande>`, sans script `mellifera` installé.
    - Code de sortie 4 ajouté : erreur de configuration (contact vide, variable invalide, base absente, éleveur inconnu ou ambigu, aucune reine).
    - `--breeder` accepte un trigramme (casse indifférente) ou un id Index Mellifera ; il refuse un trigramme porté par plusieurs éleveurs.
    - Le service Compose `mellifera` est rangé sous le profil `mellifera` : il ne démarre que sur demande.
    - Tâches 10, 12 et 13.
15. **Pages tronquées et pannes de parseur.** Une page sans `</html>` dans ses 200 derniers octets est une erreur de parsing, car BeautifulSoup accepterait une demi-page. Toute exception d'un parseur, pas seulement `ParseError`, compte comme une erreur de page : ses écritures sont annulées (Review Focus 2). Seules les erreurs `sqlite3` interrompent le crawl. Tâches 4 et 9.
16. **Politesse précisée.**
    - `MELLIFERA_DELAY` est refusé sous 1 s.
    - `Retry-After` en secondes est respecté jusqu'à 1 h ; donné sous forme de date HTTP, il est ignoré au profit du doublement du délai.
    - Les redirections ne sont suivies que dans la liste blanche, 3 au plus.
    - Un code de reine qui commence par un point est refusé, car httpx normaliserait `/queen/../` en `/`.
    - Tâches 1 et 2.
17. **Codes hors format** (§6, « compté dans le rapport ») : `queen.code_valid = 0`, et le nombre de ces codes en base est publié dans `crawl_run.stats.nonconforming_codes`. `invalid_codes` compte à part les codes impossibles à mettre dans une URL, qui ne sont pas demandés. Tâche 9.
18. **Pas de `import_raw` ni d'`insemination_comment`.** Ni la colonne Comment des tableaux ni le commentaire d'import ne sont stockés, car ils peuvent contenir des noms. Tâches 6 et 7.
19. **Points de §11 tranchés.**
    - Langue : cookie `django_language=en` et `Accept-Language: en` sont tous deux envoyés.
    - Échelle des notes : [0, 6].
    - Années des listes : celles de la navigation de la page.
    - Codes d'accouplement : lus dans la légende `#insemCodeLegend` des pages de saison et sur les fiches.
    - Reste ouvert : les valeurs du graphique sont-elles la moyenne ou la dernière cotation de l'année ? D'où l'exclusion du graphique par défaut à l'export.
    - Tâches 2, 6 et 10.
20. **Spec corrigée** et commitée avec ce plan.
    - §2.1, volume observé : 32 pays, dont « Unknown » (id 0, sans code ISO), et 754 éleveurs.
    - §6 : l'exemple `itle=` est retiré (artefact de lecture, pas un défaut du site).

## Global Constraints

- Accès (§4.3) : « GET uniquement. Liste blanche de chemins : `/all/`, `/countrybyyear/`, `/breederbyyear/`, `/breeder_Year/`, `/queen/`. Toute autre URL (`/rating/`, `/accounts/`, `/exportYearQueens/`, `/api/`…) est refusée **avant** l'envoi. » Jamais de connexion, de POST ni de cookie de session. Ne jamais suivre `pedigreeapis.org`, qui redirige vers un site de jeux d'argent.
- Politesse (§4.3) : « délai de 2 s entre deux requêtes (`MELLIFERA_DELAY`) », jamais moins d'1 s ; « le délai double à chaque erreur (plafond 5 min), `Retry-After` est respecté ; arrêt propre après 5 erreurs consécutives (`MELLIFERA_MAX_ERRORS`). Un 404 n'est pas une erreur ».
- Identification : « `User-Agent: BeekubeCrawler/1.0 (+<MELLIFERA_CONTACT>)` ; le crawl refuse de démarrer si `MELLIFERA_CONTACT` est vide ». Aucune adresse réelle dans le code ni dans les tests : utiliser `c@example.org` et `contact@example.org`.
- Données personnelles (§2.3) : « la base ne garde de l'éleveur que son id Index Mellifera, son trigramme et son pays. Ni nom, ni adresse, ni auteur de cotation. » Ni commentaire d'import, ni colonne Comment. Le cache brut (`data/mellifera/`) et les pages brutes (`tests/mellifera/fixtures/raw/`) ne sont jamais commités ni publiés ; les fixtures commitées sont anonymisées.
- Tests (§10) : « Aucun test n'appelle le vrai site. »
- Parsing (§6) : « Une page qui ne se parse pas : `page.parse_error` renseigné, le crawl continue […] au-delà de 5 % sur un type, la commande sort en erreur ».
- Codes de sortie (§9) : « 0 succès ; 1 arrêt sur erreurs réseau consécutives ; 2 seuil d'échec de parsing dépassé ; 3 verrou déjà pris », plus 4 pour une erreur de configuration (écart 14).
- Export (§7) : « `exploitation` : 0 » ; « `percentage` : valeur / 100, échelle [0, 1] » ; « `harvest` (`honeyAmt`…) : **exclu** » ; « Le fichier de sortie est validé par `BlupV2Input` avant écriture ».
- `blup2/` n'est pas modifié. `mellifera/` n'importe de `blup2` que `schemas` et `pedigree` : l'image `mellifera` n'a ni numpy ni scipy.
- Dépendances épinglées : `pydantic==2.12.3`, `httpx==0.28.1`, `beautifulsoup4==4.14.2`.
- Tests : fichiers `tests/mellifera/test_mellifera_*.py`, dont les noms doivent rester uniques dans le dépôt (pytest importe les modules de test par leur nom de base). Jamais de `__init__.py` dans `tests/` ni dans `tests/mellifera/`. Les helpers s'importent par `from tests.mellifera.fake_site import …`.
- Commandes lancées depuis la racine du worktree : `.venv/bin/pytest …`, `.venv/bin/python …`.
- Git :
  - travailler dans le worktree `/Users/tups/Sites/beekube-blup-mellifera`, branche `index-mellifera` ;
  - **toujours** `git add <chemins listés>`, jamais `git add -A` ni `git add .` ;
  - ne jamais ajouter `tests/mellifera/fixtures/raw/` ni `data/mellifera/` ;
  - pas de push ;
  - messages en français (« Ajoute … »), terminés par `Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>` puis `Claude-Session: https://claude.ai/code/session_015J8yR1HR7mLSRtwnvEYap6`.
- Ne rien modifier dans le checkout `/Users/tups/Sites/beekube-blup` (branche `main`) : `ROADMAP.fr.md` et `docs/` y portent des modifications sans rapport.
- Langue : journaux, messages d'erreur, commentaires et README FR en français, avec les accents ; README EN en anglais.

## Review Focus

1. **Trigramme porté par deux éleveurs de pays différents, ou saisi en minuscules.** Avec environ 750 éleveurs, les collisions de trigrammes sont probables. Attendu :
   - `--breeder CBP` refuse de fusionner deux programmes et liste les candidats (id, pays) ;
   - `--breeder cbp` équivaut à `CBP` ;
   - `--breeder 145` désigne l'éleveur sans ambiguïté.

   Tests : tâche 10, `test_trigram_shared_by_two_breeders_is_refused` et `test_breeder_by_trigram_any_case_or_by_id`.
2. **Page qui fait échouer un parseur autrement que par `ParseError`** (structure inattendue) pendant le crawl initial, qui dure une trentaine d'heures. Attendu : la page est notée en échec (`parse_error` « AttributeError: … »), ses écritures sont annulées et le crawl continue. Sans cela, la reprise relirait la même page en cache et buterait au même endroit, indéfiniment. Test : tâche 9, `test_unexpected_parser_exception_does_not_block_the_crawl`.
3. **`status` ou `export-blup` pendant le crawl initial**, quand le verrou est pris pendant des heures. Attendu : réponse immédiate (lecture seule, WAL). `reparse` et un second crawl sortent en code 3. Test : tâche 12, `test_lock_held_by_a_crawl`.
4. **Variable d'environnement mal saisie** (`MELLIFERA_DELAY=2s`, `MELLIFERA_DELAY=0.2`, `MELLIFERA_TZ=Europe/Pariss`). Attendu : code 4 et un message qui nomme la variable. Jamais de crawl sans délai, ni de trace Python (`ZoneInfoNotFoundError` hérite de `KeyError`). Tests : tâche 1, `test_invalid_settings_name_the_variable` ; tâche 12, `test_invalid_environment_is_a_configuration_error`.
5. **Crawl hebdomadaire en échec dans le service `schedule`** (base verrouillée, bogue). Attendu : l'erreur est journalisée et un nouvel essai a lieu la semaine suivante, sans arrêt de la boucle. Test : tâche 11, `test_run_schedule_survives_a_failing_job`.

---

## Carte des fichiers

| Fichier | Rôle | Tâche |
|---|---|---|
| `requirements-mellifera.txt` | dépendances du crawler (image `mellifera`) | 1 |
| `requirements-dev.txt` | + `-r requirements-mellifera.txt` | 1 |
| `.gitignore`, `.dockerignore` | + `data/mellifera/` (et pages brutes hors des images) | 1 |
| `mellifera/__init__.py` | paquet | 1 |
| `mellifera/config.py` | `Settings` : variables `MELLIFERA_*`, validées | 1 |
| `mellifera/urls.py` | liste blanche, type de page, construction des chemins | 1 |
| `mellifera/cache.py` | cache brut gzip, une version par URL | 1 |
| `mellifera/client.py` | client HTTP poli (délai, ralentissement, redirections) | 2 |
| `scripts/make_mellifera_fixtures.py` | anonymisation des pages brutes en fixtures | 3 |
| `tests/mellifera/fixtures/*.html` | 8 pages réelles anonymisées | 3 |
| `mellifera/parse/__init__.py`, `parse/common.py` | outils communs des parseurs | 4 |
| `mellifera/parse/countries.py`, `parse/country.py` | `/all/`, `/countrybyyear/` | 4 |
| `mellifera/parse/queens_list.py`, `parse/queen.py` | `/breederbyyear/`, `/queen/` | 5 |
| `mellifera/parse/season.py` | `/breeder_Year/` | 6 |
| `mellifera/db.py`, `traits.csv`, `criteria.csv` | schéma SQLite, référentiel, connexions | 7 |
| `mellifera/ingest.py` | écritures idempotentes, identification des reines évaluées | 8 |
| `mellifera/crawl.py` | complet, incrémental, fermeture, reprise, reparse | 9 |
| `mellifera/export.py` | export `BlupV2Input` | 10 |
| `mellifera/schedule.py` | planification hebdomadaire | 11 |
| `mellifera/cli.py`, `mellifera/__main__.py` | `python -m mellifera …`, verrou, codes de sortie | 12 |
| `Dockerfile-Mellifera`, `docker-compose.yml` | image légère, service `mellifera` (profil) | 13 |
| `.github/workflows/tests.yml` | + construction et contrôle de l'image `mellifera` | 13 |
| `README.md`, `README.fr.md` | section « crawler Index Mellifera » | 13 |
| `tests/mellifera/fake_site.py` | faux site (`httpx.MockTransport`), horloge, pages synthétiques, site de démonstration | 2, 4, 5, 6, 8, 9 |
| `tests/mellifera/test_mellifera_config.py` | réglages, chemins, cache | 1 |
| `tests/mellifera/test_mellifera_client.py` | client | 2 |
| `tests/mellifera/test_mellifera_fixtures.py` | anonymisation des fixtures | 3 |
| `tests/mellifera/test_mellifera_parse.py` | outils communs, pays, éleveurs | 4 |
| `tests/mellifera/test_mellifera_parse_queen.py` | listes et fiches | 5 |
| `tests/mellifera/test_mellifera_parse_season.py` | saisons, robustesse de tous les parseurs | 6 |
| `tests/mellifera/test_mellifera_db.py` | schéma, référentiel, connexions | 7 |
| `tests/mellifera/test_mellifera_ingest.py` | ingestion, identification (§6) | 8 |
| `tests/mellifera/test_mellifera_crawl.py` | parcours | 9 |
| `tests/mellifera/test_mellifera_export.py` | export | 10 |
| `tests/mellifera/test_mellifera_schedule.py` | planification | 11 |
| `tests/mellifera/test_mellifera_cli.py` | commandes de bout en bout | 12 |

Aucun fichier de `blup2/`, `app.py`, `openapi/` ni `Dockerfile-Python` n'est modifié.

Nombre de tests attendus (suite par défaut, `.venv/bin/pytest -q`), à partir de `210 passed, 12 deselected` :

| Tâche | Ajout | Total |
|---|---|---|
| 1 | 9 | 219 |
| 2 | 23 | 242 |
| 3 | 14 | 256 |
| 4 | 15 | 271 |
| 5 | 6 | 277 |
| 6 | 25 | 302 |
| 7 | 4 | 306 |
| 8 | 16 | 322 |
| 9 | 7 | 329 |
| 10 | 10 | 339 |
| 11 | 6 | 345 |
| 12 | 9 | 354 |

Total final : **354 passed, 12 deselected**. Si la suite de départ ne compte pas exactement 210 tests, c'est l'ajout par tâche qui fait foi.

---

### Task 1 : Socle du paquet (dépendances, réglages, chemins autorisés, cache brut)

Spec §3, §4.3 et §4.4 ; écarts 12, 13 et 16 ; Review Focus 4.

Rôle de chaque module :
- `Settings` lit les variables `MELLIFERA_*` et les valide dès la construction : `ValueError` qui nomme la variable fautive, délai d'au moins 1 s, fuseau horaire connu.
- `urls.py` est la seule définition de la liste blanche : le client (tâche 2) refuse tout chemin pour lequel `kind_of` renvoie `None`. `KINDS` donne aussi l'ordre de relecture de `reparse` : parents avant enfants.
- Le cache garde une seule version par URL, sous `<sha1[:2]>/<sha1>.html.gz`. L'écriture est atomique (fichier temporaire puis `os.replace`) et le gzip a un `mtime` nul.

**Files:**
- Create: `requirements-mellifera.txt`, `mellifera/__init__.py`, `mellifera/config.py`, `mellifera/urls.py`, `mellifera/cache.py`
- Modify: `requirements-dev.txt`, `.gitignore`, `.dockerignore`
- Test: `tests/mellifera/test_mellifera_config.py`

**Interfaces:**
- Consumes : rien (premier module).
- Produces :
  - `mellifera.config` :
    - `BASE_URL = "https://index-mellifera.org"` et `MIN_DELAY = 1.0` ;
    - `Settings(data_dir: Path, contact: str, delay: float = 2.0, max_errors: int = 5, schedule: str = "sun 02:00", timezone: str = "Europe/Paris", base_url: str = BASE_URL)`, gelée, qui lève `ValueError` si une valeur est invalide ;
    - propriétés `db_path` (`mellifera.sqlite`), `cache_dir` (`cache`) et `lock_path` (`crawl.lock`) ;
    - `Settings.from_env(env: Mapping[str, str] = os.environ) -> Settings`, qui lit `MELLIFERA_DATA` (défaut `data/mellifera`), `MELLIFERA_CONTACT` (espaces retirés), `MELLIFERA_DELAY`, `MELLIFERA_MAX_ERRORS`, `MELLIFERA_SCHEDULE` et `MELLIFERA_TZ`.
  - `mellifera.urls` :
    - `KINDS = ("countries", "country", "queens_list", "season", "queen")` ;
    - `kind_of(path: str) -> Optional[str]` ;
    - `path_ids(path: str) -> tuple[str, ...]` (lève `ValueError` hors de la liste blanche) ;
    - `countries_path()`, `country_path(country_id)`, `queens_list_path(breeder_id, year)` et `season_path(breeder_id, year)` ;
    - `queen_path(code) -> Optional[str]` (`None` si le code n'est pas sûr dans une URL).
  - `mellifera.cache` : `RawCache(root: Path)`, avec `.relative_path(url) -> str`, `.put(url, body: bytes) -> tuple[str, str]` (chemin relatif, sha256 du contenu) et `.get(rel) -> bytes`.

- [ ] **Step 1 : Vérifier le point de départ**

```bash
cd /Users/tups/Sites/beekube-blup-mellifera
git branch --show-current
git status --short
ls tests/mellifera/fixtures/raw/index-mellifera | wc -l
git check-ignore -q tests/mellifera/fixtures/raw/index-mellifera/all_.html && echo "pages brutes ignorées"
```

Expected : `index-mellifera`, aucune ligne de `git status`, `14`, `pages brutes ignorées`. Sinon, s'arrêter et le signaler.

- [ ] **Step 2 : Déclarer les dépendances et les exclusions**

`requirements-mellifera.txt` (nouveau) :

```text
pydantic==2.12.3
httpx==0.28.1
beautifulsoup4==4.14.2
```

`requirements-dev.txt` (remplacer tout le contenu) :

```text
-r requirements.txt
-r requirements-mellifera.txt
pytest==8.3.3
```

`.gitignore` : ajouter à la fin.

```text
# Base et cache du crawler Index Mellifera (le cache contient des noms) : jamais commités
/data/mellifera/
```

`.dockerignore` : ajouter à la fin. L'image `api` copie `tests/`, et les pages brutes, qui contiennent des noms, ne doivent entrer dans aucune image.

```text
data/mellifera/
tests/mellifera/fixtures/raw/
```

- [ ] **Step 3 : Créer l'environnement et vérifier la suite existante**

```bash
/opt/homebrew/bin/python3.14 -m venv .venv
.venv/bin/pip install -q -r requirements-dev.txt
.venv/bin/python -c "import httpx, bs4; print(httpx.__version__, bs4.__version__)"
.venv/bin/pytest -q | tail -1
```

Expected : `0.28.1 4.14.2`, puis `210 passed, 12 deselected in …`. Si `/opt/homebrew/bin/python3.14` n'existe pas, prendre un autre Python ≥ 3.11.

- [ ] **Step 4 : Écrire les tests (en échec)**

`tests/mellifera/test_mellifera_config.py` :

```python
import hashlib

import pytest

from mellifera.cache import RawCache
from mellifera.config import Settings
from mellifera.urls import kind_of, path_ids, queen_path, queens_list_path, season_path


# ------------------------------------------------------------------ réglages
def test_settings_from_env():
    s = Settings.from_env({"MELLIFERA_DATA": "/x", "MELLIFERA_CONTACT": " me@example.org ", "MELLIFERA_DELAY": "3"})
    assert (str(s.db_path), s.contact, s.delay, s.max_errors, s.schedule) == (
        "/x/mellifera.sqlite", "me@example.org", 3.0, 5, "sun 02:00")
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
```

- [ ] **Step 5 : Vérifier l'échec**

Run : `.venv/bin/pytest tests/mellifera/test_mellifera_config.py -q`
Expected : erreur de collecte `ModuleNotFoundError: No module named 'mellifera'`.

- [ ] **Step 6 : Écrire le paquet, les réglages, les chemins et le cache**

`mellifera/__init__.py` :

```python
"""Crawler d'Index Mellifera (index-mellifera.org) : pages publiques → SQLite → export BlupV2Input."""
```

`mellifera/config.py` :

```python
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
```

`mellifera/urls.py` :

```python
"""Chemins des pages publiques lues par le crawler : construction, type de page, liste blanche."""
from __future__ import annotations

import re
from typing import Optional

# Liste blanche : seuls ces chemins sont demandés au site (spec §4.3). Tout le reste est refusé avant l'envoi.
_KINDS = (
    ("countries", re.compile(r"^/all/$")),
    ("country", re.compile(r"^/countrybyyear/(\d+)/$")),
    ("queens_list", re.compile(r"^/breederbyyear/(\d+)/(\d{4})/$")),
    ("season", re.compile(r"^/breeder_Year/(\d+)/(\d{4})/$")),
    ("queen", re.compile(r"^/queen/([A-Za-z0-9()_-][A-Za-z0-9()._-]*)/$")),
)
KINDS = tuple(k for k, _ in _KINDS)        # ordre de relecture (reparse) : parents avant enfants
# Caractères admis dans un code mis en URL ; pas de point en tête (« .. » serait normalisé en « / » par httpx).
_SAFE_CODE = re.compile(r"^[A-Za-z0-9()_-][A-Za-z0-9()._-]*$")


def kind_of(path: str) -> Optional[str]:
    """Type de page d'un chemin autorisé, None s'il est hors liste blanche."""
    for kind, rx in _KINDS:
        if rx.match(path):
            return kind
    return None


def path_ids(path: str) -> tuple[str, ...]:
    """Identifiants contenus dans le chemin (« /breeder_Year/145/2022/ » → ('145', '2022'))."""
    for _, rx in _KINDS:
        m = rx.match(path)
        if m:
            return m.groups()
    raise ValueError(f"chemin hors liste blanche : {path}")


def countries_path() -> str:
    return "/all/"


def country_path(country_id: int) -> str:
    return f"/countrybyyear/{country_id}/"


def queens_list_path(breeder_id: int, year: int) -> str:
    return f"/breederbyyear/{breeder_id}/{year}/"


def season_path(breeder_id: int, year: int) -> str:
    return f"/breeder_Year/{breeder_id}/{year}/"


def queen_path(code: str) -> Optional[str]:
    """None si le code contient des caractères qu'on ne met pas dans une URL (compté par le crawler)."""
    return f"/queen/{code}/" if _SAFE_CODE.match(code or "") else None
```

`mellifera/cache.py` :

```python
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
```

- [ ] **Step 7 : Vérifier que les tests passent**

Run : `.venv/bin/pytest tests/mellifera/test_mellifera_config.py -q && .venv/bin/pytest -q | tail -1`
Expected : `9 passed`, puis `219 passed, 12 deselected`.

- [ ] **Step 8 : Commit**

```bash
git add requirements-mellifera.txt requirements-dev.txt .gitignore .dockerignore mellifera/__init__.py mellifera/config.py mellifera/urls.py mellifera/cache.py tests/mellifera/test_mellifera_config.py
git status --short
git commit -F - <<'EOF'
Ajoute le socle du crawler Index Mellifera : réglages, chemins autorisés, cache brut

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015J8yR1HR7mLSRtwnvEYap6
EOF
```

Expected : `git status --short` ne montre, avant le commit, que les 9 fichiers ajoutés (`A`/`M`), et en particulier rien sous `tests/mellifera/fixtures/raw/`.

---

### Task 2 : Client HTTP poli

Spec §4.3 ; écarts 16 et 19.

Le client fait des GET, un à la fois, avec au moins `delay` secondes entre deux requêtes. Un chemin hors liste blanche lève `ForbiddenPath` avant tout envoi.

| Situation | Comportement |
|---|---|
| 200 | `ok` |
| 404, 410 | `missing` (n'est pas une erreur) |
| 401, 403 | `forbidden` |
| redirection vers `/accounts/login/`, vers un autre hôte ou hors liste blanche | `forbidden`, non suivie |
| redirection dans la liste blanche | suivie, 3 fois au plus |
| 429, 5xx, délai dépassé, connexion refusée | erreur : le délai double (au plus 300 s) ; `Retry-After` en secondes est respecté jusqu'à 1 h ; après `max_errors` erreurs consécutives : `TooManyErrors` |
| autre 4xx | `error` |

Tout succès remet le délai à sa valeur de départ.

Les tests utilisent `httpx.MockTransport` et une horloge simulée (`FakeClock`), qui n'attend jamais vraiment. `tests/mellifera/fake_site.py` est créé ici. Les tâches 4, 5, 6, 8 et 9 lui ajoutent des sections à la fin. Ses imports couvrent déjà ces ajouts : `dataclass` et `field` servent dès la tâche 5.

**Files:**
- Create: `mellifera/client.py`, `tests/mellifera/fake_site.py`
- Test: `tests/mellifera/test_mellifera_client.py`

**Interfaces:**
- Consumes : `mellifera.urls.kind_of` (tâche 1).
- Produces :
  - `mellifera.client` :
    - constantes `MAX_BACKOFF = 300.0`, `MAX_RETRY_AFTER = 3600.0`, `MAX_REDIRECTS = 3` et `LOGIN_PATH = "/accounts/login/"` ;
    - exceptions `MissingContact(ValueError)`, `ForbiddenPath(ValueError)` et `TooManyErrors(RuntimeError)` ;
    - `Response(path, status: Literal["ok", "missing", "forbidden", "error"], http_status: Optional[int], body: bytes = b"")`, gelée ;
    - `PoliteClient(base_url, contact, delay=2.0, max_errors=5, *, transport=None, sleep=time.sleep, clock=time.monotonic, timeout=30.0)`, avec `.get(path) -> Response`, `.close()` et un gestionnaire de contexte.
  - `tests/mellifera/fake_site.py` :
    - `FIXTURES` (dossier des fixtures) et `fixture(name) -> bytes` ;
    - `FakeSite(pages: Optional[dict[str, bytes]] = None)`, dont les attributs `failures[chemin]` (liste de codes HTTP, ou de couples (code, Retry-After), renvoyés avant la vraie page), `redirects[chemin] = cible` et `requests` sont modifiables, avec la propriété `paths` et la méthode `transport() -> httpx.MockTransport` ;
    - `FakeClock`, avec les attributs `now` et `sleeps` et les méthodes `clock()` et `sleep(s)`.

- [ ] **Step 1 : Créer le faux site et l'horloge simulée**

`tests/mellifera/fake_site.py` :

```python
"""Pages synthétiques au balisage d'Index Mellifera et faux site (httpx.MockTransport) pour les tests."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import httpx

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(name: str) -> bytes:
    return (FIXTURES / name).read_bytes()


# ------------------------------------------------------------------ faux site et horloge
class FakeSite:
    """Sert des pages par chemin. `failures[chemin]` : réponses d'erreur renvoyées avant la vraie page, code HTTP
    ou (code, Retry-After)."""

    def __init__(self, pages: Optional[dict[str, bytes]] = None):
        self.pages: dict[str, bytes] = dict(pages or {})
        self.failures: dict[str, list] = {}
        self.redirects: dict[str, str] = {}
        self.requests: list[httpx.Request] = []

    @property
    def paths(self) -> list[str]:
        return [r.url.path for r in self.requests]

    def handler(self, request: httpx.Request) -> httpx.Response:
        self.requests.append(request)
        path = request.url.path
        if request.method != "GET":
            return httpx.Response(405)
        pending = self.failures.get(path)
        if pending:
            status, retry_after = (pending.pop(0), None) if isinstance(pending[0], int) else pending.pop(0)
            return httpx.Response(status, headers={"Retry-After": retry_after} if retry_after else {})
        if path in self.redirects:
            return httpx.Response(302, headers={"Location": self.redirects[path]})
        if path in self.pages:
            return httpx.Response(200, content=self.pages[path])
        return httpx.Response(404)

    def transport(self) -> httpx.MockTransport:
        return httpx.MockTransport(self.handler)


class FakeClock:
    """Horloge simulée : sleep() avance le temps et note chaque attente."""

    def __init__(self):
        self.now = 1000.0
        self.sleeps: list[float] = []

    def clock(self) -> float:
        return self.now

    def sleep(self, seconds: float) -> None:
        self.sleeps.append(round(seconds, 6))
        self.now += seconds
```

- [ ] **Step 2 : Écrire les tests (en échec)**

`tests/mellifera/test_mellifera_client.py` :

```python
import httpx
import pytest

from mellifera.client import ForbiddenPath, MissingContact, PoliteClient, TooManyErrors
from tests.mellifera.fake_site import FakeClock, FakeSite

BASE = "https://index-mellifera.org"


def make_client(site: FakeSite, clock: FakeClock, **kw) -> PoliteClient:
    return PoliteClient(BASE, "contact@example.org", kw.pop("delay", 2.0), kw.pop("max_errors", 5),
                        transport=site.transport(), sleep=clock.sleep, clock=clock.clock, **kw)


def test_refuses_to_start_without_contact():
    with pytest.raises(MissingContact):
        PoliteClient(BASE, "  ")


def test_user_agent_language_cookie_and_get_only():
    site = FakeSite({"/all/": b"<html></html>"})
    resp = make_client(site, FakeClock()).get("/all/")
    assert (resp.status, resp.http_status, resp.body) == ("ok", 200, b"<html></html>")
    req = site.requests[0]
    assert req.method == "GET"
    assert req.headers["user-agent"] == "BeekubeCrawler/1.0 (+contact@example.org)"
    assert "django_language=en" in req.headers["cookie"]


@pytest.mark.parametrize("path", ["/rating/2759", "/accounts/login/", "/exportYearQueens/145/2022/", "/api/queens/",
                                  "/queen/../accounts/", "/queen/../", "/queenTree/FR-CBP-B-0276-2022/3/", "/all",
                                  "/welcome/"])
def test_paths_outside_whitelist_are_refused_before_sending(path):
    site = FakeSite()
    with pytest.raises(ForbiddenPath):
        make_client(site, FakeClock()).get(path)
    assert site.requests == []


def test_delay_between_requests():
    site = FakeSite({"/all/": b"x", "/countrybyyear/8/": b"y"})
    clock = FakeClock()
    client = make_client(site, clock)
    client.get("/all/")
    clock.now += 0.5                        # 0,5 s de traitement entre les deux requêtes
    client.get("/countrybyyear/8/")
    assert clock.sleeps == [1.5]


def test_backoff_doubles_on_429_and_5xx_and_respects_retry_after():
    site = FakeSite({"/all/": b"x"})
    site.failures["/all/"] = [429, 503, 500]
    clock = FakeClock()
    client = make_client(site, clock)
    assert client.get("/all/").status == "ok"
    assert clock.sleeps == [4.0, 8.0, 16.0]     # délai doublé à chaque erreur
    clock.now += 100
    site.pages["/countrybyyear/8/"] = b"y"
    client.get("/countrybyyear/8/")
    assert clock.sleeps[-1:] == [16.0] and len(clock.sleeps) == 3   # succès : délai revenu à 2 s


def test_retry_after_longer_than_backoff_is_respected():
    site = FakeSite({"/all/": b"x"})
    site.failures["/all/"] = [(429, "120")]
    clock = FakeClock()
    make_client(site, clock).get("/all/")
    assert clock.sleeps == [120.0]


def test_retry_after_as_http_date_falls_back_to_doubling():
    site = FakeSite({"/all/": b"x"})
    site.failures["/all/"] = [(503, "Wed, 21 Oct 2026 07:28:00 GMT")]
    clock = FakeClock()
    make_client(site, clock).get("/all/")
    assert clock.sleeps == [4.0]


def test_backoff_capped_at_five_minutes():
    site = FakeSite({"/all/": b"x"})
    site.failures["/all/"] = [503] * 9
    clock = FakeClock()
    make_client(site, clock, max_errors=20).get("/all/")
    assert max(clock.sleeps) == 300.0


def test_stops_after_max_consecutive_errors():
    site = FakeSite({"/all/": b"x"})
    site.failures["/all/"] = [503] * 10
    with pytest.raises(TooManyErrors):
        make_client(site, FakeClock(), max_errors=5).get("/all/")
    assert len(site.requests) == 5


def test_timeouts_count_as_errors():
    clock = FakeClock()

    def boom(request):
        raise httpx.ReadTimeout("lent", request=request)
    client = PoliteClient(BASE, "c", 2.0, 3, transport=httpx.MockTransport(boom), sleep=clock.sleep, clock=clock.clock)
    with pytest.raises(TooManyErrors):
        client.get("/all/")


def test_404_is_missing_not_an_error():
    site = FakeSite()
    client = make_client(site, FakeClock(), max_errors=2)
    for _ in range(3):
        assert client.get("/queen/FR-XX-B-0001-2020/").status == "missing"


def test_redirect_to_login_is_forbidden_and_not_followed():
    site = FakeSite()
    site.redirects["/breederbyyear/145/2022/"] = f"{BASE}/accounts/login/?next=/breederbyyear/145/2022/"
    resp = make_client(site, FakeClock()).get("/breederbyyear/145/2022/")
    assert (resp.status, resp.http_status) == ("forbidden", 302)
    assert site.paths == ["/breederbyyear/145/2022/"]


def test_redirect_inside_whitelist_is_followed():
    site = FakeSite({"/queen/FR-CBP-B-0276-2022/": b"fiche"})
    site.redirects["/queen/fr-cbp-b-0276-2022/"] = "/queen/FR-CBP-B-0276-2022/"
    resp = make_client(site, FakeClock()).get("/queen/fr-cbp-b-0276-2022/")
    assert (resp.path, resp.status, resp.body) == ("/queen/fr-cbp-b-0276-2022/", "ok", b"fiche")


def test_redirect_to_other_host_is_forbidden():
    site = FakeSite()
    site.redirects["/all/"] = "https://pedigreeapis.org/all/"
    assert make_client(site, FakeClock()).get("/all/").status == "forbidden"
    assert len(site.requests) == 1


def test_redirect_chain_is_cut_after_three_hops():
    site = FakeSite()
    for a, b in (("A", "B"), ("B", "C"), ("C", "D"), ("D", "E")):
        site.redirects[f"/queen/{a}/"] = f"/queen/{b}/"
    resp = make_client(site, FakeClock()).get("/queen/A/")
    assert (resp.status, resp.http_status) == ("forbidden", 302)
    assert site.paths == ["/queen/A/", "/queen/B/", "/queen/C/", "/queen/D/"]
```

- [ ] **Step 3 : Vérifier l'échec**

Run : `.venv/bin/pytest tests/mellifera/test_mellifera_client.py -q`
Expected : erreur de collecte `ModuleNotFoundError: No module named 'mellifera.client'`.

- [ ] **Step 4 : Écrire le client**

`mellifera/client.py` :

```python
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
            headers={"User-Agent": self.user_agent, "Accept-Language": "en"},
            cookies={"django_language": "en"},
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
        self._last = self._clock()

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
```

- [ ] **Step 5 : Vérifier que les tests passent**

Run : `.venv/bin/pytest tests/mellifera/test_mellifera_client.py -q && .venv/bin/pytest -q | tail -1`
Expected : `23 passed`, puis `242 passed, 12 deselected`.

- [ ] **Step 6 : Commit**

```bash
git add mellifera/client.py tests/mellifera/fake_site.py tests/mellifera/test_mellifera_client.py
git commit -F - <<'EOF'
Ajoute le client HTTP poli du crawler Index Mellifera

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015J8yR1HR7mLSRtwnvEYap6
EOF
```

---

### Task 3 : Fixtures anonymisées

Spec §2.3 et §10 (« pages réelles téléchargées pendant l'étude, anonymisées »).

Les pages brutes de `tests/mellifera/fixtures/raw/index-mellifera/` contiennent des noms de personnes en trois endroits :
- les auteurs des cotations, dans les liens `/rating/<id>` des fiches ;
- le commentaire d'import des fiches ;
- la colonne Comment des tableaux d'évaluations.

Le script neutralise ces trois zones par expressions régulières et garde le reste de la page octet pour octet, pour que les parseurs soient testés sur le vrai balisage du site. Les tests vérifient la forme de l'anonymisation, sans lister de vrais noms. Les sommes SHA-256 attendues permettent de vérifier que les fixtures régénérées sont exactement celles du prototype.

**Files:**
- Create: `scripts/make_mellifera_fixtures.py`, puis, générés par le script : `tests/mellifera/fixtures/all.html`, `countrybyyear_8.html`, `breederbyyear_145_2022.html`, `queen_FR-CBP-B-0276-2022.html`, `queen_FR-NG-B-4536-2024.html`, `breeder_Year_145_2021.html`, `breeder_Year_145_2022.html` et `breeder_Year_984_2026.html`
- Test: `tests/mellifera/test_mellifera_fixtures.py`

**Interfaces:**
- Consumes : `tests.mellifera.fake_site.FIXTURES` (tâche 2).
- Produces :
  - les 8 fixtures ci-dessus, lues par `fixture(name)` dans les tâches 4 à 8 ;
  - `scripts/make_mellifera_fixtures.py` : `PAGES` (page brute → fixture), `anonymize(html: str) -> str` et `main(argv) -> int`.

- [ ] **Step 1 : Écrire les tests (en échec)**

`tests/mellifera/test_mellifera_fixtures.py` :

```python
"""Les fixtures commitées sont anonymisées : vérifications de forme (aucune liste de vrais noms ici)."""
import re

import pytest
from bs4 import BeautifulSoup

from tests.mellifera.fake_site import FIXTURES

ALL = sorted(p.name for p in FIXTURES.glob("*.html"))


def test_fixtures_present():
    assert ALL == ["all.html", "breeder_Year_145_2021.html", "breeder_Year_145_2022.html",
                   "breeder_Year_984_2026.html", "breederbyyear_145_2022.html", "countrybyyear_8.html",
                   "queen_FR-CBP-B-0276-2022.html", "queen_FR-NG-B-4536-2024.html"]


@pytest.mark.parametrize("name", ALL)
def test_rating_authors_are_anonymous(name):
    soup = BeautifulSoup((FIXTURES / name).read_bytes(), "html.parser")
    for a in soup.find_all("a", href=re.compile(r"^/rating/\d+$")):
        assert a.get_text(" ", strip=True).endswith(" - Testeur anonyme")


@pytest.mark.parametrize("name", [n for n in ALL if n.startswith("queen_")])
def test_import_comment_is_neutral(name):
    html = (FIXTURES / name).read_text(encoding="utf-8")
    m = re.search(r"<strong>Import : </strong>(?:(?!<strong>).)*<strong> -- Comment : </strong>\s*([^<]*?)\s*<",
                  html, re.S)
    assert m is not None and m.group(1) == "No"


@pytest.mark.parametrize("name", [n for n in ALL if n.startswith("breeder_Year_")])
def test_evaluation_comments_are_empty(name):
    soup = BeautifulSoup((FIXTURES / name).read_bytes(), "html.parser")
    tables = soup.find_all("table")
    assert tables
    for table in tables:
        rows = table.find_all("tr")
        header = [c.get_text(strip=True) for c in rows[0].find_all("td")]
        if header[-1:] != ["Comment"]:
            continue
        for tr in rows[1:]:
            assert tr.find_all("td")[-1].get_text(strip=True) == ""
```

- [ ] **Step 2 : Vérifier l'échec**

Run : `.venv/bin/pytest tests/mellifera/test_mellifera_fixtures.py -q`
Expected : `1 failed, 3 skipped` (`test_fixtures_present` échoue sur une liste vide ; les tests paramétrés n'ont encore aucun fichier).

- [ ] **Step 3 : Écrire le script d'anonymisation**

`scripts/make_mellifera_fixtures.py` :

```python
"""Produit les fixtures de test anonymisées à partir des pages brutes d'Index Mellifera.

    .venv/bin/python scripts/make_mellifera_fixtures.py tests/mellifera/fixtures/raw/index-mellifera tests/mellifera/fixtures

Les pages brutes (ignorées par git) contiennent des noms de personnes : auteurs des cotations, commentaire
d'import, commentaires des tableaux d'évaluations. Ces trois zones sont neutralisées ; le reste de la page est
conservé octet pour octet, pour que les tests portent sur le vrai balisage du site.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

PAGES = {   # page brute → fixture commitée
    "all_.html": "all.html",
    "countrybyyear_8_.html": "countrybyyear_8.html",
    "breederbyyear_145_2022.html": "breederbyyear_145_2022.html",
    "queen_FR-CBP-B-0276-2022_.html": "queen_FR-CBP-B-0276-2022.html",
    "queen_FR-NG-B-4536-2024_.html": "queen_FR-NG-B-4536-2024.html",
    "breeder_Year_145_2021.html": "breeder_Year_145_2021.html",
    "breeder_Year_145_2022.html": "breeder_Year_145_2022.html",
    "breeder_Year_984_2026_.html": "breeder_Year_984_2026.html",
}
ANONYMOUS = "Testeur anonyme"
_RATING = re.compile(r'(<a href="/rating/\d+"[^>]*>)((?:(?!</a>).)*)(</a>)', re.S)
_IMPORT_COMMENT = re.compile(r"(<strong>Import : </strong>(?:(?!<strong>).)*<strong> -- Comment : </strong>)"
                             r"((?:(?!<strong>|</p>).)*)", re.S)
_LAST_CELL = re.compile(r'(<td class="td\d+">)((?:(?!</td>).)*)(</td></tr>)', re.S)


def _rating(m: re.Match) -> str:
    parts = [p.strip() for p in re.sub(r"\s+", " ", m.group(2)).split(" - ")]
    return f"{m.group(1)}{' - '.join(parts[:2] + [ANONYMOUS])}{m.group(3)}"


def anonymize(html: str) -> str:
    html = _RATING.sub(_rating, html)
    html = _IMPORT_COMMENT.sub(lambda m: f"{m.group(1)}\n    No\n    ", html)
    return _LAST_CELL.sub(lambda m: f"{m.group(1)}{m.group(3)}", html)   # colonne Comment des tableaux


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print(__doc__)
        return 2
    raw, out = Path(argv[1]), Path(argv[2])
    out.mkdir(parents=True, exist_ok=True)
    for src, dst in PAGES.items():
        (out / dst).write_text(anonymize((raw / src).read_text(encoding="utf-8")), encoding="utf-8")
        print(f"{src} → {out / dst}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
```

- [ ] **Step 4 : Générer les fixtures et vérifier leurs empreintes**

```bash
.venv/bin/python -I scripts/make_mellifera_fixtures.py tests/mellifera/fixtures/raw/index-mellifera tests/mellifera/fixtures
(cd tests/mellifera/fixtures && shasum -a 256 -c - <<'EOF'
1b3fdb7cdb3283c746114eb794605aa1c34d6496eb32fee827e4566279dbe506  all.html
57ee78eba2b31ae8a57d6dd30b83102b7705b885db4ae20011195254888f5d81  breeder_Year_145_2021.html
57443bccba064387bcff4b91f14635bf878b78d42c02860456b591b227379794  breeder_Year_145_2022.html
17d7483921e9b94e3f4282c17f2f95a7c635b06b98d293825cafc9419df5cec8  breeder_Year_984_2026.html
2024771c38240ad6c1abc6ab590541be9f37b7effdf1e33222b46a69f5c9a56a  breederbyyear_145_2022.html
cf469ccd2660b150d4fd90a74574228419d5602fd7a7ff6e209798b426c1bf38  countrybyyear_8.html
b42a2aea9a8e3e1a25c626b17bb46b826cb0c969901b414ad34a82cf0ced67cf  queen_FR-CBP-B-0276-2022.html
da8bf0bb57df874470db9c9beb6d1542eee2ca0baed29f97cde70d01a24b22d5  queen_FR-NG-B-4536-2024.html
EOF
)
```

Expected : 8 lignes `<page brute> → tests/mellifera/fixtures/<fixture>`, puis 8 lignes `<fixture>: OK`.

Une ligne `FAILED` signifie que la page brute n'est pas celle de l'étude : s'arrêter et le signaler. Ne jamais éditer une fixture à la main.

- [ ] **Step 5 : Vérifier que les tests passent et que rien de brut n'est suivi**

```bash
.venv/bin/pytest tests/mellifera/test_mellifera_fixtures.py -q
.venv/bin/pytest -q | tail -1
git status --short --untracked-files=all tests/mellifera scripts
```

Expected : `14 passed`, puis `256 passed, 12 deselected` ; `git status` liste 10 fichiers (`??`) : le script, le test et les 8 fixtures, et **rien** sous `tests/mellifera/fixtures/raw/`.

- [ ] **Step 6 : Commit**

```bash
git add scripts/make_mellifera_fixtures.py tests/mellifera/test_mellifera_fixtures.py tests/mellifera/fixtures/all.html tests/mellifera/fixtures/countrybyyear_8.html tests/mellifera/fixtures/breederbyyear_145_2022.html tests/mellifera/fixtures/queen_FR-CBP-B-0276-2022.html tests/mellifera/fixtures/queen_FR-NG-B-4536-2024.html tests/mellifera/fixtures/breeder_Year_145_2021.html tests/mellifera/fixtures/breeder_Year_145_2022.html tests/mellifera/fixtures/breeder_Year_984_2026.html
git commit -F - <<'EOF'
Ajoute les fixtures anonymisées d'Index Mellifera et leur script

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015J8yR1HR7mLSRtwnvEYap6
EOF
```

---

### Task 4 : Parseurs, outils communs, pays et éleveurs d'un pays

Spec §6 (« fonction pure », « HTML approximatif », repérage par libellés, dates Django) ; écart 15.

Les parseurs ne font ni réseau ni accès à la base. Chacun prend les octets d'une page et rend des dataclasses gelées, ou lève `ParseError`. `make_soup` refuse les pages vides ou tronquées. `labelled_fields` découpe un bloc en champs « `<strong>Libellé</strong> valeur » dans l'ordre du document, quelle que soit l'imbrication : le site laisse des `<p>` non fermés. Un titre `h1` à `h6` clôt un champ, et le texte des `<script>` et `<style>` est ignoré. `field_after(fields, ancre, libellé)` lève l'ambiguïté des libellés répétés (« Comment », « Description », « Number »).

- `/all/` : un bloc par pays, avec lien `/countrybyyear/<id>/`, drapeau `/flags/XX…` et deux `<h3>` (nom, puis « n Breeders »). Le pays « Unknown » (id 0) n'a pas de code ISO.
- `/countrybyyear/<id>/` : une ligne `li.list-group-item` par éleveur, avec lien `/breeder/<id>/` (le trigramme) et liens de saison `/breeder_Year/<id>/<année>/`. L'id du pays vient de l'URL (tâche 8).

**Files:**
- Create: `mellifera/parse/__init__.py`, `mellifera/parse/common.py`, `mellifera/parse/countries.py`, `mellifera/parse/country.py`
- Modify: `tests/mellifera/fake_site.py` (ajout en fin de fichier)
- Test: `tests/mellifera/test_mellifera_parse.py`

**Interfaces:**
- Consumes : fixtures `all.html` et `countrybyyear_8.html` (tâche 3) ; `fixture` (tâche 2).
- Produces :
  - `mellifera.parse.common` :
    - `PARSER_VERSION = 1`, `CODE_RE` et `ParseError(ValueError)` ;
    - `make_soup(body: bytes) -> BeautifulSoup`, `norm(s) -> str` et `node_text(node) -> str` ;
    - `queen_code_from_href(href) -> Optional[str]` et `is_code(s) -> bool` ;
    - `parse_django_date(s) -> Optional[dt.date]`, `parse_number(s) -> Optional[float]` (sans `%`, virgule décimale acceptée) et `parse_int(s) -> Optional[int]` ;
    - `LabelledField(label, text, hrefs: tuple[str, ...])`, `labelled_fields(container) -> list[LabelledField]` et `field_after(fields, anchor, label, stop=None) -> Optional[LabelledField]`.
  - `mellifera.parse.countries` : `CountryEntry(id: int, name: str, iso: Optional[str], n_breeders: Optional[int])` et `parse_countries(body) -> tuple[CountryEntry, ...]`.
  - `mellifera.parse.country` : `BreederEntry(id: int, trigram: str, season_years: tuple[int, ...])` (années décroissantes) et `parse_country(body) -> tuple[BreederEntry, ...]`.
  - `fake_site` : `countries_page([(id, nom, iso ou None, nombre d'éleveurs)]) -> bytes` et `country_page([(id, trigramme, [années de saison])]) -> bytes`.

- [ ] **Step 1 : Ajouter les pages synthétiques pays et éleveurs au faux site**

Ajouter à la fin de `tests/mellifera/fake_site.py`, précédé de deux lignes vides :

```python
# ------------------------------------------------------------------ pages synthétiques (balisage du site)
def _page(body: str) -> bytes:
    return f"<!DOCTYPE html><html><head><title>Index Mellifera</title></head><body>{body}</body></html>".encode()


def countries_page(countries: list[tuple[int, str, Optional[str], int]]) -> bytes:
    """[(id, nom, ISO ou None, nombre d'éleveurs)]"""
    blocks = "".join(
        f'<div class="col-lg-2"><a href="/countrybyyear/{cid}/">'
        f'<img src="/media/flags/{iso or "XX_abc"}.png" alt="flag"></a><h3>{name}</h3><h3>{n} Breeders</h3></div>'
        for cid, name, iso, n in countries)
    return _page(f'<div class="row">{blocks}</div>')


def _season_links(breeder_id: int, years: list[int]) -> str:
    return "".join(f'<a href="/breeder_Year/{breeder_id}/{y}/">{y},</a>' for y in years)


def country_page(breeders: list[tuple[int, str, list[int]]]) -> bytes:
    """[(id éleveur, trigramme, années de saison)]"""
    rows = "".join(
        f'<li class="list-group-item"><div class="row"><div class="col-3"><a href="/breeder/{bid}/">{tri}</a></div>'
        f'<div class="col-7">{_season_links(bid, years)}</div>'
        f'<div class="col-1"><a href="/profile/{bid}/">Profil</a></div></div></li>'
        for bid, tri, years in breeders)
    head = '<li class="list-group-item bg-light"><div class="row"><div class="col-3"><span>Breeder</span></div></div></li>'
    return _page(f'<ul class="list-group">{head}{rows}</ul>')
```

- [ ] **Step 2 : Écrire les tests (en échec)**

`tests/mellifera/test_mellifera_parse.py` :

```python
import datetime as dt

import pytest

from mellifera.parse.common import is_code, labelled_fields, make_soup, parse_django_date
from mellifera.parse.countries import parse_countries
from mellifera.parse.country import parse_country
from tests.mellifera.fake_site import countries_page, country_page, fixture


# ------------------------------------------------------------------ outils communs
@pytest.mark.parametrize("text, expected", [
    ("Sept. 7, 2020", dt.date(2020, 9, 7)), ("June 6, 2024", dt.date(2024, 6, 6)),
    ("Nov. 5, 2021", dt.date(2021, 11, 5)), ("May 3, 2021", dt.date(2021, 5, 3)),
    ("-", None), ("", None), ("Feb. 30, 2021", None),
])
def test_parse_django_date(text, expected):
    assert parse_django_date(text) == expected


@pytest.mark.parametrize("code, ok", [
    ("FR-CBP-B-0276-2022", True), ("LU-PJ-V-005H-2015", True), ("B4536(NG)", False), ("A087", False),
])
def test_is_code(code, ok):
    assert is_code(code) is ok


def test_labelled_fields_ignore_unclosed_paragraphs_and_stop_at_headings():
    soup = make_soup(b"<html><body><div><p><strong>A : </strong>1<p><strong> -- B : </strong> deux "
                     b"<a href='/x/'>lien</a><h1>Details</h1>parasite</div></body></html>")
    fields = labelled_fields(soup.find("div"))
    assert [(f.label, f.text, f.hrefs) for f in fields] == [("A", "1", ()), ("B", "deux lien", ("/x/",))]


# ------------------------------------------------------------------ /all/
def test_countries():
    countries = parse_countries(fixture("all.html"))
    assert len(countries) == 32
    assert sum(c.n_breeders for c in countries) == 754
    by_id = {c.id: c for c in countries}
    assert (by_id[8].name, by_id[8].iso, by_id[8].n_breeders) == ("France", "FR", 85)
    assert by_id[35].iso == "RE"                    # drapeau « RE_HNhsPa6.png »
    assert (by_id[0].name, by_id[0].iso) == ("Unknown", None)


# ------------------------------------------------------------------ /countrybyyear/
def test_country():
    breeders = parse_country(fixture("countrybyyear_8.html"))
    assert len(breeders) == 85
    first = breeders[0]
    assert (first.id, first.trigram, first.season_years) == (1140, "AAD", (2024, 2022, 2021))
    cbp = next(b for b in breeders if b.trigram == "CBP")
    assert cbp.id == 145 and 2022 in cbp.season_years


def test_synthetic_country_pages_parse_like_the_site():
    countries = parse_countries(countries_page([(8, "France", "FR", 2), (0, "Unknown", None, 0)]))
    assert [(c.id, c.name, c.n_breeders) for c in countries] == [(8, "France", 2), (0, "Unknown", 0)]
    breeders = parse_country(country_page([(145, "CBP", [2022, 2021]), (200, "NG", [])]))
    assert [(b.id, b.trigram, b.season_years) for b in breeders] == [(145, "CBP", (2022, 2021)), (200, "NG", ())]
    assert parse_country(country_page([])) == ()           # pays sans éleveur (« Unknown ») : vide, pas une erreur
```

- [ ] **Step 3 : Vérifier l'échec**

Run : `.venv/bin/pytest tests/mellifera/test_mellifera_parse.py -q`
Expected : erreur de collecte `ModuleNotFoundError: No module named 'mellifera.parse'`.

- [ ] **Step 4 : Écrire les outils communs et les deux parseurs**

`mellifera/parse/__init__.py` :

```python
"""Parseurs des pages d'Index Mellifera : fonctions pures, octets HTML → structures."""
```

`mellifera/parse/common.py` :

```python
"""Outils communs aux parseurs : erreur, champs étiquetés, codes de reine, dates et nombres."""
from __future__ import annotations

import datetime as dt
import re
from dataclasses import dataclass
from typing import Optional

from bs4 import BeautifulSoup, Comment, NavigableString, Tag

PARSER_VERSION = 1

# Code international : PAYS-ÉLEVEUR-RACE-NUMÉRO-ANNÉE (FR-CBP-B-0276-2022, LU-PJ-V-005H-2015).
CODE_RE = re.compile(r"^[A-Z]{2}-[A-Z0-9]{1,6}-[A-Z0-9]-[0-9A-Z]{1,6}-\d{4}$")
_QUEEN_HREF = re.compile(r"^/queen/([^/]+)/$")
_HEADINGS = {"h1", "h2", "h3", "h4", "h5", "h6"}
_LABEL_TRIM = re.compile(r"^[\s:\-–]+|[\s:\-–]+$")
_MONTHS = {"jan": 1, "feb": 2, "mar": 3, "march": 3, "apr": 4, "april": 4, "may": 5, "jun": 6, "june": 6,
           "jul": 7, "july": 7, "aug": 8, "sep": 9, "sept": 9, "oct": 10, "nov": 11, "dec": 12}


class ParseError(ValueError):
    """La page ne contient pas la structure attendue."""


def make_soup(body: bytes) -> BeautifulSoup:
    """Page complète uniquement : vide ou sans « </html> » final (réponse tronquée) → ParseError."""
    if not body or not body.strip():
        raise ParseError("page vide")
    if b"</html>" not in body[-200:].lower():
        raise ParseError("page tronquée (pas de </html> final)")
    return BeautifulSoup(body, "html.parser")


def norm(s: Optional[str]) -> str:
    """Espaces (y compris insécables) réduits à un seul, bords retirés."""
    return re.sub(r"\s+", " ", s or "").strip()


def node_text(node: Tag) -> str:
    return norm(node.get_text(" "))


def queen_code_from_href(href: Optional[str]) -> Optional[str]:
    m = _QUEEN_HREF.match(href or "")
    return m.group(1) if m else None


def is_code(s: Optional[str]) -> bool:
    return bool(CODE_RE.match(s or ""))


def parse_django_date(s: Optional[str]) -> Optional[dt.date]:
    """« Sept. 7, 2020 », « June 6, 2024 », « May 3, 2021 » → date ; sinon None."""
    m = re.fullmatch(r"([A-Za-z]+)\.?\s+(\d{1,2}),\s+(\d{4})", norm(s))
    if not m:
        return None
    month = _MONTHS.get(m.group(1).lower())
    if month is None:
        return None
    try:
        return dt.date(int(m.group(3)), month, int(m.group(2)))
    except ValueError:
        return None


def parse_number(s: Optional[str]) -> Optional[float]:
    t = norm(s).replace(",", ".").rstrip("%")
    if t in ("", "-", "None"):
        return None
    try:
        return float(t)
    except ValueError:
        return None


def parse_int(s: Optional[str]) -> Optional[int]:
    v = parse_number(s)
    return int(v) if v is not None and v.is_integer() else None


@dataclass(frozen=True)
class LabelledField:
    """Un <strong>libellé</strong> suivi de sa valeur (texte et liens jusqu'au <strong> suivant)."""

    label: str
    text: str
    hrefs: tuple[str, ...]


def labelled_fields(container: Tag) -> list[LabelledField]:
    """Champs « <strong>Libellé : </strong> valeur » d'un bloc, dans l'ordre du document.

    Indépendant de l'imbrication : les <p> mal fermés du site ne changent pas le découpage.
    Le texte des <script> et <style> est ignoré.
    """
    elems = list(container.descendants)
    starts = [i for i, e in enumerate(elems) if isinstance(e, Tag) and e.name == "strong"]
    out: list[LabelledField] = []
    for k, i in enumerate(starts):
        strong = elems[i]
        own = {id(x) for x in strong.descendants}
        end = starts[k + 1] if k + 1 < len(starts) else len(elems)
        parts: list[str] = []
        hrefs: list[str] = []
        for e in elems[i + 1:end]:
            if isinstance(e, Tag) and e.name in _HEADINGS:
                break                      # un titre de section clôt le champ (« Details », « Cotations »)
            if id(e) in own or isinstance(e, Comment):
                continue
            if isinstance(e, NavigableString):
                if e.parent is not None and e.parent.name in ("script", "style"):
                    continue
                parts.append(str(e))
            elif e.name == "a" and e.get("href"):
                hrefs.append(e["href"])
        label = _LABEL_TRIM.sub("", norm(strong.get_text(" ")))
        out.append(LabelledField(label, _LABEL_TRIM.sub("", norm(" ".join(parts))), tuple(hrefs)))
    return out


def field_after(fields: list[LabelledField], anchor: str, label: str,
                stop: Optional[str] = None) -> Optional[LabelledField]:
    """Premier champ `label` situé après le champ `anchor` (et avant `stop` s'il est donné).

    `label == anchor` renvoie le champ `anchor` lui-même.
    """
    seen = False
    for f in fields:
        if not seen:
            if f.label == anchor:
                seen = True
                if label == anchor:
                    return f
            continue
        if stop is not None and f.label == stop:
            return None
        if f.label == label:
            return f
    return None
```

`mellifera/parse/countries.py` :

```python
"""Page `/all/` : liste des pays."""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Optional

from .common import ParseError, make_soup, node_text, parse_int

_HREF = re.compile(r"^/countrybyyear/(\d+)/$")
_FLAG = re.compile(r"/flags/([A-Za-z]{2})(?:_\w+)?\.\w+$")   # « DE.png », « RE_HNhsPa6.png »


@dataclass(frozen=True)
class CountryEntry:
    id: int
    name: str
    iso: Optional[str]                  # None pour « Unknown » (id 0), sans drapeau
    n_breeders: Optional[int]


def parse_countries(body: bytes) -> tuple[CountryEntry, ...]:
    soup = make_soup(body)
    out: dict[int, CountryEntry] = {}
    for a in soup.find_all("a", href=_HREF):
        cid = int(_HREF.match(a["href"]).group(1))
        box = a.parent
        img = a.find("img")
        m = _FLAG.search(img.get("src", "")) if img is not None else None
        h3 = box.find_all("h3") if box is not None else []
        name = node_text(h3[0]) if h3 else ""
        count = parse_int(node_text(h3[1]).split(" ")[0]) if len(h3) > 1 else None
        out.setdefault(cid, CountryEntry(cid, name, m.group(1).upper() if m else None, count))
    if not out:
        raise ParseError("aucun pays")
    return tuple(out.values())
```

`mellifera/parse/country.py` :

```python
"""Page `/countrybyyear/<pays>/` : éleveurs du pays et leurs saisons de pedigree."""
from __future__ import annotations

import re
from dataclasses import dataclass

from .common import ParseError, make_soup, node_text

_BREEDER = re.compile(r"^/breeder/(\d+)/$")
_SEASON = re.compile(r"^/breeder_Year/(\d+)/(\d{4})/$")


@dataclass(frozen=True)
class BreederEntry:
    id: int
    trigram: str
    season_years: tuple[int, ...]       # années des pages /breeder_Year/<id>/<année>/, décroissantes


def parse_country(body: bytes) -> tuple[BreederEntry, ...]:
    soup = make_soup(body)
    if soup.find("ul") is None and soup.find("li") is None:
        raise ParseError("liste des éleveurs absente")
    out: list[BreederEntry] = []
    for li in soup.find_all("li", class_="list-group-item"):
        a = li.find("a", href=_BREEDER)
        if a is None:
            continue
        years = {int(_SEASON.match(s["href"]).group(2)) for s in li.find_all("a", href=_SEASON)}
        out.append(BreederEntry(int(_BREEDER.match(a["href"]).group(1)), node_text(a), tuple(sorted(years, reverse=True))))
    return tuple(out)
```

- [ ] **Step 5 : Vérifier que les tests passent**

Run : `.venv/bin/pytest tests/mellifera/test_mellifera_parse.py -q && .venv/bin/pytest -q | tail -1`
Expected : `15 passed`, puis `271 passed, 12 deselected`.

- [ ] **Step 6 : Commit**

```bash
git add mellifera/parse/__init__.py mellifera/parse/common.py mellifera/parse/countries.py mellifera/parse/country.py tests/mellifera/fake_site.py tests/mellifera/test_mellifera_parse.py
git commit -F - <<'EOF'
Ajoute les parseurs des pays et des éleveurs d'Index Mellifera

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015J8yR1HR7mLSRtwnvEYap6
EOF
```

---

### Task 5 : Parseurs des listes de reines et des fiches

Spec §2.1 et §2.2 (sémantique du côté mâle, `1dr`), §6 (graphique des cotations) ; écarts 6, 10 et 19.

- **« Queens list »** (`/breederbyyear/<éleveur>/<année>/`) : une ligne `li.list-group-item` par reine, à 7 colonnes (nom, année, mère, Drones, cotations, série, code). L'en-tête `bg-light` est ignoré. Le code est celui du **lien** `/queen/<code>/` : le nom affiché peut différer (FR-CBP-B-0293-2022 affiché, lien vers FR-CBP-B-0193-2022). La page donne aussi l'éleveur (fil d'Ariane `/breeder/<id>/`), son pays et la navigation vers les autres années.
- **Fiche** (`/queen/<code>/`, bloc `div.cadree`) : champs repérés par libellé. « Drones : X -- Number : n » se lit `drone_mother_code = X` et `n_drone_queens = n` (convention Buckfast : n filles de X produisent les mâles). La catégorie d'accouplement est l'attribut `title` du paragraphe « Insemination Code ». La fiabilité du COI est convertie en fraction. Le nombre de cotations compte les liens `/rating/<id>`, qui ne sont jamais suivis. Le graphique vient du JavaScript intégré : `labels: [...]`, puis un jeu `label: 'AAAA' … data: [...]` par année. Les valeurs vides du site (`""`, `No`, `None`, `Unknown`) donnent `None`. Les erreurs de saisie (lien mère 0098, commentaire 0298) sont gardées telles quelles.

**Files:**
- Create: `mellifera/parse/queens_list.py`, `mellifera/parse/queen.py`
- Modify: `tests/mellifera/fake_site.py` (ajout en fin de fichier)
- Test: `tests/mellifera/test_mellifera_parse_queen.py`

**Interfaces:**
- Consumes : `ParseError`, `LabelledField`, `field_after`, `labelled_fields`, `make_soup`, `node_text`, `norm`, `parse_django_date`, `parse_int`, `parse_number` et `queen_code_from_href` (tâche 4) ; fixtures `breederbyyear_145_2022.html`, `queen_FR-CBP-B-0276-2022.html` et `queen_FR-NG-B-4536-2024.html`.
- Produces :
  - `mellifera.parse.queens_list` :
    - `QueenListRow(code, name, birth_year, mother_code, drone_mother_code, cotation_count, serial_raw)` ;
    - `QueensListPage(breeder_id: int, trigram: str, country_id: Optional[int], years: tuple[int, ...], rows: tuple[QueenListRow, ...])` ;
    - `parse_queens_list(body) -> QueensListPage`.
  - `mellifera.parse.queen` :
    - `RadarSeries(year: int, values: dict[str, float])` ;
    - `QueenPage(code, name, birth_year, breeder_id, mother_code, mother_comment, drone_mother_code, n_drone_queens, drone_comment, single_drone: bool, generation_raw, coi, coi_reliability, not_found_raw, serial_raw, insemination_date, mating_code, mating_description, mating_category_raw, race_code, race_description, lineage, origin_raw, data_origin_raw, cotation_count: int, radar: tuple[RadarSeries, ...])` ;
    - `parse_queen(body) -> QueenPage`.
  - `fake_site` :
    - `Q(code, year, breeder_id, name=None, mother=None, drones=None, number=None, cotations=0, mating=("ins", "Instrumental insemination", "Instrumental insemination (ins)"), radar={})`, avec la propriété `label` ;
    - `queens_list_page(breeder_id, trigram, country_id, years, queens) -> bytes` et `queen_page(q: Q) -> bytes`.

- [ ] **Step 1 : Ajouter les reines synthétiques au faux site**

Ajouter à la fin de `tests/mellifera/fake_site.py`, précédé de deux lignes vides :

```python
@dataclass
class Q:
    """Reine synthétique : sert à la fois aux listes et aux fiches."""
    code: str
    year: int
    breeder_id: int
    name: Optional[str] = None
    mother: Optional[str] = None
    drones: Optional[str] = None
    number: Optional[int] = None
    cotations: int = 0
    mating: tuple[str, str, str] = ("ins", "Instrumental insemination", "Instrumental insemination (ins)")
    radar: dict[int, list[float]] = field(default_factory=dict)   # année → 5 valeurs

    @property
    def label(self) -> str:
        return self.name or self.code


def _link(code: Optional[str], empty: str) -> str:
    return f'<a href="/queen/{code}/">{code}</a>' if code else empty


def queens_list_page(breeder_id: int, trigram: str, country_id: int, years: list[int], queens: list[Q]) -> bytes:
    nav = "".join(f'<a href="/breederbyyear/{breeder_id}/{y}/">{y}</a> ' for y in years)
    head = ('<li class="list-group-item bg-light"><div class="row">' + "".join(
        f'<div class="col">{h}</div>' for h in ("Name", "Birth", "Mother", "Drones", "Cotation", "Serial", "Queen"))
        + "</div></li>")
    rows = "".join(
        '<li class="list-group-item"><div class="row">'
        f'<div class="col"><a href="/queen/{q.code}/">{q.label}</a></div><div class="col">{q.year}</div>'
        f'<div class="col">{_link(q.mother, "Unknowne")}</div><div class="col">{_link(q.drones, "Unknown")}</div>'
        f'<div class="col">{q.cotations}</div><div class="col">None</div>'
        f'<div class="col"><a href="/queen/{q.code}/">{q.code}</a></div></div></li>'
        for q in queens)
    crumbs = f'<a href="/countrybyyear/{country_id}/">FR</a> / <a href="/breeder/{breeder_id}/">{trigram}</a>'
    return _page(f"{crumbs}<div>{nav}</div><ul>{head}{rows}</ul>")


def queen_page(q: Q) -> bytes:
    radar = ""
    if q.radar:
        sets = "".join(f"{{ label: '{y}', data: [{', '.join(f'{v:g}' for v in vals)}] }}," for y, vals in q.radar.items())
        radar = ("<script>var marksData = { labels: ['swarming', 'mildness', 'frameHolding', 'wintering', "
                 f"'broodVitality'], datasets: [{sets}] }};</script>")
    code, desc, title = q.mating
    ratings = "".join(f'<a href="/rating/{i + 1}" class="list-group-item">{q.code} - 01/07/2023 - Testeur anonyme</a>'
                      for i in range(q.cotations))
    return _page(
        f'<div class="cadree"><h1>{q.label}</h1>'
        f"<p><strong>International Code : </strong>{q.code}</p>"
        f"<p><strong>Birth Year : </strong>{q.year}</p>"
        f"<p><strong>Mother : </strong>{_link(q.mother, 'Unknown')}</p>"
        f"<p><strong>Drones : </strong>{_link(q.drones, 'Unknown')}"
        f"<strong> -- Number : </strong>{q.number if q.number is not None else 'None'}</p>"
        f'<p><strong>Breeder : </strong><a href="/breederbyyear/{q.breeder_id}/{q.year}/">Queens list {q.year}.</a></p>'
        f"<h1>Details</h1>"
        f'<p title="{title}"><strong>Insemination Code : </strong>{code}<strong> -- Description </strong>{desc}</p>'
        f"<p><strong>Declared Race : </strong>No<strong> -- Race Code : </strong>B"
        f"<strong> -- Lineage : </strong>Unknown<strong> -- Description </strong>buckfast</p>"
        f"<p><strong>Import : </strong>no<strong> -- Comment : </strong>No</p>{radar}</div>"
        f'<div class="list-group">{ratings}</div>')
```

- [ ] **Step 2 : Écrire les tests (en échec)**

`tests/mellifera/test_mellifera_parse_queen.py` :

```python
import datetime as dt

import pytest

from mellifera.parse.common import ParseError
from mellifera.parse.queen import parse_queen
from mellifera.parse.queens_list import parse_queens_list
from tests.mellifera.fake_site import Q, fixture, queen_page


# ------------------------------------------------------------------ /breederbyyear/
def test_queens_list():
    page = parse_queens_list(fixture("breederbyyear_145_2022.html"))
    assert (page.breeder_id, page.trigram, page.country_id) == (145, "CBP", 8)
    assert page.years == tuple(range(2024, 2012, -1))
    assert len(page.rows) == 16
    first = page.rows[0]
    assert (first.code, first.birth_year, first.mother_code, first.drone_mother_code, first.cotation_count) == (
        "FR-CBP-A-0242-2022", 2022, "FR-CBP-A-0203-2020", "FR-SPD-B-0109-2018", 1)


def test_queens_list_code_comes_from_the_link_not_the_displayed_name():
    page = parse_queens_list(fixture("breederbyyear_145_2022.html"))
    row = next(r for r in page.rows if r.name == "FR-CBP-B-0293-2022")
    assert row.code == "FR-CBP-B-0193-2022"


# ------------------------------------------------------------------ /queen/
def test_queen_with_cotation_and_radar():
    q = parse_queen(fixture("queen_FR-CBP-B-0276-2022.html"))
    assert (q.code, q.name, q.birth_year, q.breeder_id) == ("FR-CBP-B-0276-2022", "FR-CBP-B-0276-2022", 2022, 145)
    assert (q.mother_code, q.drone_mother_code, q.n_drone_queens) == ("FR-CBP-B-0006-2020", "FR-CBP-B-0943-2019", 9)
    assert q.coi == pytest.approx(0.0078125) and q.coi_reliability == pytest.approx(0.754)
    assert q.not_found_raw == "31/126"
    assert q.insemination_date == dt.date(2022, 6, 15)
    assert (q.mating_code, q.mating_description, q.mating_category_raw) == (
        "ins", "Instrumental insemination", "Instrumental insemination (ins)")
    assert (q.race_code, q.race_description, q.lineage) == ("B", "buckfast", None)
    assert q.cotation_count == 1
    assert [(s.year, s.values) for s in q.radar] == [(2023, {
        "swarming": 6.0, "mildness": 6.0, "frameHolding": 6.0, "wintering": 5.0, "broodVitality": 6.0})]


def test_queen_without_cotation():
    q = parse_queen(fixture("queen_FR-NG-B-4536-2024.html"))
    assert (q.code, q.name, q.breeder_id, q.single_drone) == ("FR-NG-B-4536-2024", "B4536(NG)", 984, False)
    # données du site conservées telles quelles : lien mère 0098, commentaire 0298 (spec §2.2)
    assert (q.mother_code, q.mother_comment) == ("FR-NG-B-0098-2022", "FR-NG-B-0298-2022")
    assert (q.drone_mother_code, q.n_drone_queens) == ("FR-CBP-B-0104-2021", 3)
    assert (q.mating_code, q.mating_category_raw) == ("inng", "Instrumental insemination (FR)")
    assert q.lineage == "PrimGreen00"
    assert q.cotation_count == 0 and q.radar == ()


def test_queen_unknown_mother_single_drone_and_multi_year_radar():
    q = parse_queen(queen_page(Q("FR-PJ-B-0562-2020", 2020, 7, name="B562(PJ)1dr", drones="FR-PJ-B-0001-2018",
                                 number=1, radar={2021: [5, 4, 5, 5, 6], 2022: [6, 5, 5, 4, 6]})))
    assert q.mother_code is None and q.single_drone is True and q.n_drone_queens == 1
    assert [s.year for s in q.radar] == [2021, 2022]
    assert q.radar[1].values["wintering"] == 4.0


def test_queen_page_without_card_is_a_parse_error():
    with pytest.raises(ParseError):
        parse_queen(b"<html><body><h1>Page not found</h1></body></html>")
```

- [ ] **Step 3 : Vérifier l'échec**

Run : `.venv/bin/pytest tests/mellifera/test_mellifera_parse_queen.py -q`
Expected : erreur de collecte `ModuleNotFoundError: No module named 'mellifera.parse.queen'`.

- [ ] **Step 4 : Écrire les deux parseurs**

`mellifera/parse/queens_list.py` :

```python
"""Page `/breederbyyear/<éleveur>/<année>/` (« Queens list ») : reines nées cette année chez l'éleveur."""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Optional

from bs4 import Tag

from .common import ParseError, make_soup, node_text, parse_int, queen_code_from_href

_BREEDER = re.compile(r"^/breeder/(\d+)/$")
_COUNTRY = re.compile(r"^/countrybyyear/(\d+)/$")
_YEAR = re.compile(r"^/breederbyyear/(\d+)/(\d{4})/$")


@dataclass(frozen=True)
class QueenListRow:
    code: str                           # code de l'URL /queen/<code>/ (le nom affiché peut différer)
    name: str
    birth_year: Optional[int]
    mother_code: Optional[str]
    drone_mother_code: Optional[str]    # colonne « Drones (RAM Mother) »
    cotation_count: Optional[int]
    serial_raw: Optional[str]


@dataclass(frozen=True)
class QueensListPage:
    breeder_id: int
    trigram: str
    country_id: Optional[int]
    years: tuple[int, ...]              # années proposées par la navigation de la page
    rows: tuple[QueenListRow, ...]


def _code(col: Tag) -> Optional[str]:
    a = col.find("a", href=True)
    return queen_code_from_href(a["href"]) if a is not None else None


def parse_queens_list(body: bytes) -> QueensListPage:
    soup = make_soup(body)
    breeder = soup.find("a", href=_BREEDER)
    if breeder is None:
        raise ParseError("éleveur absent du fil d'Ariane")
    country = soup.find("a", href=_COUNTRY)
    years = sorted({int(_YEAR.match(a["href"]).group(2)) for a in soup.find_all("a", href=_YEAR)}, reverse=True)
    rows: list[QueenListRow] = []
    for li in soup.find_all("li", class_="list-group-item"):
        if "bg-light" in (li.get("class") or []):
            continue                    # ligne d'en-tête
        row = li.find("div", class_="row")
        cols = row.find_all("div", recursive=False) if row is not None else []
        if len(cols) < 7:
            continue
        code = _code(cols[0])
        if code is None:
            continue
        serial = node_text(cols[5])
        rows.append(QueenListRow(code=code, name=node_text(cols[0]) or code, birth_year=parse_int(node_text(cols[1])),
                                 mother_code=_code(cols[2]), drone_mother_code=_code(cols[3]),
                                 cotation_count=parse_int(node_text(cols[4])),
                                 serial_raw=None if serial in ("", "None") else serial))
    return QueensListPage(breeder_id=int(_BREEDER.match(breeder["href"]).group(1)), trigram=node_text(breeder),
                          country_id=int(_COUNTRY.match(country["href"]).group(1)) if country is not None else None,
                          years=tuple(years), rows=tuple(rows))
```

`mellifera/parse/queen.py` :

```python
"""Page `/queen/<code>/` : fiche d'une reine (pedigree, accouplement, consanguinité, graphique des cotations)."""
from __future__ import annotations

import datetime as dt
import re
from dataclasses import dataclass
from typing import Optional

from .common import (LabelledField, ParseError, field_after, labelled_fields, make_soup, node_text, norm,
                     parse_django_date, parse_int, parse_number, queen_code_from_href)

_BREEDER = re.compile(r"^/breederbyyear/(\d+)/\d{4}/$")
_RATING = re.compile(r"^/rating/\d+$")
_RADAR_LABELS = re.compile(r"labels:\s*\[([^\]]*)\]")
_RADAR_SET = re.compile(r"label:\s*'(\d{4})'.*?data:\s*\[([^\]]*)\]", re.S)
_EMPTY = ("", "No", "None", "Unknown")


@dataclass(frozen=True)
class RadarSeries:
    year: int
    values: dict[str, float]            # libellé du graphique (« frameHolding ») → valeur


@dataclass(frozen=True)
class QueenPage:
    code: str
    name: str
    birth_year: Optional[int]
    breeder_id: Optional[int]
    mother_code: Optional[str]
    mother_comment: Optional[str]
    drone_mother_code: Optional[str]    # « Drones : X » : mère des reines productrices de mâles
    n_drone_queens: Optional[int]       # « Number : n »
    drone_comment: Optional[str]
    single_drone: bool                  # suffixe « 1dr » du nom
    generation_raw: Optional[str]
    coi: Optional[float]
    coi_reliability: Optional[float]    # fraction (« 70.6% » → 0.706)
    not_found_raw: Optional[str]
    serial_raw: Optional[str]
    insemination_date: Optional[dt.date]
    mating_code: Optional[str]
    mating_description: Optional[str]
    mating_category_raw: Optional[str]  # title du paragraphe : « Instrumental insemination (FR) »
    race_code: Optional[str]
    race_description: Optional[str]
    lineage: Optional[str]
    origin_raw: Optional[str]
    data_origin_raw: Optional[str]
    cotation_count: int
    radar: tuple[RadarSeries, ...]


def _value(f: Optional[LabelledField]) -> Optional[str]:
    return None if f is None or f.text in _EMPTY else f.text


def _queen_code(f: Optional[LabelledField]) -> Optional[str]:
    for h in f.hrefs if f is not None else ():
        code = queen_code_from_href(h)
        if code:
            return code
    return None


def _radar(soup) -> tuple[RadarSeries, ...]:
    for script in soup.find_all("script"):
        src = script.string or ""
        m = _RADAR_LABELS.search(src)
        if m is None or "datasets" not in src:
            continue
        labels = re.findall(r"'([^']+)'", m.group(1))
        out = []
        for year, data in _RADAR_SET.findall(src):
            values = {lab: v for lab, v in zip(labels, (parse_number(x) for x in data.split(","))) if v is not None}
            if values:
                out.append(RadarSeries(int(year), values))
        return tuple(out)
    return ()


def _mating_category(card) -> Optional[str]:
    for strong in card.find_all("strong"):
        if norm(strong.get_text(" ")).rstrip(" :") == "Insemination Code":
            p = strong.find_parent("p")
            return norm(p.get("title")) or None if p is not None else None
    return None


def parse_queen(body: bytes) -> QueenPage:
    soup = make_soup(body)
    card = soup.find("div", class_="cadree")
    if card is None:
        raise ParseError("fiche absente")
    fields = labelled_fields(card)
    code_f = field_after(fields, "International Code", "International Code")
    if code_f is None or not code_f.text:
        raise ParseError("code international absent")
    h1 = card.find("h1")
    name = node_text(h1) if h1 is not None else ""
    name = name or code_f.text
    breeder_f = field_after(fields, "Breeder", "Breeder")
    breeder_id = next((int(m.group(1)) for h in (breeder_f.hrefs if breeder_f else ()) if (m := _BREEDER.match(h))),
                      None)
    reliability = parse_number(_value(field_after(fields, "COI", "Fiability")))
    n_dq = parse_int(_value(field_after(fields, "Drones", "Number", stop="Breeder")))
    return QueenPage(
        code=code_f.text,
        name=name,
        birth_year=parse_int(_value(field_after(fields, "Birth Year", "Birth Year"))),
        breeder_id=breeder_id,
        mother_code=_queen_code(field_after(fields, "Mother", "Mother")),
        mother_comment=_value(field_after(fields, "Mother", "Comment", stop="Drones")),
        drone_mother_code=_queen_code(field_after(fields, "Drones", "Drones")),
        n_drone_queens=n_dq if n_dq is not None and n_dq >= 1 else None,
        drone_comment=_value(field_after(fields, "Drones", "Comment", stop="Breeder")),
        single_drone=name.endswith("1dr"),
        generation_raw=_value(field_after(fields, "Generation", "Generation")),
        coi=parse_number(_value(field_after(fields, "COI", "COI"))),
        coi_reliability=round(reliability / 100, 6) if reliability is not None else None,
        not_found_raw=_value(field_after(fields, "COI", "Not found")),
        serial_raw=_value(field_after(fields, "General comment", "In a serial")),
        insemination_date=parse_django_date(_value(field_after(fields, "Insemination Date", "Insemination Date"))),
        mating_code=_value(field_after(fields, "Insemination Code", "Insemination Code")),
        mating_description=_value(field_after(fields, "Insemination Code", "Description", stop="Declared Race")),
        mating_category_raw=_mating_category(card),
        race_code=_value(field_after(fields, "Declared Race", "Race Code")),
        race_description=_value(field_after(fields, "Lineage", "Description", stop="Import")),
        lineage=_value(field_after(fields, "Declared Race", "Lineage")),
        origin_raw=_value(field_after(fields, "Import", "Geographic")),
        data_origin_raw=_value(field_after(fields, "Data origin", "Data origin")),
        cotation_count=len(soup.find_all("a", href=_RATING)),
        radar=_radar(soup),
    )
```

- [ ] **Step 5 : Vérifier que les tests passent**

Run : `.venv/bin/pytest tests/mellifera/test_mellifera_parse_queen.py -q && .venv/bin/pytest -q | tail -1`
Expected : `6 passed`, puis `277 passed, 12 deselected`.

- [ ] **Step 6 : Commit**

```bash
git add mellifera/parse/queens_list.py mellifera/parse/queen.py tests/mellifera/fake_site.py tests/mellifera/test_mellifera_parse_queen.py
git commit -F - <<'EOF'
Ajoute les parseurs des listes de reines et des fiches

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015J8yR1HR7mLSRtwnvEYap6
EOF
```

---

### Task 6 : Parseur des pages de saison et robustesse de tous les parseurs

Spec §2.1, §6 (tableaux d'évaluations, identification de la reine d'une ligne : étapes 1 et 2) ; écarts 1, 2, 18 et 19.

Une page `/breeder_Year/<éleveur>/<année>/` contient un tableau par protocole, chacun sous un titre `<h5>mode : evjp</h5>`. Les colonnes sont lues dans l'en-tête de chaque tableau : deux tableaux de la même page n'ont pas les mêmes. Les colonnes descriptives sont Queen, BirthYear, Tester, YearEva, tester (id du testeur), year (date Django) et Comment (non conservée) ; toutes les autres sont des critères. Une cellule `-` ou vide est absente. Une ligne dont le nombre de cellules diffère de l'en-tête est une `ParseError` : un décalage de colonnes produirait des notes fausses.

La reine d'une ligne est résolue **sur la page**, dans cet ordre (écart 1) :
1. un lien `/queen/<code>/` unique de même texte ;
2. sinon, le texte lui-même s'il a le format d'un code et n'a aucun lien.

Les étapes 3 et 4 (recherche en base) relèvent de l'ingestion (tâche 8). La légende des codes de fécondation (`#insemCodeLegend`) donne le code, la description et la catégorie (`title`).

Les tests de robustesse vérifient enfin que les cinq parseurs lèvent `ParseError` sur une page vide, blanche ou tronquée, synthétique ou réelle.

**Files:**
- Create: `mellifera/parse/season.py`
- Modify: `tests/mellifera/fake_site.py` (ajout en fin de fichier)
- Test: `tests/mellifera/test_mellifera_parse_season.py`

**Interfaces:**
- Consumes : outils de `mellifera.parse.common` (tâche 4) ; `parse_countries` et `parse_country` (tâche 4) ; `parse_queens_list` et `parse_queen` (tâche 5) ; fixtures `breeder_Year_145_2021.html` et `breeder_Year_145_2022.html`.
- Produces :
  - `mellifera.parse.season` :
    - `FIXED_COLUMNS` ;
    - `SeasonEvaluation(queen_name_raw, queen_birth_year, queen_code, mode, tester, tester_breeder_id, eval_year, eval_date, values: dict[str, str])` ;
    - `MatingCodeEntry(code, description, category_raw)` et `SeasonPage(evaluations: tuple[SeasonEvaluation, ...], mating_codes: tuple[MatingCodeEntry, ...])` ;
    - `parse_season(body) -> SeasonPage`.
  - `fake_site` :
    - `E(queen, year, values, mode="evjp", tester="CBP", tester_id=145, date="July 6, 2023")` ;
    - `season_page(evaluations, links=(), legend=()) -> bytes`, où `links` = [(code, texte du lien)] et `legend` = [(code, description, title)].

- [ ] **Step 1 : Ajouter les pages de saison synthétiques au faux site**

Ajouter à la fin de `tests/mellifera/fake_site.py`, précédé de deux lignes vides :

```python
@dataclass
class E:
    """Ligne d'un tableau d'évaluations."""
    queen: str                          # texte de la cellule Queen
    year: int
    values: dict[str, str]
    mode: str = "evjp"
    tester: str = "CBP"
    tester_id: int = 145
    date: str = "July 6, 2023"


def season_page(evaluations: list[E], links: list[tuple[str, str]] = (),
                legend: list[tuple[str, str, str]] = ()) -> bytes:
    """links : [(code, texte du lien)] présents sur la page ; legend : [(code, description, title)]."""
    anchors = "".join(f'<a href="/queen/{c}/">{t}</a> ' for c, t in links)
    tables = ""
    for mode in dict.fromkeys(e.mode for e in evaluations):
        rows = [e for e in evaluations if e.mode == mode]
        crit = list(dict.fromkeys(k for e in rows for k in e.values))
        head = "".join(f'<td class="tdh01">{h}</td>' for h in
                       ["Queen", "BirthYear", "Tester", "YearEva", "tester", "year", *crit, "Comment"])
        body = "".join(
            "<tr>" + "".join(f'<td class="td01">{v}</td>' for v in
                             [e.queen, e.year, e.tester, e.date[-4:], e.tester_id, e.date,
                              *(e.values.get(c, "-") for c in crit), ""]) + "</tr>"
            for e in rows)
        tables += f"<h5>mode : {mode}</h5><table><tr>{head}</tr>{body}</table>"
    leg = "".join(f'<div class="row"><div class="col-2" title="{t}">{c}</div><div class="col-10">{d}</div></div>'
                  for c, d, t in legend)
    return _page(f'{anchors}<h3>Evaluations</h3>{tables}<div id="insemCodeLegend">{leg}</div>')
```

- [ ] **Step 2 : Écrire les tests (en échec)**

`tests/mellifera/test_mellifera_parse_season.py` :

```python
import datetime as dt

import pytest

from mellifera.parse.common import ParseError
from mellifera.parse.countries import parse_countries
from mellifera.parse.country import parse_country
from mellifera.parse.queen import parse_queen
from mellifera.parse.queens_list import parse_queens_list
from mellifera.parse.season import parse_season
from tests.mellifera.fake_site import E, fixture, season_page

PARSERS = [parse_countries, parse_country, parse_queens_list, parse_queen, parse_season]


# ------------------------------------------------------------------ /breeder_Year/
def test_season_two_tables_with_different_columns():
    page = parse_season(fixture("breeder_Year_145_2022.html"))
    assert {e.mode for e in page.evaluations} == {"evjp", "ceta"}
    assert all(e.queen_code is not None for e in page.evaluations)
    ev = next(e for e in page.evaluations if e.mode == "evjp" and e.queen_code == "FR-CBP-A-0203-2020")
    assert (ev.tester, ev.tester_breeder_id, ev.eval_year, ev.eval_date) == ("CBP", 145, 2023, dt.date(2023, 7, 6))
    assert ev.values["mildness"] == "6" and "hyg24h" not in ev.values      # « - » = absent
    ceta = next(e for e in page.evaluations if e.mode == "ceta" and e.queen_code == "FR-CBP-A-0087-2020")
    assert (ceta.values["honeyAmt"], ceta.values["hyg24h"]) == ("48", "100")
    assert {m.code: m.category_raw for m in page.mating_codes}["mdm"] == "Mating station (FR)"


def test_season_names_resolved_by_page_links():
    page = parse_season(fixture("breeder_Year_145_2021.html"))
    ev = page.evaluations[0]
    assert (ev.queen_name_raw, ev.queen_code, ev.mode) == ("D992(CBP)", "FR-CBP-D-0992-2019", "evjp")
    assert {e.mode for e in page.evaluations} == {"evjp", "evdk"}


def test_season_link_wins_over_code_looking_name():
    body = season_page([E("FR-CBP-B-0293-2022", 2022, {"mildness": "5"})],
                       links=[("FR-CBP-B-0193-2022", "FR-CBP-B-0293-2022")])
    assert parse_season(body).evaluations[0].queen_code == "FR-CBP-B-0193-2022"


def test_season_ambiguous_or_unlinked_name_left_unresolved():
    body = season_page([E("A087", 2020, {"mildness": "5"}), E("FR-CBP-A-0041-2020", 2020, {"mildness": "4"})],
                       links=[("FR-CBP-A-0087-2020", "A087"), ("FR-XX-A-0087-2020", "A087")])
    evs = parse_season(body).evaluations
    assert evs[0].queen_code is None                    # deux liens différents pour « A087 »
    assert evs[1].queen_code == "FR-CBP-A-0041-2020"    # pas de lien : texte au format code


def test_season_row_with_missing_cells_is_a_parse_error():
    body = season_page([E("A087", 2020, {"mildness": "5"})]).replace(b'<td class="td01">5</td>', b"")
    with pytest.raises(ParseError):
        parse_season(body)


# ------------------------------------------------------------------ robustesse de tous les parseurs
@pytest.mark.parametrize("parse", PARSERS)
@pytest.mark.parametrize("body", [b"", b"   ", b"<html><body><div class='cadree'><h1>FR-"])
def test_empty_or_truncated_page_is_a_parse_error(parse, body):
    with pytest.raises(ParseError):
        parse(body)


@pytest.mark.parametrize("parse, name", [
    (parse_countries, "all.html"), (parse_country, "countrybyyear_8.html"),
    (parse_queens_list, "breederbyyear_145_2022.html"), (parse_queen, "queen_FR-CBP-B-0276-2022.html"),
    (parse_season, "breeder_Year_145_2022.html"),
])
def test_truncated_fixture_is_a_parse_error(parse, name):
    body = fixture(name)
    with pytest.raises(ParseError):
        parse(body[: len(body) // 2])
```

- [ ] **Step 3 : Vérifier l'échec**

Run : `.venv/bin/pytest tests/mellifera/test_mellifera_parse_season.py -q`
Expected : erreur de collecte `ModuleNotFoundError: No module named 'mellifera.parse.season'`.

- [ ] **Step 4 : Écrire le parseur**

`mellifera/parse/season.py` :

```python
"""Page `/breeder_Year/<éleveur>/<année>/` : tableaux d'évaluations et légende des codes de fécondation."""
from __future__ import annotations

import datetime as dt
from collections import defaultdict
from dataclasses import dataclass
from typing import Optional

from .common import ParseError, is_code, make_soup, node_text, norm, parse_django_date, parse_int, queen_code_from_href

# Colonnes descriptives ; toutes les autres sont des critères. « Comment » n'est pas conservée.
FIXED_COLUMNS = frozenset({"Queen", "BirthYear", "Tester", "YearEva", "tester", "year", "Comment"})


@dataclass(frozen=True)
class SeasonEvaluation:
    queen_name_raw: str
    queen_birth_year: Optional[int]
    queen_code: Optional[str]           # résolu sur la page (lien de même texte, sinon texte au format code)
    mode: str                           # protocole du tableau : « evjp », « ceta », « evdk »…
    tester: str                         # trigramme du testeur
    tester_breeder_id: Optional[int]
    eval_year: Optional[int]
    eval_date: Optional[dt.date]
    values: dict[str, str]              # critère → valeur brute ; cellules « - » et vides exclues


@dataclass(frozen=True)
class MatingCodeEntry:
    code: str
    description: str
    category_raw: Optional[str]         # title de la légende : « Mating station (FR) »


@dataclass(frozen=True)
class SeasonPage:
    evaluations: tuple[SeasonEvaluation, ...]
    mating_codes: tuple[MatingCodeEntry, ...]


def _link_codes(soup) -> dict[str, set[str]]:
    found: dict[str, set[str]] = defaultdict(set)
    for a in soup.find_all("a", href=True):
        code = queen_code_from_href(a["href"])
        if code:
            found[node_text(a)].add(code)
    return found


def _resolve(name: str, links: dict[str, set[str]]) -> Optional[str]:
    codes = links.get(name, set())
    if len(codes) == 1:
        return next(iter(codes))        # 1. lien de la page au texte identique (prime sur le format)
    return name if is_code(name) and not codes else None   # 2. texte au format code


def parse_season(body: bytes) -> SeasonPage:
    soup = make_soup(body)
    links = _link_codes(soup)
    evaluations: list[SeasonEvaluation] = []
    for h5 in soup.find_all("h5"):
        title = node_text(h5)
        if not title.startswith("mode") or ":" not in title:
            continue
        mode = norm(title.split(":", 1)[1])
        table = h5.find_next("table")
        rows = table.find_all("tr") if table is not None else []
        if not rows:
            continue
        header = [node_text(c) for c in rows[0].find_all(["td", "th"])]
        if "Queen" not in header:
            raise ParseError(f"tableau « {mode} » sans colonne Queen")
        for tr in rows[1:]:
            cells = [node_text(c) for c in tr.find_all(["td", "th"])]
            if len(cells) != len(header):
                raise ParseError(f"tableau « {mode} » : {len(cells)} cellules pour {len(header)} colonnes")
            row = dict(zip(header, cells))
            name = row["Queen"]
            if not name:
                continue
            evaluations.append(SeasonEvaluation(
                queen_name_raw=name,
                queen_birth_year=parse_int(row.get("BirthYear")),
                queen_code=_resolve(name, links),
                mode=mode,
                tester=row.get("Tester", ""),
                tester_breeder_id=parse_int(row.get("tester")),
                eval_year=parse_int(row.get("YearEva")),
                eval_date=parse_django_date(row.get("year")),
                values={k: v for k, v in row.items() if k not in FIXED_COLUMNS and v not in ("", "-")},
            ))
    legend = soup.find(id="insemCodeLegend")
    codes: list[MatingCodeEntry] = []
    for row in legend.find_all("div", class_="row") if legend is not None else []:
        cols = row.find_all("div", recursive=False)
        if len(cols) >= 2 and node_text(cols[0]):
            codes.append(MatingCodeEntry(node_text(cols[0]), node_text(cols[1]), norm(cols[0].get("title")) or None))
    return SeasonPage(tuple(evaluations), tuple(codes))
```

- [ ] **Step 5 : Vérifier que les tests passent**

Run : `.venv/bin/pytest tests/mellifera/test_mellifera_parse_season.py -q && .venv/bin/pytest -q | tail -1`
Expected : `25 passed`, puis `302 passed, 12 deselected`.

- [ ] **Step 6 : Commit**

```bash
git add mellifera/parse/season.py tests/mellifera/fake_site.py tests/mellifera/test_mellifera_parse_season.py
git commit -F - <<'EOF'
Ajoute le parseur des pages de saison

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015J8yR1HR7mLSRtwnvEYap6
EOF
```

---

### Task 7 : Schéma SQLite et référentiel des critères

Spec §5 ; écarts 2 à 6, 11 et 18.

La base est en mode WAL avec `busy_timeout`, si bien que `status` et `export-blup` la lisent en lecture seule (`mode=ro`) pendant qu'un crawl écrit. `open_db` crée le dossier et le schéma (`CREATE … IF NOT EXISTS`, rejouable), puis recharge le référentiel. `traits.csv` est mis à jour en place, par id. Les liens de `criterion` vers les traits sont remis à zéro puis rechargés depuis `criteria.csv` : les critères inconnus créés au parsing restent, sans trait. Une ligne de `criteria.csv` qui cite un trait inconnu lève `ValueError`.

`now_iso()` donne l'heure UTC à la microseconde (écart 11). `mother_code` et `drone_mother_code` sont des codes sans clé étrangère, car un parent peut ne pas être en base.

**Files:**
- Create: `mellifera/db.py`, `mellifera/traits.csv`, `mellifera/criteria.csv`
- Test: `tests/mellifera/test_mellifera_db.py`

**Interfaces:**
- Consumes : rien.
- Produces :
  - `mellifera.db` :
    - `SCHEMA`, `TRAITS_CSV` et `CRITERIA_CSV` ;
    - `now_iso() -> str` ;
    - `connect(path) -> sqlite3.Connection` (`row_factory = sqlite3.Row`, WAL) ;
    - `connect_readonly(path) -> sqlite3.Connection`, qui lève `FileNotFoundError` (« base introuvable : … ») si le fichier manque ;
    - `load_reference(conn)` et `open_db(path) -> sqlite3.Connection`.
  - Tables :
    - `country(id, name, iso)` ;
    - `breeder(id, trigram, country_id, first_seen, last_seen)` ;
    - `queen(id AUTOINCREMENT, code UNIQUE, code_valid, name, birth_year, breeder_id, mother_code, drone_mother_code, n_drone_queens, mating_code, single_drone, insemination_date, race_code, race_description, lineage, generation_raw, coi, coi_reliability, not_found_raw, cotation_count, list_cotation_count, serial_raw, origin_raw, data_origin_raw, mother_comment_raw, drone_comment_raw, parse_version, fetched_at, parsed_at)` ;
    - `mating_code(code, description, category_raw, mating_type, mating_type_override)` ;
    - `trait(id, name, label, kind, scale_min, scale_max, direction)` ;
    - `criterion(id, mode, code, trait_id)` ;
    - `evaluation(id, queen_name_raw, queen_birth_year_raw, queen_code, queen_key, page_breeder_id, tester, tester_breeder_id, eval_date, eval_year, mode, criterion_code, value, value_raw, source, url)`, unique sur `(queen_key, source, mode, tester, eval_date, criterion_code)` ;
    - `page(url, kind, status, http_status, fetched_at, sha256, cache_path, parse_error, parse_version, run_id)` ;
    - `crawl_run(id, mode, started_at, finished_at, stats)`.
  - Vue `queen_mating` : toutes les colonnes de `queen`, plus `mating_type`.

- [ ] **Step 1 : Écrire les tests (en échec)**

`tests/mellifera/test_mellifera_db.py` :

```python
import re
import sqlite3

import pytest

from mellifera.db import connect_readonly, now_iso, open_db


@pytest.fixture
def conn(tmp_path):
    c = open_db(tmp_path / "m.sqlite")
    yield c
    c.close()


def count(conn, sql: str, *args) -> int:
    return conn.execute(f"SELECT COUNT(*) FROM {sql}", args).fetchone()[0]


# ------------------------------------------------------------------ schéma et référentiel
def test_schema_wal_and_reference_tables(conn):
    assert conn.execute("PRAGMA journal_mode").fetchone()[0] == "wal"
    assert count(conn, "trait") == 22
    assert count(conn, "criterion WHERE trait_id IS NULL") == 0
    mapped = dict(conn.execute("SELECT c.mode || '.' || c.code, t.name FROM criterion c JOIN trait t "
                               "ON t.id = c.trait_id").fetchall())
    assert mapped["radar.frameHolding"] == mapped["evjp.framHold"] == mapped["ceta.framHold"] == "framHold"


def test_reopening_keeps_unknown_criteria_and_reloads_mapping(tmp_path):
    c = open_db(tmp_path / "m.sqlite")
    c.execute("INSERT INTO criterion (mode, code) VALUES ('evjp', 'smr')")
    c.commit()
    c.close()
    c = open_db(tmp_path / "m.sqlite")
    assert count(c, "criterion WHERE mode = 'evjp' AND code = 'smr' AND trait_id IS NULL") == 1
    assert count(c, "criterion WHERE trait_id IS NOT NULL") == 48
    c.close()


def test_readonly_connection_for_status_and_export(tmp_path):
    with pytest.raises(FileNotFoundError, match="base introuvable"):
        connect_readonly(tmp_path / "absente.sqlite")
    open_db(tmp_path / "m.sqlite").close()
    ro = connect_readonly(tmp_path / "m.sqlite")
    assert ro.execute("SELECT COUNT(*) FROM trait").fetchone()[0] == 22
    with pytest.raises(sqlite3.OperationalError):
        ro.execute("DELETE FROM trait")
    ro.close()


def test_timestamps_are_utc_to_the_microsecond():
    # « éleveur vu pour la première fois pendant ce passage » se décide en comparant à crawl_run.started_at
    assert re.fullmatch(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d\.\d{6}\+00:00", now_iso())
```

- [ ] **Step 2 : Vérifier l'échec**

Run : `.venv/bin/pytest tests/mellifera/test_mellifera_db.py -q`
Expected : erreur de collecte `ModuleNotFoundError: No module named 'mellifera.db'`.

- [ ] **Step 3 : Écrire le référentiel**

`mellifera/traits.csv` :

```text
id,name,label,kind,scale_min,scale_max,direction
1,swarming,Essaimage,note,0,6,1
2,mildness,Douceur,note,0,6,1
3,framHold,Tenue au cadre,note,0,6,1
4,winter,Hivernage,note,0,6,1
5,broVital,Vitalité du couvain,note,0,6,1
6,beeVital,Vitalité des abeilles,note,0,6,1
7,fertil,Fécondité,note,0,6,1
8,sprCrop,Récolte de printemps (note),note,0,6,1
9,midCrop,Récolte de mi-saison (note),note,0,6,1
10,sumCrop,Récolte d'été (note),note,0,6,1
11,falCrop,Récolte d'automne (note),note,0,6,1
12,polCrop,Récolte de pollen (note),note,0,6,1
13,bldCombs,Bâtisse,note,0,6,1
14,propolis,Propolis,note,0,6,1
15,honeyIdx,Indice miel,note,0,6,1
16,aveCropIdx,Indice de récolte moyen,note,0,6,1
17,hygBehaIdx,Indice de comportement hygiénique,note,0,6,1
18,varIdx,Indice varroa,note,0,6,1
19,hyg24h,Couvain tué nettoyé en 24 h (%),percentage,,,1
20,hyg48h,Couvain tué nettoyé en 48 h (%),percentage,,,1
21,honeyAmt,Quantité de miel (kg),harvest,,,1
22,aveHonAmt,Quantité de miel moyenne (kg),harvest,,,1
```

`mellifera/criteria.csv` :

```text
mode,code,trait
evjp,beeVital,beeVital
evjp,broVital,broVital
evjp,mildness,mildness
evjp,framHold,framHold
evjp,swarming,swarming
evjp,fertil,fertil
evjp,sprCrop,sprCrop
evjp,sumCrop,sumCrop
evjp,falCrop,falCrop
evjp,polCrop,polCrop
evjp,winter,winter
evjp,bldCombs,bldCombs
evjp,propolis,propolis
evjp,hyg24h,hyg24h
ceta,swarming,swarming
ceta,mildness,mildness
ceta,framHold,framHold
ceta,fertil,fertil
ceta,winter,winter
ceta,bldCombs,bldCombs
ceta,propolis,propolis
ceta,honeyIdx,honeyIdx
ceta,honeyAmt,honeyAmt
ceta,aveHonAmt,aveHonAmt
ceta,aveCropIdx,aveCropIdx
ceta,sprCrop,sprCrop
ceta,midCrop,midCrop
ceta,sumCrop,sumCrop
ceta,falCrop,falCrop
ceta,polCrop,polCrop
ceta,broVital,broVital
ceta,beeVital,beeVital
ceta,varIdx,varIdx
ceta,hygBehaIdx,hygBehaIdx
ceta,hyg24h,hyg24h
ceta,hyg48h,hyg48h
evdk,swarming,swarming
evdk,mildness,mildness
evdk,framHold,framHold
evdk,honeyIdx,honeyIdx
evdk,honeyAmt,honeyAmt
evdk,hygBehaIdx,hygBehaIdx
evdk,hyg24h,hyg24h
radar,swarming,swarming
radar,mildness,mildness
radar,frameHolding,framHold
radar,wintering,winter
radar,broodVitality,broVital
```

- [ ] **Step 4 : Écrire le schéma et les connexions**

`mellifera/db.py` :

```python
"""Base SQLite (mode WAL) : schéma, référentiel des critères, connexions."""
from __future__ import annotations

import csv
import datetime as dt
import sqlite3
from pathlib import Path

HERE = Path(__file__).parent
TRAITS_CSV = HERE / "traits.csv"
CRITERIA_CSV = HERE / "criteria.csv"

SCHEMA = """
CREATE TABLE IF NOT EXISTS country (
    id INTEGER PRIMARY KEY,
    name TEXT NOT NULL,
    iso TEXT
);
CREATE TABLE IF NOT EXISTS breeder (
    id INTEGER PRIMARY KEY,
    trigram TEXT NOT NULL,
    country_id INTEGER,
    first_seen TEXT NOT NULL,
    last_seen TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS breeder_trigram ON breeder (trigram);
CREATE TABLE IF NOT EXISTS queen (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    code TEXT NOT NULL UNIQUE,
    code_valid INTEGER NOT NULL,
    name TEXT,
    birth_year INTEGER,
    breeder_id INTEGER,
    mother_code TEXT,
    drone_mother_code TEXT,
    n_drone_queens INTEGER,
    mating_code TEXT,
    single_drone INTEGER NOT NULL DEFAULT 0,
    insemination_date TEXT,
    race_code TEXT,
    race_description TEXT,
    lineage TEXT,
    generation_raw TEXT,
    coi REAL,
    coi_reliability REAL,
    not_found_raw TEXT,
    cotation_count INTEGER,
    list_cotation_count INTEGER,
    serial_raw TEXT,
    origin_raw TEXT,
    data_origin_raw TEXT,
    mother_comment_raw TEXT,
    drone_comment_raw TEXT,
    parse_version INTEGER,
    fetched_at TEXT,
    parsed_at TEXT
);
CREATE INDEX IF NOT EXISTS queen_name ON queen (name, birth_year);
CREATE TABLE IF NOT EXISTS mating_code (
    code TEXT PRIMARY KEY,
    description TEXT,
    category_raw TEXT,
    mating_type TEXT NOT NULL CHECK (mating_type IN ('natural', 'station', 'insemination', 'unknown')),
    mating_type_override TEXT CHECK (mating_type_override IN ('natural', 'station', 'insemination', 'unknown'))
);
CREATE VIEW IF NOT EXISTS queen_mating AS
    SELECT q.*, COALESCE(m.mating_type_override, m.mating_type, 'unknown') AS mating_type
    FROM queen q LEFT JOIN mating_code m ON m.code = q.mating_code;
CREATE TABLE IF NOT EXISTS trait (
    id INTEGER PRIMARY KEY,
    name TEXT NOT NULL UNIQUE,
    label TEXT NOT NULL,
    kind TEXT NOT NULL CHECK (kind IN ('note', 'percentage', 'harvest')),
    scale_min REAL,
    scale_max REAL,
    direction INTEGER NOT NULL CHECK (direction IN (1, -1)),
    CHECK (kind != 'note' OR scale_max > scale_min)
);
CREATE TABLE IF NOT EXISTS criterion (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    mode TEXT NOT NULL,
    code TEXT NOT NULL,
    trait_id INTEGER REFERENCES trait (id),
    UNIQUE (mode, code)
);
CREATE TABLE IF NOT EXISTS evaluation (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    queen_name_raw TEXT NOT NULL,
    queen_birth_year_raw INTEGER,
    queen_code TEXT,
    queen_key TEXT NOT NULL,
    page_breeder_id INTEGER,
    tester TEXT NOT NULL,
    tester_breeder_id INTEGER,
    eval_date TEXT NOT NULL,
    eval_year INTEGER,
    mode TEXT NOT NULL,
    criterion_code TEXT NOT NULL,
    value REAL,
    value_raw TEXT NOT NULL,
    source TEXT NOT NULL CHECK (source IN ('season_table', 'radar')),
    url TEXT NOT NULL,
    UNIQUE (queen_key, source, mode, tester, eval_date, criterion_code)
);
CREATE INDEX IF NOT EXISTS evaluation_queen ON evaluation (queen_code);
CREATE TABLE IF NOT EXISTS page (
    url TEXT PRIMARY KEY,
    kind TEXT NOT NULL,
    status TEXT NOT NULL CHECK (status IN ('ok', 'missing', 'forbidden', 'error')),
    http_status INTEGER,
    fetched_at TEXT NOT NULL,
    sha256 TEXT,
    cache_path TEXT,
    parse_error TEXT,
    parse_version INTEGER,
    run_id INTEGER
);
CREATE TABLE IF NOT EXISTS crawl_run (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    mode TEXT NOT NULL CHECK (mode IN ('full', 'incremental', 'reparse')),
    started_at TEXT NOT NULL,
    finished_at TEXT,
    stats TEXT
);
"""


def now_iso() -> str:
    """Horodatage UTC à la microseconde : « vu pour la première fois pendant ce passage » se décide par
    comparaison avec crawl_run.started_at, sans égalité possible entre deux passages."""
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="microseconds")


def connect(path: Path) -> sqlite3.Connection:
    conn = sqlite3.connect(path)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA busy_timeout=5000")
    return conn


def connect_readonly(path: Path) -> sqlite3.Connection:
    """Lecture seule (export, moteur) : n'écrit rien, même pendant un crawl (WAL)."""
    if not Path(path).is_file():
        raise FileNotFoundError(f"base introuvable : {path} (lancez d'abord « python -m mellifera crawl --full »)")
    conn = sqlite3.connect(f"file:{Path(path).resolve()}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    return conn


def _opt_float(s: str):
    return float(s) if s.strip() else None


def load_reference(conn: sqlite3.Connection) -> None:
    """Recharge traits.csv et criteria.csv (sources de vérité) ; les critères inconnus restent sans trait."""
    with TRAITS_CSV.open(encoding="utf-8") as f:
        for r in csv.DictReader(f):
            conn.execute(
                "INSERT INTO trait (id, name, label, kind, scale_min, scale_max, direction) VALUES (?, ?, ?, ?, ?, ?, ?) "
                "ON CONFLICT (id) DO UPDATE SET name = excluded.name, label = excluded.label, kind = excluded.kind, "
                "scale_min = excluded.scale_min, scale_max = excluded.scale_max, direction = excluded.direction",
                (int(r["id"]), r["name"], r["label"], r["kind"], _opt_float(r["scale_min"]),
                 _opt_float(r["scale_max"]), int(r["direction"])))
    conn.execute("UPDATE criterion SET trait_id = NULL")
    with CRITERIA_CSV.open(encoding="utf-8") as f:
        for r in csv.DictReader(f):
            trait = conn.execute("SELECT id FROM trait WHERE name = ?", (r["trait"],)).fetchone()
            if trait is None:
                raise ValueError(f"criteria.csv : trait inconnu « {r['trait']} » ({r['mode']}, {r['code']})")
            conn.execute("INSERT INTO criterion (mode, code, trait_id) VALUES (?, ?, ?) "
                         "ON CONFLICT (mode, code) DO UPDATE SET trait_id = excluded.trait_id",
                         (r["mode"], r["code"], trait["id"]))
    conn.commit()


def open_db(path: Path) -> sqlite3.Connection:
    """Base en écriture (crawler) : crée le dossier, le schéma et recharge le référentiel."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    conn = connect(path)
    conn.executescript(SCHEMA)
    load_reference(conn)
    return conn
```

- [ ] **Step 5 : Vérifier que les tests passent**

Run : `.venv/bin/pytest tests/mellifera/test_mellifera_db.py -q && .venv/bin/pytest -q | tail -1`
Expected : `4 passed`, puis `306 passed, 12 deselected`.

- [ ] **Step 6 : Commit**

```bash
git add mellifera/db.py mellifera/traits.csv mellifera/criteria.csv tests/mellifera/test_mellifera_db.py
git commit -F - <<'EOF'
Ajoute le schéma SQLite et le référentiel des critères

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015J8yR1HR7mLSRtwnvEYap6
EOF
```

---

### Task 8 : Ingestion des pages et identification des reines évaluées

Spec §5 (« écritures idempotentes ») et §6 (identification, étapes 3 à 5, et reprise de la résolution en fin de crawl) ; écarts 1, 3, 6, 7 et 10.

`ingest(conn, path, body, seen_at)` choisit le parseur d'après le chemin et écrit le résultat, sans commit (l'appelant valide). Toutes les écritures sont des upserts, si bien qu'ingérer deux fois les mêmes pages donne la même base (`db_snapshot`).

- Une liste crée l'éleveur et les ébauches de reines. La fiche remplace l'ébauche et ajoute ses valeurs du graphique (`source = 'radar'`, date au 1er janvier de l'année du jeu).
- Une ligne de tableau dont le code n'a pas été résolu sur la page est cherchée en base :
  1. par (nom, année, éleveur), l'éleveur étant le trigramme entre parenthèses s'il y en a un, sinon l'éleveur de la page ;
  2. sinon par (nom, année), si le résultat est unique ;
  3. sinon, la ligne reste non résolue, sous la clé `name:<nom>|<année>|<éleveur de la page>`.
- `resolve_pending` rejoue ces deux recherches. Une ligne qui se résout prend son code pour clé : `UPDATE OR IGNORE`, puis suppression du reste, qui doublonnait une ligne déjà enregistrée sous ce code.

**Files:**
- Create: `mellifera/ingest.py`
- Modify: `tests/mellifera/fake_site.py` (ajout en fin de fichier)
- Test: `tests/mellifera/test_mellifera_ingest.py`

**Interfaces:**
- Consumes :
  - `now_iso` et `open_db` (tâche 7) ;
  - `PARSER_VERSION`, `is_code` et `parse_number` (tâche 4) ;
  - les cinq parseurs et leurs dataclasses (tâches 4 à 6) ;
  - `kind_of` et `path_ids` (tâche 1).
- Produces :
  - `mellifera.ingest` :
    - `mating_type_for(category_raw, description) -> str` ;
    - `store_countries(conn, entries)` et `store_country(conn, country_id, breeders, seen_at)` ;
    - `store_queens_list(conn, page, seen_at)` et `store_queen(conn, q, path, seen_at)` ;
    - `lookup_queen(conn, name, birth_year, page_breeder_id) -> Optional[str]` et `unresolved_key(name, birth_year, page_breeder_id) -> str` ;
    - `store_season(conn, page, path, page_breeder_id)` et `resolve_pending(conn) -> int` (nombre de reines résolues) ;
    - `ingest(conn, path, body, seen_at)`, qui rend la structure parsée et laisse passer `ParseError`.
  - `fake_site` : `db_snapshot(conn) -> dict[str, list[tuple]]`, le contenu comparable de la base, sans id d'évaluation ni dates de parsing.

- [ ] **Step 1 : Ajouter l'instantané de base au faux site**

Ajouter à la fin de `tests/mellifera/fake_site.py`, précédé de deux lignes vides :

```python
# ------------------------------------------------------------------ base
def db_snapshot(conn) -> dict[str, list[tuple]]:
    """Contenu comparable de la base : sans id d'évaluation ni dates de parsing."""
    tables = {
        "country": "SELECT * FROM country ORDER BY id",
        "breeder": "SELECT * FROM breeder ORDER BY id",
        "queen": "SELECT * FROM queen ORDER BY code",
        "mating_code": "SELECT * FROM mating_code ORDER BY code",
        "criterion": "SELECT mode, code, trait_id FROM criterion ORDER BY mode, code",
        "evaluation": "SELECT * FROM evaluation ORDER BY queen_key, source, mode, tester, eval_date, criterion_code",
    }
    out = {}
    for name, sql in tables.items():
        rows = conn.execute(sql).fetchall()
        drop = {"parsed_at", "id"} if name in ("queen", "evaluation") else set()
        out[name] = [tuple(v for k, v in zip(r.keys(), r) if k not in drop) for r in rows]
    return out
```

- [ ] **Step 2 : Écrire les tests (en échec)**

`tests/mellifera/test_mellifera_ingest.py` :

```python
import pytest

from mellifera.db import open_db
from mellifera.ingest import ingest, mating_type_for, resolve_pending
from tests.mellifera.fake_site import E, Q, db_snapshot, fixture, queen_page, queens_list_page, season_page

SEEN = "2026-10-13T10:00:00+00:00"


@pytest.fixture
def conn(tmp_path):
    c = open_db(tmp_path / "m.sqlite")
    yield c
    c.close()


def load(conn, pages: dict[str, bytes]) -> None:
    for path, body in pages.items():
        ingest(conn, path, body, SEEN)
    conn.commit()


def count(conn, sql: str, *args) -> int:
    return conn.execute(f"SELECT COUNT(*) FROM {sql}", args).fetchone()[0]


FIXTURE_PAGES = {
    "/all/": fixture("all.html"),
    "/countrybyyear/8/": fixture("countrybyyear_8.html"),
    "/breederbyyear/145/2022/": fixture("breederbyyear_145_2022.html"),
    "/queen/FR-CBP-B-0276-2022/": fixture("queen_FR-CBP-B-0276-2022.html"),
    "/queen/FR-NG-B-4536-2024/": fixture("queen_FR-NG-B-4536-2024.html"),
    "/breeder_Year/145/2021/": fixture("breeder_Year_145_2021.html"),
    "/breeder_Year/145/2022/": fixture("breeder_Year_145_2022.html"),
    "/breeder_Year/984/2026/": fixture("breeder_Year_984_2026.html"),
}


# ------------------------------------------------------------------ type d'accouplement
@pytest.mark.parametrize("category, description, expected", [
    ("Instrumental insemination (FR)", "Thorigné (FR) instrumental insemination", "insemination"),
    ("Mating station (FR)", "Madame island (FR) mating station", "station"),
    (None, "Belegstelle Norderney", "station"),
    (None, "Instrumental insemination", "insemination"),
    ("Free mating", "Open mating", "unknown"),
    (None, None, "unknown"),
])
def test_mating_type_for(category, description, expected):
    assert mating_type_for(category, description) == expected


# ------------------------------------------------------------------ ingestion des fixtures
def test_ingest_fixtures(conn):
    load(conn, FIXTURE_PAGES)
    assert count(conn, "country") == 32
    assert count(conn, "breeder") == 85
    q = conn.execute("SELECT * FROM queen_mating WHERE code = 'FR-CBP-B-0276-2022'").fetchone()
    assert (q["mother_code"], q["drone_mother_code"], q["n_drone_queens"], q["mating_type"]) == (
        "FR-CBP-B-0006-2020", "FR-CBP-B-0943-2019", 9, "insemination")
    stub = conn.execute("SELECT * FROM queen WHERE code = 'FR-CBP-A-0242-2022'").fetchone()
    assert (stub["mother_code"], stub["parsed_at"]) == ("FR-CBP-A-0203-2020", None)   # ébauche issue de la liste
    radar = conn.execute("SELECT eval_date, value FROM evaluation WHERE source = 'radar' "
                         "AND criterion_code = 'wintering'").fetchall()
    assert [tuple(r) for r in radar] == [("2023-01-01", 5.0)]
    assert count(conn, "evaluation WHERE queen_code IS NULL") == 0
    assert conn.execute("SELECT mating_type FROM mating_code WHERE code = 'mdm'").fetchone()[0] == "station"


def test_ingest_is_idempotent(conn):
    load(conn, FIXTURE_PAGES)
    first = db_snapshot(conn)
    load(conn, FIXTURE_PAGES)
    assert db_snapshot(conn) == first


def test_same_evaluation_on_two_season_pages_is_stored_once(conn):
    rows = [E("FR-CBP-B-0104-2021", 2021, {"mildness": "5", "swarming": "6"})]
    load(conn, {"/breeder_Year/145/2022/": season_page(rows), "/breeder_Year/984/2026/": season_page(rows)})
    assert count(conn, "evaluation") == 2


def test_unknown_mode_and_column_create_unmapped_criteria(conn):
    load(conn, {"/breeder_Year/145/2022/": season_page([E("FR-CBP-B-0104-2021", 2021, {"smr": "40", "mildness": "5"},
                                                         mode="xyz")])})
    assert count(conn, "criterion WHERE mode = 'xyz' AND trait_id IS NULL") == 2
    assert count(conn, "evaluation WHERE mode = 'xyz'") == 2


def test_non_numeric_value_kept_raw(conn):
    load(conn, {"/breeder_Year/145/2022/": season_page([E("FR-CBP-B-0104-2021", 2021, {"mildness": "bof"})])})
    r = conn.execute("SELECT value, value_raw FROM evaluation").fetchone()
    assert (r["value"], r["value_raw"]) == (None, "bof")


# ------------------------------------------------------------------ identification des reines évaluées (spec §6)
def _list(*queens: Q) -> bytes:
    return queens_list_page(145, "CBP", 8, [2020], list(queens))


def test_row_by_code(conn):
    load(conn, {"/breeder_Year/145/2022/": season_page([E("FR-CBP-A-0041-2020", 2020, {"mildness": "4"})])})
    assert conn.execute("SELECT queen_code, queen_key FROM evaluation").fetchone()[:] == (
        "FR-CBP-A-0041-2020", "FR-CBP-A-0041-2020")


def test_row_by_name_with_trigram_resolved_in_db(conn):
    load(conn, {"/queen/FR-SPD-B-0131-2019/": queen_page(Q("FR-SPD-B-0131-2019", 2019, 12, name="B131(SPD)")),
                "/breederbyyear/12/2019/": queens_list_page(12, "SPD", 8, [2019], [])})
    load(conn, {"/breeder_Year/145/2022/": season_page([E("B131(SPD)", 2019, {"mildness": "4"})])})
    assert conn.execute("SELECT queen_code FROM evaluation").fetchone()[0] == "FR-SPD-B-0131-2019"


def test_row_by_bare_name_uses_page_breeder(conn):
    load(conn, {"/breederbyyear/145/2020/": _list(Q("FR-CBP-A-0087-2020", 2020, 145, name="A087")),
                "/breederbyyear/12/2020/": queens_list_page(12, "SPD", 8, [2020],
                                                            [Q("FR-SPD-A-0087-2020", 2020, 12, name="A087")])})
    load(conn, {"/breeder_Year/145/2022/": season_page([E("A087", 2020, {"mildness": "4"})])})
    assert conn.execute("SELECT queen_code FROM evaluation").fetchone()[0] == "FR-CBP-A-0087-2020"


def test_ambiguous_bare_name_left_unresolved(conn):
    load(conn, {"/breederbyyear/11/2020/": queens_list_page(11, "AAA", 8, [2020], [Q("FR-AAA-A-0087-2020", 2020, 11,
                                                                                       name="A087")]),
                "/breederbyyear/12/2020/": queens_list_page(12, "SPD", 8, [2020],
                                                            [Q("FR-SPD-A-0087-2020", 2020, 12, name="A087")])})
    load(conn, {"/breeder_Year/145/2022/": season_page([E("A087", 2020, {"mildness": "4"})])})
    r = conn.execute("SELECT queen_code, queen_key FROM evaluation").fetchone()
    assert (r["queen_code"], r["queen_key"]) == (None, "name:A087|2020|145")


def test_unresolved_row_resolved_after_queen_arrives_and_merged(conn):
    row = E("A087", 2020, {"mildness": "4"})
    load(conn, {"/breeder_Year/145/2022/": season_page([row])})
    assert count(conn, "evaluation WHERE queen_code IS NULL") == 1
    # la même évaluation, vue ailleurs sous son code (page d'un autre éleveur)
    load(conn, {"/breeder_Year/984/2026/": season_page([E("FR-CBP-A-0087-2020", 2020, {"mildness": "4"})])})
    load(conn, {"/breederbyyear/145/2020/": _list(Q("FR-CBP-A-0087-2020", 2020, 145, name="A087"))})
    assert resolve_pending(conn) == 1
    rows = conn.execute("SELECT queen_code, queen_key FROM evaluation").fetchall()
    assert [tuple(r) for r in rows] == [("FR-CBP-A-0087-2020", "FR-CBP-A-0087-2020")]
```

- [ ] **Step 3 : Vérifier l'échec**

Run : `.venv/bin/pytest tests/mellifera/test_mellifera_ingest.py -q`
Expected : erreur de collecte `ModuleNotFoundError: No module named 'mellifera.ingest'`.

- [ ] **Step 4 : Écrire l'ingestion**

`mellifera/ingest.py` :

```python
"""Écriture des pages parsées dans la base (idempotente) et identification des reines évaluées (spec §6)."""
from __future__ import annotations

import re
import sqlite3
from typing import Optional

from .db import now_iso
from .parse.common import PARSER_VERSION, is_code, parse_number
from .parse.countries import CountryEntry, parse_countries
from .parse.country import BreederEntry, parse_country
from .parse.queen import QueenPage, parse_queen
from .parse.queens_list import QueensListPage, parse_queens_list
from .parse.season import SeasonPage, parse_season
from .urls import kind_of, path_ids

_TRIGRAM = re.compile(r"\(([A-Z0-9]{2,6})\)")


def mating_type_for(category_raw: Optional[str], description: Optional[str]) -> str:
    """Type d'accouplement d'un code : catégorie de la légende d'abord, description en repli."""
    for text in (category_raw, description):
        t = (text or "").lower()
        if "insemination" in t:
            return "insemination"
        if "mating station" in t or "belegstelle" in t:
            return "station"
    return "unknown"


def _upsert_mating_code(conn: sqlite3.Connection, code: Optional[str], description: Optional[str],
                        category_raw: Optional[str]) -> None:
    if not code:
        return
    conn.execute(
        "INSERT INTO mating_code (code, description, category_raw, mating_type) VALUES (?, ?, ?, ?) "
        "ON CONFLICT (code) DO UPDATE SET description = COALESCE(excluded.description, description), "
        "category_raw = COALESCE(excluded.category_raw, category_raw), mating_type = excluded.mating_type",
        (code, description, category_raw, mating_type_for(category_raw, description)))


def _ensure_criterion(conn: sqlite3.Connection, mode: str, code: str) -> None:
    """Colonne ou mode inconnu : critère créé sans trait (ignoré par l'export tant qu'il n'est pas mappé)."""
    conn.execute("INSERT OR IGNORE INTO criterion (mode, code) VALUES (?, ?)", (mode, code))


def store_countries(conn: sqlite3.Connection, entries: tuple[CountryEntry, ...]) -> None:
    for c in entries:
        conn.execute("INSERT INTO country (id, name, iso) VALUES (?, ?, ?) "
                     "ON CONFLICT (id) DO UPDATE SET name = excluded.name, iso = excluded.iso", (c.id, c.name, c.iso))


def _upsert_breeder(conn: sqlite3.Connection, breeder_id: int, trigram: str, country_id: Optional[int],
                    seen_at: str) -> None:
    conn.execute(
        "INSERT INTO breeder (id, trigram, country_id, first_seen, last_seen) VALUES (?, ?, ?, ?, ?) "
        "ON CONFLICT (id) DO UPDATE SET trigram = excluded.trigram, "
        "country_id = COALESCE(excluded.country_id, country_id), last_seen = MAX(last_seen, excluded.last_seen)",
        (breeder_id, trigram, country_id, seen_at, seen_at))


def store_country(conn: sqlite3.Connection, country_id: int, breeders: tuple[BreederEntry, ...], seen_at: str) -> None:
    for b in breeders:
        _upsert_breeder(conn, b.id, b.trigram, country_id, seen_at)


def store_queens_list(conn: sqlite3.Connection, page: QueensListPage, seen_at: str) -> None:
    """Éleveur + une ébauche par reine inconnue (la fiche, lue ensuite, fait foi et complète la ligne)."""
    _upsert_breeder(conn, page.breeder_id, page.trigram, page.country_id, seen_at)
    for r in page.rows:
        conn.execute(
            "INSERT OR IGNORE INTO queen (code, code_valid, name, birth_year, breeder_id, mother_code, "
            "drone_mother_code, serial_raw) VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            (r.code, int(is_code(r.code)), r.name, r.birth_year, page.breeder_id, r.mother_code,
             r.drone_mother_code, r.serial_raw))


def store_queen(conn: sqlite3.Connection, q: QueenPage, path: str, seen_at: str) -> None:
    _upsert_mating_code(conn, q.mating_code, q.mating_description, q.mating_category_raw)
    values = dict(
        code=q.code, code_valid=int(is_code(q.code)), name=q.name, birth_year=q.birth_year, breeder_id=q.breeder_id,
        mother_code=q.mother_code, drone_mother_code=q.drone_mother_code, n_drone_queens=q.n_drone_queens,
        mating_code=q.mating_code, single_drone=int(q.single_drone),
        insemination_date=q.insemination_date.isoformat() if q.insemination_date else None,
        race_code=q.race_code, race_description=q.race_description, lineage=q.lineage,
        generation_raw=q.generation_raw, coi=q.coi, coi_reliability=q.coi_reliability,
        not_found_raw=q.not_found_raw, cotation_count=q.cotation_count, serial_raw=q.serial_raw,
        origin_raw=q.origin_raw, data_origin_raw=q.data_origin_raw, mother_comment_raw=q.mother_comment,
        drone_comment_raw=q.drone_comment, parse_version=PARSER_VERSION, fetched_at=seen_at, parsed_at=now_iso(),
    )
    cols = ", ".join(values)
    marks = ", ".join("?" for _ in values)
    updates = ", ".join(f"{c} = excluded.{c}" for c in values if c != "code")
    conn.execute(f"INSERT INTO queen ({cols}) VALUES ({marks}) ON CONFLICT (code) DO UPDATE SET {updates}",
                 tuple(values.values()))
    for series in q.radar:
        for label, value in series.values.items():
            _ensure_criterion(conn, "radar", label)
            _upsert_evaluation(conn, dict(
                queen_name_raw=q.name, queen_birth_year_raw=q.birth_year, queen_code=q.code, queen_key=q.code,
                page_breeder_id=q.breeder_id, tester="", tester_breeder_id=None,
                eval_date=f"{series.year}-01-01", eval_year=series.year, mode="radar", criterion_code=label,
                value=value, value_raw=f"{value:g}", source="radar", url=path))


def _upsert_evaluation(conn: sqlite3.Connection, row: dict) -> None:
    cols = ", ".join(row)
    marks = ", ".join("?" for _ in row)
    conn.execute(f"INSERT INTO evaluation ({cols}) VALUES ({marks}) "
                 "ON CONFLICT (queen_key, source, mode, tester, eval_date, criterion_code) DO UPDATE SET "
                 "value = excluded.value, value_raw = excluded.value_raw, url = excluded.url, "
                 "tester_breeder_id = excluded.tester_breeder_id, eval_year = excluded.eval_year",
                 tuple(row.values()))


def lookup_queen(conn: sqlite3.Connection, name: str, birth_year: Optional[int],
                 page_breeder_id: Optional[int]) -> Optional[str]:
    """Étapes 3 et 4 de la spec §6 : (nom, année, éleveur), puis (nom, année) si unique en base."""
    if birth_year is None:
        return None
    m = _TRIGRAM.search(name)
    if m:
        rows = conn.execute("SELECT q.code FROM queen q JOIN breeder b ON b.id = q.breeder_id "
                            "WHERE q.name = ? AND q.birth_year = ? AND b.trigram = ?",
                            (name, birth_year, m.group(1))).fetchall()
    else:
        rows = conn.execute("SELECT code FROM queen WHERE name = ? AND birth_year = ? AND breeder_id = ?",
                            (name, birth_year, page_breeder_id)).fetchall()
    if len(rows) == 1:
        return rows[0]["code"]
    rows = conn.execute("SELECT code FROM queen WHERE name = ? AND birth_year = ?", (name, birth_year)).fetchall()
    return rows[0]["code"] if len(rows) == 1 else None


def unresolved_key(name: str, birth_year: Optional[int], page_breeder_id: Optional[int]) -> str:
    """Clé d'une ligne non résolue ; l'éleveur de la page évite de fusionner deux « A087 » d'élevages différents."""
    return f"name:{name}|{'' if birth_year is None else birth_year}|{'' if page_breeder_id is None else page_breeder_id}"


def store_season(conn: sqlite3.Connection, page: SeasonPage, path: str, page_breeder_id: int) -> None:
    for mc in page.mating_codes:
        _upsert_mating_code(conn, mc.code, mc.description, mc.category_raw)
    for e in page.evaluations:
        code = e.queen_code or lookup_queen(conn, e.queen_name_raw, e.queen_birth_year, page_breeder_id)
        for criterion, raw in e.values.items():
            _ensure_criterion(conn, e.mode, criterion)
            _upsert_evaluation(conn, dict(
                queen_name_raw=e.queen_name_raw, queen_birth_year_raw=e.queen_birth_year, queen_code=code,
                queen_key=code or unresolved_key(e.queen_name_raw, e.queen_birth_year, page_breeder_id),
                page_breeder_id=page_breeder_id, tester=e.tester, tester_breeder_id=e.tester_breeder_id,
                eval_date=e.eval_date.isoformat() if e.eval_date else "", eval_year=e.eval_year, mode=e.mode,
                criterion_code=criterion, value=parse_number(raw), value_raw=raw, source="season_table", url=path))


def resolve_pending(conn: sqlite3.Connection) -> int:
    """Relance les étapes 3 et 4 sur les évaluations non résolues ; une ligne résolue prend son code pour clé
    et fusionne avec la ligne identique déjà enregistrée sous ce code. Retourne le nombre de reines résolues."""
    pending = conn.execute("SELECT DISTINCT queen_key, queen_name_raw, queen_birth_year_raw, page_breeder_id "
                           "FROM evaluation WHERE queen_code IS NULL").fetchall()
    resolved = 0
    for p in pending:
        code = lookup_queen(conn, p["queen_name_raw"], p["queen_birth_year_raw"], p["page_breeder_id"])
        if code is None:
            continue
        conn.execute("UPDATE OR IGNORE evaluation SET queen_code = ?, queen_key = ? WHERE queen_key = ?",
                     (code, code, p["queen_key"]))
        # Ce qui reste sous l'ancienne clé doublonne une ligne déjà enregistrée sous le code.
        conn.execute("DELETE FROM evaluation WHERE queen_key = ?", (p["queen_key"],))
        resolved += 1
    return resolved


def ingest(conn: sqlite3.Connection, path: str, body: bytes, seen_at: str):
    """Parse la page selon son type et l'enregistre. Lève ParseError si la page ne se parse pas."""
    kind = kind_of(path)
    ids = path_ids(path)
    if kind == "countries":
        parsed = parse_countries(body)
        store_countries(conn, parsed)
    elif kind == "country":
        parsed = parse_country(body)
        store_country(conn, int(ids[0]), parsed, seen_at)
    elif kind == "queens_list":
        parsed = parse_queens_list(body)
        store_queens_list(conn, parsed, seen_at)
    elif kind == "queen":
        parsed = parse_queen(body)
        store_queen(conn, parsed, path, seen_at)
    elif kind == "season":
        parsed = parse_season(body)
        store_season(conn, parsed, path, int(ids[0]))
    else:
        raise ValueError(f"chemin hors liste blanche : {path}")
    return parsed
```

- [ ] **Step 5 : Vérifier que les tests passent**

Run : `.venv/bin/pytest tests/mellifera/test_mellifera_ingest.py -q && .venv/bin/pytest -q | tail -1`
Expected : `16 passed`, puis `322 passed, 12 deselected`.

- [ ] **Step 6 : Commit**

```bash
git add mellifera/ingest.py tests/mellifera/fake_site.py tests/mellifera/test_mellifera_ingest.py
git commit -F - <<'EOF'
Ajoute l'ingestion des pages et l'identification des reines évaluées

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015J8yR1HR7mLSRtwnvEYap6
EOF
```

---

### Task 9 : Parcours du site (complet, incrémental, fermeture, reprise, relecture du cache)

Spec §4.1, §4.2, §4.4 et §6 (échecs de parsing, seuil de 5 %) ; écarts 10, 11, 15 et 17 ; Review Focus 2.

- **Complet.** Lecture de `/all/`, puis de chaque pays. Pour chaque éleveur, la liste de sa dernière saison donne la navigation vers toutes ses années ; on lit toutes les listes, toutes les fiches, puis toutes les saisons. Vient ensuite la **fermeture** : toute reine citée comme mère ou comme « Drones », sans fiche lue et pas encore demandée dans ce passage, est lue à son tour, jusqu'à épuisement.
- **Incrémental.** Mêmes pages pays. Un éleveur nouveau est parcouru en entier. Pour les autres, seules les listes et saisons des années ≥ Y − 2 sont lues, et une fiche n'est relue que si la reine est nouvelle ou si son nombre de cotations a changé dans la liste. La fermeture suit.
- **Reprise.** Chaque page est journalisée dans `page` avec le `run_id`. Un passage arrêté par `TooManyErrors` reste inachevé (`finished_at` nul). Le passage suivant de même mode le reprend et relit depuis le cache les pages déjà lues dans ce passage, sans requête.
- **Relecture.** `reparse([kind])` reconstruit la base depuis le cache, sans client. Les évaluations des types relus sont effacées puis recréées ; le reste est mis à jour en place, et les id des reines sont conservés.
- **Échecs de parsing.** Chaque page est ingérée dans un `SAVEPOINT`. Toute exception, hors `sqlite3.Error`, annule ses écritures et se note dans `page.parse_error` ; le crawl continue (écart 15). Au-delà de 5 % d'échecs sur un type de page, le rapport le signale (`failing_kinds`).
- **Statistiques** (`crawl_run.stats`, JSON) : `pages` (type → statut → nombre), `parse_errors`, `queens_added`, `queens_updated`, `resumed_pages`, `invalid_codes`, `nonconforming_codes`, `evaluations_resolved` et `evaluations_unresolved`.

**Files:**
- Create: `mellifera/crawl.py`
- Modify: `tests/mellifera/fake_site.py` (ajout en fin de fichier)
- Test: `tests/mellifera/test_mellifera_crawl.py`

**Interfaces:**
- Consumes :
  - `RawCache` (tâche 1) ; `PoliteClient` et `TooManyErrors` (tâche 2) ;
  - `now_iso` et `open_db` (tâche 7) ; `ingest` et `resolve_pending` (tâche 8) ;
  - `PARSER_VERSION`, `ParseError` et `BreederEntry` (tâche 4) ;
  - `KINDS`, `countries_path`, `country_path`, `kind_of`, `queen_path`, `queens_list_path` et `season_path` (tâche 1).
- Produces :
  - `mellifera.crawl` :
    - constantes `MAX_PARSE_FAILURE = 0.05` et `RECENT_YEARS = 2` ;
    - `RunReport(run_id: int, mode: str, stats: dict, stopped: bool = False, failing_kinds: dict[str, float])` ;
    - `Crawler(conn, client: Optional[PoliteClient], cache: RawCache, *, today: Callable[[], dt.date] = dt.date.today)`, avec `.full() -> RunReport`, `.incremental() -> RunReport` et `.reparse(kind: Optional[str] = None) -> RunReport`.
  - `fake_site` :
    - `DEMO` (dict des reines A, B, C, D, M, G) ;
    - `demo_site() -> FakeSite`, avec deux éleveurs français (CBP 145, NG 200). M et G ne figurent dans aucune liste : seule la fermeture les lit. FR-OLD-B-0001-2019 n'existe pas.

- [ ] **Step 1 : Ajouter le site de démonstration au faux site**

Ajouter à la fin de `tests/mellifera/fake_site.py`, précédé de deux lignes vides :

```python
# ------------------------------------------------------------------ site de démonstration (crawl, CLI)
DEMO = {
    "A": Q("FR-CBP-B-0001-2022", 2022, 145, name="B001(CBP)", mother="FR-CBP-B-0100-2020",
           drones="FR-OLD-B-0001-2019", number=8, cotations=1),
    "B": Q("FR-CBP-B-0002-2022", 2022, 145, mother="FR-CBP-B-0100-2020", cotations=1),
    "C": Q("FR-CBP-B-0003-2021", 2021, 145, name="C003"),
    "D": Q("FR-NG-B-0001-2024", 2024, 200, mother="FR-CBP-B-0001-2022", drones="FR-CBP-B-0003-2021", number=3),
    "M": Q("FR-CBP-B-0100-2020", 2020, 145, name="M100", mother="FR-CBP-B-0050-2018"),
    "G": Q("FR-CBP-B-0050-2018", 2018, 145, name="G050"),
}


def demo_site() -> FakeSite:
    """Deux éleveurs français (CBP 145, NG 200). M et G ne sont dans aucune liste : seule la fermeture les lit.
    FR-OLD-B-0001-2019 (côté mâle de A) n'existe pas sur le site."""
    a, b, c, d = DEMO["A"], DEMO["B"], DEMO["C"], DEMO["D"]
    return FakeSite({
        "/all/": countries_page([(8, "France", "FR", 2), (0, "Unknown", None, 0)]),
        "/countrybyyear/8/": country_page([(145, "CBP", [2022]), (200, "NG", [2024])]),
        "/countrybyyear/0/": country_page([]),
        "/breederbyyear/145/2022/": queens_list_page(145, "CBP", 8, [2022, 2021], [a, b]),
        "/breederbyyear/145/2021/": queens_list_page(145, "CBP", 8, [2022, 2021], [c]),
        "/breederbyyear/200/2024/": queens_list_page(200, "NG", 8, [2024], [d]),
        **{f"/queen/{q.code}/": queen_page(q) for q in DEMO.values()},
        "/breeder_Year/145/2022/": season_page(
            [E("B001(CBP)", 2022, {"mildness": "5"}), E("FR-CBP-B-0002-2022", 2022, {"mildness": "4"}),
             E("C003", 2021, {"mildness": "3"}), E("M100", 2020, {"mildness": "6"})],
            links=[(a.code, "B001(CBP)")]),
        "/breeder_Year/200/2024/": season_page([E("FR-CBP-B-0001-2022", 2022, {"mildness": "5"})]),
    })
```

- [ ] **Step 2 : Écrire les tests (en échec)**

`tests/mellifera/test_mellifera_crawl.py` :

```python
import datetime as dt

import pytest

import mellifera.ingest
from mellifera.cache import RawCache
from mellifera.client import PoliteClient
from mellifera.crawl import Crawler
from mellifera.db import open_db
from mellifera.urls import kind_of, queen_path
from tests.mellifera.fake_site import (DEMO, E, FakeClock, FakeSite, Q, country_page, db_snapshot, demo_site,
                                       queen_page, queens_list_page, season_page)

A, B, C = DEMO["A"], DEMO["B"], DEMO["C"]
QUEENS = tuple(DEMO.values())

TODAY = dt.date(2024, 6, 1)


def make_crawler(tmp_path, site: FakeSite, max_errors: int = 5):
    conn = open_db(tmp_path / "m.sqlite")
    clock = FakeClock()
    client = PoliteClient("https://index-mellifera.org", "c@example.org", 2.0, max_errors,
                          transport=site.transport(), sleep=clock.sleep, clock=clock.clock)
    return Crawler(conn, client, RawCache(tmp_path / "cache"), today=lambda: TODAY), conn


def parsed_codes(conn) -> set[str]:
    return {r[0] for r in conn.execute("SELECT code FROM queen WHERE parsed_at IS NOT NULL")}


def test_full_crawl_with_closure(tmp_path):
    site = demo_site()
    crawler, conn = make_crawler(tmp_path, site)
    report = crawler.full()
    assert not report.stopped and report.failing_kinds == {}
    assert parsed_codes(conn) == {q.code for q in QUEENS}
    assert conn.execute("SELECT status FROM page WHERE url = '/queen/FR-OLD-B-0001-2019/'").fetchone()[0] == "missing"
    assert all(kind_of(p) for p in site.paths)
    assert len(site.paths) == len(set(site.paths))              # aucune page demandée deux fois
    assert report.stats["evaluations_unresolved"] == 0          # « M100 » résolue après la fermeture
    assert report.stats["evaluations_resolved"] == 1
    assert conn.execute("SELECT COUNT(*) FROM evaluation WHERE queen_code = ?", (A.code,)).fetchone()[0] == 1
    assert conn.execute("SELECT finished_at IS NOT NULL FROM crawl_run").fetchone()[0] == 1


def test_resume_after_network_stop(tmp_path):
    site = demo_site()
    site.failures[queen_path(B.code)] = [503] * 5
    crawler, conn = make_crawler(tmp_path, site, max_errors=5)
    first = crawler.full()
    assert first.stopped
    assert conn.execute("SELECT finished_at FROM crawl_run").fetchone()[0] is None
    before = set(site.paths)
    conn.close()

    crawler, conn = make_crawler(tmp_path, site)
    site.requests.clear()
    second = crawler.full()
    assert second.run_id == first.run_id and not second.stopped
    already = before - {queen_path(B.code)}
    assert not already & set(site.paths)                        # rien de ce qui était lu n'est redemandé
    assert parsed_codes(conn) == {q.code for q in QUEENS}

    fresh, fresh_conn = make_crawler(tmp_path / "fresh", demo_site())
    fresh.full()
    assert db_snapshot(conn)["evaluation"] == db_snapshot(fresh_conn)["evaluation"]


def test_incremental_reads_only_new_or_changed(tmp_path):
    site = demo_site()
    crawler, conn = make_crawler(tmp_path, site)
    crawler.full()
    conn.close()

    b2 = Q(B.code, 2022, 145, mother=B.mother, cotations=2)           # nouvelle cotation
    n = Q("FR-CBP-B-0004-2022", 2022, 145)                          # nouvelle reine
    o = Q("FR-ABC-B-0001-2015", 2015, 300)                          # reine d'un nouvel éleveur
    site.pages.update({
        "/countrybyyear/8/": country_page([(145, "CBP", [2022]), (200, "NG", [2024]), (300, "ABC", [2015])]),
        "/breederbyyear/145/2022/": queens_list_page(145, "CBP", 8, [2022, 2021], [A, b2, n]),
        queen_path(b2.code): queen_page(b2), queen_path(n.code): queen_page(n), queen_path(o.code): queen_page(o),
        "/breederbyyear/300/2015/": queens_list_page(300, "ABC", 8, [2015], [o]),
        "/breeder_Year/300/2015/": season_page([E(o.code, 2015, {"mildness": "5"}, tester="ABC", tester_id=300)]),
    })
    site.requests.clear()
    crawler, conn = make_crawler(tmp_path, site)
    report = crawler.incremental()
    assert not report.stopped
    queens = {p for p in site.paths if kind_of(p) == "queen"}
    assert queens == {queen_path(c) for c in (B.code, n.code, o.code, "FR-OLD-B-0001-2019")}
    assert "/breederbyyear/145/2021/" not in site.paths             # année < Y − 2
    assert {"/breeder_Year/145/2022/", "/breeder_Year/200/2024/", "/breederbyyear/300/2015/",
            "/breeder_Year/300/2015/"} <= set(site.paths)           # nouvel éleveur : parcouru en entier
    assert conn.execute("SELECT cotation_count FROM queen WHERE code = ?", (B.code,)).fetchone()[0] == 2
    assert n.code in parsed_codes(conn) and o.code in parsed_codes(conn)


def test_reparse_rebuilds_the_same_database_without_network(tmp_path):
    site = demo_site()
    crawler, conn = make_crawler(tmp_path, site)
    crawler.full()
    before = db_snapshot(conn)
    report = Crawler(conn, None, RawCache(tmp_path / "cache"), today=lambda: TODAY).reparse()
    assert not report.stopped and report.failing_kinds == {}
    assert db_snapshot(conn) == before


def test_parse_failures_are_recorded_and_reported(tmp_path):
    site = demo_site()
    site.pages[queen_path(C.code)] = b"<html><body><h1>Server error page</h1></body></html>"
    crawler, conn = make_crawler(tmp_path, site)
    report = crawler.full()
    assert report.failing_kinds == {"queen": pytest.approx(1 / 6, abs=1e-3)}
    err = conn.execute("SELECT parse_error FROM page WHERE url = ?", (queen_path(C.code),)).fetchone()[0]
    assert "fiche absente" in err
    assert len(parsed_codes(conn)) == 5                             # le crawl a continué


def test_unexpected_parser_exception_does_not_block_the_crawl(tmp_path, monkeypatch):
    real = mellifera.ingest.store_queen

    def buggy(conn, q, path, seen_at):
        real(conn, q, path, seen_at)
        if q.code == C.code:
            raise AttributeError("'NoneType' object has no attribute 'get'")   # bogue face à une structure inattendue
    monkeypatch.setattr(mellifera.ingest, "store_queen", buggy)
    crawler, conn = make_crawler(tmp_path, demo_site())
    report = crawler.full()
    assert not report.stopped                                       # passage terminé : la reprise ne bute pas dessus
    assert report.failing_kinds == {"queen": pytest.approx(1 / 6, abs=1e-3)}
    err = conn.execute("SELECT parse_error FROM page WHERE url = ?", (queen_path(C.code),)).fetchone()[0]
    assert err.startswith("AttributeError")
    assert C.code not in parsed_codes(conn)                         # écritures de la page annulées


def test_nonconforming_codes_are_kept_and_counted(tmp_path):
    odd = Q("B4536(NG)", 2024, 200)                                 # code mis en URL, hors format international
    site = demo_site()
    site.pages["/breederbyyear/200/2024/"] = queens_list_page(200, "NG", 8, [2024], [DEMO["D"], odd])
    site.pages[queen_path(odd.code)] = queen_page(odd)
    crawler, conn = make_crawler(tmp_path, site)
    report = crawler.full()
    assert report.stats["nonconforming_codes"] == 1
    assert conn.execute("SELECT code_valid FROM queen WHERE code = ?", (odd.code,)).fetchone()[0] == 0
```

- [ ] **Step 3 : Vérifier l'échec**

Run : `.venv/bin/pytest tests/mellifera/test_mellifera_crawl.py -q`
Expected : erreur de collecte `ModuleNotFoundError: No module named 'mellifera.crawl'`.

- [ ] **Step 4 : Écrire le crawler**

`mellifera/crawl.py` :

```python
"""Parcours du site : crawl complet, incrémental, fermeture sur les parents, reprise, relecture du cache."""
from __future__ import annotations

import datetime as dt
import json
import logging
import sqlite3
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Callable, Optional

from .cache import RawCache
from .client import PoliteClient, TooManyErrors
from .db import now_iso
from .ingest import ingest, resolve_pending
from .parse.common import PARSER_VERSION, ParseError
from .parse.country import BreederEntry
from .urls import KINDS, countries_path, country_path, kind_of, queen_path, queens_list_path, season_path

log = logging.getLogger("mellifera")

MAX_PARSE_FAILURE = 0.05        # au-delà de 5 % d'échecs de parsing sur un type de page : sortie en erreur
RECENT_YEARS = 2                # incrémental : années ≥ Y − 2


@dataclass
class RunReport:
    run_id: int
    mode: str
    stats: dict
    stopped: bool = False                                        # arrêt sur erreurs réseau consécutives
    failing_kinds: dict[str, float] = field(default_factory=dict)  # type → taux d'échec de parsing (> 5 %)


class Crawler:
    def __init__(self, conn: sqlite3.Connection, client: Optional[PoliteClient], cache: RawCache, *,
                 today: Callable[[], dt.date] = dt.date.today):
        self.conn, self.client, self.cache = conn, client, cache
        self.today = today
        self.run_id = 0
        self.started_at = ""
        self.pages: dict[str, Counter] = defaultdict(Counter)    # type → statut → nombre
        self.parse_errors: Counter = Counter()
        self.parsed_ok: Counter = Counter()
        self.counters: Counter = Counter()

    # ------------------------------------------------------------------ passage
    def _start(self, mode: str) -> None:
        row = self.conn.execute("SELECT id, started_at FROM crawl_run WHERE mode = ? AND finished_at IS NULL "
                                "ORDER BY id DESC LIMIT 1", (mode,)).fetchone()
        if row is not None and mode != "reparse":
            self.run_id, self.started_at = row["id"], row["started_at"]
            log.info("reprise du passage %s n° %s", mode, self.run_id)
        else:
            self.started_at = now_iso()
            cur = self.conn.execute("INSERT INTO crawl_run (mode, started_at) VALUES (?, ?)", (mode, self.started_at))
            self.run_id = cur.lastrowid
        self.conn.commit()

    def _stats(self) -> dict:
        return {
            "pages": {k: dict(v) for k, v in sorted(self.pages.items())},
            "parse_errors": dict(self.parse_errors),
            "queens_added": self.counters["queens_added"],
            "queens_updated": self.counters["queens_updated"],
            "resumed_pages": self.counters["resumed_pages"],
            "invalid_codes": self.counters["invalid_codes"],          # codes impossibles à mettre dans une URL
            "nonconforming_codes": self.conn.execute(                 # codes hors format international (gardés)
                "SELECT COUNT(*) FROM queen WHERE code_valid = 0").fetchone()[0],
            "evaluations_resolved": self.counters["evaluations_resolved"],
            "evaluations_unresolved": self.conn.execute(
                "SELECT COUNT(DISTINCT queen_key) FROM evaluation WHERE queen_code IS NULL").fetchone()[0],
        }

    def _finish(self, mode: str, stopped: bool) -> RunReport:
        if not stopped:
            self.counters["evaluations_resolved"] += resolve_pending(self.conn)
        stats = self._stats()
        failing = {k: round(self.parse_errors[k] / (self.parse_errors[k] + self.parsed_ok[k]), 3)
                   for k in self.parse_errors
                   if self.parse_errors[k] / (self.parse_errors[k] + self.parsed_ok[k]) > MAX_PARSE_FAILURE}
        self.conn.execute("UPDATE crawl_run SET finished_at = ?, stats = ? WHERE id = ?",
                          (None if stopped else now_iso(), json.dumps(stats, sort_keys=True), self.run_id))
        self.conn.commit()
        return RunReport(self.run_id, mode, stats, stopped, failing)

    # ------------------------------------------------------------------ pages
    def _fetch(self, path: str) -> tuple[Optional[bytes], str]:
        """Corps de la page (None si absente ou refusée) et date de lecture. Une page déjà lue dans ce passage
        (reprise) est relue depuis le cache, sans requête."""
        kind = kind_of(path)
        row = self.conn.execute("SELECT status, cache_path, fetched_at, run_id FROM page WHERE url = ?",
                                (path,)).fetchone()
        if row is not None and row["run_id"] == self.run_id and row["status"] != "error":
            self.counters["resumed_pages"] += 1
            self.pages[kind][row["status"]] += 1
            body = self.cache.get(row["cache_path"]) if row["status"] == "ok" else None
            return body, row["fetched_at"]
        resp = self.client.get(path)
        fetched_at = now_iso()
        rel = sha = None
        if resp.status == "ok":
            rel, sha = self.cache.put(path, resp.body)
        self.conn.execute(
            "INSERT INTO page (url, kind, status, http_status, fetched_at, sha256, cache_path, run_id) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?) ON CONFLICT (url) DO UPDATE SET status = excluded.status, "
            "http_status = excluded.http_status, fetched_at = excluded.fetched_at, "
            "sha256 = COALESCE(excluded.sha256, sha256), cache_path = COALESCE(excluded.cache_path, cache_path), "
            "parse_error = NULL, run_id = excluded.run_id",
            (path, kind, resp.status, resp.http_status, fetched_at, sha, rel, self.run_id))
        self.pages[kind][resp.status] += 1
        return (resp.body if resp.status == "ok" else None), fetched_at

    def _parse(self, path: str, body: bytes, seen_at: str):
        """Parse et enregistre la page. Un échec (ParseError, ou autre exception d'un parseur devant une structure
        inattendue) annule les écritures de la page, est noté dans page.parse_error et n'arrête pas le crawl :
        sinon la reprise relirait la même page en cache et s'arrêterait au même endroit, indéfiniment."""
        kind = kind_of(path)
        self.conn.execute("SAVEPOINT ingest")
        try:
            parsed = ingest(self.conn, path, body, seen_at)
        except sqlite3.Error:
            raise                                   # base indisponible : ce n'est pas un défaut de la page
        except Exception as exc:
            self.conn.execute("ROLLBACK TO ingest")
            error = str(exc) if isinstance(exc, ParseError) else f"{type(exc).__name__}: {exc}"
            self.conn.execute("UPDATE page SET parse_error = ?, parse_version = ? WHERE url = ?",
                              (error, PARSER_VERSION, path))
            self.parse_errors[kind] += 1
            log.warning("échec de parsing %s : %s", path, error, exc_info=not isinstance(exc, ParseError))
            parsed = None
        else:
            self.conn.execute("UPDATE page SET parse_error = NULL, parse_version = ? WHERE url = ?",
                              (PARSER_VERSION, path))
            self.parsed_ok[kind] += 1
        self.conn.execute("RELEASE ingest")
        self.conn.commit()
        return parsed

    def _visit(self, path: str):
        body, seen_at = self._fetch(path)
        if body is None:
            self.conn.commit()
            return None
        return self._parse(path, body, seen_at)

    def _visit_queen(self, code: str, list_count: Optional[int] = None) -> None:
        path = queen_path(code)
        if path is None:
            self.counters["invalid_codes"] += 1
            return
        before = self.conn.execute("SELECT parsed_at FROM queen WHERE code = ?", (code,)).fetchone()
        if self._visit(path) is None:
            return
        self.counters["queens_updated" if before is not None and before["parsed_at"] else "queens_added"] += 1
        if list_count is not None:
            self.conn.execute("UPDATE queen SET list_cotation_count = ? WHERE code = ?", (list_count, code))
            self.conn.commit()

    # ------------------------------------------------------------------ parcours
    def _breeder(self, b: BreederEntry, min_year: Optional[int]) -> None:
        """Listes de reines (toutes les années, ou années ≥ min_year), fiches à lire, saisons."""
        if min_year is None:        # parcours complet : dernière saison connue (page qui existe)
            first = b.season_years[0] if b.season_years else self.today().year
        else:                       # incrémental : année courante ; sa navigation donne les autres années
            first = self.today().year
        first_page = self._visit(queens_list_path(b.id, first))
        years = set(first_page.years) if first_page is not None else set(b.season_years)
        lists = [first_page] if first_page is not None else []
        for y in sorted(years - {first}, reverse=True):
            if min_year is None or y >= min_year:
                page = self._visit(queens_list_path(b.id, y))
                if page is not None:
                    lists.append(page)
        for page in lists:
            for row in page.rows:
                known = self.conn.execute("SELECT parsed_at, list_cotation_count FROM queen WHERE code = ?",
                                          (row.code,)).fetchone()
                stale = (known is None or known["parsed_at"] is None
                         or known["list_cotation_count"] != row.cotation_count)
                if min_year is None or stale:
                    self._visit_queen(row.code, row.cotation_count)
        for y in b.season_years:
            if min_year is None or y >= min_year:
                self._visit(season_path(b.id, y))

    def _closure(self) -> None:
        """Lit toute reine citée comme mère ou « Drones » et absente de la base, jusqu'à épuisement."""
        while True:
            refs = [r[0] for r in self.conn.execute(
                "SELECT ref FROM (SELECT mother_code AS ref FROM queen UNION SELECT drone_mother_code FROM queen) "
                "WHERE ref IS NOT NULL AND ref NOT IN (SELECT code FROM queen WHERE parsed_at IS NOT NULL) "
                "ORDER BY ref")]
            todo = []
            for code in refs:
                path = queen_path(code)
                if path is None:
                    continue
                seen = self.conn.execute("SELECT run_id FROM page WHERE url = ?", (path,)).fetchone()
                if seen is None or seen["run_id"] != self.run_id:
                    todo.append(code)
            if not todo:
                return
            for code in todo:
                self._visit_queen(code)

    def _countries(self, new_breeder: Callable[[BreederEntry], bool], min_year: Optional[int]) -> None:
        for c in self._visit(countries_path()) or ():
            for b in self._visit(country_path(c.id)) or ():
                self._breeder(b, None if new_breeder(b) else min_year)

    def _run(self, mode: str, body: Callable[[], None]) -> RunReport:
        self._start(mode)
        try:
            body()
        except TooManyErrors as exc:
            log.error("arrêt : %s ; le prochain passage %s reprendra ici", exc, mode)
            self.conn.commit()
            return self._finish(mode, stopped=True)
        return self._finish(mode, stopped=False)

    def full(self) -> RunReport:
        def body():
            self._countries(lambda b: True, None)
            self._closure()
        return self._run("full", body)

    def incremental(self) -> RunReport:
        def is_new(b: BreederEntry) -> bool:
            row = self.conn.execute("SELECT first_seen FROM breeder WHERE id = ?", (b.id,)).fetchone()
            return row is None or row["first_seen"] >= self.started_at

        def body():
            self._countries(is_new, self.today().year - RECENT_YEARS)
            self._closure()
        return self._run("incremental", body)

    def reparse(self, kind: Optional[str] = None) -> RunReport:
        """Reconstruit la base depuis le cache (aucune requête). Les évaluations des types relus sont
        effacées puis recréées ; les autres tables sont mises à jour en place (id des reines conservés)."""
        def body():
            kinds = [k for k in KINDS if kind is None or k == kind]
            if "season" in kinds:
                self.conn.execute("DELETE FROM evaluation WHERE source = 'season_table'")
            if "queen" in kinds:
                self.conn.execute("DELETE FROM evaluation WHERE source = 'radar'")
            for k in kinds:
                for row in self.conn.execute("SELECT url, cache_path, fetched_at FROM page "
                                             "WHERE kind = ? AND status = 'ok' ORDER BY url", (k,)).fetchall():
                    self.pages[k]["ok"] += 1
                    self._parse(row["url"], self.cache.get(row["cache_path"]), row["fetched_at"])
        return self._run("reparse", body)
```

- [ ] **Step 5 : Vérifier que les tests passent**

Run : `.venv/bin/pytest tests/mellifera/test_mellifera_crawl.py -q && .venv/bin/pytest -q | tail -1`
Expected : `7 passed`, puis `329 passed, 12 deselected`.

- [ ] **Step 6 : Commit**

```bash
git add mellifera/crawl.py tests/mellifera/fake_site.py tests/mellifera/test_mellifera_crawl.py
git commit -F - <<'EOF'
Ajoute le parcours du site : complet, incrémental, fermeture, reprise

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015J8yR1HR7mLSRtwnvEYap6
EOF
```

---

### Task 10 : Export `BlupV2Input` d'une sélection d'éleveurs

Spec §7 ; écarts 4, 7, 8, 9 et 14 ; Review Focus 1.

**Éleveurs.** Chaque `--breeder` est un trigramme (casse indifférente) ou un id Index Mellifera. Un trigramme porté par plusieurs éleveurs est refusé, avec la liste des candidats. Sont sélectionnées les reines de ces éleveurs nées à partir de `from_year`, puis toute leur ascendance présente en base (mères et « Drones »), quel que soit l'éleveur. Les ascendants hors sélection ont `available_as_mother = available_as_drone_source = false`.

**Reines.**
- `id` = `queen.id` (entier interne positif) ; `born` = `birth_year`, et une année hors [1900, Y+1] exclut la reine (comptée) ; `variety` = `race_description`.
- `mother` = id de la mère si elle est exportée, sinon `mother_ref_ext` = son code.
- `mating` :
  - `type` vient de la vue `queen_mating`, mais un type `unknown` avec côté mâle connu devient `insemination` (compté) ;
  - `drone_queens_mother` = id de la mère des RPM si elle est exportée, sinon `external_ref` = son code ;
  - `n_drone_queens` est repris tel quel ; `n_drones = 1` si `single_drone` ; `year` = `birth_year`.
- Les liens erronés sont coupés et comptés : auto-références, parent plus jeune, cycles (écart 9).

**Enregistrements.**
- Ne sont exportés que les traits mappés de kind `note` ou `percentage`.
- Sont écartés, et comptés par motif : critère non mappé, récolte, valeur non numérique, absence de date, valeur hors échelle, absence de rucher.
- Un pourcentage est divisé par 100 et ramené à l'échelle [0, 1].
- `apiary` = `tester_breeder_id`, ou l'éleveur de la reine pour le graphique ; `evaluator` = `tester_breeder_id` ; `hive_class = "hive"`.
- Doublons et graphique : écart 8.

**Rapport et validation.** Le résultat est validé par `BlupV2Input.model_validate` ; `ExportReport.lines()` produit le rapport imprimé par la CLI. Le module n'importe de `blup2` que `schemas` et `pedigree` : un test le vérifie dans un processus séparé, car l'image `mellifera` n'a pas numpy.

**Files:**
- Create: `mellifera/export.py`
- Test: `tests/mellifera/test_mellifera_export.py`

**Interfaces:**
- Consumes :
  - `blup2.schemas.BORN_MIN`, `BlupV2Input` et `QueenIn` ; `blup2.pedigree.find_cycle(queens: list[QueenIn]) -> list[int]` (moteur v2, inchangés) ;
  - `blup2.engine.run_blup`, dans les tests seulement ;
  - vue `queen_mating`, tables `breeder`, `country`, `evaluation`, `criterion` et `trait` (tâche 7) ;
  - `ingest` et `open_db` (tâches 7 et 8) ; `queen_path` (tâche 1) ; `E`, `Q`, `queen_page`, `queens_list_page` et `season_page` (tâches 5 et 6).
- Produces :
  - `mellifera.export` :
    - `SOURCES = ("season_table", "radar")` ;
    - `LABELS` (motif → libellé français) ;
    - `ExportReport(breeders: list[str], n_queens: int, n_selected: int, n_records: int, counts: Counter)`, avec `.lines() -> list[str]`, dont la première ligne est « éleveurs : CBP (145, France) » ;
    - `export_blup(conn, breeders: Iterable[str], *, from_year: Optional[int] = None, sources: Iterable[str] = ("season_table",), profile: str = "standard") -> tuple[BlupV2Input, ExportReport]`, qui lève `ValueError` pour un éleveur inconnu ou ambigu, une source inconnue ou une sélection vide.

- [ ] **Step 1 : Écrire les tests (en échec)**

`tests/mellifera/test_mellifera_export.py` :

```python
import subprocess
import sys
from pathlib import Path

import pytest

from blup2.engine import run_blup
from blup2.schemas import BlupV2Input
from mellifera.db import open_db
from mellifera.export import export_blup
from mellifera.ingest import ingest
from mellifera.urls import queen_path
from tests.mellifera.fake_site import E, Q, queen_page, queens_list_page, season_page

ROOT = Path(__file__).resolve().parents[2]
SEEN = "2026-10-13T10:00:00.000000+00:00"
OTHER = ("xx", "Some place", "Other")         # code d'accouplement sans type reconnu
H = Q("FR-NG-B-0009-2016", 2016, 200, name="H009")                                     # autre éleveur (ascendance)
G = Q("FR-CBP-B-0050-2018", 2018, 145, mother="FR-ZZZ-B-0001-2010", drones=H.code, number=6)
M = Q("FR-CBP-B-0100-2020", 2020, 145, mother=G.code)
A = Q("FR-CBP-B-0001-2022", 2022, 145, mother=M.code, drones="FR-OLD-B-0001-2019", number=8,
      radar={2023: [6, 5, 6, 5, 6]})
B = Q("FR-CBP-B-0002-2022", 2022, 145, name="B002(CBP)1dr", mother=M.code, drones=G.code, mating=OTHER)
X = Q("FR-CBP-B-0010-2019", 2019, 145, mother="FR-CBP-B-0011-2019")                   # cycle X ↔ Y
Y = Q("FR-CBP-B-0011-2019", 2019, 145, mother="FR-CBP-B-0010-2019")
Z = Q("FR-CBP-B-0012-2019", 2019, 145, mother="FR-CBP-B-0012-2019")                   # sa propre mère
P = Q("FR-CBP-B-0013-2018", 2018, 145, mother=A.code)                                 # mère plus jeune
W = Q("FR-CBP-B-0014-2019", 2019, 145)                    # sa mère n'a pas d'année valide (voir la fixture)
D = Q("FR-NG-B-0001-2024", 2024, 200, mother=A.code)                                  # hors sélection
QUEENS = (H, G, M, A, B, X, Y, Z, P, W, D)


@pytest.fixture
def conn(tmp_path):
    c = open_db(tmp_path / "m.sqlite")
    pages = {queen_path(q.code): queen_page(q) for q in QUEENS}
    pages["/breederbyyear/145/2022/"] = queens_list_page(145, "CBP", 8, [2022], [A, B])
    pages["/breederbyyear/200/2024/"] = queens_list_page(200, "NG", 8, [2024], [D])
    pages["/breeder_Year/145/2023/"] = season_page([
        E(A.code, 2022, {"mildness": "5", "hyg24h": "90", "smr": "3"}),
        E(A.code, 2022, {"mildness": "4", "honeyAmt": "48"}, mode="ceta"),
        E(B.code, 2022, {"mildness": "9", "swarming": "-", "framHold": "bof"}),
        E(M.code, 2020, {"mildness": "6"}, date="-"),
    ])
    for path, body in pages.items():
        ingest(c, path, body, SEEN)
    # reine sans année valide (ligne de liste illisible) citée comme mère
    c.execute("INSERT INTO queen (code, code_valid, birth_year, breeder_id) VALUES ('FR-CBP-B-0099-2017', 1, NULL, 145)")
    c.execute("UPDATE queen SET mother_code = 'FR-CBP-B-0099-2017' WHERE code = ?", (W.code,))
    c.commit()
    yield c
    c.close()


def by_code(conn, inp: BlupV2Input) -> dict:
    ids = dict(conn.execute("SELECT id, code FROM queen").fetchall())
    return {ids[q.id]: q for q in inp.queens}


def test_export_is_valid_and_runs(conn):
    inp, report = export_blup(conn, ["CBP"])
    assert isinstance(inp, BlupV2Input) and inp.exploitation == 0
    out = run_blup(inp)
    assert {q.id for q in out.queens} == {q.id for q in inp.queens}
    assert all(q.id > 0 for q in inp.queens)


def test_selection_and_ancestry(conn):
    inp, report = export_blup(conn, ["CBP"])
    queens = by_code(conn, inp)
    assert D.code not in queens                                   # autre éleveur, pas un ancêtre
    assert H.code in queens                                       # ancêtre d'un autre éleveur
    assert (queens[H.code].available_as_mother, queens[H.code].available_as_drone_source) == (False, False)
    assert queens[A.code].available_as_mother is True
    assert report.n_selected == 9 and report.n_queens == 10


def test_pedigree_links(conn):
    inp, report = export_blup(conn, ["CBP"])
    queens = by_code(conn, inp)
    a, b, g = queens[A.code], queens[B.code], queens[G.code]
    assert a.mother == queens[M.code].id
    assert (a.mating.type, a.mating.external_ref, a.mating.n_drone_queens, a.mating.year) == (
        "insemination", "FR-OLD-B-0001-2019", 8, 2022)                # côté mâle absent de la base
    assert (g.mother, g.mother_ref_ext) == (None, "FR-ZZZ-B-0001-2010")
    assert g.mating.drone_queens_mother == queens[H.code].id
    assert (b.mating.type, b.mating.drone_queens_mother, b.mating.n_drones) == ("insemination", g.id, 1)
    assert report.counts["mating_type_inferred"] >= 1
    assert queens[W.code].mother_ref_ext == "FR-CBP-B-0099-2017"  # mère sans année valide
    assert report.counts["invalid_birth_year"] == 1


def test_errors_in_source_data_are_cut(conn):
    inp, report = export_blup(conn, ["CBP"])
    queens = by_code(conn, inp)
    assert queens[Z.code].mother is None and report.counts["self_reference_cut"] == 1
    assert queens[P.code].mother is None and report.counts["younger_parent_cut"] == 1
    assert report.counts["cycle_link_cut"] == 1
    assert (queens[X.code].mother is None) != (queens[Y.code].mother is None)   # un seul des deux liens coupé


def test_records_scales_and_exclusions(conn):
    inp, report = export_blup(conn, ["CBP"])
    names = {c.id: c.name for c in inp.criteria}
    assert set(names.values()) == {"mildness", "hyg24h"}
    kinds = {c.name: c.kind for c in inp.criteria}
    assert kinds == {"mildness": "note", "hyg24h": "percentage"}
    a_id = by_code(conn, inp)[A.code].id
    recs = {(names[r.criterion], r.queen): r for r in inp.records}
    mild = recs[("mildness", a_id)]
    assert (mild.value, mild.scale_min, mild.scale_max) == (4.5, 0.0, 6.0)   # evjp 5 et ceta 4 : moyenne
    assert (mild.apiary, mild.evaluator, mild.hive_class, str(mild.date)) == (145, 145, "hive", "2023-07-06")
    assert recs[("hyg24h", a_id)].value == 0.9
    c = report.counts
    assert (c["merged_duplicates"], c["harvest"], c["unmapped_criterion"], c["out_of_scale"], c["non_numeric"],
            c["no_date"]) == (1, 1, 1, 1, 1, 1)


def test_radar_source(conn):
    inp, report = export_blup(conn, ["CBP"], sources=["season_table", "radar"])
    names = {c.id: c.name for c in inp.criteria}
    a_id = by_code(conn, inp)[A.code].id
    radar = [r for r in inp.records if r.queen == a_id and str(r.date) == "2023-01-01"]
    assert {names[r.criterion] for r in radar} == {"swarming", "framHold", "winter", "broVital"}
    assert all(r.apiary == 145 and r.evaluator is None for r in radar)   # éleveur de la reine
    assert report.counts["radar_covered"] == 1                           # douceur déjà dans le tableau 2023


def test_from_year_filter_and_unknown_breeder(conn):
    inp, _ = export_blup(conn, ["CBP"], from_year=2022)
    queens = by_code(conn, inp)
    assert queens[A.code].available_as_mother and not queens[M.code].available_as_mother
    assert X.code not in queens
    with pytest.raises(ValueError, match="éleveur inconnu"):
        export_blup(conn, ["XYZ"])
    with pytest.raises(ValueError, match="sources inconnues"):
        export_blup(conn, ["CBP"], sources=["tables"])
    with pytest.raises(ValueError, match="aucune reine"):
        export_blup(conn, ["CBP"], from_year=2030)


def test_breeder_by_trigram_any_case_or_by_id(conn):
    reference, report = export_blup(conn, ["CBP"])
    assert report.lines()[0] == "éleveurs : CBP (145, pays inconnu)"      # pays absent de cette mini-base
    for ref in (["cbp"], ["145"], [" CBP ", "145"]):
        assert export_blup(conn, ref)[0] == reference


def test_trigram_shared_by_two_breeders_is_refused(conn):
    conn.execute("INSERT INTO country (id, name, iso) VALUES (3, 'Germany', 'DE')")
    conn.execute("INSERT INTO breeder (id, trigram, country_id, first_seen, last_seen) VALUES (999, 'CBP', 3, ?, ?)",
                 (SEEN, SEEN))
    with pytest.raises(ValueError, match=r"plusieurs éleveurs : 145 \(pays inconnu\), 999 \(Germany\)"):
        export_blup(conn, ["CBP"])
    inp, report = export_blup(conn, ["145"])
    assert report.n_selected == 9


def test_export_does_not_import_numpy():
    """L'image du service mellifera n'installe ni numpy ni scipy : l'export n'importe que blup2.schemas et
    blup2.pedigree."""
    code = ("import sys, mellifera.export; "
            "heavy = sorted(m for m in sys.modules if m.split('.')[0] in ('numpy', 'scipy')); "
            "assert not heavy, heavy")
    subprocess.run([sys.executable, "-c", code], cwd=ROOT, check=True)
```

- [ ] **Step 2 : Vérifier l'échec**

Run : `.venv/bin/pytest tests/mellifera/test_mellifera_export.py -q`
Expected : erreur de collecte `ModuleNotFoundError: No module named 'mellifera.export'`.

- [ ] **Step 3 : Écrire l'export**

`mellifera/export.py` :

```python
"""Export d'une sélection de la base au format BlupV2Input (spec §7)."""
from __future__ import annotations

import datetime as dt
import sqlite3
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Iterable, Optional

from blup2.pedigree import find_cycle
from blup2.schemas import BORN_MIN, BlupV2Input, QueenIn

SOURCES = ("season_table", "radar")
LABELS = {
    "invalid_birth_year": "reines sans année de naissance valide (exclues)",
    "self_reference_cut": "liens auto-référents coupés",
    "younger_parent_cut": "liens vers un parent plus jeune que la reine coupés",
    "cycle_link_cut": "liens coupés pour rompre un cycle",
    "mating_type_inferred": "accouplements « unknown » à côté mâle connu passés en insemination",
    "unmapped_criterion": "valeurs de critères non mappés (écartées)",
    "harvest": "récoltes en kg (écartées)",
    "non_numeric": "valeurs non numériques (écartées)",
    "no_date": "évaluations sans date (écartées)",
    "out_of_scale": "valeurs hors échelle (écartées)",
    "no_apiary": "évaluations sans testeur ni éleveur (écartées)",
    "radar_covered": "valeurs du graphique déjà couvertes par un tableau (écartées)",
    "merged_duplicates": "doublons fusionnés par moyenne",
}


@dataclass
class ExportReport:
    breeders: list[str] = field(default_factory=list)   # « CBP (145, France) »
    n_queens: int = 0
    n_selected: int = 0
    n_records: int = 0
    counts: Counter = field(default_factory=Counter)   # motifs d'exclusion et corrections, par nom

    def lines(self) -> list[str]:
        out = [f"éleveurs : {', '.join(self.breeders)}",
               f"reines : {self.n_queens} (dont {self.n_selected} sélectionnées)", f"enregistrements : {self.n_records}"]
        out += [f"  {LABELS.get(k, k)} : {v}" for k, v in sorted(self.counts.items())]
        return out


def _valid_year(y: Optional[int]) -> bool:
    return y is not None and BORN_MIN <= y <= dt.date.today().year + 1


def _breeders(conn: sqlite3.Connection, refs: Iterable[str]) -> list[sqlite3.Row]:
    """Éleveurs désignés par trigramme (sans tenir compte de la casse) ou par id Index Mellifera. Un trigramme
    porté par plusieurs éleveurs (pays différents) est refusé : sinon deux programmes seraient fusionnés."""
    sql = ("SELECT b.id, b.trigram, c.name AS country FROM breeder b LEFT JOIN country c ON c.id = b.country_id "
           "WHERE {} ORDER BY b.id")
    out: dict[int, sqlite3.Row] = {}
    for ref in (r.strip() for r in refs):
        if ref.isdigit():
            rows = conn.execute(sql.format("b.id = ?"), (int(ref),)).fetchall()
        else:
            rows = conn.execute(sql.format("UPPER(b.trigram) = UPPER(?)"), (ref,)).fetchall()
        if not rows:
            raise ValueError(f"éleveur inconnu : « {ref} »")
        if len(rows) > 1:
            listed = ", ".join(f"{r['id']} ({r['country'] or 'pays inconnu'})" for r in rows)
            raise ValueError(f"trigramme « {ref} » porté par plusieurs éleveurs : {listed} ; "
                             f"désignez l'éleveur par son id (--breeder {rows[0]['id']})")
        out[rows[0]["id"]] = rows[0]
    return list(out.values())


def _select(conn: sqlite3.Connection, breeder_ids: list[int], from_year: Optional[int]) -> dict[str, sqlite3.Row]:
    """Reines des éleveurs choisis, nées à partir de from_year."""
    marks = ", ".join("?" for _ in breeder_ids)
    sql = f"SELECT * FROM queen_mating WHERE breeder_id IN ({marks})"
    args: list = list(breeder_ids)
    if from_year is not None:
        sql += " AND birth_year >= ?"
        args.append(from_year)
    return {r["code"]: r for r in conn.execute(sql, args)}


def _ancestry(conn: sqlite3.Connection, selected: dict[str, sqlite3.Row]) -> dict[str, sqlite3.Row]:
    """Sélection + toute son ascendance présente en base (mères et « Drones »), quel que soit l'éleveur."""
    rows = dict(selected)
    frontier = list(selected.values())
    while frontier:
        refs = {c for r in frontier for c in (r["mother_code"], r["drone_mother_code"]) if c and c not in rows}
        frontier = []
        for code in sorted(refs):
            r = conn.execute("SELECT * FROM queen_mating WHERE code = ?", (code,)).fetchone()
            if r is not None:
                rows[code] = r
                frontier.append(r)
    return rows


def _queens(rows: dict[str, sqlite3.Row], selected: set[str], report: ExportReport) -> list[dict]:
    kept = {}
    for code, r in rows.items():
        if _valid_year(r["birth_year"]):
            kept[code] = r
        else:
            report.counts["invalid_birth_year"] += 1
    ids = {code: r["id"] for code, r in kept.items()}
    queens = []
    for code, r in sorted(kept.items(), key=lambda kv: kv[1]["id"]):
        mother = r["mother_code"]
        drone = r["drone_mother_code"]
        if mother == code:
            mother = None
            report.counts["self_reference_cut"] += 1
        if drone == code:
            drone = None
            report.counts["self_reference_cut"] += 1
        q: dict = {"id": r["id"], "born": r["birth_year"], "variety": r["race_description"]}
        if mother in ids and kept[mother]["birth_year"] > r["birth_year"]:
            mother = None                                   # mère plus jeune que sa fille : saisie erronée
            report.counts["younger_parent_cut"] += 1
        if mother in ids:
            q["mother"] = ids[mother]
        elif mother:
            q["mother_ref_ext"] = mother
        mtype = r["mating_type"]
        if drone and mtype == "unknown":
            mtype = "insemination"
            report.counts["mating_type_inferred"] += 1
        if drone in ids and kept[drone]["birth_year"] > r["birth_year"]:
            drone = None
            report.counts["younger_parent_cut"] += 1
        mating: dict = {"type": mtype, "year": r["birth_year"]}
        if drone in ids:
            mating["drone_queens_mother"] = ids[drone]
        elif drone:
            mating["external_ref"] = drone
        if r["n_drone_queens"]:
            mating["n_drone_queens"] = r["n_drone_queens"]
        if r["single_drone"]:
            mating["n_drones"] = 1
        q["mating"] = mating
        if code not in selected:
            q["available_as_mother"] = q["available_as_drone_source"] = False
        queens.append(q)
    _cut_cycles(queens, report)
    return queens


def _cut_cycles(queens: list[dict], report: ExportReport) -> None:
    """Coupe, tant qu'il en reste, un lien de la reine d'id maximal prise dans un cycle (mère d'abord)."""
    by_id = {q["id"]: q for q in queens}
    while True:
        cycle = find_cycle([QueenIn.model_validate(q) for q in queens])
        if not cycle:
            return
        members = set(cycle)
        for qid in sorted(cycle, reverse=True):
            q = by_id[qid]
            if q.get("mother") in members:
                q.pop("mother")
            elif q["mating"].get("drone_queens_mother") in members:
                q["mating"].pop("drone_queens_mother")
            else:
                continue
            report.counts["cycle_link_cut"] += 1
            break
        else:
            raise RuntimeError(f"cycle sans lien à couper : {sorted(cycle)}")


def _records(conn: sqlite3.Connection, codes: dict[int, str], breeders: dict[int, Optional[int]],
             sources: Iterable[str], report: ExportReport) -> tuple[list[dict], list[dict]]:
    marks = ", ".join("?" for _ in sources)
    rows = conn.execute(
        "SELECT e.*, t.id AS trait_id, t.name AS trait, t.kind, t.scale_min, t.scale_max, t.direction "
        "FROM evaluation e LEFT JOIN criterion c ON c.mode = e.mode AND c.code = e.criterion_code "
        f"LEFT JOIN trait t ON t.id = c.trait_id WHERE e.queen_code IS NOT NULL AND e.source IN ({marks}) ORDER BY e.id",
        list(sources)).fetchall()
    id_of = {code: qid for qid, code in codes.items()}
    cells: dict[tuple, list[float]] = defaultdict(list)
    traits: dict[int, sqlite3.Row] = {}
    for e in rows:
        qid = id_of.get(e["queen_code"])
        if qid is None:
            continue
        if e["trait_id"] is None:
            report.counts["unmapped_criterion"] += 1
            continue
        if e["kind"] == "harvest":
            report.counts["harvest"] += 1
            continue
        if e["value"] is None:
            report.counts["non_numeric"] += 1
            continue
        if not e["eval_date"]:
            report.counts["no_date"] += 1
            continue
        lo, hi = (e["scale_min"], e["scale_max"]) if e["kind"] == "note" else (0.0, 1.0)
        value = e["value"] / 100 if e["kind"] == "percentage" else e["value"]
        if not lo <= value <= hi:
            report.counts["out_of_scale"] += 1
            continue
        apiary = e["tester_breeder_id"] if e["source"] == "season_table" else breeders.get(qid)
        if apiary is None:
            report.counts["no_apiary"] += 1
            continue
        traits[e["trait_id"]] = e
        cells[(qid, e["trait_id"], apiary, e["eval_date"], e["source"], e["tester_breeder_id"])].append(value)
    season_years = {(k[0], k[1], k[3][:4]) for k in cells if k[4] == "season_table"}
    records = []
    for (qid, tid, apiary, date, source, evaluator), values in sorted(cells.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][3], kv[0][4], kv[0][2])):
        if source == "radar" and (qid, tid, date[:4]) in season_years:
            report.counts["radar_covered"] += 1
            continue
        if len(values) > 1:
            report.counts["merged_duplicates"] += len(values) - 1
        t = traits[tid]
        lo, hi = (t["scale_min"], t["scale_max"]) if t["kind"] == "note" else (0.0, 1.0)
        records.append({"queen": qid, "criterion": tid, "value": round(sum(values) / len(values), 4),
                        "scale_min": lo, "scale_max": hi, "date": date, "evaluator": evaluator, "apiary": apiary,
                        "hive_class": "hive"})
    criteria = [{"id": tid, "name": t["trait"], "kind": t["kind"], "direction": t["direction"]}
                for tid, t in sorted(traits.items())]
    return criteria, records


def export_blup(conn: sqlite3.Connection, breeders: Iterable[str], *, from_year: Optional[int] = None,
                sources: Iterable[str] = ("season_table",), profile: str = "standard") -> tuple[BlupV2Input, ExportReport]:
    """`breeders` : trigrammes ou id Index Mellifera. ValueError : éleveur inconnu ou ambigu, source inconnue,
    aucune reine sélectionnée."""
    sources = tuple(sources)
    bad = [s for s in sources if s not in SOURCES]
    if bad or not sources:
        raise ValueError(f"sources inconnues : {bad} (attendu : {', '.join(SOURCES)})")
    chosen = _breeders(conn, breeders)
    report = ExportReport(breeders=[f"{b['trigram']} ({b['id']}, {b['country'] or 'pays inconnu'})" for b in chosen])
    selected = _select(conn, [b["id"] for b in chosen], from_year)
    if not selected:
        raise ValueError(f"aucune reine en base pour {', '.join(report.breeders)}"
                         + (f" née à partir de {from_year}" if from_year is not None else ""))
    rows = _ancestry(conn, selected)
    queens = _queens(rows, set(selected), report)
    codes = {r["id"]: code for code, r in rows.items()}
    breeders = {r["id"]: r["breeder_id"] for r in rows.values()}
    kept = {q["id"] for q in queens}
    criteria, records = _records(conn, {i: c for i, c in codes.items() if i in kept}, breeders, sources, report)
    report.n_queens = len(queens)
    report.n_selected = sum(1 for q in queens if q.get("available_as_mother", True))
    report.n_records = len(records)
    inp = BlupV2Input.model_validate({"exploitation": 0, "profile": profile, "criteria": criteria,
                                      "queens": queens, "records": records})
    return inp, report
```

- [ ] **Step 4 : Vérifier que les tests passent**

Run : `.venv/bin/pytest tests/mellifera/test_mellifera_export.py -q && .venv/bin/pytest -q | tail -1`
Expected : `10 passed`, puis `339 passed, 12 deselected`.

- [ ] **Step 5 : Commit**

```bash
git add mellifera/export.py tests/mellifera/test_mellifera_export.py
git commit -F - <<'EOF'
Ajoute l'export BlupV2Input d'une sélection d'éleveurs

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015J8yR1HR7mLSRtwnvEYap6
EOF
```

---

### Task 11 : Planification hebdomadaire du crawl incrémental

Spec §9 (« `mellifera schedule` = incrémental chaque dimanche à 2 h (Europe/Paris, `MELLIFERA_SCHEDULE`) ») ; Review Focus 5.

- `parse_schedule("sun 02:00")` rend (jour 0–6, heure, minute) et lève `ValueError` sur un format invalide.
- `next_run(now, spec)` rend la prochaine échéance **strictement** après `now`, dans le fuseau de `now`. L'arithmétique se fait en heure locale : les jours de changement d'heure, l'échéance garde son heure murale.
- `run_schedule` attend par tranches d'1 h au plus, si bien qu'une mise en veille ou un changement d'heure ne fait pas manquer l'échéance. Il lance ensuite le job et recommence. Une exception du job est journalisée et la boucle continue (Review Focus 5). `iterations` borne la boucle pour les tests.

**Files:**
- Create: `mellifera/schedule.py`
- Test: `tests/mellifera/test_mellifera_schedule.py`

**Interfaces:**
- Consumes : rien (`zoneinfo` de la bibliothèque standard).
- Produces : `mellifera.schedule` :
  - `DAYS` ;
  - `parse_schedule(spec: str) -> tuple[int, int, int]` ;
  - `next_run(now: dt.datetime, spec: str) -> dt.datetime` ;
  - `run_schedule(spec: str, timezone: str, job: Callable[[], int], *, now: Optional[Callable[[], dt.datetime]] = None, sleep: Callable[[float], None] = time.sleep, iterations: Optional[int] = None) -> None`.

- [ ] **Step 1 : Écrire les tests (en échec)**

`tests/mellifera/test_mellifera_schedule.py` :

```python
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
```

- [ ] **Step 2 : Vérifier l'échec**

Run : `.venv/bin/pytest tests/mellifera/test_mellifera_schedule.py -q`
Expected : erreur de collecte `ModuleNotFoundError: No module named 'mellifera.schedule'`.

- [ ] **Step 3 : Écrire la planification**

`mellifera/schedule.py` :

```python
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
```

- [ ] **Step 4 : Vérifier que les tests passent**

Run : `.venv/bin/pytest tests/mellifera/test_mellifera_schedule.py -q && .venv/bin/pytest -q | tail -1`
Expected : `6 passed`, puis `345 passed, 12 deselected`.

- [ ] **Step 5 : Commit**

```bash
git add mellifera/schedule.py tests/mellifera/test_mellifera_schedule.py
git commit -F - <<'EOF'
Ajoute la planification hebdomadaire du crawl incrémental

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015J8yR1HR7mLSRtwnvEYap6
EOF
```

---

### Task 12 : Commandes `python -m mellifera`

Spec §9 ; écart 14 ; Review Focus 3 et 4.

| Commande | Rôle |
|---|---|
| `crawl --full` / `crawl --incremental` | crawl, sous verrou |
| `reparse [--kind <type>]` | relecture du cache, sous verrou, sans réseau |
| `status` | dernier passage, volumes par table, pages par type et statut, critères non mappés (lecture seule) |
| `export-blup --breeder … [--from] [--sources] [--profile] -o fichier.json` | export en lecture seule ; le JSON est écrit après validation (`model_dump_json(indent=1, exclude_none=True)`) et le rapport est imprimé |
| `schedule` | incrémental planifié |

- **Verrou.** `fcntl.flock` exclusif et non bloquant sur `data/mellifera/crawl.lock` : deux crawls, ou un crawl et un `reparse`, ne tournent jamais ensemble. `status` et `export-blup` ne prennent pas le verrou.
- **Configuration.** Les réglages sont lus dans le `try` : une variable invalide donne le code 4, pas une trace Python. Un contact vide donne aussi le code 4, avant toute création de base.
- **Codes de sortie.** 0 ; 1 si le passage s'est arrêté sur erreurs réseau ; 2 si un type de page dépasse 5 % d'échecs de parsing ; 3 si le verrou est pris ; 4 pour une erreur de configuration. Dans ce dernier cas, le message est journalisé (`MissingContact`, `FileNotFoundError`, `ValueError`).
- **Image légère.** L'export et la planification ne sont importés qu'à la demande. Un test vérifie, dans un processus séparé, que `import mellifera.cli` ne charge ni numpy ni scipy ; le test équivalent pour `mellifera.export` est à la tâche 10.

**Files:**
- Create: `mellifera/cli.py`, `mellifera/__main__.py`
- Test: `tests/mellifera/test_mellifera_cli.py`

**Interfaces:**
- Consumes :
  - `Settings` (tâche 1) ; `RawCache` (tâche 1) ; `MissingContact` et `PoliteClient` (tâche 2) ;
  - `Crawler` et `RunReport` (tâche 9) ; `connect_readonly` et `open_db` (tâche 7) ;
  - `export_blup` (tâche 10, importé à la demande) ; `run_schedule` (tâche 11, importé à la demande) ; `KINDS` (tâche 1) ;
  - `DEMO`, `FakeClock` et `demo_site` (tâches 2 et 9).
- Produces : `mellifera.cli` :
  - `EXIT_OK, EXIT_NETWORK, EXIT_PARSE, EXIT_LOCKED, EXIT_CONFIG = 0, 1, 2, 3, 4` ;
  - `LockBusy` et `crawl_lock(path: Path)` (gestionnaire de contexte) ;
  - `run_crawl(settings, mode, *, transport=None, sleep=None, clock=None) -> int` et `run_reparse(settings, kind) -> int` ;
  - `status_lines(settings) -> list[str]` et `run_export(settings, args) -> int` ;
  - `build_parser() -> argparse.ArgumentParser` et `main(argv: Optional[list[str]] = None, settings: Optional[Settings] = None) -> int`.

- [ ] **Step 1 : Écrire les tests (en échec)**

`tests/mellifera/test_mellifera_cli.py` :

```python
import json
import subprocess
import sys
from pathlib import Path

import pytest

from blup2.engine import run_blup
from blup2.schemas import BlupV2Input
from mellifera.cli import (EXIT_CONFIG, EXIT_LOCKED, EXIT_NETWORK, EXIT_OK, EXIT_PARSE, crawl_lock, main, run_crawl,
                           status_lines)
from mellifera.config import Settings
from tests.mellifera.fake_site import DEMO, FakeClock, demo_site

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture
def settings(tmp_path):
    return Settings(data_dir=tmp_path / "mellifera", contact="c@example.org")


def crawl(settings, site, mode="full"):
    clock = FakeClock()
    return run_crawl(settings, mode, transport=site.transport(), sleep=clock.sleep, clock=clock.clock)


def test_crawl_status_export_end_to_end(settings, tmp_path, capsys):
    assert crawl(settings, demo_site()) == EXIT_OK
    lines = status_lines(settings)
    assert lines[0].startswith("dernier passage : n° 1 full") and "terminé" in lines[0]
    assert "queen : 6" in lines and "evaluation : 4" in lines
    out = tmp_path / "cbp.json"
    assert main(["export-blup", "--breeder", "CBP", "-o", str(out)], settings=settings) == EXIT_OK
    inp = BlupV2Input.model_validate(json.loads(out.read_text()))
    assert len(inp.queens) == 5 and len(inp.records) == 4
    run_blup(inp)
    assert "reines : 5 (dont 5 sélectionnées)" in capsys.readouterr().out


def test_parse_failure_threshold_exit_code(settings):
    site = demo_site()
    site.pages[f"/queen/{DEMO['C'].code}/"] = b"<html><body>oops</body></html>"
    assert crawl(settings, site) == EXIT_PARSE


def test_network_stop_exit_code(settings):
    site = demo_site()
    site.failures["/all/"] = [503] * 5
    assert crawl(settings, site) == EXIT_NETWORK


def test_lock_held_by_a_crawl(settings, tmp_path):
    assert crawl(settings, demo_site()) == EXIT_OK
    with crawl_lock(settings.lock_path):                    # un crawl en cours (des heures pour le premier)
        assert crawl(settings, demo_site(), "incremental") == EXIT_LOCKED
        assert main(["reparse"], settings=settings) == EXIT_LOCKED
        assert main(["status"], settings=settings) == EXIT_OK          # lecture seule : jamais bloquée
        out = tmp_path / "cbp.json"
        assert main(["export-blup", "--breeder", "CBP", "-o", str(out)], settings=settings) == EXIT_OK
    assert crawl(settings, demo_site(), "incremental") == EXIT_OK


def test_missing_contact_refuses_to_start(tmp_path):
    s = Settings(data_dir=tmp_path / "m", contact="")
    assert main(["crawl", "--full"], settings=s) == EXIT_CONFIG
    assert main(["schedule"], settings=s) == EXIT_CONFIG
    assert not s.db_path.exists()


def test_invalid_environment_is_a_configuration_error(settings, monkeypatch):
    assert crawl(settings, demo_site()) == EXIT_OK
    monkeypatch.setenv("MELLIFERA_DATA", str(settings.data_dir))
    assert main(["status"]) == EXIT_OK                      # environnement valide : la base répond
    monkeypatch.setenv("MELLIFERA_TZ", "Mars/Olympus")
    assert main(["status"]) == EXIT_CONFIG
    monkeypatch.setenv("MELLIFERA_TZ", "Europe/Paris")
    monkeypatch.setenv("MELLIFERA_DELAY", "2s")
    assert main(["status"]) == EXIT_CONFIG


def test_status_and_export_without_database(settings, tmp_path):
    assert main(["status"], settings=settings) == EXIT_CONFIG
    assert main(["export-blup", "--breeder", "CBP", "-o", str(tmp_path / "x.json")], settings=settings) == EXIT_CONFIG


def test_reparse_command(settings):
    assert crawl(settings, demo_site()) == EXIT_OK
    assert main(["reparse"], settings=settings) == EXIT_OK
    assert main(["reparse", "--kind", "season"], settings=settings) == EXIT_OK


def test_cli_does_not_import_numpy():
    """L'image du service mellifera n'installe ni numpy ni scipy : la CLI se charge sans eux (l'export est
    vérifié par test_mellifera_export)."""
    code = ("import sys, mellifera.cli; "
            "heavy = sorted(m for m in sys.modules if m.split('.')[0] in ('numpy', 'scipy')); "
            "assert not heavy, heavy")
    subprocess.run([sys.executable, "-c", code], cwd=ROOT, check=True)
```

- [ ] **Step 2 : Vérifier l'échec**

Run : `.venv/bin/pytest tests/mellifera/test_mellifera_cli.py -q`
Expected : erreur de collecte `ModuleNotFoundError: No module named 'mellifera.cli'`.

- [ ] **Step 3 : Écrire les commandes**

`mellifera/cli.py` :

```python
"""Commandes : python -m mellifera crawl | reparse | status | export-blup | schedule."""
from __future__ import annotations

import argparse
import contextlib
import fcntl
import json
import logging
import time
from pathlib import Path
from typing import Iterator, Optional

from .cache import RawCache
from .client import MissingContact, PoliteClient
from .config import Settings
from .crawl import Crawler, RunReport
from .db import connect_readonly, open_db
from .urls import KINDS

EXIT_OK, EXIT_NETWORK, EXIT_PARSE, EXIT_LOCKED, EXIT_CONFIG = 0, 1, 2, 3, 4
log = logging.getLogger("mellifera")


class LockBusy(RuntimeError):
    pass


@contextlib.contextmanager
def crawl_lock(path: Path) -> Iterator[None]:
    """Verrou exclusif non bloquant : deux crawls (ou crawl et reparse) ne tournent jamais ensemble."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a") as f:
        try:
            fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise LockBusy(f"un crawl est déjà en cours ({path})") from exc
        try:
            yield
        finally:
            fcntl.flock(f, fcntl.LOCK_UN)


def _exit_code(report: RunReport) -> int:
    for kind, rate in sorted(report.failing_kinds.items()):
        log.error("échecs de parsing « %s » : %.1f %% (seuil 5 %%) : structure du site modifiée ?", kind, 100 * rate)
    log.info("passage %s n° %s : %s", report.mode, report.run_id, json.dumps(report.stats, sort_keys=True))
    if report.stopped:
        return EXIT_NETWORK
    return EXIT_PARSE if report.failing_kinds else EXIT_OK


def run_crawl(settings: Settings, mode: str, *, transport=None, sleep=None, clock=None) -> int:
    """Crawl complet ou incrémental sous verrou. transport/sleep/clock : injection pour les tests."""
    try:
        with crawl_lock(settings.lock_path):
            with PoliteClient(settings.base_url, settings.contact, settings.delay, settings.max_errors,
                              transport=transport, sleep=sleep or time.sleep, clock=clock or time.monotonic) as client:
                conn = open_db(settings.db_path)
                try:
                    crawler = Crawler(conn, client, RawCache(settings.cache_dir))
                    report = crawler.full() if mode == "full" else crawler.incremental()
                finally:
                    conn.close()
    except LockBusy as exc:
        log.error("%s", exc)
        return EXIT_LOCKED
    return _exit_code(report)


def run_reparse(settings: Settings, kind: Optional[str]) -> int:
    try:
        with crawl_lock(settings.lock_path):
            conn = open_db(settings.db_path)
            try:
                report = Crawler(conn, None, RawCache(settings.cache_dir)).reparse(kind)
            finally:
                conn.close()
    except LockBusy as exc:
        log.error("%s", exc)
        return EXIT_LOCKED
    return _exit_code(report)


def status_lines(settings: Settings) -> list[str]:
    conn = connect_readonly(settings.db_path)
    try:
        out = []
        run = conn.execute("SELECT * FROM crawl_run ORDER BY id DESC LIMIT 1").fetchone()
        if run is None:
            out.append("aucun passage")
        else:
            state = f"terminé le {run['finished_at']}" if run["finished_at"] else "inachevé (reprise au prochain passage)"
            out.append(f"dernier passage : n° {run['id']} {run['mode']}, commencé le {run['started_at']}, {state}")
        for table in ("country", "breeder", "queen", "evaluation", "mating_code", "criterion", "page"):
            out.append(f"{table} : {conn.execute(f'SELECT COUNT(*) FROM {table}').fetchone()[0]}")
        for r in conn.execute("SELECT kind, status, COUNT(*) AS n, SUM(parse_error IS NOT NULL) AS errors "
                              "FROM page GROUP BY kind, status ORDER BY kind, status"):
            out.append(f"pages {r['kind']} {r['status']} : {r['n']} (échecs de parsing : {r['errors']})")
        unmapped = conn.execute("SELECT COUNT(*) FROM criterion WHERE trait_id IS NULL").fetchone()[0]
        out.append(f"critères non mappés : {unmapped}")
        return out
    finally:
        conn.close()


def run_export(settings: Settings, args: argparse.Namespace) -> int:
    from .export import export_blup
    conn = connect_readonly(settings.db_path)
    try:
        inp, report = export_blup(conn, args.breeder, from_year=args.from_year,
                                  sources=[s.strip() for s in args.sources.split(",") if s.strip()],
                                  profile=args.profile)
    finally:
        conn.close()
    Path(args.output).write_text(inp.model_dump_json(indent=1, exclude_none=True), encoding="utf-8")
    for line in report.lines():
        print(line)
    return EXIT_OK


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="python -m mellifera", description="Crawler d'Index Mellifera")
    sub = p.add_subparsers(dest="command", required=True)
    crawl = sub.add_parser("crawl", help="lit le site (GET poli, liste blanche)")
    mode = crawl.add_mutually_exclusive_group(required=True)
    mode.add_argument("--full", action="store_const", const="full", dest="mode")
    mode.add_argument("--incremental", action="store_const", const="incremental", dest="mode")
    reparse = sub.add_parser("reparse", help="reconstruit la base depuis le cache, sans réseau")
    reparse.add_argument("--kind", choices=KINDS)
    sub.add_parser("status", help="dernier passage et volumes")
    exp = sub.add_parser("export-blup", help="exporte une sélection au format BlupV2Input")
    exp.add_argument("--breeder", action="append", required=True,
                     help="trigramme ou id Index Mellifera de l'éleveur (répétable)")
    exp.add_argument("--from", dest="from_year", type=int, help="année de naissance minimale des reines choisies")
    exp.add_argument("--sources", default="season_table", help="season_table, radar (séparés par des virgules)")
    exp.add_argument("--profile", choices=("standard", "breeder"), default="standard")
    exp.add_argument("-o", "--output", required=True)
    sub.add_parser("schedule", help="crawl incrémental hebdomadaire (MELLIFERA_SCHEDULE)")
    return p


def main(argv: Optional[list[str]] = None, settings: Optional[Settings] = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    args = build_parser().parse_args(argv)
    try:
        settings = settings or Settings.from_env()
        if args.command == "crawl":
            return run_crawl(settings, args.mode)
        if args.command == "reparse":
            return run_reparse(settings, args.kind)
        if args.command == "status":
            print("\n".join(status_lines(settings)))
            return EXIT_OK
        if args.command == "export-blup":
            return run_export(settings, args)
        from .schedule import run_schedule
        if not settings.contact:
            raise MissingContact("MELLIFERA_CONTACT est vide : renseignez un contact avant de planifier le crawl")
        run_schedule(settings.schedule, settings.timezone, lambda: run_crawl(settings, "incremental"))
        return EXIT_OK
    except (MissingContact, FileNotFoundError, ValueError) as exc:
        log.error("%s", exc)
        return EXIT_CONFIG
```

`mellifera/__main__.py` :

```python
import sys

from .cli import main

sys.exit(main())
```

- [ ] **Step 4 : Vérifier que les tests passent, puis les commandes à la main**

```bash
.venv/bin/pytest tests/mellifera/test_mellifera_cli.py -q
.venv/bin/pytest -q | tail -1
.venv/bin/python -m mellifera --help
CHECK="${TMPDIR:-/tmp}/mellifera-check"; rm -rf "$CHECK"
MELLIFERA_CONTACT= MELLIFERA_DATA="$CHECK" .venv/bin/python -m mellifera crawl --full; echo "code $?"
MELLIFERA_DATA="$CHECK" .venv/bin/python -m mellifera status; echo "code $?"
test ! -e "$CHECK/mellifera.sqlite" && echo "aucune base créée"
rm -rf "$CHECK"
```

Expected :
- `9 passed`, puis `354 passed, 12 deselected` ;
- l'aide liste `crawl`, `reparse`, `status`, `export-blup` et `schedule` ;
- `crawl --full` avec un contact vide journalise « MELLIFERA_CONTACT est vide : renseignez un contact (e-mail ou URL) avant de crawler » et affiche `code 4`, **sans aucune requête** ;
- `status` journalise « base introuvable : … » et affiche `code 4` ;
- puis `aucune base créée`.

Le `MELLIFERA_CONTACT=` explicite protège d'un contact déjà exporté dans le shell. Ne **jamais** lancer `crawl` avec un contact renseigné pendant l'exécution du plan : ce serait un vrai crawl du site.

- [ ] **Step 5 : Commit**

```bash
git add mellifera/cli.py mellifera/__main__.py tests/mellifera/test_mellifera_cli.py
git commit -F - <<'EOF'
Ajoute les commandes python -m mellifera

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015J8yR1HR7mLSRtwnvEYap6
EOF
```

---

### Task 13 : Image, service Compose, CI et documentation

Spec §3 (image Python légère, sans R) et §9 (service `mellifera`, volume `./data`, crawl initial lancé à la main) ; écarts 13 et 14.

- **Image.** python:3.12-slim avec `tzdata`, `TZ=Europe/Paris` et `MELLIFERA_DATA=/app/data/mellifera`. Elle ne copie que `blup2/` et `mellifera/`, et lance `python -m mellifera schedule`.
- **Service Compose.** Le service est rangé sous le profil `mellifera` : `docker compose up` ne le démarre pas, et un service visé explicitement (`build`, `run`) est activé. Il a `init: true`, pour que `docker stop` interrompe proprement Python lancé en PID 1, et `restart: unless-stopped`. `MELLIFERA_CONTACT` est lu dans `.env.production`, comme les variables de l'API.
- **CI.** Le job `docker` construit l'image et vérifie qu'elle importe la CLI et l'export sans numpy. Le job `python` installe déjà `requirements-dev.txt`, qui inclut désormais `requirements-mellifera.txt` depuis la tâche 1.

**Files:**
- Create: `Dockerfile-Mellifera`
- Modify: `docker-compose.yml`, `.github/workflows/tests.yml`, `README.md`, `README.fr.md`

**Interfaces:**
- Consumes : `python -m mellifera` (tâche 12) ; `requirements-mellifera.txt` (tâche 1).
- Produces : le service Compose `mellifera` (profil `mellifera`) et l'étape CI de construction.

- [ ] **Step 1 : Écrire l'image**

`Dockerfile-Mellifera` :

```dockerfile
# Crawler Index Mellifera : image Python légère, sans R ni numpy (le moteur n'y tourne pas).
FROM python:3.12-slim

RUN apt-get update && apt-get install -y --no-install-recommends tzdata \
    && rm -rf /var/lib/apt/lists/*

ENV TZ=Europe/Paris \
    PYTHONUNBUFFERED=1 \
    MELLIFERA_DATA=/app/data/mellifera

WORKDIR /app

COPY requirements-mellifera.txt .
RUN pip install --no-cache-dir -r requirements-mellifera.txt

# blup2 : seuls schemas.py et pedigree.py sont importés (validation de l'export BlupV2Input)
COPY blup2 ./blup2
COPY mellifera ./mellifera

CMD ["python", "-m", "mellifera", "schedule"]
```

- [ ] **Step 2 : Ajouter le service Compose**

`docker-compose.yml` : ajouter, sous le service `api` existant et à la même indentation, après une ligne vide :

```yaml
  # Crawler Index Mellifera : démarré seulement sur demande (profil), MELLIFERA_CONTACT requis.
  mellifera:
    profiles: ["mellifera"]
    build:
      context: .
      dockerfile: Dockerfile-Mellifera
    init: true
    restart: unless-stopped
    volumes:
      - ./data:/app/data
    env_file:
      - path: .env.production
        required: false
```

- [ ] **Step 3 : Ajouter la construction de l'image à la CI**

`.github/workflows/tests.yml`, job `docker` : insérer, juste après la ligne `      - run: docker compose build api` :

```yaml
      - run: docker compose build mellifera
      - name: Image du crawler sans numpy
        run: docker compose run --rm --no-deps mellifera python -c "import sys, mellifera.cli, mellifera.export; assert not {'numpy', 'scipy'} & set(sys.modules)"
```

- [ ] **Step 4 : Documenter (EN et FR)**

`README.md` : insérer la section suivante juste avant `## Project Structure`.

````markdown
## Index Mellifera crawler

The `mellifera/` package reads the public pages of [Index Mellifera](https://index-mellifera.org) (countries, breeders, queen lists, queen pages, yearly evaluation tables) with a polite crawl, stores them in a SQLite database (`data/mellifera/mellifera.sqlite`) and exports the queens of chosen breeders, with their whole ancestry, as a `BlupV2Input` file for the v2 engine. Design (in French): [docs/superpowers/specs/2026-10-13-index-mellifera-design.md](docs/superpowers/specs/2026-10-13-index-mellifera-design.md).

The crawler only sends GET requests to an allow-list of public paths, waits at least 1 s between two requests (2 s by default), slows down on errors and identifies itself with a contact address. The database keeps only the Index Mellifera id, trigram and country of each breeder. The raw page cache (`data/mellifera/cache/`) contains names: never commit or publish it.

| Variable | Default | Role |
|---|---|---|
| `MELLIFERA_CONTACT` | (empty) | e-mail or URL in the `User-Agent`; **required**: the crawl refuses to start without it |
| `MELLIFERA_DATA` | `data/mellifera` (`/app/data/mellifera` in the image) | database, raw cache and lock file |
| `MELLIFERA_DELAY` | `2.0` | seconds between two requests (at least 1) |
| `MELLIFERA_MAX_ERRORS` | `5` | consecutive network errors before a clean stop |
| `MELLIFERA_SCHEDULE` | `sun 02:00` | day and time of the weekly incremental crawl |
| `MELLIFERA_TZ` | `Europe/Paris` | time zone of the schedule |

With Docker, the `mellifera` service belongs to the `mellifera` profile: `docker compose up` does not start it. Put `MELLIFERA_CONTACT` in `.env.production`, then:

```bash
docker compose run --rm mellifera python -m mellifera crawl --full   # initial crawl: more than a day, resumed if interrupted
docker compose up -d mellifera                                       # weekly incremental crawl
docker compose run --rm mellifera python -m mellifera status
```

Commands (`python -m mellifera …`, or `.venv/bin/python -m mellifera …` outside Docker):

- `crawl --full` / `crawl --incremental`: crawl the site; an interrupted run is resumed by the next run of the same mode;
- `reparse [--kind countries|country|queens_list|season|queen]`: rebuild the database from the raw cache, without network;
- `status`: last run, row counts, pages by type and status, unmapped criteria;
- `export-blup --breeder CBP --breeder 200 [--from 2018] [--sources season_table,radar] [--profile breeder] -o data/mellifera/export.json`: export the queens of these breeders and their ancestry. `--breeder` takes a trigram (any case) or an Index Mellifera id; a trigram shared by several breeders is refused. Only the evaluation tables are exported by default; `radar` adds the yearly values of the queen page chart;
- `schedule`: weekly incremental crawl (`MELLIFERA_SCHEDULE`).

Exit codes: `0` success, `1` stopped after consecutive network errors, `2` more than 5 % of the pages of one type failed to parse (site structure changed?), `3` another crawl holds the lock, `4` configuration error (empty contact, invalid variable, missing database, unknown or ambiguous breeder).
````

`README.md`, liste « Project Structure » : ajouter après la ligne `` - `docker-compose.yml` : Configuration Docker Compose `` :

```markdown
- `Dockerfile-Mellifera` : Dockerfile of the Index Mellifera crawler image (no R, no numpy)
- `mellifera/` : Index Mellifera crawler (SQLite database, `BlupV2Input` export)
```

`README.md`, section « Tests » : ajouter après le paragraphe qui commence par « `pytest -m docker` builds » :

```markdown
The crawler tests (`tests/mellifera/`) never contact Index Mellifera: they use a fake site (`httpx.MockTransport`) and anonymized real pages (`tests/mellifera/fixtures/`). The raw pages these fixtures come from (`tests/mellifera/fixtures/raw/`, ignored by git) contain names and stay on the machine that downloaded them; `scripts/make_mellifera_fixtures.py` regenerates the fixtures from them.
```

`README.fr.md` : insérer la section suivante juste avant `## Structure du projet`.

````markdown
## Crawler Index Mellifera

Le paquet `mellifera/` lit les pages publiques d'[Index Mellifera](https://index-mellifera.org) (pays, éleveurs, listes de reines, fiches, tableaux d'évaluations de saison) par un crawl poli, les stocke dans une base SQLite (`data/mellifera/mellifera.sqlite`) et exporte les reines d'éleveurs choisis, avec toute leur ascendance, en fichier `BlupV2Input` pour le moteur v2. Conception : [docs/superpowers/specs/2026-10-13-index-mellifera-design.md](docs/superpowers/specs/2026-10-13-index-mellifera-design.md).

Le crawler n'envoie que des GET vers une liste blanche de chemins publics, attend au moins 1 s entre deux requêtes (2 s par défaut), ralentit en cas d'erreur et s'identifie par une adresse de contact. La base ne garde de chaque éleveur que son id Index Mellifera, son trigramme et son pays. Le cache brut des pages (`data/mellifera/cache/`) contient des noms : ne jamais le commiter ni le publier.

| Variable | Défaut | Rôle |
|---|---|---|
| `MELLIFERA_CONTACT` | (vide) | e-mail ou URL dans le `User-Agent` ; **obligatoire** : le crawl refuse de démarrer sans |
| `MELLIFERA_DATA` | `data/mellifera` (`/app/data/mellifera` dans l'image) | base, cache brut et fichier de verrou |
| `MELLIFERA_DELAY` | `2.0` | secondes entre deux requêtes (au moins 1) |
| `MELLIFERA_MAX_ERRORS` | `5` | erreurs réseau consécutives avant un arrêt propre |
| `MELLIFERA_SCHEDULE` | `sun 02:00` | jour et heure du crawl incrémental hebdomadaire |
| `MELLIFERA_TZ` | `Europe/Paris` | fuseau horaire de la planification |

Avec Docker, le service `mellifera` appartient au profil `mellifera` : `docker compose up` ne le démarre pas. Renseignez `MELLIFERA_CONTACT` dans `.env.production`, puis :

```bash
docker compose run --rm mellifera python -m mellifera crawl --full   # crawl initial : plus d'une journée, repris s'il est interrompu
docker compose up -d mellifera                                       # crawl incrémental hebdomadaire
docker compose run --rm mellifera python -m mellifera status
```

Commandes (`python -m mellifera …`, ou `.venv/bin/python -m mellifera …` hors Docker) :

- `crawl --full` / `crawl --incremental` : parcourt le site ; un passage interrompu est repris par le passage suivant de même mode ;
- `reparse [--kind countries|country|queens_list|season|queen]` : reconstruit la base depuis le cache brut, sans réseau ;
- `status` : dernier passage, volumes, pages par type et statut, critères non mappés ;
- `export-blup --breeder CBP --breeder 200 [--from 2018] [--sources season_table,radar] [--profile breeder] -o data/mellifera/export.json` : exporte les reines de ces éleveurs et leur ascendance. `--breeder` accepte un trigramme (casse indifférente) ou un id Index Mellifera ; un trigramme porté par plusieurs éleveurs est refusé. Seuls les tableaux d'évaluations sont exportés par défaut ; `radar` ajoute les valeurs annuelles du graphique des fiches ;
- `schedule` : crawl incrémental hebdomadaire (`MELLIFERA_SCHEDULE`).

Codes de sortie : `0` succès, `1` arrêt après des erreurs réseau consécutives, `2` plus de 5 % d'échecs de parsing sur un type de page (structure du site modifiée ?), `3` un autre crawl tient le verrou, `4` erreur de configuration (contact vide, variable invalide, base absente, éleveur inconnu ou ambigu).
````

`README.fr.md`, liste « Structure du projet » : ajouter après la ligne `` - `docker-compose.yml` : Configuration Docker Compose `` :

```markdown
- `Dockerfile-Mellifera` : Dockerfile de l'image du crawler Index Mellifera (sans R ni numpy)
- `mellifera/` : crawler Index Mellifera (base SQLite, export `BlupV2Input`)
```

`README.fr.md`, section « Tests » : ajouter après le paragraphe qui commence par « `pytest -m docker` construit » :

```markdown
Les tests du crawler (`tests/mellifera/`) ne contactent jamais Index Mellifera : ils utilisent un faux site (`httpx.MockTransport`) et des pages réelles anonymisées (`tests/mellifera/fixtures/`). Les pages brutes d'où viennent ces fixtures (`tests/mellifera/fixtures/raw/`, ignorées par git) contiennent des noms et restent sur la machine qui les a téléchargées ; `scripts/make_mellifera_fixtures.py` régénère les fixtures à partir d'elles.
```

- [ ] **Step 5 : Construire et contrôler l'image**

```bash
export COMPOSE_PROJECT_NAME=mellifera-check
docker compose config --services
docker compose build mellifera
docker compose run --rm --no-deps mellifera python -c "import sys, mellifera.cli, mellifera.export; assert not {'numpy', 'scipy'} & set(sys.modules); print('imports OK')"
docker compose run --rm --no-deps -e MELLIFERA_CONTACT= mellifera python -m mellifera schedule; echo "code $?"
docker compose run --rm --no-deps -e MELLIFERA_TZ=Mars/Olympus mellifera python -m mellifera status; echo "code $?"
COMPOSE_PROFILES=mellifera docker compose down --rmi local
unset COMPOSE_PROJECT_NAME
```

Expected :
- `docker compose config --services` n'affiche que `api` : le service `mellifera` est caché par son profil ;
- la construction réussit ;
- `imports OK` ;
- `schedule` sans contact journalise « MELLIFERA_CONTACT est vide … » et affiche `code 4` ;
- `status` avec un fuseau inconnu journalise « MELLIFERA_TZ inconnu : « Mars/Olympus » » et affiche `code 4` ;
- l'image `mellifera-check-mellifera` est supprimée.

Les deux `run` passent une variable explicite : sinon, un `.env.production` local pourrait fournir un contact et lancer un vrai crawl. Si Docker répond `open ~/.docker/buildx/current: permission denied` (fichier créé par un `sudo docker`), préfixer les commandes par `BUILDX_CONFIG="$TMPDIR/buildx"`, sans modifier `~/.docker`.

- [ ] **Step 6 : Vérification d'ensemble**

```bash
.venv/bin/pytest -q | tail -1
/opt/homebrew/bin/python3.12 -m venv "$TMPDIR/venv312" && "$TMPDIR/venv312/bin/pip" install -q -r requirements-dev.txt && "$TMPDIR/venv312/bin/python" -m pytest -q | tail -1
git ls-files tests/mellifera/fixtures/raw data/mellifera | wc -l
git status --short
```

Expected :
- `354 passed, 12 deselected`, deux fois (Python 3.14, puis 3.12) ;
- `0` fichier suivi sous `raw/` et `data/mellifera/` ;
- `git status` ne montre que les 5 fichiers de cette tâche.

Si `python3.12` est absent, la CI couvre 3.12.

- [ ] **Step 7 : Commit**

```bash
git add Dockerfile-Mellifera docker-compose.yml .github/workflows/tests.yml README.md README.fr.md
git commit -F - <<'EOF'
Ajoute l'image et le service Compose du crawler, la CI et la documentation

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015J8yR1HR7mLSRtwnvEYap6
EOF
```

Le premier crawl réel (`docker compose run --rm mellifera python -m mellifera crawl --full`) ne fait pas partie de ce plan. L'utilisateur le lance lui-même, après avoir renseigné `MELLIFERA_CONTACT`.
