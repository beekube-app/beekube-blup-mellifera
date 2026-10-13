# Index Mellifera — spécification (crawler, base, export BLUP, résolution des lignées)

Date : 13/10/2026. Statut : validée (13/10/2026).
Entrées : [`docs/superpowers/specs/2026-10-12-blup-v2-design.md`](2026-10-12-blup-v2-design.md)
(contrat v2, sémantique du pedigree §1) ; code du moteur v2 sur `main` (5f2170a).

Décisions prises (13/10/2026) :
- source : **Index Mellifera** (index-mellifera.org) par **crawl poli** des pages publiques ; la
  base Van Dyck (perso.unamur.be) est hors périmètre ;
- usage Beekube : Beekube ne garde que le **code** de la lignée (`parentRefExt`, `droneRefExt`) ;
  au calcul, le **moteur** ajoute les ancêtres Index Mellifera au pedigree comme reines de base.
  Aucune copie dans l'exploitation, aucune note Index Mellifera dans le BLUP ;
- fraîcheur : **crawl complet une fois, puis incrémental hebdomadaire** ; un code inconnu est
  ignoré (avertissement) jusqu'au passage suivant, pas de récupération à la demande ;
- architecture : **tout dans ce dépôt**, paquet Python `mellifera/` + base **SQLite** dans le
  volume `data/`, lue par le moteur ;
- deux volets, deux plans : **volet 1** crawler, base, export BLUP ; **volet 2** résolution des
  références externes dans `blup2`. Les deux partent de `main` (moteur v2 fusionné).

---

## 1. Objectifs et critères de réussite

1. **Tester le BLUP** sur des données réelles d'éleveurs Buckfast : pedigrees profonds,
   accouplements contrôlés, évaluations multi-critères.
2. **Alimenter Beekube** : quand un apiculteur indique une lignée connue (mère ou côté mâle par
   son code Index Mellifera), sa parenté, sa consanguinité et les suggestions d'accouplement
   tiennent compte de l'ascendance publiée.

Critères de réussite :
- le crawl initial va au bout (reprise comprise), avec un taux d'échec de parsing < 5 % par type
  de page ;
- la base contient toutes les reines des « Queens list » et la fermeture de leurs parents ;
- `mellifera export-blup --breeder CBP` produit un JSON qui passe `BlupV2Input` et que
  `blup2.engine.run_blup` traite sans erreur ;
- une entrée Beekube dont une reine porte `mother_ref_ext = "FR-CBP-B-0276-2022"` reçoit
  l'ascendance de cette reine dans son pedigree ; l'ensemble des reines en sortie, la base 100
  et les suggestions d'accouplement ne contiennent aucun ancêtre Index Mellifera.

## 2. La source

### 2.1 Pages publiques (sans compte)

| URL | Contenu |
|---|---|
| `/all/` | pays (id, nom, nombre d'éleveurs) |
| `/countrybyyear/<pays>/` | éleveurs (id, trigramme) et leurs saisons de pedigree (`/breeder_Year/<éleveur>/<année>/`) |
| `/breederbyyear/<éleveur>/<année>/` | « Queens list » : toutes les reines nées cette année chez l'éleveur — nom, année, mère, « Drones (RAM Mother) », nombre de cotations, série, code ; navigation vers les autres années |
| `/queen/<code>/` | fiche : code international, nom, année, mère, côté mâle (« Drones : X — Number : n »), éleveur, génération, consanguinité (COI, fiabilité, ancêtres introuvables), série, date et code d'insémination, race, lignée, origine, import, origine de la donnée, graphique des cotations (5 critères par année, JavaScript intégré) |
| `/breeder_Year/<éleveur>/<année>/` | pedigree de saison des reines reproductrices et tableaux d'évaluations détaillées (jusqu'à ~45 colonnes ; deux tableaux d'une même page n'ont pas les mêmes colonnes) |
| `/queenTree/<code>/<3..7>/` | arbre JSON (non utilisé : la fiche suffit) |

Non publiques (redirection vers `/accounts/login/`) : `/rating/<id>` (détail d'une cotation),
`/exportYearQueens/<éleveur>/<année>/` (export CSV).

Volume observé le 13/10/2026 : 32 pays (dont « Unknown », sans code ISO), 754 éleveurs ;
France = 85 éleveurs, 412 pages de saison ; ~56 reines par page. Estimation : ~3 700 listes,
~3 500 pages de saison, quelques dizaines de milliers de fiches.

### 2.2 Sémantique des champs

- **Côté mâle (convention Buckfast)** : « Drones : X — Number : n » sur la fiche de Q signifie
  que Q a été accouplée avec des mâles issus de **n filles de X**. Donc, pour Q :
  `drone_queens_mother = X`, `n_drone_queens = n` — pas `n_drones`. Vérifié : la fiche de
  FR-CBP-B-0276-2022 (« Number : 9 ») et son arbre (« ins, (9 C.) ») concordent. C'est
  l'accouplement **de Q** (père de ses filles et de ses ouvrières), conforme à la spec v2 §1.
- Suffixe `1dr` dans un nom (`B562(PJ)1dr`) : insémination mono-mâle.
- Code d'insémination (`ins`, `inng`, `insng`…) avec une description lisible sur la fiche
  (« Thorigné (FR) instrumental insemination »).
- Mère « Unknown » : reine de base.
- Code international : `PAYS-ÉLEVEUR-RACE-NUMÉRO-ANNÉE` (`FR-CBP-B-0276-2022`,
  `LU-PJ-V-005H-2015`).
- Évaluations détaillées : publiques seulement pour les reines présentes dans un pedigree de
  saison, c'est-à-dire les **reproductrices** ; les autres reines notées n'ont que le graphique
  (5 critères par année). **Biais de sélection** à garder en tête pour toute étude BLUP sur ces
  données.
- Qualité : les données sont saisies par les éleveurs (ex. lien mère `FR-NG-B-0098-2022` et
  commentaire `FR-NG-B-0298-2022` sur la même fiche). Rien n'est corrigé au crawl.

### 2.3 Cadre d'accès

`robots.txt` d'Index Mellifera : `Disallow: *` (seul `/welcome/` autorisé). Conditions
d'utilisation : consultation « gratuite et ouverte à tous », pedigrees « libre d'accès »,
réutilisation « à des fins de recherche sur les abeilles sous licence CETA Mellifera ».
L'utilisateur a choisi le crawl poli en connaissance de cette tension (13/10/2026). Le crawler
s'y tient strictement (§4.3).

Données personnelles : la base ne garde de l'éleveur que son **id Index Mellifera, son
trigramme et son pays**. Ni nom, ni adresse, ni auteur de cotation. Le cache brut contient ces
informations : il reste local au volume `data/`, n'est jamais commité ni publié. Les fixtures de
test sont anonymisées.

## 3. Architecture

```
mellifera/                 volet 1
  http.py                  client poli : délai, ralentissement, User-Agent, URL interdites, cache brut
  parse/                   un module par type de page ; fonctions pures octets → structures
    countries.py  country.py  queens_list.py  queen.py  season.py  common.py
  db.py                    schéma SQLite, écritures idempotentes
  crawl.py                 parcours complet, incrémental, fermeture, reprise
  export.py                export BlupV2Input
  cli.py                   mellifera crawl | reparse | status | export-blup | schedule
blup2/external.py          volet 2 : résolution des références externes
tests/mellifera/           tests du volet 1 (+ fixtures anonymisées)
tests/test_external.py     tests du volet 2
data/mellifera/            mellifera.sqlite, cache/   (volume, non commité)
```

Dépendances ajoutées à `requirements.txt` : `httpx`, `beautifulsoup4`. Le service `mellifera`
utilise une image Python légère (sans R).

Flux :

```
index-mellifera.org ──http.py──▶ cache brut ──parse/──▶ db.py ──▶ mellifera.sqlite
                                                                    │        │
                                                     export.py ◀────┘        └────▶ blup2/external.py
                                                  (BlupV2Input JSON)          (lecture seule, au calcul)
```

## 4. Crawl (volet 1)

### 4.1 Parcours complet

1. `/all/` → pays.
2. `/countrybyyear/<pays>/` → éleveurs et leurs saisons de pedigree.
3. Pour chaque éleveur, une première `/breederbyyear/<éleveur>/<année>/` (année tirée de ses
   saisons) donne la liste de ses années ; on lit toutes les années.
4. `/queen/<code>/` pour chaque reine des listes.
5. `/breeder_Year/<éleveur>/<année>/` pour chaque saison (évaluations détaillées).
6. **Fermeture** : toute reine citée comme mère ou comme « Drones » et absente de la base est
   lue à son tour, jusqu'à ce qu'il n'en reste plus.

Le crawl est **reprenable** : chaque page lue est journalisée (`page`) ; un redémarrage saute ce
qui est déjà en cache pour ce passage.

### 4.2 Incrémental (hebdomadaire)

Avec Y = année courante :
- toutes les pages pays (nouveaux éleveurs, nouvelles saisons) ;
- un nouvel éleveur est parcouru en entier (§4.1, étapes 3 à 5) ;
- `breederbyyear` et `breeder_Year` des années ≥ Y − 2 ;
- la fiche d'une reine est relue seulement si elle est nouvelle ou si son nombre de cotations a
  changé dans la liste ;
- fermeture.

Les autres pages ne sont jamais relues.

### 4.3 Politesse

- Une seule connexion ; délai de 2 s entre deux requêtes (`MELLIFERA_DELAY`).
- `User-Agent: BeekubeCrawler/1.0 (+<MELLIFERA_CONTACT>)` ; le crawl refuse de démarrer si
  `MELLIFERA_CONTACT` est vide.
- 429, 5xx ou délai dépassé : le délai double à chaque erreur (plafond 5 min), `Retry-After` est
  respecté ; arrêt propre après 5 erreurs consécutives (`MELLIFERA_MAX_ERRORS`). Un 404 n'est
  pas une erreur : la page est marquée `missing`.
- GET uniquement. Liste blanche de chemins : `/all/`, `/countrybyyear/`, `/breederbyyear/`,
  `/breeder_Year/`, `/queen/`. Toute autre URL (`/rating/`, `/accounts/`, `/exportYearQueens/`,
  `/api/`…) est refusée **avant** l'envoi. Une redirection vers `/accounts/login/` est traitée
  comme `forbidden` et n'est pas suivie.
- Langue forcée en anglais (cookie de langue) pour des libellés stables.

### 4.4 Cache brut

Chaque réponse est stockée compressée (`cache/<2 premiers caractères>/<sha1(url)>.html.gz`) ;
`page` garde l'URL, le type, le statut HTTP, la date, l'empreinte SHA-256 et le chemin. Seule la
dernière version est conservée. Le parseur ne lit que le cache : corriger un parseur, c'est
lancer `mellifera reparse`, pas re-crawler.

## 5. Modèle de données (SQLite, mode WAL)

Chaque champ normalisé est stocké à côté de sa valeur brute (`*_raw`).

| Table | Colonnes principales |
|---|---|
| `country` | `id` (id Index Mellifera), `name` |
| `breeder` | `id` (id Index Mellifera), `trigram`, `country_id`, `first_seen`, `last_seen` |
| `queen` | `id` (entier interne, autoincrément), `code` (unique), `name` (`B4536(NG)`), `birth_year`, `breeder_id`, `mother_code`, `drone_mother_code`, `n_drone_queens`, `mating_code`, `mating_type`, `single_drone`, `insemination_date`, `race_code`, `race_description`, `lineage`, `generation_raw`, `coi`, `coi_reliability`, `not_found_raw`, `cotation_count`, `serial_raw`, `origin_raw`, `data_origin_raw`, `mother_comment_raw`, `drone_comment_raw`, `parse_version`, `fetched_at`, `parsed_at` |
| `mating_code` | `code`, `description`, `mating_type` (déduit), `mating_type_override` (manuel) |
| `criterion` | `id` (entier), `code` (`mildness`, `honeyAmt`, `hyg24h`…), `label`, `kind` (`note`, `percentage`, `harvest`, `unknown`), `scale_min`, `scale_max`, `direction` |
| `evaluation` | `id`, `queen_name_raw`, `queen_birth_year_raw` (toujours remplis), `queen_code` (nullable : résolution §6), `queen_key` (non nul : `queen_code` si résolu, sinon `name:<nom>\|<année>`), `tester` (trigramme), `tester_breeder_id`, `eval_date`, `eval_year`, `criterion_code`, `value`, `value_raw`, `source` (`season_table`, `radar`), `url` ; unique (`queen_key`, `tester`, `eval_date`, `criterion_code`, `source`) — `queen_key` est non nul parce que SQLite laisse passer les doublons d'une contrainte d'unicité contenant NULL |
| `page` | `url` (clé), `kind`, `status` (`ok`, `missing`, `forbidden`, `error`), `http_status`, `fetched_at`, `sha256`, `cache_path`, `parse_error`, `run_id` |
| `crawl_run` | `id`, `mode` (`full`, `incremental`, `reparse`), `started_at`, `finished_at`, `stats` (JSON) |

- `mother_code` et `drone_mother_code` sont des codes **sans clé étrangère** : un parent peut ne
  pas (encore) être en base.
- `criterion` est initialisée à partir d'un fichier de correspondance versionné
  (`mellifera/criteria.csv`) ; une colonne inconnue rencontrée au parsing crée un critère
  `kind = unknown`, ignoré par l'export tant qu'il n'est pas mappé.
- Le moteur ouvre la base en lecture seule pendant que le crawler écrit (WAL).

## 6. Parsing et normalisation (volet 1)

- Un parseur par type de page, fonction pure : octets → structures, sans réseau ni base.
  `PARSER_VERSION` est enregistrée sur chaque ligne ; `mellifera reparse [--kind <type>]`
  reconstruit la base depuis le cache.
- HTML approximatif (`<p>` non fermés) : parseur tolérant
  (BeautifulSoup) ; repérage par **libellés** (« International Code », « Mother », « Drones »,
  « Number »…), jamais par position.
- Tableaux d'évaluations : colonnes lues d'après l'en-tête de chaque tableau ; `-` → absent ;
  dates Django (`Sept. 7, 2020`) → ISO.
- **Identification de la reine d'une ligne d'évaluation.** La cellule « Queen » est du texte
  sans lien : tantôt un code (`FR-CBP-A-0203-2020`), tantôt un nom affiché (`B943(CBP)`,
  `A087`). Résolution, dans cet ordre :
  1. la cellule correspond à l'expression régulière du code → ce code ;
  2. un lien `/queen/<code>/` **de la même page** a ce texte exact et un seul `href` → ce code
     (vérifié sur `breeder_Year/145/2021` et `2022` : toutes les lignes se résolvent ainsi) ;
  3. recherche dans `queen` par (`name`, `birth_year`, éleveur) — éleveur = trigramme entre
     parenthèses s'il y en a (`B131(SPD)` → SPD), sinon éleveur de la page ;
  4. recherche par (`name`, `birth_year`) si le résultat est unique en base ;
  5. sinon non résolue : `queen_code` nul, comptée dans le rapport, écartée de l'export.

  Les étapes 3 et 4 sont refaites à chaque `reparse` et en fin de crawl, car la fermeture peut
  faire entrer la reine en base plus tard ; une ligne qui se résout prend alors pour
  `queen_key` son code, et fusionne avec une éventuelle ligne identique déjà enregistrée sous ce
  code.
- Graphique des cotations : extrait du JavaScript intégré (`labels`, un jeu de valeurs par
  année) → `evaluation` avec `source = radar`, `eval_date` = 1er janvier de l'année du jeu
  (seule l'année est connue).
- Accouplement :
  - « Drones : X — Number : n » → `drone_mother_code = X`, `n_drone_queens = n` ;
  - `mating_code` → `mating_type` via la table `mating_code` : description contenant
    « insemination » → `insemination` ; « mating station » ou « Belegstelle » → `station` ;
    sinon `unknown` ; `mating_type_override` prime ;
  - suffixe `1dr` → `single_drone = 1`.
- Code international vérifié par expression régulière ; non conforme → gardé brut, compté dans
  le rapport.
- Une page qui ne se parse pas : `page.parse_error` renseigné, le crawl continue. Le rapport de
  fin compte les échecs par type ; au-delà de 5 % sur un type, la commande sort en erreur
  (changement probable de structure du site).

## 7. Export BLUP (volet 1)

`mellifera export-blup --breeder CBP [--breeder NG …] [--from 2015] [--sources season_table,radar] -o fichier.json`

Produit un `BlupV2Input` (profil `standard` par défaut, `--profile breeder` possible) :

- `exploitation` : 0.
- `queens` : reines des éleveurs choisis (nées à partir de `--from`) **et** toute leur
  ascendance présente en base, quel que soit l'éleveur.
  - `id` = `queen.id` ; `born` = `birth_year` (reine sans année valide : exclue, comptée) ;
  - `mother` = id de la mère si elle est en base, sinon `null` et `mother_ref_ext` = son code ;
  - `mating` : `type` = `mating_type`, sauf si `drone_mother_code` est connu et `mating_type`
    vaut `unknown` : `type` = `insemination` (déduit, compté dans le rapport ; sinon
    `_sire_spec` ignorerait l'accouplement) ; `drone_queens_mother` = id de
    `drone_mother_code` s'il est en base, sinon `external_ref` = ce code ; `n_drone_queens` ;
    `n_drones` = 1 si `single_drone`, sinon `null` ; `year` = `birth_year` ;
  - liens auto-référents ou cycliques dans les données Index Mellifera : coupés à l'export
    (comptés), pour que le fichier passe la validation `BlupV2Input` ;
  - `variety` = `race_description` ; reines de l'ascendance hors sélection :
    `available_as_mother = available_as_drone_source = false`.
- `criteria` : uniquement les critères mappés (`kind` ≠ `unknown`) et non `harvest` ; `id` =
  `criterion.id`, `direction` et `name` depuis la table.
- `records` :
  - `note` : valeur et échelle de `criterion` ;
  - `percentage` : valeur / 100, échelle [0, 1] ;
  - `harvest` (`honeyAmt`…) : **exclu** dans cette version (le contrat impose [0, 1], la mise à
    l'échelle est une décision statistique séparée) ;
  - `date` = `eval_date` ; `evaluator` = `tester_breeder_id` ; `apiary` = `tester_breeder_id`
    (pseudo-rucher : le site ne publie pas de rucher, d'où un groupe de contemporaines
    testeur × année) ; `hive_class = hive`.
  - valeur hors échelle : enregistrement écarté, compté.
- Le fichier de sortie est validé par `BlupV2Input` avant écriture ; un rapport (reines,
  enregistrements écartés par motif) est imprimé.

## 8. Résolution des références externes dans le moteur (volet 2)

### 8.1 Activation

- Variable `MELLIFERA_DB` (chemin de `mellifera.sqlite`) dans le service `api`, base ouverte en
  **lecture seule** (`mode=ro`).
- Nouvelle option `options.external_pedigree_depth` (entier ≥ 0, défaut 7 ; 0 = désactivé).
- Variable absente ou fichier illisible : aucune résolution ; si l'entrée contient des
  références, avertissement `external_pedigree_unavailable`.

### 8.2 Règles

Nouvelle étape `expand_external(inp)` appelée dans `engine._evaluate` **avant**
`build_pedigree`. Elle ne modifie pas `inp` : elle renvoie une **copie** des reines dont les
liens sont complétés, la liste des ancêtres ajoutés et des avertissements.

- Recherche d'un code : code international exact après normalisation (majuscules, espaces
  retirés) ; à défaut, nom affiché (`B4536(NG)`) s'il est **unique** en base ; sinon inconnu.
- **Mère** : reine avec `mother = null` et `mother_ref_ext` trouvé → `mother` = id de l'ancêtre.
- **Côté mâle** : reine dont `mating.drone_queens` est vide et `mating.external_ref` trouvé →
  interprété selon la **convention Index Mellifera** (§2.2) : `drone_queens_mother` = id de
  l'ancêtre ; `type` passe de `unknown` à `insemination` (sinon `_sire_spec` ignore
  l'accouplement) ; `external_ref` est conservé, ce qui fait partager le même nœud sire aux
  reines accouplées avec la même lignée ; `n_drone_queens` reste `null` (le moteur prend alors
  S = 1, règle v2 §2.4). L'aide de Beekube devra préciser que le code saisi côté mâle est celui
  de la **mère des reines productrices de mâles**, comme sur Index Mellifera.
- **Ascendance** : à partir de chaque ancêtre trouvé, on remonte `mother_code` et
  `drone_mother_code` jusqu'à `external_pedigree_depth` générations. Chaque ancêtre devient un
  `QueenIn` (`born` = `birth_year`, `mating` construit comme en §7) avec un **id négatif**
  (−1, −2, …) qui ne peut pas entrer en collision avec les id Beekube. Ancêtre sans année
  valide : la remontée s'arrête là (compté). Les liens vers des reines au-delà de la profondeur
  sont coupés (reines de base). Le `type` d'un ancêtre suit la même déduction qu'en §7
  (côté mâle connu et type `unknown` → `insemination`).
- Une même reine Index Mellifera atteinte par plusieurs chemins n'est ajoutée qu'une fois.
- **Les erreurs de saisie d'Index Mellifera ne font jamais échouer le calcul** : un lien
  auto-référent ou qui fermerait un cycle (parmi les ancêtres, ou entre un ancêtre et une reine
  Beekube) est coupé et compté dans `external_ref_resolved`. La validation `BlupV2Input` a lieu
  avant l'expansion ; sans cette règle, `build_pedigree` lèverait « pedigree cyclique » pour une
  erreur qui n'est pas dans les données de l'utilisateur.
- Les ancêtres reçoivent `available_as_mother = available_as_drone_source = false`.
- Les évaluations Index Mellifera ne sont **jamais** ajoutées.

### 8.3 Isolation des ancêtres

`base_birth_years`, `compute_indexes`, `queen_outputs`, `suggest_matings`, `genetic_trend` et
`connectivity` parcourent `inp.queens`. Les ancêtres ne doivent pas y entrer. Ils sont donc
passés **uniquement** à `build_pedigree` (nouveau paramètre `extra_queens`), qui les place dans
le graphe comme les reines de base (`added_queen_ids`), sans les compter dans
`ped.queen_ids`. Conséquences garanties par les tests :
- mêmes reines en sortie, avec ou sans base Index Mellifera ;
- même `base.birth_years` et même `base.n_queens` ;
- aucun id négatif dans `matings`, `queens` ou `trend` ;
- l'avertissement `unknown_parent_added` ne compte pas les ancêtres Index Mellifera (ils ont
  leur propre avertissement).

### 8.4 Avertissements et contrat

Nouveaux codes dans `WarningOut.code` : `external_ref_resolved` (nombre de références
résolues, nombre d'ancêtres ajoutés), `external_ref_unknown` (codes introuvables, au plus
20 cités), `external_pedigree_unavailable`.

`WarningOut.code` est une liste fermée publiée dans `openapi/openapi.json` : l'ajout de codes et
de l'option `external_pedigree_depth` modifie le contrat. Le snapshot OpenAPI est régénéré, et
**`MELLIFERA_DB` n'est activé en production qu'après régénération et publication du SDK PHP**,
pour qu'un SDK qui valide cette liste ne rejette pas une réponse.

## 9. Erreurs, journalisation, exploitation

- Commandes : `mellifera crawl --full | --incremental`, `mellifera reparse [--kind]`,
  `mellifera status` (dernier passage, erreurs, volumes par table), `mellifera export-blup`,
  `mellifera schedule`.
- Verrou fichier (`data/mellifera/crawl.lock`) : deux crawls ne tournent jamais ensemble.
- Journaux sur la sortie standard ; statistiques de chaque passage dans `crawl_run.stats`
  (pages par type et statut, reines ajoutées et modifiées, échecs de parsing).
- Codes de sortie : 0 succès ; 1 arrêt sur erreurs réseau consécutives ; 2 seuil d'échec de
  parsing dépassé ; 3 verrou déjà pris.
- `docker-compose.yml` : service `mellifera` (image Python légère, volume `./data`), commande
  `mellifera schedule` = incrémental chaque dimanche à 2 h (Europe/Paris, `MELLIFERA_SCHEDULE`).
  Le crawl initial se lance à la main : `docker compose run --rm mellifera mellifera crawl --full`.
- Service `api` : `MELLIFERA_DB=/app/data/mellifera/mellifera.sqlite`, activé selon §8.4.

## 10. Tests

Aucun test n'appelle le vrai site.

Volet 1 (`tests/mellifera/`) :
- **Parseurs**, sur les pages réelles téléchargées pendant l'étude, anonymisées (noms
  d'éleveurs et d'auteurs remplacés) : fiche avec et sans cotations (graphique compris) ; page
  de saison à deux tableaux de colonnes différentes ; « Queens list » ; page pays ; liste des
  pays ; HTML tronqué ou vide → erreur de parsing, pas d'exception non gérée.
- **Client HTTP** (`httpx.MockTransport`, horloge simulée) : délai entre requêtes ;
  ralentissement sur 429 et 5xx et respect de `Retry-After` ; arrêt après N erreurs ;
  `User-Agent` ; refus des URL hors liste blanche avant envoi ; redirection vers la connexion →
  `forbidden` ; démarrage refusé sans `MELLIFERA_CONTACT`.
- **Crawl** (site simulé à partir des fixtures) : reprise après interruption ; sélection
  incrémentale (reine relue seulement si nouvelle ou si son nombre de cotations change) ;
  fermeture sur les parents absents ; verrou.
- **Base** : écritures idempotentes ; dédoublonnage des évaluations entre pages de saison ;
  `reparse` donne la même base.
- **Identification des reines évaluées** (§6) : ligne par code ; ligne par nom résolue par un
  lien de la page ; ligne par nom avec trigramme résolue en base ; ligne par nom nu ambigu →
  non résolue ; une reine non résolue au premier passage est résolue après la fermeture ; la
  même évaluation vue sous son nom puis sous son code donne une seule ligne.
- **Export** : la sortie passe `BlupV2Input` et `run_blup` sur une mini-base issue des
  fixtures ; échelles (notes, pourcentages / 100, récoltes exclues) ; ascendance incluse ;
  mère absente → `mother_ref_ext` ; côté mâle connu et type inconnu → `insemination` ; cycle
  dans les données → coupé, export valide.

Volet 2 (`tests/test_external.py`, mini-base SQLite construite par le test) :
- résolution par code, par nom unique, nom ambigu → inconnu ;
- côté mâle : `type` `unknown` → `insemination`, `drone_queens_mother` renseigné, nœud sire
  partagé entre reines de même `external_ref` ;
- profondeur respectée, id négatifs, ancêtre partagé ajouté une fois, ancêtre sans année ;
- cycle et auto-référence dans la mini-base → liens coupés, calcul réussi ;
- isolation (§8.3) : sorties, base 100 et suggestions identiques en ensemble de reines avec et
  sans base ; la parenté entre deux reines Beekube rattachées à des lignées sœurs augmente ;
- base absente ou illisible → `external_pedigree_unavailable`, calcul inchangé ;
- snapshot OpenAPI à jour (`tests/test_openapi.py`).

## 11. Points à vérifier à l'implémentation

- Nom et valeur du cookie de langue (Django : `django_language=en` probable).
- Les valeurs du graphique d'une année sont-elles une moyenne de plusieurs cotations ou une
  seule ? (Comparer une reine à plusieurs cotations la même année avec la page de saison.)
- Échelle des notes : 1–6 ou 0–6 (le graphique va de 0 à 6) ; à fixer dans `criteria.csv`
  d'après les valeurs observées.
- Les années listées dans `/breederbyyear/` couvrent-elles toutes les reines d'un éleveur, y
  compris les années sans saison de pedigree ?
- Codes d'insémination et de station présents en base, pour compléter `mating_code`.

## 12. Hors périmètre

Base Van Dyck (source complémentaire possible pour les ancêtres anciens) ; récupération à la
demande d'un code inconnu ; mise à l'échelle des quantités récoltées ; affichage de l'ascendance
Index Mellifera dans Beekube et autocomplétion des codes ; injection des évaluations Index
Mellifera dans le BLUP (évaluation inter-exploitations, spec v2 §7) ; détail des cotations
(réservé aux membres).
