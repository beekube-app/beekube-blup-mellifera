# Index Mellifera, volet 3 — notes des ancêtres dans le BLUP et catalogue des caractères

Date : 15/10/2026. Statut : validée, puis amendée le 15/10/2026 après prototype (contemporaines, §4 et §5).
Entrées : [`2026-10-13-index-mellifera-design.md`](2026-10-13-index-mellifera-design.md) (crawler, base,
export, résolution des lignées : volets 1 et 2, fusionnés dans `main` à 632045e) ;
[`2026-10-12-blup-v2-design.md`](2026-10-12-blup-v2-design.md) (contrat v2, modèle) ; étude des sources
externes du 14/10/2026 (`docs/analyses/2026-10-14-etude-sources-externes.md`, hors git).

Décisions prises (15/10/2026) :
- objectif : **de meilleures estimations des reines de l'exploitation**, en empruntant l'information des
  parentes notées sur Index Mellifera. Les reines externes restent invisibles dans les sorties (isolement
  du volet 2 conservé). La comparaison des reines de l'exploitation aux lignées externes sur une même
  échelle est hors périmètre ;
- correspondance des critères : **déclarée par Beekube, critère par critère**, par un champ `trait` qui
  prend une valeur d'un **catalogue de caractères standard** figé dans le contrat (donc dans le SDK). Le
  catalogue unifie les caractères connus d'Index Mellifera, de BeeBreed et du registre Karl Kehrle, pour
  que Beekube puisse plus tard exporter vers ces outils ;
- reines qui apportent leurs notes : **les ancêtres résolus par le volet 2** (les notes
  détaillées du site ne sont publiques que pour les reproductrices, qui sont précisément les reines des
  pedigrees), **et leurs contemporaines** (amendement, voir ci-dessous) ;
- intégration : **au moment du calcul, dans le moteur** (`blup2`), qui lit les notes dans la base du
  crawler comme il lit déjà l'ascendance ; aucune logique de fusion côté Beekube ;
- valeurs d'élevage BeeBreed : jamais injectées dans le BLUP (une valeur d'élevage ne se convertit pas en
  note ; voir l'étude du 14/10). Les sources BeeBreed et Karl Kehrle n'entrent ici que par le catalogue.

**Amendement du 15/10/2026 (prototype).** Les notes des ancêtres seuls n'apportent rien : chaque note
tombe seule dans son groupe testeur × année, l'effet du groupe l'absorbe et le gain de fiabilité est nul.
Le moteur charge donc aussi les **contemporaines** : les autres reines notées dans le même groupe (même
testeur, ou même éleveur pour le radar, même année, même caractère). Elles entrent au pedigree comme les
ancêtres, ne sortent jamais, et rendent l'estimation peu sensible à la sélection des reproductrices
[Henderson 1975]. Mesures du prototype (§7) : fiabilité des reines issues des lignées + 0,068 en moyenne,
corrélation avec la vraie valeur 0,472 → 0,517. Deux règles de conception en découlent : aucun effet fixe
partagé entre notes externes et notes de l'exploitation (§5), et une méthode de fiabilité qui ne bascule
jamais en approché à cause des seules notes externes (§5).

---

## 1. Objectifs et critères de réussite

1. Une exploitation Beekube dont des reines citent des lignées Index Mellifera (`mother_ref_ext`,
   `mating.external_ref`) et dont des critères déclarent un `trait` voit la fiabilité de ses reines
   apparentées augmenter, sans qu'aucune reine externe n'apparaisse dans ses résultats.
2. Une exploitation sans `trait`, ou sans référence externe, ou sans `MELLIFERA_DB`, obtient des résultats
   **strictement identiques** à la version 2.1.0 (test bit à bit).
3. Aucune donnée de la base du crawler ne fait échouer un calcul : au pire, un avertissement et le calcul
   sans les notes externes.
4. Le contrat 2.2.0 expose le catalogue des caractères comme une énumération, avec ses correspondances
   dans la documentation, de sorte que le SDK PHP régénéré la porte.
5. Une simulation sur le banc de tests lents montre le gain attendu et l'absence d'effet sur les reines
   non apparentées, y compris quand les notes externes portent un biais de sélection.

## 2. Catalogue des caractères standard

Identifiant stable en anglais (`snake_case`), libellé français, type au sens du contrat (`note`,
`percentage`, `harvest`), et correspondance **une pour une** avec chaque source (une case vide : la
source ne publie pas ce caractère). Les index composites (`Gesamtzuchtwert` BeeBreed, `combination`
Karl Kehrle) n'en font pas partie : ce sont des index, pas des caractères.

| Standard | Libellé | Type | Index Mellifera (`trait.name`) | BeeBreed | Karl Kehrle |
|---|---|---|---|---|---|
| `honey_yield` | Récolte de miel (kg) | harvest | `honeyAmt` | Honey yield | `honeykgpromean` |
| `honey_yield_average` | Récolte de miel moyenne (kg) | harvest | `aveHonAmt` | | |
| `honey_score` | Indice miel (note) | note | `honeyIdx` | | `honeyindex` |
| `average_crop_score` | Indice de récolte moyen (note) | note | `aveCropIdx` | | |
| `spring_crop` | Récolte de printemps (note) | note | `sprCrop` | | `cropearly` |
| `mid_season_crop` | Récolte de mi-saison (note) | note | `midCrop` | | |
| `summer_crop` | Récolte d'été (note) | note | `sumCrop` | | `cropsumm` |
| `autumn_crop` | Récolte d'automne (note) | note | `falCrop` | | `croplate` |
| `pollen_crop` | Récolte de pollen (note) | note | `polCrop` | | |
| `gentleness` | Douceur | note | `mildness` | Defensive behavior | `temperofbees` |
| `calmness` | Tenue au cadre | note | `framHold` | Steadiness on comb | `calmbehaviour` |
| `swarming` | Essaimage (note haute = essaime peu) | note | `swarming` | Swarming drive | `swarmtendency` |
| `overwintering` | Hivernage | note | `winter` | Robustness in winter | `overwintering` |
| `spring_development` | Développement de printemps | note | | Development in Spring | |
| `colony_strength` | Force de la colonie | note | | Colony strength | |
| `brood_vitality` | Vitalité du couvain | note | `broVital` | | `vitalityofbrood` |
| `bee_vitality` | Vitalité des abeilles | note | `beeVital` | | `vitalityofbees` |
| `fertility` | Fécondité | note | `fertil` | | `fertility` |
| `comb_building` | Bâtisse | note | `bldCombs` | | `buildingcombs` |
| `propolis` | Propolis | note | `propolis` | | `buildingpropol` |
| `hygienic_24h` | Test d'aiguille, couvain nettoyé à 24 h (%) | percentage | `hyg24h` | rate of opened cells | `hygprozafter24h` |
| `hygienic_48h` | Test d'aiguille, couvain nettoyé à 48 h (%) | percentage | `hyg48h` | | `hygienicprozafter48h` |
| `hygienic_score` | Indice de comportement hygiénique (note) | note | `hygBehaIdx` | | `hygienicindex` |
| `varroa_score` | Indice varroa (note) | note | `varIdx` | Varroa index | |
| `varroa_infestation` | Infestation varroa (mesure) | percentage | | Varroa infestation development | `varroanumber` |
| `smr` | SMR / VSH (%) | percentage | | SMR | `varroavsh` |
| `recapping` | Recapping (%) | percentage | | Recapping | |
| `chalkbrood` | Couvain plâtré (note haute = résistant) | note | | Chalkbrood | |
| `nosema` | Nosémose (note haute = résistant) | note | | Nosemosis | `nosemaindex` |
| `cbpv` | Paralysie chronique (note haute = résistant) | note | | Chronic Bee Paralysis | |

Règles :
- le catalogue vit dans `blup2/traits.py` (`STANDARD_TRAITS`, dictionnaire ordonné identifiant → fiche
  `Trait` figée : libellé, type, nom Index Mellifera, nom BeeBreed, champ Karl Kehrle, sens `direction`).
  Le sens vaut 1 (valeur haute = mieux) sauf pour `varroa_infestation` (−1) ; la description OpenAPI le
  mentionne (« valeur haute = moins bien »). Les noms Index Mellifera sont
  ceux de `mellifera/traits.csv` ; un test vérifie que chaque nom cité existe dans ce fichier et que
  chaque nom du fichier est cité au plus une fois ;
- le moteur ne lit que la colonne Index Mellifera. Les colonnes BeeBreed et Karl Kehrle sont
  documentaires (description de l'énumération dans l'OpenAPI, `METHODOLOGIE.md`), pour les exports futurs
  de Beekube ;
- les caractères de type `harvest` sont dans le catalogue pour Beekube, mais leurs notes Index Mellifera
  (kilos) **ne sont pas** chargées : le contrat impose des valeurs dans l'échelle déclarée, et la mise à
  l'échelle des récoltes reste une décision statistique séparée (hors périmètre, comme au volet 1) ;
- le sens des notes Index Mellifera (note haute = mieux, y compris pour l'essaimage) est celui de
  `traits.csv` (`direction = 1`), jamais vérifié sur le site : point ouvert hérité de la spec précédente
  (§11), à confirmer avant la mise en production.

## 3. Contrat (version 2.2.0)

- `CriterionIn.trait : Optional[Literal[<identifiants du catalogue>]] = None`. Description OpenAPI :
  « caractère standard du critère ; sert au moteur à rattacher les notes Index Mellifera des ancêtres
  résolus, et à Beekube pour ses exports ». L'énumération est régénérée dans `openapi/openapi.json`.
- Validation de `BlupV2Input` : deux critères qui déclarent le même `trait` sont refusés (« critères N et
  M déclarent le même caractère standard `x` »), pour ne jamais compter deux fois une note. Aucune
  contrainte entre `kind` du critère et type du caractère : chaque enregistrement porte son échelle, et le
  moteur sait déjà mélanger des échelles sur un critère (`scale_mixed`).
- Aucune option nouvelle. `external_pedigree_depth` borne l'ascendance, donc les notes ; `0` désactive
  tout ; `MELLIFERA_DB` absente : rien ne change et aucun avertissement (écart 1 du volet 2, confirmé).
- Avertissements nouveaux :
  - `external_records_added` : « N note(s) Index Mellifera ajoutée(s) sur K critère(s) : A ancêtre(s) et
    C contemporaine(s) », suivi des notes écartées par motif quand il y en a, dans cet ordre (« ; écartées :
    a en kg, b sans date, c hors échelle, d non numériques, e sans rucher, f antérieures à la naissance,
    g sans fiche de reine, h sans contemporaine »), puis « ; plafond de 1000 contemporaines atteint » le cas
    échéant ;
  - `external_records_unavailable` : « notes Index Mellifera illisibles : calcul avec l'ascendance mais
    sans leurs notes » ; émis si la lecture des notes échoue alors que l'ascendance a été résolue.
- Sorties inchangées en forme. Les statistiques par critère (`n_records`, `n_queens`,
  `n_contemporary_groups`, par reine `n_records`, `n_evaluators`, `n_years`) **ne comptent que
  l'exploitation** : les notes externes sont comptées dans l'avertissement, pas dans les sorties.
- `ENGINE_VERSION`, `app.py`, `METHODOLOGIE.md` (§ contrat et § lignées externes), `README.fr.md`,
  `README.md` passent à 2.2.0. Le SDK PHP est à régénérer pour que Beekube puisse envoyer `trait` ; sans
  lui, rien ne casse (le champ est optionnel).

## 4. Lecture des notes dans le moteur

Nouveau module `blup2/external_records.py` (lecture et nettoyage), appelé par `blup2/external.py` après la
résolution de l'ascendance (volet 2) :

1. **Condition** : au moins un critère de l'entrée a un `trait` dont le catalogue donne un nom Index
   Mellifera **et** au moins une note de l'exploitation (un critère sans note de l'exploitation ne reçoit
   rien : les notes externes n'y classeraient que des reines invisibles), et au moins un ancêtre a été
   résolu. Sinon rien n'est lu.
2. **Requêtes** (lecture seule, même connexion que l'ascendance), sur `evaluation` jointe au critère du
   site (`criterion (mode, code)`), au caractère (`trait`) et, en jointure externe, à la fiche de la reine
   (`queen.breeder_id`, `queen.birth_year`), limitées aux noms Index Mellifera demandés et aux sources
   `season_table` et `radar` : (a) les notes des codes des ancêtres résolus ; (b) les notes des groupes de
   ces ancêtres (testeurs pour `season_table`, éleveurs pour `radar`), filtrées ensuite sur l'année et le
   critère. Le moteur garde sa propre lecture SQL et n'importe rien du paquet `mellifera/` (contrainte
   d'architecture du volet 2) ; le test de contrat avec la vraie base couvre ces tables.
3. **Nettoyage**, mêmes règles que l'export du volet 1 (`mellifera/export.py::_records`), chaque motif
   compté : type `harvest` écarté (« en kg ») ; `value` non numérique (« non numériques ») ; `eval_date`
   vide, non ISO ou invalide (« sans date ») ; `percentage` divisé par 100 ; échelle invalide ou valeur hors
   de l'échelle (« hors échelle ») ; rucher = testeur (`tester_breeder_id`) pour `season_table`, éleveur de
   la reine (`queen.breeder_id`) pour `radar`, sinon « sans rucher » ; note antérieure à l'année de
   naissance (« antérieures à la naissance ») ; contemporaine sans code ou d'année de naissance invalide
   (« sans fiche de reine ») ; `radar` ignoré quand `season_table` couvre la même reine, le même caractère
   et la même année ; plusieurs valeurs pour une même clé (reine, critère, date, source, rucher)
   moyennées. La valeur est ramenée sur [0, 1] avec l'échelle du caractère, puis inversée (`1 − v`) quand
   le sens du critère de l'exploitation (`direction`) diffère de celui du caractère dans la base.
4. **Contemporaines** : une reine non ancêtre notée dans un groupe (critère, source, rucher, année) qui
   contient une note d'ancêtre. Au plus `MAX_CONTEMPORARIES = 1000` par calcul, retenues par code
   croissant. Un groupe n'est gardé que s'il compte au moins deux reines dont un ancêtre ; les notes
   d'ancêtres d'un groupe écarté sont comptées « sans contemporaine ». Les contemporaines reçoivent les id
   négatifs qui suivent ceux des ancêtres (par code croissant) et sont construites comme eux (`_ancestor` :
   liées à leur mère et à leur côté mâle quand ils sont dans l'ensemble, cycles coupés).
5. **Identifiants** : rucher et évaluateur d'une note externe = un même entier négatif par (source,
   testeur ou éleveur), attribués dans l'ordre trié sous `min(−1, plus petit rucher et plus petit
   évaluateur de l'exploitation)` (−1 est la valeur de l'évaluateur inconnu dans `records.py`). Ils ne
   peuvent entrer en collision ni entre eux ni avec l'exploitation.
6. **Enregistrements** : un `RecordIn` par clé, `queen` = id négatif de la reine, `criterion` = id du
   critère de l'exploitation qui porte le caractère, `value` moyennée et bornée à [0, 1], échelle 0–1,
   `date` = `eval_date`, `apiary` = `evaluator` = l'entier du point 5, `hive_class = "hive"`. Ces
   enregistrements ne passent pas par la validation de `BlupV2Input` (leurs reines n'y sont pas).
7. **Résultat** : `Expansion` gagne `contemporaries: list[QueenIn]` et `records: list[RecordIn]` ; le
   message `external_records_added` est composé ici.

## 5. Ajustement et isolement

- `prepare_criterion(inp, criterion_id, f_workers, extra_queens=(), extra_records=())` : `born` et
  `f_workers` s'étendent aux reines externes (`f_workers` est calculé sur toutes les reines du pedigree) ;
  les enregistrements externes du critère suivent ceux de l'exploitation. **Aucun effet fixe partagé** :
  rucher × année et évaluateur sont séparés par les identifiants négatifs (§4.5), la classe de ruche et la
  classe d'âge des notes externes reçoivent des niveaux propres (préfixe `x`), la consanguinité des
  ouvrières des notes externes est une covariable propre (`f_workers_ext`). Les notes externes ne se
  connectent à l'exploitation que par la parenté.
- `CriterionData` gagne `n_farm_obs` et `farm_X` (plan des effets fixes de l'exploitation seule). Les
  statistiques du critère et par reine, les a priori de variance (variance phénotypique, ā_WW et, en
  breeder, ā_QQ et ā_WQ), le choix de lancer le REML (`should_estimate`), `few_records`, `n_queens` et
  `known_mating_share` portent sur l'exploitation seule. Quand le REML tourne, il voit toutes les données
  du critère, notes externes comprises.
- **Méthode de fiabilité** : elle se décide sur les équations de l'exploitation et de son ascendance. Les
  équations dues aux seules notes externes (effets fixes propres, environnement permanent des reines
  externes notées, nœuds des contemporaines, comptés par une seconde construction du pedigree sans
  elles) ne comptent pas dans `EXACT_LIMIT`. Si le total dépasse la limite à cause d'elles, les PEV des
  reines de l'exploitation restent exactes, lues par résolutions du système déjà factorisé
  (`mme.inverse_entries`), sans inverser la matrice entière. Sans cette règle, le prototype passait en
  fiabilité approchée sur la base de 9 000 reines (8 647 équations), avec un écart allant jusqu'à 0,35
  pour une reine.
- Isolement (spec précédente §8.3) conservé et testé : `queen_outputs`, `compute_indexes` (base 100),
  `suggest_matings`, `evaluate_pairs`, `genetic_trend`, `connectivity` ne voient ni les ancêtres ni leurs
  notes. `connectivity` continue de ne recevoir que `inp`.
- Repli du volet 2 (pedigree invalide avec les ancêtres → calcul sans eux) : il écarte aussi leurs notes.

## 6. Erreurs et journalisation

- Toute erreur pendant la lecture ou la construction des notes (SQL, type inattendu, bogue) :
  `logger.exception`, avertissement `external_records_unavailable`, calcul avec l'ascendance mais sans
  notes. Jamais de 500.
- Un `trait` sans nom Index Mellifera (caractère BeeBreed seul), dont le nom n'existe pas dans la base
  (vieille base, `traits.csv` plus récent), ou dont le critère n'a aucune note de l'exploitation : rien
  d'ajouté pour ce critère, sans avertissement.
- Tables de notes absentes (base du volet 2) : `external_records_unavailable` seulement si un critère
  déclare un caractère lisible ; sinon aucun avertissement.
- Biais de sélection : les notes publiques sont celles des reproductrices, donc des meilleures reines de
  chaque testeur. Les contemporaines remettent dans le calcul les reines du même groupe, donc une bonne
  part de celles parmi lesquelles le choix a été fait ; le reste du biais passe par la parenté. Mesure du
  §7 (c) documentée dans `METHODOLOGIE.md`.

## 7. Validation statistique (porte avant fusion)

Simulation sur le banc de tests lents (`pytest -m slow`, `tests/simdata.py`) :
- une population « externe » notée (plusieurs testeurs × années), une exploitation dont une partie des
  reines descend de reines externes ;
- (a) avec variances fixées : les filles de reines externes notées gagnent en fiabilité par rapport au
  calcul sans notes externes, et la corrélation de leurs EBV avec la vraie valeur ne baisse pas ;
- (b) les reines de l'exploitation **sans lien avec les lignées externes ni groupe commun** avec leurs
  filles gardent exactement leurs EBV (à la tolérance numérique du solveur) ; une reine non apparentée
  mais notée dans le même rucher × année qu'une fille de lignée externe peut bouger, comme dans tout BLUP ;
- (c) biais de sélection : en ne gardant que la meilleure moitié des reines externes notées, la
  corrélation de rang des reines de l'exploitation ne se dégrade pas par rapport au calcul sans notes
  externes.

Seuils fixés d'après le prototype (20 populations, `tests/test_external_calibration.py`), consignés dans
`METHODOLOGIE.md` §10 : (a) gain de fiabilité positif sur chaque population et ≥ 0,05 en moyenne (mesuré
0,068, minimum 0,047), corrélation + 0,02 au moins (mesuré 0,472 → 0,517), méthode de fiabilité exacte
avec et sans notes ; (b) écart ≤ 1e-9 (mesuré 0) ; (c) Spearman avec notes ≥ sans notes (mesuré 0,452 →
0,472, gain de fiabilité 0,033).

Hors de ce volet, après le crawl complet : prendre un éleveur de la base comme exploitation fictive et
le reste comme source externe, et mesurer le gain prédictif par validation croisée avec les scripts de
`docs/analyses/reel/`.

## 8. Tests

- `tests/test_traits.py` : cohérence du catalogue avec `traits.csv`, unicité des correspondances,
  énumération du contrat égale aux clés du catalogue.
- `tests/test_external_records.py` (mini-base de `tests/external_base.py` étendue aux tables
  `trait`, `criterion`, `evaluation`) : rattachement caractère par caractère ; règles de nettoyage par
  motif ; radar ignoré quand la saison couvre ; doublons moyennés ; échelles ; identifiants négatifs sous
  ceux de l'exploitation ; caractère inconnu de la base ; critère sans `trait` ; message de
  `external_records_added`.
- `tests/test_schemas.py` : `trait` accepté, valeur hors catalogue refusée,
  deux critères sur un même caractère refusés.
- `tests/test_external_engine.py` : sorties identiques bit à bit sans `trait`, sans référence ou sans
  base ; fiabilité en hausse pour une fille de reine externe notée ; isolement des sorties, de la base
  100, des accouplements, de la tendance et de la connectivité ; statistiques du critère inchangées ;
  repli `external_records_unavailable` quand la lecture des notes échoue après l'ascendance.
- `tests/mellifera/test_mellifera_external_contract.py` : la requête du moteur fonctionne sur une base
  créée par `mellifera/db.py` et peuplée par l'ingestion des fixtures (fiches avec radar, page de saison,
  une contemporaine).
- `tests/test_records.py` : niveaux propres aux notes externes, notes ignorées sans note de
  l'exploitation, covariable externe. `tests/test_mme.py` : `inverse_entries` égale l'inverse dense.
  `tests/test_external_engine.py` : la méthode de fiabilité ne bascule pas à cause des notes externes
  (profils `standard` et `breeder`).
- Snapshot OpenAPI régénéré et comparé ; version 2.2.0 partout.
- Simulation lente du §7 et mesure de coût du §9, marquées `slow` (`tests/test_external_calibration.py`).

## 9. Points à vérifier à l'implémentation

- `fit.n_records` et les autres compteurs de `CriterionFit` : s'assurer qu'ils sont bien calculés sur
  les enregistrements de l'exploitation seulement, ou filtrés avant les sorties.
- Classe d'âge d'un ancêtre : `age_class(année d'évaluation, born)` avec `born` de l'ancêtre ; vérifier
  qu'une note antérieure à la naissance (erreur de saisie du site) est écartée et comptée.
- Collisions d'identifiants de ruchers négatifs si Beekube en envoyait : couverts par la règle du §4.4.
- Coût, mesuré sur le prototype (base synthétique de 9 000 reines, 300 reines de l'exploitation citant
  100 lignées) : 1 789 notes ajoutées (789 ancêtres, 1 000 contemporaines, plafond atteint), calcul exact
  en 0,45 à 0,5 s contre 0,7 à 0,8 s sans notes. Le test lent borne à 30 s et exige la méthode exacte.
- Sens des notes Index Mellifera (§2) avant toute activation en production.

## 10. Hors périmètre

Comparaison des reines de l'exploitation aux lignées externes sur une même échelle (sortie des reines
externes) ; notes des parentes non ancêtres (sœurs, filles) hors des groupes des ancêtres ; mise à
l'échelle des récoltes en kg ;
injection de valeurs d'élevage BeeBreed ou Karl Kehrle ; groupes génétiques ; crawl ou demande de données
BeeBreed et Karl Kehrle (sujets séparés, voir l'étude du 14/10) ; exports Beekube vers ces outils
(Beekube) ; affichage des notes externes dans Beekube.
