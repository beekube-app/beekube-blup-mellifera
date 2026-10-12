# Analyse de la littérature — fiabiliser le BLUP Beekube

Date : 11/10/2026. Sources analysées (dans `docs/`) :

- **[BB]** Brascamp E.W., Bijma P. (2014). *Methods to estimate breeding values in honey bees*.
  Genetics Selection Evolution 46:53. — `12711_2014_Article_53.pdf`
- **[BER08]** Bienefeld K., Ehrhardt K., Reinhardt F. (2008). *Noticeable success in honey bee
  selection after the introduction of genetic evaluation using BLUP*. American Bee Journal, août
  2008, p. 739–742 (trad. K. S. Traynor). — `BienefeldEhrhardtandReinhardtarticlefromAugust08ABJ.pdf`

Méthode : quatre analyses indépendantes (deux par article, une revue croisée code ↔ littérature, une
analyse numérique en Python), puis vérification dans le code des points retenus. Les numéros de page
renvoient à la pagination du PDF (BB : p. 1–15 ; BER08 : p. 1–4 = p. 739–742 de la revue).
« blup.r:N » renvoie au code au commit `bfe8634`.

---

## 1. Résumé

1. **L'index publié n'est pas encore un BLUP génétique.** Les deux articles décrivent un *modèle
   animal* où la parenté (A⁻¹) entre dans les équations du modèle mixte, de sorte que toutes les
   apparentées (mère, sœurs, filles) contribuent à l'index (BER08 p. 2 ; BB p. 4, éq. 4 et 8).
   Dans `blup.r`, la matrice A ne sert que de **poids** via sa diagonale (blup.r:400) : l'index est
   un effet aléatoire reine indépendant (modèle de répétabilité). Déjà l'item 12 de la feuille de
   route — les articles le confirment et le précisent.
2. **Le poids `diag(A)` est à contresens.** Une reine consanguine (1+F > 1) reçoit *plus* de poids,
   alors que BER08 (p. 1–2) corrige la consanguinité de la reine et des ouvrières comme un **effet**
   défavorable. Sans consanguinité, tous les poids valent 1,0208 et n'ont aucun effet. → à supprimer
   (nouveau P0).
3. **Les constantes de parenté « abeille » de `calculate_bee_kinship` ne correspondent à aucune
   équation de BB** : facteur global `A × (1 + 1/48)` (blup.r:233) et sœurs fixées à
   `0,25 + 0,5/12 = 0,2917` (blup.r:245). La parenté exacte entre pleines sœurs vaut
   ¼ + 1/(2D) + 1/(4S) (BB éq. 24a), soit **0,3229** pour S = 8 et D = 12, et **0,5208** si la
   productrice de mâles est une reine unique (S = 1). Ces constantes sont de toute façon sans effet
   aujourd'hui (seule la diagonale est lue), mais elles coûtent O(n²) (39 s pour 1 200 reines).
4. **La valeur à publier est celle d'un accouplement**, pas d'une reine seule : valeur attendue d'une
   fille moyenne = effet ouvrières + effet reine du groupe d'ouvrières (BB éq. 2). Les deux
   composantes (reine / ouvrières) sont publiées séparément en Allemagne et sont **génétiquement
   corrélées négativement** (BER08 p. 2 ; BB p. 10 : rG = −0,5 en simulation).
5. **Présentation** : en Allemagne, les index sont publiés en % (100 = moyenne des reines des
   5 dernières années, base mobile) depuis 1997, précisément parce que les unités brutes étaient
   illisibles (BER08 p. 2). Les index Beekube sont en unités brutes, centrés sur les reines de
   l'appel, donc non comparables d'un calcul à l'autre.
6. **Plusieurs bugs nouveaux** ont été trouvés et vérifiés (section 4) : paramètres publiés issus d'un
   autre modèle que les index, étiquette `methods` réduite à la première reine, accouplements filtrés
   sur le seul premier critère, fichiers d'échange partagés entre requêtes, deux échelles de notes
   mélangées dans un même critère.
7. **La simulation chiffre le gain** (section 5.2). Le modèle actuel atteint une précision de 0,515.
   Avec la parenté maternelle déjà disponible, elle monte à 0,564. Avec une parenté abeille complète,
   accouplements connus, elle atteint 0,635 et la réponse à la sélection gagne 25 %. Le modèle donne
   aussi un index aux reines sans note : sur l'exemple, 180 des 259 reines non notées ont une
   apparentée notée. L'index actuel est mal calibré : pente de 0,41, pour 1 attendu.
8. **Ordre de mise en œuvre conseillé** : mettre A dans le modèle (item 12) d'abord, puis enregistrer
   les accouplements (items 16, 31), puis passer aux deux effets reine / ouvrières (item 15).
   Raffiner la D bloc-diagonale (item 14) ne vient qu'ensuite. Selon BB (p. 11–12), elle n'apporte
   que +4 à 5 % de réponse et classe les reines à peu près comme l'approximation de Bienefeld et
   al. 2007. Dans la simulation, les hypothèses D et S ne changent pas la précision.

---

## 2. Ce que disent les articles

### 2.1 Bienefeld, Ehrhardt & Reinhardt 2008 (vulgarisation, Allemagne)

- **Standardiser les tests est nécessaire mais insuffisant** (Ruttner 1972) : miellée et normes de
  notation du comportement (douceur, tenue au cadre) varient tellement d'un rucher et d'un apiculteur
  à l'autre que la comparaison n'a de sens qu'à l'intérieur d'un rucher (p. 1).
- **Tester des sœurs dans un même rucher** atténue les biais d'environnement, mais l'écart à la
  moyenne du rucher reste biaisé par le **niveau génétique des colonies voisines** (exemple du sprint
  régional vs finale olympique, p. 1). Le BLUP corrige ce biais car il connaît la valeur des
  comparantes (p. 2). C'est l'argument central contre un modèle sans parenté.
- **Modèle animal abeille** (Hohen Neuendorf, 1992–1994) : toute la parenté, deux valeurs d'élevage
  par caractère (reine et ouvrières), parenté entre reines, entre ouvrières moyennes et entre reines
  et ouvrières ; base Willham 1963, adaptation Bienefeld et al. 2007 (p. 2).
- **Corrections préalables** : consanguinité de la reine et des ouvrières (la performance est corrigée
  « comme si F = 0 » ; Bienefeld et al. 1989) et nombre variable de colonies comparées par rucher
  (formule non détaillée) (p. 1–2).
- **Corrélation génétique négative** entre effets reine et ouvrières (Bienefeld & Pirchner) : l'ignorer
  a fait échouer la sélection dans d'autres espèces (p. 2).
- **Définition de l'index** : écart à la moyenne de la population en %, la descendance attendue étant
  la moyenne des parents (reine 120 % × mâles moyens 100 % → 110 %) (p. 1).
- **Base 100 mobile sur 5 ans**, critère par critère, depuis 1997 ; les anciens index « baissent »
  quand la population progresse (p. 2, fig. 3–4).
- **Résultats** : progrès génétique sur le miel de +0,11 %/an avant 1994 à +0,59 %/an après le BLUP
  (+0,51 → +0,72 kg/colonie/an) ; douceur de +0,01 à +0,41 %/an (la fig. 4 indique 0,44) (p. 3).
- **Piège** : la période 1989–1994 (réunification, peu de données : 4 058 colonies contre 8 000 à
  16 000 ailleurs) crée une fausse chute de la tendance génétique ; l'inclure ramène le progrès
  d'avant 1994 à 0,003 % (p. 2–3, fig. 2). → Les tendances calculées sur peu de colonies sont
  trompeuses.
- **Transparence** : beebreed.eu publie index et consanguinité (reine, ouvrières) et propose un
  planificateur d'accouplement donnant l'index et la consanguinité attendus de la descendance (p. 3).
- **Ce que l'article ne contient pas** : aucune valeur d'héritabilité, aucun nombre minimal de sœurs,
  rien sur la fiabilité des index, ni sur la structure exacte des effets d'environnement.

### 2.2 Brascamp & Bijma 2014 (méthode)

- **Trois types de nœuds dans le pedigree** (p. 3) : reines ; « sires » = groupe des S reines
  productrices de mâles (RPM), en pratique des sœurs, traité comme un seul individu ; groupes
  d'ouvrières (la colonie). Chaque nœud a deux parents (une reine, un sire). Les mâles d'une RPM ne
  portent que les gènes de la RPM (haploïdie).
- **Modèle** (éq. 1–5) : y = Xb + Z_W a_W + Z_Q a_Q + e ; effet ouvrières porté par le groupe
  d'ouvrières, effet reine par la reine ; même A pour les deux effets ; covariance génétique
  rG·σ_AW·σ_AQ. Critère de sélection d'une reine fécondée : Â_W + Â_Q de son groupe d'ouvrières
  (éq. 2) — on sélectionne un **accouplement**.
- **Inverse de A** : A⁻¹ = (I−M)ᵀ D⁻¹ (I−M) (éq. 8), comme chez les autres espèces, mais **D n'est
  pas diagonale** : deux pleines sœurs peuvent porter le même gamète paternel (un mâle = un clone de
  spermatozoïdes). D est bloc-diagonale, un bloc par couple reine × sire, facile à inverser
  (p. 5–8).
- **Formules clés** : variance de ségrégation d'une reine (éq. 13) ; d'un sire (éq. 17) ; d'un
  groupe d'ouvrières = cov(δ_FS) (éq. 19) ; sous hypothèse de Poisson cov(δ_FS) = (2 − a_ss)/(4D)·σ²_A
  (éq. 21a) ; parenté entre RPM sœurs a_ss (éq. 24a) ; consanguinité F = ½ a(mère, sire) (éq. 25).
  Six cas de parents inconnus (annexe 1) : un sire inconnu n'est **pas** un père fantôme commun,
  la fille a simplement var(δ) = ¾ − ¼F_mère.
- **Paramètres** : D = 12 mâles par reine (valeur Bienefeld), consensus ≈ 16, mais utiliser un
  **nombre effectif plus petit** car les contributions ne suivent pas Poisson ; S = 8 ; variances
  de simulation σ²_AW = 1, σ²_AQ = 0,5, σ²_E = 2, rG = −0,5 (p. 10, 12).
- **Résultats** (tables 2–3) : la méthode exacte est sans biais (régression valeur vraie / index ≈ 1).
  L'approximation de Bienefeld et al. 2007 surestime la variance des colonies (régression **0,42**
  sur l'effet ouvrières des colonies), mais classe à peu près pareil : corrélation 0,308 contre
  0,278, réponse +4 à 5 %. Le gain principal est pour **l'estimation des paramètres génétiques**.
- **Si les RPM ne sont pas des sœurs** (hors Beebreed), les mettre individuellement dans le pedigree
  avec des fractions 1/(2S) (p. 12) — cas le plus probable chez les utilisateurs Beekube.
- **Limites pour Beekube** : structure Beebreed supposée (stations, IA, 6 000 colonies/an) ;
  fécondation libre = sire inconnu, sans modélisation ; variances supposées connues ; une
  observation par colonie, un seul effet fixe, un seul caractère, pas de notes ordinales ni
  d'évaluateurs.

### 2.3 Référence « Bienefeld et al. 2007 »

La référence de la feuille de route est correcte et précisée par BER08 (p. 3–4) :
Bienefeld K., Ehrhardt K., Reinhardt F. (2007). *Genetic evaluation in the honey bee considering
queen and worker effects — a BLUP-animal model approach*. Apidologie 38:77–85.
Elle porte sur le **modèle reine / ouvrières** (item 15) ; les **groupes de reines productrices de
mâles** et la D bloc-diagonale relèvent de Brascamp & Bijma 2014 (item 14).

---

## 3. Code actuel face à la littérature

| Étape | Code actuel | Littérature | Écart / risque |
|---|---|---|---|
| Pedigree | `sire = drone_parent` traité comme un mâle diploïde (blup.r:322–354) | Sire = groupe de S RPM ; 3 types de nœuds (BB p. 3) | Sens de `drone_parent` non défini : RPM elle-même ou mère du groupe ? Père de la reine ou de ses ouvrières ? |
| Parenté | `getA` × (1 + 1/48) ; sœurs forcées à 0,2917 (blup.r:230–258) | D bloc-diagonale, éq. 13–24a (BB p. 4–8) | Constantes sans fondement ; écrasement a posteriori qui ne se propage pas aux descendantes ; O(n²) |
| Poids | `weights = diag(A_subset)` (blup.r:400) | A⁻¹ dans les équations (BB éq. 4) ; F en effet (BER08 p. 1–2) | Contresens : reines consanguines surpondérées ; matrice dense n_notes² (3,2 Go pour 20 000 notes) |
| Modèle | `note ~ (1|queenbee) [+ (1|drone_parent)] + (1|user) + (1|apiary) + (1|beehiveType) + (1|month) + (1|year)` | Effets reine + ouvrières liés à A, rG, environnement en effets fixes (BB éq. 3–5) | Modèle de répétabilité ; notes répétées traitées comme indépendantes |
| Paramètres | Second modèle (blup.r:69–81), h² = (v_reine + v_drone)/total, SE = √(2/n) | σ²_AW, σ²_AQ, rG estimés sur le modèle à deux effets (BB éq. 5) | Les paramètres publiés ne viennent pas du modèle qui produit les index |
| Valeur publiée | `ranef(queenbee) + ranef(drone_parent)` (blup.r:405–417) | Â_W + Â_Q du groupe d'ouvrières (BB éq. 2) ; deux composantes publiées (BER08 p. 2) | Un seul nombre mêlant génétique et environnement permanent ; reines non notées = `null` même avec des filles notées |
| Échelle | Unités brutes, centrées sur l'appel | % base 100, moyenne mobile 5 ans (BER08 p. 2) | Non comparable entre calculs |
| Accouplements | Ratio d'ancêtres communs sur 3 générations ; 0,7 × BLUP + 0,3 × diversité | Descendance = ½(mère + mâles) ; F = ½ a(mère, sire) (BER08 p. 1, 3 ; BB éq. 25) | La diversité pilote le classement (terme BLUP ≈ 0,2 d'amplitude contre 0,3) |

---

## 4. Bugs nouveaux, vérifiés dans le code

L'agent d'audit a exécuté `blup.r` (R 4.6.1) sur les données d'exemple converties au format de
l'API ; les mesures chiffrées en viennent (temps de calcul, 4 reines exclues, échec d'un critère).
J'ai confirmé les causes en relisant le code pour B2, B3, B5, B6 et pour les poids `diag(A)`.

| # | Constat | Lieu | Effet |
|---|---|---|---|
| B1 | Les variances et l'« héritabilité » publiées sont estimées par un **autre modèle** que celui des index (moins d'effets, pas de poids) ; h² additionne une variance entre reines et une variance entre groupes paternels sans facteur ×4 | blup.r:69–99 | Les paramètres affichés ne décrivent pas l'index publié ; ni h² ni répétabilité |
| B2 | `methods[[critere]]` reçoit un scalaire (`"Not enough data"`) puis est lu à l'indice `[i]` | blup.r:372 (aussi 435, 468), 642 | Seule la première reine porte l'étiquette ; une reine non notée et un échec renvoient tous deux `null` |
| B3 | Les suggestions d'accouplement ne retiennent que les reines ayant un BLUP sur le **premier** critère | blup.r:553, 563 | Sur l'exemple, 4 des 69 reines indexées sont exclues comme reine et comme partenaire |
| B4 | `check.nobs.vs.nlev = "ignore"` n'est passé qu'au modèle d'héritabilité, pas au `lmer` des index | blup.r:78–81 vs 400 | Un critère avec une note par reine échoue ; sa clé disparaît du JSON (complète l'item 1) |
| B5 | Boucle sur toutes les paires de reines et matrice dense n_notes × n_notes, pour un résultat dont seule la diagonale est lue | blup.r:241–258, 397 | 0,8 s pour 300 reines, 39 s pour 1 200 (au-delà du timeout gunicorn de 30 s par défaut) |
| B6 | Fichiers d'échange fixes `data/input.json` / `resultats_index.json` ; code retour de Rscript ignoré ; avertissements R perdus | app.py:125–157, blup.r:19–20 | Avec plusieurs workers, une exploitation peut recevoir les résultats d'une autre ; échec R renvoyé en HTTP 200 |
| B7 | `data/input-example.json` n'a pas le schéma de l'API (rejeté par pydantic) ; un champ `apiary` dans la reine produit `apiary.x` / `apiary.y` au `left_join` et l'effet rucher disparaît | data/, blup.r:367, 390 | Pas de jeu de test utilisable ; perte silencieuse de l'effet rucher |
| B8 | `results: Union[BLUPResultOutput, dict]` : la validation de sortie retombe sur `dict` dès qu'une valeur ne valide pas | app.py:84 | Sortie non validée |
| B9 | `brms` masque `ranef` / `ngrps` ; l'extraction dépend de l'ordre des `library()` | blup.r:1–8, 405 | Fragilité (complète l'item 6) |

### Données d'exemple (`data/input-example.json`)

- 328 reines, dont 69 notées ; toutes les notes datent de 2024 ; 5 évaluateurs, 6 ruchers, 3 types
  de ruche ; `drone_parent` absent pour 100 % des reines → **la branche « abeille » du code n'a
  jamais tourné sur ces données**.
- **Deux échelles de notes coexistent dans un même critère** : {0,25 ; 0,5 ; 0,75 ; 1} de mars à mai,
  {0,2 ; 0,4 ; 0,6 ; 0,8 ; 1} de juin à septembre, pour presque tous les évaluateurs (critères 44,
  45, 46). Le changement d'échelle est confondu avec l'effet mois. Aucune note ne vaut 0.
- **Ruchers non connectés** : l'évaluateur 848 note seul trois ruchers, les quatre autres évaluateurs
  notent les trois autres ; une seule reine change de rucher. Les deux blocs ne sont reliés que par le
  pedigree, que le modèle n'utilise pas : leurs index ne sont pas comparables.
- **Information perdue** : 23 mères non notées ont des filles notées (ex. une reine avec 31 filles dont
  6 notées) et reçoivent `null` ; inversement, des dizaines de filles de mères notées n'ont aucun index.
  Un vrai modèle animal leur donnerait un index (avec une fiabilité).
- 232 combinaisons reine × critère × mois ont plusieurs notes, traitées comme indépendantes
  (pseudo-réplication).
- Les BLUP vont de −0,13 à +0,14 : `blup_rescaled` reste entre 0,37 et 0,64, donc la diversité
  pèse davantage que l'index dans le score d'accouplement.

---

## 5. Analyses numériques

### 5.1 Parenté abeille calculée selon Brascamp & Bijma

Le papier ne donne aucune matrice A chiffrée. Les valeurs ci-dessous sont dérivées de ses équations
(calcul exact en fractions, vérifications : A = (I−M)⁻¹ D (I−M)ᵀ⁻¹, A·A⁻¹ = I, diagonale = 1 + F,
sœurs = éq. 24a).

| Relation | blup.r actuel | `getA` standard | BB, S = 8, D = 12 | BB, S = 8, D = 16 | BB, S = 1, D = 12 |
|---|---|---|---|---|---|
| Pleines sœurs (même reine, même sire) | 0,2917 | 0,5 | **0,3229** | 0,3125 | **0,5208** |
| Reine – son groupe d'ouvrières | — | — | 0,5 | 0,5 | 0,5 |
| Variance du groupe d'ouvrières (diag. A) | — | — | 0,3229 | 0,3125 | 0,5208 |
| Sire de base (diag. A) | 1,0208 | 1 | 0,125 = 1/S | 0,125 | 1 |
| Reine de base (diag. A) | 1,0208 | 1 | 1 | 1 | 1 |

Lecture :
- 0,2917 = ¼ + 1/(2D) est la limite de l'éq. 24a quand le nombre de RPM tend vers l'infini : elle
  sous-estime la parenté réelle, surtout quand la productrice de mâles est une reine unique (cas
  probable de `drone_parent` aujourd'hui). Le `getA` standard (0,5) est même plus proche du vrai que
  la valeur « corrigée » dans ce cas.
- Le facteur global (1 + 1/48) gonfle aussi les reines de base, ce qu'aucune équation ne justifie.
- **Point d'implémentation découvert en vérifiant** : avec une D non diagonale, le calcul ligne par
  ligne de A exige un tri **par familles de pleines sœurs** (blocs consécutifs, en ordre
  topologique). Un simple tri « parents avant enfants » donne une A fausse quand une nièce précède sa
  tante (écart mesuré 1/48 sur l'exemple). Le calcul de A⁻¹ par l'éq. 8 n'est pas concerné. La
  contrainte disparaît si l'on ajoute à la récurrence le terme Σ_k D_ik T_jk (T = (I−M)⁻¹), ou si
  l'on calcule A = (I−M)⁻¹ D (I−M)⁻ᵀ directement. C'est ce que fait `docs/analyses/bee_kinship.py`,
  dont le résultat ne dépend pas de l'ordre (vérifié par `test_order.py`).

Comparaison sur un pedigree de type Beebreed (S = 8, D = 12) des coefficients produits par chaque
méthode :

| Couple | A abeille S = 8 | A abeille S = 1 | blup.r | `getA` |
|---|---|---|---|---|
| Mère – fille | 0,5000 | 0,5000 | 0,5104 | 0,5000 |
| Pleines sœurs, parents apparentés | 0,4025 | 0,6458 | 0,2917 | 0,6250 |
| Demi-sœurs maternelles (sires différents) | 0,2500 | 0,2500 | **0,2917** | 0,2500 |
| Demi-sœurs paternelles (même sire) | 0,0312 | 0,2500 | 0,2552 | 0,2500 |
| Reine – son sire | 0,0625 | 0,5000 | 0,5104 | 0,5000 |
| Plus petite valeur propre de A | 0,0103 | > 0 | **−0,205** | > 0 |

- blup.r attribue 0,2917 à **toutes** les reines de même mère, y compris des demi-sœurs de sires
  différents (le test blup.r:244 ne compare que les mères).
- **La matrice de blup.r n'est pas définie positive** (valeur propre −0,205) : même si on la
  branchait comme covariance dans le modèle, celui-ci serait invalide. Elle ne peut pas servir de
  base à l'item 12 : il faut repartir de `getA` (ou de la construction de BB).
- Selon que `drone_parent` désigne une RPM unique (S = 1) ou un groupe de 8 sœurs, les parentés
  paternelles varient d'un facteur 8.
- Sur `data/input-example.json` (aucun `drone_parent`), la matrice de blup.r est identique à `getA`
  « mère seule » et tous les poids valent 1 : **la parenté n'a aucun effet sur les index publiés**.
  Pourtant, **180 des 259 reines sans note ont au moins une apparentée notée** et pourraient recevoir
  un index avec un modèle animal.

### 5.2 Simulation : précision des différents modèles

Protocole : pedigree simulé sur 5 générations (≈ 270 reines, 30 sires de 8 RPM, D = 12), vraies
valeurs tirées selon la covariance abeille de BB ; 70 % des reines notées, 2 à 5 notes, effets
rucher et environnement permanent ; h² par note ≈ 0,11, répétabilité ≈ 0,28. Équations du modèle
mixte résolues à variances connues, 200 répétitions (erreur-type ≈ 0,005).

| Modèle | Précision, reines notées | Pente (calibration, idéal 1) | Précision, reines sans note | Valeur vraie des 10 % meilleures (σ_A) |
|---|---|---|---|---|
| (a) BLUP avec A abeille, accouplements connus | **0,635** | **1,00** | **0,501** | **0,665** |
| (a') A abeille, accouplements inconnus (« sire fantôme ») | 0,564 | 1,08 | 0,325 | 0,582 |
| (b) Reine iid sans parenté = **blup.r actuel** | 0,515 | 0,41 | aucun index | 0,530 |
| (b2) blup.r avec `drone_parent` iid | 0,521 | 0,42 | 0,138 | 0,536 |
| (c) Modèle animal diploïde, `drone_parent` = père | 0,554 | 0,93 | 0,311 | 0,572 |
| (d0) Modèle animal diploïde, mère seule | 0,564 | 0,92 | 0,325 | 0,583 |

Précision = corrélation entre l'index et la valeur génétique vraie de la colonie.

Scénario à héritabilité plus faible (σ²_e = 4, 100 répétitions) : (a) 0,545, (b) 0,413, (c) 0,457,
(d0) 0,467. La réponse à la sélection des 10 % meilleures passe de 0,464 (b) à 0,618 (a), soit +33 %.

**Sensibilité aux hypothèses** (vrais S = 8, D = 12) : supposer D entre 4 et 48 ou S entre 4 et 100
ne change pas la précision (0,635 ± 0,001) ; seule la pente varie (0,92 à 1,04). Supposer S = 1 à
tort coûte 0,02 de précision.

Conclusions :

1. **Mettre la parenté dans le modèle est le gain principal et il est disponible tout de suite** :
   avec le seul pedigree maternel déjà transmis (d0), la précision passe de 0,515 à 0,564 (+10 %),
   l'index est correctement calibré (pente 0,92 au lieu de 0,41) et les reines sans note en
   reçoivent un. Rien ne bloque cette étape côté données.
2. **Connaître l'accouplement de chaque reine** apporte le reste du gain : 0,564 → 0,635 sur les
   reines notées, 0,325 → 0,501 sur les reines sans note. C'est l'enjeu des items 16 et 31 côté
   Beekube.
3. **Traiter `drone_parent` comme un père diploïde (c) ne fait pas mieux que l'ignorer (d0).** La
   structure abeille (sire = groupe, groupe d'ouvrières) compte quand l'accouplement est connu.
4. **L'index actuel est surdispersé** (pente 0,41) : il mêle génétique et environnement permanent.
   Un écart de 1 sur l'index correspond en moyenne à 0,41 de valeur génétique.
5. **D et S ont peu d'importance pratique** pour le classement : prendre D = 12 à 16 et S = 4 à 8
   par défaut. L'item 36 est un réglage de calibration, pas une priorité.

Limites : un seul caractère (effet ouvrières), variances connues (pas d'estimation REML), pedigree
sans erreur, accouplements aléatoires, hypothèse de Poisson. Les modèles (b), (c) et (d0) reçoivent
la vraie variance des colonies, ce qui les avantage par rapport à une estimation sur petits
effectifs. `getA` et `calculate_bee_kinship` ont été répliqués en Python d'après le code R, sans
être exécutés en R.

### 5.3 Reproduire les analyses

Scripts dans [`docs/analyses/`](analyses/) (Python 3, `numpy`) :

| Script | Rôle |
|---|---|
| `bee_kinship.py` | Construction de A abeille et de A⁻¹ (BB), réplique de `getA` et de `calculate_bee_kinship` |
| `verify_paper.py` | Vérification des équations de BB, des identités matricielles et par tirage de gènes (*gene dropping*) |
| `test_order.py` | Indépendance de A vis-à-vis de l'ordre du pedigree |
| `compare_r.py` | Écarts blup.r / A abeille sur un pedigree type et sur `data/input-example.json` |
| `simulate.py` | Simulation de la section 5.2 (sorties `sim_main.txt`, `sim_main_100.txt`, `sim_lowh2.txt`) |

```bash
cd docs/analyses && python verify_paper.py && python compare_r.py && python simulate.py
```

---

## 6. Données nécessaires côté Beekube

Pour un BLUP abeille, l'API devrait recevoir :

**Par reine**
- `born` complet (au minimum année sur 4 chiffres), date d'introduction dans la colonie.
- Fécondation, en deux objets distincts (BB fig. 1 distingue le père de la reine et le père de ses
  ouvrières) :
  - `sireGroup` : les mâles qui ont fécondé la **mère** de la reine ;
  - `matingGroup` : les mâles qui ont fécondé la **reine elle-même** ;
  - chacun : mère des RPM, liste des RPM si connue, S, station / IA / libre, année, référence externe
    (`droneRefExt`). En fécondation libre : **aucun** sire (et non `"0"`).
- Origine / fournisseur et année d'origine (groupes génétiques).
- Date de mort, élimination et motif, disponibilité comme mère ou productrice de mâles.

**Par note**
- Date (au jour) ou identifiant de visite, au lieu de `month` / `year` sur 2 chiffres.
- Échelle de la note (min, max, pas, ou version de grille) et sens de 0.
- Colonie physique, type de ruche normalisé (ruchette / ruche / divisée), rucher, évaluateur.

**Par critère** (objet au lieu de la liste d'entiers `evaluate`)
- Identifiant, nom, sens (+1 / −1), type (ordinal / continu / binaire), unité, poids dans l'index de
  l'utilisateur, éliminatoire ou non, saison attendue, héritabilité a priori facultative.

**Par calcul** : version des paramètres et années de base, pour que deux calculs soient comparables.

---

## 7. Spécification : construire A et A⁻¹ abeille

À partir des équations de BB (détail dans la section 2.2) :

1. **Nœuds** : une reine Q par `queenbee` ; un sire S par groupe de mâles distinct (type « groupe »,
   avec S RPM) ou, si les RPM sont connues individuellement, pas de nœud sire mais des fractions
   1/(2S) vers chaque RPM (BB p. 12) ; un groupe d'ouvrières W par reine fécondée (sire connu).
2. **Tri** : par familles de pleines sœurs (même mère × même sire), blocs en ordre topologique.
3. **Récurrence** : F_i = ½ a(mère, sire) (éq. 25) ; variance de ségrégation selon le type de nœud
   et les parents connus (éq. 13, 17, 19 ; annexe 1, cas 1 à 6) ; covariance entre pleines sœurs
   cov(δ_FS) (éq. 21 / 21a) dans le bloc ; a_ss du sire par l'éq. 24 ; puis
   A[i,j] = ½ A[mère,j] + ½ A[sire,j] + D[i,j].
4. **Inverse** : A⁻¹ = (I−M)ᵀ · blocdiag(D_B⁻¹) · (I−M), M creuse avec ½ (ou 1/(2S)) sur les
   parents.
5. **Modèle** : effets W et Q sur la même A, G = [[σ²_AW, rG σ_AW σ_AQ], [·, σ²_AQ]] ⊗ A ; valeur
   publiée = Â_W + Â_Q du groupe d'ouvrières (éq. 2), plus les deux composantes ; consanguinité d'un
   accouplement planifié = ½ A[reine, sire] (éq. 25).
6. **Tests d'invariants** : diagonale des reines = 1 + F ; sœurs de parents de base =
   ¼ + 1/(2D) + 1/(4S) ; reine – son groupe d'ouvrières = 0,5 pour des parents non apparentés ;
   A·A⁻¹ = I ; sire inconnu ⇒ var(δ) = ¾ − ¼F_mère ; permutation nièce / tante sans effet.

Outils R possibles (pedigreemm ne sait pas représenter une D non diagonale) : `sommer`, ou BLUPF90 /
ASReml avec A⁻¹ fournie.

---

## 8. Limites

- BER08 est un article de vulgarisation : pas d'équations, pas d'héritabilités, pas de structure
  précise des effets d'environnement. Les chiffres de tendance sont lus en partie sur des figures.
- BB suppose la structure Beebreed (stations, insémination, RPM sœurs, gros effectifs), une
  observation par colonie et des variances connues. Chez les utilisateurs Beekube (fécondation libre
  fréquente, petites exploitations, notes ordinales répétées), son apport est surtout **conceptuel**
  (nœuds, valeur publiée, parents inconnus) ; le gain chiffré de la D exacte (+4 à 5 %) sera noyé dans
  le bruit des petits effectifs.
- Les variances ne seront pas estimables sur quelques dizaines de colonies : utiliser des valeurs
  publiées (BB p. 10 ; Bienefeld & Pirchner 1990 ; Bernstein et al. 2023).

---

## 9. Références

- Bienefeld K., Ehrhardt K., Reinhardt F. (2007). Genetic evaluation in the honey bee considering
  queen and worker effects — a BLUP-animal model approach. *Apidologie* 38:77–85.
- Bienefeld K., Ehrhardt K., Reinhardt F. (2008). Noticeable success in honey bee selection after
  the introduction of genetic evaluation using BLUP. *American Bee Journal*, août 2008, 739–742.
- Bienefeld K., Pirchner F. (1990). Heritabilities for several colony traits in the honeybee
  (*Apis mellifera carnica*). *Apidologie* 21:175–183.
- Bienefeld K., Pirchner F. (1991). Genetic correlations among several colony characters in the
  honey bee taking queen and worker effects into account. *Annals of the Entomological Society of
  America* 84:324–331.
- Bienefeld K., Reinhardt F., Pirchner F. (1989). Inbreeding effects of queen and workers on colony
  traits in the honeybee. *Apidologie* 20:439–450.
- Brascamp E.W., Bijma P. (2014). Methods to estimate breeding values in honey bees. *Genetics
  Selection Evolution* 46:53.
- Ruttner H. (1972). Technische Empfehlung zur Methodik der Leistungsprüfung von Bienenvölkern.
  Apimondia, Bucarest, 103–107.
- Willham R.L. (1963). The covariance between relatives for characters composed of components
  contributed by related individuals. *Biometrics* 19:18–27.
