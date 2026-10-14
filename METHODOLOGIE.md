# Méthodologie du moteur BLUP v2 (Beekube)

Ce document décrit ce que calcule réellement le moteur v2 (`POST /v2/blup`, paquet `blup2/`),
les hypothèses retenues et les limites connues. Il s'adresse aux apiculteurs et sélectionneurs qui
lisent les index dans Beekube, et aux personnes qui relisent le code. Les références entre crochets
renvoient à la section « Références » ; [BB] = Brascamp & Bijma 2014, [BER08] = Bienefeld,
Ehrhardt & Reinhardt 2008. Version du moteur : `2.1.0`, profils `standard` et `breeder`.

## 1. Pourquoi un BLUP avec parenté

Comparer une colonie à la moyenne de son rucher ne suffit pas : l'écart est biaisé par le niveau
génétique des colonies voisines. Une bonne colonie entourée de sœurs excellentes paraît moyenne ;
une colonie moyenne entourée de colonies faibles paraît bonne (l'exemple du sprint régional contre
la finale olympique, [BER08 p. 1]). Le BLUP corrige ce biais parce qu'il connaît, par la parenté,
la valeur génétique des colonies comparées ; il fait aussi profiter chaque reine des notes de ses
apparentées (mère, sœurs, filles) [BER08 p. 2 ; BB p. 4].

En Allemagne, le passage à une évaluation BLUP de ce type a fait passer le progrès génétique sur
la production de miel de **+0,11 % par an à +0,59 % par an** (et la douceur de +0,01 à +0,41 % par
an) [BER08 p. 3].

## 2. Données utilisées

- **Reines** : identifiant, année de naissance, mère, accouplement (type, reines productrices de
  mâles « RPM », mère des RPM, nombre de RPM S, nombre de mâles D, station, année, référence
  externe), statut (morte, éliminée, disponible comme mère ou comme source de mâles).
- **Notes** : reine, critère, valeur et échelle (`scale_min`, `scale_max`), date, évaluateur,
  rucher, classe de ruche (`hive`, `nuc`, `split`), visite.
- **Critères** : sens (+1 si « plus = mieux », −1 sinon), type (`note`, `percentage`, `harvest`),
  poids dans l'index total, caractère éliminatoire, h² et répétabilité facultatives (si une seule
  des deux est fournie, l'autre prend sa valeur par défaut ; répétabilité seule : h² = min(0,25,
  répétabilité)), et, pour le profil `breeder` seulement, `queen_ratio` (ρ) et `rg` (rG)
  facultatifs (§6).
- **Profil** : champ `profile` de la requête, `"standard"` (défaut) ou `"breeder"` (§5, §5 bis).

Préparation des notes, critère par critère :
- une note est ramenée sur [0, 1] : (valeur − `scale_min`) / (`scale_max` − `scale_min`) ; les
  pourcentages et récoltes sont déjà sur [0, 1] ;
- les notes d'une même reine prises lors d'une même visite (même `visit`, même rucher, même année ;
  à défaut même date et même rucher) sont **moyennées** et forment une observation ;
- plusieurs échelles pour un même critère déclenchent l'avertissement `scale_mixed`.

## 3. Pedigree abeille

Le pedigree comporte trois sortes de nœuds [BB p. 3] :
- **Q** : la reine ;
- **S** : le « sire » d'un accouplement, c'est-à-dire le groupe des S reines productrices de mâles
  (RPM, sœurs entre elles) dont les mâles ont fécondé la reine ; avec une seule RPM connue, le sire
  est la RPM elle-même ;
- **W** : le groupe d'ouvrières d'une colonie, dont les parents sont la reine et son accouplement.

Règles :
- la mère d'une reine est sa mère ; **son père est l'accouplement de sa mère** (les filles et les
  ouvrières d'une reine ont le même père). Le moteur v1 décalait le pedigree d'une génération en
  prenant l'accouplement de la reine pour son propre père ;
- fécondation naturelle ou inconnue : **pas de nœud sire** (père inconnu, aucun père fictif
  commun) ; il en va de même d'un accouplement par station ou insémination sans aucune information
  (ni RPM, ni mère des RPM, ni station, ni référence externe) ;
- un parent cité mais absent des données est ajouté comme reine de base (`unknown_parent_added`) ;
- une reine citée comme sa propre mère, sa propre RPM ou la mère de ses RPM (erreur de saisie) :
  ce lien est retiré et l'avertissement `self_reference_removed` liste les reines concernées ;
- un pedigree cyclique (une reine ancêtre d'elle-même par au moins une autre reine) est refusé
  (HTTP 422).

## 4. Parenté

La matrice de parenté A entre tous les nœuds suit [BB] : chaque nœud reçoit la moitié de la valeur
de chacun de ses parents plus un aléa de ségrégation dont la variance dépend de la consanguinité
des parents, de S (nombre de RPM) et de D (nombre de mâles, option `default_n_drones`, 12 par défaut). La covariance
particulière entre pleines sœurs (mêmes mère et sire) est prise en compte [BB éq. 21–24].
A⁻¹ est creuse et construite directement à partir des parents de chaque nœud, sans inverser A.

Propriétés vérifiées par les tests (`tests/test_kinship.py`) :
- identité avec l'implémentation de référence `docs/analyses/bee_kinship.py` (écart < 1e-12) ;
- diagonale d'une reine = 1 + F ;
- pleines sœurs issues de parents de base : ¼ + 1/(2D) + 1/(4S) (0,3229 pour S = 8, D = 12) ;
- A · A⁻¹ = I ; le résultat ne dépend pas de l'ordre des reines dans l'entrée.

Consanguinité publiée : F de la reine (nœud Q) et F des ouvrières (nœud W) = ½ a(reine, sire)
[BB éq. 25].

Père inconnu : la reine reçoit la variance de ségrégation ¾ − ¼ F_mère [BB annexe 1, cas 3] ; le
groupe d'ouvrières reçoit 1/(2D) (moyenne de D gamètes paternels de colonies non apparentées).

## 5. Modèle (profil `standard`)

Pour chaque critère :

```
y = Xb + Z_W a + Z_pe pe + e
a  ~ N(0, A σ²_A)    valeur génétique, portée par le nœud W de la colonie
pe ~ N(0, I σ²_pe)   environnement permanent de la colonie (une colonie = une reine)
e  ~ N(0, I σ²_e)
```

Effets fixes, chacun gardé seulement s'il a au moins deux niveaux avec données :
- groupe de contemporaines **rucher × année** d'évaluation ;
- **évaluateur** (une note sans évaluateur forme son propre niveau) ;
- **classe de ruche** (`hive`, `nuc`, `split`) ;
- **âge de la reine** à l'évaluation (0, 1, 2 ans et plus) ;
- **consanguinité des ouvrières** F_W en covariable (correction « comme si F = 0 »,
  [BER08 p. 1–2]).
Les colonnes redondantes (par exemple un évaluateur qui note un seul rucher) sont retirées.

Résolution : équations du modèle mixte de Henderson, creuses, factorisées par `scipy` (SuperLU).

## 5 bis. Modèle (profil `breeder`)

Profil demandé par `profile = "breeder"`. La note d'une colonie dépend de deux effets génétiques
distincts [BB éq. 1–5 ; Bienefeld et al. 2007] : l'effet des **ouvrières** (comportement du
groupe d'ouvrières, porté par le nœud W de la colonie) et l'effet de la **reine** (ponte,
phéromones, porté par le nœud Q de sa reine).

```
y = Xb + Z_W a_W + Z_Q a_Q + Z_pe pe + e
[a_W ; a_Q] ~ N(0, G0 ⊗ A)    G0 = [ σ²_AW            rG σ_AW σ_AQ ]
                                   [ rG σ_AW σ_AQ     σ²_AQ        ]
```

Effets fixes, environnement permanent et matrice de parenté A sont ceux du profil `standard` ;
chaque nœud porte les deux effets. La corrélation génétique rG entre effet reine et effet
ouvrières est en général **négative** [BER08 p. 2 ; BB p. 10] : une reine favorable par son effet
reine transmet souvent des ouvrières moins favorables, ce que le profil `standard` ne voit pas.

Valeurs publiées pour une reine (T = a_W + a_Q, valeur totale d'un nœud) :
- `ebv` et `index` : T̂ du groupe d'ouvrières de sa colonie (Â_W + Â_Q du nœud W), valeur
  génétique moyenne d'une fille de l'accouplement [BB éq. 2], comparable à l'index `standard` ;
- `index_worker` : composante ouvrières Â_W du même nœud ;
- `index_queen` : T̂ du nœud Q, valeur reproductrice totale de la reine ; c'est elle qui sert aux
  accouplements (§9).

Écarts-types des index : √(σ²_T ā_WW) pour `index`, √(σ²_AW ā_WW) pour `index_worker`, √σ²_T
pour `index_queen`, avec σ²_T = σ²_AW + σ²_AQ + 2 rG σ_AW σ_AQ et ā_WW = parenté moyenne a(W, W)
des colonies notées (§6).

**Repli** : séparer les deux effets exige de connaître le père (l'accouplement). Si moins de 50 %
des reines notées d'un critère ont un accouplement connu (accouplement connu : accouplement qui a un
nœud sire, §3 ; une fécondation naturelle ou inconnue ne compte pas), ce critère est calculé avec le
profil `standard` et l'avertissement `breeder_fallback` en donne la proportion ; la part est aussi
publiée dans `parameters.known_mating_share`. Sur les données Beekube actuelles, c'est le cas
courant : parmi les 170 critères comptant au moins 10 reines notées, un seul atteint 50 %
d'accouplements connus, avec 14 reines notées (les fécondations sont rarement saisies).

## 6. Paramètres de variance

Par critère, dans cet ordre :
1. h² et répétabilité fournies par Beekube (`source = "input"` ; si une seule est fournie, l'autre
   prend sa valeur par défaut ; répétabilité seule : h² = min(0,25, répétabilité)) ;
2. estimation **REML** (`source = "reml"`) si l'option `estimate_variances` vaut `"always"`, ou `"auto"`
   (défaut) avec au
   moins 100 reines notées, 300 notes et 2 groupes de contemporaines. L'algorithme est un REML sans
   dérivées (DF-REML : vraisemblance restreinte calculée sur la factorisation creuse, maximisée par
   Nelder-Mead [Graser et al. 1987 ; Meyer 1989]). `"never"` désactive le REML. Une variance d'environnement permanent nulle
   (répétabilité = h²) est une estimation valide. Si l'optimiseur échoue ou si une variance devient
   infinie, les valeurs a priori sont gardées (`reml_not_converged`) ; si la variance génétique
   estimée est nulle (héritabilité nulle sur ces données), aussi, avec l'avertissement
   `low_heritability` ;
3. valeurs a priori (`source = "prior"`) : h² = 0,25, répétabilité = 0,30. Sur les données réelles
   de Beekube, la vraisemblance de l'exploitation la plus riche rejette une répétabilité de 0,40
   (intervalle de confiance d'environ 0,1 à moins de 0,4) ; le passage à 0,30 ne change pas les classements.

Conversion : σ²_P = variance des notes après effets fixes ; σ²_a = h² σ²_P (variance génétique de
la colonie) ; σ²_pe = (répétabilité − h²) σ²_P ; σ²_e = (1 − répétabilité) σ²_P ;
σ²_A = σ²_a / (a_WW moyen des colonies notées).

Les variances publiées dans `parameters` (`sigma2_a`, `sigma2_pe`, `sigma2_e` et, pour un critère
calculé en profil `breeder`, `sigma2_aq`) sont dans l'unité des notes ramenées sur [0, 1] (§2) :
pour un critère de type `note`, la variance sur l'échelle brute divisée par (`scale_max` −
`scale_min`)² (le diviseur est celui de chaque note si les échelles diffèrent, `scale_mixed`) ; les
pourcentages et récoltes, déjà sur [0, 1], sont inchangés. Les rapports (h², répétabilité, ρ) ne
dépendent pas de cette unité.

Profil `breeder` : mêmes règles, avec deux paramètres de plus, `queen_ratio` ρ = σ²_AQ / σ²_AW
(défaut 0,5) et `rg` (défaut −0,5, ramené dans [−0,95 ; 0,95]) [BB p. 10]. Si l'un des quatre
paramètres (h², répétabilité, ρ, rG) est fourni, `source = "input"` et les autres prennent leur
valeur par défaut ; quand seule la répétabilité est fournie, h² = min(0,25 ; répétabilité), dans
les deux profils. h² reste l'héritabilité de la **note de colonie** :
σ²_colonie = h² σ²_P = σ²_AW (ā_WW + ρ ā_QQ + 2 rG √ρ ā_WQ), où ā_WW, ā_QQ et ā_WQ sont les
moyennes, sur les colonies notées, des parentés a(W, W), a(Q, Q) et a(W, Q).

REML du profil `breeder` : en `"auto"` (mêmes seuils), seuls σ²_AW et σ²_pe sont estimés, ρ et rG
restent a priori (sur 6 populations simulées, ce REML partiel converge 5 fois sur 6, avec un h²
moyen de 0,18 pour une vraie valeur de 0,22 ; estimations de 0,03 à 0,42 selon le tirage, sur
environ 360 colonies notées). En `"always"`, σ²_AW, σ²_pe, ρ et rG sont tous estimés, mais ils sont
peu identifiables : sur des populations simulées d'environ 360 colonies, le REML complet ne converge
que 2 fois sur 10 ; sur environ 970 colonies, il converge mais ρ estimé varie de 0,10 à 1,29 pour
une vraie valeur de 0,5, et rG de −0,52 à −0,81 pour −0,5. Dans `parameters`, pour un critère
calculé en profil `breeder` (`profile_used = "breeder"`), `sigma2_a` vaut σ²_AW (variance de base de
l'effet ouvrières, et non la variance génétique de la colonie comme en `standard`), `sigma2_aq` =
σ²_AQ et `rg` = rG utilisé ; un critère en repli publie `profile_used = "standard"`, `sigma2_a` est
alors la variance génétique de la colonie et `sigma2_aq` est absent (`null`). Validité du REML :
mêmes règles qu'en `standard` (σ²_pe nulle valide ; σ²_AW et σ²_AQ nulles → a priori et
`low_heritability`) ; en plus, |rG| ≥ 0,94, ou ρ à une borne en `"always"`, rend l'estimation
inutilisable (a priori et `reml_not_converged`). Ces règles sont appliquées par cohérence avec le
profil `standard`, sans validation sur données réelles : aucune exploitation n'a encore 100 reines
notées dont au moins la moitié à accouplement connu (§11).

Critère sans variation : si les notes d'un critère ne varient pas une fois retirés les effets fixes
(par exemple toutes identiques), le critère ne classe pas les reines : pas de REML, index 100 pour
toutes et avertissement `no_variation`.

## 7. Fiabilité

r² = 1 − PEV / (a_ii σ²_A), où PEV est la variance d'erreur de prédiction lue sur l'inverse des
équations du modèle mixte et a_ii σ²_A la variance génétique a priori du nœud. Calcul **exact**
jusqu'à 6 000 équations ; au-delà, **approximation** (`reliability_method = "approx"`,
avertissement `reliability_approx`) : information propre de chaque colonie après absorption des
effets fixes, propagée dans le pedigree [Misztal & Wiggans 1988]. Sur simulation, l'approximation
suit le classement exact (corrélation 0,92 à 0,99) mais **surestime la fiabilité d'environ 0,1**.

Profil `breeder` : r² = 1 − PEV(T) / (a_ii σ²_T) pour la valeur totale T du nœud W publié (σ²_T :
§5 bis) ; PEV(T) = PEV(a_W) + PEV(a_Q) + 2 cov(a_W, a_Q), lues sur le bloc 2 × 2 de l'inverse.
Chaque nœud portant deux effets, il faut environ 5 équations par reine (contre 3 en `standard`) : la
limite de 6 000 équations est atteinte dès 1 200 reines environ (2 000 en `standard`). Au-delà,
l'approximation ci-dessus est appliquée à un modèle à un seul effet équivalent (même variance
génétique d'une note, σ²_A = σ²_colonie / ā_WW), puis la fiabilité obtenue au nœud W est multipliée
par corr²(T_W, g), où g = a_W(W) + a_Q(Q) (effet ouvrières du groupe et effet reine de sa reine) est
la valeur génétique que mesure la note ; sur simulation, corrélation 0,82 à 0,88 avec l'exacte et
surestimation moyenne de 0,04 à 0,05.

## 8. Index

- **Statut** par reine et critère : `evaluated` (notes propres), `pedigree` (pas de note mais une
  apparentée notée : fiabilité ≥ 0,001), `none` (aucune information : pas d'index), `failed`
  (calcul du critère en échec, avertissement `criterion_failed`).
- **Index base 100** : `index = 100 + 10 × sens × (EBV − m) / σ`, m = moyenne des EBV des reines
  de la **base** (reines informées nées dans les `base_years` (5 par défaut) dernières années de
  naissance présentes, base mobile [BER08 p. 2]), σ = écart-type génétique du nœud publié (profil
  `standard` : √σ²_a pour la colonie, √σ²_A pour `index_queen` ; profil `breeder` :
  √(σ²_T ā_WW) pour `index`, √(σ²_AW ā_WW) pour `index_worker`, √σ²_T pour `index_queen`, §5 bis).
  **10 points = 1 écart-type génétique.** Comme en Allemagne, les anciens index « baissent » quand
  la population progresse.
- **Index total** : moyenne des index pondérée par les poids des critères non éliminatoires
  (poids > 0) où la reine a un index ; fiabilité totale = moyenne pondérée des r².
- **Tendance génétique** : moyenne des index et effectif par année de naissance. Le moteur
  ne filtre pas les années : une année avec très peu de reines n'est pas interprétable
  [BER08 p. 2–3].
- **Connexion** : ruchers × années reliés par une reine commune ou des reines apparentées
  (a ≥ 0,125) ; les index ne sont pas comparables entre blocs non reliés. `connectivity.minor_share`
  donne la part des reines notées hors du bloc principal ; l'avertissement `disconnected_groups`
  n'est émis que si elle atteint 10 % (un petit rucher isolé n'avertit pas toute l'exploitation).

## 9. Accouplements

Pour une mère M et une reine productrice de mâles P (vivantes, non éliminées, disponibles, P ≠ M) :
- index attendu des filles par critère = 100 + ½ (index_Q(M) − 100) + ½ (index_Q(P) − 100), à
  partir de l'index de la **reine** (nœud Q) ; une reine sans index compte pour 100 ;
- consanguinité attendue des filles F = ½ a(M, P) (un couple au pedigree inconnu est noté F = 0) ;
- score = index total attendu − `inbreeding_penalty` × F (100 par défaut : 1 % de consanguinité
  coûte 1 point) ; avertissement
  `high_inbreeding` si F ≥ `max_inbreeding` (0,0625 par défaut) ;
- les `n_suggestions` (10 par défaut) meilleures sources de mâles par mère sont proposées ; `POST /v2/mating` évalue des couples
  choisis.

## 10. Validation

- Tests unitaires (`pytest`) : schémas d'entrée, pedigree, parenté (invariants ci-dessus),
  préparation des notes, solveur comparé aux formules de Henderson denses, REML comparé à la
  vraisemblance dense, index, accouplements, API.
- Calibration par simulation (`pytest -m slow`, 20 populations simulées de 4 générations, sires
  de 8 RPM, D = 12). Valeurs mesurées avec ce moteur : corrélation index / vraie valeur 0,571
  (seuil ≥ 0,55 ; corrélation des reines sans note : 0,440), pente de la vraie valeur sur l'index 1,009 (attendu 1 [BB p. 9–10], seuil
  0,8–1,2), fiabilité moyenne annoncée 0,308 pour une corrélation au carré de 0,333 (écart toléré
  0,08). Le REML sous-estime légèrement l'héritabilité simulée : h² moyen 0,144 pour une vraie valeur
  d'environ 0,17 (soit environ 15 % de moins), avec une forte dispersion selon les tirages
  (estimations de 0,05 à 0,33) ; il converge 10 fois sur 10. La corrélation est inférieure à celle
  du prototype d'analyse (environ 0,61) parce que le groupe de contemporaines du moteur (rucher ×
  année) est, dans cette simulation, confondu avec la génération ; un solveur de référence des
  équations du modèle mixte avec le même modèle donne 0,572.
- Profil `breeder` : comparaison aux formules de Henderson denses (EBV et PEV de T), calcul à la
  main sur deux colonies, convergence vers le profil `standard` quand σ²_AQ → 0. Calibration par
  simulation (`pytest -m slow`, 20 populations simulées au schéma de BB p. 10 : σ²_AW = 1,
  σ²_AQ = 0,5, rG = −0,5, 5 générations, 450 reines, accouplements tous connus, paramètres vrais
  fournis) : corrélation index / vraie valeur totale 0,546 (seuil ≥ 0,45), pente de la vraie valeur
  sur l'index 1,082 pour la valeur totale et 1,050 pour la valeur reine (nœud Q) (seuils
  0,8–1,2), 1,120 pour la composante ouvrières (plus bruitée, seuils 0,7–1,4), fiabilité moyenne
  annoncée 0,252 pour une corrélation au carré de 0,307 (écart toléré 0,08). Sur ces données, le
  classement des colonies est très proche de celui du profil `standard` (EBV breeder / standard
  corrélées à 0,98 ; corrélation avec la vraie valeur totale 0,543 en standard contre 0,546 en
  breeder) : l'apport du profil est la séparation reine / ouvrières et l'index reine utilisé aux
  accouplements. Le REML complet (`"always"`) converge sur deux populations simulées d'environ
  1 000 colonies notées, mais ρ et rG y sont peu précis et le REML complet est instable sur de
  plus petites populations (§6).
- Performance : 3 000 reines et 30 000 notes en moins de 60 s pour un critère (environ 24 s
  mesurées, pic mémoire d'environ 640 Mio) ; le REML coûte environ 15 s par critère éligible à
  3 000 reines (5 critères : environ 70 s ; 10 : environ 99 s). Au-delà d'environ 12 critères
  éligibles à cette taille, le délai de 180 s du serveur est dépassé ; fournir h² / répétabilité
  évite le REML.
- Performance du profil `breeder` : 3 000 reines et 30 000 notes en moins de 60 s pour un critère
  (environ 30 s mesurées, dont environ 27 s de REML ; 15 380 équations, fiabilité approchée), avec
  au moins 50 % d'accouplements connus (sinon le critère est calculé en profil `standard`). Le REML
  coûte donc environ 27 s par critère éligible à 3 000 reines : au-delà d'environ 6 critères
  `breeder` éligibles à cette taille, le délai de 180 s du serveur est dépassé (environ 12 en
  profil `standard`) ; fournir h², répétabilité, ρ ou rG évite le REML. En `"always"` (4 paramètres),
  le REML coûte environ 65 s par critère à 3 000 reines, jusqu'à environ 105 s au plafond de 400
  évaluations : à éviter au-delà d'un critère à cette taille.
- Non-régression sur le jeu d'exemple anonymisé (`data/input-example.json`).

## 11. Hypothèses et limites (ce qui n'est pas fait)

- **Évaluation par exploitation** : les index ne sont comparables qu'à l'intérieur d'une
  exploitation, et seulement entre blocs connectés.
- **Notes traitées comme continues** : pas de modèle à seuil pour les notes ordinales.
- **Pas de groupes génétiques** pour les parents inconnus : toutes les reines de base sont
  supposées issues d'une même population.
- **Une seule famille notée** : si toutes les reines notées sont sœurs, la moyenne de la famille
  se confond avec la moyenne générale ; la mère ne reçoit alors aucune information (statut `none`).
- **Père inconnu** : pas de covariance « super-sœurs » entre filles d'une même mère fécondée
  naturellement (hypothèse prudente : moins d'information transmise).
- **h² et répétabilité a priori** (0,25 et 0,30) tant que les effectifs sont faibles ; ce sont
  des valeurs de travail, à affiner par critère [Bienefeld & Pirchner 1990].
- **D et S** : 12 mâles par défaut ; la simulation montre qu'ils changent peu le classement.
- **Fiabilité approchée** au-delà de 6 000 équations : surestimée d'environ 0,1 (profil
  `standard`).
- **Profil `breeder`** (§5 bis) : ρ et rG ne sont estimés qu'en `"always"` et restent imprécis
  jusqu'à un millier de colonies au moins ; par défaut ce sont les valeurs de BB p. 10 (ρ = 0,5,
  rG = −0,5), à confirmer par critère. La fiabilité approchée (> 6 000 équations) surestime
  d'environ 0,05. La composante reine seule (Â_Q du nœud W) n'a pas d'index propre. Repli sur le
  profil `standard`, critère par critère, sous 50 % d'accouplements connus. Le profil n'a pas
  été validé sur données réelles : aucune exploitation n'a assez de reines notées à accouplement
  connu (il en faut au moins la moitié des reines notées d'un critère, et 100 reines notées pour
  le REML) ; le seul critère à 50 % ne compte que 14 reines notées. Ses réglages viennent de la
  simulation et du profil `standard`.
- **Évaluation univariée** : chaque critère est ajusté séparément ; l'index total est une moyenne
  pondérée des index par critère, pas un index de sélection économique. Aucune corrélation
  génétique entre critères ; la covariance reine / ouvrières n'est modélisée que par le profil
  `breeder` (§5 bis).
- **Groupe de contemporaines fixe** (rucher × année) : il peut absorber des différences
  génétiques quand un rucher-année ne contient qu'une famille ou qu'une génération ; c'est ce qui
  abaisse la corrélation mesurée en simulation (§10). Les avertissements de connexion aident à
  repérer ce cas.
- La classe de ruche est la seule information de ruche utilisée (`hive_type` est reçu mais non
  modélisé).

## Références

- [BB] Brascamp E.W., Bijma P. (2014). Methods to estimate breeding values in honey bees.
  *Genetics Selection Evolution* 46:53.
- [BER08] Bienefeld K., Ehrhardt K., Reinhardt F. (2008). Noticeable success in honey bee
  selection after the introduction of genetic evaluation using BLUP. *American Bee Journal*,
  août 2008, 739–742.
- Bienefeld K., Ehrhardt K., Reinhardt F. (2007). Genetic evaluation in the honey bee considering
  queen and worker effects — a BLUP-animal model approach. *Apidologie* 38:77–85.
- Bienefeld K., Pirchner F. (1990). Heritabilities for several colony traits in the honeybee
  (*Apis mellifera carnica*). *Apidologie* 21:175–183.
- Bienefeld K., Reinhardt F., Pirchner F. (1989). Inbreeding effects of queen and workers on colony
  traits in the honeybee. *Apidologie* 20:439–450.
- Graser H.-U., Smith S.P., Tier B. (1987). A derivative-free approach for estimating variance
  components in animal models by restricted maximum likelihood. *Journal of Animal Science*
  64:1362–1370.
- Henderson C.R. (1975). Best linear unbiased estimation and prediction under a selection model.
  *Biometrics* 31:423–447.
- Meyer K. (1989). Restricted maximum likelihood to estimate variance components for animal models
  with several random effects using a derivative-free algorithm. *Genetics Selection Evolution*
  21:317–340.
- Misztal I., Wiggans G.R. (1988). Approximation of prediction error variance in large-scale
  animal models. *Journal of Dairy Science* 71 (suppl. 2):27–32.
- Willham R.L. (1963). The covariance between relatives for characters composed of components
  contributed by related individuals. *Biometrics* 19:18–27.

Analyse détaillée de la littérature et des écarts avec le moteur v1 :
[`docs/analyse-litterature-blup.fr.md`](docs/analyse-litterature-blup.fr.md).
