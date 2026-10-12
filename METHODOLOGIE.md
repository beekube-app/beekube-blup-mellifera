# Méthodologie du moteur BLUP v2 (Beekube)

Ce document décrit ce que calcule réellement le moteur v2 (`POST /v2/blup`, paquet `blup2/`),
les hypothèses retenues et les limites connues. Il s'adresse aux apiculteurs et sélectionneurs qui
lisent les index dans Beekube, et aux personnes qui relisent le code. Les références entre crochets
renvoient à la section « Références » ; [BB] = Brascamp & Bijma 2014, [BER08] = Bienefeld,
Ehrhardt & Reinhardt 2008. Version du moteur : `2.0.0`, profil `standard`.

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
  répétabilité)).

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

## 8. Index

- **Statut** par reine et critère : `evaluated` (notes propres), `pedigree` (pas de note mais une
  apparentée notée : fiabilité ≥ 0,001), `none` (aucune information : pas d'index), `failed`
  (calcul du critère en échec, avertissement `criterion_failed`).
- **Index base 100** : `index = 100 + 10 × sens × (EBV − m) / σ`, m = moyenne des EBV des reines
  de la **base** (reines informées nées dans les `base_years` (5 par défaut) dernières années de
  naissance présentes, base
  mobile [BER08 p. 2]), σ = écart-type génétique du nœud publié (√σ²_a pour la colonie,
  √σ²_A pour `index_queen`). **10 points = 1 écart-type génétique.** Comme en Allemagne, les
  anciens index « baissent » quand la population progresse.
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
- Performance : 3 000 reines et 30 000 notes en moins de 60 s pour un critère (environ 24 s
  mesurées, pic mémoire d'environ 640 Mio) ; le REML coûte environ 15 s par critère éligible à
  3 000 reines (5 critères : environ 70 s ; 10 : environ 99 s). Au-delà d'environ 12 critères
  éligibles à cette taille, le délai de 180 s du serveur est dépassé ; fournir h² / répétabilité
  évite le REML.
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
- **Fiabilité approchée** au-delà de 6 000 équations : surestimée d'environ 0,1.
- **Profil `breeder`** (effets reine et ouvrières séparés, corrélation génétique négative
  [BER08 p. 2 ; BB p. 10]) : décrit dans un document séparé ; tant qu'il n'est pas disponible, le
  moteur applique le profil `standard` et renvoie `breeder_fallback`.
- **Évaluation univariée** : chaque critère est ajusté séparément ; l'index total est une moyenne
  pondérée des index par critère, pas un index de sélection économique. Aucune corrélation
  génétique entre critères, ni covariance reine / ouvrières (profil `breeder`, plan ultérieur).
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
