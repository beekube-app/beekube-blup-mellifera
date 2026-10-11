# Feuille de route — fiabilité du BLUP

Objectif : que l'index Beekube soit **interprétable, vérifiable et fiable**, au point que les
sélectionneurs s'en servent pour décider. Cette liste est issue des questions d'un utilisateur
sélectionneur (Discord #evaluation-selection, 11/10/2026) et d'une relecture de `blup.r`.

Repères : **[Qn]** = question correspondante du message (1 groupe de mâles, 2 valeur publiée,
3 variances, 4 héritabilités, 5 inter-exploitations, 6 fiabilité, 7 notes répétées,
8 accouplements). Côté **[BLUP]** = ce dépôt, côté **[Beekube]** = application.

## Déjà corrigé

- Effet mâle additionné par position au lieu d'être associé à chaque reine : les index publiés
  étaient faux dès qu'un `drone_parent` était renseigné (corrélation de rang 0,54 avec la bonne
  valeur sur données simulées, 1,00 après correction). [Q2]
- `(1|drone_parent)` présent deux fois dans la formule du modèle BLUP. [Q2]
- Diversité : les parents inconnus (`"0"`) comptaient comme un ancêtre commun (deux reines
  d'origine inconnue avaient une diversité de 0), et une mère avec sa fille obtenait une
  diversité plus élevée que deux reines sans lien. [Q8]
- Docker : `pydantic` 2.9.2 → 2.12.3 (l'image `r-base` embarque Python 3.14).

## P0 — Bugs et honnêteté des sorties (rapide, à faire d'abord)

1. **Échecs silencieux** [BLUP] : dans les `error = function(e)` de `calculate_blup`, les
   affectations `blups[[critere]] <- …` sont locales au gestionnaire (portée R) : un critère en
   échec disparaît du résultat sans trace. Utiliser `<<-` ou retourner une valeur.
   La méthode de repli `lm` n'est jamais exécutée (`lmer` n'est pas dans un `try`).
2. **SE de l'héritabilité fictive** [BLUP] [Q4] : `sqrt(2/n)` n'est pas une erreur standard.
   La remplacer par une vraie estimation (profil de vraisemblance ou bootstrap paramétrique
   `bootMer`) ou ne rien renvoyer.
3. **Noms trompeurs** [BLUP] [Q3] : `v_colony` est la variance rucher ; `v_additive_queen` est
   une variance entre reines (génétique + environnement permanent) ; l'« héritabilité » est une
   répétabilité. Renommer (avec compatibilité du SDK) ou documenter.
4. **Inconnus regroupés dans le modèle** [BLUP] [Q1] : toutes les reines sans `drone_parent`
   partagent le niveau `"0"` de l'effet mâle, comme si elles avaient le même père. Les exclure
   de cet effet (ou créer des groupes génétiques explicites : origine, année, fournisseur).
5. **Critères éliminatoires** [BLUP] : `evaluate_elimination` est reçu mais jamais utilisé.
   Exclure les reines éliminées des suggestions d'accouplement, ou retirer le champ.
6. **Dépendances inutiles** [BLUP] : `brms`, `rstan`, `lm.beta` chargés sans être utilisés
   (image lourde, démarrage lent).

## P1 — Fiabilité par reine (ce qui rend l'index utilisable)

7. **Fiabilité / précision par reine** [BLUP] [Q6] : calculer la PEV (`ranef(condVar = TRUE)`)
   et renvoyer pour chaque reine et critère `r² = 1 − PEV / σ²` ainsi que le nombre de notes,
   d'évaluateurs et d'années utilisés.
8. **Affichage** [Beekube] [Q6] : montrer la fiabilité et le nombre de notes à côté de chaque
   BLUP ; griser ou masquer les index sous un seuil (par ex. r² < 0,3) ; ne pas classer des
   reines dont l'index repose sur une seule note.
9. **Exposer les paramètres génétiques** [Beekube] [Q3][Q4] : les héritabilités et les
   variances sont renvoyées par l'API mais ni stockées ni affichées. Les stocker par calcul et
   les montrer (avec effectifs et SE) dans une page « Méthode / paramètres ».
10. **Documentation publique de la méthode** [BLUP] : `METHODOLOGIE.md` décrivant le modèle
    réel, ses hypothèses et ses limites (ce qui est fait et ce qui ne l'est pas), avec des
    références. Lien depuis Beekube.
11. **Tests de validation** [BLUP] : intégrer les simulations utilisées pour le diagnostic (vraie
    valeur connue → corrélation index/valeur vraie, calibration de la fiabilité, scénarios avec
    parents inconnus) et les lancer en CI.

## P2 — Modèle génétique (passer d'un modèle de répétabilité à un vrai BLUP)

12. **Utiliser la parenté dans le modèle** [BLUP] [Q2][Q5] : aujourd'hui la matrice A ne sert
    qu'en poids (sa diagonale) ; les apparentées ne contribuent pas à l'index. Passer à un modèle
    animal avec matrice de parenté (`pedigreemm`, `sommer`, ou BLUPF90 pour la montée en charge).
13. **Séparer génétique et environnement permanent** [BLUP] [Q3][Q7] : ajouter un effet
    colonie/environnement permanent distinct de l'effet génétique, sinon les répétitions gonflent
    la « variance génétique ».
14. **Parenté spécifique à l'abeille** [BLUP] [Q1][Q2] : matrice de parenté haplo-diploïde avec
    groupes de reines productrices de mâles (Bienefeld et al. 2007 ; Brascamp & Bijma 2014),
    plutôt que de traiter la reine productrice de mâles comme un « père » diploïde.
15. **Effets reine / ouvrières** [BLUP] [Q2][Q3] : modèle à deux effets (maternel = reine,
    direct = ouvrières) avec covariance estimée, et publication de la valeur totale
    (éq. 2 de Brascamp & Bijma) ainsi que des deux composantes. Nécessite des données suffisantes ;
    sinon utiliser des paramètres de la littérature.
16. **Groupes de mâles dans les données** [Beekube] [Q1] : modéliser un « groupe de reines
    productrices de mâles » (reines sœurs + leur mère, rucher de fécondation, année), au lieu d'une
    seule reine `drone`. Transmettre aussi les mâles en référence externe (`droneRefExt`) comme
    groupe identifié plutôt que comme inconnu.
17. **Paramètres de variance robustes** [BLUP] [Q4] : quand les effectifs sont faibles, utiliser
    des héritabilités publiées par critère (par ex. Bernstein et al. 2023) comme valeurs fixées
    ou a priori, et n'estimer que lorsqu'un seuil de colonies est atteint (seuil à établir par
    simulation et à documenter).
18. **Effets d'environnement** [BLUP] [Q7] : groupes de contemporaines rucher × année (× saison)
    plutôt que mois/année/rucher séparés ; évaluateur en effet fixe quand il y en a peu.
19. **Notes d'échelle** [BLUP] [Q7] : modèle à seuil (`ordinal::clmm` ou `brms` cumulatif) pour
    les critères notés sur échelle, ou au minimum vérification que le modèle linéaire reste
    acceptable sur les distributions réelles.

## P3 — Comparabilité et accouplements

20. **Évaluation inter-exploitations** [Beekube][BLUP] [Q5] : évaluation commune facultative
    (consentement explicite) entre exploitations qui échangent des reines, avec contrôle de la
    connexion génétique ; les sœurs évaluées ailleurs contribuent alors à l'index.
21. **Diversité fondée sur la parenté** [BLUP] [Q8] : remplacer le ratio d'ancêtres communs par
    la coparenté issue de la matrice A (consanguinité attendue de la descendance), voie paternelle
    comprise. Signaler les couples dont le pedigree est inconnu au lieu de les noter « diversité
    maximale ».
22. **Score total cohérent** [BLUP][Beekube] [Q8] : utiliser l'index pondéré de l'utilisateur
    (comme le « Score BLUP » affiché) au lieu de la moyenne non pondérée des critères ; ne plus
    supposer une plage de BLUP de −0,5 à 0,5 (standardiser) ; rendre le poids 0,7 / 0,3
    configurable et l'afficher.
23. **Sens des critères** [BLUP] : vérifier que tous les critères sont orientés « plus = mieux »
    avant de les agréger (agressivité, essaimage…).
