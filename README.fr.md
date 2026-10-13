# Beekube BLUP Mellifera API

## Description
Beekube BLUP Mellifera API est une application Flask qui fournit une interface API pour effectuer des calculs BLUP (Best Linear Unbiased Prediction) sur des données d'élevage d'abeilles. Cette API a été conçue pour être utilisée avec le système de gestion apicole Beekube, mais le modèle de données peut être utilisée pour d'autres fins.


*Lisez ceci dans d'autres langues : [English](README.md), [Français](README.fr.md)*

## Fonctionnalités
- Calcul BLUP pour l'évaluation génétique des reines d'abeilles
- Gestion des données de pedigree et d'évaluation
- Interface API RESTful avec documentation OpenAPI

## Prérequis
- Docker
- Docker Compose ≥ 2.24 (la syntaxe `env_file` avec `required` n'est pas prise en charge par les versions antérieures)

## Installation et démarrage
1. Clonez ce dépôt :
   ```
   git clone beekube-app/beekube-blup-mellifera
   cd beekube-blup-mellifera
   ```

2. Lancez l'application avec Docker Compose :
   ```
   docker compose up --build
   ```

3. L'API sera accessible à l'adresse : `http://localhost:8081`

## Utilisation
L'API expose les endpoints suivants :

- `GET /` : Page d'accueil de l'API
- `POST /blup` : Effectue le calcul BLUP (moteur v1, R)
- `POST /v2/blup` : Évaluation génétique v2 (parenté abeille, index base 100, accouplements) — méthode décrite dans [METHODOLOGIE.md](METHODOLOGIE.md)
- `POST /v2/mating` : Index et consanguinité attendus pour des couples mère × source de mâles choisis

Pour plus de détails sur l'utilisation de l'API, consultez la documentation OpenAPI disponible à l'adresse de l'API une fois lancée :
`http://localhost:8081/openapi`

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

- `crawl --full` / `crawl --incremental` : parcourt le site ; un passage interrompu est repris par le passage suivant de même mode ; tant qu'aucun crawl complet n'a abouti, `--incremental` reprend le complet ; si la liste des pays (`/all/`) est illisible, le passage sort en erreur (code 1) ;
- `reparse [--kind countries|country|queens_list|season|queen]` : reconstruit la base depuis le cache brut, sans réseau ;
- `status` : dernier passage, volumes, pages par type et statut, critères non mappés ;
- `export-blup --breeder CBP --breeder 200 [--from 2018] [--sources season_table,radar] [--profile breeder] -o data/mellifera/export.json` : exporte les reines de ces éleveurs et leur ascendance. `--breeder` accepte un trigramme (casse indifférente) ou un id Index Mellifera ; un trigramme porté par plusieurs éleveurs est refusé. Seuls les tableaux d'évaluations sont exportés par défaut ; `radar` ajoute les valeurs annuelles du graphique des fiches ;
- `schedule` : crawl incrémental hebdomadaire (`MELLIFERA_SCHEDULE`).

Codes de sortie : `0` succès, `1` arrêt après des erreurs réseau consécutives, `2` plus de 5 % d'échecs de parsing sur un type de page (structure du site modifiée ?), `3` un autre crawl tient le verrou, `4` erreur de configuration (contact vide, variable invalide, base absente, éleveur inconnu ou ambigu).

## Structure du projet
- `app.py` : Point d'entrée de l'application Flask
- `blup.r` : Script R pour les calculs BLUP
- `Dockerfile-Python` : Dockerfile pour construire l'image Python
- `docker-compose.yml` : Configuration Docker Compose
- `Dockerfile-Mellifera` : Dockerfile de l'image du crawler Index Mellifera (sans R ni numpy)
- `mellifera/` : crawler Index Mellifera (base SQLite, export `BlupV2Input`)
- `data/` : Répertoire pour les fichiers de données temporaires

## Développement
Pour contribuer au projet :

1. Forkez le dépôt
2. Créez votre branche de fonctionnalité (`git checkout -b feature/AmazingFeature`)
3. Committez vos changements (`git commit -m 'Add some AmazingFeature'`)
4. Poussez vers la branche (`git push origin feature/AmazingFeature`)
5. Ouvrez une Pull Request

## Tests

```bash
python3 --version                        # doit être ≥ 3.11 ; sinon appeler explicitement, par ex., python3.12
python3 -m venv .venv && .venv/bin/pip install -r requirements-dev.txt
.venv/bin/pytest -v                      # API Flask (Rscript simulé)
docker compose build api                 # obligatoire après toute modification de blup.r ou de tests/
docker compose run --rm api Rscript tests/r/test_blup.R   # tests R
.venv/bin/pytest -m docker -v            # bout en bout sur data/input-example.json
.venv/bin/pytest -m slow -v              # moteur v2 : calibration par simulation et performance
.venv/bin/python scripts/export_openapi.py                # après un changement de modèle pydantic
```

`pytest -m docker` construit et démarre le service `api` dans son propre projet Compose, `blup-test` par défaut (modifiable avec `COMPOSE_PROJECT_NAME`, sans le pointer vers votre projet habituel), puis le supprime, sans interférer avec votre projet Compose habituel. Le test échoue immédiatement, avant toute commande Compose, si un service écoute déjà sur le port de test (8081 par défaut). Sur un poste où l'API tourne déjà sur le port 8081, utilisez un port libre : `BLUP_API_PORT=8082 BLUP_API_URL=http://localhost:8082 .venv/bin/pytest -m docker -v`.

Les tests du crawler (`tests/mellifera/`) ne contactent jamais Index Mellifera : ils utilisent un faux site (`httpx.MockTransport`) et des pages réelles anonymisées (`tests/mellifera/fixtures/`). Les pages brutes d'où viennent ces fixtures (`tests/mellifera/fixtures/raw/`, ignorées par git) contiennent des noms et restent sur la machine qui les a téléchargées ; `scripts/make_mellifera_fixtures.py` régénère les fixtures à partir d'elles.

## Licence
MIT Licence

## Contact
Beekube - Damien TUPINIER
https://www.beekube.com