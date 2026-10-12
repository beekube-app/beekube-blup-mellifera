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
- `POST /blup` : Effectue le calcul BLUP

Pour plus de détails sur l'utilisation de l'API, consultez la documentation OpenAPI disponible à l'adresse de l'API une fois lancée :
`http://localhost:8081/openapi`

## Structure du projet
- `app.py` : Point d'entrée de l'application Flask
- `blup.r` : Script R pour les calculs BLUP
- `Dockerfile-Python` : Dockerfile pour construire l'image Python
- `docker-compose.yml` : Configuration Docker Compose
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
python3 -m venv .venv && .venv/bin/pip install -r requirements-dev.txt
.venv/bin/pytest -v                      # API Flask (Rscript simulé)
docker compose build api                 # obligatoire après toute modification de blup.r ou de tests/
docker compose run --rm api Rscript tests/r/test_blup.R   # tests R
.venv/bin/pytest -m docker -v            # bout en bout sur data/input-example.json
.venv/bin/python scripts/export_openapi.py                # après un changement de modèle pydantic
```

`pytest -m docker` construit et démarre le service `api` dans son propre projet Compose, `blup-test` par défaut (modifiable avec `COMPOSE_PROJECT_NAME`, sans le pointer vers votre projet habituel), puis le supprime, sans interférer avec votre projet Compose habituel. Le test échoue immédiatement, avant toute commande Compose, si un service écoute déjà sur le port de test (8081 par défaut). Sur un poste où l'API tourne déjà sur le port 8081, utilisez un port libre : `BLUP_API_PORT=8082 BLUP_API_URL=http://localhost:8082 .venv/bin/pytest -m docker -v`.

## Licence
MIT Licence

## Contact
Beekube - Damien TUPINIER
https://www.beekube.com