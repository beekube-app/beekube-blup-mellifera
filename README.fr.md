# Beekube BLUP Melifera API

## Description
Beekube BLUP Melifera API est une application Flask qui fournit une interface API pour effectuer des calculs BLUP (Best Linear Unbiased Prediction) sur des données d'élevage d'abeilles. Cette API a été conçue pour être utilisée avec le système de gestion apicole Beekube, mais le modèle de données peut-elle être utilisée pour d'autres fins.


*Lisez ceci dans d'autres langues : [English](README.md), [Français](README.fr.md)*

## Fonctionnalités
- Calcul BLUP pour l'évaluation génétique des reines d'abeilles
- Gestion des données de pedigree et d'évaluation
- Interface API RESTful avec documentation OpenAPI

## Prérequis
- Docker
- Docker Compose

## Installation et démarrage
1. Clonez ce dépôt :
   ```
   git clone beekube-app/beekube-blup-melifera
   cd beekube-blup-melifera
   ```

2. Lancez l'application avec Docker Compose :
   ```
   docker-compose up --build
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

## Licence
MIT Licence

## Contact
Beekube - Damien TUPINIER
https://www.beekube.com