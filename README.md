# Beekube BLUP Mellifera API

## Description
Beekube BLUP Mellifera API is a Flask application that provides an API interface to perform BLUP (Best Linear Unbiased Prediction) calculations on bee breeding data. This API was designed to be used with the Beekube beekeeping management system, but the data model can be used for other purposes.


*Read this in other languages: [English](README.md), [Français](README.fr.md)*


## Features
- BLUP calculation for the genetic evaluation of queen bees
- Management of pedigree and evaluation data
- RESTful API interface with OpenAPI documentation
- 
## Prerequisites
- Docker
- Docker Compose ≥ 2.24 (the `env_file` syntax with `required` is not supported by older versions)

## Installation and start-up
1. Clone this repository:
   ```
   git clone beekube-app/beekube-blup-mellifera
   cd beekube-blup-mellifera
   ```

2. Start the application with Docker Compose:
   ```
   docker compose up --build
   ```

3. The API will be accessible at the address: `http://localhost:8081`

## Utilization
The API exposes the following endpoints:

- `GET /` Homepage of the API
- `POST /blup` Perform the calculation BLUP (v1 engine, R)
- `POST /v2/blup` Genetic evaluation v2 (honey bee relationships, base-100 indexes, matings) — method described (in French) in [METHODOLOGIE.md](METHODOLOGIE.md)
- `POST /v2/mating` Expected indexes and inbreeding for chosen mother × drone source pairs

For more details on the use of the API, refer to the OpenAPI documentation available at the API address once it is launched:
`http://localhost:8081/openapi`

## Project Structure
- `app.py` : Entry point of the Flask application
- `blup.r` : R script for BLUP calculations
- `Dockerfile-Python` Dockerfile to build the Python image
- `docker-compose.yml` : Configuration Docker Compose
- `data/` : Directory for temporary data files

## Development
To contribute to the project:

1. Clone the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## Tests

```bash
python3 --version                        # must be >= 3.11; otherwise call e.g. python3.12 explicitly
python3 -m venv .venv && .venv/bin/pip install -r requirements-dev.txt
.venv/bin/pytest -v                      # Flask API (Rscript mocked)
docker compose build api                 # required after any change to blup.r or tests/
docker compose run --rm api Rscript tests/r/test_blup.R   # R tests
.venv/bin/pytest -m docker -v            # end-to-end on data/input-example.json
.venv/bin/pytest -m slow -v              # v2 engine: simulation calibration and performance
.venv/bin/python scripts/export_openapi.py                # after any pydantic model change
```

`pytest -m docker` builds and starts the `api` service in its own Compose project, `blup-test` by default (override with `COMPOSE_PROJECT_NAME`, but do not point it at your usual project), then removes it, so it does not interfere with your usual Compose project. It fails right away, before running any Compose command, if a service already listens on the test port (8081 by default). On a machine where the API already runs on 8081, point the test to a free port: `BLUP_API_PORT=8082 BLUP_API_URL=http://localhost:8082 .venv/bin/pytest -m docker -v`.

## Licence
MIT Licence

## Contact
Beekube - Damien TUPINIER
https://www.beekube.com