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

## Index Mellifera crawler

The `mellifera/` package reads the public pages of [Index Mellifera](https://index-mellifera.org) (countries, breeders, queen lists, queen pages, yearly evaluation tables) with a polite crawl, stores them in a SQLite database (`data/mellifera/mellifera.sqlite`) and exports the queens of chosen breeders, with their whole ancestry, as a `BlupV2Input` file for the v2 engine. Design (in French): [docs/superpowers/specs/2026-10-13-index-mellifera-design.md](docs/superpowers/specs/2026-10-13-index-mellifera-design.md).

The crawler only sends GET requests to an allow-list of public paths, waits at least 1 s between two requests (2 s by default), slows down on errors and identifies itself with a contact address. The database keeps only the Index Mellifera id, trigram and country of each breeder. The raw page cache (`data/mellifera/cache/`) contains names: never commit or publish it.

| Variable | Default | Role |
|---|---|---|
| `MELLIFERA_CONTACT` | (empty) | e-mail or URL in the `User-Agent`; **required**: the crawl refuses to start without it |
| `MELLIFERA_DATA` | `data/mellifera` (`/app/data/mellifera` in the image) | database, raw cache and lock file |
| `MELLIFERA_DELAY` | `2.0` | seconds between two requests (at least 1) |
| `MELLIFERA_MAX_ERRORS` | `10` | consecutive network errors before a clean stop (with the delay doubling up to 5 min, an outage of about twenty minutes is tolerated). A page that exceeds the read timeout (30 s) three times is skipped for this pass, logged and recorded as an error; the crawl goes on |
| `MELLIFERA_SCHEDULE` | `sun 02:00` | day and time of the weekly incremental crawl |
| `MELLIFERA_TZ` | `Europe/Paris` | time zone of the schedule |

With Docker, the `mellifera` service belongs to the `mellifera` profile: `docker compose up` does not start it. Put `MELLIFERA_CONTACT` in `.env.production`, then:

```bash
docker compose run --rm mellifera python -m mellifera crawl --full   # initial crawl: more than a day, resumed if interrupted
docker compose up -d mellifera                                       # weekly incremental crawl
docker compose run --rm mellifera python -m mellifera status
```

Commands (`python -m mellifera …`, or `.venv/bin/python -m mellifera …` outside Docker):

- `crawl --full` / `crawl --incremental`: crawl the site; an interrupted run is resumed by the next run of the same mode; until a full crawl has completed, `--incremental` resumes the full crawl; if the country list (`/all/`) cannot be read, the run exits with an error (code 1);
- `reparse [--kind countries|country|queens_list|season|queen]`: rebuild the database from the raw cache, without network;
- `status`: last run, row counts, pages by type and status, unmapped criteria;
- `export-blup --breeder CBP --breeder 200 [--from 2018] [--sources season_table,radar] [--profile breeder] -o data/mellifera/export.json`: export the queens of these breeders and their ancestry. `--breeder` takes a trigram (any case) or an Index Mellifera id; a trigram shared by several breeders is refused. Only the evaluation tables are exported by default; `radar` adds the yearly values of the queen page chart;
- `schedule`: weekly incremental crawl (`MELLIFERA_SCHEDULE`).

Exit codes: `0` success, `1` stopped after consecutive network errors, `2` more than 5 % of the pages of one type failed to parse (site structure changed?), `3` another crawl holds the lock, `4` configuration error (empty contact, invalid variable, missing database, unknown or ambiguous breeder).

### Index Mellifera lineages in the engine

When a queen cites an Index Mellifera lineage by its code (`mother_ref_ext` for her mother, `mating.external_ref` for the mother of the drone-producing queens), the v2 engine reads the crawler database read-only and adds that lineage's ancestry to the pedigree (at most `options.external_pedigree_depth` generations, 7 by default, 0 = disabled, and at most 1,000 ancestors in total). These ancestors only feed relationship and inbreeding; they never appear in any output. Since version 2.2.0, a criterion that declares its standard trait (`criteria[].trait`, a list frozen in the contract) also receives the Index Mellifera records of these ancestors and of the queens recorded in the same tester × year groups (at most 1,000), with no fixed effect shared with the farm's records (see `METHODOLOGIE.md`, §3 bis). A code is looked up as given, then without spaces and case-insensitively (case matters on the site: `DE-ABg-B-0003-2017`; two codes differing only by case make the input ambiguous), then by displayed name (`B4536(NG)`) when unique; an unknown code is ignored (`external_ref_unknown` warning).

Activation: `MELLIFERA_DB=/app/data/mellifera/mellifera.sqlite` in the `api` service's `.env.production`. Without this variable, resolution is disabled and no new warning is emitted; if the database is missing or unreadable, an input citing lineages gets the `external_pedigree_unavailable` warning. Version 2.1.0 adds the `external_pedigree_depth` option and three warning codes (`external_ref_resolved`, `external_ref_unknown`, `external_pedigree_unavailable`) to the contract, version 2.2.0 the `criteria[].trait` field and two codes (`external_records_added`, `external_records_unavailable`): **enable `MELLIFERA_DB` in production only after regenerating and publishing the PHP SDK** from `openapi/openapi.json`.

## Project Structure
- `app.py` : Entry point of the Flask application
- `blup.r` : R script for BLUP calculations
- `Dockerfile-Python` Dockerfile to build the Python image
- `docker-compose.yml` : Configuration Docker Compose
- `Dockerfile-Mellifera` : Dockerfile of the Index Mellifera crawler image (no R, no numpy)
- `mellifera/` : Index Mellifera crawler (SQLite database, `BlupV2Input` export)
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

The crawler tests (`tests/mellifera/`) never contact Index Mellifera: they use a fake site (`httpx.MockTransport`) and anonymized real pages (`tests/mellifera/fixtures/`). The raw pages these fixtures come from (`tests/mellifera/fixtures/raw/`, ignored by git) contain names and stay on the machine that downloaded them; `scripts/make_mellifera_fixtures.py` regenerates the fixtures from them.

## Licence
MIT Licence

## Contact
Beekube - Damien TUPINIER
https://www.beekube.com