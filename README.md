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
- Docker Compose

## Installation and start-up
1. Clone this repository:
   ```
   git clone beekube-app/beekube-blup-mellifera
   cd beekube-blup-mellifera
   ```

2. Start the application with Docker Compose:
   ```
   docker-compose up --build
   ```

3. The API will be accessible at the address: `http://localhost:8081`

## Utilization
The API exposes the following endpoints:

- `GET /` Homepage of the API
- `POST /blup` Perform the calculation BLUP

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

## Licence
MIT Licence

## Contact
Beekube - Damien TUPINIER
https://www.beekube.com