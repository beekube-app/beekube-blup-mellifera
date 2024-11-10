from flask import jsonify
from flask_openapi3 import OpenAPI, Info, Tag, Server
from pydantic import BaseModel, Field, RootModel
import subprocess
import json
import os
from typing import List, Optional, Union, Dict, Any

# Récupérer l'URL de base depuis les variables d'environnement
API_BASE_URL = os.getenv('API_BASE_URL', 'http://localhost:8081')

info = Info(title="Beekube BLUP Mellifera API", version="1.0.0")
servers = [
    Server(url=API_BASE_URL),
]
app = OpenAPI(__name__, info=info, servers=servers)

blup_tag = Tag(name="blup", description="BLUP operations")

# Input Models
class Evaluation(BaseModel):
    note: float
    user: int
    apiary: int
    beehiveType: str
    month: str
    year: str

class QueenBeeEvaluations(RootModel):
    root: Dict[str, List[Evaluation]]

class QueenBeeInput(BaseModel):
    queenbee: int
    queenbee_parent: Optional[int] = None
    drone_parent: Optional[int] = None
    apiary_default: Optional[int] = None
    born: str
    evaluate: QueenBeeEvaluations = Field(default_factory=dict)

class BLUPInput(BaseModel):
    exploitation: int
    evaluate: List[int]
    evaluate_elimination: List[int]
    data: List[QueenBeeInput]

# Output Models
class HeritabilityStats(BaseModel):
    heritability: float
    se: float
    v_additive_queen: float
    v_additive_drone: float
    v_colony: float
    v_residual: float

class QueenBeeOutput(BaseModel):
    queenbee: str
    queenbee_parent: str
    drone_parent: str
    apiary_default: Optional[int]
    born: str
    blups: Dict[str, Optional[float]]
    methods: Dict[str, Optional[str]]

class BLUPResultOutput(BaseModel):
    blup: List[QueenBeeOutput]
    heritabilities: Dict[str, HeritabilityStats]

class BLUPOutput(BaseModel):
    status: str
    results: Union[BLUPResultOutput, dict]

class WelcomeResponse(BaseModel):
    status: str
    results: str

class ErrorResponse(BaseModel):
    status: str
    message: str
    output: str
    error: str
    results: dict

@app.get("/",
         tags=[blup_tag],
         responses={
             200: WelcomeResponse,
         })
def welcome():
    """Welcome endpoint"""
    return jsonify({
        'status': 'success',
        'results': "Beekube BLUP Mellifera"
    })

@app.post("/blup",
          tags=[blup_tag],
          responses={
              200: BLUPOutput,
              500: ErrorResponse
          })
def calulate(body: BLUPInput):
    """
    Perform BLUP calculation
    """
    try:
        # Make sure that the directory "data" exists
        os.makedirs('./data', exist_ok=True)

        # Clean up any existing files
        for filename in ['input.json', 'resultats_index.json']:
            filepath = f'./data/{filename}'
            if os.path.exists(filepath):
                os.remove(filepath)

        # Save input data
        with open('./data/input.json', 'w') as f:
            json.dump(body.dict(), f)

        # Run R script
        result = subprocess.run(['Rscript', 'blup.r'], capture_output=True, text=True)

        # Read results
        if os.path.exists('./data/resultats_index.json'):
            with open('./data/resultats_index.json', 'r') as f:
                results = json.load(f)

            # Parse results into our output model
            response = BLUPOutput(
                status='success',
                results=results
            )
        else:
            response = BLUPOutput(
                status='error',
                results={"error": "Results file not found"}
            )

        # Clean up
        for filename in ['input.json', 'resultats_index.json']:
            filepath = f'./data/{filename}'
            if os.path.exists(filepath):
                os.remove(filepath)

        return jsonify(response.model_dump())

    except Exception as e:
        error_response = ErrorResponse(
            status='error',
            message=str(e),
            output=result.stdout if 'result' in locals() else '',
            error=result.stderr if 'result' in locals() else '',
            results={"error": "An unexpected error occurred"}
        )
        return jsonify(error_response.model_dump()), 500

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8081)