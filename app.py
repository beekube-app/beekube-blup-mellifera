from flask import jsonify
from flask_openapi3 import OpenAPI, Info, Tag
from pydantic import BaseModel, Field, RootModel
import subprocess
import json
import os
from typing import List, Optional, Union, Dict

info = Info(title="Beekube BLUP Melifera API", version="1.0.0")
app = OpenAPI(__name__, info=info)

blup_tag = Tag(name="blup", description="BLUP operations")


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
    apiary: Optional[int] = None
    born: str
    evaluate: QueenBeeEvaluations = Field(default_factory=dict)

class QueenBeeOutput(BaseModel):
    queenbee: int
    queenbee_parent: Optional[int]
    drone_parent: Optional[int]
    apiary: Optional[int]
    born: str
    blups: dict
    methods: dict


class BLUPInput(BaseModel):
    exploitation: int
    evaluate: List[int]
    evaluate_elimination: List[int]
    data: List[QueenBeeInput]


class BLUPOutputDebug(BaseModel):
    status: str
    output: str
    error: str
    results: Union[List[QueenBeeOutput], dict]


class BLUPOutput(BaseModel):
    status: str
    results: Union[List[QueenBeeOutput], dict]


@app.get("/", tags=[blup_tag])
def welcome():
    """Welcome endpoint"""
    return jsonify({
        'status': 'success',
        'results': "Beekube BLUP Melifera"
    })


@app.post("/blup", tags=[blup_tag])

def blup(body: BLUPInput):
    try:
        """
        Perform BLUP calculation
        """
        # Make sure that the directory "data" exists
        os.makedirs('./data', exist_ok=True)

        # Delete the file input.json
        if os.path.exists('./data/input.json'):
            os.remove('./data/input.json')

        # Delete the file results_index.json
        if os.path.exists('./data/resultats_index.json'):
            os.remove('./data/resultats_index.json')

        # Save the JSON data in a file
        with open('./data/input.json', 'w') as f:
            json.dump(body.dict(), f)

        # Run the R script
        result = subprocess.run(['Rscript', 'blup.r'], capture_output=True, text=True)

        # Read the results
        if os.path.exists('./data/resultats_index.json'):
            with open('./data/resultats_index.json', 'r') as f:
                results = json.load(f)
        else:
            results = {"error": "Results file not found"}

        # Prepare the response
        response = BLUPOutput(
            status='success' if result.returncode == 0 else 'error',
            results=results
        )

        # Delete the file input.json
        if os.path.exists('./data/input.json'):
            os.remove('./data/input.json')

        # Delete the file resultats_index.json.
        if os.path.exists('./data/resultats_index.json'):
            os.remove('./data/resultats_index.json')

        return jsonify(response.model_dump())
    except Exception as e:
        return jsonify({
            'status': 'error',
            'message': str(e),
            'output': result.stdout if 'result' in locals() else '',
            'error': result.stderr if 'result' in locals() else '',
            'results': {"error": "An unexpected error occurred"}
        }), 500


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8081)
