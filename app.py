from flask import jsonify, make_response, request
from flask_openapi3 import OpenAPI, Info, Tag, Server
from flask_openapi3.utils import make_validation_error_response
from pydantic import BaseModel, Field, RootModel, ValidationError
import subprocess
import tempfile
import json
import os
from typing import Any, Dict, List, Optional

# Récupérer l'URL de base depuis les variables d'environnement
API_BASE_URL = os.getenv('API_BASE_URL', 'http://localhost:8081')
APP_DIR = os.path.dirname(os.path.abspath(__file__))
BLUP_SCRIPT = os.path.join(APP_DIR, 'blup.r')
# Doit rester strictement inférieur au --timeout de gunicorn (180 s, Dockerfile-Python)
RSCRIPT_TIMEOUT = int(os.getenv('RSCRIPT_TIMEOUT', '170'))

info = Info(title="Beekube BLUP Mellifera API", version="2.2.0", description="SDK PHP pour l'API BLUP de Beekube - Calcul génétique des abeilles mellifères")
servers = [
    Server(url=API_BASE_URL),
]
def validation_error(e):
    """422 au format v2 ({"status": "error", "errors": [...]}) sur /v2/*, format v1 inchangé ailleurs."""
    if not request.path.startswith("/v2/"):
        return make_validation_error_response(e)
    response = make_response(jsonify({"status": "error", "errors": json.loads(e.json(include_url=False, include_input=False))}))
    response.status_code = 422
    return response


app = OpenAPI(__name__, info=info, servers=servers, validation_error_callback=validation_error)
# Une exploitation réelle pèse environ 1 Ko par reine : 50 Mo laissent une marge très large
# tout en bornant ce qu'un corps de requête peut allouer en mémoire (413 au-delà).
app.config['MAX_CONTENT_LENGTH'] = 50 * 1024 * 1024

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
    variety: str
    died: bool = False

class BLUPInput(BaseModel):
    exploitation: int
    evaluate: List[int]
    evaluate_elimination: List[int]
    data: List[QueenBeeInput]

# Output Models
class HeritabilityStats(BaseModel):
    # Nullables : NA côté R (pas assez de données, échec du modèle ; se n'est jamais estimée)
    heritability: Optional[float]
    se: Optional[float]
    v_additive_queen: Optional[float]
    v_additive_drone: Optional[float]
    v_colony: Optional[float]
    v_residual: Optional[float]

class MatingSuggestion(BaseModel):
    blup_score: float
    diversity_score: float
    mate_blup: float
    queen: str
    queen_blup: float
    total_score: float

class QueenBeeOutput(BaseModel):
    queenbee: str
    queenbee_parent: str
    drone_parent: str
    apiary_default: Optional[int]
    born: str
    blups: Dict[str, Optional[float]] = Field(description="BLUP values for different criteria", json_schema_extra={'additionalProperties': 'true'})
    methods: Dict[str, Optional[str]] = Field(description="Methods used for calculations")
    variety: str
    died: bool = False
    mating_suggestions: Optional[List[MatingSuggestion]] = None


class BLUPResultOutput(BaseModel):
    blup: List[QueenBeeOutput]
    heritabilities: Dict[str, HeritabilityStats]

class BLUPOutput(BaseModel):
    status: str
    results: BLUPResultOutput


class WelcomeResponse(BaseModel):
    status: str
    results: str

class ErrorResponse(BaseModel):
    status: str
    message: str
    output: str
    error: str
    results: dict


class BlupScriptError(Exception):
    """Échec de blup.r : code retour non nul, timeout ou fichier de sortie absent."""

    def __init__(self, message: str, stdout: str = '', stderr: str = ''):
        super().__init__(message)
        self.stdout = stdout
        self.stderr = stderr


def _as_text(value: Any) -> str:
    if value is None:
        return ''
    if isinstance(value, bytes):
        return value.decode('utf-8', errors='replace')
    return value


def run_blup_script(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Lance blup.r dans un répertoire temporaire propre à la requête et renvoie son JSON."""
    with tempfile.TemporaryDirectory(prefix='blup-') as workdir:
        input_path = os.path.join(workdir, 'input.json')
        output_path = os.path.join(workdir, 'output.json')
        with open(input_path, 'w') as f:
            json.dump(payload, f)

        try:
            result = subprocess.run(
                ['Rscript', BLUP_SCRIPT, input_path, output_path],
                capture_output=True, text=True, timeout=RSCRIPT_TIMEOUT, cwd=APP_DIR,
            )
        except subprocess.TimeoutExpired as e:
            raise BlupScriptError(f'Rscript timed out after {RSCRIPT_TIMEOUT} s',
                                  _as_text(e.stdout), _as_text(e.stderr))

        if result.returncode != 0:
            raise BlupScriptError(f'Rscript exited with code {result.returncode}',
                                  result.stdout, result.stderr)
        if not os.path.exists(output_path):
            raise BlupScriptError('Results file not found', result.stdout, result.stderr)
        if result.stderr:
            app.logger.warning('blup.r stderr (last 2000 chars): %s', result.stderr[-2000:])

        with open(output_path) as f:
            return json.load(f)


def error_response(message: str, stdout: str = '', stderr: str = ''):
    body = ErrorResponse(status='error', message=message, output=stdout or '',
                         error=stderr or '', results={'error': message})
    return jsonify(body.model_dump()), 500


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
    if not body.data:
        empty = BLUPOutput(status='success', results=BLUPResultOutput(blup=[], heritabilities={}))
        return jsonify(empty.model_dump())

    try:
        raw_results = run_blup_script(body.model_dump())
    except BlupScriptError as e:
        return error_response(str(e), e.stdout, e.stderr)
    except Exception as e:
        return error_response(f'Unexpected error: {e}')

    try:
        results = BLUPResultOutput.model_validate(raw_results)
    except ValidationError as e:
        return error_response(f'Invalid output from blup.r: {e}')

    return jsonify(BLUPOutput(status='success', results=results).model_dump())

# ------------------------------------------------------------------ v2 (moteur Python blup2)
from blup2.engine import run_blup, run_matings  # noqa: E402
from blup2.schemas import BlupV2Input, BlupV2Output, ErrorV2Out, MatingV2Input, MatingV2Output  # noqa: E402

# Même tag que v1 (`blup_tag`) : le SDK PHP génère une seule classe BlupApi (plan 04, tâche A1).
@app.post("/v2/blup", tags=[blup_tag], operation_id="blup_v2_post", responses={200: BlupV2Output, 422: ErrorV2Out, 500: ErrorV2Out})
def blup_v2(body: BlupV2Input):
    """Évaluation génétique v2 : index par critère, index total, accouplements suggérés."""
    try:
        return jsonify(run_blup(body).model_dump(mode="json"))
    except Exception as e:  # filet de sécurité : jamais de 200 avec status = error (spec §5.2)
        app.logger.exception("erreur inattendue dans /v2/blup")
        return jsonify(ErrorV2Out(message=f"{type(e).__name__}: {e}").model_dump()), 500


@app.post("/v2/mating", tags=[blup_tag], operation_id="mating_v2_post", responses={200: MatingV2Output, 422: ErrorV2Out, 500: ErrorV2Out})
def mating_v2(body: MatingV2Input):
    """Évalue des couples mère × source de mâles choisis (index attendu, consanguinité attendue)."""
    try:
        return jsonify(run_matings(body).model_dump(mode="json"))
    except Exception as e:  # filet de sécurité
        app.logger.exception("erreur inattendue dans /v2/mating")
        return jsonify(ErrorV2Out(message=f"{type(e).__name__}: {e}").model_dump()), 500


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8081)
