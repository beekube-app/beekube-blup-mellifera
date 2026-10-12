import importlib.util
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OPENAPI = os.path.join(ROOT, "openapi", "openapi.json")


def load_exporter():
    spec = importlib.util.spec_from_file_location(
        "export_openapi", os.path.join(ROOT, "scripts", "export_openapi.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def committed_text():
    with open(OPENAPI, encoding="utf-8") as f:
        return f.read()


def schemas():
    return json.loads(committed_text())["components"]["schemas"]


def test_committed_openapi_is_up_to_date():
    assert committed_text() == load_exporter().render_spec(), \
        "openapi/openapi.json est périmé : lancer .venv/bin/python scripts/export_openapi.py"


def test_v1_contract_field_names_unchanged():
    expected = {
        "BLUPInput": {"exploitation", "evaluate", "evaluate_elimination", "data"},
        "QueenBeeInput": {"queenbee", "queenbee_parent", "drone_parent", "apiary_default",
                          "born", "evaluate", "variety", "died"},
        "Evaluation": {"note", "user", "apiary", "beehiveType", "month", "year"},
        "BLUPOutput": {"status", "results"},
        "BLUPResultOutput": {"blup", "heritabilities"},
        "QueenBeeOutput": {"queenbee", "queenbee_parent", "drone_parent", "apiary_default", "born",
                           "blups", "methods", "variety", "died", "mating_suggestions"},
        "HeritabilityStats": {"heritability", "se", "v_additive_queen", "v_additive_drone",
                              "v_colony", "v_residual"},
        "MatingSuggestion": {"blup_score", "diversity_score", "mate_blup", "queen",
                             "queen_blup", "total_score"},
        "ErrorResponse": {"status", "message", "output", "error", "results"},
    }
    s = schemas()
    for name, fields in expected.items():
        assert set(s[name]["properties"]) == fields, name


def test_heritability_stats_are_nullable():
    for field, schema in schemas()["HeritabilityStats"]["properties"].items():
        assert {"type": "null"} in schema.get("anyOf", []), field


def test_results_is_not_a_free_dict():
    assert schemas()["BLUPOutput"]["properties"]["results"] == {
        "$ref": "#/components/schemas/BLUPResultOutput"}


def test_blups_values_stay_untyped_for_sdk():
    # additionalProperties: true (et non un schéma typé) évite le modèle `BlupsValue`
    # dans le SDK PHP, retiré en 1.0.7 : les valeurs de `blups` restent un `mixed`.
    assert schemas()["QueenBeeOutput"]["properties"]["blups"]["additionalProperties"] is True


def test_operation_id_kept_for_sdk():
    spec = json.loads(committed_text())
    assert spec["paths"]["/blup"]["post"]["operationId"] == "calulate_blup_post"
    assert spec["servers"] == [{"url": "https://mellifera.beekube.com"}]


def test_all_422_responses_use_the_python_312_label():
    # http.HTTPStatus libelle 422 « Unprocessable Content » dès Python 3.13 ; le fichier versionné
    # doit être identique sous 3.12 et 3.14 (matrice CI)
    responses = [operation["responses"]["422"]
                 for path_item in json.loads(committed_text())["paths"].values()
                 for operation in path_item.values() if "422" in operation.get("responses", {})]
    assert responses, "aucune réponse 422 dans openapi.json"
    assert {r["description"] for r in responses} == {"Unprocessable Entity"}
