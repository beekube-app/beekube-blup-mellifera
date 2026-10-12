"""Entrées et sorties minimales du contrat v1 (POST /blup), partagées par les tests."""


def _record(note):
    return {"note": note, "user": 7, "apiary": 100, "beehiveType": "dadant-10", "month": "04", "year": "24"}


def minimal_input():
    """Entrée valide : 2 reines, critère 10 noté, critère éliminatoire 20 vide (comme Beekube l'envoie)."""
    return {
        "exploitation": 1,
        "evaluate": [10],
        "evaluate_elimination": [20],
        "data": [
            {"queenbee": 1, "queenbee_parent": None, "drone_parent": None, "apiary_default": 100,
             "born": "2023", "variety": "Buckfast", "died": False,
             "evaluate": {"10": [_record(0.5)], "20": []}},
            {"queenbee": 2, "queenbee_parent": 1, "drone_parent": None, "apiary_default": 100,
             "born": "2023", "variety": "Buckfast", "died": False,
             "evaluate": {"10": [_record(0.75)], "20": []}},
        ],
    }


def valid_output():
    """Sortie de blup.r conforme à BLUPResultOutput (se = None : item 2 de la ROADMAP)."""
    return {
        "blup": [
            {"queenbee": "1", "queenbee_parent": "0", "drone_parent": "0", "apiary_default": 100,
             "born": "2023", "blups": {"10": 0.0125}, "methods": {"10": "BLUP"},
             "variety": "Buckfast", "died": False,
             "mating_suggestions": [
                 {"queen": "2", "total_score": 0.61, "blup_score": 0.0, "diversity_score": 0.5,
                  "queen_blup": 0.0125, "mate_blup": -0.0125},
             ]},
            {"queenbee": "2", "queenbee_parent": "1", "drone_parent": "0", "apiary_default": 100,
             "born": "2023", "blups": {"10": -0.0125}, "methods": {"10": "BLUP"},
             "variety": "Buckfast", "died": False, "mating_suggestions": []},
        ],
        "heritabilities": {
            "10": {"heritability": 0.31, "se": None, "v_additive_queen": 0.01,
                   "v_additive_drone": 0.0, "v_colony": 0.002, "v_residual": 0.02},
        },
    }
