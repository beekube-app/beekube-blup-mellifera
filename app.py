from flask import Flask, request, jsonify
import subprocess
import json
import os

app = Flask(__name__)

@app.get('/')
def welcome():
    # Préparer la réponse
    response = {
        'status': 'success',
        'results': "Beekube BLUP Melifera"
    }
    return jsonify(response)


@app.route('/blup', methods=['POST'])
def blup():
    data = request.json

    # Assurez-vous que le répertoire data existe
    os.makedirs('./data', exist_ok=True)

    # Sauvegarder les données JSON dans un fichier
    with open('/app/data/input.json', 'w') as f:
        json.dump(data, f)

    # Exécuter le script R
    result = subprocess.run(['Rscript', 'blup.r'], capture_output=True, text=True)

    # Lire les résultats
    if os.path.exists('/app/data/resultats_index.json'):
        with open('/app/data/resultats_index.json', 'r') as f:
            results = json.load(f)
    else:
        results = {"error": "Fichier de résultats non trouvé"}

    # Préparer la réponse
    response = {
        'status': 'success' if result.returncode == 0 else 'error',
        'output': result.stdout,
        'error': result.stderr,
        'results': results
    }

    return jsonify(response)


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8081)