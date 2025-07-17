from flask import Flask, request, render_template, jsonify
from flask_cors import CORS
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address
from Lindel.Predictor import gen_prediction, format_predictions
import pickle as pkl
import os
import re

app = Flask(__name__)
CORS(app, supports_credentials=True, resources={r"/*": {"origins": "*"}})

# Set up rate limiting
limiter = Limiter(
    get_remote_address,
    app=app,
    default_limits=["200 per day", "50 per hour"]
)

# Define the path to the Lindel folder
lindel_folder = os.path.join(os.path.dirname(__file__), 'Lindel')

# Load model weights and prerequisites
try:
    with open(os.path.join(lindel_folder, "Model_weights.pkl"), 'rb') as weights_file:
        weights = pkl.load(weights_file)
    with open(os.path.join(lindel_folder, 'model_prereq.pkl'), 'rb') as prereq_file:
        prerequesites = pkl.load(prereq_file)
except FileNotFoundError as e:
    raise RuntimeError(f"Required file not found: {e}")

@app.route("/predict", methods=['POST'])
@limiter.limit("10 per minute")
def predict():
    """Handles prediction requests with rate limiting and input validation.""" 
    seq = request.form.get('string')
    if not seq:
        return jsonify({'error': 'No sequence provided'}), 400
    if not re.match(r"^[ATCG]{60}$", seq):
        return jsonify({'error': 'Invalid sequence. Must be 60 characters and only contain A, T, C, or G.'}), 400

    try:
        y_hat, fs = gen_prediction(seq, weights, prerequesites)
    except ValueError as e:
        return jsonify({'error': str(e)}), 400

    rev_index = prerequesites[1]
    pred_freq = {rev_index[i]: y_hat[i] for i in range(len(y_hat)) if y_hat[i] != 0}
    pred_sorted = sorted(pred_freq.items(), key=lambda kv: kv[1], reverse=True)

    return jsonify(format_predictions(seq, pred_sorted, pred_freq, output_type='json'))

@app.route('/')
def static_page():
    """Serves the main application page."""
    return render_template('index_Lindel.html')

if __name__ == "__main__":
    app.run(port=5000)