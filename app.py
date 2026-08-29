import logging
import os
import time
import uuid
from collections import defaultdict

from flask import Flask, jsonify, render_template, request
from werkzeug.utils import secure_filename

from utils.sequence_analyzer import SequenceAnalyzer
from utils.neutralization_predictor import neutralization_predictor

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
UPLOAD_FOLDER = os.path.join(BASE_DIR, "uploads")
ALLOWED_EXTENSIONS = {"fasta", "fa", "txt", "seq"}
MAX_CONTENT_LENGTH = 16 * 1024 * 1024  # 16 MB
UPLOAD_RATE_LIMIT = int(os.environ.get("UPLOAD_RATE_LIMIT", "20"))
UPLOAD_RATE_WINDOW_SEC = int(os.environ.get("UPLOAD_RATE_WINDOW_SEC", "60"))

app = Flask(__name__)
app.config["UPLOAD_FOLDER"] = UPLOAD_FOLDER
app.config["MAX_CONTENT_LENGTH"] = MAX_CONTENT_LENGTH
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

_analyzer = None
_upload_timestamps = defaultdict(list)


def get_analyzer():
    global _analyzer
    if _analyzer is None:
        _analyzer = SequenceAnalyzer()
    return _analyzer


def allowed_file(filename):
    return "." in filename and filename.rsplit(".", 1)[1].lower() in ALLOWED_EXTENSIONS


def client_ip():
    forwarded = request.headers.get("X-Forwarded-For", "")
    if forwarded:
        return forwarded.split(",")[0].strip()
    return request.remote_addr or "unknown"


def rate_limit_ok(ip_address):
    now = time.time()
    recent = [stamp for stamp in _upload_timestamps[ip_address] if now - stamp < UPLOAD_RATE_WINDOW_SEC]
    _upload_timestamps[ip_address] = recent
    if len(recent) >= UPLOAD_RATE_LIMIT:
        return False
    _upload_timestamps[ip_address].append(now)
    return True


def build_response(result, detected_mutations):
    raw_prediction = neutralization_predictor.predict_variant_neutralization(detected_mutations)
    details = result.get("details", {})

    payload = {
        "success": True,
        "similarity_score": float(result.get("similarity_score", 100)),
        "sequence_info": {
            "length": int(details.get("sequence_length", 0)),
            "gc_content": float(details.get("gc_content", 0)),
        },
        "mutations": {
            "list": detected_mutations,
            "has_d614g": bool(result.get("has_d614g", False)),
        },
        "variant": {
            "name": result.get("variant"),
            "confidence": result.get("variant_confidence"),
        },
        "neutralization_results": raw_prediction,
    }
    if result.get("warning"):
        payload["warning"] = result["warning"]
    return payload


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/upload", methods=["POST"])
def upload_file():
    if not rate_limit_ok(client_ip()):
        return jsonify({"error": "Too many upload requests. Please wait and try again."}), 429

    if "file" not in request.files:
        return jsonify({"error": "No file selected"}), 400

    file = request.files["file"]
    if not file.filename:
        return jsonify({"error": "No file selected"}), 400

    if not allowed_file(file.filename):
        return jsonify({"error": "File type not allowed. Use .fasta, .fa, .txt, or .seq"}), 400

    safe_name = secure_filename(file.filename)
    filename = os.path.join(app.config["UPLOAD_FOLDER"], f"{uuid.uuid4().hex}_{safe_name}")

    try:
        file.save(filename)
        result = get_analyzer().analyze_sequence_file(filename)

        if not isinstance(result, dict):
            return jsonify({"error": "Analyzer returned invalid data"}), 500

        if not result.get("success", False):
            return jsonify({
                "success": False,
                "error": result.get("error", "Analysis failed"),
                "similarity_score": result.get("similarity_score", 0),
            }), 400

        detected = result.get("detected_mutations", {})
        if isinstance(detected, dict):
            detected_mutations = list(detected.keys())
        elif isinstance(detected, list):
            detected_mutations = detected
        else:
            detected_mutations = []

        logger.info("Detected mutations: %s", detected_mutations)
        return jsonify(build_response(result, detected_mutations))

    except Exception:
        logger.exception("Upload analysis failed")
        return jsonify({"error": "Analysis failed. Please check your FASTA file and try again."}), 500

    finally:
        if os.path.exists(filename):
            os.remove(filename)


@app.route("/health")
def health():
    return jsonify({"status": "ok"}), 200


@app.route("/ready")
def ready():
    analyzer = get_analyzer()
    return jsonify({
        "status": "ok",
        "predictor_ready": neutralization_predictor.ready,
        "references_loaded": len(analyzer.reference_genome or ""),
    })


if __name__ == "__main__":
    debug_mode = os.environ.get("FLASK_DEBUG", "").lower() in {"1", "true", "yes"}
    app.run(host="127.0.0.1", port=5000, debug=debug_mode)
