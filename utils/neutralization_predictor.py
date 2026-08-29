import json
import logging
import os

import joblib
import numpy as np

from utils.paths import COV_UNIBIND_DIR, MODEL_DIR

logger = logging.getLogger(__name__)

# Wild-type IC50 baselines (μg/ml) from CoV-UniBind data
DEFAULT_WILDTYPE_IC50 = {
    "Casirivimab": 0.011725,
    "Imdevimab": 0.015525,
}


class NeutralizationPredictor:
    def __init__(self, model_dir=None):
        if model_dir is None:
            model_dir = MODEL_DIR

        self.model_dir = model_dir
        self.target_antibodies = ["Casirivimab", "Imdevimab"]
        self.models = {}
        self.feature_columns = []
        self.wildtype_ic50 = dict(DEFAULT_WILDTYPE_IC50)
        self.ready = False

        self._load_wildtype_baselines()
        self._load_models()

    def _load_wildtype_baselines(self):
        wt_path = os.path.join(COV_UNIBIND_DIR, "wildtype_ic50.csv")

        if not os.path.exists(wt_path):
            return

        try:
            import csv

            with open(wt_path, newline="", encoding="utf-8") as f:
                for row in csv.DictReader(f):
                    name = row.get("antibody_id", "").strip()
                    ic50 = float(row["ic50_ug_ml"])
                    if name in self.target_antibodies:
                        self.wildtype_ic50[name] = ic50
        except Exception as exc:
            logger.warning("Failed to load wildtype IC50 baselines: %s", exc)

    def _normalize_input_mutation(self, mutation):
        """E484K / S:E484K -> 484K to match training feature names."""
        m = str(mutation).strip().upper().replace("S:", "")
        if len(m) >= 3 and m[0].isalpha() and any(ch.isdigit() for ch in m[1:]):
            m = m[1:]
        return m

    def _load_models(self):
        feature_path = os.path.join(self.model_dir, "feature_columns.json")
        if not os.path.exists(feature_path):
            logger.warning("feature_columns.json not found in %s", self.model_dir)
            return

        try:
            with open(feature_path, "r", encoding="utf-8") as f:
                self.feature_columns = json.load(f)
        except Exception as exc:
            logger.error("Failed to load feature columns: %s", exc)
            return

        for ab in self.target_antibodies:
            model_path = os.path.join(self.model_dir, f"{ab.lower()}_model.pkl")
            if not os.path.exists(model_path):
                logger.warning("Model file not found for %s: %s", ab, model_path)
                continue

            try:
                self.models[ab] = joblib.load(model_path)
                logger.info("Loaded ML model for %s", ab)
            except Exception as exc:
                logger.error("Failed to load model for %s: %s", ab, exc)

        self.ready = len(self.models) > 0 and len(self.feature_columns) > 0
        if self.ready:
            logger.info(
                "ML predictor ready with %d features and %d models",
                len(self.feature_columns),
                len(self.models),
            )
        else:
            logger.warning("ML predictor is not ready")

    @staticmethod
    def models_available(model_dir=None):
        model_dir = model_dir or MODEL_DIR
        feature_path = os.path.join(model_dir, "feature_columns.json")
        if not os.path.exists(feature_path):
            return False
        return all(
            os.path.exists(os.path.join(model_dir, f"{ab.lower()}_model.pkl"))
            for ab in ["Casirivimab", "Imdevimab"]
        )

    def _mutations_to_feature_vector(self, mutations):
        normalized = [self._normalize_input_mutation(m) for m in mutations]
        normalized_set = set(normalized)
        vec = [1.0 if feat in normalized_set else 0.0 for feat in self.feature_columns]
        return np.array(vec, dtype=float).reshape(1, -1), normalized

    def _risk_from_ic50(self, ic50_value):
        if ic50_value >= 0.1:
            return "High"
        if ic50_value >= 0.03:
            return "Moderate"
        return "Low"

    def predict_variant_neutralization(self, mutations):
        if not self.ready:
            return {
                "individual_analysis": {},
                "cocktail_prediction": 0.0,
                "summary_risk": "Unavailable",
                "matched_mutation_count": 0,
                "input_mutations": mutations,
                "error": "ML predictor not initialized. Train models first.",
            }

        X, normalized_mutations = self._mutations_to_feature_vector(mutations)
        known_mutations = [m for m in normalized_mutations if m in set(self.feature_columns)]

        analysis = {}
        predicted_ic50s = []

        for ab in self.target_antibodies:
            model = self.models.get(ab)
            if model is None:
                continue

            wt_ic50 = self.wildtype_ic50.get(ab, 0.01)

            try:
                # Model predicts log10(fold_change) relative to wild-type
                pred_log10_fold = float(model.predict(X)[0])
                pred_log10_fold = float(np.clip(pred_log10_fold, -2.0, 4.0))
                fold_change = float(10 ** pred_log10_fold)
                pred_ic50 = float(wt_ic50 * fold_change)
            except Exception as exc:
                logger.error("Prediction failed for %s: %s", ab, exc)
                fold_change = 1.0
                pred_ic50 = wt_ic50

            predicted_ic50s.append(pred_ic50)

            analysis[ab] = {
                "predicted_ic50_ug_ml": round(pred_ic50, 4),
                "fold_change": round(fold_change, 2),
                "wildtype_ic50_ug_ml": round(wt_ic50, 4),
                "model_type": type(model).__name__,
                "risk_level": self._risk_from_ic50(pred_ic50),
                "recognized_mutations": known_mutations,
            }

        cocktail_prediction = round(min(predicted_ic50s), 4) if predicted_ic50s else 0.0

        if cocktail_prediction >= 0.1:
            summary_risk = "High"
        elif cocktail_prediction >= 0.03:
            summary_risk = "Moderate"
        else:
            summary_risk = "Potentially Effective"

        return {
            "individual_analysis": analysis,
            "cocktail_prediction": cocktail_prediction,
            "summary_risk": summary_risk,
            "matched_mutation_count": len(set(known_mutations)),
            "input_mutations": [str(m) for m in mutations],
            "model_based": True,
        }


_instance = None


def get_neutralization_predictor():
    """Lazy singleton — sklearn model load must happen after SequenceAnalyzer init."""
    global _instance
    if _instance is None:
        _instance = NeutralizationPredictor()
    return _instance


class _LazyPredictorProxy:
    @property
    def ready(self):
        if _instance is not None:
            return _instance.ready
        return NeutralizationPredictor.models_available()

    def predict_variant_neutralization(self, mutations):
        return get_neutralization_predictor().predict_variant_neutralization(mutations)


neutralization_predictor = _LazyPredictorProxy()
