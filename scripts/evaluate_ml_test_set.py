"""
Evaluate Ridge IC50 models on the held-out ML test set (ml_test.csv).
Reports MAE / RMSE on log10(fold_change) for Casirivimab and Imdevimab.
"""

import json
import os
import sys

import joblib
import numpy as np
import pandas as pd

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

from utils.paths import EVALUATION_DIR, ML_TEST_CSV, MODEL_DIR

TEST_PATH = ML_TEST_CSV
OUT_PATH = os.path.join(EVALUATION_DIR, "ml_test_results.json")

TARGET_ANTIBODIES = ["Casirivimab", "Imdevimab"]


def load_feature_columns():
    path = os.path.join(MODEL_DIR, "feature_columns.json")
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def metrics(y_true, y_pred):
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    mae = float(np.mean(np.abs(y_true - y_pred)))
    rmse = float(np.sqrt(np.mean((y_true - y_pred) ** 2)))
    if len(y_true) > 1 and np.std(y_true) > 0:
        corr = float(np.corrcoef(y_true, y_pred)[0, 1])
    else:
        corr = None
    return {"mae": round(mae, 4), "rmse": round(rmse, 4), "r": round(corr, 4) if corr is not None else None}


def main():
    if not os.path.exists(TEST_PATH):
        print(f"Missing {TEST_PATH}. Run: python scripts/split_ml_dataset.py", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(TEST_PATH)
    mutations = load_feature_columns()
    meta_cols = {"rx_name", "iso_name", "target_y", "ref_name", "source_file", "n_measurements"}

    rows_out = []
    summary_by_ab = {}

    for ab in TARGET_ANTIBODIES:
        model_path = os.path.join(MODEL_DIR, f"{ab.lower()}_model.pkl")
        if not os.path.exists(model_path):
            print(f"Missing model: {model_path}", file=sys.stderr)
            continue

        model = joblib.load(model_path)
        ab_df = df[df["rx_name"] == ab]
        if ab_df.empty:
            continue

        feature_cols = [c for c in mutations if c in ab_df.columns]
        missing = [c for c in mutations if c not in ab_df.columns]
        if missing:
            for col in missing:
                ab_df[col] = 0.0
        X = ab_df[mutations].values
        y_true = ab_df["target_y"].values
        y_pred = model.predict(X)

        for (_, row), pred in zip(ab_df.iterrows(), y_pred):
            rows_out.append({
                "iso_name": row["iso_name"],
                "antibody": ab,
                "target_log10_fold": round(float(row["target_y"]), 4),
                "predicted_log10_fold": round(float(pred), 4),
                "error": round(float(pred - row["target_y"]), 4),
            })

        summary_by_ab[ab] = {
            "n": int(len(ab_df)),
            **metrics(y_true, y_pred),
        }

    overall_true = [r["target_log10_fold"] for r in rows_out]
    overall_pred = [r["predicted_log10_fold"] for r in rows_out]

    result = {
        "description": "Held-out ML test evaluation (iso_names excluded from training)",
        "test_file": os.path.relpath(TEST_PATH, PROJECT_ROOT),
        "model_dir": os.path.relpath(MODEL_DIR, PROJECT_ROOT),
        "per_antibody": summary_by_ab,
        "overall": {"n": len(rows_out), **metrics(overall_true, overall_pred)},
        "predictions": rows_out,
    }

    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    with open(OUT_PATH, "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)

    print("ML test set evaluation (log10 fold-change)")
    print("-" * 50)
    for ab, s in summary_by_ab.items():
        print(f"{ab:14s} n={s['n']:2d}  MAE={s['mae']:.4f}  RMSE={s['rmse']:.4f}  r={s['r']}")
    print(f"Overall       n={result['overall']['n']:2d}  MAE={result['overall']['mae']:.4f}  RMSE={result['overall']['rmse']:.4f}")
    print(f"\nSaved {OUT_PATH}")


if __name__ == "__main__":
    main()
