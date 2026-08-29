import json
import os
import sys

import joblib
import pandas as pd
from sklearn.linear_model import Ridge

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

from utils.paths import ML_TRAIN_CSV, MODEL_DIR, PROCESSED_FEATURES_CSV


def main():
    train_path = ML_TRAIN_CSV
    data_path = train_path if os.path.exists(train_path) else PROCESSED_FEATURES_CSV
    os.makedirs(MODEL_DIR, exist_ok=True)

    if not os.path.exists(data_path):
        raise FileNotFoundError(
            f"Training file not found: {data_path}. "
            "Run scripts/split_ml_dataset.py first to generate ml_train.csv."
        )

    using_holdout = os.path.basename(data_path) == "ml_train.csv"
    print(f"Training data: {data_path} ({'hold-out split' if using_holdout else 'FULL dataset — not recommended'})")

    df = pd.read_csv(data_path)

    meta_cols = ["rx_name", "iso_name", "target_y", "ref_name", "source_file", "n_measurements"]
    mutations = sorted([c for c in df.columns if c not in meta_cols])
    antibodies = sorted(df["rx_name"].unique().tolist())

    print(f"Found {len(antibodies)} antibodies: {antibodies}")
    print(f"Found {len(mutations)} mutation features")

    feature_path = os.path.join(MODEL_DIR, "feature_columns.json")
    with open(feature_path, "w", encoding="utf-8") as f:
        json.dump(mutations, f, indent=2)

    for ab in antibodies:
        ab_df = df[df["rx_name"] == ab]
        X = ab_df[mutations].values
        y = ab_df["target_y"].values

        if len(X) == 0:
            print(f"Skip {ab}: no data")
            continue

        model = Ridge(alpha=1.0)
        model.fit(X, y)

        safe_ab_name = ab.replace("/", "_").replace(" ", "_").replace("–", "_").lower()
        model_path = os.path.join(MODEL_DIR, f"{safe_ab_name}_model.pkl")
        joblib.dump(model, model_path)

        print(f"\nTrained model: {ab}")
        print(f"  samples: {len(X)}")
        print(f"  saved: {model_path}")

        coef_dict = dict(zip(mutations, model.coef_))
        sorted_muts = sorted(coef_dict.items(), key=lambda item: abs(item[1]), reverse=True)
        print("  top mutation weights:")
        for m, weight in sorted_muts[:3]:
            print(f"    - {m}: {weight:+.4f}")

    print("\nAll models saved.")


if __name__ == "__main__":
    main()
