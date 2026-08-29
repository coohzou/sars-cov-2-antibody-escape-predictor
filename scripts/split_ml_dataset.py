"""
Split CoV-UniBind ML features into train / test sets.
Test iso_names are defined in data/evaluation/split_config.json and must
never be used during model training.
"""

import json
import os
import sys

import pandas as pd

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

from utils.paths import (
    EVALUATION_DIR,
    ML_TEST_CSV,
    ML_TRAIN_CSV,
    PROCESSED_FEATURES_CSV,
    SPLIT_CONFIG_PATH,
)

META_OUT = os.path.join(EVALUATION_DIR, "split_summary.json")


def main():
    with open(SPLIT_CONFIG_PATH, encoding="utf-8") as f:
        config = json.load(f)

    test_isos = set(config["ml_test_iso_names"])
    df = pd.read_csv(PROCESSED_FEATURES_CSV)

    train_df = df[~df["iso_name"].isin(test_isos)].copy()
    test_df = df[df["iso_name"].isin(test_isos)].copy()

    if test_df.empty:
        print("ERROR: test set is empty. Check split_config.json iso_names.", file=sys.stderr)
        sys.exit(1)
    if train_df.empty:
        print("ERROR: train set is empty.", file=sys.stderr)
        sys.exit(1)

    train_df.to_csv(ML_TRAIN_CSV, index=False)
    test_df.to_csv(ML_TEST_CSV, index=False)

    summary = {
        "source_rows": len(df),
        "train_rows": len(train_df),
        "test_rows": len(test_df),
        "train_iso_names": sorted(train_df["iso_name"].unique().tolist()),
        "test_iso_names": sorted(test_df["iso_name"].unique().tolist()),
        "train_path": os.path.relpath(ML_TRAIN_CSV, PROJECT_ROOT),
        "test_path": os.path.relpath(ML_TEST_CSV, PROJECT_ROOT),
    }
    with open(META_OUT, "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    print(f"Train: {len(train_df)} rows -> {ML_TRAIN_CSV}")
    print(f"Test:  {len(test_df)} rows -> {ML_TEST_CSV}")
    print(f"Held-out iso_names: {sorted(test_df['iso_name'].unique())}")


if __name__ == "__main__":
    main()
