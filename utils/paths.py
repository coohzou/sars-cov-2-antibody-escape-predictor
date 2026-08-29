"""Central path constants for data, training, prediction, and evaluation."""

import os

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

DATA_DIR = os.path.join(PROJECT_ROOT, "data")
TRAINING_DATA_DIR = os.path.join(DATA_DIR, "training")
PREDICTION_DATA_DIR = os.path.join(DATA_DIR, "prediction")
EVALUATION_DIR = os.path.join(DATA_DIR, "evaluation")
EVALUATION_INTERNAL_DIR = os.path.join(EVALUATION_DIR, "internal")

COV_UNIBIND_DIR = os.path.join(TRAINING_DATA_DIR, "cov_unibind")
PROCESSED_FEATURES_CSV = os.path.join(TRAINING_DATA_DIR, "processed_ml_features.csv")
ML_TRAIN_CSV = os.path.join(TRAINING_DATA_DIR, "ml_train.csv")
ML_TEST_CSV = os.path.join(TRAINING_DATA_DIR, "ml_test.csv")
MODEL_DIR = os.path.join(TRAINING_DATA_DIR, "models")

MANIFEST_PATH = os.path.join(PREDICTION_DATA_DIR, "manifest.json")
SPLIT_CONFIG_PATH = os.path.join(EVALUATION_DIR, "split_config.json")


def prediction_fasta(filename: str) -> str:
    return os.path.join(PREDICTION_DATA_DIR, filename)
