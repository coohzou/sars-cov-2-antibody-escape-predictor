"""
Run full unbiased evaluation workflow:
1. Rebuild IC50 features from CoV-DRDB
2. Split ML data into train / test
3. Retrain models on train only
4. Evaluate ML, pipeline, and external validation sets
"""

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PY = Path(sys.executable)

STEPS = [
    ("Rebuild IC50 features from CoV-DRDB", ["scripts/expand_ic50_from_drdb.py"]),
    ("Split ML dataset", ["scripts/split_ml_dataset.py"]),
    ("Train on train split only", ["scripts/train_neutralization_model.py"]),
    ("Evaluate ML test set", ["scripts/evaluate_ml_test_set.py"]),
    ("Evaluate pipeline test set", ["scripts/evaluate_pipeline_test_set.py"]),
    ("Evaluate external NCBI + IC50 validation", ["scripts/evaluate_pipeline_external_test.py"]),
]


def main():
    for title, args in STEPS:
        print(f"\n=== {title} ===")
        cmd = [str(PY), *[str(ROOT / a) for a in args]]
        result = subprocess.run(cmd, cwd=ROOT)
        if result.returncode != 0:
            sys.exit(result.returncode)
    print("\nDone. See data/evaluation/ for ml_test_results.json and pipeline_test_results.json")


if __name__ == "__main__":
    main()
