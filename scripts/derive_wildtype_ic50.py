"""Regenerate wildtype_ic50.csv from CoV-DRDB B.1 Spike median potencies."""

import os
import sys

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

from utils.wildtype_baseline import (
    EXPANDED_POTENCY_CSV,
    WILDTYPE_IC50_CSV,
    WILDTYPE_IC50_JSON,
    compute_b1_spike_medians_ng_ml,
    write_wildtype_files,
)


def main():
    if not os.path.exists(EXPANDED_POTENCY_CSV):
        raise SystemExit(f"Missing {EXPANDED_POTENCY_CSV}. Run expand_ic50_from_drdb.py first.")

    ng_ml = write_wildtype_files()
    print("Wild-type IC50 baselines (B.1 Spike medians, ng/ml):")
    for ab, val in ng_ml.items():
        print(f"  {ab}: {val} ng/ml ({val / 1000:.6f} ug/ml)")
    print(f"\nWrote {WILDTYPE_IC50_CSV}")
    print(f"Wrote {WILDTYPE_IC50_JSON}")


if __name__ == "__main__":
    main()
