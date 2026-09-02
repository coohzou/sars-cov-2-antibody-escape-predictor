"""Shared wild-type IC50 baselines derived from CoV-DRDB B.1 Spike medians."""

from __future__ import annotations

import csv
import json
import os
from statistics import median
from typing import Dict

from utils.paths import COV_UNIBIND_DIR

WILDTYPE_IC50_CSV = os.path.join(COV_UNIBIND_DIR, "wildtype_ic50.csv")
WILDTYPE_IC50_JSON = os.path.join(COV_UNIBIND_DIR, "wildtype_ic50_provenance.json")
EXPANDED_POTENCY_CSV = os.path.join(COV_UNIBIND_DIR, "expanded_raw_potency.csv")

TARGET_ANTIBODIES = ("Casirivimab", "Imdevimab")
REFERENCE_ISOFORM = "B.1 Spike"

# Fallback matches median B.1 Spike potencies in expanded_raw_potency.csv (ng/ml).
DEFAULT_WT_NG_ML: Dict[str, float] = {
    "Casirivimab": 8.35,
    "Imdevimab": 7.5,
}


def compute_b1_spike_medians_ng_ml(potency_csv: str | None = None) -> Dict[str, float]:
    """Median cleaned IC50 (ng/ml) per antibody for iso_name B.1 Spike."""
    path = potency_csv or EXPANDED_POTENCY_CSV
    if not os.path.exists(path):
        return dict(DEFAULT_WT_NG_ML)

    values: Dict[str, list[float]] = {ab: [] for ab in TARGET_ANTIBODIES}
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row.get("iso_name") != REFERENCE_ISOFORM:
                continue
            ab = row.get("rx_name", "")
            if ab not in TARGET_ANTIBODIES:
                continue
            values[ab].append(float(row["cleaned_potency"]))

    out: Dict[str, float] = {}
    for ab in TARGET_ANTIBODIES:
        if values[ab]:
            out[ab] = round(median(values[ab]), 4)
        else:
            out[ab] = DEFAULT_WT_NG_ML[ab]
    return out


def write_wildtype_files(ng_ml: Dict[str, float] | None = None) -> Dict[str, float]:
    """Write wildtype_ic50.csv and provenance JSON; return ng/ml baselines."""
    ng_ml = ng_ml or compute_b1_spike_medians_ng_ml()
    os.makedirs(COV_UNIBIND_DIR, exist_ok=True)

    with open(WILDTYPE_IC50_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["antibody_id", "ic50_ug_ml"])
        writer.writeheader()
        for ab in TARGET_ANTIBODIES:
            writer.writerow({
                "antibody_id": ab,
                "ic50_ug_ml": round(ng_ml[ab] / 1000.0, 6),
            })

    provenance = {
        "description": (
            "Wild-type IC50 baselines for inference and external validation. "
            "Median cleaned potency (ng/ml) for CoV-DRDB isoform B.1 Spike per antibody."
        ),
        "source_file": EXPANDED_POTENCY_CSV,
        "iso_name": REFERENCE_ISOFORM,
        "ic50_ng_ml": ng_ml,
        "ic50_ug_ml": {ab: round(v / 1000.0, 6) for ab, v in ng_ml.items()},
    }
    with open(WILDTYPE_IC50_JSON, "w", encoding="utf-8") as f:
        json.dump(provenance, f, indent=2)
        f.write("\n")

    return ng_ml


def load_wildtype_ic50_ng_ml() -> Dict[str, float]:
    """Load canonical wild-type IC50 baselines in ng/ml."""
    if os.path.exists(WILDTYPE_IC50_CSV):
        out: Dict[str, float] = {}
        with open(WILDTYPE_IC50_CSV, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                name = row.get("antibody_id", "").strip()
                if name in TARGET_ANTIBODIES:
                    out[name] = round(float(row["ic50_ug_ml"]) * 1000.0, 4)
        if len(out) == len(TARGET_ANTIBODIES):
            return out

    return dict(DEFAULT_WT_NG_ML)


def load_wildtype_ic50_ug_ml() -> Dict[str, float]:
    ng = load_wildtype_ic50_ng_ml()
    return {ab: round(v / 1000.0, 6) for ab, v in ng.items()}
