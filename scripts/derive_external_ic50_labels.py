"""
Derive external-validation IC50 labels from committed CoV-DRDB export.

Reads data/training/cov_unibind/expanded_raw_potency.csv, computes per-isoform
median IC50 (ng/ml) for Casirivimab and Imdevimab, and writes provenance metadata
into data/evaluation/pipeline_external_test.json.

Method matches scripts/expand_ic50_from_drdb.py (median of cleaned_potency per
rx_name + iso_name). Wild-type reference IC50 for log10 fold comparisons uses
B.1 Spike medians from the same table.
"""

from __future__ import annotations

import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
POTENCY_CSV = ROOT / "data" / "training" / "cov_unibind" / "expanded_raw_potency.csv"
CONFIG_PATH = ROOT / "data" / "evaluation" / "pipeline_external_test.json"
PROVENANCE_PATH = ROOT / "data" / "evaluation" / "pipeline_external_test_provenance.json"

TARGET_ANTIBODIES = ("Casirivimab", "Imdevimab")
WT_ISO = "B.1 Spike"
CENSOR_FLOOR_NG_ML = 10000.0


def load_measurements():
    groups: dict[tuple[str, str], list[dict]] = defaultdict(list)
    with POTENCY_CSV.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            rx = row["rx_name"].strip()
            if rx not in TARGET_ANTIBODIES:
                continue
            iso = row["iso_name"].strip()
            val = float(row["cleaned_potency"])
            groups[(rx, iso)].append(
                {
                    "value_ng_ml": val,
                    "ref_name": row["ref_name"].strip(),
                    "source_file": row["source_file"].strip(),
                    "section": row.get("section", "").strip(),
                    "assay_name": row.get("assay_name", "").strip(),
                }
            )
    return groups


def summarize_group(entries: list[dict]) -> dict:
    values = [e["value_ng_ml"] for e in entries]
    median = float(statistics.median(values))
    refs = sorted({e["ref_name"] for e in entries})
    sources = sorted({e["source_file"] for e in entries})
    return {
        "value_ng_ml": round(median, 4),
        "censored": median >= CENSOR_FLOOR_NG_ML,
        "n_measurements": len(entries),
        "drdb_references": refs,
        "drdb_source_files": sources,
        "derivation": "median cleaned IC50 (ng/ml) from expanded_raw_potency.csv",
        "source_table": "data/training/cov_unibind/expanded_raw_potency.csv",
        "drdb_payload": "https://github.com/hivdb/covid-drdb-payload",
    }


def main() -> None:
    groups = load_measurements()
    config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))

    wt_block = {}
    provenance = {
        "description": (
            "CoV-DRDB IC50 labels for external validation. "
            "Medians computed from expanded_raw_potency.csv (CoV-DRDB rx_potency export)."
        ),
        "method": "median(cleaned_potency) grouped by rx_name and iso_name",
        "wild_type_iso_name": WT_ISO,
        "source_table": "data/training/cov_unibind/expanded_raw_potency.csv",
        "drdb_repository": "https://github.com/hivdb/covid-drdb-payload",
        "cases": {},
    }

    for antibody in TARGET_ANTIBODIES:
        wt = summarize_group(groups[(antibody, WT_ISO)])
        wt_block[antibody] = wt["value_ng_ml"]
        provenance["wild_type_ic50_ng_ml"] = provenance.get("wild_type_ic50_ng_ml", {})
        provenance["wild_type_ic50_ng_ml"][antibody] = wt

    config["wt_ic50_ng_ml"] = wt_block
    config["wt_ic50_derivation"] = {
        "iso_name": WT_ISO,
        "method": provenance["method"],
        "source_table": provenance["source_table"],
    }

    for case in config["cases"]:
        iso = case["drdb_iso_name"]
        case_prov = {"drdb_iso_name": iso, "antibodies": {}}
        literature = {}
        for antibody in TARGET_ANTIBODIES:
            key = (antibody, iso)
            if key not in groups:
                raise SystemExit(f"No DRDB measurements for {antibody} / {iso}")
            summary = summarize_group(groups[key])
            literature[antibody] = {
                "value": summary["value_ng_ml"],
                "censored": summary["censored"],
            }
            case_prov["antibodies"][antibody] = summary
        case["literature_ic50_ng_ml"] = literature
        provenance["cases"][case["name"]] = case_prov

    CONFIG_PATH.write_text(json.dumps(config, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    PROVENANCE_PATH.write_text(json.dumps(provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Updated {CONFIG_PATH}")
    print(f"Wrote {PROVENANCE_PATH}")
    print("Wild-type (B.1 Spike) IC50 ng/ml:", wt_block)


if __name__ == "__main__":
    main()
