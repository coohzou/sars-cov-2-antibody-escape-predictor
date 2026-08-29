"""
Evaluate variant identification + mutation detection on held-out NCBI references.

These genomes are NOT used in CoV-UniBind ML training (different lineages).
Use this output for the report — not batch_test_variants.py on BA.1/BQ.1.1 etc.
"""

import json
import os
import sys
import time

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

from utils.paths import EVALUATION_DIR, MANIFEST_PATH, PREDICTION_DATA_DIR, SPLIT_CONFIG_PATH
from utils.sequence_analyzer import SequenceAnalyzer
from utils.neutralization_predictor import neutralization_predictor

SPLIT_CONFIG = SPLIT_CONFIG_PATH
MANIFEST_PATH = MANIFEST_PATH
OUT_JSON = os.path.join(EVALUATION_DIR, "pipeline_test_results.json")
OUT_CSV = os.path.join(EVALUATION_DIR, "pipeline_test_summary.csv")

EXPECTED_MAJOR = {
    "Alpha": "Alpha",
    "Beta": "Beta",
    "Gamma": "Gamma",
    "Delta": "Delta",
    "Omicron": "Omicron",
    "BA.4": "Omicron",
    "JN.1": "Omicron",
}

EXPECTED_MUTATIONS = {
    "Alpha": {"D614G", "N501Y", "P681H"},
    "Beta": {"K417N", "E484K", "N501Y", "D614G"},
    "Gamma": {"D614G", "N501Y", "E484K", "K417T", "H655Y"},
    "Delta": {"D614G", "L452R", "T478K", "P681R"},
    "Omicron": {"D614G", "N501Y", "K417N", "P681H"},
    "BA.4": {"D614G", "N501Y", "K417N", "L452R", "F486V", "P681H"},
    "JN.1": {"D614G", "N501Y", "K417N", "R346T"},
}


def load_json(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def resolve_path(item):
    return os.path.join(PREDICTION_DATA_DIR, item["file"])


def mutation_recall(detected, expected):
    if not expected:
        return 1.0
    hit = len(set(detected) & expected)
    return round(hit / len(expected), 3)


def main():
    config = load_json(SPLIT_CONFIG)
    manifest = load_json(MANIFEST_PATH)
    test_names = set(config["pipeline_test_variants"])
    items = [v for v in manifest["variants"] if v["name"] in test_names]

    analyzer = SequenceAnalyzer()
    rows = []

    print("Pipeline held-out test set (variant ID + mutation detection)\n")
    print(f"{'Variant':<12} {'Match':<12} {'Sim%':>6} {'Mut#':>5} {'Recall':>7} Time")
    print("-" * 60)

    for item in items:
        path = resolve_path(item)
        start = time.time()
        result = analyzer.analyze_sequence_file(path)
        elapsed = round(time.time() - start, 1)

        if not result.get("success"):
            print(f"{item['name']:<12} FAILED")
            continue

        detected = list(result.get("detected_mutations", {}).keys())
        expected = EXPECTED_MUTATIONS.get(item["name"], set())
        recall = mutation_recall(detected, expected)
        prediction = neutralization_predictor.predict_variant_neutralization(detected)

        expected_major = EXPECTED_MAJOR.get(item["name"], item["name"])
        matched = result.get("variant")
        id_correct = matched == expected_major

        row = {
            "name": item["name"],
            "lineage": item.get("lineage", ""),
            "accession": item.get("accession", ""),
            "expected_major_variant": expected_major,
            "matched_variant": matched,
            "identification_correct": id_correct,
            "similarity": result.get("similarity_score"),
            "mutation_count": len(detected),
            "detected_mutations": detected,
            "expected_key_mutations": sorted(expected),
            "mutation_recall": recall,
            "cocktail_ic50": prediction.get("cocktail_prediction"),
            "summary_risk": prediction.get("summary_risk"),
            "elapsed_sec": elapsed,
            "evaluation_split": "pipeline_test",
        }
        rows.append(row)

        ok = "OK" if row["identification_correct"] else "MISS"
        print(
            f"{item['name']:<12} {str(result.get('variant')):<12} "
            f"{result.get('similarity_score', 0):>6.1f} {len(detected):>5} "
            f"{recall:>6.0%} {elapsed}s  [{ok}]"
        )

    id_acc = sum(1 for r in rows if r["identification_correct"]) / len(rows) if rows else 0
    mean_recall = sum(r["mutation_recall"] for r in rows) / len(rows) if rows else 0

    payload = {
        "description": "Held-out pipeline test (lineages not used as CoV-UniBind ML training iso_names)",
        "identification_accuracy": round(id_acc * 100, 1),
        "mean_mutation_recall": round(mean_recall, 3),
        "n": len(rows),
        "results": rows,
    }

    os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
    with open(OUT_JSON, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)

    import csv
    with open(OUT_CSV, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "variant", "lineage", "matched_variant", "identification_correct", "similarity_pct",
            "mutation_count", "key_mutation_recall", "detected_mutations",
            "cocktail_ic50_ug_ml", "risk", "elapsed_sec",
        ])
        for r in rows:
            w.writerow([
                r["name"], r["lineage"], r["matched_variant"], r["identification_correct"],
                r["similarity"], r["mutation_count"], r["mutation_recall"],
                "; ".join(r["detected_mutations"]), r["cocktail_ic50"],
                r["summary_risk"], r["elapsed_sec"],
            ])

    print(f"\nIdentification accuracy: {id_acc:.0%} ({sum(r['identification_correct'] for r in rows)}/{len(rows)})")
    print(f"Mean key-mutation recall: {mean_recall:.0%}")
    print(f"Saved {OUT_JSON}")
    print(f"Saved {OUT_CSV}")


if __name__ == "__main__":
    main()
