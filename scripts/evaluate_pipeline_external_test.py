"""
Evaluate the full pipeline on NCBI genomes that were never used as website references.

Downloads FASTAs (if missing), runs variant ID + mutation detection + IC50 prediction,
then compares predictions to CoV-DRDB literature IC50 (median ng/ml).
"""

import csv
import json
import math
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

from utils.paths import EVALUATION_DIR
from utils.sequence_analyzer import SequenceAnalyzer
from utils.neutralization_predictor import neutralization_predictor
from utils.wildtype_baseline import load_wildtype_ic50_ng_ml

CONFIG_PATH = os.path.join(EVALUATION_DIR, "pipeline_external_test.json")
EXTERNAL_DIR = os.path.join(EVALUATION_DIR, "external")
OUT_JSON = os.path.join(EVALUATION_DIR, "pipeline_external_test_results.json")
OUT_CSV = os.path.join(EVALUATION_DIR, "pipeline_external_test_summary.csv")
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def load_json(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def fetch_fasta(accession):
    params = urllib.parse.urlencode({
        "db": "nuccore",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text",
    })
    url = f"{EFETCH_URL}?{params}"
    with urllib.request.urlopen(url, timeout=60) as resp:
        return resp.read().decode("utf-8")


def normalize_fasta(text, accession, lineage):
    lines = [line.strip() for line in text.strip().splitlines() if line.strip()]
    if not lines:
        raise ValueError("Empty FASTA response")

    header = lines[0] if lines[0].startswith(">") else f">{accession} SARS-CoV-2 {lineage} complete genome"
    if "complete genome" not in header.lower():
        raise ValueError(f"Refusing partial record (not a complete genome): {header[:120]}")
    sequence = "".join(line for line in lines[1:] if not line.startswith(">")).upper()
    if len(sequence) < 29000:
        raise ValueError(f"Sequence too short ({len(sequence)} bp)")

    wrapped = [header]
    for i in range(0, len(sequence), 70):
        wrapped.append(sequence[i:i + 70])
    return "\n".join(wrapped) + "\n"


def ensure_fasta(case):
    out_path = os.path.join(EXTERNAL_DIR, case["file"])
    if os.path.exists(out_path):
        return out_path

    os.makedirs(EXTERNAL_DIR, exist_ok=True)
    print(f"  Downloading {case['name']} ({case['accession']})...")
    raw = fetch_fasta(case["accession"])
    fasta = normalize_fasta(raw, case["accession"], case["lineage"])
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(fasta)
    time.sleep(0.8)
    return out_path


def log10_fold_from_ng_ml(obs_ng_ml, wt_ng_ml):
    return math.log10(obs_ng_ml) - math.log10(wt_ng_ml)


def ic50_metrics(rows):
    """MAE on log10 fold for uncensored antibody predictions."""
    errors = []
    for row in rows:
        for ab in ("Casirivimab", "Imdevimab"):
            comp = row.get("ic50_comparison", {}).get(ab)
            if not comp or comp.get("censored"):
                continue
            if comp.get("log10_fold_error") is not None:
                errors.append(abs(comp["log10_fold_error"]))
    if not errors:
        return {"n": 0, "mae_log10_fold": None}
    return {"n": len(errors), "mae_log10_fold": round(sum(errors) / len(errors), 4)}


def main():
    config = load_json(CONFIG_PATH)
    wt = load_wildtype_ic50_ng_ml()
    analyzer = SequenceAnalyzer()
    rows = []

    print("External pipeline validation (new NCBI genomes + DRDB IC50 labels)\n")
    print(
        f"{'Case':<12} {'Match':<10} {'Sim%':>6}  "
        f"{'Casi lit':>10} {'Casi pred':>10}  {'Imd lit':>10} {'Imd pred':>10}  ID"
    )
    print("-" * 95)

    for case in config["cases"]:
        path = ensure_fasta(case)
        start = time.time()
        result = analyzer.analyze_sequence_file(path)
        elapsed = round(time.time() - start, 1)

        if not result.get("success"):
            print(f"{case['name']:<12} FAILED")
            continue

        detected = list(result.get("detected_mutations", {}).keys())
        prediction = neutralization_predictor.predict_variant_neutralization(detected)
        matched = result.get("variant")
        id_correct = matched == case["expected_major_variant"]

        ic50_comparison = {}
        for ab in ("Casirivimab", "Imdevimab"):
            lit = case["literature_ic50_ng_ml"][ab]
            obs_ng = float(lit["value"])
            censored = bool(lit.get("censored", False))
            obs_log_fold = log10_fold_from_ng_ml(obs_ng, wt[ab])

            ab_pred = prediction.get("individual_analysis", {}).get(ab, {})
            pred_ug = ab_pred.get("predicted_ic50_ug_ml")
            pred_ng = pred_ug * 1000.0 if pred_ug is not None else None
            pred_log_fold = math.log10(ab_pred["fold_change"]) if ab_pred.get("fold_change") else None

            log_err = None
            if pred_log_fold is not None and not censored:
                log_err = round(pred_log_fold - obs_log_fold, 4)

            ic50_comparison[ab] = {
                "literature_ng_ml": obs_ng,
                "literature_log10_fold": round(obs_log_fold, 4),
                "predicted_ug_ml": pred_ug,
                "predicted_ng_ml": round(pred_ng, 2) if pred_ng is not None else None,
                "predicted_log10_fold": round(pred_log_fold, 4) if pred_log_fold is not None else None,
                "log10_fold_error": log_err,
                "censored": censored,
                "direction_correct": (
                    (pred_ng >= obs_ng) if censored and pred_ng is not None else None
                ),
            }

        lit_casi = case["literature_ic50_ng_ml"]["Casirivimab"]["value"]
        lit_imd = case["literature_ic50_ng_ml"]["Imdevimab"]["value"]
        pred_casi = ic50_comparison["Casirivimab"].get("predicted_ng_ml")
        pred_imd = ic50_comparison["Imdevimab"].get("predicted_ng_ml")
        lit_cocktail = (lit_casi + lit_imd) / 2.0
        pred_cocktail = prediction.get("cocktail_prediction", 0) * 1000

        row = {
            "name": case["name"],
            "lineage": case["lineage"],
            "accession": case["accession"],
            "drdb_iso_name": case["drdb_iso_name"],
            "expected_major_variant": case["expected_major_variant"],
            "matched_variant": matched,
            "identification_correct": id_correct,
            "similarity": result.get("similarity_score"),
            "mutation_count": len(detected),
            "detected_mutations": detected,
            "cocktail_ic50_literature_ng_ml": lit_cocktail,
            "cocktail_ic50_predicted_ng_ml": round(pred_cocktail, 2),
            "summary_risk": prediction.get("summary_risk"),
            "ic50_comparison": ic50_comparison,
            "elapsed_sec": elapsed,
            "evaluation_split": "pipeline_external",
            "notes": case.get("notes", ""),
        }
        rows.append(row)

        ok = "OK" if id_correct else "MISS"
        casi_lit_s = f">{lit_casi:.0f}" if case["literature_ic50_ng_ml"]["Casirivimab"].get("censored") else f"{lit_casi:.1f}"
        imd_lit_s = f">{lit_imd:.0f}" if case["literature_ic50_ng_ml"]["Imdevimab"].get("censored") else f"{lit_imd:.1f}"
        print(
            f"{case['name']:<12} {str(matched):<10} {result.get('similarity_score', 0):>6.1f}  "
            f"{casi_lit_s:>10} {pred_casi or 0:>10.1f}  "
            f"{imd_lit_s:>10} {pred_imd or 0:>10.1f}  [{ok}]"
        )

    id_acc = sum(1 for r in rows if r["identification_correct"]) / len(rows) if rows else 0
    ic50_stats = ic50_metrics(rows)
    censored_n = sum(
        1
        for r in rows
        for ab in ("Casirivimab", "Imdevimab")
        if r["ic50_comparison"][ab]["censored"]
    )
    direction_ok = sum(
        1
        for r in rows
        for ab in ("Casirivimab", "Imdevimab")
        if r["ic50_comparison"][ab]["censored"] and r["ic50_comparison"][ab]["direction_correct"]
    )

    payload = {
        "description": config["description"],
        "identification_accuracy": round(id_acc * 100, 1),
        "ic50_mae_log10_fold_uncensored": ic50_stats,
        "censored_comparisons": {"n": censored_n, "predicted_at_or_above_limit": direction_ok},
        "n": len(rows),
        "results": rows,
    }

    os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
    with open(OUT_JSON, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)

    with open(OUT_CSV, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "case", "lineage", "accession", "expected_major", "matched_variant",
            "identification_correct", "similarity_pct", "mutation_count",
            "casi_lit_ng_ml", "casi_pred_ng_ml", "casi_log10_fold_err",
            "imd_lit_ng_ml", "imd_pred_ng_ml", "imd_log10_fold_err",
            "cocktail_lit_ng_ml", "cocktail_pred_ng_ml", "risk", "elapsed_sec",
        ])
        for r in rows:
            c = r["ic50_comparison"]
            w.writerow([
                r["name"], r["lineage"], r["accession"], r["expected_major_variant"],
                r["matched_variant"], r["identification_correct"], r["similarity"],
                r["mutation_count"],
                c["Casirivimab"]["literature_ng_ml"], c["Casirivimab"]["predicted_ng_ml"],
                c["Casirivimab"]["log10_fold_error"],
                c["Imdevimab"]["literature_ng_ml"], c["Imdevimab"]["predicted_ng_ml"],
                c["Imdevimab"]["log10_fold_error"],
                r["cocktail_ic50_literature_ng_ml"], r["cocktail_ic50_predicted_ng_ml"],
                r["summary_risk"], r["elapsed_sec"],
            ])

    print(f"\nIdentification accuracy: {id_acc:.0%} ({sum(r['identification_correct'] for r in rows)}/{len(rows)})")
    if ic50_stats["mae_log10_fold"] is not None:
        print(
            f"IC50 MAE (log10 fold, uncensored): {ic50_stats['mae_log10_fold']:.3f} "
            f"(n={ic50_stats['n']} antibody measurements)"
        )
    if censored_n:
        print(f"Censored IC50 (at assay limit): {direction_ok}/{censored_n} predicted at/above literature floor")
    print(f"Saved {OUT_JSON}")
    print(f"Saved {OUT_CSV}")


if __name__ == "__main__":
    main()
