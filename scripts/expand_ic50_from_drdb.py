"""
Build processed_ml_features.csv from Stanford CoV-DRDB + Arora baseline.

Uses the website's 95-mutation panel, median-aggregates duplicate isoform
measurements, and writes to data/training/processed_ml_features.csv.

DRDB payload (one-time):
  git clone --depth 1 https://github.com/hivdb/covid-drdb-payload.git \\
    data/training/local_expansion/drdb_payload
"""

from __future__ import annotations

import json
import os
import sys
from glob import glob

import numpy as np
import pandas as pd

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

from scripts.build_ml_features import VARIANT_MUTATIONS, parse_iso_to_mutations
from utils.paths import (
    COV_UNIBIND_DIR,
    MODEL_DIR,
    PROCESSED_FEATURES_CSV,
    PROJECT_ROOT as ROOT,
)

DRDB_ROOT = os.path.join(ROOT, "data", "training", "local_expansion", "drdb_payload")
RX_POTENCY_DIR = os.path.join(DRDB_ROOT, "tables", "rx_potency")
BASELINE_CSV = os.path.join(COV_UNIBIND_DIR, "raw", "arora.txt")
ISOLATES_CSV = os.path.join(DRDB_ROOT, "tables", "isolates.csv")
VARIANT_CONSENSUS_CSV = os.path.join(DRDB_ROOT, "tables", "variant_consensus.csv")
ISOLATE_MUTATIONS_CSV = os.path.join(DRDB_ROOT, "tables", "isolate_mutations.csv")
PANEL_PATH = os.path.join(MODEL_DIR, "feature_columns.json")
OUT_RAW = os.path.join(COV_UNIBIND_DIR, "expanded_raw_potency.csv")

TARGET_ANTIBODIES = {"Casirivimab", "Imdevimab"}
COCKTAIL_RX = {"Casirivimab+Imdevimab"}

RX_NORMALIZE = {
    "Casirivimab/Imdevimab": "Casirivimab+Imdevimab",
    "Casirivimab_Imdevimab": "Casirivimab+Imdevimab",
    "Casirivimab–imdevimab": "Casirivimab+Imdevimab",
    "Casirivimab-imdevimab": "Casirivimab+Imdevimab",
    "Casirivimab + Imdevimab": "Casirivimab+Imdevimab",
    "REGN10933": "Casirivimab",
    "REGN10987": "Imdevimab",
}


def spike_token(position, amino_acid: str) -> str | None:
    aa = str(amino_acid).strip()
    if not aa or aa.lower() == "nan":
        return None
    if aa == "del":
        return f"{int(position)}del"
    if len(aa) == 1 and aa.isalpha():
        return f"{int(position)}{aa.upper()}"
    return None


def build_variant_consensus_map() -> dict[str, list[str]]:
    if not os.path.exists(VARIANT_CONSENSUS_CSV):
        return {}
    vc = pd.read_csv(VARIANT_CONSENSUS_CSV)
    vc = vc[vc["gene"] == "S"]
    out: dict[str, list[str]] = {}
    for var_name, grp in vc.groupby("var_name"):
        tokens = []
        for _, row in grp.iterrows():
            tok = spike_token(row["position"], row["amino_acid"])
            if tok:
                tokens.append(tok)
        out[str(var_name)] = sorted(set(tokens))
    return out


def build_isolate_lookup() -> dict[str, str]:
    iso = pd.read_csv(ISOLATES_CSV)
    return dict(zip(iso["iso_name"].astype(str), iso["var_name"].astype(str)))


def build_isolate_mutation_map() -> dict[str, list[str]]:
    if not os.path.exists(ISOLATE_MUTATIONS_CSV):
        return {}
    im = pd.read_csv(ISOLATE_MUTATIONS_CSV)
    im = im[im["gene"] == "S"]
    out: dict[str, list[str]] = {}
    for iso_name, grp in im.groupby("iso_name"):
        tokens = []
        for _, row in grp.iterrows():
            tok = spike_token(row["position"], row["amino_acid"])
            if tok:
                tokens.append(tok)
        out[str(iso_name)] = sorted(set(tokens))
    return out


def extend_variant_mutations(consensus: dict[str, list[str]], isolates: dict[str, str]) -> dict[str, list[str]]:
    extended = {k: list(v) for k, v in VARIANT_MUTATIONS.items()}
    alias = {
        "B.1.351": "Beta",
        "Omicron/BA.1": "BA.1",
        "Omicron/BA.2": "BA.2",
        "Omicron/BQ.1.1": "BQ.1.1",
        "Omicron/XBB.1": "XBB.1",
        "Omicron/XBB.1.5": "XBB.1.5",
        "Omicron/CH.1.1": "CH.1.1",
        "Omicron/EG.5": "EG.5",
        "KP.2": "KP.2",
        "KP.3": "KP.3",
        "KP.3.1.1": "KP.3.1.1",
        "JN.1": "JN.1",
    }
    for var_name, short in alias.items():
        if var_name in consensus and short not in extended:
            extended[short] = consensus[var_name]
    for iso_name, var_name in isolates.items():
        if " Spike" not in iso_name:
            continue
        base = iso_name.replace(" Spike", "").split(":")[0].strip()
        if base not in extended and var_name in consensus:
            extended[base] = consensus[var_name]
    return extended


def iso_to_mutations(iso_name: str, variant_map: dict[str, list[str]], isolate_mut_map: dict[str, list[str]]) -> list[str]:
    name = str(iso_name).strip()
    if name in isolate_mut_map and isolate_mut_map[name]:
        return isolate_mut_map[name]
    parsed = parse_iso_to_mutations(name)
    if parsed:
        return parsed
    base = name.replace(" Spike", "").split(":")[0].strip()
    if base in variant_map:
        return list(variant_map[base])
    if name.endswith(" Spike"):
        base2 = name[:-6].strip()
        if base2 in variant_map:
            return list(variant_map[base2])
    return []


def load_drdb_potency() -> pd.DataFrame:
    paths = glob(os.path.join(RX_POTENCY_DIR, "*.csv"))
    frames = []
    for path in paths:
        try:
            df = pd.read_csv(path)
        except Exception:
            continue
        if df.empty or "rx_name" not in df.columns:
            continue
        df["source_file"] = os.path.basename(path)
        frames.append(df)
    if not frames:
        raise FileNotFoundError(f"No rx_potency CSV files under {RX_POTENCY_DIR}")
    all_df = pd.concat(frames, ignore_index=True)
    all_df["rx_name"] = all_df["rx_name"].replace(RX_NORMALIZE)
    all_df = all_df[all_df["potency_type"].astype(str).str.upper() == "IC50"]
    all_df = all_df[all_df["potency_unit"].astype(str).str.lower() == "ng/ml"]
    mask = all_df["rx_name"].isin(TARGET_ANTIBODIES | COCKTAIL_RX)
    return all_df[mask].copy()


def load_baseline_potency() -> pd.DataFrame:
    df = pd.read_csv(BASELINE_CSV)
    df["rx_name"] = df["rx_name"].replace(RX_NORMALIZE)
    df["ref_name"] = "AroraEmbedded"
    df["source_file"] = "arora.txt"
    df["potency_type"] = "IC50"
    return df


def clean_potency(row) -> float:
    val = float(row["potency"])
    upper = float(row["potency_upper_limit"])
    return upper if val >= upper else val


def aggregate_by_iso(combined: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (_, _), grp in combined.groupby(["rx_name", "iso_name"], sort=False):
        row = grp.iloc[0].copy()
        row["cleaned_potency"] = float(np.median(grp["cleaned_potency"].astype(float)))
        row["ref_name"] = ";".join(sorted(set(grp["ref_name"].astype(str))))
        row["source_file"] = ";".join(sorted(set(grp["source_file"].astype(str))))
        row["n_measurements"] = len(grp)
        rows.append(row)
    return pd.DataFrame(rows)


def build_feature_matrix(combined: pd.DataFrame, panel: set[str]) -> pd.DataFrame:
    combined = combined.copy()
    combined["mutation_list"] = combined["mutation_list"].apply(lambda muts: [m for m in muts if m in panel])
    combined = combined[combined["mutation_list"].map(len) > 0]
    combined = aggregate_by_iso(combined)

    mutations = sorted(panel)
    rows = []
    for _, row in combined.iterrows():
        feature_dict = {mut: 0 for mut in mutations}
        for m in row["mutation_list"]:
            feature_dict[m] = 1

        rx = row["rx_name"]
        wt_rows = combined[
            (combined["rx_name"] == rx)
            & (combined["iso_name"].isin(["Wuhan-Hu-1", "Wuhan-Hu-1 Spike", "B.1 Spike"]))
        ]
        wt_potency = wt_rows["cleaned_potency"].mean() if not wt_rows.empty else 1.0

        feature_dict["target_y"] = np.log10(row["cleaned_potency"]) - np.log10(wt_potency)
        feature_dict["rx_name"] = rx
        feature_dict["iso_name"] = row["iso_name"]
        rows.append(feature_dict)

    return pd.DataFrame(rows)


def main():
    if not os.path.isdir(RX_POTENCY_DIR):
        raise SystemExit(
            "Missing DRDB payload. Run:\n"
            "  git clone --depth 1 https://github.com/hivdb/covid-drdb-payload.git "
            "data/training/local_expansion/drdb_payload"
        )
    if not os.path.exists(PANEL_PATH):
        raise SystemExit(f"Missing mutation panel: {PANEL_PATH}")

    meta_names = {"ref_name", "source_file", "n_measurements", "target_y", "rx_name", "iso_name"}
    panel = {f for f in json.load(open(PANEL_PATH, encoding="utf-8")) if f not in meta_names}
    consensus = build_variant_consensus_map()
    isolates = build_isolate_lookup()
    isolate_mut = build_isolate_mutation_map()
    variant_map = extend_variant_mutations(consensus, isolates)

    drdb = load_drdb_potency()
    base = load_baseline_potency()
    combined = pd.concat([drdb, base], ignore_index=True, sort=False)
    combined["rx_name"] = combined["rx_name"].replace(RX_NORMALIZE)
    combined["cleaned_potency"] = combined.apply(clean_potency, axis=1)
    combined["mutation_list"] = combined["iso_name"].apply(
        lambda x: iso_to_mutations(x, variant_map, isolate_mut)
    )
    combined = combined[combined["mutation_list"].map(len) > 0].copy()
    combined["dedupe_key"] = combined.apply(
        lambda r: (str(r.get("ref_name", "")), str(r["rx_name"]), str(r["iso_name"]), round(float(r["cleaned_potency"]), 4)),
        axis=1,
    )
    combined = combined.drop_duplicates(subset=["dedupe_key"], keep="first")
    combined.to_csv(OUT_RAW, index=False)

    ml_df = build_feature_matrix(combined, panel)
    ml_df.to_csv(PROCESSED_FEATURES_CSV, index=False)

    print("IC50 feature matrix rebuilt from CoV-DRDB")
    print(f"  Raw potency rows: {len(combined)}")
    print(f"  Feature matrix:   {len(ml_df)} rows -> {PROCESSED_FEATURES_CSV}")
    print(f"  Unique isoforms:  {ml_df['iso_name'].nunique()}")
    print(f"  By antibody:      {ml_df['rx_name'].value_counts().to_dict()}")
    print(f"  Mutation panel:   {len(panel)} features")


if __name__ == "__main__":
    main()
