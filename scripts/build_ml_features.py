"""
Rebuild ML feature matrix from CoV-UniBind / CoV-DRDB sources.

Preferred entry point:
  python scripts/expand_ic50_from_drdb.py

This legacy script still exports the original embedded Arora table (~83 rows)
if you need to reproduce the first report baseline.
"""

import io
import json
import os
import re
import sys

import joblib
import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge

raw_data = """rx_name,iso_name,potency,potency_upper_limit,potency_lower_limit,potency_unit
Casirivimab,Wuhan-Hu-1,0.65,10000,0.01,ng/ml
Casirivimab,S:19R+142D+156del+157del+158G+452R+478K+614G+681R+950N,0.48,10000,0.01,ng/ml
Imdevimab,Wuhan-Hu-1,8.35,10000,0.01,ng/ml
Imdevimab,S:19R+142D+156del+157del+158G+452R+478K+614G+681R+950N,4.35,10000,0.01,ng/ml
Casirivimab,Wuhan-Hu-1,22.8,10000,0.01,ng/ml
Casirivimab,S:19R+142D+156del+157del+158G+452R+478K+614G+681R+950N,11.5,10000,0.01,ng/ml
Casirivimab,S:19R+142D+156del+157del+158G+417N+452R+478K+614G+681R+950N,9.2,10000,0.01,ng/ml
Imdevimab,Wuhan-Hu-1,22.7,10000,0.01,ng/ml
Imdevimab,S:19R+142D+156del+157del+158G+452R+478K+614G+681R+950N,33.5,10000,0.01,ng/ml
Imdevimab,S:19R+142D+156del+157del+158G+417N+452R+478K+614G+681R+950N,13.7,10000,0.01,ng/ml
Casirivimab,S:614G,19.15,10000,1,ng/ml
Imdevimab,S:614G,23.73,10000,1,ng/ml
Casirivimab/Imdevimab,S:614G,7.605,10000,1,ng/ml
Casirivimab,BA.1 Spike,2266,10000,1,ng/ml
Imdevimab,BA.1 Spike,10000,10000,1,ng/ml
Casirivimab/Imdevimab,BA.1 Spike,2372,10000,1,ng/ml
Casirivimab,BA.2 Spike,1998,10000,1,ng/ml
Imdevimab,BA.2 Spike,1166,10000,1,ng/ml
Casirivimab/Imdevimab,BA.2 Spike,485.9,10000,1,ng/ml
Casirivimab,BA.3 Spike,3475,10000,1,ng/ml
Imdevimab,BA.3 Spike,10000,10000,1,ng/ml
Casirivimab/Imdevimab,BA.3 Spike,5092,10000,1,ng/ml
Casirivimab,B.1 Spike,9.7,2000,1,ng/ml
Imdevimab,B.1 Spike,10.5,2000,1,ng/ml
Casirivimab_Imdevimab,B.1 Spike,11.1,2000,1,ng/ml
Casirivimab,C.1.2 Spike:-144del-243del-244del-859N,71.6,2000,1,ng/ml
Imdevimab,C.1.2 Spike:-144del-243del-244del-859N,8.6,2000,1,ng/ml
Casirivimab_Imdevimab,C.1.2 Spike:-144del-243del-244del-859N,9.0,2000,1,ng/ml
Casirivimab,C.1.2 Spike:25L+152R+478K+879T,29.5,2000,1,ng/ml
Imdevimab,C.1.2 Spike:25L+152R+478K+879T,7.3,2000,1,ng/ml
Casirivimab_Imdevimab,C.1.2 Spike:25L+152R+478K+879T,6.7,2000,1,ng/ml
Casirivimab,S:95I+144S+145N+346K+484K+501Y+614G+681H+950N,58.6,2000,1,ng/ml
Imdevimab,S:95I+144S+145N+346K+484K+501Y+614G+681H+950N,10.1,2000,1,ng/ml
Casirivimab_Imdevimab,S:95I+144S+145N+346K+484K+501Y+614G+681H+950N,18.7,2000,1,ng/ml
Casirivimab,B.1.351 Spike,220.6,2000,1,ng/ml
Imdevimab,B.1.351 Spike,5.8,2000,1,ng/ml
Casirivimab_Imdevimab,B.1.351 Spike,17.8,2000,1,ng/ml
Casirivimab,S:614G,3.6,10000,1,ng/ml
Casirivimab,B.1.640.2 Spike,10.0,10000,1,ng/ml
Imdevimab,S:614G,4.9,10000,1,ng/ml
Imdevimab,B.1.640.2 Spike,2.1,10000,1,ng/ml
Casirivimab,S:614G,12.1,5000,1,ng/ml
Imdevimab,S:614G,11.8,5000,1,ng/ml
Casirivimab+Imdevimab,S:614G,6.6,5000,1,ng/ml
Casirivimab,BA.1 Spike,1456,5000,1,ng/ml
Imdevimab,BA.1 Spike,5000,5000,1,ng/ml
Casirivimab+Imdevimab,BA.1 Spike,1469,5000,1,ng/ml
Casirivimab,BA.2 Spike,938,5000,1,ng/ml
Imdevimab,BA.2 Spike,415,5000,1,ng/ml
Casirivimab+Imdevimab,BA.2 Spike,267,5000,1,ng/ml
Casirivimab,BA.2.12.1 Spike,480,5000,1,ng/ml
Imdevimab,BA.2.12.1 Spike,345,5000,1,ng/ml
Casirivimab+Imdevimab,BA.2.12.1 Spike,293,5000,1,ng/ml
Casirivimab,BA.4/5 Spike,5000,5000,1,ng/ml
Imdevimab,BA.4/5 Spike,844,5000,1,ng/ml
Casirivimab+Imdevimab,BA.4/5 Spike,966,5000,1,ng/ml
Casirivimab,B.1 Spike,21,50000,1,ng/ml
Imdevimab,B.1 Spike,19,50000,1,ng/ml
Casirivimab–imdevimab,B.1 Spike,9,50000,1,ng/ml
Casirivimab,BA.1 Spike,1890,50000,1,ng/ml
Imdevimab,BA.1 Spike,50000,50000,1,ng/ml
Casirivimab–imdevimab,BA.1 Spike,3642,50000,1,ng/ml
Casirivimab,BA.4/5 Spike,50000,50000,1,ng/ml
Imdevimab,BA.4/5 Spike,994,50000,1,ng/ml
Casirivimab–imdevimab,BA.4/5 Spike,2611,50000,1,ng/ml
Casirivimab,BA.4.6 Spike,50000,50000,1,ng/ml
Imdevimab,BA.4.6 Spike,2109,50000,1,ng/ml
Casirivimab–imdevimab,BA.4.6 Spike,5395,50000,1,ng/ml
Casirivimab,BA.2.75.2 Spike,50000,50000,1,ng/ml
Imdevimab,BA.2.75.2 Spike,50000,50000,1,ng/ml
Casirivimab–imdevimab,BA.2.75.2 Spike,50000,50000,1,ng/ml
Casirivimab,BJ.1 Spike,880,50000,1,ng/ml
Imdevimab,BJ.1 Spike,50000,50000,1,ng/ml
Casirivimab–imdevimab,BJ.1 Spike,2456,50000,1,ng/ml
Casirivimab,BQ.1.1 Spike,50000,50000,1,ng/ml
Imdevimab,BQ.1.1 Spike,50000,50000,1,ng/ml
Casirivimab–imdevimab,BQ.1.1 Spike,50000,50000,1,ng/ml
Casirivimab,B.1 Spike,7,50000,1,ng/ml
Imdevimab,B.1 Spike,7,50000,1,ng/ml
Casirivimab-imdevimab,B.1 Spike,6,50000,1,ng/ml
Casirivimab,XBB.1 Spike,50000,50000,1,ng/ml
Imdevimab,XBB.1 Spike,50000,50000,1,ng/ml
Casirivimab-imdevimab,XBB.1 Spike,50000,50000,1,ng/ml
"""

# Variant isoform to mutation token mapping
VARIANT_MUTATIONS = {
    "B.1": ["614G"], "S:614G": ["614G"],
    "B.1.351": ["18F", "80A", "215G", "241del", "242del", "243del", "417N", "484K", "501Y", "614G", "701V"],
    "B.1.640.2": ["42V", "94G", "144del", "190S", "446S", "490S", "501Y", "572I", "614G", "655Y", "679K", "681H",
                  "859N", "1139H"],
    "C.1.2": ["9L", "215G", "449H", "484K", "501Y", "614G", "655Y", "716I", "859N"],
    "BA.1": ["67V", "69del", "70del", "95I", "142D", "143del", "144del", "145del", "211del", "212I", "214EPE", "339D",
             "371L", "373P", "375F", "417N", "440K", "446S", "477N", "478K", "484A", "493R", "496S", "498R", "501Y",
             "505H", "547K", "614G", "655Y", "679K", "681H", "764K", "796Y", "856K", "954H", "969K", "981F"],
    "BA.2": ["19I", "24del", "25del", "26del", "27S", "142D", "213G", "339D", "371F", "373P", "375F", "376A", "405N",
             "408S", "417N", "440K", "477N", "478K", "484A", "493R", "498R", "501Y", "505H", "614G", "655Y", "679K",
             "681H", "764K", "796Y", "954H", "969K"],
    "BA.3": ["67V", "69del", "70del", "95I", "142D", "143del", "144del", "145del", "211del", "212I", "339D", "371F",
             "373P", "375F", "405N", "417N", "440K", "446S", "477N", "478K", "484A", "493R", "498R", "501Y", "505H",
             "614G", "655Y", "679K", "681H", "764K", "796Y", "954H", "969K"],
    "BA.2.12.1": ["19I", "24del", "25del", "26del", "27S", "142D", "213G", "339D", "371F", "373P", "375F", "376A",
                  "405N", "408S", "417N", "440K", "452Q", "477N", "478K", "484A", "493R", "498R", "501Y", "505H",
                  "614G", "655Y", "679K", "681H", "704L", "764K", "796Y", "954H", "969K"],
    "BA.4/5": ["19I", "24del", "25del", "26del", "27S", "69del", "70del", "142D", "213G", "339D", "371F", "373P",
               "375F", "376A", "405N", "408S", "417N", "440K", "452R", "477N", "478K", "484A", "486V", "498R", "501Y",
               "505H", "614G", "655Y", "679K", "681H", "764K", "796Y", "954H", "969K"],
    "BA.4.6": ["19I", "24del", "25del", "26del", "27S", "69del", "70del", "142D", "213G", "346T", "339D", "371F",
               "373P", "375F", "376A", "405N", "408S", "417N", "440K", "452R", "477N", "478K", "484A", "486V", "498R",
               "501Y", "505H", "614G", "655Y", "679K", "681H", "764K", "796Y", "954H", "969K"],
    "BA.2.75.2": ["19I", "24del", "25del", "26del", "27S", "147E", "152R", "157L", "210V", "213G", "257S", "346T",
                  "339D", "371F", "373P", "375F", "376A", "405N", "408S", "417N", "440K", "446S", "460K", "477N",
                  "478K", "484A", "486S", "498R", "501Y", "505H", "614G", "655Y", "679K", "681H", "764K", "796Y",
                  "954H", "969K", "1199N"],
    "BJ.1": ["19I", "24del", "25del", "26del", "27S", "83A", "142D", "147E", "152R", "157L", "210V", "213G", "257S",
             "346T", "445P", "446S", "460K", "484V", "614G", "655Y", "679K", "681H", "764K", "796Y", "954H", "969K"],
    "BQ.1.1": ["19I", "24del", "25del", "26del", "27S", "69del", "70del", "142D", "213G", "346T", "339D", "371F",
               "373P", "375F", "376A", "405N", "408S", "417N", "440K", "444T", "452R", "460K", "477N", "478K", "484A",
               "486V", "498R", "501Y", "505H", "614G", "655Y", "679K", "681H", "764K", "796Y", "954H", "969K"],
    "XBB.1": ["19I", "24del", "25del", "26del", "27S", "142D", "213G", "346T", "368I", "339D", "371F", "373P", "375F",
              "376A", "405N", "408S", "417N", "440K", "445P", "446S", "460K", "477N", "478K", "484A", "486S", "490S",
              "498R", "501Y", "505H", "614G", "655Y", "679K", "681H", "764K", "796Y", "954H", "969K"]
}


def parse_iso_to_mutations(iso_name):
    name_str = str(iso_name).strip()
    if name_str in ["Wuhan-Hu-1", "Wuhan-Hu-1 Spike"]: return []
    mut_set = set()
    base_part, inline_part = name_str, ""
    if ":" in name_str: base_part, inline_part = name_str.split(":", 1)

    clean_base = base_part.replace(" Spike", "").strip()
    if clean_base in VARIANT_MUTATIONS:
        mut_set.update(VARIANT_MUTATIONS[clean_base])
    elif "S:" in base_part:
        inline_part = base_part

    target_text = inline_part if inline_part else name_str
    if "+" in target_text or ("-" in target_text and "Spike" not in target_text):
        target_text = target_text.replace("S:", "")
        for t in re.split(r'[+\-]', target_text):
            t = t.strip()
            if t and t not in ["Spike", "del"]: mut_set.add(t)
    return list(mut_set)


def main():
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.insert(0, project_root)
    from utils.paths import MODEL_DIR, PROCESSED_FEATURES_CSV

    os.makedirs(MODEL_DIR, exist_ok=True)

    # Clean potency values and build feature matrix
    df = pd.read_csv(io.StringIO(raw_data))

    # Normalise antibody names
    df['rx_name'] = df['rx_name'].replace({
        'Casirivimab/Imdevimab': 'Casirivimab+Imdevimab',
        'Casirivimab_Imdevimab': 'Casirivimab+Imdevimab',
        'Casirivimab–imdevimab': 'Casirivimab+Imdevimab',
        'Casirivimab-imdevimab': 'Casirivimab+Imdevimab'
    })

    def clean_potency(row):
        val = str(row['potency']).strip()
        if val in ["10000", "50000"] or float(row['potency']) >= float(row['potency_upper_limit']):
            return float(row['potency_upper_limit'])
        return float(val)

    df['cleaned_potency'] = df.apply(clean_potency, axis=1)
    df['mutation_list'] = df['iso_name'].apply(parse_iso_to_mutations)
    all_mutations = sorted(list(set([m for sublist in df['mutation_list'] for m in sublist])))

    rows = []
    for _, row in df.iterrows():
        feature_dict = {mut: 0 for mut in all_mutations}
        for m in row['mutation_list']:
            if m in feature_dict: feature_dict[m] = 1

        rx = row['rx_name']
        wt_rows = df[(df['rx_name'] == rx) & (df['iso_name'].isin(["Wuhan-Hu-1", "Wuhan-Hu-1 Spike", "B.1 Spike"]))]
        wt_potency = wt_rows['cleaned_potency'].mean() if not wt_rows.empty else 1.0

        feature_dict['target_y'] = np.log10(row['cleaned_potency']) - np.log10(wt_potency)
        feature_dict['rx_name'] = rx
        feature_dict['iso_name'] = row['iso_name']
        rows.append(feature_dict)

    ml_df = pd.DataFrame(rows)

    csv_path = PROCESSED_FEATURES_CSV
    ml_df.to_csv(csv_path, index=False)

    antibodies = sorted(ml_df['rx_name'].unique())

    print("Feature matrix export")
    print(f"  Saved: {csv_path}")
    print(f"  Mutation features: {len(all_mutations)}")
    print(f"  Antibodies: {antibodies}\n")

    feature_json_path = os.path.join(MODEL_DIR, "feature_columns.json")
    with open(feature_json_path, "w", encoding="utf-8") as f:
        json.dump(all_mutations, f, indent=2)
    print(f"  Feature columns: {feature_json_path}\n")

    print("Training Ridge models")
    for ab in antibodies:
        ab_df = ml_df[ml_df['rx_name'] == ab]
        X = ab_df[all_mutations].values
        y = ab_df['target_y'].values

        if len(X) == 0:
            print(f"  Skip {ab}: no samples")
            continue

        model = Ridge(alpha=1.0)
        model.fit(X, y)

        safe_ab_name = ab.replace("/", "_").replace("+", "_").replace(" ", "_").replace("–", "_").lower()
        model_path = os.path.join(MODEL_DIR, f"{safe_ab_name}_model.pkl")

        joblib.dump(model, model_path)

        print(f"\n  {ab}: {len(X)} samples -> {model_path}")

        coef_dict = dict(zip(all_mutations, model.coef_))
        sorted_muts = sorted(coef_dict.items(), key=lambda item: abs(item[1]), reverse=True)
        print("  Top weights:")
        for m, weight in sorted_muts[:3]:
            if abs(weight) > 1e-4:
                print(f"    {m}: {weight:+.4f}")

    print("\nDone.")


if __name__ == "__main__":
    main()