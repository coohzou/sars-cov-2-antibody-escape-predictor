# Reproducibility guide

This repository supports two workflows: **verify cached metrics** (fast) and **full rebuild**
(requires downloading CoV-DRDB).

## Environment

```bash
python -m venv .venv
# Windows
.\.venv\Scripts\pip install -r requirements.txt
# macOS / Linux
# source .venv/bin/activate && pip install -r requirements.txt
```

Python 3.12 is used in Docker and Render deployments.

## Workflow A — verify paper metrics (recommended)

Uses committed CSVs, models, and evaluation configs. No external DRDB clone required.

```bash
python scripts/split_ml_dataset.py
python scripts/train_neutralization_model.py
python scripts/evaluate_ml_test_set.py
python scripts/evaluate_pipeline_test_set.py
python scripts/evaluate_pipeline_external_test.py
```

Expected headline results:

- Pipeline test: 7/7 major-category identification, 100% key-mutation recall
- External validation: 7/7 identification; MAE on log10 fold-change ≈ 0.45 (10 uncensored)
- ML test: MAE ≈ 0.50, Pearson r ≈ 0.87 on 40 antibody measurements (20 isoforms)

Compare outputs under `data/evaluation/` with the values reported in the manuscript.

## Workflow B — full rebuild from CoV-DRDB

```bash
git clone --depth 1 https://github.com/hivdb/covid-drdb-payload.git data/training/drdb_payload
python scripts/run_evaluation.py
```

`run_evaluation.py` runs, in order:

1. `expand_ic50_from_drdb.py`
2. `split_ml_dataset.py`
3. `train_neutralization_model.py`
4. `evaluate_ml_test_set.py`
5. `evaluate_pipeline_test_set.py`
6. `evaluate_pipeline_external_test.py`

Step 1 requires the DRDB payload directory above. The repository ships pre-built
`processed_ml_features.csv`, train/test splits, and trained models so Workflow A works
without cloning DRDB.

## Web application smoke test

```bash
python app.py
```

Then upload `data/prediction/gamma_complete.fasta` via the UI or POST to `/upload`.
Check `/ready` returns `predictor_ready: true`.

## Unit tests

```bash
python -m unittest discover -s tests -p "test_*.py"
```

## Refresh reference genomes

```bash
python scripts/download_variants.py
```

Requires network access to NCBI E-utilities.
