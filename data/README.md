# Data directory guide

Project data are organised into three folders by purpose.

## 1. `data/training/` — machine learning

Used to train Casirivimab / Imdevimab IC50 Ridge models (offline only).

| Path | Description |
|------|-------------|
| `cov_unibind/` | Processed IC50 tables (`mutation_ic50.csv`, `wildtype_ic50.csv`) |
| `cov_unibind/raw/` | Archived Arora source CSVs |
| `processed_ml_features.csv` | DRDB-expanded feature matrix (560 rows, 95 features) |
| `cov_unibind/expanded_raw_potency.csv` | Raw Casi/Imdev IC50 rows from CoV-DRDB + Arora |
| `ml_train.csv` | Training split (507 rows) |
| `ml_test.csv` | Held-out ML test split (53 rows; 20 isoforms) |
| `models/` | Trained Ridge models and `feature_columns.json` |
| `drdb_payload/` | Optional local clone of CoV-DRDB (not committed) |

Scripts: `expand_ic50_from_drdb.py`, `split_ml_dataset.py`, `train_neutralization_model.py`.

## 2. `data/prediction/` — web runtime

Fourteen NCBI reference genomes loaded by the web app for variant matching and spike
mutation calling.

| Path | Description |
|------|-------------|
| `manifest.json` | Variant metadata (name, lineage, accession, FASTA filename) |
| `*_complete.fasta` | Full genomes including Wuhan-Hu-1 |

Used by `utils/sequence_comparator.py` and `app.py`. Refresh with `scripts/download_variants.py`.

## 3. `data/evaluation/` — test configuration and results

Held-out evaluation outputs (separate from training data).

| Path | Description |
|------|-------------|
| `split_config.json` | ML / pipeline test split rules |
| `split_summary.json` | Split statistics |
| `ml_test_results.json` | ML test metrics (MAE, RMSE, r) |
| `pipeline_test_results.json` | Pipeline test results (7 strains) |
| `pipeline_external_test.json` | External validation case definitions |
| `external/*.fasta` | Seven external NCBI genomes |

Run: `python scripts/run_evaluation.py` (see `REPRODUCIBILITY.md`).

## Quick reference

| Task | Folder |
|------|--------|
| Train / retrain models | `data/training/` |
| Upload analysis (web app) | `data/prediction/` |
| Report test metrics | `data/evaluation/` |
