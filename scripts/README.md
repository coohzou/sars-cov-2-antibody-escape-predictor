# Scripts

| Script | Purpose |
|--------|---------|
| `run_evaluation.py` | Full workflow: expand → split → train → three evaluations |
| `split_ml_dataset.py` | Split `processed_ml_features.csv` into train/test |
| `train_neutralization_model.py` | Train Ridge models on `ml_train.csv` only |
| `evaluate_ml_test_set.py` | IC50 metrics on held-out `ml_test.csv` |
| `evaluate_pipeline_test_set.py` | Variant ID + mutation recall on 7 held-out NCBI strains |
| `evaluate_pipeline_external_test.py` | External NCBI genomes + DRDB IC50 comparison |
| `expand_ic50_from_drdb.py` | Rebuild `processed_ml_features.csv` from CoV-DRDB |
| `build_ml_features.py` | Legacy Arora-only rebuild (~83 rows) |
| `download_variants.py` | Re-download reference FASTA from NCBI |
| `derive_external_ic50_labels.py` | Recompute external-validation IC50 medians from CoV-DRDB export |

## Typical workflows

```bash
# Verify paper metrics (no DRDB clone)
python scripts/evaluate_ml_test_set.py
python scripts/evaluate_pipeline_test_set.py
python scripts/evaluate_pipeline_external_test.py

# Full rebuild (requires data/training/local_expansion/drdb_payload/)
python scripts/run_evaluation.py

# Start web app
python app.py
```
