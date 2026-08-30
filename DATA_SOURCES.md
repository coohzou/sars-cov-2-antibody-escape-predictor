# Data sources and attribution

## Reference genomes (NCBI)

Fourteen complete SARS-CoV-2 genomes under `data/prediction/` are listed in
`manifest.json`. External validation genomes under `data/evaluation/external/` were
downloaded separately and are not part of the web reference panel.

When citing sequence accessions, refer to NCBI GenBank records linked in
`manifest.json` and `data/evaluation/pipeline_external_test.json`.

## CoV-DRDB (Stanford)

Expanded IC50 training data were built from the Stanford Coronavirus Antiviral and
Resistance Database payload:

- Repository: https://github.com/hivdb/covid-drdb-payload
- Processed tables: `data/training/cov_unibind/expanded_raw_potency.csv`
- Feature matrix: `data/training/processed_ml_features.csv`

To rebuild features locally, clone the payload into `data/training/local_expansion/drdb_payload/` and
run `scripts/expand_ic50_from_drdb.py`. See `data/training/local_expansion/README.md`.

## CoV-UniBind / Arora baseline

Baseline potency tables and archived source CSVs live under `data/training/cov_unibind/`.
These complement DRDB measurements for model training.

## Evaluation labels

External validation IC50 literature values in `data/evaluation/pipeline_external_test.json`
were taken from CoV-DRDB medians for the corresponding spike isoforms.

## Models

Ridge regression models (`data/training/models/*.pkl`) were trained on
`data/training/ml_train.csv` using scikit-learn with fixed hyperparameters documented in
the manuscript Methods section.
