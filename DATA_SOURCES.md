# Data sources and attribution

## Reference genomes (NCBI)

Fourteen NCBI reference assemblies under `data/prediction/` are listed in
`manifest.json`. Five accessions (BA.2, BA.4.6/OR325409.1, BQ.1.1, XBB.1.5, JN.1)
are partial genome records with complete spike coding sequence. External validation
genomes under `data/evaluation/external/` were downloaded separately and are not part
of the web reference panel.

When citing sequence accessions, refer to NCBI GenBank records linked in
`manifest.json` and `data/evaluation/pipeline_external_test.json`.

## CoV-DRDB (Stanford)

Expanded IC50 training data were built from the Stanford Coronavirus Antiviral and
Resistance Database payload:

- Repository: https://github.com/hivdb/covid-drdb-payload
- Raw export: `data/training/cov_unibind/expanded_raw_potency.csv`
- Feature matrix: `data/training/processed_ml_features.csv`

To rebuild features locally, clone the payload into
`data/training/local_expansion/drdb_payload/` and run
`scripts/expand_ic50_from_drdb.py`. See `data/training/local_expansion/README.md`.

## CoV-UniBind / Arora baseline

Baseline potency tables and archived source CSVs under `data/training/cov_unibind/`
(e.g. `mutation_ic50.csv`, `wildtype_ic50.csv`, `raw/arora*.csv`) complement DRDB
measurements during feature expansion.

## External validation IC50 labels

Literature comparison values in `data/evaluation/pipeline_external_test.json` are
**recomputed CoV-DRDB medians**: for each spike isoform and antibody, the median of
`cleaned_potency` (ng/ml) in `expanded_raw_potency.csv`. Values at the assay upper
limit (≥10,000 ng/ml) are marked censored.

- Regenerate labels: `python scripts/derive_external_ic50_labels.py`
- Full provenance (study IDs, source files, measurement counts):
  `data/evaluation/pipeline_external_test_provenance.json`
- Wild-type reference for log10 fold-change in external evaluation: B.1 Spike
  medians (Casirivimab 8.35 ng/ml, Imdevimab 7.5 ng/ml from the same table)

## Models

Ridge regression models (`data/training/models/*.pkl`) were trained on
`data/training/ml_train.csv` using scikit-learn with fixed hyperparameters documented in
the manuscript Methods section.
