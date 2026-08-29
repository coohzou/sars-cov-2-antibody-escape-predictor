# Evaluation data

## Train / test splits

- **ML test** (`data/training/ml_test.csv`): 53 rows, 20 held-out isoforms; excluded from training.
- **Pipeline test** (7 variants in `split_config.json`): Alpha, Beta, Gamma, Delta, Omicron, BA.4, JN.1.
- **External pipeline test** (`pipeline_external_test.json`): 7 NCBI complete genomes not in
  `manifest.json`, with CoV-DRDB IC50 labels for prediction comparison.

See `split_config.json` for rules.

## Output files

```
data/evaluation/split_config.json
data/evaluation/split_summary.json
data/evaluation/ml_test_results.json
data/evaluation/pipeline_test_results.json
data/evaluation/pipeline_test_summary.csv
data/evaluation/pipeline_external_test_results.json
data/evaluation/pipeline_external_test_summary.csv
data/evaluation/external/*.fasta
```

## Run

```bash
python scripts/run_evaluation.py
```

For a fast verification path without rebuilding DRDB features, see `REPRODUCIBILITY.md`.
