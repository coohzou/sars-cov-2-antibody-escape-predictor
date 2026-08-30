# SARS-CoV-2 Antibody Escape Predictor

Flask web application for SARS-CoV-2 variant identification, spike mutation detection, and Casirivimab/Imdevimab cocktail neutralization (IC50) prediction.

## Quick start

```bash
git clone https://github.com/coohzou/sars-cov-2-antibody-escape-predictor.git
cd sars-cov-2-antibody-escape-predictor
python -m venv .venv
# Windows
.\.venv\Scripts\pip install -r requirements.txt
.\.venv\Scripts\python.exe app.py
# macOS / Linux
# source .venv/bin/activate && pip install -r requirements.txt && python app.py
```

Open http://127.0.0.1:5000 and upload a `.fasta` file.

## Live demo (Render)

Production deployment uses [Render Standard (`1c-2g`)](https://render.com) for always-on hosting with sufficient RAM for genome alignment.
After connecting the GitHub repo in the Render dashboard, the service URL is typically:

`https://sars-cov-2-antibody-escape-predictor.onrender.com`

Step-by-step setup: [`docs/DEPLOY_RENDER.md`](docs/DEPLOY_RENDER.md).

Pre-trained Ridge models are included under `data/training/models/`. If `/ready` reports
`predictor_ready=false`, run:

```bash
python scripts/train_neutralization_model.py
```

## Project structure

```
sars-cov-2-antibody-escape-predictor/
├── app.py
├── requirements.txt
├── templates/index.html
├── utils/
├── data/
│   ├── training/       # ML features, splits, models
│   ├── prediction/     # 14 NCBI reference genomes
│   └── evaluation/     # Held-out test configs and results
├── scripts/
└── tests/
```

See `data/README.md` for the data layout and `REPRODUCIBILITY.md` for evaluation commands.

## Reproduce paper metrics

```bash
python scripts/evaluate_ml_test_set.py
python scripts/evaluate_pipeline_test_set.py
python scripts/evaluate_pipeline_external_test.py
```

Or run the full workflow (requires a local CoV-DRDB clone for step 1):

```bash
python scripts/run_evaluation.py
```

| Test set | Samples | Metric | Result |
|----------|---------|--------|--------|
| Pipeline test | 7 NCBI strains | Major-category ID + mutation recall | 7/7, 100% |
| External validation | 7 NCBI genomes | Major-category ID | 7/7, 100% |
| ML test | 40 rows (20 held-out isoforms) | log10(fold) MAE / r | MAE 0.50, r 0.87 |

Outputs: `data/evaluation/ml_test_results.json`, `pipeline_test_results.json`,
`pipeline_external_test_results.json`.

## Recommended test files

| Variant | Path |
|---------|------|
| Wild Type | `data/prediction/wuhan_hu1_complete.fasta` |
| Gamma | `data/prediction/gamma_complete.fasta` |
| Delta | `data/prediction/delta_complete.fasta` |
| JN.1 | `data/prediction/jn1_complete.fasta` |

## API

| Route | Method | Description |
|-------|--------|-------------|
| `/` | GET | Main page |
| `/upload` | POST | Upload FASTA, returns JSON analysis |
| `/health` | GET | Liveness check |
| `/ready` | GET | Model and reference readiness |

Run tests:

```bash
python -m unittest discover -s tests -p "test_*.py"
```

## Dependencies

Python 3.12 recommended (see `Dockerfile` and `render.yaml`). Core packages: Flask, Biopython,
pandas, scikit-learn, joblib, numpy.

## License

MIT. See `LICENSE`. Training data sources are documented in `DATA_SOURCES.md`.
