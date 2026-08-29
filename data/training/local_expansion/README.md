# Local IC50 expansion workspace

The main pipeline now uses **`scripts/expand_ic50_from_drdb.py`**, which reads DRDB data from
`drdb_payload/` in this folder.

To refresh DRDB sources:

```powershell
git -C data/training/local_expansion/drdb_payload pull
python scripts/expand_ic50_from_drdb.py
python scripts/run_evaluation.py
```

Legacy exploratory outputs (`models_panel/`, `expansion_summary.json`, etc.) may remain here from
the first local experiment; the canonical training artifacts live under `data/training/`.
