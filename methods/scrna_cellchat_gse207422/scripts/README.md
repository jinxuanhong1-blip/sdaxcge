# scripts

| Script | Role |
|---|---|
| `00_download.py` | GEO processed UMI + metadata; CellPhoneDB v5 CSVs |
| `01_build_pairs.py` | Parse CellPhoneDB `interactors` → `resources/cellphonedb_v5_lr_pairs.tsv` |
| `02_run_ccc.py` | A3 annotation, documented score, optional LIANA, figures, `RESULTS.md` |
