# B1_CHOL scripts

Public-data analog of the B1 TACSTD2–CLDN4 surface-rank / co-expression
readout, run in TCGA-CHOL.

```bash
pip install pandas numpy scipy matplotlib openpyxl
python3 scripts/w200/B1_CHOL/download_data.py   # → /tmp/b1_chol_data (override: B1_CHOL_DATA)
python3 scripts/w200/B1_CHOL/run_analysis.py    # → results/w200/B1_CHOL/
```

Raw matrices are not committed. Every number under `results/w200/B1_CHOL/`
is produced by `run_analysis.py` from those downloads.
