# A11 galectin analysis

Run from the repository root:

```bash
python3 scripts/w200/A11_galectin/analyze.py
```

The script retrieves public lung RNA through the cBioPortal API:

- Primary bulk: TCGA LUAD and LUSC PanCancer Atlas
- Independent bulk: OncoSG LUAD, CPTAC LUAD, CPTAC LUSC, CAS LUAD
- Tumor-cell check: CCLE / DepMap Broad 2025 lung and NSCLC cell lines
- Purity check: CPTAC ESTIMATE and CAS pathologist purity partial Spearman

Outputs are written to `results/w200/A11_galectin/`.

Runtime dependencies: Python 3, NumPy, SciPy, and Requests.
