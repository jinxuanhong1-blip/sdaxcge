# A11 galectin analysis

Run from the repository root:

```bash
python3 scripts/w200/A11_galectin/analyze.py
```

The script retrieves public TCGA LUAD and LUSC PanCancer Atlas primary-tumor
RNA-seq values through the cBioPortal API. It writes the concise interpretation,
gene-level statistics, sample-level expression extract, and provenance metadata
to `results/w200/A11_galectin/`.

Runtime dependencies: Python 3, NumPy, SciPy, and Requests.
