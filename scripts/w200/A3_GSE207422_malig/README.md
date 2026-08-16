# W200-A3 · GSE207422 malignant TACSTD2

Recompute Hu et al. 2023 (GEO GSE207422) **malignant-only TACSTD2**: NMPR vs MPR, and vs T/NK fraction.

Processed UMI matrix is 175.5 MB gzip (< 2 GB budget) → **recompute**, do not catalog-only.

Author per-cell labels are **not public** (see `results/w200/A3_GSE207422_malig/label_search.json`).
Compartments are marker-based, aligned to the paper's published canonical markers.

```bash
python3 scripts/w200/A3_GSE207422_malig/download.py
python3 scripts/w200/A3_GSE207422_malig/extract.py
python3 scripts/w200/A3_GSE207422_malig/analyze.py
```
