# TROP2–CLDN4 protein: histology and partial correlations

Recompute TACSTD2 (TROP2) versus CLDN4 on the public Gygi/Nusinow CCLE protein table. The cohort list is fixed in `analyze.py`. The reported maximum is the largest Spearman in that list with n ≥ 15 (partials also need residual df ≥ 10). No cell line is dropped to raise ρ, and no table is built to equal n = 118.

```bash
python3 scripts/depmap_trop2_cldn4_maxrho/download.py
python3 scripts/depmap_trop2_cldn4_maxrho/analyze.py
```

Raw downloads stay in `data/` (gitignored). Tables and figures land in `results/depmap_trop2_cldn4_maxrho/`.
