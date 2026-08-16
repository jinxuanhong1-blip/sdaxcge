# B2 lung-only TACSTD2–CLDN4 ρ (DepMap / CCLE)

Recomputes Spearman and Pearson correlation between **TACSTD2** and **CLDN4** on **DepMap Public 24Q4** lung cell lines only.

The user-claimed value is **ρ = 0.69**. This pipeline does not tune lineage, histology, or correlation method to hit that number.

```bash
python3 scripts/w200/B2_lung/download.py
python3 scripts/w200/B2_lung/analyze.py
```

Outputs land in `results/w200/B2_lung/`.
