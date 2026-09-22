# Fig. 1 — CosMx cold neighborhoods

Section-level CosMx NSCLC (He et al. 2022) panels. Python and matplotlib only. No methods schematic.

```bash
python3 nature_figures/fig1_cosmx_cold/plot_fig1_cosmx_cold.py
```

Writes `fig1_cosmx_cold.pdf`, `.svg`, and `.png`, plus `source_stats.json`.

## What is drawn

| Panel | Contrast | Source table |
|---|---|---|
| a, c (TACSTD2) | Section-median TACSTD2, CD8+NK neighbor counts at 10/20/50/100 µm | `source/marginal_cd8nk_by_section.csv` |
| b | CLDN4 Q4 vs Q1, CD8+NK counts at 20/40 µm (60/80 in the same file) | `source/cldn4_q4q1_cd8nk_by_section.csv` |
| c (CLDN4) | Section-median CLDN4, same count definition as panel a | `source/marginal_cd8nk_by_section.csv` |
| d | GZMB, PRF1, NKG7, IFNG CPM in nearby CD8 or NK cells | `source/cldn4_nearby_effector_cpm.csv` |

Provenance of the copied rows:

- Marginal TACSTD2/CLDN4 counts: commit `b79c19a`, `results/cosmx_tacstd2_cldn4/tables/marginal_by_section.csv` (median cut, `cd8nk_n`).
- Quartile counts: commit `9e6fedf`, `results/cosmx_nsclc_cldn4_nk_cd8/tables/section_paired_cd8nk_counts.csv`.
- Effector CPM: commit `db7689d`, `results/cosmx_squidpy_fov/tables/muzzling_ratios.csv`.

`source/cldn4_median_squidpy_cd8nk_by_section.csv` is the Squidpy median-split count table at 50/100 µm (same commit). It is not drawn. Its ratios agree with the median-split fade in panel c (about 0.88 and 0.99 as median section ratios; not 8/8).

## What is not drawn

A handoff summary of cytotoxic neighbor ratios 0.36 at 50 µm and 0.52 at 100 µm (8/8, sign P = 0.031) is cited in later notes and has no section-level rows in this repository. Those two numbers are not plotted. Full cell matrices and spatial coordinates are not in the checkout and were not downloaded.
