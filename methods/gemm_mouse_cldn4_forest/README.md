# Mouse-level Cldn4 %pos vs immune / IFN

Public autochthonous GEMM lung scRNA only. The private 8 KL matrices are not used.

Primary question, fixed before the pooled p-value: among biological mice, is the percent of epithelial cells with Cldn4 detected associated with the T/NK fraction (unsorted digests) or with an IFN-only score inside those epithelial cells?

The pool is a random-effects meta-analysis of within-study Spearman correlations (Fisher z, REML, Hartung-Knapp). A pooled p-value is reported only with at least 3 studies of at least 4 mice. The sensitivity grid is written out in full. The primary rows are not replaced by the spec with the smallest p-value.

## Run

```bash
# optional: rescore GEO matrices into tables/new_library_scores.tsv
python3 methods/gemm_mouse_cldn4_forest/score_matrices.py
python3 methods/gemm_mouse_cldn4_forest/score_matrices.py --remaining

python3 methods/gemm_mouse_cldn4_forest/analyze.py
```

`analyze.py` reads the vendored mouse tables in `data/published/` plus `tables/new_library_scores.tsv`. Raw count matrices are downloaded from GEO and are not stored in this repository.

## Outputs

- `FINDING.md` — primary result and the sensitivity note
- `tables/mouse_level.tsv` — one row per mouse
- `tables/primary_study_effects.tsv` — within-study Spearman
- `tables/sensitivity_grid.tsv` — every pre-listed spec
- `figures/forest_primary.png` — T/NK and IFN forests
- `figures/sensitivity_pooled.png` — pooled estimate under each spec
