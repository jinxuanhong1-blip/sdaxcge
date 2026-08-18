# GSE154989 mouse-level Cldn4 scores

Public KP/K plate scRNA (Marjanovic et al., processed GEO h5). Not KL.

- `mouse_level_scores.tsv` — one row per biological mouse (tumor IDs collapsed). `pass_min_cells` marks mice with ≥10 cells (n=28 scored).
- `high_vs_low_summary.tsv` — Spearman and median-split tests for Cldn4 vs T/NK and IFN/MHC.
- `genes_used.tsv` — Cldn4-only plus which T/NK and IFN/MHC symbols were present.
- `cell_level_scores.tsv` — per-cell values used for the mouse means.

See `methods/mouse_kl_cldn4_hunt/FINDING.md`.
