# scripts

| Script | What it does |
| --- | --- |
| `lib.py` | Histology map, Spearman, Fisher-z meta |
| `extract_gse241934.py` | Stream TACSTD2/CLDN4 from GEO MTX; author Epi vs T+NK |
| `01_compute_cohort_effects.py` | Histology-stratified ρ from harvested tables |
| `02_meta_and_forest.py` | Within-histology meta + forest |
| `selftest.py` | Recompute a few locked numbers |

GSE241934 MTX stays in `/tmp/gse241934/` (not committed). The per-sample extract is in `harvested/`.
