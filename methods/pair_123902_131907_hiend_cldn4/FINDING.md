# FINDING — CLDN4-only high-end CellChat / LIANA on GSE123902 + GSE131907

ADDITIVE. **CLDN4 only. No dual-high.** Patient/donor is the unit
(GSE123902 donor; GSE131907 sample). Do **not** add GSE148071. Do
**not** pile the seven-cohort pool.

The pair %pos Spearman is **taken as given** from PR #459 and is **not
re-audited**:

| combo | score | N | ρ (p, I²) | Q4 vs Q1 r (p) |
|---|---|---:|---|---|
| GSE123902+GSE131907 | %pos | 34 | −0.575 (0.001, 0%) | −0.700 (0.012; 10 vs 8) |

Singles (given, same PR): GSE123902 n=13 ρ=−0.659; GSE131907 n=21 ρ=−0.522.

LR tables are written by `scripts/analyze.py` after the public GEO
matrices are scored. Primary files:

- `results/lr_table.tsv` — CellChat-style outgoing (n_patients ≥ 8)
- `results/lr_table_liana.tsv` — LIANA-style outgoing (n_patients ≥ 8)

This stub is replaced when the analyzer finishes.
