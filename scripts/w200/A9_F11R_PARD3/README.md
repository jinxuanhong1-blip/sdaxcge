# A9 — F11R and PARD3 with TACSTD2-high public lung

Honest public-data test of whether **F11R** (JAM-A) and **PARD3** track
TACSTD2-high lung tumors. Outputs land in `results/w200/A9_F11R_PARD3/`.

This does **not** tune lineage, split, or correlation method to force a
positive call. Pooled LUAD+LUSC is a sensitivity only: TACSTD2-high tumors
are LUSC-enriched.

```bash
python3 scripts/w200/A9_F11R_PARD3/download.py
python3 scripts/w200/A9_F11R_PARD3/analyze.py
# optional stroma-free slice
python3 scripts/w200/A9_F11R_PARD3/download.py --with-depmap
python3 scripts/w200/A9_F11R_PARD3/analyze.py
```

## Pre-specified call

A claim gene is `up_in_tacstd2_high` only if, **within one histology**:

1. Spearman ρ vs continuous TACSTD2 is > 0
2. Tertile high-minus-low log2FC is > 0
3. Welch t-test BH-FDR (across F11R and PARD3 in that cohort) is < 0.05

Primary cohorts are TCGA-LUAD and TCGA-LUSC primary tumors, separately.
OncoSG LUAD 2020 is the independent LUAD check. DepMap 24Q4 lung cell lines
are a no-stroma sensitivity.

## Sources

All public / open access. Raw Xena matrices are not committed (see
`results/w200/A9_F11R_PARD3/data/.gitignore`).
