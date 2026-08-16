# A10 results: GRHL1 vs TACSTD2 / CLDN4 (public TCGA lung)

**Verdict on the pre-specified primary cohort (mixed NSCLC primary tumors, n=1017): `PARTIAL_ONLY_ONE_TARGET`.**

GRHL1 tracks TACSTD2. It does not track CLDN4 once LUAD and LUSC are pooled.

## Claim pairs (Spearman)

| slice | n | GRHL1–TACSTD2 | GRHL1–CLDN4 | both ≥ 0.30 and p ≤ 0.05 |
|---|---:|---|---|---|
| NSCLC primary (primary) | 1017 | ρ = 0.445, p = 1.6e-50 | ρ = 0.021, p = 0.50 | no |
| LUAD primary (sensitivity) | 515 | ρ = 0.397, p = 7.3e-21 | ρ = 0.391, p = 3.2e-20 | yes |
| LUSC primary (sensitivity) | 502 | ρ = 0.387, p = 2.4e-19 | ρ = 0.245, p = 2.8e-8 | no (CLDN4 is only weak) |

The mixed-NSCLC CLDN4 number is not “no biology.” Pooling LUAD and LUSC can cancel a within-LUAD GRHL1–CLDN4 correlation. That is reported as a sensitivity, not used to change the primary score.

## What this is not

- Not evidence that GRHL1 directly transcribes TACSTD2 or CLDN4.
- Not a DepMap / cell-line result. This is bulk TCGA tumor RNA only.
- GRHL2 / GRHL3 rows in `correlations.csv` are context. They do not rescue the GRHL1 claim.

## Files

- `expression_tcga_lung.csv` — tidy extract (GRHL1/2/3, TACSTD2, CLDN3/4/7, ELF3, KLF5)
- `correlations.csv` — all pairs × slices
- `summary.json` — machine-readable verdict
- `fig_grhl1_vs_targets_nsclc.png` — scatter, LUAD vs LUSC colored
- `fig_grhl1_rho_by_histology.png` — ρ bars vs the 0.30 support line
- `download_manifest.json` — Xena URLs and SHA-256 of the source matrices
