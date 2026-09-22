# CLDN4 exclusion sensitivity (PPT later slide)

Locks stay labeled. Numbers below are either the published lock or a value written by the scripts from public inputs. The searched minimum is descriptive. It is not a confirmatory p-value.

## Locks

| Role | Estimand | Value |
|---|---|---|
| LOCK | He 2022 CosMx, cytotoxic neighbors, 50 µm | ratio **0.36**, 8/8 sections, 5/5 patients, sign P = 0.031 |
| LOCK | Same, 100 µm | ratio **0.52**, same signs |
| LOCK | Concordant-4 malignant CLDN4 % positive (UMI > 0) vs T/NK fraction | **ρ = −0.531**, N = 65, I² = 0%, p = 1.65×10⁻⁵, Cliff δ = −0.724 |

The CosMx ratios are the published lock (exclusion, not muzzling). This note does not recompute them. The concordant-4 ρ is recomputed from `unit_scores.tsv` and matches the published value.

## CosMx sensitivity

Object: figshare 25976224, `cosmx_human_nsclc_clustered.h5ad`, 765,771 cells, 8 sections, 5 patients. Index cells are the patient-matched tumor label. Centroids are global pixels × 0.18 µm. The index cell is excluded.

Fraction grid: radii 9, 10, 12, 15, 20, 25, 40, 50, 100 µm; rings 10–25, 20–50, and 50–100 µm; classes broad immune, CD8+NK, and CD8; cuts count ≥ 1, 2, 5, 10, 20 and the top 50/25/10/5% of positive counts, each versus count 0; denominator all other cells or non-tumor cells. 408 specs had both arms at n ≥ 30 in every section.

A fraction spec is eligible when the high arm is lower in 8/8 sections and 5/5 patients and every absent-arm section mean is ≥ 0.005. The headline also requires a positive high-arm mean in every section and the same direction among cells that have a neighbor.

Calibration, same code: 10 µm, count ≥ 1 vs 0, broad immune, ratio **0.377**. That matches the previous short-range recomputation.

**Headline fraction.** Broad immune, 9 µm, CLDN4 at or above the section 75th percentile of positive counts versus count 0. Means 0.0112 / 0.0490, ratio **0.229**. Contact fraction among cells with a neighbor: 0.313. Quietest absent section mean 0.0053. 8/8 and 5/5. Rings, the non-tumor denominator, the top 5% cut, and radii through 100 µm did not yield a smaller ratio under the same rule. This is the same specification as the previous 9 µm search.

**Stable count.** Every section low-arm mean ≥ 0.05, every high-arm mean > 0, 8/8 and 5/5. Five specs pass. The smallest ratio is broad immune counts at **15 µm**, positive-count median versus absent: 0.153 / 0.361, ratio **0.423**. Section ratios run from 0.141 (LUAD-9 R2) to 0.900 (LUAD-5 R2). In LUSC-6 the positive-count median is 1, so that section’s cut is detected versus absent. This cell class is broader than the locked cytotoxic count, and 0.423 is a weaker ratio than the locked 0.36.

**Cytotoxic class on this grid.** No cytotoxic fraction spec and no cytotoxic count spec cleared its floor at 8/8 and 5/5. The smallest cytotoxic count ratio with a positive mean in every section is 0.217 at 10 µm (top 25% of positive counts; 0.0028 / 0.0129). The quietest absent section mean is 0.0014, so it stays off the stable call. Detected-versus-absent cytotoxic counts at 50 µm are 0.776 (5/8, 3/5) and at 100 µm are 0.895 (3/8, 2/5). Those rows are a detected-versus-absent count contrast. The locked 0.36 / 0.52 summary remains the cytotoxic lock.

## Concordant-4 sensitivity

Cohorts: GSE123902, GSE131907, GSE205335, GSE189357. Unit is the patient. Grid: seven CLDN4 scores, four T/NK outcomes (fraction, logit, T/NK per malignant cell, T/NK among non-malignant cells), EPCAM and log malignant-cell adjustments, and filters all / primary / drop AIS. A panel is concordant when all four cohorts have a negative Spearman and the stacked Q4 vs Q1 Cliff delta on the raw outcome is negative. Pooled ρ is DerSimonian–Laird on Fisher z.

**Locked outcome.** Malignant CLDN4 % positive versus the T/NK fraction remains the most negative I² = 0 panel on that outcome (ρ = −0.531, Cliff δ = −0.724, N = 65, 95% interval −0.697 to −0.312). The same panel’s median of the four within-cohort Q4/Q1 median ratios is **0.349** (0.471, 0.124, 0.226, 0.517).

**Strongest I² = 0 fold on the T/NK fraction.** Percent of malignant cells with CLDN4 UMI ≥ 2. Median cohort ratio **0.297** (0.242, 0.202, 0.352, 0.539). ρ = −0.455 (p = 3.5×10⁻⁴, I² = 0%, interval −0.641 to −0.218), Cliff δ = −0.671, N = 65. The fold is larger than 0.349 and |ρ| is smaller than 0.531.

**EPCAM partial of the locked score.** Partial Spearman on malignant EPCAM (log1p CP10K): ρ = −0.421, I² = 0%, Cliff δ = −0.678, p = 0.0017, N = 65. The association stays negative and is smaller in magnitude than the lock.

**T/NK per malignant cell is not used.** Unadjusted ρ = −0.536 (I² = 0%, Cliff δ = −0.743). CLDN4 % positive tracks malignant-cell number (Spearman 0.44), and T/NK per malignant cell tracks the reciprocal (Spearman −0.64). After a partial on log malignant-cell number the ρ is −0.198 (I² = 0%). Dropping AIS gives ρ = −0.600 with I² = 32%. Primary-only could not be scored: GSE131907 has fewer than 6 primary units in this table.

## Reproduce

```bash
python3 scripts/concordant4_outcome_sensitivity.py
python3 scripts/cosmx_cytotoxic_fraction_sensitivity.py
python3 -m pytest tests/test_concordant4_sensitivity.py -q
```

The CosMx script expects `data/cosmx_nsclc/cosmx_human_nsclc_clustered.h5ad` (figshare file 46841842). That file is gitignored. Dependencies: `scripts/requirements-cldn4-sensitivity.txt`.
