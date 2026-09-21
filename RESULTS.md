# TISMO ICB: Cldn4 versus NHEJ and IFN gene scores

Public recompute on the PR #542 pairing. Tacstd2 **49/64** is unchanged. Every number below is in `tables/`.

## Paper sentence

On the locked 64 TISMO ICB slices (Tacstd2 still 49/64 up, Wilcoxon p = 5.84×10⁻⁵; Cldn4 still 34/64 up, 28 down, 2 ties, Wilcoxon p = 0.0767), the ICB change in Cldn4 is uncorrelated with the ICB change in Hallmark IFN-γ (ρ = -0.03 (95% bootstrap CI -0.29 to +0.25, p = 0.828, q = 0.828, n = 64)) and with the GO:0006303 NHEJ change (ρ = +0.05 (95% bootstrap CI -0.19 to +0.30, p = 0.68, q = 0.828, n = 64)). Averaging slices to 17 cell lines leaves those change-correlations null (IFN-γ ρ = -0.16, p = 0.541; NHEJ ρ = -0.10, p = 0.708). IFN-γ rose in 49/64 slices (Wilcoxon p = 1.20×10⁻⁶) and NHEJ fell in 46/64 (Wilcoxon p = 6.36×10⁻⁵); Cldn4 does not follow either shift. A slice-level baseline correlation of Cldn4 with IFN-γ (ρ = +0.39 (95% bootstrap CI +0.13 to +0.60, p = 0.00128, q = 0.00768, n = 64)) does not remain after cell-line averaging (ρ = +0.08, p = 0.772, n = 17). Baseline Cldn4 versus NHEJ is ρ = +0.11 (95% bootstrap CI -0.15 to +0.36, p = 0.375, q = 0.751, n = 64) at slice level and ρ = +0.17 (p = 0.516) across cell lines. After ICB, slice-level Cldn4 versus IFN-γ is ρ = +0.17 (95% bootstrap CI -0.10 to +0.41, p = 0.184, q = 0.551, n = 64) and versus NHEJ is ρ = +0.04 (95% bootstrap CI -0.24 to +0.30, p = 0.754, q = 0.828, n = 64); cell-line means are IFN-γ ρ = -0.26 (p = 0.305) and NHEJ ρ = +0.11 (p = 0.68). Lung ICB in this pairing is LLC only (2/64), not KL or KP. Scores are bulk-tumor RNA.

## What was held fixed

Live TISMO Gene-module export matched the archived 64-slice means (Tacstd2 max |Δ difference| = 9.65e-16; Cldn4 max |Δ difference| = 1.72e-15). Naive is `Baseline == 1`. ICB is `Baseline == 0` (responders and non-responders pooled; anti-PD-1, anti-PD-L1, anti-PD-L2, anti-CTLA4, and combos). The paired location is the mean, the same rule that produces Tacstd2 49/64. Scores use only samples present in the Tacstd2 export. The current TISMO export keeps each sample once. Two slices had ICB rows and no isotype rows; their baselines are the same sample IDs the PR #542 export repeated onto that arm (4T1_GSE130472_old_antiCTLA4 from 4T1_GSE130472_old_antiPDL1 (n=7); 4T1_GSE130472_young_antiCTLA4 from 4T1_GSE130472_young_antiPDL1 (n=7)).

TISMO has no KL or KP lung ICB model. The two lung slices are LLC (GSE155972), which is Lewis lung carcinoma, not KL.

## Scores

Each score is the mean of per-gene z-scores. Z-scores are fit once on the samples that belong to the 64 stems (population sd, ddof = 0). A gene enters a score only if it is non-missing in at least 80% of those samples and has non-zero variance. A sample is scored only if at least 80% of the retained genes are present.

Primary IFN is MSigDB mouse Hallmark interferon-gamma response (185/188 genes used; absent in TISMO: Cxcl11, Mx2, Tmt1b). Primary NHEJ is MSigDB mouse GO:0006303, double-strand break repair via nonhomologous end joining (78/79 genes used; absent in TISMO: Bend2). TISMO still uses the previous MGI symbols Wars and Ddx58; those columns fill Hallmark members Wars1 and Rigi. Reactome NHEJ (R-MMU-5693571) is secondary because 33 of its 67 genes are histones; a histone-stripped Reactome score and a 12-gene c-NHEJ machinery score are reported beside it. Hallmark interferon-alpha and the 7-gene WikiPathways NHEJ set are secondary.

These are bulk syngeneic tumors. The IFN-γ score mixes immune infiltrate with any tumor-cell interferon response. It is not a malignant-cell program.

## Primary Spearman tests (BH q across these six)

| Contrast | n | ρ | 95% CI | p | q |
|---|---:|---:|---|---:|---:|
| delta Cldn4 vs IFN_gamma | 64 | -0.03 | -0.29 to +0.25 | 0.828 | 0.828 |
| delta Cldn4 vs NHEJ_GO | 64 | +0.05 | -0.19 to +0.30 | 0.68 | 0.828 |
| baseline Cldn4 vs IFN_gamma | 64 | +0.39 | +0.13 to +0.60 | 0.00128 | 0.00768 |
| baseline Cldn4 vs NHEJ_GO | 64 | +0.11 | -0.15 to +0.36 | 0.375 | 0.751 |
| icb Cldn4 vs IFN_gamma | 64 | +0.17 | -0.10 to +0.41 | 0.184 | 0.551 |
| icb Cldn4 vs NHEJ_GO | 64 | +0.04 | -0.24 to +0.30 | 0.754 | 0.828 |

Partial Spearman (Pearson correlation of rank residuals), not in the six-test family:

| Contrast | n | ρ | p |
|---|---:|---:|---:|
| delta Cldn4 vs IFN_gamma given NHEJ_GO | 64 | -0.02 | 0.863 |
| delta Cldn4 vs NHEJ_GO given IFN_gamma | 64 | +0.05 | 0.696 |
| baseline Cldn4 vs IFN_gamma given NHEJ_GO | 64 | +0.38 | 0.00188 |
| baseline Cldn4 vs NHEJ_GO given IFN_gamma | 64 | -0.03 | 0.787 |

Holding NHEJ constant, the partial Spearman of the Cldn4 ICB change with IFN-γ is -0.02 (p = 0.863). Holding IFN-γ constant, the partial Spearman with NHEJ is +0.05 (p = 0.696).

## Do the scores themselves move after ICB?

| Score | Up | Down | Tie | Mean Δ | Wilcoxon p | Binomial p |
|---|---:|---:|---:|---:|---:|---:|
| IFN_gamma | 49 | 15 | 0 | +0.159 | 1.20×10⁻⁶ | 2.44×10⁻⁵ |
| IFN_alpha | 48 | 16 | 0 | +0.160 | 8.97×10⁻⁶ | 7.73×10⁻⁵ |
| NHEJ_GO | 18 | 46 | 0 | -0.109 | 6.36×10⁻⁵ | 0.000617 |
| NHEJ_Reactome | 16 | 46 | 0 | -0.121 | 1.99×10⁻⁵ | 0.000176 |
| NHEJ_Reactome_no_histone | 17 | 47 | 0 | -0.146 | 1.05×10⁻⁵ | 0.000227 |
| NHEJ_WP | 17 | 47 | 0 | -0.215 | 9.55×10⁻⁶ | 0.000227 |
| NHEJ_core | 20 | 44 | 0 | -0.110 | 0.000124 | 0.00369 |

Cldn4’s own paired test is unchanged from PR #542 (34 up, 28 down, 2 ties). IFN-γ is also up in 49/64 slices. The overlap with the Tacstd2-up slices is 38/49 (about 37.5 expected if the two up-sets were independent), so the matching counts are not one shared set of slices and they do not revise the Tacstd2 lock.

## Sign concordance with Cldn4

Slices with a zero Cldn4 delta or a zero score delta are omitted.

| Score | Same direction | Cldn4 up, score down | Cldn4 down, score up | n |
|---|---:|---:|---:|---:|
| IFN_gamma | 35 | 7 | 20 | 62 |
| NHEJ_GO | 26 | 26 | 10 | 62 |
| IFN_alpha | 34 | 8 | 20 | 62 |
| NHEJ_core | 31 | 23 | 8 | 62 |
| NHEJ_Reactome | 24 | 26 | 10 | 62 |
| NHEJ_Reactome_no_histone | 26 | 27 | 9 | 62 |
| NHEJ_WP | 28 | 26 | 8 | 62 |

## Dependence checks

The 64 slices are not 64 independent tumors. They are 17 cell lines and 22 studies, and one study contributes many slices. The cell-line check averages slice deltas within each cell line and correlates those means. The leave-one-study check drops the study with the most slices.

| Check | n | IFN-γ ρ | IFN-γ p | NHEJ ρ | NHEJ p |
|---|---:|---:|---:|---:|---:|
| 64 slices (primary unit) | 64 | -0.03 | 0.828 | +0.05 | 0.68 |
| cell-line mean of slice deltas | 17 | -0.16 | 0.541 | -0.10 | 0.708 |
| drop GSE124821 (19 slices) | 45 | +0.04 | 0.817 | -0.07 | 0.628 |

The same slice-versus-cell-line check for Cldn4 levels, rather than changes:

| When | Unit | n | IFN-γ ρ | IFN-γ p | NHEJ ρ | NHEJ p |
|---|---|---:|---:|---:|---:|---:|
| baseline | 64 slices | 64 | +0.39 | 0.00128 | +0.11 | 0.375 |
| baseline | cell-line means | 17 | +0.08 | 0.772 | +0.17 | 0.516 |
| after ICB | 64 slices | 64 | +0.17 | 0.184 | +0.04 | 0.754 |
| after ICB | cell-line means | 17 | -0.26 | 0.305 | +0.11 | 0.68 |

Cancer groups with at least 8 slices, correlation of the ICB changes:

| Group | n | IFN-γ ρ | IFN-γ p | NHEJ ρ | NHEJ p |
|---|---:|---:|---:|---:|---:|
| Colorectal | 10 | +0.27 | 0.446 | -0.03 | 0.934 |
| Mammary | 29 | -0.05 | 0.784 | +0.15 | 0.452 |
| Melanoma | 14 | +0.15 | 0.62 | -0.46 | 0.0962 |

Baseline levels inside those same groups:

| Group | n | IFN-γ ρ | IFN-γ p | NHEJ ρ | NHEJ p |
|---|---:|---:|---:|---:|---:|
| Colorectal | 10 | +0.59 | 0.0739 | +0.48 | 0.162 |
| Mammary | 29 | +0.11 | 0.565 | -0.09 | 0.647 |
| Melanoma | 14 | +0.66 | 0.0104 | +0.11 | 0.696 |

## Lung

Both lung slices are LLC from GSE155972. With n = 2 there is no Wilcoxon claim.

| Slice | Cldn4 Δ | IFN-γ Δ | NHEJ Δ |
|---|---:|---:|---:|
| LLC_GSE155972_Setdb1_KO_antiCTLA4&antiPD1 | +0.334 | +0.204 | -0.378 |
| LLC_GSE155972_antiCTLA4&antiPD1 | +0.081 | +0.304 | -0.060 |

## Secondary scores

Same 64-slice Spearman of the Cldn4 mean-delta against each secondary score delta. These p values are not in the six-test BH family.

| Score | Genes used | Slices | ρ | p | Score up |
|---|---:|---:|---:|---:|---:|
| IFN_alpha | 92 | 64 | -0.01 | 0.929 | 48/64 |
| NHEJ_Reactome | 59 | 62 | +0.15 | 0.257 | 16/62 |
| NHEJ_Reactome_no_histone | 34 | 64 | +0.14 | 0.279 | 17/64 |
| NHEJ_core | 12 | 64 | +0.08 | 0.516 | 20/64 |
| NHEJ_WP | 7 | 64 | +0.14 | 0.285 | 17/64 |

Full Reactome NHEJ is unscored when histone genes are missing and coverage falls under 80%. Unscored slices (2): EMT6_GSE107801_antiPDL1, EMT6_GSE107801_antiTGFb_trap_antiPDL1. The histone-stripped Reactome score and the c-NHEJ machinery score still cover those slices.

## Figures

- `figures/fig1_delta_cldn4_vs_scores.png` — Cldn4 change versus IFN-γ change and versus NHEJ change. Black edges mark LLC.
- `figures/fig2_score_waterfalls.png` — paired score changes across the 64 slices.
- `figures/fig3_levels_pre_post.png` — slice-level Cldn4 versus each score before ICB and after ICB. The baseline IFN-γ panel is the ρ = +0.39 association that shrinks to ρ = +0.08 across 17 cell lines.

## Sources

- TISMO Gene module, `POST https://tismo.pku-genomics.org/rtismo/gene/downVivoExprn`, type=3, all ICB treatments × all tumors. Download date is in `tables/download_manifest.json`.
- Lock: PR #542 (`results/tismo/tables/pairs_merged.tsv` on `cursor/tismo-icb-tacstd2-cldn4-tj-aa3d`), copied to `methods/tismo_cldn4_nhej_ifn/locked_pairs_pr542.tsv`.
- Mouse MSigDB gene sets (MGI symbols): Hallmark IFN-γ MM3878, Hallmark IFN-α MM3877, GO:0006303 MM4659, Reactome NHEJ MM15296 (R-MMU-5693571), WikiPathways NHEJ MM15989 (WP1242). Lists are in `methods/tismo_cldn4_nhej_ifn/gene_sets.json`.

## Reproduce

```bash
pip install -r methods/tismo_cldn4_nhej_ifn/requirements.txt
python3 methods/tismo_cldn4_nhej_ifn/download.py
python3 methods/tismo_cldn4_nhej_ifn/analyze.py
```

The analyzer refuses to write results if the recomputed Tacstd2 or Cldn4 mean-deltas disagree with the lock.
