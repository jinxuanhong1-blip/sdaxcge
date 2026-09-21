# Concordant-4: CNV+ epithelium vs annotation-only malignant

Additive. **CLDN4-only.** The locked annotation result is recomputed on the
same 65 units and left as the annotation number. CNV does not replace it.

Patient / donor / sample is the unit. Do not quote cell counts as n.
GSE123902 = donor. GSE131907 = sample. GSE205335 = patient. GSE189357 = patient.
GSE189357 has no normal epithelium on GEO; its CNV null is the GSE123902
normal-lung profile. p-values are descriptive.

## Annotation-only (same definition as the lock)

Malignant CLDN4 **%pos** vs full-unit T/NK fraction. Marker gate on
GSE123902 and GSE189357. Author malignant on GSE131907 and GSE205335.
DerSimonian–Laird on the four Fisher-z values. Recompute matches the lock
on every unit (max |Δ %pos| = 0, max |Δ T/NK| = 0).

| cohort | n | ρ | p |
|---|---:|---:|---:|
| GSE123902 | 13 | −0.659 | 0.0142 |
| GSE131907 | 21 | −0.522 | 0.0152 |
| GSE205335 | 22 | −0.435 | 0.0429 |
| GSE189357 | 9 | −0.600 | 0.0876 |
| DL meta | 65 | −0.531 | 1.65×10⁻⁵ |

I² = 0%. 95% CI −0.697 to −0.312. Stacked within-cohort Q4 vs Q1
rank-biserial **r = −0.724** (19/16, p = 2.88×10⁻⁴).

## CopyKAT-aneuploid epithelium

Epithelial cells whose CopyKAT residual (normal-epithelium corrected, chr6
out, 5 Mb bins) is aneuploid by the v1.2 Ward / Wasserstein rule and above
the normal-epithelium 95th percentile. A unit is in the Spearman only when
that set has ≥20 cells. GSE131907 `NS_12` has 22 T/NK cells, under the
30-cell reference floor, so it has no CNV call.

| cohort | n | ρ | p |
|---|---:|---:|---:|
| GSE123902 | 9 | +0.259 | 0.50 |
| GSE131907 | 12 | −0.138 | 0.67 |
| GSE205335 | 15 | −0.036 | 0.90 |
| GSE189357 | 9 | −0.033 | 0.93 |
| DL meta | 45 | −0.009 | 0.96 |

I² = 0%. 95% CI −0.336 to +0.321. Stacked Q4 vs Q1 **r = +0.114**
(10/7, p = 0.73).

On these same 45 units the annotation score is still negative:
DL **ρ = −0.382** (p = 0.021, I² = 0%, CI −0.631 to −0.061); Q4 vs Q1
r = −0.552 (13/11, p = 0.024).

Tumor units are called aneuploid more often than normal epithelium
(median fraction of drawn epithelial cells 0.20 vs 0.03). In the two
author-labeled cohorts the non-malignant epithelium in the draw is usually
below the threshold (median call fraction 0), while the median
author-malignant call fraction is 0.32 (GSE131907) and 0.27 (GSE205335).

## InferCNV+ epithelium

Same normal-epithelium 95th percentile, on the InferCNV preliminary
residual (log2 CP10k, T/NK mean removed, window 101, HMM off).

Units with ≥20 InferCNV+ cells: 2 + 3 + 12 + 2 = 19. Only GSE205335 has
n > 3, so a four-cohort meta is not formed. That cohort is
**ρ = −0.028** (n = 12, p = 0.93).

Median fraction of drawn tumor epithelium above the threshold is 0.019.
Normal-epithelium controls sit on the percentile they define
(mean call fraction 0.043, 0.044, and 0.050 in GSE123902, GSE131907, and
GSE205335).

## Both callers

CopyKAT-aneuploid and InferCNV+ together reach 20 cells in **2 / 65**
units (one GSE123902 donor, one GSE205335 patient). No cohort Spearman.

## Scope

- Draw caps: 800 epithelial / 250 T/NK per unit (GSE131907: 500 / 150).
  Annotation %pos uses every malignant cell. CNV %pos uses the called
  cells inside the draw.
- On GSE123902 and GSE189357 the epithelial gate is the marker-malignant
  definition, so CNV filters that set. On GSE131907 and GSE205335
  epithelium is wider than author-malignant.
- Chromosome 6 is excluded. The R `copykat` MCMC sampler and inferCNV HMM
  were not run. Bins are fixed 5 Mb windows.
- Not merged: GSE148071, GSE127465, GSE207422, GSE154826, GSE200563,
  E-MTAB-13526.

Figures: `results/figures/fig_cldn4_vs_tnk.png`,
`results/figures/fig_rho_compare.png`.
