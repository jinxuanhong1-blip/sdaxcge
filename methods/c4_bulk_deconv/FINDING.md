# Concordant-4 signatures on TCGA-LUAD

Additive. Public data only. CLDN4 versus an estimated CD8 fraction, compared
with keratin-adjusted bulk CLDN4 versus CD8A on the same matrix.

This is not the Stanford CIBERSORTx S-mode container and not the BayesPrism
Gibbs sampler. ν-SVR is the Newman 2015 fractions engine. Bayes MAP is a
simplex mixture with a flat Dirichlet prior and fixed scRNA profiles.
Rules: `PROTOCOL.md`.

## Cohort

TCGA-LUAD STAR TPM (Xena GDC hub), primary tumor, one sample per patient.
**n = 516.** Treatment-naive. Not an ICI endpoint.

Reference cells (class n), not the test n:

| Dataset | malignant | cd8 | cd4 | nk | b | myeloid | endothelial | fibroblast |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| GSE123902 tumor donors (normals out) | 991 | 4441 | 1812 | 4585 | 2278 | 5727 | 277 | 1029 |
| GSE131907 tLung, author subtype | 6352 | 3134 | 7307 | 500 | 5312 | 8794 | 655 | 1718 |
| GSE205335 author lineage | 28512 | 13452 | 21433 | 3880 | 9333 | 10511 | 737 | 2688 |
| GSE189357 marker gate, 9 resections | 4644 | 12564 | 8513 | 22424 | 13829 | 31483 | 2545 | 2865 |

GSE131907 malignant cells are tS1/tS2/tS3 only. GSE189357 NK/myeloid counts
are high because the marker gate calls NKG7/GNLY-positive, CD3-negative cells
as NK and LYZ/CD68/CD14-positive cells as myeloid. Profiles are an unweighted
mean of the four datasets. In that mean, CD8A peaks in cd8, MS4A1 in b, EPCAM
and KRT19 in malignant, PECAM1 in endothelial, COL1A1 in fibroblast. CLDN4
also peaks in malignant and is held out of every signature.

## Estimator check (before reading the CLDN4 claim)

Synthetic mixtures were recovered (Pearson of true vs estimated θ: ν-SVR 0.990,
NNLS 1.000, Bayes MAP 1.000). On TCGA the picture splits.

Primary signature, 867 genes, median reconstruction r = 0.796 (ν-SVR) and
0.877 (NNLS):

| Engine | CD8A vs CD8 fraction | KRT8 vs malignant fraction |
|---|---|---|
| ν-SVR | ρ = −0.086, p = 0.052 | ρ = 0.005, p = 0.91 |
| NNLS | ρ = −0.111, p = 0.012 | ρ = −0.159, p = 2.9×10⁻⁴ |
| Bayes MAP | ρ = 0.738, p = 1.0×10⁻⁸⁹ | ρ = 0.290, p = 1.8×10⁻¹¹ |

ν-SVR and NNLS do not measure CD8 here. Their CLDN4 coefficients are not a
CD8 result. Bayes MAP CD8 ranks with CD8A. The absolute Bayes CD8 fraction is
compressed (mean 0.003, 90th percentile 0.008), so it is a rank, not a tissue
percentage. ν-SVR and Bayes CD8 agree poorly (ρ = 0.029).

Author-labeled malignant cells only (GSE131907 tS1–tS3 and GSE205335), same
immune profiles: Bayes CD8A control ρ = 0.756, KRT8 control ρ = 0.565.
Holding KRT8/KRT18/KRT19/KRT7/EPCAM out of the feature list leaves the Bayes
CD8A control at ρ = 0.708.

## CLDN4 vs CD8

Same 516 tumors. Spearman, two-sided. Partial Spearman residualizes midranks.
BH q is across the pre-specified biological family (bulk, all three engines,
both sensitivities). Intervals are Fisher z.

| Contrast | ρ | 95% interval | p | q |
|---|---:|---|---:|---:|
| Bulk CLDN4 vs CD8A | −0.098 | −0.183 to −0.012 | 0.026 | 0.33 |
| Bulk CLDN4 vs CD8A given KRT8+KRT18+KRT19 | −0.087 | −0.172 to −0.0003 | 0.049 | 0.33 |
| Bayes CD8 vs CLDN4 | −0.044 | −0.129 to 0.043 | 0.32 | 0.78 |
| Bayes CD8 vs CLDN4 given malignant fraction | 0.002 | −0.085 to 0.088 | 0.97 | 0.99 |
| Bayes CD8 vs CLDN4 given KRT8+KRT18+KRT19 | −0.056 | −0.142 to 0.031 | 0.21 | 0.70 |
| Author-malignant Bayes CD8 vs CLDN4 | −0.087 | −0.172 to −0.001 | 0.049 | 0.33 |
| Author-malignant Bayes, given malignant fraction | −0.035 | — | 0.42 | 0.86 |
| Author-malignant Bayes, given KRT8+KRT18+KRT19 | −0.041 | — | 0.35 | 0.78 |
| CLDN4 vs KRT8 given only KRT18+KRT19 | 0.075 | −0.012 to 0.160 | 0.090 | 0.40 |

Keratin genes barely move the bulk correlation (−0.098 to −0.087). The CD8
fraction that actually tracks CD8A does not reproduce that correlation on the
primary signature. The author-malignant sensitivity matches the bulk number
(ρ = −0.087) and loses it after the malignant fraction or the keratins are
partialled out. No CLDN4 test in this family has q below 0.3.

ν-SVR CLDN4 vs CD8 is ρ = 0.007. That number is withheld as a biological
result because the CD8A control failed.

## What this does not say

It does not replace the locked seven-cohort keratin table, and it does not
re-estimate that q. The LUAD specificity contrast here (CLDN4 vs KRT8, both
adjusted only for KRT18 and KRT19) is ρ = 0.075, p = 0.090 on STAR TPM.
TCGA is not ICI. GSE205335 contributes cells from lymph-node and liver
metastases. Marker-gated malignant cells in GSE123902 and GSE189357 include
epithelium that was not called malignant by an author. Bayes CD8 fractions
are not cell percentages.
