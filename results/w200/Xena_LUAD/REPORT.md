# UCSC Xena TCGA-LUAD: TACSTD2 / CLDN4 vs immune after ESTIMATE and ABSOLUTE purity

**Cohort:** TCGA-LUAD primary tumors (`-01`). HiSeqV2 RNA-seq n = 515. **Primary analysis set n = 502** with HiSeqV2 + MD Anderson ESTIMATE + PanCanAtlas ABSOLUTE. 500 also have methylation leukocyte fraction; 444 have Wolf signature scores. STAR TPM sensitivity: 516 primaries, 502 with both purities.

**Question:** Are TACSTD2 (TROP2) and CLDN4 associated with immune features in UCSC Xena TCGA-LUAD after adjusting for tumor purity, using both ESTIMATE (RNA) and ABSOLUTE (DNA)?

This is not a re-run of `results/w200/A1_LUAD/` (TACSTD2 only, ABSOLUTE only). It adds CLDN4, ESTIMATE vs ABSOLUTE side-by-side, and a GDC STAR TPM sensitivity matrix.

## Verdict: WEAK / PARTIAL SUPPORT — not “immune cold”

Both genes show a **small** cytotoxicity-poor, Treg / TGF-beta-high skew that survives either purity method. Effect sizes are |rho| ≤ 0.22. There is **no association with total immune infiltration** on the primary HiSeqV2 matrix (methylation leukocyte fraction: TACSTD2 |ABS rho = 0.03, FDR = 0.57; CLDN4 |ABS rho = −0.09, FDR = 0.10). ESTIMATE Immune_score, GEP18, CD274/PD-L1, FOXP3, and CXCL9/10 are null for TACSTD2 after ABSOLUTE adjustment.

Calling TACSTD2-high or CLDN4-high LUAD “immune excluded” or “immune cold” from these data would be overstating it.

Purity adjustment is almost a no-op for both genes: **TACSTD2 is uncorrelated with purity** (vs ESTIMATE rho = −0.006, p = 0.90; vs ABSOLUTE rho = 0.007, p = 0.87). **CLDN4 is only weakly related** (vs ESTIMATE rho = 0.082, p = 0.064; vs ABSOLUTE rho = 0.048, p = 0.29). The immune features themselves *are* strongly purity-associated (leukocyte fraction vs ABSOLUTE rho = −0.77; CYT vs ABSOLUTE rho = −0.47). The confounder is real; these two genes are not it.

## Purity machinery works; these genes are orthogonal to it

| Pair (HiSeqV2) | n | Spearman rho | p |
|---|---:|---:|---:|
| ESTIMATE purity vs ABSOLUTE purity | 502 | 0.66 | 8e−64 |
| MDACC ESTIMATE purity vs Aran ESTIMATE | 511 | 1.00 | 0 |
| ABSOLUTE vs Aran ABSOLUTE | 352 | 0.92 | 3e−147 |
| ESTIMATE purity vs leukocyte fraction | 513 | −0.73 | 3e−85 |
| ABSOLUTE purity vs leukocyte fraction | 500 | −0.77 | 4e−100 |
| ESTIMATE Immune_score vs ESTIMATE purity | 515 | **−0.93** | ~0 |
| ESTIMATE Immune_score vs ABSOLUTE purity | 502 | −0.62 | 8e−54 |
| TACSTD2 vs ESTIMATE / ABSOLUTE purity | 515 / 502 | −0.006 / 0.007 | 0.90 / 0.87 |
| CLDN4 vs ESTIMATE / ABSOLUTE purity | 515 / 502 | 0.082 / 0.048 | 0.064 / 0.29 |
| TACSTD2 vs CLDN4 | 515 | **0.46** | 3e−28 |

The Immune_score vs ESTIMATE-purity rho of −0.93 is the pre-declared circularity: Immune_score is a term in ESTIMATE_score, and purity is a cosine transform of that score. Partialling ESTIMATE purity out of Immune_score is not a valid test. ABSOLUTE (DNA) and methylation leukocyte fraction are the non-circular readouts.

## Key results — TACSTD2, partial Spearman, HiSeqV2

Significant after ABSOLUTE adjustment (BH-FDR < 0.05 across 37 features):

| Feature | Unadj. rho | \| ESTIMATE rho | \| ABSOLUTE rho | ABS FDR |
|---|---:|---:|---:|---:|
| Wolf TGF-beta response | 0.18 | 0.20 | **0.22** | 1.7e−04 |
| GZMB | −0.14 | −0.16 | **−0.16** | 7.4e−03 |
| Wolf IFN-gamma signature | 0.12 | 0.13 | **0.15** | 0.022 |
| GZMA | −0.12 | −0.15 | **−0.14** | 0.022 |
| LAG3 | −0.11 | −0.13 | **−0.12** | 0.041 |
| CIBERSORT Tregs | 0.12 | 0.12 | **0.12** | 0.041 |
| CIBERSORT CD8 T | −0.12 | −0.12 | **−0.12** | 0.041 |
| CYT (mean GZMA/PRF1) | −0.11 | −0.14 | **−0.12** | 0.041 |
| CD8 score (CD8A/CD8B) | −0.11 | −0.13 | **−0.12** | 0.041 |

Notably **not** significant after ABSOLUTE: leukocyte fraction (0.03), ESTIMATE Immune_score (0.07), GEP18 (−0.03), CD274 (0.05), PDCD1, FOXP3, CXCL9/10, IDO1, CD3E. CD8A, IFNG, PRF1, CTLA4 sit at FDR 0.058.

ESTIMATE adjustment and ABSOLUTE adjustment agree to two decimal places for this gene, as expected once TACSTD2 ⟂ both purities. Dual adjustment (both purities at once) does not change the picture (TGF-beta 0.22; GZMB −0.17; leukocyte fraction 0.02).

## Key results — CLDN4, partial Spearman, HiSeqV2

Significant after ABSOLUTE (FDR < 0.05):

| Feature | Unadj. rho | \| ESTIMATE rho | \| ABSOLUTE rho | ABS FDR |
|---|---:|---:|---:|---:|
| GZMA | −0.17 | −0.17 | **−0.17** | 0.003 |
| GZMB | −0.17 | −0.16 | **−0.17** | 0.003 |
| CYT | −0.16 | −0.15 | **−0.15** | 0.007 |
| CIBERSORT M2 macrophages | 0.15 | 0.15 | **0.15** | 0.007 |
| IFNG gene | −0.15 | −0.14 | **−0.15** | 0.007 |
| CIBERSORT Tregs | 0.13 | 0.15 | **0.14** | 0.014 |
| PRF1 | −0.13 | −0.11 | **−0.12** | 0.027 |
| CD3E | −0.13 | −0.12 | **−0.12** | 0.027 |
| Wolf TGF-beta | 0.09 | 0.15 | **0.13** | 0.027 |
| NKG7 | −0.13 | −0.11 | **−0.12** | 0.028 |
| LAG3 | −0.13 | −0.10 | **−0.12** | 0.028 |
| Wolf wound-healing CSR | 0.12 | 0.12 | **0.12** | 0.028 |
| Wolf lymphocyte infiltration | −0.13 | −0.10 | **−0.12** | 0.042 |

Not significant after ABSOLUTE: leukocyte fraction (−0.09, FDR 0.10), CD8A / CD8 score / CIBERSORT CD8, GEP18, CD274, ESTIMATE Immune_score. Same small-effect ceiling as TACSTD2 (max |rho| = 0.17).

## STAR TPM sensitivity — direction holds, FDR set does not

TACSTD2 vs CLDN4 is a bit stronger on STAR (rho = 0.54). Partial-rho vectors for TACSTD2 | ABSOLUTE agree in sign for 32/37 features (Spearman of the two rho vectors = 0.92). That is the stable part.

The FDR-significant set moves:

- STAR TACSTD2 | ABSOLUTE picks up Wolf lymphocyte infiltration (−0.18), CIBERSORT M2 (0.17), CD274 (0.16), HAVCR2 (0.13) and drops CYT / CD8 score / Tregs / LAG3 below FDR 0.05 (they stay ~−0.10).
- STAR CLDN4 | ABSOLUTE **does** reach FDR < 0.05 for leukocyte fraction (rho = −0.18). HiSeqV2 does not (−0.09, FDR 0.10). CLDN4 also tracks ESTIMATE purity more on STAR (rho = 0.14, p = 0.001) than on HiSeqV2 (0.08, p = 0.064).

Honest reading: a small cytotoxic-low / TGF-beta-high skew is reproducible across quantification pipelines. Which exact gene or score crosses FDR 0.05 is not. Do not quote the STAR leukocyte-fraction hit for CLDN4 without the HiSeqV2 miss, or vice versa.

## Honest interpretation

1. **This is not a purity artifact of TACSTD2 or CLDN4.** Both genes are essentially orthogonal to ESTIMATE and ABSOLUTE in LUAD. Adjusting for either purity barely moves the estimates. Immune *features* are purity-associated; these two *predictors* are not.
2. **The residual association is real, small, and only partly coherent.** Cytotoxic effectors (GZMA/GZMB/CYT) and Tregs / TGF-beta point the same way for both genes. Total immune mass (leukocyte fraction, ESTIMATE Immune_score, GEP18) does not. A cytotoxicity-poor, TGF-beta-high *skew* at |rho| ≈ 0.12–0.22 is the most the data support. That is <5% of rank variance.
3. **Wolf IFN-gamma signature vs IFNG gene still disagree** (TACSTD2 |ABS: signature +0.15, gene −0.10). Same inconsistency as A1. Different measurements, both small; do not cherry-pick.
4. **ESTIMATE Immune_score is not an independent immune readout here.** It is −0.93 with ESTIMATE purity by construction. After ABSOLUTE adjustment it is null vs both genes. Use leukocyte fraction or gene-level CYT/CD8 if the question is “how much immune.”
5. **CLDN4 is not a cleaner immune-exclusion marker than TACSTD2.** Co-expression is moderate (0.46–0.54). CLDN4’s cytotoxic-gene rhos are slightly larger; its leukocyte-fraction result is pipeline-dependent. Neither gene is a strong anti-immune correlate.

## Methods

- Expression, primary: UCSC Xena `TCGA.LUAD.sampleMap/HiSeqV2`, log2(RSEM normalized count + 1); `-01` primaries.
- Expression, sensitivity: UCSC Xena GDC hub `TCGA-LUAD.star_tpm.tsv.gz`, log2(TPM+1), GENCODE v36 symbols, aliquots averaged to patient `-01`.
- ESTIMATE purity: MD Anderson LUAD RNAseqV2 scores; `purity = cos(0.6049872018 + 0.0001467884 · ESTIMATE_score)` (Yoshihara 2013). Formula matches Aran 2015 ESTIMATE column at rho = 1.00.
- ABSOLUTE purity: PanCanAtlas `TCGA_mastercalls.abs_tables_JSedit.fixed.txt`.
- Features (37): methylation leukocyte fraction; ESTIMATE Immune/Stromal scores; 5 Wolf signatures; 5 CIBERSORT relative fractions; CYT, CD8 score, Ayers GEP18; 21 marker genes.
- Unadjusted: Spearman. Adjusted: partial Spearman (Pearson on rank residuals), t-test with n−2−k df. BH-FDR within each (matrix × predictor × purity method) block. Dual-purity model is sensitivity only.
- Reproduce: `bash scripts/xena_luad_download.sh && python3 scripts/xena_luad_tacstd2_cldn4_immune.py`

## Caveats

- Bulk-RNA observational association; not causal, not protein, not ICI outcome.
- ESTIMATE purity is Affymetrix-calibrated and RNA-derived. It is a valid *comparison* confounder, not an independent one.
- CIBERSORT fractions are relative to leukocyte content, unfiltered by CIBERSORT p-value.
- Wolf signatures come from a different pan-cancer expression freeze (n drops to 444).
- Single-covariate rank residual cannot remove all compositional confounding. Given TACSTD2/CLDN4 ⟂ purity, residual purity confounding is unlikely to be the story for these two genes.

Tables: `tables/correlations.tsv`, `tables/purity_context.tsv`, per-sample matrices in `tables/`. Figures: `figures/`. Provenance (URLs, sha256): `provenance.json`.
