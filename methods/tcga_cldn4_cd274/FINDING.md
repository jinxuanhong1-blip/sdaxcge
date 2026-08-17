# Finding — TCGA-LUAD / TCGA-LUSC RNA: CLDN4 vs CD274 and HLA-A/B/C, partial on ESTIMATE ImmuneScore

**Additive public RNA only.** UCSC Xena legacy **HiSeqV2** log2(RSEM norm_count+1) for **TCGA-LUAD** and **TCGA-LUSC** primary tumors (`-01`). Official MD Anderson **ESTIMATE RNAseqV2 ImmuneScore** (Yoshihara 2013). This folder does not re-audit prior TACSTD2 / purity PRs.

Primary question: **CLDN4 vs CD274 (PD-L1)** and **CLDN4 vs HLA-A / HLA-B / HLA-C**, unadjusted and after residualizing on **ESTIMATE ImmuneScore**. MHC-I = mean(HLA-A, HLA-B, HLA-C) is a same-run companion.

## Honest n

| Filter | TCGA-LUAD | TCGA-LUSC |
|---|---:|---:|
| HiSeqV2 primary tumors (`-01`, 15-char) | 515 | 502 |
| Official ESTIMATE RNAseqV2 primaries | 515 | 501 |
| **Complete: CLDN4 + CD274 + HLA-A/B/C + ImmuneScore** | **515** | **501** |

Primary tests use the complete-case n. No imputation. One row per 15-character barcode; replicate `-01` aliquots averaged.

## Verdict table

| Cohort | Endpoint | n | Unadj ρ | Unadj p | Partial ρ \| ImmuneScore | Partial p |
|---|---|---:|---:|---:|---:|---:|
| TCGA-LUAD | CD274 | 515 | +0.045 | 0.312 | +0.112 | 0.011 |
| TCGA-LUAD | HLA-A | 515 | +0.059 | 0.180 | +0.109 | 0.014 |
| TCGA-LUAD | HLA-B † | 515 | +0.035 | 0.424 | +0.114 | 0.010 |
| TCGA-LUAD | HLA-C | 515 | +0.023 | 0.596 | +0.076 | 0.084 |
| TCGA-LUAD | MHC_I † | 515 | +0.042 | 0.347 | +0.108 | 0.014 |
| TCGA-LUSC | CD274 | 501 | +0.023 | 0.608 | +0.027 | 0.554 |
| TCGA-LUSC | HLA-A | 501 | +0.032 | 0.474 | +0.046 | 0.309 |
| TCGA-LUSC | HLA-B † | 501 | +0.007 | 0.880 | +0.015 | 0.733 |
| TCGA-LUSC | HLA-C | 501 | +0.023 | 0.610 | +0.035 | 0.431 |
| TCGA-LUSC | MHC_I † | 501 | +0.019 | 0.676 | +0.032 | 0.476 |

† **HLA-B is in Yoshihara Immune141.** HLA-A, HLA-C, CD274, and CLDN4 are not. MHC-I includes HLA-B, so that companion partial is partly circular. The HLA-B row is reported because it was requested; it is not an independent infiltrate control.

## What holds / what does not

**Unadjusted: null in both histologies.** CLDN4 is not associated with CD274 or HLA-A/B/C on the raw HiSeqV2 matrix. LUAD |ρ| ≤ 0.059 (all p ≥ 0.18). LUSC |ρ| ≤ 0.032 (all p ≥ 0.47).

**After ImmuneScore: small LUAD residual; LUSC stays null.** LUAD n=515: CLDN4 vs CD274 +0.045 (p=0.312) → partial +0.112 (p=0.011); vs HLA-A +0.059 (p=0.180) → partial +0.109 (p=0.014); vs HLA-C +0.023 (p=0.596) → partial +0.076 (p=0.084). CD274 and HLA-A partials survive BH q=0.036 across the 8 primary tests. HLA-C does not (q=0.168). Effect size is |ρ| ≈ 0.11 (~1% of rank variance). LUSC n=501 remains null on every endpoint (all partial p ≥ 0.31).

The LUAD residual is a **suppressor** pattern, not an immune-cold claim. CLDN4 vs ImmuneScore is weakly negative and non-significant (ρ=−0.065, p=0.14). CD274 and HLA track ImmuneScore strongly (ρ=+0.52 to +0.70). Removing infiltrate unmasks a small positive CLDN4–CD274 / CLDN4–HLA-A association that the unadjusted test does not show. Do not quote the partial without the unadjusted null.

**Do not over-read HLA-B / MHC-I partials.** HLA-B ∈ Immune141. Residualizing ImmuneScore out of HLA-B subtracts a score that already contains HLA-B. HLA-A, HLA-C, and CD274 are the non-circular endpoints.

LUSC dropped **1** HiSeqV2 primary that has no official ESTIMATE row (502 → 501). That is the honest analysis n, not 502.

## Context: each gene vs ImmuneScore

| Cohort | Pair | n | ρ | p |
|---|---|---:|---:|---:|
| TCGA-LUAD | CLDN4 vs ImmuneScore | 515 | -0.065 | 0.143 |
| TCGA-LUAD | CD274 vs ImmuneScore | 515 | +0.641 | 5.79e-61 |
| TCGA-LUAD | HLA-A vs ImmuneScore | 515 | +0.519 | 8.28e-37 |
| TCGA-LUAD | HLA-B vs ImmuneScore | 515 | +0.703 | 5.27e-78 |
| TCGA-LUAD | HLA-C vs ImmuneScore | 515 | +0.590 | 1.40e-49 |
| TCGA-LUAD | MHC_I vs ImmuneScore | 515 | +0.641 | 5.61e-61 |
| TCGA-LUSC | CLDN4 vs ImmuneScore | 501 | -0.004 | 0.925 |
| TCGA-LUSC | CD274 vs ImmuneScore | 501 | +0.380 | 1.28e-18 |
| TCGA-LUSC | HLA-A vs ImmuneScore | 501 | +0.645 | 3.43e-60 |
| TCGA-LUSC | HLA-B vs ImmuneScore | 501 | +0.759 | 3.66e-95 |
| TCGA-LUSC | HLA-C vs ImmuneScore | 501 | +0.685 | 1.37e-70 |
| TCGA-LUSC | MHC_I vs ImmuneScore | 501 | +0.730 | 1.46e-84 |

## Methods

- **Matrix:** UCSC Xena `TCGA.{LUAD,LUSC}.sampleMap/HiSeqV2`, log2(RSEM normalized count + 1). Primary tumors only (`sample type 01`). Aliquots collapsed to the 15-character barcode by mean.
- **ESTIMATE:** official MD Anderson RNAseqV2 tables (`Immune_score`). Not recomputed. Matched to the same RNAseqV2 freeze as HiSeqV2.
- **Genes:** CLDN4, CD274, HLA-A, HLA-B, HLA-C — all present on HiSeqV2. MHC-I companion = unweighted mean of the three HLA log2 values (not a gene-set ssGSEA).
- **Unadjusted:** Spearman. **Partial:** first-order partial Spearman (algebraic) of CLDN4 vs endpoint controlling for ImmuneScore; df = n − 3. Rank-residual Pearson is a sensitivity column in `tables/correlations.tsv` and agrees to three decimals here.
- **CI:** Fisher z, variance 1/(n−3) unadjusted and 1/(n−4) partial.
- **FDR:** BH across the 8 primary tests (2 cohorts × CD274/HLA-A/HLA-B/HLA-C). MHC-I is not in the FDR set.

| Test | LUAD q unadj / q partial | LUSC q unadj / q partial |
|---|---|---|
| CLDN4 vs CD274 | 0.697 / 0.036 | 0.697 / 0.633 |
| CLDN4 vs HLA-A | 0.697 / 0.036 | 0.697 / 0.495 |
| CLDN4 vs HLA-B | 0.697 / 0.036 | 0.880 / 0.733 |
| CLDN4 vs HLA-C | 0.697 / 0.168 | 0.697 / 0.574 |

## What is not done

- No TACSTD2 restatement. No ICI outcome. No protein / IHC. No STAR-TPM sensitivity (different freeze than official ESTIMATE RNAseqV2).
- No ABSOLUTE / methylation leukocyte-fraction residual (those are prior Xena PRs). The requested covariate is **ImmuneScore**.
- ImmuneScore is RNA-derived. Partialling it does not equal adjusting for DNA purity.

## Files

- `tables/correlations.tsv` — all n / ρ / p / CI
- `tables/counts.tsv` — honest-n filters
- `tables/samples.tsv` — per-sample genes + ImmuneScore
- `tables/provenance.json` — URLs and sha256
- `figures/forest_partial.png` — unadj vs ImmuneScore-partial ρ
- Reproduce: `python3 methods/tcga_cldn4_cd274/analyze.py`

