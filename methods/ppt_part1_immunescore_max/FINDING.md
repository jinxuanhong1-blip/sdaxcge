# PPT Part 1 — maximum |ρ| vs ImmuneScore

Public LUAD bulk only. The 221 genes and their order are the PR 590 CLDN4-high malignant signature (`methods/cldn4_high_malignant_signature/tables/signature_genes.tsv`, SHA-256 `9d398e71f8a2569d3c5c7bb940577a24eaca1f19db0a1d589c84e13a0bac6a3a`). CLDN4 and TACSTD2 are held out of that list. Prefix length is the only size knob. No gene was added or removed because of ImmuneScore. GSE10072, GSE11969, and GSE248378 stay closed. Bulk ρ is not a spatial exclusion result and does not replace the locked CosMx ratios or the concordant-4 malignant CLDN4 % vs T/NK result.

ImmuneScore is the mean of within-cohort z-scores of CD8A, GZMA, GZMB, IFNG, EOMES, CXCL9, CXCL10, and TBX21. None of these genes are in the 221. ESTIMATE immune ssGSEA (α=0.25, common-gene filter, same engine as the stromal covariate) is reported separately and is not called ImmuneScore.

## Objective

Primary number: absolute DerSimonian–Laird meta-analytic Spearman versus ImmuneScore. Studies are OncoSG, GSE273377 (discovery and validation inverse-variance combined first), GSE282774, and GSE233774 tumors. Eligible specs have a finite correlation in all five strata.

Three feature classes are maximized separately, then the largest of the three is the Part 1 maximum:

1. Signature score (z-mean or ssGSEA; α in {0, 0.25, 0.75, 1}; prefix size 5–221).
2. TACSTD2 itself. Spearman is invariant to a monotone transform, so the gene has no scale knob.
3. Equal-weight mean of the within-cohort z-score of a signature score and the within-cohort z-score of TACSTD2. Weights are not fit to ImmuneScore.

Covariate modes inside each class: unadjusted; partial on published PURITY or ESTIMATE StromalScore; partial on KRT8, KRT18, and KRT19; partial on purity/stroma plus those keratins; gene-level residual on purity/stroma. The p-value on a selected row is the p-value of a maximized |ρ|.

## Check against published cells

All 7 locked cells matched, including OncoSG z-mean size 221 versus CD8A at the stored PR 590 precision and the PR 693 ImmuneScore cohort rhos at ssGSEA α=0.75, size 163 (tolerance 0.0015).

ImmuneScore gene coverage: OncoSG uses 8/8; GSE273377 discovery uses 8/8; GSE273377 validation uses 8/8; GSE233774 tumor uses 7/8 (IFNG absent); GSE282774 uses 8/8.

## Part 1 maximum

Largest |meta ρ| versus ImmuneScore is **-0.585** (fixed, p=1.78e-40, I²=0%, n sum=420, k=4). Spec: **ssGSEA α=0.75, size 163, partial | KRT8/18/19, vs ImmuneScore**.

OncoSG -0.626 (p=1.77e-19, n=169); GSE273377 discovery -0.563 (p=1.06e-09, n=103); GSE273377 validation -0.476 (p=0.000185, n=60); GSE282774 -0.599 (p=1.32e-06, n=58); GSE233774 tumor -0.578 (p=0.0016, n=30).

## Signature alone

Signature-only maximum: **-0.585** (fixed, p=1.78e-40, I²=0%, n sum=420). Spec: ssGSEA α=0.75, size 163, partial | KRT8/18/19, vs ImmuneScore.

OncoSG -0.626 (p=1.77e-19, n=169); GSE273377 discovery -0.563 (p=1.06e-09, n=103); GSE273377 validation -0.476 (p=0.000185, n=60); GSE282774 -0.599 (p=1.32e-06, n=58); GSE233774 tumor -0.578 (p=0.0016, n=30).

Pre-specified baseline (z-mean, size 221, unadjusted): meta ρ -0.460 (DL random, p=7.99e-07, I²=71%, n sum=420).

Same signature score across covariate modes:

keratin: -0.585 (fixed, I²=0%); none: -0.574 (fixed, I²=0%); purity_keratin: -0.475 (fixed, I²=0%); purity: -0.462 (fixed, I²=0%); residual: -0.356 (fixed, I²=0%).

Unadjusted correlation at this same score is meta ρ -0.574 (fixed, I²=0%). That unadjusted score is the PR 693 ImmuneScore maximum. Partial correlation on KRT8/18/19 is the covariate mode that increases |ρ| over that unadjusted score. Partial correlation on purity/stroma decreases |ρ|.

## TACSTD2 alone

TACSTD2 maximum: **-0.278** (DL random, p=7.30e-05, I²=42%, n sum=420). Spec: TACSTD2, partial | KRT8/18/19, vs ImmuneScore.

OncoSG -0.385 (p=2.92e-07, n=169); GSE273377 discovery -0.011 (p=0.912, n=103); GSE273377 validation -0.351 (p=0.00749, n=60); GSE282774 -0.184 (p=0.178, n=58); GSE233774 tumor -0.398 (p=0.0399, n=30).

Pre-specified TACSTD2 unadjusted ImmuneScore meta ρ is -0.268 (DL random, p=0.0029, I²=65%).

At the TACSTD2 maximum, 3/5 strata are inverse. Flat strata: GSE273377 discovery -0.011 (p=0.912), GSE282774 -0.184 (p=0.178). Flat and positive strata stay in the meta. TACSTD2 is not a concordant cold result on these five strata.

## Equal-weight signature + TACSTD2

Equal-weight maximum: **-0.542** (DL random, p=4.48e-21, I²=30%, n sum=420). Spec: equal-weight signature + TACSTD2, ssGSEA α=1, size 205, partial | KRT8/18/19, vs ImmuneScore.

OncoSG -0.625 (p=2.26e-19, n=169); GSE273377 discovery -0.431 (p=7.67e-06, n=103); GSE273377 validation -0.568 (p=3.97e-06, n=60); GSE282774 -0.449 (p=0.000576, n=58); GSE233774 tumor -0.553 (p=0.00277, n=30).

Equal-weight with TACSTD2 does not raise |meta ρ| above the signature-only maximum (-0.542 vs -0.585).

## Held-out OncoSG

The GEO-only maximum inside the ImmuneScore search (OncoSG not used to pick the spec) is ssGSEA α=1, size 157, unadjusted, vs ImmuneScore. GEO meta ρ -0.564 (fixed, p=2.96e-23, I²=0%, n sum=251).

Applied to OncoSG: ρ=-0.583 (p=8.58e-17, n=169).

Frozen GEO spec, all five strata (OncoSG was not used to choose it): OncoSG -0.583 (p=8.58e-17, n=169); GSE273377 discovery -0.571 (p=2.98e-10, n=103); GSE273377 validation -0.472 (p=0.000142, n=60); GSE282774 -0.624 (p=1.66e-07, n=58); GSE233774 tumor -0.592 (p=0.000569, n=30).

That frozen spec has four-study meta ρ -0.572 (fixed, p=2.02e-39, I²=0%, n sum=420).

## ESTIMATE immune ssGSEA

The largest |meta ρ| versus ESTIMATE immune ssGSEA is -0.544 (DL random, p=2.44e-21, I²=32%, n sum=420). Spec: ssGSEA α=0.25, size 10, unadjusted, vs ESTIMATE immune ssGSEA. This row is not relabeled ImmuneScore.

## What the 221 genes are

The list is the CLDN4-tracking malignant program: concordant-4 CLDN4-high versus CLDN4-low genes that still track CLDN4 in TCGA-LUAD after KRT8, KRT18, KRT19, and ABSOLUTE purity. Top gene CLDN3 (family TJ, TCGA partial ρ vs CLDN4 +0.600). Genes labeled TJ in the signature table: 4/221 (CLDN3, DLG3, TBCD, MAP2K7). CLDN4 in the list: no. TACSTD2 in the list: no. Overlap with ImmuneScore genes: none. Overlap with the ESTIMATE immune set: none. Classical junction-panel genes present: CLDN3. This score is not a KEGG tight-junction module.

## Plateau

Largest eligible ImmuneScore |meta ρ| values:

| feature | spec | covariate | meta ρ | I² | model |
|---|---|---|---:|---:|---|
| signature | ssGSEA α=0.75, size 163 | keratin | -0.585 | 0% | fixed |
| signature | ssGSEA α=0.75, size 157 | keratin | -0.584 | 0% | fixed |
| signature | ssGSEA α=0.75, size 156 | keratin | -0.584 | 0% | fixed |
| signature | ssGSEA α=0.25, size 159 | keratin | -0.584 | 10% | DL random |
| signature | ssGSEA α=0.25, size 161 | keratin | -0.584 | 14% | DL random |
| signature | ssGSEA α=0.25, size 162 | keratin | -0.584 | 12% | DL random |
| signature | ssGSEA α=0.75, size 159 | keratin | -0.584 | 0% | fixed |
| signature | ssGSEA α=0.25, size 160 | keratin | -0.584 | 8% | DL random |

## Reproduction

```bash
python3 methods/ppt_part1_immunescore_max/analyze.py
```

Matrices download to `/tmp/cldn4sig` and are not committed. Purity covariate: OncoSG uses the cBioPortal published PURITY column. GEO uses ESTIMATE stromal ssGSEA (α=0.25). Partial correlation is Pearson of rank residuals.
