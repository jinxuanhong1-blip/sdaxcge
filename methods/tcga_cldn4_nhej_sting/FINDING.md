# Finding — TCGA-LUAD / TCGA-LUSC: CLDN4-low vs ssGSEA NHEJ, IFN-γ, STING, and APM

**Additive public RNA.** UCSC Xena legacy **HiSeqV2** log2(RSEM norm_count+1) for primary tumors (`-01`). Official MD Anderson **ESTIMATE RNAseqV2 ImmuneScore**. Single-sample GSEA is gseapy 1.3.1 (rank normalization, weight 0.25, no gene-set permutation). This folder does not restate prior TACSTD2, purity, or keratin-rank results.

Questions, pre-specified: in **CLDN4-low** tumors (bottom quartile, Q1) versus **CLDN4-high** (top quartile, Q4), is the **NHEJ** ssGSEA score lower, and are **IFN-γ / STING / MHC-I APM** scores higher? Separately, Spearman of **CLDN4 vs PRKDC, LIG4, STING1, and HLA-A/B/C**.

## Honest n

| Filter | TCGA-LUAD | TCGA-LUSC |
|---|---:|---:|
| HiSeqV2 primary tumors (`-01`, 15-char) | 515 | 502 |
| Official ESTIMATE RNAseqV2 primaries | 515 | 501 |
| **Complete: CLDN4 + scores + PRKDC/LIG4/TMEM173/HLA + ImmuneScore** | **515** | **501** |
| Q1 / Q4 (barcode order breaks CLDN4 ties) | 129 / 129 | 126 / 125 |

Quartile sizes differ by one tumor when n is not divisible by 4 (LUAD 129/129/128/129; LUSC 126/125/125/125). LUSC drops the one HiSeqV2 primary with no official ESTIMATE row (502 → 501). That is the same complete-case n as the earlier CLDN4–HLA RNA table. A regression check reproduces those unadjusted HLA Spearman values (LUAD HLA-A ρ=+0.059, HLA-B ρ=+0.035, HLA-C ρ=+0.023; LUSC HLA-A ρ=+0.032, HLA-B ρ=+0.007, HLA-C ρ=+0.023).

## Answer

**CLDN4-low tumors do not have a lower NHEJ score, and they do not have higher IFN-γ, STING, or APM scores.** NHEJ moves the other way in both histologies: the bottom CLDN4 quartile has a higher NHEJ score than the top quartile. IFN-γ, the STING pathway score, and the MHC-I APM score are null on the unadjusted test. Absolute ssGSEA NES is not an enrichment p-value, and NES levels are not comparable across gene sets; the contrasts below are within a cohort.

| Question | TCGA-LUAD (n=515) | TCGA-LUSC (n=501) |
|---|---|---|
| NHEJ lower in CLDN4-low? | No. Score is higher. ρ=-0.177 (q=4.20e-04); Q1−Q4 Δ=+0.020, Cliff's δ=+0.28, MW q=7.05e-04 | No. Score is higher. ρ=-0.141 (q=0.006); Q1−Q4 Δ=+0.012, Cliff's δ=+0.19, MW q=0.036 |
| IFN-γ higher in CLDN4-low? | No. ρ=+0.020 (q=0.838); Q1−Q4 Δ=-0.003, Cliff's δ=+0.02, MW q=0.831 | No. ρ=+0.004 (q=0.925); Q1−Q4 Δ=+0.017, Cliff's δ=+0.05, MW q=0.551 |
| STING score higher in CLDN4-low? | No. ρ=+0.090 (q=0.110); Q1−Q4 Δ=-0.017, Cliff's δ=-0.08, MW q=0.491 | No. ρ=+0.071 (q=0.221); Q1−Q4 Δ=-0.011, Cliff's δ=-0.06, MW q=0.528 |
| APM score higher in CLDN4-low? | No. ρ=-0.042 (q=0.544); Q1−Q4 Δ=+0.010, Cliff's δ=+0.08, MW q=0.491 | No. ρ=-0.015 (q=0.838); Q1−Q4 Δ=+0.019, Cliff's δ=+0.06, MW q=0.528 |

**NHEJ is higher in CLDN4-low, not lower.** LUAD Spearman ρ=-0.177 (q=4.20e-04); LUSC ρ=-0.141 (q=0.006). ImmuneScore partials stay negative (LUAD -0.209, LUSC -0.147), and so do partials on mean KRT8/KRT18/KRT19 (LUAD -0.182, p=3.24e-05; LUSC -0.149, p=8.46e-04). MKI67 does not explain the LUAD result (partial ρ=-0.130, p=0.003). In LUSC the score association shrinks after MKI67 (partial ρ=-0.058, p=0.198), so the LUSC NHEJ-score shift is partly shared with a proliferation marker. The NES scale is tight (within-cohort SD ≈ 0.05), so the median shifts (+0.020 LUAD, +0.012 LUSC) are modest rank effects (Cliff's δ = +0.28 in LUAD and +0.19 in LUSC), not a separated population.

**IFN-γ is not higher in CLDN4-low.** Unadjusted ρ=+0.020 in LUAD and +0.004 in LUSC; both quartile tests are null. The LUAD ImmuneScore partial is a small positive residual (ρ=+0.130, q=0.006): after removing infiltrate, higher CLDN4, not lower CLDN4, sits slightly higher on IFN-γ. The unadjusted test is null, and IFN-γ vs ImmuneScore itself is ρ=+0.826, so that residual is a suppressor, not an IFN-high CLDN4-low result. The LUSC partial is null (ρ=+0.017).

**The STING pathway score is not higher in CLDN4-low.** Unadjusted LUAD ρ=+0.090 (q=0.110) and LUSC ρ=+0.071 (q=0.221) do not survive BH, and the Q1 vs Q4 tests are null. The LUAD ImmuneScore partial is positive (ρ=+0.162, q=9.44e-04), which is the other direction. The 11-gene set mixes TMEM173/TBK1/IRF3/IFI16 with negative regulators (TREX1, NLRC3, NLRP4, DTX4), so it is a membership score, not an induced-ISG score. The direct STING1 readout is the TMEM173 row below.

**A sensitivity, not a primary endpoint:** the 8 APM genes that are outside Hallmark IFN-γ (HLA-C, HLA-E, TAP2, ERAP1, ERAP2, CALR, PDIA3, CANX) correlate modestly negatively with CLDN4 (LUAD ρ=-0.115, p=0.009; LUSC ρ=-0.107, p=0.017), i.e. slightly higher in CLDN4-low. In LUSC that association does not survive the keratin partial (ρ=-0.015, p=0.730). The primary 19-gene APM score, which shares its interferon-inducible members with Hallmark IFN-γ, stays null. This sensitivity was defined from set overlap before looking at CLDN4, and it is not in the BH family.

BH q for quartile tests is across the 8 primary contrasts (2 cohorts × NHEJ / IFN-γ / STING / APM). BH q for Spearman is separate for the 8 primary score correlations and for the 12 primary gene correlations (2 cohorts × PRKDC / LIG4 / STING1 / HLA-A / HLA-B / HLA-C). APM-outside-IFN-γ and the MHC-I mean are companions.

## Gene correlations — CLDN4 vs PRKDC, LIG4, STING1, HLA

STING1 is **TMEM173** on this HiSeqV2 freeze. There is no `STING1` or `CGAS` row. MHC-I = mean(HLA-A, HLA-B, HLA-C) on the log2 scale.

| Cohort | Gene | n | Unadj ρ | Unadj p | Unadj q | Partial ρ \| ImmuneScore | Partial p | Partial q |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| TCGA-LUAD | PRKDC | 515 | -0.126 | 0.004 | 0.013 | -0.142 | 0.001 | 0.004 |
| TCGA-LUSC | PRKDC | 501 | -0.229 | 2.09e-07 | 2.51e-06 | -0.237 | 8.30e-08 | 9.79e-07 |
| TCGA-LUAD | LIG4 | 515 | -0.168 | 1.28e-04 | 5.13e-04 | -0.163 | 2.08e-04 | 8.32e-04 |
| TCGA-LUSC | LIG4 | 501 | -0.082 | 0.067 | 0.161 | -0.083 | 0.064 | 0.110 |
| TCGA-LUAD | STING1 | 515 | +0.193 | 1.04e-05 | 6.23e-05 | +0.228 | 1.63e-07 | 9.79e-07 |
| TCGA-LUSC | STING1 | 501 | -0.001 | 0.978 | 0.978 | +0.001 | 0.981 | 0.981 |
| TCGA-LUAD | HLA-A | 515 | +0.059 | 0.180 | 0.359 | +0.109 | 0.014 | 0.027 |
| TCGA-LUSC | HLA-A | 501 | +0.032 | 0.474 | 0.710 | +0.046 | 0.309 | 0.413 |
| TCGA-LUAD | HLA-B † | 515 | +0.035 | 0.424 | 0.710 | +0.114 | 0.010 | 0.024 |
| TCGA-LUSC | HLA-B † | 501 | +0.007 | 0.880 | 0.960 | +0.015 | 0.733 | 0.799 |
| TCGA-LUAD | HLA-C | 515 | +0.023 | 0.596 | 0.732 | +0.076 | 0.084 | 0.126 |
| TCGA-LUSC | HLA-C | 501 | +0.023 | 0.610 | 0.732 | +0.035 | 0.431 | 0.517 |
| TCGA-LUAD | MHC-I | 515 | +0.042 | 0.347 | — | +0.108 | 0.014 | — |
| TCGA-LUSC | MHC-I | 501 | +0.019 | 0.676 | — | +0.032 | 0.476 | — |

† **HLA-B is in Yoshihara Immune141.** Partialling ImmuneScore out of HLA-B is not an independent infiltrate control. HLA-A, HLA-C, PRKDC, LIG4, and TMEM173 are the non-circular gene endpoints. NHEJ ssGSEA **includes** PRKDC and LIG4, so the NHEJ score correlation is not a separate replication of those two genes.

**PRKDC and LIG4 agree with the higher-NHEJ direction in CLDN4-low tumors, where they are significant.** PRKDC ρ=-0.126 in LUAD (q=0.013) and -0.229 in LUSC (q=2.51e-06). The LUSC PRKDC association is the largest gene effect here and remains after keratin (partial ρ=-0.213, p=1.47e-06) and after MKI67 (partial ρ=-0.148, p=9.24e-04). LIG4 is negative in LUAD (ρ=-0.168, q=5.13e-04) and does not clear BH in LUSC (ρ=-0.082, q=0.161).

**STING1 is higher with higher CLDN4 in LUAD, and is null in LUSC.** TMEM173 ρ=+0.193 (q=6.23e-05), ImmuneScore partial +0.228 (q=9.79e-07). That LUAD association remains after keratin (partial ρ=+0.243, p=2.32e-08) and is smaller after MKI67 (partial ρ=+0.153, p=4.85e-04). LUSC TMEM173 ρ=-0.001. This is not a CLDN4-low STING-high result, and it does not replicate in squamous tumors.

**HLA-A/B/C stay unadjusted-null in both histologies**, at the same coefficients as the earlier TCGA RNA table. LUAD partials on ImmuneScore are small and positive for HLA-A and HLA-B (q≈0.02–0.03): higher CLDN4, not lower, with slightly higher HLA after the infiltrate score is removed. HLA-C does not clear BH. LUSC partials stay null. Do not quote a LUAD partial without the unadjusted null, and do not read HLA-B's partial as independent of ImmuneScore.

## Score correlations and the ImmuneScore confound

| Cohort | Score | n | Unadj ρ | p | q | Partial ρ \| ImmuneScore | p | q | ρ(score, ImmuneScore) | Partial ρ \| KRT8/18/19 (p) | Partial ρ \| MKI67 (p) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| TCGA-LUAD | NHEJ | 515 | -0.177 | 5.25e-05 | 4.20e-04 | -0.209 | 1.70e-06 | 1.36e-05 | -0.320 | -0.182 (3.24e-05) | -0.130 (0.003) |
| TCGA-LUSC | NHEJ | 501 | -0.141 | 0.002 | 0.006 | -0.147 | 9.86e-04 | 0.003 | -0.250 | -0.149 (8.46e-04) | -0.058 (0.198) |
| TCGA-LUAD | IFN-γ | 515 | +0.020 | 0.653 | 0.838 | +0.130 | 0.003 | 0.006 | +0.826 | +0.015 (0.729) | +0.034 (0.447) |
| TCGA-LUSC | IFN-γ | 501 | +0.004 | 0.925 | 0.925 | +0.017 | 0.712 | 0.814 | +0.878 | +0.093 (0.037) | -0.025 (0.571) |
| TCGA-LUAD | STING | 515 | +0.090 | 0.041 | 0.110 | +0.162 | 2.36e-04 | 9.44e-04 | +0.601 | +0.103 (0.020) | +0.064 (0.145) |
| TCGA-LUSC | STING | 501 | +0.071 | 0.111 | 0.221 | +0.094 | 0.036 | 0.057 | +0.617 | +0.118 (0.009) | +0.037 (0.410) |
| TCGA-LUAD | APM | 515 | -0.042 | 0.340 | 0.544 | +0.006 | 0.896 | 0.896 | +0.715 | -0.055 (0.216) | -0.026 (0.557) |
| TCGA-LUSC | APM | 501 | -0.015 | 0.733 | 0.838 | -0.018 | 0.688 | 0.814 | +0.736 | +0.059 (0.189) | -0.027 (0.547) |
| TCGA-LUAD | APM outside IFN-γ | 515 | -0.115 | 0.009 | — | -0.095 | 0.032 | — | +0.545 | -0.107 (0.015) | -0.104 (0.018) |
| TCGA-LUSC | APM outside IFN-γ | 501 | -0.107 | 0.017 | — | -0.130 | 0.004 | — | +0.596 | -0.015 (0.730) | -0.107 (0.017) |

IFN-γ and APM track ImmuneScore, as they should: both gene sets are full of infiltrate and interferon genes. CLDN4 itself is only weakly related to ImmuneScore. A partial correlation can therefore move even when the unadjusted association is null. Quote the unadjusted number next to any partial. Keratin (mean of KRT8/KRT18/KRT19) and MKI67 partials are sensitivities, not extra primary tests.

## Q1 versus the other three quartiles

The primary contrast is Q1 vs Q4, matching earlier TCGA quartile cuts. Q1 versus everyone else is the same direction test with a larger control arm.

| Cohort | Score | n Q1 | n rest | Δ median (Q1−rest) | MW p |
|---|---|---:|---:|---:|---:|
| TCGA-LUAD | NHEJ | 129 | 386 | +0.010 | 0.006 |
| TCGA-LUSC | NHEJ | 126 | 375 | +0.005 | 0.032 |
| TCGA-LUAD | IFN-γ | 129 | 386 | +0.017 | 0.262 |
| TCGA-LUSC | IFN-γ | 126 | 375 | +0.007 | 0.455 |
| TCGA-LUAD | STING | 129 | 386 | -0.016 | 0.141 |
| TCGA-LUSC | STING | 126 | 375 | -0.006 | 0.407 |
| TCGA-LUAD | APM | 129 | 386 | +0.015 | 0.108 |
| TCGA-LUSC | APM | 126 | 375 | +0.014 | 0.335 |

## What the gene sets are, and what they share

- **NHEJ** (13 genes), KEGG hsa03450: RAD50, DNTT, FEN1, XRCC6, POLL, POLM, LIG4, MRE11A, PRKDC, DCLRE1C, XRCC4, XRCC5, NHEJ1. `MRE11` is `MRE11A` on HiSeqV2.
- **IFN-γ** (200 genes), Hallmark interferon gamma response. Four symbols were mapped onto this freeze: WARS1→WARS, MARCHF1→MARCH1, HELZ2→PRIC285, CMTR1→FTSJD2. Genes absent after that map are listed in `tables/genesets.tsv` only if dropped; the scored set has 200 genes.
- **STING** (11 genes), Reactome R-HSA-1834941 after removing KEGG-NHEJ members (PRKDC, XRCC5, XRCC6, MRE11) and CGAS (not on HiSeqV2): DDX41, DTX4, IFI16, IRF3, NLRC3, NLRP4, STAT6, TMEM173, TBK1, TREX1, TRIM21. `STING1` is `TMEM173`. The set mixes positive components (TMEM173, TBK1, IRF3, IFI16) with negative regulators (TREX1, NLRC3, NLRP4, DTX4), so it is a pathway membership score, not a pure induced-ISG score. TRIM21 is also in Hallmark IFN-γ.
- **APM** (19 genes): HLA-A, HLA-B, HLA-C, HLA-E, B2M, TAP1, TAP2, TAPBP, PSMB8, PSMB9, PSMB10, PSME1, PSME2, ERAP1, ERAP2, NLRC5, CALR, PDIA3, CANX.
- **APM outside IFN-γ** (8 genes), sensitivity: HLA-C, HLA-E, TAP2, ERAP1, ERAP2, CALR, PDIA3, CANX. The other APM genes sit inside the Hallmark IFN-γ set, so the primary APM score is not an independent IFN test.

## What this does not say

- Bulk primary-tumor RNA. It is not a malignant-cell NHEJ or STING measurement, and it is not an ICI-response result.
- A null quartile test is not evidence that NHEJ or cGAS–STING is irrelevant in CLDN4-high cells. It is evidence that this bulk contrast does not show the pre-specified shift.
- ImmuneScore is RNA-derived. Partialling it is not an ABSOLUTE purity adjustment. Hallmark IFN-γ and the APM set overlap Immune141-style genes, so those partials are partly circular. HLA-B is the gene-level case that is explicitly inside Immune141.
- No STAR-TPM rerun, no protein, no CPTAC, no single-cell split.

## Methods

- **Matrix.** UCSC Xena `TCGA.{LUAD,LUSC}.sampleMap/HiSeqV2`. Primary tumors only (sample type 01). Replicate aliquots averaged to the 15-character barcode.
- **ESTIMATE.** Official RNAseqV2 `Immune_score`, not recomputed.
- **ssGSEA.** gseapy single-sample GSEA. Sample normalization `rank`, weight α=0.25, `permutation_num=0` (the NES column is gseapy's normalized enrichment score, not a permutation p-value). Minimum set size 5. A synthetic check in this environment confirms that a sample with the NHEJ genes pinned to the top of the ranking receives a positive NES. Absolute NES is not a significance call, and NES cannot be compared across gene sets of different sizes. Every test is a within-cohort contrast.
- **CLDN4-low / high.** Quartiles of CLDN4 within each complete-case cohort. `rank(method='first')` after sorting barcodes, then `qcut` into four equal groups. Q1 is the lowest CLDN4.
- **Quartile test.** Two-sided Mann–Whitney Q1 vs Q4. Effect sizes: median(Q1)−median(Q4) and Cliff's delta (positive means Q1 values are larger). Q1 vs Q2–Q4 is secondary and is not in the BH family.
- **Correlation.** Spearman. Partial = first-order partial Spearman given ImmuneScore (df = n−3). Fisher z intervals use 1/(n−3) unadjusted and 1/(n−4) partial. Keratin mean and MKI67 partials use the same formula and are labeled sensitivities.
- **FDR.** Benjamini–Hochberg within the families named above. Families are not pooled with each other.

## Specification sweep

An exploratory grid over cutoffs, gene sets, ssGSEA / GSVA / AUCell / z-mean, covariates, and cohort strata is in `SWEEP.md`. It does not replace the pre-specified tests above. No quadruple in that grid is significant in the thesis direction for NHEJ together with IFN, STING, and APM.

## Files

- `tables/correlations.tsv` — Spearman, partials, CI, q
- `tables/quartiles.tsv` — Q1 vs Q4 and Q1 vs rest
- `tables/counts.tsv` — honest-n filters
- `tables/samples.tsv` — per-tumor CLDN4, quartiles, scores, genes, ImmuneScore
- `tables/genesets.tsv` — genes that entered each ssGSEA set
- `tables/provenance.json` — URLs, sha256, gseapy version
- `figures/q1q4_ssgsea.png` — quartile NES
- `figures/forest_partial.png` — unadjusted vs ImmuneScore-partial ρ
- `figures/scatter_genes.png` — CLDN4 vs PRKDC, LIG4, STING1, MHC-I
- `figures/scatter_scores.png` — CLDN4 vs the four ssGSEA scores
- Reproduce: `python3 methods/tcga_cldn4_nhej_sting/analyze.py`

