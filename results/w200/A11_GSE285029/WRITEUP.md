# GSE285029: TACSTD2 / CLDN4 vs PD-L1, IFN/MHC-I, IL-6/STAT3, CD8/GEP

> Honest public-matrix report. Numbers are computed from the GEO author file.
> A11 (galectin / nectin / TGF-β / CD47 with TROP2) and C (SKB264 raises PD-L1,
> so ADC+ICI is required) are taken as given. This cohort is extra
> **baseline co-expression** evidence only.
> GEO has no RECIST, histology, PD-L1 IHC, TMB, or purity.

## TL;DR

**TROP2-high: no. CLDN4-high: IFN/MHC-I yes, CD274-Q4 no.**

The two anchors co-express (ρ = +0.526, p = 5.01e-18, n = 234) and then split.

- **TACSTD2 vs CD274** ρ = +0.157 (p = 0.016, BH q = 0.112, n = 234).
  vs IFN-compact ρ = +0.097 (p = 0.139).
  vs GEP18 ρ = +0.008 (p = 0.905).
  vs CD8A ρ = -0.084 (p = 0.199).
  TACSTD2 Q4 ∩ CD274 Q4 = 25.4% (expected 25%; OR = 1.01, p = 1.000).
  TACSTD2 Q4 is **depleted** for GEP18 Q4 (13.6%, OR = 0.38, p = 0.023)
  and CD8A Q4 (11.9%, OR = 0.32, p = 0.006).
- **CLDN4 vs IFN-compact** ρ = +0.204 (p = 0.002, BH q = 0.006).
  vs MHC-I ρ = +0.263 (p = 4.62e-05).
  vs Hallmark IFN-γ ρ = +0.301 (p = 2.65e-06).
  vs CD274 ρ = +0.167 (p = 0.011).
  CLDN4 Q4 ∩ IFN-compact Q4 = 42.4% (OR = 3.05, p = 8.54e-04).
  CLDN4 Q4 ∩ CD274 Q4 = 32.2% (OR = 1.60, p = 0.167).
  CD8A and IL6 stay null.

TROP2-high sit-on-program: **NO — Q4 overlap with high CD274 / high IFN is consistent with chance.**
CLDN4-high sit-on-program: **PARTIAL — Q4 is enriched for a high IFN program, not for high CD274.**

Extra combination-rationale only (A11 and C stay given): TROP2-high bulk tumors here are not already PD-L1-high or IFN-high. CLDN4-high tumors sit on an IFN / MHC-I transcriptional program without being CD274-Q4-high or CD8-high. The ADC-target (TROP2) side of this matrix does not supply a baseline “already hot / already PD-L1-high” argument.

## Cohort

| Item | Value |
|---|---|
| GEO | [GSE285029](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285029) |
| Paper | Koh et al., *J Immunother Cancer* 2025; PMID [40050048](https://pubmed.ncbi.nlm.nih.gov/40050048/) |
| Design | Pre-ICI NSCLC tumor WTS (PD-1 or PD-L1 blockade), Illumina HiSeq 2500 |
| Public matrix | `GSE285029_WTS_expr_count_235_032820.txt.gz` |
| Samples | **n = 234** (Case1–Case234) |
| Genes | 20069 symbols |
| GEO phenotype | tissue = Lung; cell type = cancer; genotype = wt. No response / histology / purity |
| Author n in paper | 234 (ICI–RNA-seq cohort) |

The filename says “count” and 235 columns (gene id + 234 samples). Values are
continuous with **13.6% negatives** (global min
-6339.0, max 252628.2, median
column sum 590632). That is not raw integer counts
and not TPM (TPM column sums would be ~1e6). It is the author-processed WTS
matrix released on GEO. Primary transform: `log2(pmax(x,0)+1)`. Sensitivity:
Spearman on the raw author values (negatives kept). Rank correlations are the
claim; the clip is only to put the matrix on a standard log scale.

## Pre-specified design

- Anchors: `TACSTD2` (TROP2), `CLDN4` (junction partner, not an immune gene).
- Primary partners: `CD274`; IFN-compact (`IFNG STAT1 IRF1 CXCL9 CXCL10 CXCL11 IDO1 GBP1`);
  MHC-I cassette (`HLA-A/B/C B2M TAP1 TAP2 NLRC5 PSMB8 PSMB9 TAPBP`);
  `IL6`; IL6/STAT3-compact (`IL6 IL6R IL6ST JAK1 JAK2 STAT3 SOCS3`);
  `CD8A`; Ayers 2017 GEP18 (unweighted z-mean; not NanoString TIS weights).
- Sensitivity: MSigDB Hallmark IFN-γ / IFN-α / IL6-JAK-STAT3 (2024.1.Hs symbols,
  bundled), A1 8-gene immune, CYT (`GZMA+PRF1`), single-gene `IFNG` and `STAT3`.
- Score = unweighted mean of per-gene z-scores on the log2(clip0+1) matrix.
- Primary statistic: Spearman ρ, two-sided, Fisher-z 95% CI, n = complete pairs.
- TACSTD2-high / CLDN4-high: Q4 vs Q1 Mann-Whitney U. Median split is sensitivity.
- “Sit on a high program”: Q4×Q4 overlap vs 25% independence, Fisher exact.
- FDR: BH within each anchor across the 7 primary partners only. Sensitivity p-values are descriptive.
- Epithelial z-mean (`EPCAM KRT8 KRT18 KRT19`) is a purity-like **sensitivity**
  covariate, not published ABSOLUTE / ESTIMATE purity.
- No filter was tuned to produce a positive PD-L1 / IFN class effect.

## Primary result — TACSTD2

| feature | role | ρ | 95% CI | p | BH q | n | label | epi-partial ρ | epi-partial p |
|---|---|---:|---|---:|---:|---:|---|---:|---:|
| CD274 | primary | +0.157 | +0.030 to +0.280 | 0.016 | 0.112 | 234 | WEAK_POSITIVE | +0.102 | 0.119 |
| IFN_compact | primary | +0.097 | -0.032 to +0.223 | 0.139 | 0.447 | 234 | NULL | +0.001 | 0.989 |
| MHC1 | primary | +0.075 | -0.054 to +0.201 | 0.256 | 0.447 | 234 | NULL | -0.075 | 0.253 |
| IL6 | primary | -0.046 | -0.173 to +0.083 | 0.484 | 0.637 | 234 | NULL | -0.018 | 0.786 |
| IL6_STAT3_compact | primary | +0.040 | -0.089 to +0.167 | 0.546 | 0.637 | 234 | NULL | -0.042 | 0.526 |
| CD8A | primary | -0.084 | -0.210 to +0.044 | 0.199 | 0.447 | 234 | NULL | -0.090 | 0.170 |
| GEP18 | primary | +0.008 | -0.121 to +0.136 | 0.905 | 0.905 | 234 | NULL | -0.052 | 0.427 |
| HALLMARK_IFNG | sensitivity | +0.129 | +0.001 to +0.253 | 0.049 | — | 234 | WEAK_POSITIVE | -0.023 | 0.726 |
| HALLMARK_IFNA | sensitivity | +0.151 | +0.023 to +0.274 | 0.021 | — | 234 | WEAK_POSITIVE | -0.035 | 0.599 |
| HALLMARK_IL6_JAK_STAT3 | sensitivity | +0.120 | -0.008 to +0.245 | 0.066 | — | 234 | NULL | -0.002 | 0.973 |
| immune8 | sensitivity | -0.043 | -0.170 to +0.086 | 0.516 | — | 234 | NULL | -0.056 | 0.392 |
| CYT | sensitivity | -0.045 | -0.172 to +0.084 | 0.495 | — | 234 | NULL | -0.049 | 0.460 |
| STAT3 | sensitivity | +0.180 | +0.053 to +0.301 | 0.006 | — | 234 | WEAK_POSITIVE | +0.031 | 0.634 |
| IFNG | sensitivity | +0.043 | -0.085 to +0.171 | 0.509 | — | 234 | NULL | +0.033 | 0.611 |

TACSTD2 vs CLDN4 (same matrix): ρ = +0.526 (p = 5.01e-18, n = 234).

## Primary result — CLDN4

| feature | role | ρ | 95% CI | p | BH q | n | label | epi-partial ρ | epi-partial p |
|---|---|---:|---|---:|---:|---:|---|---:|---:|
| CD274 | primary | +0.167 | +0.040 to +0.289 | 0.011 | 0.015 | 234 | WEAK_POSITIVE | +0.101 | 0.124 |
| IFN_compact | primary | +0.204 | +0.077 to +0.323 | 0.002 | 0.006 | 234 | ASSOCIATED_POSITIVE | +0.102 | 0.120 |
| MHC1 | primary | +0.263 | +0.139 to +0.379 | 4.62e-05 | 3.23e-04 | 234 | ASSOCIATED_POSITIVE | +0.112 | 0.088 |
| IL6 | primary | -0.035 | -0.163 to +0.093 | 0.591 | 0.591 | 234 | NULL | +0.005 | 0.945 |
| IL6_STAT3_compact | primary | +0.191 | +0.064 to +0.312 | 0.003 | 0.006 | 234 | WEAK_POSITIVE | +0.119 | 0.071 |
| CD8A | primary | +0.067 | -0.061 to +0.194 | 0.305 | 0.356 | 234 | NULL | +0.093 | 0.159 |
| GEP18 | primary | +0.190 | +0.063 to +0.310 | 0.004 | 0.006 | 234 | WEAK_POSITIVE | +0.153 | 0.020 |
| HALLMARK_IFNG | sensitivity | +0.301 | +0.180 to +0.414 | 2.65e-06 | — | 234 | ASSOCIATED_POSITIVE | +0.147 | 0.025 |
| HALLMARK_IFNA | sensitivity | +0.342 | +0.224 to +0.451 | 7.97e-08 | — | 234 | ASSOCIATED_POSITIVE | +0.152 | 0.020 |
| HALLMARK_IL6_JAK_STAT3 | sensitivity | +0.267 | +0.144 to +0.382 | 3.54e-05 | — | 234 | ASSOCIATED_POSITIVE | +0.142 | 0.030 |
| immune8 | sensitivity | +0.080 | -0.048 to +0.207 | 0.220 | — | 234 | NULL | +0.090 | 0.172 |
| CYT | sensitivity | +0.095 | -0.034 to +0.220 | 0.148 | — | 234 | NULL | +0.122 | 0.062 |
| STAT3 | sensitivity | +0.307 | +0.186 to +0.418 | 1.73e-06 | — | 234 | ASSOCIATED_POSITIVE | +0.144 | 0.028 |
| IFNG | sensitivity | +0.112 | -0.016 to +0.237 | 0.087 | — | 234 | NULL | +0.118 | 0.071 |

## TACSTD2-high vs TACSTD2-low (Q4 vs Q1)

| feature | median Q4 | median Q1 | Δmedian | MWU p | n Q4 | n Q1 |
|---|---:|---:|---:|---:|---:|---:|
| CD274 | +2.464 | +1.947 | +0.517 | 0.017 | 59 | 59 |
| IFN_compact | -0.224 | -0.416 | +0.192 | 0.202 | 59 | 59 |
| MHC1 | +0.003 | -0.058 | +0.061 | 0.327 | 59 | 59 |
| IL6 | +1.384 | +1.638 | -0.254 | 0.477 | 59 | 59 |
| IL6_STAT3_compact | +0.128 | +0.062 | +0.066 | 0.813 | 59 | 59 |
| CD8A | +1.811 | +2.084 | -0.273 | 0.276 | 59 | 59 |
| GEP18 | -0.149 | -0.167 | +0.019 | 0.953 | 59 | 59 |
| HALLMARK_IFNG | +0.026 | -0.021 | +0.046 | 0.177 | 59 | 59 |
| HALLMARK_IL6_JAK_STAT3 | +0.098 | -0.001 | +0.099 | 0.152 | 59 | 59 |

## CLDN4-high vs CLDN4-low (Q4 vs Q1)

| feature | median Q4 | median Q1 | Δmedian | MWU p | n Q4 | n Q1 |
|---|---:|---:|---:|---:|---:|---:|
| CD274 | +2.838 | +1.947 | +0.891 | 0.003 | 59 | 59 |
| IFN_compact | +0.256 | -0.539 | +0.796 | 8.63e-04 | 59 | 59 |
| MHC1 | +0.425 | -0.147 | +0.572 | 1.22e-04 | 59 | 59 |
| IL6 | +1.658 | +1.741 | -0.083 | 0.953 | 59 | 59 |
| IL6_STAT3_compact | +0.329 | +0.062 | +0.267 | 0.003 | 59 | 59 |
| CD8A | +2.538 | +2.184 | +0.354 | 0.100 | 59 | 59 |
| GEP18 | +0.236 | -0.182 | +0.418 | 0.002 | 59 | 59 |
| HALLMARK_IFNG | +0.320 | -0.141 | +0.460 | 4.19e-06 | 59 | 59 |
| HALLMARK_IL6_JAK_STAT3 | +0.222 | -0.001 | +0.223 | 5.42e-05 | 59 | 59 |

## Do they sit on a high PD-L1 or high IFN program?

Independence expectation for two Q4 calls is 25%. Enrichment is the extra
combination-rationale test.

### TACSTD2 Q4

| feature | both Q4 | anchor Q4 | frac in feature Q4 | expected if independent | OR | Fisher p |
|---|---:|---:|---:|---:|---:|---:|
| CD274 | 15 | 59 | 0.254 | 0.250 | 1.01 | 1.000 |
| IFN_compact | 11 | 59 | 0.186 | 0.250 | 0.61 | 0.225 |
| HALLMARK_IFNG | 12 | 59 | 0.203 | 0.250 | 0.70 | 0.387 |
| GEP18 | 8 | 59 | 0.136 | 0.250 | 0.38 | 0.023 |
| IL6 | 15 | 59 | 0.254 | 0.250 | 1.01 | 1.000 |
| CD8A | 7 | 59 | 0.119 | 0.250 | 0.32 | 0.006 |

### CLDN4 Q4

| feature | both Q4 | anchor Q4 | frac in feature Q4 | expected if independent | OR | Fisher p |
|---|---:|---:|---:|---:|---:|---:|
| CD274 | 19 | 59 | 0.322 | 0.250 | 1.60 | 0.167 |
| IFN_compact | 25 | 59 | 0.424 | 0.250 | 3.05 | 8.54e-04 |
| HALLMARK_IFNG | 25 | 59 | 0.424 | 0.250 | 3.05 | 8.54e-04 |
| GEP18 | 21 | 59 | 0.356 | 0.250 | 1.99 | 0.039 |
| IL6 | 14 | 59 | 0.237 | 0.250 | 0.90 | 0.863 |
| CD8A | 17 | 59 | 0.288 | 0.250 | 1.28 | 0.490 |

**TROP2-high:** NO — Q4 overlap with high CD274 / high IFN is consistent with chance.

**CLDN4-high:** PARTIAL — Q4 is enriched for a high IFN program, not for high CD274.

## Honest reading

- **CD274:** TACSTD2 ρ = +0.157 (WEAK_POSITIVE; BH q = 0.112).
  CLDN4 ρ = +0.167 (WEAK_POSITIVE; BH q = 0.015).
  Both continuous associations are small. Neither Q4 is enriched for CD274 Q4
  (TROP2 25.4% p = 1.00; CLDN4 32.2% p = 0.17). Epithelial partial drops both
  CD274 ρ values to ~0.10 (p ≈ 0.12). This is not a high-PD-L1 class for either anchor.
- **IFN / MHC-I:** TACSTD2 vs IFN-compact ρ = +0.097 (NULL); vs MHC-I ρ = +0.075 (NULL).
  Hallmark IFN-γ / IFN-α are only WEAK_POSITIVE for TACSTD2 and go to null after the
  epithelial residual. CLDN4 vs IFN-compact ρ = +0.204; vs MHC-I ρ = +0.263;
  vs Hallmark IFN-γ ρ = +0.301. CLDN4 Q4 is enriched for IFN-compact Q4
  (42.4%, OR = 3.05, p = 8.5e-4). That is a real IFN / antigen-presentation neighborhood
  for CLDN4-high, not for TROP2-high. Epithelial partial attenuates CLDN4–IFN-compact
  to ρ = +0.102 (p = 0.120); Hallmark IFN-γ
  and GEP18 remain weakly positive after the residual.
- **CD8 / GEP:** TACSTD2 vs CD8A ρ = -0.084; vs GEP18 ρ = +0.008.
  TACSTD2 Q4 is depleted for CD8A Q4 and GEP18 Q4. CLDN4 vs CD8A ρ = +0.067 (NULL);
  vs GEP18 ρ = +0.190 (WEAK_POSITIVE). CLDN4-high is IFN/MHC-I-high without
  being CD8-high. GEP18 here is pulled by the IFN/MHC genes in the 18-gene set, not by CD8A.
- **IL-6 / STAT3:** TACSTD2 vs IL6 ρ = -0.046; vs IL6/STAT3-compact ρ = +0.040.
  CLDN4 vs IL6 ρ = -0.035 (NULL); vs IL6/STAT3-compact ρ = +0.191.
  The source paper’s PD-L1→IL-6 axis is about *CD274-high* tumors. IL6 itself is null
  for both anchors. The compact / Hallmark STAT3 scores that track CLDN4 are not IL6 mRNA.
- **Epithelial partial:** TACSTD2 and CLDN4 are epithelial genes. Residualizing on
  `EPCAM/KRT8/18/19` is a sensitivity check, not published ABSOLUTE / ESTIMATE purity.
  The TROP2–CD274 and CLDN4–IFN-compact primary ρ values both lose p < 0.05 after this
  residual. Hallmark IFN-γ / GEP18 for CLDN4 do not.
- **Combination-rationale (extra only):** A11 and C stay as given. This public slice
  does **not** add “TROP2-high already transcribes high PD-L1 or high IFN.”
  It adds a CLDN4-high IFN/MHC-I neighborhood and a TROP2-high CD8/GEP depletion.
  Any ADC+ICI PD-L1 argument for the TROP2-high class in this cohort is not a
  baseline co-expression argument.

## What this does **not** show

- It does not show ICI response, ADC response, or SKB264 on-treatment PD-L1 induction.
  GEO has no RECIST labels; those analyses are not performed.
- It does not show protein PD-L1 (IHC TPS) or serum IL-6. `CD274` mRNA ≠ TPS.
- It does not show tumor-cell-intrinsic PD-L1 signaling. Bulk `CD274` mixes tumor and immune cells.
- Hallmark IL6-JAK-STAT3 is a mixed transcriptional set, not phospho-STAT3.
- GEP18 is an unweighted z-mean, not the commercial TIS assay.
- n = 234 is one institution’s ICI-era NSCLC WTS. It is not TCGA and not a TROP2-ADC trial.

## Gene-set coverage

| score | role | present / defined | missing |
|---|---|---|---|
| TACSTD2 | anchor | 1/1 | — |
| CLDN4 | anchor | 1/1 | — |
| CD274 | primary | 1/1 | — |
| IL6 | primary | 1/1 | — |
| CD8A | primary | 1/1 | — |
| STAT3 | sensitivity | 1/1 | — |
| IFNG | sensitivity | 1/1 | — |
| IFN_compact | primary | 8/8 | — |
| MHC1 | primary | 10/10 | — |
| IL6_STAT3_compact | primary | 7/7 | — |
| GEP18 | primary | 18/18 | — |
| immune8 | sensitivity | 8/8 | — |
| CYT | sensitivity | 2/2 | — |
| epithelial | covariate | 4/4 | — |
| HALLMARK_IFNG | sensitivity | 196/200 | MARCHF1,RIGI,TMT1B,WARS1 |
| HALLMARK_IFNA | sensitivity | 95/97 | TENT5A,WARS1 |
| HALLMARK_IL6_JAK_STAT3 | sensitivity | 87/87 | — |

## Files

- `correlations.csv` — marginal and epithelial-partial Spearman (log2 and raw)
- `q4_vs_q1.csv` — Q4 vs Q1 MWU
- `median_split.csv` — median-split MWU
- `q4_overlap.csv` — Q4×Q4 Fisher exact
- `gene_coverage.csv` — genes present / missing per score
- `sample_scores.csv` — per-sample anchors and scores (log2 clip)
- `matrix_qc.json` — n, negatives, column sums
- `summary.json` — machine-readable verdict
- `fig1_rho_heatmap.png` / `fig2_scatter.png` / `fig3_q4_boxplots.png` / `fig4_q4_overlap.png`
- Code: `scripts/w200/A11_GSE285029/`
