# REWORK A1 wave 2 — TCGA-LUAD double residual (ABSOLUTE + epithelium) and ImmuneScore residual

**Self-contained. Public data only. No fabricated numbers.**

**Cohort:** TCGA-LUAD primary tumors (`-01`) only. Not LUSC. Not OncoSG. Not the pooled n ≈ 1031.

**Question:** Prior LUAD-only A1 was weak after ABSOLUTE (|ρ| ≤ 0.22 markers, PR #81; |ρ| ≤ 0.13 xCell/MCP/GEP18/ESTIMATE, PR #169). OncoSG and LUSC *did* support a negative TACSTD2–immune association after purity. Can a second residual — epithelium or ESTIMATE ImmuneScore — recover the user's implied **strong inverse** on LUAD?

## Verdict: STILL_WEAK

- One-sided double residual (requested): most negative ρ = -0.072; |ρ| stays inside the prior LUAD marker bound of 0.22 and far from OncoSG −0.30.
- Two-sided partial | ABS+KRT/EPCAM: most negative ρ = -0.082; |ρ| stays inside the prior LUAD marker bound of 0.22 and far from OncoSG −0.30.
- ImmuneScore residual, one-sided: most negative ρ = -0.140; |ρ| stays inside the prior LUAD marker bound of 0.22 and far from OncoSG −0.30.
- ImmuneScore residual, two-sided: most negative ρ = -0.211; |ρ| stays inside the prior LUAD marker bound of 0.22 and far from OncoSG −0.30.

ABSOLUTE-only two-sided partials **reproduce** the histology / deconv reworks: CD8 −0.104, CYT −0.116, GEP18 z-mean −0.013, xCell CD8 −0.126 (n=502). TACSTD2 ⟂ ABSOLUTE (ρ = 0.007, p = 0.869).

The epithelial double residual **weakens** those already-small inverses (TACSTD2 vs KRT/EPCAM ρ = 0.364). GEP18 even flips to a small positive. Epithelium was not masking a strong LUAD inverse.

Scan-wide strongest primary ρ is `ABS_IMMUNE` two-sided CYT = −0.217 (n=502). That is still `WEAK_NEGATIVE` under the locked −0.22 / −0.30 rules. It is not hidden; it is not OncoSG.

**Implied strong inverse, locked before results.** OncoSG A1 after published purity (PR #139, n=169) is the user's working example of a supported inverse: CD8A partial ρ = −0.309, GEP18 = −0.349. Thresholds:

| Call | Rule |
|---|---|
| `STRONG_INVERSE` | ρ ≤ **−0.30** (OncoSG CD8/GEP) |
| `REACHES_LUSC_BOUND` | −0.30 < ρ ≤ **−0.22** (LUSC CD8 / prior LUAD marker bound) |
| `WEAK_NEGATIVE` | −0.22 < ρ < 0 |
| `NOT_INVERSE` | ρ ≥ 0 |

n ≈ 1031 in the original claim is the **pooled TCGA LUAD+LUSC** ESTIMATE-complete set (PR #89), not LUAD-only. This file stays LUAD-only. Complete-case n is in every table cell (~502 with ABSOLUTE).

## Why one-sided and two-sided are both shown

The request is: residualize **TACSTD2** on (1) ABSOLUTE and (2) epithelial fraction / KRT/EPCAM, **then** correlate with CD8 / GEP18 / CYT / xCell CD8.

That is a **one-sided** residual. Immune scores still carry epithelial / purity composition. Because KRT/EPCAM is expected to anti-correlate with CD8 in bulk RNA, one-sided residualization can *inflate* a negative ρ. The **two-sided** rank residual (standard multi-covariate partial Spearman) is the fair test of “TACSTD2 vs immune beyond tumor content and epithelium.” Both are reported. The overall tag uses the two-sided result as the claim-level answer and states the one-sided number separately.

Residualizing TACSTD2 on ESTIMATE ImmuneScore and then correlating with CD8/GEP18/CYT/xCell CD8 is valid (different features). Correlating that residual with ImmuneScore itself would be circular and is not used as a primary readout. ImmuneScore is strongly anti-correlated with ABSOLUTE (ρ = -0.616) by construction of leukocyte content.

TIMER2 xCell has **no Epithelial-cells column** in this freeze (keratinocyte column: also absent). Epithelial fraction is therefore the KRT/EPCAM gene score and ESTIMATE TumorPurity, not an xCell epithelial fraction. That absence is reported, not filled in.

## Analysis set

| Filter | n |
|---|---:|
| Xena HiSeqV2 primary tumors (`-01`) | 515 |
| + ABSOLUTE purity | 502 |
| + KRT/EPCAM score (same matrix) | 502 |
| + official ESTIMATE ImmuneScore | 515 |
| + TIMER2 xCell CD8 | 515 |

Primary tumors only. One row per 15-character barcode. Replicate aliquots averaged.

## Definitions (locked)

### TACSTD2
UCSC Xena `TCGA.LUAD.sampleMap/HiSeqV2` **log2(RSEM normalized_count + 1)**. Gene symbol `TACSTD2`. Not protein. Not ADC. Not ICI.

### Samples
TCGA-LUAD only. Sample type `01`. No LUSC, no OncoSG, no normals.

### Covariates
1. **ABSOLUTE purity** — PanCanAtlas GDC `4f277128-f793-4354-a13d-30cc7fe9f6b5`, column `purity`.
2. **KRT/EPCAM score (primary epithelium)** — mean log2(RSEM+1) of `EPCAM, KRT7, KRT8, KRT18, KRT19`. Present: EPCAM, KRT7, KRT8, KRT18, KRT19. Missing (not imputed): none. **TACSTD2 is not in the score.**
3. **ESTIMATE TumorPurity (sensitivity epithelium)** — Yoshihara 2013 transform `cos(0.6049872018 + 0.0001467884 × ESTIMATEScore)` from the official MD Anderson LUAD RNAseqV2 table. This is RNA tumor fraction, not ABSOLUTE.
4. **ESTIMATE ImmuneScore** — official `Immune_score` from the same table. Used only as the requested TACSTD2 residual covariate, not as a primary immune *outcome*.
5. **xCell epithelial** — used only if a TIMER2 `*_XCELL` epithelial column exists; otherwise skipped (see above).

Sensitivity epithelium: `EPCAM, KRT7, KRT8, KRT18, KRT19, CDH1` mean; and mean of EPCAM + every `KRT[0-9]+` gene in HiSeqV2 (56 genes).

### Immune outcomes (primary)
| Name | Definition |
|---|---|
| CD8 | `CD8A` log2(RSEM+1) on the same Xena matrix |
| CYT | Rooney *Cell* 2015 cytolytic score = mean(`GZMA`, `PRF1`) on log2 data |
| GEP18 | Ayers *JCI* 2017 18 genes, **unweighted within-LUAD z-mean** (histology-rework definition, PR #107) |
| xCell CD8 | TIMER2.0 immunedeconv column `T cell CD8+_XCELL` |

Sensitivity outcomes (not in the verdict): `CD8B`, GEP18 ssGSEA (Barbie/GSVA tau=0.25, ranks 1..10000), GEP18 log-mean, official ESTIMATE ImmuneScore as a *Y* (only for ABS/KRT models).

### Estimators
- **onesided_rank (user-requested primary):** OLS of rank(TACSTD2) on rank(covariates); Pearson of that residual vs rank(immune). df = n − 2 − k.
- **onesided_value:** OLS of TACSTD2 on raw covariates; Spearman of residual vs raw immune. Same df.
- **twosided_rank:** both sides residualized on ranked covariates (partial Spearman). Same df.
- Unadjusted Spearman is the no-covariate baseline.
- 95% Fisher-z CI uses variance `1/(n − k − 3)`.
- BH-FDR across the 4 primary immune features, within each model × estimator.

## Covariate context (unadjusted Spearman)

| Pair | ρ | p | n |
|---|---:|---:|---:|
| TACSTD2 vs ABSOLUTE | 0.007 | 0.869 | 502 |
| TACSTD2 vs KRT/EPCAM | 0.364 | 1.47e-17 | 515 |
| TACSTD2 vs ESTIMATE ImmuneScore | 0.036 | 0.41 | 515 |
| ABSOLUTE vs KRT/EPCAM | 0.023 | 0.61 | 502 |
| ABSOLUTE vs ImmuneScore | -0.616 | 7.71e-54 | 502 |

If TACSTD2 is nearly orthogonal to ABSOLUTE (as in PR #81 / #169), purity-only adjustment cannot create a strong inverse. A second residual can only help if TACSTD2 shares variance with epithelium (or ImmuneScore) that was *masking* an inverse with CD8/GEP/CYT.


## Pointer

Numbers live in `REPORT.md` and `correlations.tsv`. This file is the locked definition list.
