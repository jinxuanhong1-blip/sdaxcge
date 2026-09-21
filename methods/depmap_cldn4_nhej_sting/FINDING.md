# DepMap: CLDN4 vs NHEJ and cGAS–STING

The match is in RNA, and it replicates across three DepMap releases. A 24Q4 search over lineages, Spearman and Pearson, Chronos and RNA (1,824 tests) picks one winner. Search-wide Benjamini–Hochberg on that grid still leaves it.

**Winner.** CLDN4 RNA vs **PAXX** RNA (NHEJ), all cell lines, Spearman.

| release | n | ρ | p | 95% CI |
|---|---:|---:|---:|---|
| **24Q4** | **1,673** | **+0.375** | **6.3×10⁻⁵⁷** | **[+0.33, +0.42]** |
| 24Q2 | 1,517 | +0.329 | 1.3×10⁻³⁹ | |
| 22Q4 | 1,408 | +0.292 | 4.5×10⁻²⁹ | |

Search-wide *q* on the 1,824-test 24Q4 grid: **1.2×10⁻⁵³**. Lineage-median residual on 24Q4: ρ=+0.201, p=1.4×10⁻¹⁶ (n=1,651). Partial Spearman given EPCAM: ρ=+0.178, p=2.1×10⁻¹³. Given MKI67: ρ=+0.360, p=3.5×10⁻⁵². Within-lineage Stouffer (20 lineages, n≥25): weighted mean ρ=+0.251, p=8.7×10⁻²⁴.

**Strongest cGAS–STING gene on the same grid.** CLDN4 RNA vs **IRF3** RNA.

| release | all-line ρ (p) | lineage-residual ρ (p) | within-lineage Stouffer p |
|---|---|---|---|
| **24Q4** | **+0.325 (1.4×10⁻⁴²)**, n=1,673, CI [+0.28, +0.37] | **+0.285 (2.5×10⁻³²)** | **3.3×10⁻⁴⁵** |
| 24Q2 | +0.278 (2.4×10⁻²⁸), n=1,517 | +0.254 (2.0×10⁻²³) | 1.7×10⁻³¹ |
| 22Q4 | +0.259 (4.5×10⁻²³), n=1,408 | +0.241 (6.5×10⁻²⁰) | 3.5×10⁻²⁵ |

24Q4 partial given EPCAM: ρ=+0.188, p=1.1×10⁻¹⁴. Given MKI67: ρ=+0.305, p=3.3×10⁻³⁷. Of 20 lineages with n≥25, **19** have positive CLDN4–IRF3 ρ and **12** have BH *q*<0.05 across those lineages. The largest single lineage is skin / melanoma: skin ρ=+0.597 (n=120, p=6.0×10⁻¹³, CI [+0.47, +0.69]); melanoma ρ=+0.618 (n=110, p=6.4×10⁻¹³). 24Q2 skin repeats (ρ=+0.521, n=95, p=6.1×10⁻⁸).

**STING1 is the lung filter.** Pan-cancer STING1 RNA is flat (24Q4 ρ=−0.040, p=0.099). Across 20 lineages, lung is the only one with STING1 BH *q*<0.05: ρ=**+0.309**, n=214, p=4.1×10⁻⁶, *q*=8.3×10⁻⁵, CI [+0.18, +0.43]. Partial given EPCAM inside lung: ρ=+0.464, p=9.4×10⁻¹³. That partial repeats on 24Q2 (ρ=+0.427, n=204, p=2.1×10⁻¹⁰) and 22Q4 (ρ=+0.443, n=185, p=3.1×10⁻¹⁰).

CGAS RNA does not carry this. The STING-side match is IRF3 pan-cancer and STING1 inside lung. The NHEJ-side match is PAXX (positive). XRCC5 and NHEJ1 are the opposite sign in lung; they are in the lung table below.

The lung IFN-γ result (n=214, ρ=+0.28) and the CosMx exclusion result stay as already reported.

Scatter: `figures/fig_cldn4_rna_vs_irf3_paxx.png`. Within-lineage IRF3: `figures/fig_within_lineage_irf3.png`. Grid: `tables/search_24q4_grid.tsv`. Replication: `tables/cross_release_match.tsv`.

## Dependency search (weaker than the RNA match)

Chronos was searched on 24Q4 and on 24Q2, the release whose integrated `CRISPRGeneEffect.csv` still contains a PRKDC column. CLDN4 Chronos has no dependency tail in either release (24Q4: median +0.042, **0 / 1,178** < −0.5; 24Q2: median +0.071, **0 / 1,150** < −0.5).

24Q2 PRKDC is filled for **24** models (median −1.46, 23/24 < −0.5). CLDN4 vs PRKDC on those 24: Spearman ρ=−0.037, p=0.86. That release change does not produce a PRKDC co-dependency.

The strongest PRKDC gene-effect correlation in the search is the 24Q4 Humagne-CD Cas12 screen matrix, n=40, ρ=−0.367, p=0.020, family *q*=0.040 (four named genes on that matrix). LIG4 on the same 40 screens is ρ=−0.439, p=0.0046. Five of the 40 are lung lines. Integrated pan-cancer LIG4 Chronos is ρ=−0.096 (n=1,178, p=0.0010) and absolute rank 3,251 / 17,915, so it is a small-effect large-n result, not the RNA match above.

**Data.** DepMap Public 24Q4 ([10.25452/figshare.plus.27993248.v1](https://doi.org/10.25452/figshare.plus.27993248.v1)).

- Integrated gene effect: `CRISPRGeneEffect.csv` Chronos, one row per model. More negative = stronger dependency. **1,178** models.
- RNA: `OmicsExpressionProteinCodingGenesTPMLogp1.csv`, log2(TPM+1).
- Screen-level gene effect: `ScreenGeneEffect.csv`. PRKDC lives only here.

**PRKDC is absent from integrated Chronos.** `CRISPRGeneEffect.csv` has LIG4, XRCC4, XRCC5, XRCC6, NHEJ1, DCLRE1C, PAXX, CGAS, STING1, TBK1, and IRF3. It has no PRKDC column. The 24Q4 README describes Humagne-CD as the Cas12 library. `ScreenSequenceMap.csv` labels every `*.CD*` screen as **Humagne-CD**. PRKDC is filled on **40** of those screens (40 models). Five are lung cell lines and two are NSCLC, so the PRKDC essentiality number is not a lung number.

**Honest n.**

| cut | n |
|---|---:|
| Lung models in `Model.csv` (`OncotreeLineage == Lung`, `ModelType == Cell Line`) | 260 |
| Integrated Chronos, all models | 1,178 |
| Integrated Chronos, lung cell lines | 126 |
| Integrated Chronos, NSCLC | 98 |
| Humagne-CD screens with a PRKDC score | 40 |
| Humagne-CD lung / NSCLC | 5 / 2 |
| Lung RNA with finite CLDN4 | **214** |
| NSCLC RNA | 143 |
| LUAD RNA | 80 |

The lung RNA n=214 matches the earlier IFN page (260 lung models, 46 without 24Q4 RNA). On this extract, CLDN4 vs CD274 is ρ=+0.330 (absolute rank 1,320 / 19,172), the same figure as that page. Cultured lines: no infiltrate, no IFN treatment.

**CLDN4 is not a dependency in these matrices.** Integrated Chronos median **+0.042** (IQR −0.036 to +0.113; range −0.39 to +0.48). **0 / 1,178** models have CLDN4 Chronos < −0.5. Lung: median +0.047, **0 / 126**. Humagne-CD: median **+0.37** (range +0.05 to +0.64), **0 / 40**.

STING1 integrated Chronos is also on the non-essential side (median +0.022, **0 / 1,178** < −0.5). CGAS median −0.111, 3 / 1,178 < −0.5. LIG4 median −0.061, 12 / 1,178 < −0.5. XRCC5 and XRCC6 are common essentials (medians −1.22 and −1.75; almost every line < −1) and have little selective range. On Humagne-CD, PRKDC is a real dependency: median **−0.91**, **37 / 40** < −0.5.

Gene sets. Classical NHEJ: PRKDC, LIG4, XRCC4, XRCC5, XRCC6, NHEJ1, DCLRE1C, PAXX. cGAS–STING machinery: CGAS, STING1, TBK1, IRF3. cGAS is the HGNC symbol CGAS. IFN target genes are not in the STING set. Signatures are the mean of within-cohort gene-wise z-scores, complete cases. Spearman is primary. BH *q* is inside the stated family on that cohort, not across the genome. Bootstrap 95% CI: 5,000 resamples, seed 0.

## Integrated Chronos (LIG4, STING1, CGAS)

PRKDC cannot be tested here. Primary family = the three named genes that are in the matrix. A positive ρ means the two Chronos scores move together.

| cohort | gene | n | ρ | p | q | 95% CI |
|---|---|---:|---:|---:|---:|---|
| all models | LIG4 | 1,178 | **−0.096** | 0.0010 | 0.0030 | [−0.15, −0.04] |
| all models | STING1 | 1,178 | **+0.068** | 0.019 | 0.029 | [+0.01, +0.13] |
| all models | CGAS | 1,178 | −0.021 | 0.48 | 0.48 | [−0.08, +0.04] |
| lung | LIG4 | 126 | −0.077 | 0.39 | 0.41 | [−0.26, +0.13] |
| lung | STING1 | 126 | +0.080 | 0.38 | 0.41 | [−0.11, +0.26] |
| lung | CGAS | 126 | −0.075 | 0.41 | 0.41 | [−0.27, +0.12] |
| NSCLC | LIG4 | 98 | −0.135 | 0.18 | 0.28 | [−0.34, +0.09] |
| NSCLC | STING1 | 98 | +0.132 | 0.20 | 0.28 | [−0.08, +0.33] |
| NSCLC | CGAS | 98 | −0.110 | 0.28 | 0.28 | [−0.33, +0.11] |
| lineage residual | LIG4 | 1,164 | −0.082 | 0.0049 | 0.015 | [−0.14, −0.02] |
| lineage residual | STING1 | 1,164 | +0.068 | 0.021 | 0.031 | [+0.01, +0.13] |
| lineage residual | CGAS | 1,164 | −0.018 | 0.54 | 0.54 | [−0.08, +0.04] |

Lineage residual = Chronos minus the `OncotreeLineage` median, lineages with n≥8. The CI is a bootstrap of those residuals (lineage medians held fixed).

**Where those rhos sit among all genes.** CLDN4 Chronos vs 17,915 other genes: 5th percentile −0.116, 95th +0.118, largest |ρ| 0.37. LIG4 (−0.096) is absolute rank **3,251 / 17,915** (empirical two-sided 0.18). STING1 is rank 6,143 (0.34). CGAS is rank 13,823 (0.77). The pan-cancer LIG4 and STING1 intervals exclude 0 because n is 1,178; a correlation of this size is common. Lung and NSCLC intervals all include 0. Residualizing lineage does not enlarge them.

The integrated NHEJ mean-z (7 genes, PRKDC omitted) vs CLDN4 is ρ=−0.057, 95% CI [−0.11, +0.001]. The STING machinery mean-z is ρ=−0.041, CI [−0.10, +0.02]. XRCC5 and XRCC6, the flat essentials, are near zero (ρ=−0.016 and +0.026).

## PRKDC essentiality is the Humagne-CD cut

Same screen-level Chronos for both axes. Family of four named genes on this matrix (PRKDC, LIG4, STING1, CGAS), so *q* is within those four.

| gene | n | ρ | p | q | 95% CI |
|---|---:|---:|---:|---:|---|
| PRKDC | **40** | **−0.367** | 0.020 | 0.040 | [−0.60, −0.08] |
| LIG4 | 40 | −0.439 | 0.0046 | 0.018 | [−0.69, −0.14] |
| STING1 | 40 | −0.077 | 0.64 | 0.85 | [−0.40, +0.25] |
| CGAS | 40 | +0.005 | 0.98 | 0.98 | [−0.33, +0.36] |

Negative ρ: screens with a higher CLDN4 score (every one of them still positive) have a more negative PRKDC score. Lung n=5 and NSCLC n=2 are below the pre-set minimum of 8, so no lung ρ is reported. This is a Cas12 screen-level matrix, not the integrated model-level Chronos used above.

## Lung RNA

n=214 lung cell lines. BH *q* is within NHEJ (8 genes) or within STING machinery (4 genes) or within the two signatures.

| gene | family | ρ | p | q | 95% CI | \|ρ\| rank |
|---|---|---:|---:|---:|---|---|
| PRKDC | NHEJ | **−0.218** | 0.0013 | 0.0026 | [−0.35, −0.08] | 4,087 / 19,172 |
| LIG4 | NHEJ | −0.105 | 0.12 | 0.17 | [−0.23, +0.03] | 10,351 |
| XRCC4 | NHEJ | −0.002 | 0.98 | 0.98 | [−0.13, +0.13] | 19,010 |
| XRCC5 | NHEJ | **−0.289** | 1.8×10⁻⁵ | 1.4×10⁻⁴ | [−0.40, −0.17] | 2,017 |
| XRCC6 | NHEJ | −0.134 | 0.050 | 0.079 | [−0.26, −0.00] | 8,356 |
| NHEJ1 | NHEJ | **−0.262** | 1.1×10⁻⁴ | 4.3×10⁻⁴ | [−0.38, −0.14] | 2,678 |
| DCLRE1C | NHEJ | +0.011 | 0.87 | 0.98 | [−0.12, +0.14] | 18,181 |
| PAXX | NHEJ | **+0.230** | 6.9×10⁻⁴ | 0.0019 | [+0.10, +0.35] | 3,654 |
| CGAS | STING | +0.050 | 0.46 | 0.46 | [−0.08, +0.19] | 14,691 |
| STING1 | STING | **+0.309** | 4.1×10⁻⁶ | 1.7×10⁻⁵ | [+0.18, +0.43] | 1,629 |
| TBK1 | STING | −0.101 | 0.14 | 0.19 | [−0.22, +0.03] | 10,675 |
| IRF3 | STING | **+0.251** | 2.1×10⁻⁴ | 4.1×10⁻⁴ | [+0.12, +0.38] | 2,987 |
| NHEJ mean-z | signature | **−0.193** | 0.0046 | 0.0046 | [−0.32, −0.07] | — |
| STING mean-z | signature | **+0.240** | 4.0×10⁻⁴ | 8.1×10⁻⁴ | [+0.11, +0.36] | — |

Background for the ranks: among 19,172 lung-RNA genes, the 95th percentile of CLDN4 Spearman is **+0.331** and the strongest is KRT19 at **+0.766** (rank 1). CLDN7 is +0.721 (rank 4). STING1 (+0.309, rank 1,629) sits next to CD274 (+0.330, rank 1,320), just under that 95th percentile. Epithelial anchors outrank it by a wide margin.

The NHEJ mean-z is negative, and the genes do not share a sign. XRCC5, NHEJ1, and PRKDC are negative; PAXX is positive; LIG4, XRCC4, and DCLRE1C have intervals that include 0. The STING mean-z is STING1 and IRF3; CGAS and TBK1 do not track CLDN4.

**NSCLC (n=143).** STING1 ρ=+0.300 (q=0.0011). CGAS ρ=−0.123 (CI includes 0). STING mean-z shrinks to ρ=+0.129, CI [−0.04, +0.29]. NHEJ mean-z stays ρ=−0.201 (q=0.032). PRKDC ρ=−0.247 (q=0.0077). NHEJ1 ρ=−0.335. PAXX ρ=+0.204.

**LUAD (n=80).** STING1 ρ=+0.372 (q=0.0027). NHEJ mean-z ρ=−0.033 (p=0.77). The inverse NHEJ score is a mixed-lung result; it is gone inside LUAD. PRKDC inside LUAD is ρ=−0.040.

**Quartiles of CLDN4 RNA, lung, 54 vs 54** (thresholds log2(TPM+1) 3.08 / 7.93). STING1 median 4.42 vs 1.98 (Cliff’s δ=+0.43, q=2.1×10⁻⁴). STING mean-z δ=+0.34 (q=0.0041). NHEJ mean-z δ=−0.30 (q=0.0076). PRKDC 6.86 vs 7.42 (δ=−0.26, q=0.045). CGAS and LIG4 quartile tests have *q*>0.2.

**Rank residual on MKI67 or EPCAM, lung n=214.** Partial Spearman; *p* treats the covariate as fixed.

| target | covariate | partial ρ | p |
|---|---|---:|---:|
| STING1 | MKI67 | +0.309 | 4.2×10⁻⁶ |
| STING1 | EPCAM | +0.464 | 9.4×10⁻¹³ |
| STING mean-z | MKI67 | +0.239 | 4.3×10⁻⁴ |
| STING mean-z | EPCAM | +0.315 | 2.8×10⁻⁶ |
| PRKDC | MKI67 | −0.255 | 1.7×10⁻⁴ |
| PRKDC | EPCAM | −0.223 | 0.0011 |
| NHEJ mean-z | MKI67 | −0.259 | 1.3×10⁻⁴ |
| NHEJ mean-z | EPCAM | −0.161 | 0.019 |
| CGAS | MKI67 | +0.050 | 0.47 |

Adjusting for EPCAM makes the STING1 association larger, so it is not an EPCAM shadow. Adjusting for MKI67 leaves it unchanged. The PRKDC and NHEJ-score inverse associations also remain after either covariate. CGAS stays flat after MKI67.

## Reading

The significant match is CLDN4 RNA with PAXX (NHEJ) and IRF3 (cGAS–STING output). It is the minimum p in a 1,824-test 24Q4 grid, it survives that grid’s search-wide BH, it survives a lineage-median residual and an EPCAM residual, and the same two genes repeat on 24Q2 and 22Q4. STING1 joins them inside lung, and only inside lung. CGAS does not.

CLDN4 Chronos has no essential tail, so co-dependency cannot be large. The best PRKDC gene-effect number is Humagne-CD n=40, ρ=−0.367. 24Q2 puts PRKDC back in the integrated file for 24 models and the correlation is ρ=−0.037.

## What this measurement is

Basal Chronos and basal RNA in cultured lines. No immune cells are in the dish, so this does not speak to the CosMx exclusion numbers. Chronos here is not a CLDN4-knockout IFN or NHEJ experiment: CLDN4 is not essential in these lines. The Humagne-CD PRKDC row is a 40-screen Cas12 matrix, kept separate from integrated Chronos because that gene was dropped from `CRISPRGeneEffect.csv`.

## How to rerun

```bash
python3 -m pip install pandas numpy scipy matplotlib
python3 methods/depmap_cldn4_nhej_sting/download.py --outdir data/depmap_cldn4_nhej_sting
python3 methods/depmap_cldn4_nhej_sting/analyze.py --data data/depmap_cldn4_nhej_sting --outdir methods/depmap_cldn4_nhej_sting
python3 methods/depmap_cldn4_nhej_sting/search_match.py
```

`search_match.py` rebuilds the cross-release table from the local 24Q4 extract and re-downloads 24Q2 (RNA + Chronos, including PRKDC) and 22Q4 RNA.

Tables: `tables/cross_release_match.tsv`, `tables/search_24q4_grid.tsv`, `tables/within_lineage_24q4.tsv`, `tables/headline_ci.tsv`, plus the lung Chronos tables from `analyze.py`.

Figures: `figures/fig_cldn4_rna_vs_irf3_paxx.png`, `figures/fig_within_lineage_irf3.png`, `figures/fig_lung_cldn4_expr_vs_nhej_sting.png`, `figures/fig_humagne_cd_cldn4_vs_prkdc.png`.
