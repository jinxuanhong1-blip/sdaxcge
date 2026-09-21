# DepMap 24Q4 lung: CLDN4 expression vs NHEJ and cGAS–STING

**Additive public evidence.** Cultured lung lines only. This page is CLDN4 RNA against four NHEJ genes (`PRKDC`, `XRCC4`, `LIG4`, `TP53BP1`) and three cGAS–STING-axis genes (`STING1`, `CGAS`, `STAT1`), plus Chronos co-dependency on the CRISPR overlap. It does not restate the already reported all-lung CLDN4–Hallmark IFN-γ Spearman.

**Data.** DepMap Public 24Q4 ([10.6084/m9.figshare.27993248](https://doi.org/10.6084/m9.figshare.27993248)). RNA = `OmicsExpressionProteinCodingGenesTPMLogp1.csv`, log2(TPM+1), the CCLE/DepMap expression matrix. Dependency = `CRISPRGeneEffect.csv` Chronos (more negative = stronger dependency). Expression columns use the HGNC symbols above (`STING1`, not `TMEM173`; `CGAS`, not `MB21D1`). Integrated Chronos is missing **PRKDC** (no `PRKDC` / `XRCC7` / Entrez 5591 column). PRKDC remains only in `ScreenGeneEffect.csv`, and only on Humagne-CD screens.

**Cohort.** `OncotreeLineage == Lung` and `ModelType == Cell Line`: **260** models in `Model.csv` → **214** with finite RNA for CLDN4 and the seven partners. NSCLC **143** (LUAD 80 / LUSC 27 / other NSCLC 36). SCLC+NET **60**. Other lung **11**. CRISPR overlap with finite integrated Chronos for CLDN4 and every partner that exists in `CRISPRGeneEffect` (CLDN4, XRCC4, LIG4, TP53BP1, STING1, CGAS, STAT1): **123** lung lines (**95** NSCLC). Humagne-CD screens with a finite PRKDC gene effect: **5** lung lines. No immune infiltrate and no IFN treatment in the dish.

Scores are the mean of gene-wise z-scores computed **inside each cut**. NHEJ = PRKDC, XRCC4, LIG4, TP53BP1. STING = STING1, CGAS, STAT1. SENSOR = STING1 and CGAS only, because STAT1 is an interferon-stimulated gene and a Hallmark IFN-γ member. Spearman is two-sided, complete cases. BH *q* is within the family on that cut (seven genes; three scores; Chronos pairs). Bootstrap 95% CI: 5000 resamples. A gene is called only when q<0.05 and |ρ|≥0.15.

## CLDN4 RNA vs each gene

Each cell is Spearman ρ (BH *q* within that cut’s seven-gene family, or within that cut’s three scores). Master file, with p and bootstrap CIs: `tables/expr_correlations.tsv`.

| gene | lung n=214 | NSCLC n=143 | LUAD n=80 | LUSC n=27 | SCLC_NET n=60 | other_NSCLC n=36 |
|---|---|---|---|---|---|---|
| PRKDC | -0.218 (0.00456) | -0.247 (0.0101) | -0.040 (0.849) | -0.374 (0.275) | +0.098 (0.53) | -0.470 (0.0268) |
| XRCC4 | -0.002 (0.979) | -0.064 (0.456) | +0.002 (0.983) | -0.085 (0.786) | -0.129 (0.53) | -0.230 (0.311) |
| LIG4 | -0.105 (0.219) | -0.063 (0.456) | -0.178 (0.342) | -0.218 (0.641) | -0.080 (0.546) | +0.257 (0.311) |
| TP53BP1 | -0.159 (0.0457) | -0.131 (0.251) | -0.043 (0.849) | -0.344 (0.275) | -0.108 (0.53) | -0.090 (0.7) |
| STING1 | +0.309 (2.9e-05) | +0.300 (0.00193) | +0.372 (0.00472) | +0.167 (0.68) | +0.353 (0.0392) | +0.240 (0.311) |
| CGAS | +0.050 (0.648) | -0.123 (0.251) | -0.164 (0.342) | -0.140 (0.68) | +0.234 (0.168) | -0.007 (0.967) |
| STAT1 | -0.008 (0.979) | -0.066 (0.456) | -0.079 (0.849) | +0.041 (0.839) | +0.255 (0.168) | -0.115 (0.7) |
| NHEJ mean-z | -0.257 (0.000426) | -0.275 (0.00267) | -0.156 (0.322) | -0.566 (0.00628) | -0.096 (0.466) | -0.273 (0.323) |
| STING mean-z | +0.170 (0.0126) | +0.077 (0.36) | +0.096 (0.398) | +0.042 (0.835) | +0.429 (0.000935) | +0.072 (0.677) |
| SENSOR mean-z | +0.239 (0.00064) | +0.125 (0.204) | +0.140 (0.322) | +0.061 (0.835) | +0.452 (0.000855) | +0.095 (0.677) |

Bootstrap 95% CIs on the two large cuts: lung PRKDC [-0.35, -0.08], lung STING1 [+0.18, +0.43], NSCLC PRKDC [-0.41, -0.08], NSCLC STING1 [+0.14, +0.44]. Score CIs: lung NHEJ ρ=-0.257, p=0.000142, q=0.000426, 95% CI [-0.38, -0.13], n=214; lung STING ρ=+0.170, p=0.0126, q=0.0126, 95% CI [+0.04, +0.30], n=214; lung SENSOR ρ=+0.239, p=0.000427, q=0.00064, 95% CI [+0.11, +0.36], n=214; NSCLC NHEJ ρ=-0.275, p=0.00089, q=0.00267, 95% CI [-0.42, -0.12], n=143; NSCLC STING ρ=+0.077, p=0.36, q=0.36, 95% CI [-0.09, +0.23], n=143; NSCLC SENSOR ρ=+0.125, p=0.136, q=0.204, 95% CI [-0.04, +0.28], n=143.

NHEJ mean-z vs STING mean-z on all lung lines: ρ=-0.010, p=0.883, q=0.883, 95% CI [-0.15, +0.13], n=214. NHEJ vs SENSOR: ρ=-0.081, p=0.236, q=0.473, 95% CI [-0.22, +0.06]. On NSCLC the same NHEJ-vs-STING score correlation is ρ=+0.192, q=0.0435.

Partial Spearman (rank-residual) on all lung lines: CLDN4 vs STING score given NHEJ score ρ=+0.174 (p=0.0109, n=214); CLDN4 vs NHEJ score given STING score ρ=-0.259 (p=0.000125, n=214); CLDN4 vs SENSOR score given NHEJ score ρ=+0.226 (p=0.000861, n=214).

CLDN4 RNA quartile split (top vs bottom quartile) on the lung RNA cohort: sig_NHEJ: Q4 median -0.220 vs Q1 +0.082 (Cliff δ=-0.34, p=0.00206, n=54/54); sig_STING: Q4 median +0.267 vs Q1 +0.053 (Cliff δ=+0.22, p=0.0533, n=54/54); sig_SENSOR: Q4 median +0.513 vs Q1 -0.040 (Cliff δ=+0.33, p=0.0028, n=54/54); PRKDC: Q4 median +6.864 vs Q1 +7.417 (Cliff δ=-0.26, p=0.0225, n=54/54); XRCC4: Q4 median +3.602 vs Q1 +3.493 (Cliff δ=+0.04, p=0.728, n=54/54); LIG4: Q4 median +2.704 vs Q1 +2.846 (Cliff δ=-0.12, p=0.286, n=54/54); TP53BP1: Q4 median +4.572 vs Q1 +4.825 (Cliff δ=-0.22, p=0.0475, n=54/54); STING1: Q4 median +4.421 vs Q1 +1.975 (Cliff δ=+0.43, p=0.000103, n=54/54); CGAS: Q4 median +3.154 vs Q1 +3.224 (Cliff δ=+0.02, p=0.866, n=54/54); STAT1: Q4 median +5.713 vs Q1 +5.762 (Cliff δ=-0.05, p=0.66, n=54/54).

`sig_*` columns in `tables/lung_lines.tsv` are z-scored on the all-lung cohort. Correlations in the table above recompute z inside each cut.

## What the expression cut shows

Calls use q<0.05 and |ρ|≥0.15 inside the family for that cut.

The gene that stays positive with CLDN4 RNA in the mixed lung set, in NSCLC, and inside LUAD is **STING1** (lung ρ=+0.309, LUAD ρ=+0.372, q=0.00472, n=80). SCLC+NET is the same direction (ρ=+0.353, q=0.0392, n=60). **CGAS** is not called in any of those cuts. **STAT1** is not called in the mixed lung set or in NSCLC. STAT1 is one Hallmark IFN-γ gene; this page does not recompute that score, and a null STAT1 rank is not a revision of it.

The SENSOR mean-z (STING1+CGAS) is called on the mixed lung set and is not called in NSCLC or LUAD. CGAS dilutes STING1. The three-gene STING score is the same story: called on all-lung, not called in NSCLC.

**PRKDC** RNA is negatively associated with CLDN4 on the mixed lung set (ρ=-0.218, q=0.00456) and on mixed NSCLC (ρ=-0.247, q=0.0101). It is not a LUAD association (ρ=-0.040, q=0.849, n=80). The NSCLC number is carried by the non-LUAD slices (other NSCLC ρ=-0.470, q=0.0268, n=36; LUSC ρ=-0.374, q=0.275, n=27). **XRCC4** and **LIG4** are not called. **TP53BP1** meets the call rule only on the mixed lung set and not in NSCLC. The NHEJ mean-z is therefore not a four-gene block. It is called on mixed lung and mixed NSCLC, not on LUAD (ρ=-0.156, q=0.322). LUSC NHEJ mean-z is ρ=-0.566, q=0.00628, n=27.

NHEJ mean-z and STING mean-z are uncorrelated on all lung lines (ρ=-0.010). On NSCLC they are weakly positive (ρ=+0.192, q=0.0435). That is not an NHEJ-up / STING-down pair.

Lung gene calls: PRKDC negative at q<0.05 (ρ=-0.218, q=0.00456); XRCC4 not q<0.05 (ρ=-0.002, q=0.979); LIG4 not q<0.05 (ρ=-0.105, q=0.219); TP53BP1 negative at q<0.05 (ρ=-0.159, q=0.0457); STING1 positive at q<0.05 (ρ=+0.309, q=2.9e-05); CGAS not q<0.05 (ρ=+0.050, q=0.648); STAT1 not q<0.05 (ρ=-0.008, q=0.979).

NSCLC gene calls: PRKDC negative at q<0.05 (ρ=-0.247, q=0.0101); XRCC4 not q<0.05 (ρ=-0.064, q=0.456); LIG4 not q<0.05 (ρ=-0.063, q=0.456); TP53BP1 not q<0.05 (ρ=-0.131, q=0.251); STING1 positive at q<0.05 (ρ=+0.300, q=0.00193); CGAS not q<0.05 (ρ=-0.123, q=0.251); STAT1 not q<0.05 (ρ=-0.066, q=0.456).

## CRISPR co-dependency

Chronos on the lung CRISPR∩panel set (n=123). Profile SD is the standard deviation of the gene effect across those lines.

| gene | n | median Chronos | SD | fraction < −0.5 |
|---|---:|---:|---:|---:|
| CLDN4 | 123 | +0.050 | 0.130 | 0.000 |
| XRCC4 | 123 | -0.236 | 0.132 | 0.041 |
| LIG4 | 123 | -0.049 | 0.136 | 0.008 |
| TP53BP1 | 123 | -0.133 | 0.191 | 0.033 |
| STING1 | 123 | +0.016 | 0.091 | 0.000 |
| CGAS | 123 | -0.104 | 0.099 | 0.000 |
| STAT1 | 123 | -0.007 | 0.086 | 0.000 |

CLDN4 Chronos median +0.050, SD 0.130, fraction < −0.5 = 0.000. That profile is too tight for a co-dependency scan to mean selective CLDN4 essentiality. The CLDN4-Chronos correlations are still reported, as the available co-dependency numbers:

- CLDN4 Chronos vs partner Chronos, lung: XRCC4 ρ=-0.024 (n=123, p=0.789, q=0.789); LIG4 ρ=-0.060 (n=123, p=0.513, q=0.633); TP53BP1 ρ=+0.059 (n=123, p=0.52, q=0.633); STING1 ρ=+0.108 (n=123, p=0.236, q=0.633); CGAS ρ=-0.057 (n=123, p=0.528, q=0.633); STAT1 ρ=-0.116 (n=123, p=0.202, q=0.633)
- CLDN4 RNA vs partner Chronos, lung (expression versus dependency, not co-dependency): XRCC4 ρ=+0.165 (n=123, p=0.0679, q=0.226); LIG4 ρ=-0.036 (n=123, p=0.691, q=0.763); TP53BP1 ρ=-0.161 (n=123, p=0.0754, q=0.226); STING1 ρ=+0.056 (n=123, p=0.535, q=0.763); CGAS ρ=+0.028 (n=123, p=0.763, q=0.763); STAT1 ρ=-0.043 (n=123, p=0.64, q=0.763)

NHEJ×STING Chronos pairs use the integrated genes that exist (XRCC4, LIG4, TP53BP1 × STING1, CGAS, STAT1; BH within that 3×3 family). On the lung CRISPR cut the minimum *q* in that 3×3 is 0.51. Table: `tables/crispr_nhej_vs_sting.tsv`. Read each ρ next to the Chronos SD above. A flat profile cannot support a selective co-dependency.

**PRKDC Chronos is not in the integrated matrix.** Humagne-CD `ScreenGeneEffect` has a finite PRKDC value for 5 lung cell lines: NCI-H841 (other_lung, Chronos -1.152); DMS 273 (SCLC_NET, Chronos -0.210); RERF-LC-Ad2 (LUAD, Chronos -0.924); NCI-H23 (LUAD, Chronos -1.750); NCI-H1882 (SCLC_NET, Chronos -1.112). n=5 is too small for a Spearman. Those screen-level values are listed in `tables/prkdc_humagne_lung.tsv` and are not mixed into the integrated co-dependency table. Guide maps: Avana and KY have no PRKDC guides; Humagne-CD has two guides marked `UsedByChronos=True`.

## What this is not

- Not a tumor immune-exclusion result. These are cell lines.
- Not an IFN-stimulated or ICI result. STAT1 here is basal RNA.
- Not a restatement of the all-lung CLDN4–Hallmark IFN-γ Spearman. STAT1 is one gene inside that axis.
- Not evidence that knocking out CLDN4 changes NHEJ or STING. CLDN4 Chronos in these lines has little spread.
- Not a DNA-PK–versus–cGAS mechanism result. PRKDC dependency is almost absent from this release, and the expression ranks are reported as ranks.

Scatter and forest: `figures/fig_cldn4_nhej_sting.png`.

## How to rerun

```bash
python3 -m pip install -r methods/depmap_nhej_sting/requirements.txt
python3 methods/depmap_nhej_sting/download.py --outdir data/depmap_nhej_sting
python3 methods/depmap_nhej_sting/analyze.py --data data/depmap_nhej_sting --outdir methods/depmap_nhej_sting
```
