# GSE10245 NSCLC array — CLDN4 vs CD8A / CD274

**Additive only.** Public Kuner resected NSCLC Affymetrix U133 Plus 2.0 series (Kuner et al., *Lung Cancer* 2009, PMID 18486272; GEO [GSE10245](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10245)). Unit is the **array**. No slide was re-scored. No ICI arm. This slice does not re-claim the paper’s junction / histology DE.

Primary question: **CLDN4 vs CD8A** and **CLDN4 vs CD274**. TACSTD2 is a same-run companion only.

## Honest n

Do not write n=48. Lung Cancer Explorer lists 48 tumors for this accession; that is an external reprocessed subset, not the GEO series matrix.

| Item | Public? | n | Note |
|---|---|---:|---|
| Arrays in the series matrix | yes | **58** | 54,675 probes × 58 GSM; 0 missing gcRMA values |
| Unique GSM / unique titles | yes | **58** | `NSCLC_AC_*` / `NSCLC_SCC_*`; all unique |
| Patients | yes (1 title = 1 array) | **58** | GEO `overall_design` is “40 AC and 18 SCC”; no duplicate-patient sentence |
| Adenocarcinoma | yes | **40** | `disease state: adenocarcinoma` |
| Squamous cell carcinoma | yes | **18** | `disease state: squamous cell carcinoma` |
| Normals / paired adjacent | no | 0 | every `source_name` is tumor tissue |
| Stage | no | 0 | not a GEO characteristic |
| OS / DSS time or event | no | 0 | junction / histology paper; labels not deposited |
| ICI / treatment | no | 0 | resected atlas, not an ICI series |
| Numeric tumor % / ABSOLUTE purity | no | 0 | protocol text says tumor/stromal content was estimated then macro-dissected; values not deposited |
| CLDN4 finite (`201428_at`) | yes | **58** | max-mean chose the named probe; `1569421_at` is background (mean 2.22) |
| CD8A finite (`205758_at`) | yes | **58** | only CD8A probe on GPL570 |
| CD274 finite (`227458_at`) | yes | **58** | max-mean of two clean probes (`223834_at`, `227458_at`) |
| **Primary pairwise n (all NSCLC)** | yes | **58** | complete-case CLDN4 + CD8A + CD274 |
| Primary pairwise n (ADC) | yes | **40** | HOLDS-rule eligible (n≥40) |
| Primary pairwise n (SCC) | yes | **18** | **UNDERPOWERED** for HOLDS (n≥40) |

The computable public n is **58 arrays**. Use 40 / 18 when the row is histology-stratified. Do not pool a silent n=48.

## One-row table

| dataset | histology | platform | n arrays | CLDN4 | CD8A | CD274 | CLDN4–CD8A ρ (p) | adj ρ (p) | CLDN4–CD274 ρ (p) | adj ρ (p) | verdict | OS / ICI |
|---|---|---|---:|---|---|---|---|---|---|---|---|---|
| GSE10245 Kuner | NSCLC 40 ADC + 18 SCC | GPL570 U133 Plus 2.0 gcRMA | **58** | `201428_at` | `205758_at` | `227458_at` | **−0.115 (0.39)** | +0.006 (0.96) epi; −0.162 (0.23) histo | **−0.105 (0.43)** | −0.049 (0.72) epi; −0.004 (0.97) histo | **NO_EVIDENCE** | not deposited |

Full numbers: `tables/spearman.tsv`, `tables/label_inventory.tsv`.

## CLDN4 vs CD8A / CD274 (primary)

gcRMA as deposited. Spearman is rank-based. Bootstrap 95% CI, 2,000 resamples, seed `20260817`. Partial Spearman residualises on the pan-epithelial mean-z (EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7; **6/6** present) or on a binary ADC indicator (pooled NSCLC only). HOLDS rule (same as PR 229 / GSE4573): n≥40, ρ_adj<0, p_adj<0.05.

| pair | subset | n | ρ | 95% CI | p | ρ_adj epi (p) | ρ_adj histo (p) | verdict |
|---|---|---:|---:|---|---:|---|---|---|
| CLDN4 vs **CD8A** | all NSCLC | **58** | **−0.115** | −0.362 to +0.130 | 0.39 | +0.006 (0.96) | −0.162 (0.23) | **NO_EVIDENCE** |
| CLDN4 vs **CD274** | all NSCLC | **58** | **−0.105** | −0.384 to +0.182 | 0.43 | −0.049 (0.72) | −0.004 (0.97) | **NO_EVIDENCE** |
| CLDN4 vs **CD8A** | ADC | **40** | −0.170 | −0.434 to +0.129 | 0.29 | −0.047 (0.78) | — | **NO_EVIDENCE** |
| CLDN4 vs **CD274** | ADC | **40** | +0.089 | −0.277 to +0.419 | 0.58 | +0.129 (0.43) | — | **NO_EVIDENCE** |
| CLDN4 vs **CD8A** | SCC | **18** | −0.079 | −0.614 to +0.429 | 0.75 | −0.056 (0.83) | — | **UNDERPOWERED** |
| CLDN4 vs **CD274** | SCC | **18** | −0.003 | −0.496 to +0.478 | 0.99 | −0.040 (0.88) | — | **UNDERPOWERED** |

CLDN4-high is not CD8A-low or CD274-low on this mixed NSCLC array. Every HOLDS-eligible CI includes both a modest negative and a modest positive effect.

Named-probe CD8A (`201428_at` vs `205758_at`) is identical to the collapsed row (single CD8A probe). Named-probe CD274 `223834_at` is the same sign and still null (all NSCLC ρ=−0.151, p=0.26; epi adj −0.159, p=0.24).

### Q4 vs Q1 (descriptive)

| subset | endpoint | n Q4 vs Q1 | MWU p | rank-biserial |
|---|---|---|---:|---:|
| all NSCLC | CD8A | 15 vs 15 | 0.32 | −0.22 |
| all NSCLC | CD274 | 15 vs 15 | 0.48 | −0.16 |
| ADC | CD8A | 10 vs 10 | 0.16 | −0.38 |
| ADC | CD274 | 10 vs 10 | 0.62 | +0.14 |

SCC quartiles are 5 vs 5 and are not a claim.

## Context (not the claim)

The deposited series is a histology contrast. CLDN4 is higher in ADC than SCC on this matrix (median 9.40 vs 7.92; MWU p=5.2e-4; n=40 vs 18). That matches the Kuner junction paper and is **not** an immune finding. CD8A does not differ by histology (p=0.64). CD274 is not significant (p=0.098).

CLDN4 tracks the epithelial mean-z (all NSCLC ρ=+0.544, p=1.0e-5, n=58). After that residual the CD8A and CD274 associations go to zero.

## TACSTD2 companion (not the claim)

| pair | subset | n | ρ (p) | ρ_adj epi (p) | verdict |
|---|---|---:|---|---|---|
| TACSTD2 vs CD8A | all NSCLC | 58 | −0.222 (0.094) | −0.185 (0.17) | NO_EVIDENCE |
| TACSTD2 vs CD274 | all NSCLC | 58 | −0.030 (0.82) | −0.006 (0.97) | NO_EVIDENCE |
| CLDN4 vs TACSTD2 | all NSCLC | 58 | +0.225 (0.090) | +0.137 (0.31) | NO_EVIDENCE |
| CLDN4 vs TACSTD2 | ADC | 40 | +0.388 (0.013) | +0.223 (0.17) | NO_EVIDENCE |

Do not treat CLDN4 as interchangeable with TACSTD2. Neither gene has a HOLDS residual vs CD8A or CD274 here.

## Methods

- **Matrix:** GEO `GSE10245_series_matrix.txt.gz` (author gcRMA, Bioconductor). Probe × sample values used as published.
- **Annotation:** official GPL570 `Gene symbol` (NCBI platform table, first `///` token). **Max-mean** probe collapse to HUGO.
- **CD8** = `CD8A`. **CD274** = collapsed max-mean (`227458_at`); `223834_at` is a named-probe sensitivity row.
- **Epithelial residual:** unweighted mean of gene-wise z for EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7.
- **Histology residual:** binary ADC indicator on the pooled 58 only. Not applied inside ADC or SCC.
- **Partial Spearman:** Pearson of rank residuals; df = n − 2 − k.
- **Not done:** no ESTIMATE (no numeric purity deposited; samples were macro-dissected). No ICI / OS model. No gene-set fishing.

## Files

- `tables/spearman.tsv` — all n / ρ / p / CI / residuals
- `tables/label_inventory.tsv` — honest n
- `tables/sample_annotation.tsv` — per-array genes + histology
- `tables/probe_confirm.tsv` / `gene_coverage.tsv` / `histology_contrast.tsv` / `highlow_cldn4.tsv` / `summary.json`
- `figures/fig1_cldn4_vs_cd8a_cd274.png` / `fig2_forest.png`

## Reproduce

```bash
python3 -m pip install -r methods/gse10245_cldn4/requirements.txt
export GSE10245_CLDN4_DATA=/tmp/gse10245_cldn4
python3 methods/gse10245_cldn4/analyze.py
```
