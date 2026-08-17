# GSE10072 — Landi LUAD array, CLDN4-only vs CD8A / CD274 / ImmuneScore

**Additive only. CLDN4-only. No dual-high.** Public Landi EAGLE LUAD Affymetrix HG-U133A series (Landi et al., *PLoS ONE* 2008, PMID 18297132; GEO [GSE10072](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10072); GPL96). Unit is the **array**. Tumor only. No slide was re-scored. No ICI arm. TACSTD2 is not analysed.

Primary tests: **CLDN4 vs CD8A**, **CLDN4 vs CD274**, **CLDN4 vs ESTIMATE ImmuneScore**, each after an **epithelial residual**.

## Honest n

Do not write n=135 or n=122. The series summary still says 135 tissues; the GEO `overall_design` then drops low-tumour-cell arrays and averages QC duplicates. The deposited series matrix is **107 arrays**.

| item | public? | n | rule |
|---|---|---:|---|
| selected tissues (paper / design text) | text only | 180 | includes QC duplicates / triplicates from 14 subjects |
| high-quality RNA | text only | 148 | 180 minus insufficient RNA |
| normalized microarrays | text only | 135 | 148 minus 13 problematic assays |
| after low-tumour-cell drop | text only | 122 | 135 minus 13 low % tumor |
| after averaging 15 QC duplicates | text only | **107** | GEO `overall_design`; 58 tumor + 49 non-tumor |
| **arrays in the series matrix** | yes | **107** | 22,283 probes × 107 GSM; 0 missing RMA values |
| unique GSM / unique titles | yes | **107** | `Lung Tumor_GT*` / `Normal Lung_GT*`; all unique |
| unique subjects (`GT` id) | yes | **74** | 20 never / 26 former / 28 current (design text) |
| adenocarcinoma tumor (`source_name`) | yes | **58** | `Adenocarcinoma of the Lung`; 58 unique `GT` ids |
| non-tumor lung | yes | **49** | `Normal Lung Tissue`; **dropped** (tumor-only rule) |
| paired tumor+normal subjects | yes | 33 | same `GT` in both sources; pairing unused (tumor-only) |
| tumor-only subjects | yes | 25 | tumor array, no deposited normal |
| tumor smoking | yes | 58 | Current 24 / Former 18 / Never 16 |
| tumor stage | yes | 58 | IA 5, IB 17, IIA 3, IIB 18, IIIA 9, IIIB 3, IV 3 |
| tumor gender | yes | 58 | Male 35 / Female 23 |
| ICI / treatment / response | **no** | 0 | smoking / LUAD atlas, not an ICI series |
| OS / DSS time or event | **no** | 0 | paper discusses survival; labels **not** in the matrix |
| numeric tumor % / ABSOLUTE purity | **no** | 0 | 13 arrays dropped for low % tumor; values not deposited |
| CLDN4 finite (`201428_at`) | yes | **58** | single unique-mapped GPL96 probe |
| CD8A finite (`205758_at`) | yes | **58** | single unique-mapped GPL96 probe |
| CD274 finite | **no** | **0** | not on GPL96 HG-U133A (no CD274 / PDCD1LG1 / B7-H1 probe) |
| ESTIMATE ImmuneScore | computed | **58** | Yoshihara Immune141 ssGSEA; 138/141 genes |
| epithelial mean-z | computed | **58** | EPCAM/KRT8/KRT18/KRT19/CDH1/KRT7 (6/6) |
| **primary pairwise n (tumor LUAD)** | yes | **58** | complete-case CLDN4 + CD8A + ImmuneScore |

Primary tests use **n=58 tumor LUAD arrays**. Do not pool the 49 normals. Do not write the paper’s 135 / 122.

## One-row table

| dataset | histology rule | platform | n tumors | CLDN4 | CD8A | CD274 | ImmuneScore | CLDN4–CD8A ρ (p) | adj ρ \| epi (p) | CLDN4–ImmuneScore ρ (p) | adj ρ \| epi (p) | CLDN4–CD274 | verdict |
|---|---|---|---:|---|---|---|---|---|---|---|---|---|---|
| GSE10072 Landi | LUAD tumor only (normals dropped) | GPL96 HG-U133A RMA | **58** | `201428_at` | `205758_at` | **absent (U133A)** | ESTIMATE Immune141 (138/141) | **-0.481 (0.000132)** | -0.423 (0.00103) | **-0.177 (0.183)** | -0.124 (0.357) | not on U133A | CD8A **HOLDS**; ImmuneScore **NO_EVIDENCE** |

Full numeric rows: `tables/one_row.tsv`, `tables/stats.tsv`.

## Matrix and labels (nothing invented)

GEO series-matrix **RMA** (Bioconductor `affy`; `!Sample_data_processing`). Tests use native ranks (Spearman / ssGSEA), so the deposited log-intensity scale does not change ρ.

| field | public? | n | what is there |
|---|---|---:|---|
| histology | yes | 107 | every tumor is LUAD; every normal is non-tumor lung |
| smoking | yes | 107 | GEO `Cigarette Smoking Status` |
| stage / age / gender | yes | 107 | GEO characteristics; **not** the claim |
| ICI response | **no** | 0 | not an ICI series |
| ESTIMATE published scores | **no** | 0 | computed here from Yoshihara lists |

Max-mean unique-mapped GPL96 probe (multi-mapped `///` probes dropped for named genes):

| gene | probe | on GPL96? |
|---|---|---|
| CLDN4 | `201428_at` | yes |
| CD8A | `205758_at` | yes |
| CD274 | — | **no** (Plus 2.0 probes `223834_at` / `227458_at` are not on U133A) |

ESTIMATE gene collapse uses first-symbol mapping so Immune141 coverage is 138/141 and Stromal141 is 136/141.

## Epithelial residual and ImmuneScore

- **Epithelial residual:** unweighted mean of gene-wise z for EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7 (6/6 present: EPCAM/KRT8/KRT18/KRT19/CDH1/KRT7).
- **ImmuneScore:** ssGSEA (Barbie/GSVA, τ=0.25) on Yoshihara 2013 Immune141. R `estimate` is not required.
- **Partial Spearman:** Pearson of rank residuals on the epithelial mean-z; df = n − 3.
- **HOLDS rule** (same as PR 229 / GSE4573 / GSE10245): n≥40, ρ_adj<0, p_adj<0.05.
- **No dual-high:** no TACSTD2×CLDN4 quadrant, no dual-high vs dual-low contrast.

CD8A vs ImmuneScore (positive-control, tumor LUAD): ρ=+0.560 (p=4.96e-06, n=58).

CLDN4 vs epithelial mean-z (context, not the claim): ρ=+0.339 (p=0.00934, n=58).

## Main Spearman (primary tumor LUAD n=58)

Bootstrap 95% CI, 2,000 resamples, seed `20260817`.

| predictor | endpoint | n | ρ | 95% CI | p | ρ_adj \| epi | p_adj | verdict |
|---|---|---:|---:|---|---:|---:|---:|---|
| CLDN4 | CD8A | 58 | **-0.481** | -0.660 to -0.244 | 0.000132 | -0.423 | 0.00103 | **HOLDS** |
| CLDN4 | CD274 | 0 | — | — | — | — | — | **ABSENT_ON_GPL96** |
| CLDN4 | ImmuneScore | 58 | **-0.177** | -0.412 to +0.071 | 0.183 | -0.124 | 0.357 | **NO_EVIDENCE** |
| CD8A | ImmuneScore | 58 | +0.560 | +0.331 to +0.731 | 4.96e-06 | +0.538 | 1.61e-05 | positive control |
| CLDN4 | epithelial_z | 58 | +0.339 | +0.095 to +0.557 | 0.00934 | — | — | context |

This is **not** an ICI-response test. CD274 cannot be scored on this chip. Smoking-stratum n (Current 24 / Former 18 / Never 16) is below the HOLDS n≥40 rule and is not a claim.

## What is not done

- No dual-high TACSTD2×CLDN4 split.
- No pooling of the 49 non-tumor arrays into the primary n.
- No invented CD274 probe or PD-L1 IHC substitute.
- No ICI ORR / PFS model (labels are not deposited).
- No OS model (survival labels are not in the series matrix).
- No re-use of the paper’s 135-sample sentence as the analysis n.

## Files

- `analyze.py` — GEO download, honest-n inventory, ESTIMATE, epithelial residual, Spearman, figures
- `estimate_yoshihara_2013.tsv` — Yoshihara 2013 Supp Data 1 (141+141)
- `tables/one_row.tsv`, `stats.tsv`, `per_sample.tsv`, `probe_used.tsv`, `coverage.tsv`, `label_inventory.tsv`, `summary.json`
- `figures/fig1_cldn4_vs_cd8a.png`
- `figures/fig2_cldn4_vs_immunescore.png`
- `figures/fig3_spearman_forest.png`

```bash
python3 -m pip install -r methods/gse10072_cldn4/requirements.txt
export GSE10072_CLDN4_DATA=/tmp/gse10072_cldn4
python3 methods/gse10072_cldn4/analyze.py
```
