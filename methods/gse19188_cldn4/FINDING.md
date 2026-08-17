# GSE19188 NSCLC array — CLDN4 vs CD8A / CD274

**Additive only.** Public Hou resected early-stage NSCLC Affymetrix U133 Plus 2.0 series (Hou et al., *PLoS One* 2010, PMID 20421987; GEO [GSE19188](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19188)). Unit is the **array**. No slide was re-scored. No ICI arm.

Primary pairs are **CLDN4 vs CD8A** and **CLDN4 vs CD274** on **tumor** arrays. Adjacent-normal arrays are inventoried and are **not** mixed into those tests.

## Honest n

| Item | Public? | n | Note |
|---|---|---:|---|
| Arrays in the series matrix | yes | **156** | 54,675 probes × 156 GSM |
| Unique GSM | yes | **156** | all unique |
| Patients (GEO design text) | text | **91** | series `overall_design`: “91 tumor- and 65 adjacent normal … cohort of 91 patients” |
| Tumor arrays (`tissue type: tumor`) | yes | **91** | **PRIMARY n**; ADC 45 / SCC 27 / LCC 19 |
| Adjacent-normal arrays (`tissue type: healthy`) | yes | **65** | not mixed into the NSCLC pairwise tests |
| Title suffix T/N | misleading | 156 | e.g. GSM475659 title `2344T` is `tissue type: healthy`; use the characteristic, not the title |
| Stage | no | 0 | paper is early-stage; stage is not a GEO characteristic |
| OS months / event on tumors | partial | **82** | 9 tumors `Not available`; **not this claim** |
| ICI / treatment | no | 0 | resected atlas, not an ICI series |
| Tumor % / ABSOLUTE purity | no | 0 | only public proxy is an RNA epithelial score |
| CLDN4 finite (`201428_at`) on tumors | yes | **91** | named Plus2 probe; official GPL570 = CLDN4 / Entrez 1364 |
| CD8A finite (`205758_at`) on tumors | yes | **91** | named Plus2 probe; official GPL570 = CD8A / Entrez 925 |
| CD274 finite (`223834_at`) on tumors | yes | **91** | named Plus2 probe; official GPL570 = CD274 / Entrez 29126 |
| **Primary pairwise n (CLDN4 + CD8A tumors)** | yes | **91** | this is the n used below |
| **Primary pairwise n (CLDN4 + CD274 tumors)** | yes | **91** | this is the n used below |

Do not write n=156 for the correlations. Mixing tumor with adjacent-normal would inflate n and is not an NSCLC-tumor test. The computable public tumor n is **91 arrays**.

Deposited matrix: RMA; intensities <30 reset to 30; then log2(ratio to the probe-set geometric mean). Spearman is rank-based, so a further monotone transform would not change ρ.

## One-row table

| dataset | histology | platform | n tumor | n patients (GEO text) | CLDN4 | CD8A | CD274 | CLDN4–CD8A ρ (p) | adj ρ (p) | CLDN4–CD274 ρ (p) | adj ρ (p) | verdict | OS / ICI |
|---|---|---|---:|---:|---|---|---|---|---|---|---|---|---|
| GSE19188 Hou | NSCLC (ADC 45 / SCC 27 / LCC 19) | GPL570 U133 Plus 2.0 RMA log2-ratio | **91** | 91 | `201428_at` | `205758_at` | `223834_at` | **+0.121 (0.25)** | +0.082 (0.44) | **+0.036 (0.74)** | −0.070 (0.51) | **NO_EVIDENCE** / **NO_EVIDENCE** | OS 82/91 deposited; ICI not deposited |

Full numbers: `tables/spearman_cldn4_vs_genes.tsv`, `tables/label_inventory.tsv`.

## CLDN4 vs CD8A / CD274 (primary)

Named Plus2 probes. Bootstrap 95% CI, 2,000 resamples, seed `20260817`. Partial Spearman residualises on the pan-epithelial mean-z (EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7; 6/6 present). HOLDS rule (same as PR 229 / GSE4573): n≥40, ρ_adj<0, p_adj<0.05.

| pair | n | ρ | 95% CI | p | ρ_adj | p_adj | verdict |
|---|---:|---:|---|---:|---:|---:|---|
| CLDN4 vs **CD8A** | **91** | **+0.121** | −0.079 to +0.321 | 0.25 | +0.082 | 0.44 | **NO_EVIDENCE** |
| CLDN4 vs **CD274** | **91** | **+0.036** | −0.163 to +0.235 | 0.74 | −0.070 | 0.51 | **NO_EVIDENCE** |
| CLDN4 Q4 vs Q1 on CD8A | 23 vs 23 | — | — | MWU 0.54 | — | — | null (rank-biserial +0.11) |
| CLDN4 Q4 vs Q1 on CD274 | 23 vs 23 | — | — | MWU 0.91 | — | — | null (rank-biserial −0.02) |

CLDN4-high is not CD8A-low and is not CD274-low on this NSCLC array. Both CIs include a modest negative and a modest positive effect.

Alternate CD274 probe `227458_at` (also GPL570 = CD274; not the headline): ρ = +0.123, p = 0.25; ρ_adj = +0.059, p = 0.58. Same verdict.

## Histology split (honest n; not the primary)

| subset | pair | n | ρ | p | ρ_adj | p_adj | verdict |
|---|---|---:|---:|---:|---:|---:|---|
| ADC | CLDN4 vs CD8A | 45 | −0.158 | 0.30 | −0.230 | 0.13 | NO_EVIDENCE |
| ADC | CLDN4 vs CD274 | 45 | −0.039 | 0.80 | −0.142 | 0.36 | NO_EVIDENCE |
| SCC | CLDN4 vs CD8A | 27 | +0.043 | 0.83 | +0.117 | 0.57 | **UNDERPOWERED** |
| SCC | CLDN4 vs CD274 | 27 | −0.165 | 0.41 | −0.178 | 0.39 | **UNDERPOWERED** |
| LCC | CLDN4 vs CD8A | 19 | +0.119 | 0.63 | +0.236 | 0.35 | **UNDERPOWERED** |
| LCC | CLDN4 vs CD274 | 19 | +0.372 | 0.12 | +0.364 | 0.14 | **UNDERPOWERED** |

ADC is the only histology with n≥40. The ADC residuals are negative but not significant. Do not pool SCC/LCC into a “holds in squamous” sentence.

## CLDN4 vs immune genes (pre-specified panel)

19/19 genes present after first-symbol collapse of official GEO GPL570 annot, then named-probe override for CLDN4 / CD8A / CD274. BH *q* is across the 19 present immune genes (not TACSTD2).

| pair | n | ρ | p | q (BH) | ρ_adj | p_adj | verdict |
|---|---:|---:|---:|---:|---:|---:|---|
| CLDN4 vs CD8A | 91 | +0.121 | 0.25 | 0.40 | +0.082 | 0.44 | NO_EVIDENCE |
| CLDN4 vs CD274 | 91 | +0.036 | 0.74 | 0.74 | −0.070 | 0.51 | NO_EVIDENCE |
| CLDN4 vs CD8B | 91 | +0.047 | 0.66 | 0.73 | −0.002 | 0.99 | NO_EVIDENCE |
| CLDN4 vs PDCD1 | 91 | +0.041 | 0.70 | 0.74 | +0.001 | 0.99 | NO_EVIDENCE |
| CLDN4 vs CD2 | 91 | +0.179 | 0.089 | 0.30 | +0.150 | 0.16 | NO_EVIDENCE |
| CLDN4 vs CD3D | 91 | +0.161 | 0.13 | 0.30 | +0.129 | 0.23 | NO_EVIDENCE |
| CLDN4 vs CD3E | 91 | +0.166 | 0.12 | 0.30 | +0.149 | 0.16 | NO_EVIDENCE |
| CLDN4 vs GZMA | 91 | +0.104 | 0.32 | 0.47 | +0.095 | 0.37 | NO_EVIDENCE |
| CLDN4 vs GZMB | 91 | +0.175 | 0.098 | 0.30 | +0.095 | 0.37 | NO_EVIDENCE |
| CLDN4 vs GZMK | 91 | +0.080 | 0.45 | 0.61 | +0.101 | 0.34 | NO_EVIDENCE |
| CLDN4 vs PRF1 | 91 | +0.166 | 0.12 | 0.30 | +0.149 | 0.16 | NO_EVIDENCE |
| CLDN4 vs IFNG | 91 | +0.122 | 0.25 | 0.40 | +0.111 | 0.30 | NO_EVIDENCE |
| CLDN4 vs STAT1 | 91 | +0.170 | 0.11 | 0.30 | +0.027 | 0.80 | NO_EVIDENCE |
| CLDN4 vs CXCL9 | 91 | +0.147 | 0.17 | 0.35 | +0.078 | 0.46 | NO_EVIDENCE |
| CLDN4 vs CXCL10 | 91 | +0.052 | 0.63 | 0.73 | +0.004 | 0.97 | NO_EVIDENCE |
| CLDN4 vs IDO1 | 91 | +0.192 | 0.068 | 0.30 | +0.120 | 0.26 | NO_EVIDENCE |
| CLDN4 vs HLA-DRA | 91 | +0.219 | 0.037 | 0.30 | +0.279 | 0.0077 | OPPOSITE (adj) |
| CLDN4 vs CXCL13 | 91 | +0.062 | 0.56 | 0.71 | +0.037 | 0.73 | NO_EVIDENCE |
| CLDN4 vs MS4A1 | 91 | +0.125 | 0.24 | 0.40 | +0.088 | 0.41 | NO_EVIDENCE |

Every crude immune-gene BH *q* is ≥0.30. The one adjusted *p*<0.05 (HLA-DRA, ρ_adj = +0.28) is a **positive** residual after eating the epithelial score; it is not the primary pair and is not a CD8 or PD-L1 claim.

## CLDN4 vs locked immune signatures

Score = mean of gene-wise z. Coverage is Plus2-complete for these short lists.

| signature | n | genes used | ρ | p | ρ_adj | p_adj | verdict |
|---|---:|---:|---:|---:|---:|---:|---|
| T_cell_CD8 | 91 | 8/8 | +0.145 | 0.17 | +0.120 | 0.26 | **NO_EVIDENCE** |
| IFNg_Ayers | 91 | 6/6 | +0.163 | 0.12 | +0.115 | 0.28 | NO_EVIDENCE |

Do not call a CD8 or IFN-γ effect for CLDN4 on GSE19188.

## What CLDN4 *does* track here

| pair | n | ρ | 95% CI | p |
|---|---:|---:|---|---:|
| CLDN4 vs epithelial mean-z | 91 | **+0.614** | 0.452 to 0.735 | 9.8×10⁻¹¹ |
| CLDN4 vs TACSTD2 | 91 | **+0.621** | 0.458 to 0.752 | 4.9×10⁻¹¹ |

CLDN4 sits on the epithelial / TACSTD2 side of this array. That is compatible with a tight-junction / tumour-cell program. It is not evidence that CLDN4-high NSCLC is immune-low or PD-L1-high.

## TACSTD2 extras (companion; not this claim)

| pair | n | ρ | p | ρ_adj | p_adj | verdict |
|---|---:|---:|---:|---:|---:|---|
| TACSTD2 vs CD8A | 91 | +0.224 | 0.033 | +0.206 | 0.052 | NO_EVIDENCE |
| TACSTD2 vs CD274 | 91 | +0.156 | 0.14 | +0.102 | 0.34 | NO_EVIDENCE |

TACSTD2–CD8A is a weak crude positive that does not hold after the epithelial residual. Do not treat CLDN4 as interchangeable with TACSTD2 on this matrix; both track epithelium, and neither is CD8A-low / CD274-low after that residual.

## What this does not test

- ICI response or PFS (not on GEO).
- OS as a CLDN4 endpoint (82/91 tumors have time+event; not scored here).
- Pathologist CD8 / PD-L1 IHC.
- A histology-wide CLDN4–immune law. SCC n=27 and LCC n=19 are under the n≥40 rule.
- n=156 (tumor+normal mix).

## Reproduce

```bash
python3 -m pip install -r methods/gse19188_cldn4/requirements.txt
python3 methods/gse19188_cldn4/analyze.py
```

Downloads (not committed) go to `$GSE19188_CLDN4_DATA` (default `/tmp/gse19188_cldn4`):

- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE19nnn/GSE19188/matrix/GSE19188_series_matrix.txt.gz`
- `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570&targ=self&form=text&view=data`

## Files

- `analyze.py` — download, probe map, Spearman / partial Spearman, figures
- `tables/label_inventory.tsv` — public vs missing fields; honest n
- `tables/probe_confirm.tsv` — named CLDN4 / CD8A / CD274 / TACSTD2 probes
- `tables/histology_counts.tsv` — ADC 45 / SCC 27 / LCC 19
- `tables/spearman_cldn4_vs_genes.tsv` — CLDN4 vs CD8A, CD274, and the immune panel
- `tables/spearman_cldn4_by_histology.tsv` — ADC / SCC / LCC splits
- `tables/spearman_cldn4_vs_signatures.tsv` — CD8 / IFN-γ / epithelial
- `tables/spearman_tacstd2_companion.tsv` — TACSTD2 vs CD8A / CD274 only
- `tables/spearman_cd274_alt_probe.tsv` — `227458_at` companion
- `tables/highlow_cldn4.tsv` — CLDN4 Q4 vs Q1 on CD8A and CD274 (23 vs 23)
- `tables/sample_annotation.tsv` — 156 arrays
- `tables/summary.json`
- `figures/fig1_cldn4_vs_cd8a_cd274.png`
- `figures/fig2_cldn4_immune_forest.png`
- `figures/fig3_cldn4_vs_cd8sig.png`
