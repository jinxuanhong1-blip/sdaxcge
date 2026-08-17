# GSE4573 LUSC array — CLDN4 vs CD8A / immune genes

**Additive only.** TLS / CD8 **TACSTD2** extras on this same matrix are already in [PR 229](https://github.com/jinxuanhong1-blip/sdaxcge/pull/229) and are treated as given. This slice is **CLDN4**.

Public Raponi resected squamous lung carcinoma Affymetrix U133A series (Raponi et al., *Cancer Res* 2006, PMID 16885343; GEO [GSE4573](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4573)). Unit is the **array**. No slide was re-scored. No ICI arm.

## Honest n

| Item | Public? | n | Note |
|---|---|---:|---|
| Arrays in the series matrix | yes | **130** | 22,283 probes × 130 GSM; 0 missing MAS5 values |
| Unique GSM / unique `LS-*` titles | yes | **130** | all unique |
| Patients (GEO design text) | text only | **129** | series `overall_design`: “130 samples from 129 patients” |
| Duplicate patient pair | **no** | — | not identified in sample metadata; **not dropped** |
| Histology | series-level | 130 | every sample described as lung squamous carcinoma; no LUAD rows |
| Stage | no | 0 | not a GEO characteristic |
| OS / DSS time or event | no | 0 | prognosis paper; labels not deposited; GEO has no `suppl/` folder |
| ICI / treatment | no | 0 | resected atlas, not an ICI series |
| Tumor % / ABSOLUTE purity | no | 0 | only public proxy is an RNA epithelial score |
| CLDN4 finite (`201428_at`) | yes | **130** | named U133A probe; official GPL96 annot = CLDN4 / Entrez 1364 |
| CD8A finite (`205758_at`) | yes | **130** | named U133A probe; official GPL96 annot = CD8A / Entrez 925 |
| **Primary pairwise n (CLDN4 + CD8A)** | yes | **130** | this is the n used below |

Do not write n=129 for the correlations. The computable public n is **130 arrays**. The 129-patient sentence is GEO text only.

## One-row table

| dataset | histology | platform | n arrays | n patients (GEO text) | CLDN4 | CD8A | CD8B | CLDN4–CD8A ρ (p) | adj ρ (p) | verdict | OS / ICI |
|---|---|---|---:|---:|---|---|---|---|---|---|---|
| GSE4573 Raponi | LUSC | GPL96 U133A MAS5 | **130** | 129 (pair not ID’d) | `201428_at` | `205758_at` | **absent** (dirty `LOC100996919///CD8B`) | **−0.054 (0.54)** | −0.052 (0.56) | **NO_EVIDENCE** | not deposited |

Full numbers: `tables/spearman_cldn4_vs_genes.tsv`, `tables/label_inventory.tsv`.

## CLDN4 vs CD8A (primary)

MAS5 as deposited (linear signal). Spearman is rank-based, so a log2 transform would not change ρ. Bootstrap 95% CI, 2,000 resamples, seed `20260816`. Partial Spearman residualises on the pan-epithelial mean-z (EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7; 6/6 present). HOLDS rule (same as PR 229): n≥40, ρ_adj<0, p_adj<0.05.

| pair | n | ρ | 95% CI | p | ρ_adj | p_adj | verdict |
|---|---:|---:|---|---:|---:|---:|---|
| CLDN4 vs **CD8A** | **130** | **−0.054** | −0.227 to +0.117 | 0.54 | −0.052 | 0.56 | **NO_EVIDENCE** |
| CLDN4 Q4 vs Q1 on CD8A | 33 vs 33 | — | — | MWU 0.35 | — | — | null (rank-biserial −0.13) |

CLDN4-high is not CD8A-low on this LUSC array. The CI includes both a modest negative and a modest positive effect.

## CLDN4 vs immune genes (pre-specified panel)

21/22 genes present after first-symbol collapse of official GEO GPL96 annot (2016-08-09). **CD8B is absent**: probes `207979_s_at` and `215332_s_at` are annotated `LOC100996919///CD8B` and are not mapped to CD8B. BH *q* is across the 21 present immune genes (not TACSTD2).

| pair | n | ρ | p | q (BH) | ρ_adj | p_adj | verdict |
|---|---:|---:|---:|---:|---:|---:|---|
| CLDN4 vs CD8A | 130 | −0.054 | 0.54 | 0.98 | −0.052 | 0.56 | NO_EVIDENCE |
| CLDN4 vs CD8B | 0 | — | — | — | — | — | **ABSENT** |
| CLDN4 vs CD2 | 130 | +0.067 | 0.45 | 0.98 | +0.115 | 0.20 | NO_EVIDENCE |
| CLDN4 vs CD3D | 130 | +0.046 | 0.60 | 0.98 | +0.073 | 0.41 | NO_EVIDENCE |
| CLDN4 vs CD3E | 130 | +0.013 | 0.88 | 0.98 | +0.100 | 0.26 | NO_EVIDENCE |
| CLDN4 vs GZMA | 130 | −0.013 | 0.88 | 0.98 | −0.013 | 0.88 | NO_EVIDENCE |
| CLDN4 vs GZMB | 130 | −0.027 | 0.76 | 0.98 | +0.032 | 0.72 | NO_EVIDENCE |
| CLDN4 vs GZMK | 130 | +0.012 | 0.89 | 0.98 | +0.007 | 0.93 | NO_EVIDENCE |
| CLDN4 vs PRF1 | 130 | −0.119 | 0.18 | 0.98 | −0.079 | 0.37 | NO_EVIDENCE |
| CLDN4 vs IFNG | 130 | −0.067 | 0.45 | 0.98 | −0.031 | 0.73 | NO_EVIDENCE |
| CLDN4 vs STAT1 | 130 | +0.016 | 0.85 | 0.98 | +0.074 | 0.41 | NO_EVIDENCE |
| CLDN4 vs CXCL9 | 130 | +0.043 | 0.63 | 0.98 | +0.012 | 0.90 | NO_EVIDENCE |
| CLDN4 vs CXCL10 | 130 | +0.014 | 0.87 | 0.98 | +0.066 | 0.46 | NO_EVIDENCE |
| CLDN4 vs IDO1 | 130 | +0.068 | 0.44 | 0.98 | +0.118 | 0.18 | NO_EVIDENCE |
| CLDN4 vs HLA-DRA | 130 | +0.110 | 0.21 | 0.98 | +0.179 | 0.042 | OPPOSITE (adj only) |
| CLDN4 vs CXCL13 | 130 | −0.076 | 0.39 | 0.98 | −0.038 | 0.67 | NO_EVIDENCE |
| CLDN4 vs CCL19 | 130 | +0.039 | 0.66 | 0.98 | +0.079 | 0.38 | NO_EVIDENCE |
| CLDN4 vs CCL21 | 130 | −0.068 | 0.44 | 0.98 | −0.010 | 0.91 | NO_EVIDENCE |
| CLDN4 vs MS4A1 | 130 | +0.002 | 0.98 | 0.98 | +0.059 | 0.51 | NO_EVIDENCE |
| CLDN4 vs CD19 | 130 | −0.099 | 0.26 | 0.98 | −0.018 | 0.84 | NO_EVIDENCE |
| CLDN4 vs CD79A | 130 | −0.003 | 0.97 | 0.98 | +0.064 | 0.47 | NO_EVIDENCE |
| CLDN4 vs CD79B | 130 | −0.057 | 0.52 | 0.98 | −0.008 | 0.93 | NO_EVIDENCE |

Every crude immune-gene test has BH *q* = 0.98. The one adjusted *p*<0.05 (HLA-DRA, ρ_adj = +0.18) is a weak **positive** residual after eating the epithelial score; it is not the primary pair and is not a CD8 claim.

## CLDN4 vs locked immune signatures (PR 229 lists)

Same gene lists as `scripts/opus_tls/opus_tls_lib.py` in PR 229. Score = mean of gene-wise z. Coverage is U133A-honest (CCL3, RBP5, CD8B, BTLA, TOX2 and a few B/plasma genes are missing).

| signature | n | genes used | ρ | p | ρ_adj | p_adj | verdict |
|---|---:|---:|---:|---:|---:|---:|---|
| T_cell_CD8 | 130 | 7/8 (no CD8B) | +0.010 | 0.91 | +0.044 | 0.62 | **NO_EVIDENCE** |
| IFNg_Ayers | 130 | 6/6 | +0.059 | 0.51 | +0.114 | 0.20 | NO_EVIDENCE |
| TLS_Cabrita | 130 | 8/8 | +0.013 | 0.89 | +0.063 | 0.48 | NO_EVIDENCE |
| TLS_12chemokine | 130 | 11/12 (no CCL3) | −0.056 | 0.53 | +0.024 | 0.79 | NO_EVIDENCE |
| TLS_imprint | 130 | 8/9 (no RBP5) | +0.021 | 0.81 | +0.101 | 0.26 | NO_EVIDENCE |
| B_cell | 130 | 14/16 | −0.017 | 0.85 | +0.074 | 0.40 | NO_EVIDENCE |
| Plasma_cell | 130 | 11/14 | +0.050 | 0.57 | +0.104 | 0.24 | NO_EVIDENCE |
| Tfh | 130 | 7/9 | −0.180 | 0.040 | −0.130 | 0.14 | NO_EVIDENCE |

Tfh is the only crude *p*<0.05 among the locked signatures. It does not survive epithelial residual or BH (q = 0.32). Do not call a TLS or CD8 effect for CLDN4 on GSE4573.

## What CLDN4 *does* track here

| pair | n | ρ | 95% CI | p |
|---|---:|---:|---|---:|
| CLDN4 vs epithelial mean-z | 130 | **+0.425** | 0.264 to 0.558 | 4.8×10⁻⁷ |
| CLDN4 vs TACSTD2 | 130 | **+0.366** | 0.205 to 0.506 | 1.9×10⁻⁵ |

CLDN4 sits on the epithelial / TACSTD2 side of this array. That is compatible with a tight-junction / tumour-cell program. It is not evidence that CLDN4-high LUSC is immune-low.

## TACSTD2 extras (given; not this claim)

PR 229 already scored TACSTD2 on GSE4573 (n=130): CD8 signature / 12-chemokine / IFN-γ **HOLDS** after epithelial correction; Cabrita TLS and B-cell do **not**. This slice does not re-claim those signature tests.

A single-gene companion computed here (not in the PR 229 headline table): TACSTD2 vs **CD8A** ρ = −0.191 (CI −0.345 to −0.028), p = 0.030; ρ_adj = −0.190, p = 0.031 (HOLDS). That is the contrast: **TACSTD2–CD8A is negative; CLDN4–CD8A is not.** Do not treat CLDN4 as interchangeable with TACSTD2 on this matrix.

## What this does not test

- ICI response, PFS, or OS (not on GEO).
- Pathologist TLS or CD8 IHC.
- Patient-level n=129 (duplicate ID is not public).
- A LUSC-wide CLDN4–immune law. This is one public U133A series.

## Reproduce

```bash
python3 -m pip install -r methods/gse4573_cldn4/requirements.txt
python3 methods/gse4573_cldn4/analyze.py
```

Downloads (not committed) go to `$GSE4573_CLDN4_DATA` (default `/tmp/gse4573_cldn4`):

- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE4nnn/GSE4573/matrix/GSE4573_series_matrix.txt.gz`
- `https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz`

## Files

- `analyze.py` — download, probe map, Spearman / partial Spearman, figures
- `tables/label_inventory.tsv` — public vs missing fields; honest n
- `tables/probe_confirm.tsv` — named CLDN4 / CD8A / TACSTD2 probes; dirty CD8B
- `tables/gene_coverage.tsv` — U133A presence for locked lists
- `tables/spearman_cldn4_vs_genes.tsv` — CLDN4 vs CD8A and the immune panel
- `tables/spearman_cldn4_vs_signatures.tsv` — locked TLS / CD8 / IFN signatures
- `tables/spearman_tacstd2_companion.tsv` — TACSTD2 vs CD8A only
- `tables/highlow_cldn4_cd8a.tsv` — CLDN4 Q4 vs Q1 on CD8A (33 vs 33)
- `tables/sample_annotation.tsv` — 130 arrays
- `tables/summary.json`
- `figures/fig1_cldn4_vs_cd8a.png`
- `figures/fig2_cldn4_immune_forest.png`
- `figures/fig3_cldn4_vs_cd8sig.png`
