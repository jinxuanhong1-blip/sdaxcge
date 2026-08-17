# GSE4573 LUSC array — CLDN4 vs CD274 / HLA-A/B/C

**Additive CD274 / classical MHC-I cut.** CLDN4 vs CD8A on this same matrix is already **NS** in [PR 313](https://github.com/jinxuanhong1-blip/sdaxcge/pull/313) (n=130, ρ=−0.054, p=0.54) and is **not** the claim here. Do not headline CD8.

Public Raponi resected squamous lung carcinoma Affymetrix U133A series (Raponi et al., *Cancer Res* 2006, PMID 16885343; GEO [GSE4573](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4573)). Unit is the **array**. No slide was re-scored. No ICI arm.

## Honest n

| Item | Public? | n | Note |
|---|---|---:|---|
| Arrays in the series matrix | yes | **130** | 22,283 probes × 130 GSM; 0 missing MAS5 values |
| Unique GSM / unique `LS-*` titles | yes | **130** | all unique |
| Patients (GEO design text) | text only | **129** | series `overall_design`: “130 samples from 129 patients” |
| Duplicate patient pair | **no** | — | not identified in sample metadata; **not dropped** |
| Histology | series-level | 130 | every sample described as lung squamous carcinoma; no LUAD rows |
| Stage / OS / ICI / tumor % | no | 0 | not on GEO |
| CLDN4 finite (`201428_at`) | yes | **130** | official GPL96 = CLDN4 / Entrez 1364 |
| **CD274 finite (Entrez 29126)** | **no** | **0** | no GPL96 probe; Plus-2 `223834_at` / `227458_at` are not on U133A |
| HLA-A / HLA-B / HLA-C finite | yes | **130** | chosen max-mean probes `215313_x_at` / `209140_x_at` / `216526_x_at` |
| **Primary pairwise n (CLDN4 + CD274)** | **no** | **0** | ABSENT — do not invent a PD-L1 ρ |
| **Primary pairwise n (CLDN4 + HLA-A/B/C)** | yes | **130** | this is the n used for HLA / MHC-I |

Do not write n=129 for the correlations. The computable public n is **130 arrays**. The 129-patient sentence is GEO text only.

PDCD1LG2 (PD-L2, `220049_s_at`, Entrez 80380) **is** on U133A. It is not CD274 and is **not** substituted.

## One-row table

| dataset | histology | platform | n arrays | n patients (GEO text) | CLDN4 | CD274 | HLA-A/B/C | CLDN4–CD274 ρ (p) | CLDN4–MHC-I ρ (p) | MHC-I \| epi ρ (p) | MHC-I \| CD8A ρ (p) | verdict |
|---|---|---|---:|---:|---|---|---|---|---|---|---|---|
| GSE4573 Raponi | LUSC | GPL96 U133A MAS5 | **130** | 129 (pair not ID’d) | `201428_at` | **ABSENT** (Entrez 29126 not on U133A) | 3/3 | **ABSENT (n=0)** | **+0.026 (0.77)** | +0.124 (0.16) | +0.063 (0.48) | **NO_EVIDENCE** (HLA); **ABSENT** (CD274) |

Full numbers: `tables/one_row.tsv`, `tables/spearman_cldn4_vs_cd274_hla.tsv`, `tables/label_inventory.tsv`.

## New CD274 cut

MAS5 as deposited (linear signal). Spearman is rank-based. Bootstrap 95% CI, 2,000 resamples, seed `20260817`. MHC-I mean-z = gene-wise z of **HLA-A/B/C only** (3/3; not B2M/TAP). Partial Spearman residualises on the pan-epithelial mean-z (EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7; 6/6) or on CD8A (covariate only).

| pair | n | ρ | 95% CI | p | ρ \| epi (p) | ρ \| CD8A (p) | verdict |
|---|---:|---:|---|---:|---|---|---|
| CLDN4 vs **CD274** | **0** | — | — | — | — | — | **ABSENT** |
| CLDN4 vs HLA-A | **130** | −0.001 | −0.168 to +0.181 | 0.99 | +0.114 (0.20) | +0.025 (0.77) | NO_EVIDENCE |
| CLDN4 vs HLA-B | 130 | +0.043 | −0.137 to +0.220 | 0.63 | +0.119 (0.18) | +0.083 (0.35) | NO_EVIDENCE |
| CLDN4 vs HLA-C | 130 | +0.023 | −0.149 to +0.205 | 0.79 | +0.106 (0.23) | +0.059 (0.50) | NO_EVIDENCE |
| CLDN4 vs MHC-I (HLA-A/B/C mean-z) | 130 | **+0.026** | −0.142 to +0.200 | **0.77** | +0.124 (0.16) | +0.063 (0.48) | **NO_EVIDENCE** |

CLDN4 Q4 vs Q1 (33 vs 33) on the same endpoints:

| endpoint | n Q4 vs Q1 | rank-biserial (Q4−Q1) | MWU p | verdict |
|---|---|---:|---:|---|
| CD274 | — | — | — | **ABSENT** |
| HLA-A | 33 vs 33 | −0.052 | 0.72 | NO_EVIDENCE |
| HLA-B | 33 vs 33 | +0.017 | 0.91 | NO_EVIDENCE |
| HLA-C | 33 vs 33 | −0.023 | 0.88 | NO_EVIDENCE |
| MHC-I mean-z | 33 vs 33 | −0.005 | 0.98 | NO_EVIDENCE |

There is no public CD274 (PD-L1) measurement on this U133A series. Classical MHC-I is present and is **null** versus CLDN4 at honest n=130. Residualising on epithelium or on CD8A does not create a significant MHC-I association.

## What this does not claim

- A CLDN4–CD8 result. That pair is already NS in PR 313 and is not re-headlined.
- A PD-L1 RNA claim (CD274 is not on the chip).
- A PD-L2 (PDCD1LG2) substitute for CD274.
- ICI response, PFS, or OS (not on GEO).
- Patient-level n=129 (duplicate ID is not public).

## Reproduce

```bash
python3 -m pip install -r methods/gse4573_cldn4_cd274/requirements.txt
python3 methods/gse4573_cldn4_cd274/analyze.py
```

Downloads (not committed) go to `$GSE4573_CLDN4_CD274_DATA` (default `/tmp/gse4573_cldn4_cd274`):

- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE4nnn/GSE4573/matrix/GSE4573_series_matrix.txt.gz`
- `https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz`

## Files

- `analyze.py` — download, exact-symbol / Entrez probe map, Spearman / partial Spearman, figures
- `tables/one_row.tsv` — claim row
- `tables/label_inventory.tsv` — public vs missing fields; honest n
- `tables/probe_confirm.tsv` — CD274 absent; HLA-A/B/C probes; PDCD1LG2 noted, not used
- `tables/spearman_cldn4_vs_cd274_hla.tsv`
- `tables/highlow_cldn4_cd274_hla.tsv` — CLDN4 Q4 vs Q1
- `tables/gene_coverage.tsv` / `sample_annotation.tsv` / `summary.json`
- `figures/fig1_cldn4_vs_cd274.png` — ABSENT call
- `figures/fig2_cldn4_cd274_hla_forest.png`
- `figures/fig3_cldn4_vs_mhci.png`
