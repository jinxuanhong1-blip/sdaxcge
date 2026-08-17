# GSE33532 — LUAD array CLDN4 vs CD8A / CD274 / ImmuneScore

**Additive CLDN4-only.** No dual-high. No prior GSE33532 CLDN4–immune table is re-audited. This folder only measures **CLDN4 vs CD8A, CD274, and ESTIMATE ImmuneScore** on public LUAD **tumor** arrays, after an epithelial residual.

Public Meister / Thoraxklinik Heidelberg surgical series (Meister et al., *J Bioinformatics Research Studies* 2014, 1(1):1; GEO [GSE33532](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33532); GPL570 Affymetrix U133 Plus 2.0, RMA-PLM). Unit for independence is the **patient**. No ICI arm. No slide was re-scored.

## Honest n

Do not write n=100. The series is 20 patients × (4 tumor sub-sites + 1 matched normal). GEO `histology` is **adeno / mixed / squamous / normal lung**. Mixed is not LUAD.

| item | n | rule |
|---|---:|---|
| arrays on the series matrix | 100 | GEO `GSE33532_series_matrix.txt.gz` |
| deposited probe sets | 25906 | filtered GPL570 (full Plus2 is 54,675); not a re-annotation |
| unique GSM | 100 | all unique |
| patients (GEO design) | 20 | `overall_design`: 4 sites A–D + matched normal from 20 patients |
| matched normal lung | 20 | titles “matched normal lung”; **dropped** |
| all tumor arrays | 80 | 20 patients × 4 sites; **not** the LUAD n |
| LUAD tumor arrays (`histology: adeno`) | **40** | 10 patients × 4 sites |
| LUAD patients | **10** | **independent n** |
| mixed tumor arrays | 24 | 6 patients; not LUAD |
| squamous tumor arrays | 16 | 4 patients; not LUAD |
| ICI / response | 0 | resected early-stage atlas; not deposited |
| OS / DFS | 0 | not a GEO characteristic |
| tumor % / ABSOLUTE purity | 0 | only public proxy used here is an RNA epithelial mean-z |
| CLDN4 finite (`201428_at`) | 40 | named Plus2 probe |
| CD8A finite (`205758_at`) | 40 | named Plus2 probe |
| CD274 finite (`223834_at`) | 40 | named Plus2 probe; `227458_at` **absent** from the deposited matrix |
| ESTIMATE ImmuneScore | 40 | Yoshihara Immune141 ssGSEA (138/141 genes on this matrix) |
| epithelial mean-z | 40 | 6/6 of EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7 |

Primary tests use **10 LUAD patients** (mean of the four tumor sites). Array-level n=40 is the same 10 tumors counted four times and is **not** an independent n. HOLDS (PR 229 / GSE4573) requires n≥40, ρ_adj<0, p_adj<0.05 — the independent LUAD n is **UNDERPOWERED**.

## One-row table

| dataset | histology rule | n arrays | n patients | CLDN4–CD8A ρ (p) | adj ρ \| epi (p) | CLDN4–CD274 ρ (p) | adj ρ \| epi (p) | CLDN4–ImmuneScore ρ (p) | adj ρ \| epi (p) | verdict |
|---|---|---:|---:|---|---|---|---|---|---|---|
| GSE33532 LUAD (patient mean) | GEO `adeno` tumors only | 40 | **10** | -0.164 (0.651) | -0.172 (0.658) | +0.042 (0.907) | +0.042 (0.915) | -0.103 (0.777) | -0.106 (0.785) | **UNDERPOWERED** |
| GSE33532 LUAD (array; not independent) | same 10 tumors × 4 sites | 40 | 10 | -0.144 (0.377) | -0.113 (0.494) | +0.097 (0.55) | +0.114 (0.49) | -0.156 (0.336) | -0.119 (0.469) | NO_EVIDENCE (pseudo-replicated) |

Full numeric rows: `tables/one_row.tsv`, `tables/stats.tsv`.

## Matrix and labels (nothing invented)

Deposited series-matrix is RMA-PLM log2 intensity on a **filtered** 25,906-probe GPL570 subset. Spearman is rank-based.

| field | public? | n | what is there |
|---|---|---:|---|
| histology | yes | 100 | GEO `adeno` 40 / `mixed` 24 / `squamous` 16 / `normal lung` 20 |
| tumor stage | yes | 80 tumors | 1A / 1B / 2A / 2B; unused here |
| gender / age | yes | 100 | unused here |
| ICI response | **no** | 0 | surgical / pretreatment diagnostic series |
| ESTIMATE published scores | **no** | 0 | ImmuneScore computed here from Yoshihara lists |

Named Plus2 probes (forced; max-mean collapse is not used for the three named genes):

| gene | probe | on deposited matrix |
|---|---|---|
| CLDN4 | `201428_at` | yes |
| CD8A | `205758_at` | yes |
| CD274 | `223834_at` | yes (`227458_at` missing) |

ESTIMATE gene collapse uses first-symbol mapping so Immune141 coverage is 138/141 on this filtered matrix.

## Epithelial residual

Partial Spearman residualises ranks on the pan-epithelial mean-z (EPCAM, KRT8, KRT18, KRT19, CDH1, KRT7; **6/6** present: EPCAM, KRT8, KRT18, KRT19, CDH1, KRT7). No dual-high TACSTD2×CLDN4 split.

CLDN4 vs epithelial mean-z (LUAD arrays): ρ=+0.123 (p=0.449, n=40).

CD8A vs ImmuneScore (positive-control, LUAD arrays): ρ=+0.855 (p=2.18e-12, n=40).

## Main Spearman (primary = LUAD patient mean, n=10)

Partial = Pearson of rank residuals on epithelial mean-z; df = n − 3. Bootstrap 95% CI, 2,000 resamples, seed `20260817`.

| predictor | endpoint | unit | n | ρ | 95% CI | p | ρ_adj \| epi | p_adj | verdict |
|---|---|---|---:|---:|---|---:|---:|---:|---|
| CLDN4 | CD8A | patient | 10 | **-0.164** | -0.869 to 0.590 | 0.651 | -0.172 | 0.658 | **UNDERPOWERED** |
| CLDN4 | CD274 | patient | 10 | **0.042** | -0.648 to 0.702 | 0.907 | 0.042 | 0.915 | **UNDERPOWERED** |
| CLDN4 | ImmuneScore | patient | 10 | **-0.103** | -0.692 to 0.539 | 0.777 | -0.106 | 0.785 | **UNDERPOWERED** |
| CLDN4 | CD8A | array | 40 | -0.144 | -0.476 to 0.200 | 0.377 | -0.113 | 0.494 | NO_EVIDENCE |
| CLDN4 | CD274 | array | 40 | 0.097 | -0.202 to 0.385 | 0.55 | 0.114 | 0.49 | NO_EVIDENCE |
| CLDN4 | ImmuneScore | array | 40 | -0.156 | -0.425 to 0.117 | 0.336 | -0.119 | 0.469 | NO_EVIDENCE |

Do not treat the 40 LUAD arrays as 40 patients. The independent test is n=10 and is **UNDERPOWERED** for HOLDS. This is **not** an ICI-response test.

## Sensitivity (not the claim)

All-tumor NSCLC (LUAD+MIXED+LUSC) is 80 arrays / 20 patients. Still not n=100.

| pair | unit | n | ρ (p) | adj ρ \| epi (p) | verdict |
|---|---|---:|---|---|---|
| CLDN4 vs CD8A | patient | 20 | -0.250 (0.289) | -0.106 (0.666) | UNDERPOWERED |
| CLDN4 vs CD274 | patient | 20 | -0.041 (0.865) | +0.087 (0.724) | UNDERPOWERED |
| CLDN4 vs CD8A | array | 80 | -0.217 (0.053) | -0.102 (0.373) | NO_EVIDENCE |
| CLDN4 vs CD274 | array | 80 | -0.035 (0.756) | +0.067 (0.559) | NO_EVIDENCE |

## What is not done

- No n=100 (tumor+normal mix) and no n=80 as a LUAD n.
- No mixed or squamous arrays in the primary LUAD row.
- No TACSTD2 / CLDN4 dual-high split.
- No ICI ORR / PFS model (labels are not deposited).
- No survival model (labels are not deposited).

## Files

- `analyze.py` — GEO download, honest-n inventory, ESTIMATE ImmuneScore, epithelial residual, Spearman, figures
- `estimate_yoshihara_2013.tsv` — Yoshihara 2013 Supp Data 1 (141+141)
- `tables/one_row.tsv`, `stats.tsv`, `per_sample.tsv`, `patient_means.tsv`, `probe_used.tsv`, `coverage.tsv`, `histology_inventory.tsv`, `label_inventory.tsv`, `summary.json`
- `figures/fig1_cldn4_vs_cd8a.png`
- `figures/fig2_cldn4_vs_cd274.png`
- `figures/fig3_cldn4_vs_immunescore.png`
- `figures/fig4_spearman_forest.png`

```bash
python3 methods/gse33532_cldn4/analyze.py
```
