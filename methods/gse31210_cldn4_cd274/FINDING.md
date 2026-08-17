# GSE31210 — CLDN4 vs CD274 and HLA-A/B/C (partial on ESTIMATE)

**Additive only. New cut.** CLDN4–CD8A ρ on this series is **already known** (`methods/gse31210_cldn4_immune`) and is not re-cut. This folder only measures CLDN4 against **CD274** and **HLA-A / HLA-B / HLA-C**, with a rank residual on ESTIMATEScore.

Public Okayama / Kohno Japanese stage I–II LUAD (Okayama et al., *Cancer Res* 2012, PMID 22261853; GEO [GSE31210](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31210); GPL570 Affymetrix U133 Plus 2.0). Unit is the **primary lung tumor array**. No ICI arm. No slide was re-scored.

## Honest n

| item | n | rule |
|---|---:|---|
| arrays on the series matrix | 246 | GEO `GSE31210_series_matrix.txt.gz` |
| primary lung tumors | **226** | `tissue: primary lung tumor` |
| adjacent / non-tumor | 20 | dropped |
| CLDN4 / CD274 / HLA-A / HLA-B / HLA-C | 226 | all five genes present after collapse |
| ESTIMATEScore | 226 | Yoshihara Stromal141+Immune141 ssGSEA (136/141 + 138/141) |
| ESTIMATE TumorPurity in [0, 1] | 18 | cosine map; 208 tumors wrap outside [0, 1] |

Primary tests use **n=226** tumors. That is the same tumor rule as the sibling CD8 cut; it is not a new cohort hunt.

## One-row table

| dataset | n | CLDN4–CD274 ρ (p) | CLDN4–HLA-A ρ (p) | CLDN4–HLA-B ρ (p) | CLDN4–HLA-C ρ (p) | CLDN4–CD274 partial \| ESTIMATE (p) | CLDN4–HLA-A partial \| ESTIMATE (p) | CLDN4–HLA-B partial \| ESTIMATE (p) | CLDN4–HLA-C partial \| ESTIMATE (p) |
|---|---:|---|---|---|---|---|---|---|---|
| GSE31210 LUAD tumors | 226 | -0.284 (1.46e-05) | -0.092 (0.169) | -0.058 (0.385) | -0.138 (0.0377) | -0.087 (0.195) | -0.052 (0.434) | 0.007 (0.92) | -0.042 (0.531) |

MHC-I mean-z (HLA-A/B/C): unadj -0.100 (0.134); partial \| ESTIMATE -0.028 (0.674).

Full numeric rows: `tables/one_row.tsv`, `tables/stats.tsv`.

## Matrix and labels (nothing invented)

GEO series-matrix MAS5 linear intensity. Tests use the native ranks (Spearman / ssGSEA), so log vs linear does not change ρ. Compact extracts store **log2(MAS5+1)** for plots.

| field | public? | n | what is there |
|---|---|---:|---|
| tissue | yes | 246 | primary lung tumor 226 / other 20 |
| histology | yes | 226 | LUAD (series is stage I–II lung adenocarcinoma) |
| ICI response | **no** | 0 | treatment-naive surgical cohort |
| ESTIMATE published scores | **no** | 0 | computed here from Yoshihara lists |

Named-gene collapse (unique-mapped max-mean when a unique probe exists; otherwise first-symbol max-mean):

| gene | probe | rule |
|---|---|---|
| CLDN4 | `201428_at` | unique-mapped max-mean |
| CD274 | `227458_at` | unique-mapped max-mean |
| HLA-A | `213932_x_at` | unique-mapped max-mean |
| HLA-B | `209140_x_at` | unique-mapped max-mean |
| HLA-C | `216526_x_at` | unique-mapped max-mean |

ESTIMATE gene collapse uses first-symbol mapping (same as the sibling cut).

## ESTIMATE

R `estimate` is not required. Scores follow the sibling public implementation:

- ssGSEA (Barbie/GSVA, τ=0.25) on Yoshihara 2013 Stromal141 + Immune141
- ranks scaled to 1…10000
- `ESTIMATEScore = StromalScore + ImmuneScore`
- `TumorPurity = cos(0.6049872018 + 0.0001467884 × ESTIMATEScore)`

**Partial covariate is ESTIMATEScore**, not ImmuneScore. ImmuneScore is one addend of ESTIMATEScore and is collinear with it. On this matrix TumorPurity remains monotone with ESTIMATEScore (angles still in (0, π)), so the rank residual on TumorPurity equals the rank residual on ESTIMATEScore. The cosine is below 0 for **208 / 226** tumors and is **not** used as a 0–1 fraction. The [0, 1] subset (n=18) is only the least-impure tail and is **not** the analysis n.

CLDN4 vs ImmuneScore is reported only as a context row (same matrix as the sibling cut). It is not a new CD8 test.

## Main Spearman (n=226)

Partial = Pearson of rank residuals on ESTIMATEScore; df = n − 3.

| predictor | endpoint | n | ρ | p | partial ρ \| ESTIMATE | partial p |
|---|---|---:|---:|---:|---:|---:|
| CLDN4 | CD274 | 226 | **-0.284** | 1.46e-05 | -0.087 | 0.195 |
| CLDN4 | HLA-A | 226 | **-0.092** | 0.169 | -0.052 | 0.434 |
| CLDN4 | HLA-B | 226 | **-0.058** | 0.385 | 0.007 | 0.92 |
| CLDN4 | HLA-C | 226 | **-0.138** | 0.0377 | -0.042 | 0.531 |
| CLDN4 | MHC-I (HLA-A/B/C mean-z) | 226 | -0.100 | 0.134 | -0.028 | 0.674 |
| CLDN4 | ImmuneScore | 226 | -0.366 | 1.50e-08 | -0.031 | 0.639 |

Unadjusted CLDN4–CD274 is negative (ρ=-0.284, p=1.46e-05). HLA-A and HLA-B are null; HLA-C is weakly negative (p=0.0377). After the ESTIMATEScore residual every named endpoint is null (CD274 partial p=0.195). That is the same attenuation pattern as the already-known CD8 residual on this matrix. Antigen-presentation / PD-L1 cut only; not an ICI-response test. CLDN4–CD8A is not restated.

## What is not done

- No re-cut of CLDN4 vs CD8A / TACSTD2 (sibling folder).
- No OncoSG or CPTAC re-run.
- No ICI ORR / PFS model (labels are not deposited).
- No R `estimate` purity transform beyond the published cosine on the computed ESTIMATEScore.

## Files

- `analyze.py` — GEO download, tumor filter, ESTIMATE, Spearman, figures
- `estimate_yoshihara_2013.tsv` — Yoshihara 2013 Supp Data 1 (141+141)
- `tables/one_row.tsv`, `stats.tsv`, `per_sample.tsv`, `probe_used.tsv`, `coverage.tsv`, `summary.json`
- `figures/fig1_cldn4_vs_cd274.png`
- `figures/fig2_cldn4_vs_hla.png`
- `figures/fig3_spearman_forest.png`
- `figures/fig4_correlation_heatmap.png`

```bash
python3 methods/gse31210_cldn4_cd274/analyze.py
```
