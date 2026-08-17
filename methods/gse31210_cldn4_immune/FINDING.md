# GSE31210 — CLDN4 vs CD8A / ImmuneScore / TACSTD2 (after ESTIMATE)

**Additive only.** A9 extra-cohort 3/3 catalog (OncoSG LUAD, CPTAC LUAD RNA, GSE31210 LUAD) is **taken as given** and is not re-audited. This folder only measures CLDN4 against CD8A, ESTIMATE ImmuneScore, and TACSTD2 on the catalogued GSE31210 tumor set.

Public Okayama / Kohno Japanese stage I–II LUAD (Okayama et al., *Cancer Res* 2012, PMID 22261853; GEO [GSE31210](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31210); GPL570 Affymetrix U133 Plus 2.0). Unit is the **primary lung tumor array**. No ICI arm. No slide was re-scored.

## Honest n

| item | n | rule |
|---|---:|---|
| arrays on the series matrix | 246 | GEO `GSE31210_series_matrix.txt.gz` |
| primary lung tumors | **226** | `tissue: primary lung tumor` |
| adjacent / non-tumor | 20 | dropped |
| CLDN4 / CD8A / TACSTD2 | 226 | all three genes present after max-mean collapse |
| ESTIMATE ImmuneScore | 226 | Yoshihara Immune141 ssGSEA (138/141 genes) |
| ESTIMATE TumorPurity in [0, 1] | 18 | cosine map; 208 tumors wrap outside [0, 1] |

Primary tests use **n=226** tumors. The A9 catalog n=226 is the same tumor rule; it is not re-derived as a new cohort hunt.

## One-row table

| dataset | n | CLDN4–CD8A ρ (p) | CLDN4–ImmuneScore ρ (p) | CLDN4–TACSTD2 ρ (p) | CLDN4–CD8A partial \| purity (p) | CLDN4–ImmuneScore partial \| purity (p) |
|---|---:|---|---|---|---|---|
| GSE31210 LUAD tumors | 226 | -0.341 (1.45e-07) | -0.366 (1.50e-08) | 0.471 (7.48e-14) | -0.127 (0.0574) | -0.031 (0.639) |

Full numeric rows: `tables/one_row.tsv`, `tables/stats.tsv`.

## Matrix and labels (nothing invented)

GEO series-matrix MAS5 linear intensity. Tests use the native ranks (Spearman / ssGSEA), so log vs linear does not change ρ. Compact extracts store **log2(MAS5+1)** for plots.

| field | public? | n | what is there |
|---|---|---:|---|
| tissue | yes | 246 | primary lung tumor 226 / other 20 |
| histology | yes | 226 | LUAD (series is stage I–II lung adenocarcinoma) |
| ICI response | **no** | 0 | treatment-naive surgical cohort |
| ESTIMATE published scores | **no** | 0 | computed here from Yoshihara lists |

Max-mean unique-mapped GPL570 probe (multi-mapped `///` probes dropped for the three named genes):

| gene | probe |
|---|---|
| CLDN4 | `201428_at` |
| CD8A | `205758_at` |
| TACSTD2 | `202286_s_at` |

ESTIMATE gene collapse uses first-symbol mapping (same as A1 extra) so Immune141 coverage is 138/141 and Stromal141 is 136/141.

## ESTIMATE

R `estimate` is not required. Scores follow the A1 extra public implementation:

- ssGSEA (Barbie/GSVA, τ=0.25) on Yoshihara 2013 Stromal141 + Immune141
- ranks scaled to 1…10000
- `ESTIMATEScore = StromalScore + ImmuneScore`
- `TumorPurity = cos(0.6049872018 + 0.0001467884 × ESTIMATEScore)`

On this matrix ESTIMATEScore ranges ~5.1k–13.6k. The cosine then falls below 0 for **208 / 226** tumors. Those values are **not** usable as a 0–1 purity fraction. They remain monotone with ESTIMATEScore here (angles still in (0, π)), so the rank residual on TumorPurity equals the rank residual on ESTIMATEScore. The [0, 1] subset (n=18) is only the least-impure tail and is **not** the analysis n.

**ImmuneScore after ESTIMATEScore / TumorPurity is collinear.** ImmuneScore is one addend of ESTIMATEScore. The CLDN4–ImmuneScore partial is reported because it was requested; it is not an independent purity control. The independent residual is **CLDN4 vs CD8A | ESTIMATEScore**.

CD8A vs ImmuneScore (positive-control): ρ=0.778 (p=4.43e-47, n=226).

## Main Spearman (n=226)

Partial = Pearson of rank residuals on ESTIMATE TumorPurity; df = n − 3.

| predictor | endpoint | n | ρ | p | partial ρ \| purity | partial p |
|---|---|---:|---:|---:|---:|---:|
| CLDN4 | CD8A | 226 | **-0.341** | 1.45e-07 | -0.127 | 0.0574 |
| CLDN4 | ImmuneScore | 226 | **-0.366** | 1.50e-08 | -0.031 | 0.639 |
| CLDN4 | TACSTD2 | 226 | **0.471** | 7.48e-14 | 0.392 | 1.15e-09 |
| TACSTD2 | CD8A | 226 | -0.289 | 9.76e-06 | -0.107 | 0.108 |
| TACSTD2 | ImmuneScore | 226 | -0.306 | 2.75e-06 | -0.023 | 0.735 |
| CLDN4 | StromalScore | 226 | -0.343 | 1.18e-07 | — | — |
| TACSTD2 | StromalScore | 226 | -0.304 | 3.27e-06 | — | — |
| CLDN4 | TumorPurity | 226 | 0.395 | 7.46e-10 | — | — |
| TACSTD2 | TumorPurity | 226 | 0.332 | 3.21e-07 | — | — |
| CLDN4 | TJ (no CLDN4) | 226 | 0.614 | 7.68e-25 | 0.512 | 1.86e-16 |

CLDN4 tracks TACSTD2 and a nine-gene TJ companion, and is anti-correlated with CD8A and ImmuneScore on the unadjusted tumor matrix. After the ESTIMATE impurity axis the CD8A residual is attenuated (p=0.057). That is extra East-Asian LUAD microarray weight for **CLDN4-high with TACSTD2-high / immune-low**, not an ICI-response test.

## What is not done

- No re-audit of the A9 3/3 TJ catalog (CLDN4 / MICALL2 / PARD6B recurrence).
- No OncoSG or CPTAC re-run.
- No ICI ORR / PFS model (labels are not deposited).
- No R `estimate` purity transform beyond the published cosine on the computed ESTIMATEScore.

## Files

- `analyze.py` — GEO download, tumor filter, ESTIMATE, Spearman, figures
- `estimate_yoshihara_2013.tsv` — Yoshihara 2013 Supp Data 1 (141+141)
- `tables/one_row.tsv`, `stats.tsv`, `per_sample.tsv`, `probe_used.tsv`, `coverage.tsv`, `summary.json`
- `figures/fig1_cldn4_vs_cd8a.png`
- `figures/fig2_cldn4_vs_immunescore.png`
- `figures/fig3_cldn4_vs_tacstd2.png`
- `figures/fig4_spearman_forest.png`
- `figures/fig5_correlation_heatmap.png`
- `figures/fig6_tacstd2_vs_immune.png`

```bash
python3 methods/gse31210_cldn4_immune/analyze.py
```
