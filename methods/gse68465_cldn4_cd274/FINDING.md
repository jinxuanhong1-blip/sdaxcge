# GSE68465 LUAD — CLDN4 vs CD274 and HLA-A/B/C (partial on ESTIMATE)

**Additive only. New CD274 / classical MHC-I cut.** CLDN4 vs CD8A after ESTIMATE TumorPurity is **already known** ([PR 298](https://github.com/jinxuanhong1-blip/sdaxcge/pull/298): n=443, unadj ρ=−0.177, partial ρ=−0.092, p=0.054) and is **not** re-cut. This folder only measures CLDN4 against **CD274** and **HLA-A / HLA-B / HLA-C**.

Public Director's Challenge LUAD microarray (Shedden et al., *Nat Med* 2008; GEO [GSE68465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68465); GPL96 Affymetrix Human Genome U133A). Unit is the **tumor array**. No ICI arm. No slide was re-scored.

## Honest n

| item | n | rule |
|---|---:|---|
| arrays on the series matrix | 462 | GEO `GSE68465_series_matrix.txt.gz` (22,283 probes) |
| LUAD tumors | **443** | `disease_state` contains `Adenocarcinoma` |
| Normal arrays | 19 | dropped |
| CLDN4 finite (`201428_at`) | 443 | official GPL96 = CLDN4 / Entrez 1364 |
| **CD274 finite (Entrez 29126)** | **0** | no GPL96 probe; Plus-2 `223834_at` / `227458_at` are not on U133A |
| HLA-A / HLA-B / HLA-C finite | 443 | first-symbol max-mean; 3/3 present |
| ESTIMATEScore | 443 | Yoshihara Stromal141+Immune141 ssGSEA (136/141 + 138/141) |
| ESTIMATE TumorPurity in [0, 1] | 442 | cosine map; 1 tumor wraps outside [0, 1] |
| ICI labels | 0 | surgical / multi-site prognostic series |
| **Primary pairwise n (CLDN4 + CD274)** | **0** | ABSENT — do not invent a PD-L1 ρ |
| **Primary pairwise n (CLDN4 + HLA-A/B/C)** | **443** | this is the n used for HLA / MHC-I |

Primary HLA tests use **n=443** tumors. That is the same tumor rule as the sibling CD8 cut ([PR 298](https://github.com/jinxuanhong1-blip/sdaxcge/pull/298)); it is not a new cohort hunt.

PDCD1LG2 (PD-L2, `220049_s_at`, Entrez 80380) **is** on U133A. It is not CD274 and is **not** substituted.

## One-row table

| dataset | n tumors | CLDN4 | CD274 | HLA-A/B/C | CLDN4–CD274 ρ (p) | CLDN4–HLA-A ρ (p) | CLDN4–HLA-B ρ (p) | CLDN4–HLA-C ρ (p) | CLDN4–CD274 partial \| ESTIMATE (p) | CLDN4–HLA-A partial \| ESTIMATE (p) | CLDN4–HLA-B partial \| ESTIMATE (p) | CLDN4–HLA-C partial \| ESTIMATE (p) |
|---|---:|---|---|---|---|---|---|---|---|---|---|---|
| GSE68465 LUAD tumors | 443 | `201428_at` | **ABSENT** (Entrez 29126 not on U133A) | 3/3 | **ABSENT (n=0)** | −0.040 (0.397) | −0.025 (0.605) | −0.029 (0.536) | **ABSENT (n=0)** | +0.006 (0.897) | +0.076 (0.112) | +0.051 (0.287) |

MHC-I mean-z (HLA-A/B/C only): unadj −0.040 (0.397); partial \| ESTIMATE +0.036 (0.452).

Full numeric rows: `tables/one_row.tsv`, `tables/stats.tsv`.

## New CD274 cut

MAS5 as deposited. Tests use native ranks (Spearman / ssGSEA), so log vs linear does not change ρ. Bootstrap 95% CI, 2,000 resamples, seed `20260817`. MHC-I mean-z = gene-wise z of **HLA-A/B/C only** (3/3; not B2M/TAP). Partial = Pearson of rank residuals on ESTIMATEScore; df = n − 3. On this matrix TumorPurity is a monotone cosine of ESTIMATEScore (when the angle stays in (0, π)), so the rank residual on TumorPurity equals the rank residual on ESTIMATEScore for tumors with a defined cosine.

| pair | n | ρ | 95% CI | p | ρ \| ESTIMATE (p) | 95% CI partial | verdict |
|---|---:|---:|---|---:|---|---|---|
| CLDN4 vs **CD274** | **0** | — | — | — | — | — | **ABSENT** |
| CLDN4 vs HLA-A | **443** | −0.040 | −0.132 to +0.057 | 0.397 | +0.006 (0.897) | −0.089 to +0.102 | **NO_EVIDENCE** |
| CLDN4 vs HLA-B † | **443** | −0.025 | −0.115 to +0.068 | 0.605 | +0.076 (0.112) | −0.016 to +0.168 | **NO_EVIDENCE** |
| CLDN4 vs HLA-C | **443** | −0.029 | −0.121 to +0.069 | 0.536 | +0.051 (0.287) | −0.045 to +0.149 | **NO_EVIDENCE** |
| CLDN4 vs MHC-I (HLA-A/B/C mean-z) | **443** | −0.040 | −0.133 to +0.054 | 0.397 | +0.036 (0.452) | −0.057 to +0.133 | **NO_EVIDENCE** |

† **HLA-B is in Yoshihara Immune141.** Residualizing ESTIMATEScore (stromal + immune) out of HLA-B is partly circular. HLA-A and HLA-C are not in Immune141. The HLA-B row is reported because it was requested; it is not an independent infiltrate control.

CLDN4 Q4 vs Q1 on the same endpoints:

| endpoint | n Q4 vs Q1 | rank-biserial (Q4−Q1) | MWU p | verdict |
|---|---|---:|---:|---|
| CD274 | — | — | — | **ABSENT** |
| HLA-A | 111 vs 111 | −0.077 | 0.323 | NO_EVIDENCE |
| HLA-B | 111 vs 111 | −0.052 | 0.508 | NO_EVIDENCE |
| HLA-C | 111 vs 111 | −0.061 | 0.431 | NO_EVIDENCE |
| MHC1_HLAabc | 111 vs 111 | −0.083 | 0.285 | NO_EVIDENCE |

There is no public CD274 (PD-L1) measurement on this U133A series. Classical MHC-I is present at honest n=443. Antigen-presentation / PD-L1 cut only; not an ICI-response test. CLDN4–CD8A is not restated.

## ESTIMATE

R `estimate` is not required. Scores follow the sibling public implementation ([PR 298](https://github.com/jinxuanhong1-blip/sdaxcge/pull/298)):

- ssGSEA (Barbie/GSVA, τ=0.25) on Yoshihara 2013 Stromal141 + Immune141
- ranks scaled to 1…10000
- `ESTIMATEScore = StromalScore + ImmuneScore`
- `TumorPurity = cos(0.6049872018 + 0.0001467884 × ESTIMATEScore)`

**Partial covariate is ESTIMATEScore**, not ImmuneScore. ImmuneScore is one addend of ESTIMATEScore and is collinear with it. The cosine is below 0 for **1 / 443** tumors and is **not** used as a 0–1 fraction. The [0, 1] subset (n=442) is only the least-impure tail and is **not** the analysis n.

CLDN4 vs ImmuneScore is reported only as a context row (same matrix as the sibling cut). It is not a new CD8 test.

| predictor | endpoint | n | ρ | p | partial ρ \| ESTIMATE | partial p |
|---|---|---:|---:|---:|---:|---:|
| CLDN4 | ImmuneScore | 443 | −0.209 | 9.39e-06 | −0.091 | 0.056 |

## Given CD8 residual (not re-run)

From [PR 298](https://github.com/jinxuanhong1-blip/sdaxcge/pull/298), GSE68465 CLDN4 vs CD8A after the same ESTIMATE TumorPurity residual is **partial ρ = −0.092, p = 0.054, n = 443** (unadj −0.177, p = 1.7e-4). That CD8 table is not re-downloaded or re-fit here.

## Named-gene collapse

First-symbol max-mean (same as the sibling CD8 cut). Exact first-symbol / Entrez only — no substring aliases (`PDL1` hits `SPDL1`).

| gene | probe | rule |
|---|---|---|
| CLDN4 | `201428_at` | first-symbol max-mean |
| CD274 | **ABSENT** | no GPL96 hit for CD274 / PDCD1LG1 / Entrez 29126 |
| HLA-A | `215313_x_at` | first-symbol max-mean |
| HLA-B | `209140_x_at` | first-symbol max-mean |
| HLA-C | `216526_x_at` | first-symbol max-mean |

ESTIMATE gene collapse uses the same first-symbol mapping.

## What is not done

- No re-cut of CLDN4 vs CD8A / TACSTD2 (sibling folder / PR 298 / PR 235).
- No PD-L1 RNA claim (CD274 is not on the chip).
- No PD-L2 (PDCD1LG2) substitute for CD274.
- No ICI ORR / PFS model (labels are not deposited).
- No R `estimate` purity transform beyond the published cosine on the computed ESTIMATEScore.

## Files

- `analyze.py` — GEO download, tumor filter, ESTIMATE, Spearman, figures
- `estimate_yoshihara_2013.tsv` — Yoshihara 2013 Supp Data 1 (141+141)
- `tables/one_row.tsv`, `stats.tsv`, `per_sample.tsv`, `probe_confirm.tsv`, `coverage.tsv`, `highlow_cldn4_cd274_hla.tsv`, `summary.json`
- `figures/fig1_cldn4_vs_cd274.png` — ABSENT call
- `figures/fig2_cldn4_vs_hla.png`
- `figures/fig3_spearman_forest.png`
- `figures/fig4_correlation_heatmap.png`

```bash
python3 -m pip install -r methods/gse68465_cldn4_cd274/requirements.txt
python3 methods/gse68465_cldn4_cd274/analyze.py
```

Downloads (not committed) go to `$GSE68465_CLDN4_CD274_DATA` (default `/tmp/gse68465_cldn4_cd274`).
