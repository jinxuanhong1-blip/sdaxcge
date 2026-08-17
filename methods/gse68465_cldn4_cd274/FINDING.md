# GSE68465 LUAD — CLDN4 vs CD274 and HLA-A/B/C (partial on ESTIMATE)

**Additive only. New CD274 / classical MHC-I cut.** CLDN4 vs CD8A after ESTIMATE TumorPurity is **already known** ([PR 298](https://github.com/jinxuanhong1-blip/sdaxcge/pull/298): n=443, unadj ρ=−0.177, partial ρ=−0.092, p=0.054) and is **not** re-cut. This folder only measures CLDN4 against **CD274** and **HLA-A / HLA-B / HLA-C**.

Public Director's Challenge LUAD microarray (Shedden et al., *Nat Med* 2008; GEO [GSE68465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68465); GPL96 Affymetrix Human Genome U133A). Unit is the **tumor array**. No ICI arm.

## Honest n

| item | n | rule |
|---|---:|---|
| arrays on the series matrix | 462 | GEO `GSE68465_series_matrix.txt.gz` (22,283 probes) |
| LUAD tumors | **443** | `disease_state` contains `Adenocarcinoma` |
| Normal arrays | 19 | dropped |
| CLDN4 finite (`201428_at`) | 443 | official GPL96 = CLDN4 / Entrez 1364 |
| **CD274 finite (Entrez 29126)** | **0** | no GPL96 probe; Plus-2 `223834_at` / `227458_at` are not on U133A |
| HLA-A / HLA-B / HLA-C | 443 | first-symbol max-mean; 3/3 present on GPL96 |
| ICI labels | 0 | surgical / multi-site prognostic series |
| **Primary pairwise n (CLDN4 + CD274)** | **0** | ABSENT — do not invent a PD-L1 ρ |
| **Primary pairwise n (CLDN4 + HLA-A/B/C)** | **443** | HLA / MHC-I tests |

PDCD1LG2 (PD-L2, `220049_s_at`, Entrez 80380) **is** on U133A. It is not CD274 and is **not** substituted.

## One-row table

| dataset | n tumors | CLDN4 | CD274 | HLA-A/B/C | CLDN4–CD274 ρ (p) | CLDN4–HLA-A ρ (p) | CLDN4–HLA-B ρ (p) | CLDN4–HLA-C ρ (p) | CLDN4–CD274 partial \| ESTIMATE (p) | CLDN4–HLA-A partial \| ESTIMATE (p) | CLDN4–HLA-B partial \| ESTIMATE (p) | CLDN4–HLA-C partial \| ESTIMATE (p) |
|---|---:|---|---|---|---|---|---|---|---|---|---|---|
| GSE68465 LUAD tumors | 443 | `201428_at` | **ABSENT** (Entrez 29126 not on U133A) | 3/3 | **ABSENT (n=0)** | *run `analyze.py`* | *run `analyze.py`* | *run `analyze.py`* | **ABSENT (n=0)** | *run `analyze.py`* | *run `analyze.py`* | *run `analyze.py`* |

HLA Spearman / ESTIMATE partial numbers are written by `analyze.py` into this file and `tables/one_row.tsv`. CD274 remains **ABSENT** after the run.

## Given CD8 residual (not re-run)

From [PR 298](https://github.com/jinxuanhong1-blip/sdaxcge/pull/298), GSE68465 CLDN4 vs CD8A after ESTIMATE TumorPurity is **partial ρ = −0.092, p = 0.054, n = 443**. That CD8 table is not re-fit here.

## Reproduce

```bash
python3 -m pip install -r methods/gse68465_cldn4_cd274/requirements.txt
python3 methods/gse68465_cldn4_cd274/analyze.py
```
