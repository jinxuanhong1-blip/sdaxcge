# GSE30219 — NSCLC bulk CLDN4 vs CD8A / CD274

**Additive only.** No prior GSE30219 CLDN4–immune table is re-audited. This folder only measures **CLDN4 vs CD8A and CD274** on the public Rousseaux mixed-histology lung series.

Public Centre Léon Bérard surgical lung series (Rousseaux et al., *Sci Transl Med* 2013, PMID 23698379; GEO [GSE30219](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30219); GPL570 Affymetrix U133 Plus 2.0). Unit is the **array**. No ICI arm. No slide was re-scored.

## Honest n

GEO histology codes on this series are **not** WHO abbreviations. **`SCC` is small-cell**, not squamous. Squamous is **`SQC`**. Do not write n=293 for NSCLC.

| item | n | rule |
|---|---:|---|
| arrays on the series matrix | 307 | GEO `GSE30219_series_matrix.txt.gz` |
| GEO histology inventory | 307 | ADC 85, SQC 61, LCNE 56, BAS 39, CARCI 24, SCC 21, NTL 14, OTHER 4, LCC 3 |
| non-tumoral lung (`NTL`) | 14 | titles “Non tumoral lung”; dropped. GEO `tissue` still says lung tumour |
| all tumors (any histology) | 293 | matrix minus `NTL` = paper 293 |
| SCLC (`SCC`) | 21 | **not NSCLC**; held out |
| carcinoid (`CARCI`) | 24 | **not NSCLC**; held out |
| LCNEC (`LCNE`) | 56 | neuroendocrine; **not** in primary NSCLC n |
| other / unmapped histology | 4 | GEO `Other` unspecified tumors; held out |
| ADC | 85 | GEO `ADC` |
| SQC (squamous) | 61 | GEO `SQC` (not `SCC`) |
| LCC | 3 | GEO `LCC` |
| BAS (basaloid) | 39 | GEO `BAS` |
| **primary NSCLC** | **188** | ADC + SQC + LCC + BAS |
| CLDN4 / CD8A / CD274 finite | 188 | unique-mapped max-mean collapse |
| ESTIMATE ImmuneScore | 188 | Yoshihara Immune141 ssGSEA (138/141) |

Primary tests use **n=188** NSCLC arrays. LCNE is a sensitivity row only. Paper text “293 tumors” includes SCLC / carcinoid / LCNE and is **not** the NSCLC n.

## One-row table

| dataset | histology rule | n | CLDN4–CD8A ρ (p) | CLDN4–CD274 ρ (p) | CLDN4–CD8A partial \| purity (p) | CLDN4–CD274 partial \| purity (p) |
|---|---|---:|---|---|---|---|
| GSE30219 NSCLC | ADC+SQC+LCC+BAS | 188 | -0.074 (0.311) | -0.055 (0.45) | -0.016 (0.826) | -0.011 (0.877) |
| GSE30219 ADC | ADC only | 85 | -0.345 (0.00122) | -0.103 (0.349) | -0.141 (0.2) | 0.107 (0.334) |
| GSE30219 SQC | SQC only | 61 | -0.157 (0.225) | -0.134 (0.303) | 0.054 (0.682) | 0.015 (0.911) |
| GSE30219 +LCNE (sensitivity) | NSCLC+LCNE | 244 | -0.113 (0.077) | -0.039 (0.542) | -0.107 (0.0971) | -0.006 (0.931) |

Full numeric rows: `tables/one_row.tsv`, `tables/stats.tsv`.

## Matrix and labels (nothing invented)

GEO series-matrix MAS5 linear intensity. Tests use native ranks (Spearman / ssGSEA), so log vs linear does not change ρ. Compact extracts store **log2(MAS5+1)** for plots.

| field | public? | n | what is there |
|---|---|---:|---|
| histology | yes | 307 | GEO characteristic; codes listed above |
| ICI response | **no** | 0 | surgical / pretreatment diagnostic series |
| OS / DFS | yes | 307 | GEO follow-up / status / DFS / relapse; **not used** (not a survival claim) |
| ESTIMATE published scores | **no** | 0 | computed here from Yoshihara lists |

Max-mean unique-mapped GPL570 probe (multi-mapped `///` probes dropped for the three named genes):

| gene | probe |
|---|---|
| CLDN4 | `201428_at` |
| CD8A | `205758_at` |
| CD274 | `227458_at` |

ESTIMATE gene collapse uses first-symbol mapping (same as A1 extra / GSE31210) so Immune141 coverage is 138/141 and Stromal141 is 136/141.

## ESTIMATE

R `estimate` is not required. Scores follow the A1 extra public implementation:

- ssGSEA (Barbie/GSVA, τ=0.25) on Yoshihara 2013 Stromal141 + Immune141
- ranks scaled to 1…10000
- `ESTIMATEScore = StromalScore + ImmuneScore`
- `TumorPurity = cos(0.6049872018 + 0.0001467884 × ESTIMATEScore)`

On the NSCLC matrix ESTIMATEScore ranges -1830–13699. The cosine falls outside [0, 1] for **135 / 188** tumors. Those values are **not** usable as a 0–1 purity fraction. They remain monotone with ESTIMATEScore here, so the rank residual on TumorPurity equals the rank residual on ESTIMATEScore. The [0, 1] subset (n=53) is **not** the analysis n.

**ImmuneScore after ESTIMATEScore / TumorPurity is collinear.** The independent residual is **CLDN4 vs CD8A | ESTIMATEScore** (equals TumorPurity residual here).

CD8A vs ImmuneScore (positive-control, NSCLC): ρ=0.825 (p=6.87e-48, n=188).

## Main Spearman (primary NSCLC n=188)

Partial = Pearson of rank residuals on ESTIMATE TumorPurity; df = n − 3. Bootstrap 95% CI, 2,000 resamples, seed `20260817`.

| predictor | endpoint | n | ρ | 95% CI | p | partial ρ \| purity | partial p |
|---|---|---:|---:|---|---:|---:|---:|
| CLDN4 | CD8A | 188 | **-0.074** | -0.213 to 0.071 | 0.311 | -0.016 | 0.826 |
| CLDN4 | CD274 | 188 | **-0.055** | -0.207 to 0.104 | 0.45 | -0.011 | 0.877 |
| CLDN4 | ImmuneScore | 188 | -0.042 | -0.176 to 0.099 | 0.571 | 0.141 | 0.0543 |
| CLDN4 | TumorPurity | 188 | 0.084 | -0.060 to 0.217 | 0.252 | — | — |
| CD274 | CD8A | 188 | 0.578 | 0.450 to 0.691 | 3.71e-18 | 0.299 | 3.13e-05 |

Mixed NSCLC (ADC+SQC+LCC+BAS) is **null** for both named pairs. The ADC-only unadjusted CLDN4–CD8A inverse (see one-row table) attenuates after the ESTIMATE impurity axis, the same pattern as GSE31210 LUAD. Do not treat the mixed-NSCLC pool as evidence that CLDN4-high NSCLC is CD8-low or PD-L1-low. It is **not** an ICI-response test.

## What is not done

- No re-use of the paper’s 293-tumor “all histologies” n as NSCLC.
- No SCLC (`SCC`) or carcinoid (`CARCI`) in the primary n.
- No ICI ORR / PFS model (labels are not deposited).
- No R `estimate` purity transform beyond the published cosine on the computed ESTIMATEScore.

## Files

- `analyze.py` — GEO download, histology inventory, ESTIMATE, Spearman, figures
- `estimate_yoshihara_2013.tsv` — Yoshihara 2013 Supp Data 1 (141+141)
- `tables/one_row.tsv`, `stats.tsv`, `per_sample.tsv`, `probe_used.tsv`, `coverage.tsv`, `histology_inventory.tsv`, `summary.json`
- `figures/fig1_cldn4_vs_cd8a.png`
- `figures/fig2_cldn4_vs_cd274.png`
- `figures/fig3_cldn4_vs_immunescore.png`
- `figures/fig4_spearman_forest.png`

```bash
python3 methods/gse30219_cldn4/analyze.py
```
