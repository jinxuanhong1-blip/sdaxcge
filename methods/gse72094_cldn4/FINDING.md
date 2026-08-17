# Finding — GSE72094 LUAD microarray: CLDN4 vs CD8A after ESTIMATE

**Additive public cohort.** Schabath / Moffitt **GSE72094** (PMID 26477306): resected **lung adenocarcinoma**, **GPL15048** Rosetta/Merck HuRSTA custom Affymetrix 2.0, author IRON-normalized series matrix. All **n=442** matrix samples are `source_name = lung adenocarcinoma`. No paired normals.

Primary question: **CLDN4 vs CD8A**, unadjusted and after ESTIMATE. **TACSTD2** is a companion on the same samples, not an audit of prior TACSTD2 residual claims.

## Verdict

| Test | n | Unadj ρ | Unadj p | After ESTIMATEScore | Partial p |
|---|---:|---:|---:|---:|---:|
| **CLDN4 vs CD8A** | 442 | -0.219 | 3.35e-06 | -0.051 | 0.282 |
| TACSTD2 vs CD8A (companion) | 442 | -0.227 | 1.39e-06 | -0.085 | 0.0741 |
| CLDN4 vs GEP18 | 442 | -0.194 | 4.17e-05 | +0.023 | 0.637 |
| TACSTD2 vs GEP18 (companion) | 442 | -0.216 | 4.47e-06 | -0.046 | 0.331 |
| CLDN4 vs TACSTD2 | 442 | +0.549 | 3.33e-36 | +0.521 | 4.34e-32 |

**What holds.** Unadjusted, higher **CLDN4** tracks **lower CD8A** (n=442, ρ=-0.219, p=3.35e-06). Same sign for TACSTD2 vs CD8A (companion). CLDN4 and TACSTD2 are positively coexpressed.

**What does not hold.** After ESTIMATEScore residual, **CLDN4 vs CD8A is not significant** (partial ρ=-0.051, p=0.282). The unadjusted inverse association is largely the stromal/immune-content axis.

## Primary numbers

| Test | Covariate | n | ρ | p |
|---|---|---:|---:|---:|
| CLDN4 vs CD8A | none | 442 | -0.219 | 3.35e-06 |
| CLDN4 vs CD8A | ESTIMATEScore | 442 | -0.051 | 0.282 |
| CLDN4 vs CD8A | TumorPurity (cosine) | 442 | -0.051 | 0.282 |
| CLDN4 vs CD8A | ImmuneScore | 442 | -0.021 | 0.664 |
| TACSTD2 vs CD8A | none | 442 | -0.227 | 1.39e-06 |
| TACSTD2 vs CD8A | ESTIMATEScore | 442 | -0.085 | 0.0741 |
| TACSTD2 vs CD8A | TumorPurity (cosine) | 442 | -0.085 | 0.0741 |

Context (unadjusted vs ESTIMATE axes):

| Gene | vs ESTIMATEScore | vs ImmuneScore |
|---|---|---|
| CLDN4 | n=442 ρ=-0.252 p=7.43e-08 | n=442 ρ=-0.253 p=6.60e-08 |
| TACSTD2 | n=442 ρ=-0.233 p=7.35e-07 | n=442 ρ=-0.228 p=1.33e-06 |
| CD8A | n=442 ρ=+0.734 p=6.61e-76 | n=442 ρ=+0.819 p=4.05e-108 |

## Methods

- **Matrix:** GEO `GSE72094_series_matrix.txt.gz` (author-processed IRON + RNA-quality batch correction). Probe × sample values used as published.
- **Annotation:** GPL15048 `GeneSymbol` (NCBI platform table). First symbol if `///`. **Max-mean** probe collapse to HUGO.
- **CD8** = `CD8A`. **GEP18** = unweighted within-cohort z-mean of Ayers 2017 18-gene list (18/18 present; missing: none).
- **ESTIMATE:** ssGSEA (Barbie/GSVA τ=0.25, ranks scaled 1…10000) on Yoshihara 2013 Stromal141 + Immune141 (`data/signatures/estimate_yoshihara_2013.tsv`). Coverage **140/141** stromal, **139/141** immune. `ESTIMATEScore = Stromal + Immune`. `TumorPurity = cos(0.6049872018 + 0.0001467884 × ESTIMATEScore)`.
- **Partial Spearman:** Pearson of rank residuals; df = n − 3. Primary residual covariate = **ESTIMATEScore** (impurity axis). TumorPurity and ImmuneScore are sensitivity rows.
- **ESTIMATE scale.** ssGSEA ESTIMATEScore range 1457–14960. Yoshihara cosine TumorPurity is −0.943 to +0.683 (outside [0,1] because ssGSEA scores are on a different scale than the original ESTIMATE R package). No sample exceeds the cosine wrap threshold (~17281). TumorPurity vs ESTIMATEScore is a perfect rank anti-correlation (ρ=−1, n=442), so the TumorPurity residual equals the ESTIMATEScore residual here.

## Cohort

442 resected LUAD (Moffitt). Platform GPL15048. Genes after collapse: 22115. CLDN4, TACSTD2, and CD8A are all present.

This write-up does **not** re-open or audit prior TACSTD2 residual claims from other PRs. TACSTD2 numbers here are a same-run companion only.

## Files

- `correlations.tsv` — all n / ρ / p
- `samples.tsv` — per-sample genes + ESTIMATE + mutation annotations
- `coverage.tsv` / `provenance.json`
- Reproduce: `python scripts/gse72094_cldn4.py` after placing the GEO matrix and GPL15048 table under `data/gse72094/`
