# Finding — GSE72094 LUAD: CLDN4 vs CD274 / HLA (new CD274 cut)

**Additive public cohort.** Schabath / Moffitt **GSE72094** (PMID 26477306): resected **lung adenocarcinoma**, **GPL15048** Rosetta/Merck HuRSTA custom Affymetrix 2.0, author IRON-normalized series matrix.

**Honest n.** Series matrix has **442** arrays. All **442** have `source_name = lung adenocarcinoma`. Non-LUAD dropped: **0**. Complete-case n for CLDN4, CD274, and MHC-I genes is **442**. No ICI labels. No paired normals.

**This page is the new CD274 / HLA cut only.** Residual CLDN4 vs CD8A after ESTIMATE is already **not significant** on this matrix (PR 302: partial ρ = −0.051, p = 0.282, n = 442) and is **not** re-run or headlined here.

## Verdict (unadjusted Spearman; n=442)

| Test | n | ρ | 95% CI | p | BH q |
|---|---:|---:|---|---:|---:|
| **CLDN4 vs CD274** | 442 | -0.034 | [-0.128, +0.059] | 0.477 | 0.477 |
| **CLDN4 vs MHC-I** | 442 | -0.124 | [-0.216, -0.032] | 0.0088 | 0.0176 |
| CLDN4 vs HLA-A | 442 | -0.018 | [-0.115, +0.077] | 0.711 | — |
| CLDN4 vs HLA-B | 442 | -0.060 | [-0.153, +0.031] | 0.21 | — |
| CLDN4 vs HLA-C | 442 | -0.099 | [-0.190, -0.006] | 0.0369 | — |
| CLDN4 vs B2M | 442 | -0.241 | [-0.327, -0.149] | 3.06e-07 | — |

BH *q* is within the two new primary axes (CD274, MHC-I) only. Single HLA genes are the MHC-I decomposition, not a second family of claims. MHC-I = unweighted within-cohort z-mean of HLA-A / HLA-B / HLA-C / B2M (**4/4** present: HLA-A, HLA-B, HLA-C, B2M).

**What holds.** The new **CD274** cut is **null**: CLDN4 vs CD274 Spearman CI includes 0 (n=442, ρ=-0.034, p=0.477). CD274 Q4 vs Q1 does not separate CLDN4, and CLDN4 Q4 ∩ CD274 Q4 is chance (28/111 = 25.2%, OR=1.01, p=1). MHC-I is a weak inverse (ρ=-0.124, p=0.0088, q=0.0176); that cassette signal is **B2M-driven** (ρ=-0.241, p=3.06e-07). HLA-A and HLA-B are NS; HLA-C is weak.

**What does not hold.** CLDN4-high is not CD274-high on this LUAD microarray. Do not cite this page as a residual-vs-CD8 result.

## New CD274 cut (Q4 vs Q1)

CD274 high = ≥ cohort 75th percentile; CD274 low = ≤ 25th percentile. Ties at the percentile are kept, so arm n is reported rather than assumed n/4.

| Cut | n_high | n_low | CLDN4 median (high vs low) | Mann–Whitney p | CLDN4 Q4 ∩ CD274 Q4 | Fisher OR | Fisher p |
|---|---:|---:|---|---:|---|---:|---:|
| **CD274 Q4 vs Q1** | 111 | 111 | 10.499 vs 10.518 | 0.547 | 28 / 111 of CLDN4 Q4 (25.2%) | 1.01 | 1 |

Independence expectation for two Q4 calls is 25% of the CLDN4-Q4 arm. MHC-I Q4 overlap with CLDN4 Q4: 20 / 111 (18.0%), OR=0.58, p=0.0573.

TACSTD2 vs CD274 is a same-run companion only (n=442, ρ=-0.021, p=0.654) and is not the finding.

## What this is not

- **Not** residual CLDN4 vs CD8A. That test is already NS on this matrix (PR 302) and is not headlined.
- **Not** an ICI-response, PD-L1 IHC, or protein result. CD274 and HLA here are microarray RNA.
- **Not** a claim that CLDN4 induces PD-L1 or MHC-I. Association only.
- ESTIMATE residual for CD274 / MHC-I is in `correlations.tsv` as a sensitivity column. It is **not** the headline.

## Methods

- **Matrix:** GEO `GSE72094_series_matrix.txt.gz` (author-processed IRON + RNA-quality batch correction). Probe × sample values used as published.
- **Annotation:** GPL15048 `GeneSymbol`. First `///` token, then first whitespace token. **Max-mean** probe collapse to HUGO. HLA-A has no singleton symbol on this table; all five probes are `HLA-A LOC100507703` (Entrez 3105) and are mapped to HLA-A. One HLA-C probe is dual-annotated `HLA-C HLA-E` and is first-token mapped to HLA-C; two probes are exact `HLA-C`.
- **CD274** = collapsed `CD274` (probes `merck-NM_014143_at`, `merck-CK904742_at`, `merck2-ENST00000381577_at`; max-mean picks one).
- **MHC-I** = mean of within-cohort z-scores of HLA-A, HLA-B, HLA-C, B2M.
- **CD274 cut:** cohort quartiles on collapsed CD274. Primary contrast is Q4 vs Q1 (Mann–Whitney on CLDN4). Overlap is CLDN4 Q4 ∩ CD274 Q4 (two-sided Fisher).
- **Spearman** two-sided, complete cases. Bootstrap 95% CI: 5,000 resamples, seed 0, CLDN4 rows only.
- **ESTIMATE** (sensitivity column only): ssGSEA τ=0.25 on Yoshihara 2013 Stromal141 + Immune141. Coverage **141/141** stromal, **141/141** immune. Partial Spearman = Pearson of rank residuals; df = n − 3.

## Files

- `correlations.tsv` — n / ρ / p / CI / ESTIMATE residual (not headlined)
- `q4_cuts.tsv` / `q4_overlap.json` — new CD274 quartile cut
- `samples.tsv` — per-sample genes + MHC-I + CD274 quartile
- `coverage.tsv` / `provenance.json`

```bash
python3 -m pip install -r methods/gse72094_cldn4_cd274/requirements.txt
python3 methods/gse72094_cldn4_cd274/analyze.py
```

Place (or let the script download) `GSE72094_series_matrix.txt.gz` and `GPL15048_datatable.txt` under `data/gse72094/`.
