# FINDING — GSE19804 paired LUAD: CLDN4 Q4 vs Q1 vs CD8

**Additive only. Q4 extra.** Continuous CLDN4–CD8 Spearman and the
ESTIMATE-purity partial on this same 60-tumor table are **already known**
(PR #235 extra LUAD RNA; not a new claim). This folder only recuts tumor
CLDN4 as **Q4 vs Q1** and tests **CD8A**.

Public paired tumor / adjacent-normal arrays from never-smoking Taiwanese
women (Lu et al., *Cancer Epidemiol Biomarkers Prev* 2010, PMID 20802022;
GEO [GSE19804](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19804);
GDS3837; GPL570 U133 Plus 2.0). GDS labels the series NSCLC; the extra
LUAD table used the 60 tumors. Unit here is the **tumor array**. Paired
normals are dropped. No ICI arm.

## Honest n

| item | n | rule |
|---|---:|---|
| arrays on the GEO series | 120 | GDS3837 sample count |
| paired adjacent normals | 60 | dropped; not on the harvest |
| tumors on the harvested table | **60** | same mask as PR #235 |
| CLDN4 / CD8A non-NA | 60 / 60 | complete-case |
| Q1 / Q2 / Q3 / Q4 | 15 / 15 / 15 / 15 | `qcut(rank(method="first"), 4)` on tumor CLDN4 |
| **Q4 vs Q1 used** | **15 vs 15** | MWU n |

Do **not** write n=120 or n=60 for the quartile test. The continuous
purity-partial keeps n=60.

## Q4 vs Q1 vs CD8 (this extra)

Quartiles on raw tumor CLDN4. Primary endpoint is CD8A. Rank-biserial
\(r = 2U/(n_4 n_1)-1\) (positive = Q4 higher). 95% CI = percentile
bootstrap, 2000 resamples, seed 20260817.

| endpoint | n Q4 / Q1 | median Q4 | median Q1 | Δ (Q4−Q1) | r | 95% CI | p |
|---|---:|---:|---:|---:|---:|---|---:|
| **CD8A** | **15 / 15** | 7.56 | 8.93 | **−1.38** | **−0.72** | −0.95 to −0.43 | **7.8×10⁻⁴** |
| CD8A residual \| purity | 15 / 15 | −3.99 | +2.28 | −6.27 | −0.32 | −0.70 to +0.08 | 0.15 |

CLDN4 cuts among the 60 tumors: Q1 ≤ 8.42, Q4 ≥ 9.61 (harvested max-mean
GPL570 scale).

**What holds:** CLDN4-high tumors (Q4) have lower CD8A than Q1 at
**15 vs 15**. Same sign as the already-known continuous Spearman.

**What does not hold:** The same Q4 vs Q1 split on the **purity residual**
of CD8A is the same sign and **not significant** (p=0.15). Throwing away
Q2+Q3 (30 tumors) underpowers the residual contrast. Do not upgrade the
known n=60 partial to a residual Q4-vs-Q1 claim.

## Already known continuous (not re-claimed)

Recomputed on the harvested table; ρ matches PR #235 exactly. Partial p
in PR #235 uses the t / df = n−3 form (0.0356); Pearson p on the same
residuals is 0.0341.

| contrast | n | ρ | p | source |
|---|---:|---:|---:|---|
| CLDN4–CD8A Spearman | 60 | **−0.465** | 1.85×10⁻⁴ | PR #235 |
| CLDN4–CD8A partial \| ESTIMATE TumorPurity | 60 | **−0.274** | **0.0356** | PR #235 |
| CLDN4–GEP18 Spearman | 60 | −0.287 | 0.026 | PR #235 |
| CLDN4–GEP18 partial \| purity | 60 | +0.083 | 0.53 | PR #235 (null) |

Purity here is ESTIMATE TumorPurity (Yoshihara 141+141 ssGSEA, then the
published cosine). CD8A vs TumorPurity on this table is strong
(ρ=−0.755, p=3.3×10⁻¹², n=60), so the unadjusted CLDN4–CD8A cloud is
partly impurity. The partial is the purity-controlled continuous number
that is already on the extra LUAD table.

## Same-table extras (not requested)

Same CLDN4 Q4 vs Q1 arms (15 vs 15):

| endpoint | Δ (Q4−Q1) | r | p |
|---|---:|---:|---:|
| ESTIMATE ImmuneScore | −1585 | −0.50 | 0.020 |
| GEP18 | −0.460 | −0.51 | 0.018 |

These are the same 15 vs 15 tumors. GEP18 Q4 is negative while the
already-known GEP18 **partial** is null — do not treat GEP18 Q4 as
purity-controlled.

## Methods (this slice)

- Harvest only: `harvested/samples_GSE19804.tsv` from PR #235
  (`results/rework/A1_extra_luad/samples_GSE19804.tsv`). No new GEO
  download. No invented signature.
- Predictor: tumor CLDN4. No TACSTD2 gate.
- Quartiles: `pd.qcut(rank(method="first"), 4)` on the 60 tumors.
- Primary test: two-sided Mann–Whitney U, Q4 vs Q1, on CD8A.
- Supporting residual: CD8A rank residual on ESTIMATE TumorPurity,
  then the same Q4 vs Q1 arms (quartiles stay on raw CLDN4).
- Continuous Spearman / partial are cited from PR #235 and recomputed
  only as a checksum.

## What this does not claim

- It does not re-open the continuous purity-partial (already known).
- It does not test ICI response, PFS, or OS (labels are not deposited).
- It does not use the 60 paired normals as a CD8 contrast.
- It does not write n=120 or n=60 for the MWU.
- It does not claim that residual Q4 vs Q1 is significant.

## Files

- `harvested/samples_GSE19804.tsv`, `harvested/PROVENANCE.md`
- `tables/n_table.tsv` — honest n
- `tables/q4_vs_q1.tsv` — MWU Q4 vs Q1
- `tables/continuous_known.tsv` — checksum of the already-known continuous
- `tables/sample_scores.tsv` — per-tumor CLDN4, quartile, CD8A, residual
- `tables/summary.json`
- `figures/fig_q4_vs_q1_cd8.png`

```bash
python3 methods/gse19804_cldn4_q4/analyze.py
```
