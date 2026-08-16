# Hugo 2016 melanoma ICI: TACSTD2 and CLDN4 versus response

## Bottom line

Pretreatment tumor **TACSTD2 and CLDN4 do not separate anti-PD-1 responders
from nonresponders** in this public melanoma cohort. Both genes sit near
background in almost every sample. The only large values are a handful of
outliers, one of them a scalp subcutaneous biopsy (Pt28) with TACSTD2 FPKM
644 and CLDN4 FPKM 148.

This is a **null result**, not support for claim B5 (CLDN4-high worse ICI
response, pooled OR = 0.42). The CLDN4 point estimate here is in the
**opposite** direction (median-split response OR = 1.78) and the confidence
interval runs from 0.40 to 8.0. Melanoma is the wrong histology for these
epithelial genes; the cohort cannot test a TROP2/claudin barrier hypothesis.

## Cohort and endpoint

- Source: Hugo et al., *Cell* 2016, GEO
  [GSE78220](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78220).
- Treatment: pembrolizumab (anti-PD-1) in advanced melanoma.
- Public RNA: Cuffnorm FPKM, 25,268 genes × 28 biopsies from 27 patients.
  Every GEO sample is annotated `tissue: Melanoma biopsies`.
- Genes: **TACSTD2** (TROP2) and **CLDN4**.
- Response: GEO `Complete Response` / `Partial Response` / `Progressive
  Disease`. There is no SD in this series. Binary rule: CR+PR = responder,
  PD = nonresponder.
- Primary set: pretreatment, one row per patient.
  - Pt16 is on-treatment (PD) and is excluded.
  - Pt27 has two pretreatment CR biopsies (A/B); FPKM is averaged.
  - n = 26 patients (14 responders: 4 CR + 10 PR; 12 PD).
- Expression scale: `log2(FPKM + 1)`. Rank tests are unchanged by this
  monotonic transform.

## Results

Medians and means below are on the log2(FPKM + 1) scale. AUC is the
probability that a random responder has higher expression than a random
nonresponder (ties split by the Mann–Whitney midrank).

| Gene | Analysis | R / NR | Median R | Median NR | AUC (R higher) | 95% bootstrap CI | MW p |
|---|---|---:|---:|---:|---:|---:|---:|
| TACSTD2 | Pretreatment (primary) | 14 / 12 | 0.507 | 0.729 | 0.351 | 0.149–0.577 | 0.208 |
| TACSTD2 | Drop Pt28 | 13 / 12 | 0.494 | 0.729 | 0.301 | 0.109–0.526 | 0.097 |
| CLDN4 | Pretreatment (primary) | 14 / 12 | 0.477 | 0.369 | 0.542 | 0.310–0.762 | 0.738 |
| CLDN4 | Drop Pt28 | 13 / 12 | 0.444 | 0.369 | 0.506 | 0.276–0.737 | 0.978 |

On the raw FPKM scale the primary medians are TACSTD2 0.42 (R) vs 0.66 (NR)
and CLDN4 0.39 (R) vs 0.29 (NR). Those are background counts, not a
biomarker range.

Median-split (cohort median FPKM; high = above the median) response odds
ratio, Haldane–Anscombe 0.5-corrected Wald CI, two-sided Fisher p:

| Gene | Analysis | High R / High NR | Low R / Low NR | OR (high vs low) | 95% CI | Fisher p |
|---|---|---:|---:|---:|---:|---:|
| TACSTD2 | Primary | 6 / 7 | 8 / 5 | 0.56 | 0.12–2.52 | 0.695 |
| CLDN4 | Primary | 8 / 5 | 6 / 7 | 1.78 | 0.40–8.00 | 0.695 |

Detection (FPKM > 1) is equally uninformative: TACSTD2 3/14 vs 4/12
(Fisher p = 0.665); CLDN4 4/14 vs 3/12 (p = 1.0).

Sample-level pretreatment (both Pt27 biopsies kept) and adding Pt16 back
do not change the conclusion. Full numbers are in `summary.csv`.

TACSTD2 and CLDN4 are correlated (Spearman ρ = 0.53, p = 0.005; ρ = 0.47
after dropping Pt28). That is the expected epithelial co-expression, not
evidence of a response effect.

A three-level Kruskal–Wallis on TACSTD2 after dropping Pt28 is p = 0.044.
That comparison is an uncorrected sensitivity check in n = 25 with only
4 CR. It is driven by PR samples sitting even closer to zero (median FPKM
0.28) than CR (0.67) or PD (0.67). It is not a confirmatory finding and
it does not appear in the primary binary test.

## What the values actually look like

Primary-set FPKM percentiles:

| Gene | 25% | 50% | 75% | 90% |
|---|---:|---:|---:|---:|
| TACSTD2 | 0.29 | 0.49 | 1.12 | 7.21 |
| CLDN4 | 0.16 | 0.34 | 1.06 | 5.12 |

Only three patients exceed FPKM 5 for either gene:

| Patient | RECIST | Site | TACSTD2 | CLDN4 |
|---|---|---|---:|---:|
| Pt28 | PR | Scalp, SC | 644.3 | 147.9 |
| Pt29 | PD | Lung | 12.0 | 15.4 |
| Pt9 | CR | L forearm, SC | 10.4 | 7.0 |

Pt28 is two orders of magnitude above the rest. A scalp subcutaneous
melanoma biopsy can easily pick up keratinocytes, which express both
genes at high levels. Treating that sample as tumor TACSTD2/CLDN4 is not
justified. Rank tests are more robust than means, but the honest
description of this cohort is still: **these genes are not expressed in
the melanomas**.

## What this does and does not say about claim B5

Claim B5 is an 11-cohort ICI meta-analysis in which CLDN4-high has worse
response (user-reported OR = 0.42). This folder is the Hugo 2016 melanoma
piece only.

- CLDN4 does not differ by response (AUC 0.54, p = 0.74).
- The median-split OR is 1.78, opposite the claimed direction, with a CI
  that includes both 0.42 and large effects in the other direction.
- TACSTD2 ranks slightly lower in responders (AUC 0.35, p = 0.21). That
  trend is not significant, sits on background FPKM, and is not the B5
  claim.
- A melanoma RNA cohort is a weak place to test an epithelial tight-junction
  marker. Null here does not prove the meta-claim false in carcinoma. It
  does mean Hugo 2016 cannot be counted as supporting evidence.

## Important limitations

1. n = 26 is small. Every interval is wide. Absence of a significant p-value
   is not proof of a zero effect; it is proof that this dataset cannot
   resolve one.
2. No SD category; CR and PR are pooled. Response depth cannot be modeled
   further from the public labels.
3. GEO metadata do not include tumor purity, melanoma subtype (cutaneous vs
   occult), or a reliable contamination flag. Pt28 cannot be proven to be
   keratinocyte mix-in from these files alone; it can only be flagged as
   biologically implausible tumor expression.
4. FPKM is acceptable for a within-gene, between-sample rank test. It is
   not a raw-count differential-expression analysis.
5. Previous MAPKi is recorded but n is too small for an adjusted model.
   This analysis does not test IPRES or any other Hugo signature.
6. Association in treated patients cannot distinguish a predictive
   pembrolizumab interaction from a general prognostic effect.

## Reproduction

From this directory:

```bash
python3 -m pip install -r requirements.txt
python3 analyze.py
```

The script downloads the two public GEO inputs when absent, checks the
FPKM workbook MD5 (`548c3d627df111c107d83b3f7d33acf0`) and matrix
dimensions, and regenerates all tables and figures. `provenance.tsv`
records URLs, MD5, and SHA-256.

Outputs:

- `tacstd2_cldn4_vs_response.png` / `.svg`: primary and Pt28-dropped boxplots
- `fpkm_scatter.png` / `.svg`: raw FPKM scatter (log-like axes)
- `summary.csv`: tests, AUCs, CIs, detection, median-split ORs
- `patient_level_primary.csv`: 26-patient analysis set
- `sample_level.csv`: all 28 public biopsies
- `gene_correlation.csv`: TACSTD2–CLDN4 Spearman
- `provenance.tsv`: source URLs and checksums

## References

- Hugo W, Zaretsky JM, Sun L, et al. Genomic and Transcriptomic Features
  of Response to Anti-PD-1 Therapy in Metastatic Melanoma. *Cell*.
  2016;165:35–44.
  [doi:10.1016/j.cell.2016.02.065](https://doi.org/10.1016/j.cell.2016.02.065)
- NCBI GEO accession
  [GSE78220](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78220)
