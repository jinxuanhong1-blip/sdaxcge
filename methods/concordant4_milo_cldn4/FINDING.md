# Concordant-4 Milo: CLDN4-high-rich vs low-rich neighbourhoods

Additive. CLDN4-only. The locked patient-level result stays the
DerSimonian–Laird Spearman on malignant CLDN4 %pos vs T/NK fraction
(ρ = −0.531, P = 1.65×10⁻⁵, I² = 0%, N = 65). This note is the
neighbourhood test on the same 65 units.

The graph is a kNN on the Harmony embedding (`group.by.vars = dataset`
in the concordant-4 Seurat object). A neighbourhood is a transcriptional
state. Counts are cells inside that state, in a dissociated sample.
The object is capped at ≤350 QC cells per unit (22,653 cells in the
graph). Do not quote 22,653 as n.

## Design that was fit

High-rich vs low-rich uses the locked within-cohort quartiles:
**Q3+Q4 (31)** vs **Q1+Q2 (34)**.

| dataset | low (Q1+Q2) | high (Q3+Q4) |
|---|---:|---:|
| GSE123902 | 7 | 6 |
| GSE131907 | 11 | 10 |
| GSE189357 | 5 | 4 |
| GSE205335 | 11 | 11 |

Primary model, sample = locked unit:

```
~ dataset + cldn4_high
```

`cldn4_high` is the patient-level covariate (1 = high-rich). `dataset`
is the batch covariate. Positive logFC is higher abundance in
CLDN4-high-rich units. edgeR logFC is log2.

DA is `estimateDisp` + `glmQLFit(robust=TRUE)` + `glmQLFTest`, the same
quasi-likelihood F test `miloR::testNhoods` calls. Normalisation for the
primary model is TMM with library size = column sums of the neighbourhood
count matrix. SpatialFDR uses k-distance weights (miloR default).

Graph: Euclidean k = 30 on 30 Harmony dimensions, undirected kNN,
Milo median refinement, `prop = 0.1`, seed = 1. 1,502 neighbourhoods
(median size 69). A neighbourhood enters the primary test when cells
from **≥5 patients** fall in it: **1,287** tested, **215** held out
(211 of those 215 are malignant). Median patients per neighbourhood
among all 1,502 is 19 (range 1–45).

## Within-patient model

GEO series-matrix patient ids for the 21 GSE131907 tumour-bearing
samples are 21 distinct patients (P1006–P1058 and P3002–P3019). The
other three cohorts already use one donor or patient per row.
**0 / 65 patients have two samples.**

`~ patient_id + cldn4_high` has rank **65 / 66** and is not estimated.
There is no second sampled condition inside a patient, so a
within-patient malignant-neighbourhood contrast is not a Milo design
on this object. Per-cell CLDN4 is not stored on the Harmony embedding;
splitting one suspension into CLDN4-high and CLDN4-low pseudo-samples
would reuse the same cells. That split was not run.

## Primary result

**253 / 1,287** neighbourhoods at SpatialFDR < 0.10, and **202 / 1,287**
at SpatialFDR < 0.05. Minimum SpatialFDR = 1.7×10⁻⁸.
BH-FDR < 0.10 is 279 / 1,287.

Neighbourhoods are class-pure (median majority fraction 1).
**18 / 1,502** contain ≥3 malignant cells and ≥3 T/NK cells.

Hits at SpatialFDR < 0.10, by class of the tested neighbourhoods:

| class | tested | down in high-rich | up in high-rich | median log2FC (all tested) |
|---|---:|---:|---:|---:|
| T/NK | 536 | 36 | 1 | −0.55 |
| malignant | 266 | 26 | 118 | +1.35 |
| myeloid | 315 | 31 | 12 | −0.28 |
| B | 96 | 5 | 0 | −0.60 |
| other | 74 | 24 | 0 | −0.71 |

81% of tested T/NK neighbourhoods have log2FC < 0. 71% of tested
malignant neighbourhoods have log2FC > 0.

Hits that are also seen in ≥10 patients: **158**. Of those,
T/NK is 33 down and 1 up (hit median 16 patients), malignant is
51 up and 5 down, myeloid is 30 down and 11 up, B is 5 down, other
is 22 down. The malignant neighbourhoods that are lower in high-rich
units are mostly thinly shared (26 down-hits, median 6 patients,
5 of them in ≥10 patients).

The same direction is in the five class totals on this capped object
(log library-size offset, BH-FDR across 5 classes, not SpatialFDR):

| class | log2FC | P | BH-FDR |
|---|---:|---:|---:|
| malignant | +1.49 | 1.1×10⁻⁵ | 5.6×10⁻⁵ |
| T/NK | −0.73 | 0.017 | 0.028 |
| other | −1.35 | 0.011 | 0.027 |
| B | −1.01 | 0.052 | 0.065 |
| myeloid | −0.49 | 0.081 | 0.081 |

Group medians of the locked T/NK fraction are 0.42 (low-rich) and 0.17
(high-rich). In the capped object the embedded malignant fraction
medians are 0.09 and 0.51. The neighbourhood log fold-changes follow
that composition: high-rich units place more of the 350-cell cap in
malignant states and less in T/NK states.

## Sensitivities

| model | units | SpatialFDR < 0.10 | SpatialFDR < 0.05 | min SpatialFDR |
|---|---:|---:|---:|---:|
| Primary TMM, high vs low | 65 | 253 | 202 | 1.7×10⁻⁸ |
| logMS offset = cells/unit | 65 | 322 | 223 | 1.6×10⁻⁸ |
| Continuous CLDN4, per 10 percentage points, TMM | 65 | 257 | 202 | 3.8×10⁻¹³ |
| Primary + GSE131907 origin (`mLN`, `tL/B` vs `mBrain`) | 65 | 250 | 169 | 9.5×10⁻⁷ |
| Q4 vs Q1 only (16 vs 19) | 35 | 233 | 160 | 3.9×10⁻⁴ |
| GSE131907 only, `~ tissue + cldn4_high` | 21 | 65 | 41 | 0.018 |
| GSE123902 only | 13 | 75 | 0 | 0.056 |
| GSE131907 only, no origin term | 21 | 133 | 109 | 0.0073 |
| GSE205335 only | 22 | 126 | 91 | 0.014 |
| GSE189357 only | 9 | **0** | **0** | 0.83 |

218 of the 253 primary hits remain hits after the GSE131907 origin
terms. Continuous CLDN4 (patient-level %pos / 10) gives 58 T/NK hits,
all down, and 127 malignant hits up vs 19 down.

GSE189357 alone (9 patients, 4 vs 5) does not produce a SpatialFDR < 0.10
neighbourhood. That is the small-n behaviour of the earlier single-cohort
Milo runs. The pooled hit list is carried by the larger cohorts plus the
shared Harmony graph, under a common log fold-change and a dataset intercept.

A malignant-cell offset (library size = malignant cells in the cap,
neighbourhoods with ≥10 malignant cells) was fit as a check on
redistribution inside the malignant compartment. It was **not** limited
to the ≥5-patient rule. Its SpatialFDR < 0.10 calls sit in thinly shared
neighbourhoods (up-hits have median 2 patients). That table is stored
and is not a second result.

## Why earlier Milo tables were empty

Earlier concordant-cohort Milo runs (GSE131907 tLung, GSE205335, and the
GSE131907+GSE205335 merge) tested neighbourhood **proportions with
Spearman**. After |ρ| = 1 calls at very small n were set aside, SpatialFDR
< 0.10 was 0. Those graphs were single-cohort PCA graphs, and most
neighbourhoods were private to one sample.

This retry keeps the locked 65-unit cohort on one Harmony graph, so the
tested neighbourhoods are shared (1,287 with ≥5 patients). The DA model
is the edgeR quasi-likelihood GLM with the patient-level CLDN4 group as
the covariate and dataset as a batch term. Under that design the
composition shift is large enough to pass SpatialFDR. The GSE189357-only
fit is still empty, which is the honest small-n limit.

## What this is for

Use it as neighbourhood-resolved abundance between CLDN4-high-rich and
low-rich patients: malignant states up, T/NK states down, on dissociated
cells. The patient-level T/NK association remains the locked Spearman.
Quote **n = 65** (31 vs 34) for the primary contrast, and **n = 35**
(16 vs 19) for Q4 vs Q1.

Figures: `figures/volcano_binary_tmm.png`, `figures/logfc_by_class_binary_tmm.png`.
Primary table: `results/tables/da_binary_tmm.tsv`.
Design ranks: `results/tables/design_rank.tsv`.
