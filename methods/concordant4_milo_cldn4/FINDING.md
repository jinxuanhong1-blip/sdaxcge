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

## Sensitivity grid at SpatialFDR < 0.05

The primary model stays k = 30, d = 30, Q3+Q4 vs Q1+Q2. The grid
below was fixed before the new fits: k ∈ {15, 30, 50}, d ∈ {10, 20, 30}
(the first d Harmony dimensions), and five high-rich definitions.
A neighbourhood is tested when ≥5 patients **in that contrast**
contribute a cell. The search, among full-graph specs with at least
10 units in each arm, maximises the number of T/NK neighbourhoods
with SpatialFDR < 0.05 and log2FC < 0. Ties go to the more negative
median log2FC of those hits, then of all tested T/NK neighbourhoods.

| contrast | high | low | rule |
|---|---:|---:|---|
| Q3+Q4 vs Q1+Q2 | 31 | 34 | primary |
| Q4 vs Q1 | 16 | 19 | drop Q2 and Q3 |
| Q4 vs the rest | 16 | 49 | Q1–Q3 are the low arm |
| upper half vs lower half | 34 | 31 | within-cohort rank; the middle of an odd cohort is low |
| top vs bottom tertile | 21 | 23 | middle tertile dropped |

The count maximum is **k = 15, d = 10, Q4 vs the rest**: **65** T/NK
neighbourhoods down and **0** up, out of 725 tested T/NK neighbourhoods
(9.0%). Median log2FC of tested T/NK neighbourhoods is −0.73; of the
65 hits, −1.99. Malignant neighbourhoods in that fit are 99 up and
30 down. This row wins the pre-specified count. It uses the coarser
embedding and a low arm that includes Q2 and Q3. It does not replace
the primary.

The primary specification, refit inside the grid, is **23** T/NK
neighbourhoods down and **1** up (536 tested, 4.3%). Median log2FC
of tested T/NK neighbourhoods is −0.55; of the 23 down-hits, −2.33.
Malignant neighbourhoods are 104 up and 25 down. Minimum SpatialFDR
is 1.70×10⁻⁸, the same minimum as the primary table.

The most negative median T/NK log2FC among specs with at least 20
T/NK-down hits is **k = 30, d = 30, Q4 vs Q1** (16 vs 19): **42** down
and **0** up, median log2FC −1.06 for tested T/NK neighbourhoods and
−2.88 for the hits. That filter is applied inside the 35-unit contrast
(1,139 tested neighbourhoods). The earlier Q4-vs-Q1 table filtered on
all 65 units first, so the two rows are different fits.

The full 45-spec table is `results/tables/grid_full_ranked.tsv`.
The count heatmap is `figures/grid_tnk_down_sfdr05.png`.

## Malignant-only graph

The same k and d values were refit on the 7,424 embedded malignant
cells alone, for Q3+Q4 and for Q4 vs Q1. There are no T/NK cells on
this graph, so it cannot increase the T/NK-down count.

At the primary k and d, Q3+Q4 vs Q1+Q2: **281** tested neighbourhoods,
**71** up and **17** down at SpatialFDR < 0.05, median log2FC +1.03,
minimum SpatialFDR 1.64×10⁻⁷. The largest up count on this graph is
k = 50, d = 10, same contrast: 161 up and 2 down (499 tested, median
log2FC +1.08). These are malignant transcriptional states that are
more abundant in CLDN4-high-rich units.

## IFN, NHEJ, and STING inside differential neighbourhoods

Gene lists are in `data/gene_sets.json`.

| set | n | source |
|---|---:|---|
| IFN | 224 | MSigDB Hallmark IFNα ∪ IFNγ |
| NHEJ | 78 | GO:0006303, non-homologous end joining |
| STING | 16 | Reactome R-HSA-1834941 |
| cGAS–STING (secondary) | 29 | GO:0140896 |

Raw counts are the public matrices that built the Harmony object.
Every embedded cell matched a matrix column (4,453 + 7,350 + 3,150
+ 7,700). The score is the mean of log1p(count / nCount_RNA × 10⁴)
over the fixed gene list. A gene missing from a matrix contributes 0.
`nCount_RNA` matched the matrix column sum on the cells that were
checked (GSE123902 20/20, GSE131907 3/3, GSE189357 45/45, GSE205335
20/20).

Where the current symbol is absent and one HGNC previous symbol is
present, the count is read from that alias: STING1 = TMEM173,
CGAS = MB21D1, MRE11 = MRE11A, MARCHF1 = MARCH1, MARCHF5 = MARCH5,
WARS1 = WARS, TENT5A = FAM46A, NSD2 = WHSC1, CYREN = C7orf49,
MRNIP = C5orf45, PAXX = C9orf142, SHLD1 = C20orf196, SHLD2 = FAM35A.
MIR4691 is absent from all four matrices. SHLD3 and TREX1 are absent
from every GSE123902 donor file. SHLD3 and ATP23 are absent from
GSE131907. GSE123902 drops genes that are zero in a donor, so a few
other symbols are zero-scored in some donors only. Coverage is
`results/tables/gene_coverage.tsv`.

The rebuilt graphs match the cached neighbourhood counts, and the
recomputed SpatialFDR matches the saved values (maximum absolute gap
1×10⁻¹⁶). A cell in both an up-hit and a down-hit is kept in the up
compartment (70 such cells on the primary malignant hits). Within
each locked unit, Δ is the mean score in the focal compartment minus
the mean score in other neighbourhoods of the same class. The
Wilcoxon signed-rank is across units that contain both compartments.
The unit of that test is the locked sample.

Primary graph, malignant cells in malignant-up neighbourhoods
(SpatialFDR < 0.05 and log2FC > 0; 2,770 cells) versus other
malignant-majority neighbourhoods (3,952 cells):

| set | units | median Δ | units negative / positive | P |
|---|---:|---:|---:|---:|
| IFN | 61 | −0.0068 | 32 / 29 | 0.47 |
| NHEJ | 61 | +0.0015 | 30 / 31 | 0.59 |
| STING | 61 | +0.0025 | 30 / 31 | 0.77 |

Median IFN scores in those 61 units are 0.24 in the up compartment
and 0.26 in the other malignant compartment. The within-unit IFN
difference is centered at zero. The secondary cGAS–STING set on the
same contrast has median Δ +0.0033 (P = 0.81).

Primary graph, T/NK cells in T/NK-down neighbourhoods (997 cells)
versus other T/NK-majority neighbourhoods (7,012 cells). Both sides
are T/NK cells from the same locked unit:

| set | units | median Δ | units negative / positive | P |
|---|---:|---:|---:|---:|
| IFN | 53 | −0.036 | 38 / 15 | 0.0016 |
| NHEJ | 53 | −0.024 | 44 / 9 | 1.2×10⁻⁶ |
| STING | 53 | −0.019 | 38 / 15 | 0.0034 |

The single T/NK neighbourhood that is up in high-rich units holds
62 embedded cells. In the 22 units that also have other T/NK
neighbourhoods, its IFN score is higher (median Δ +0.069, 20/22
positive, P = 9.9×10⁻⁵).

On the count-maximum graph, malignant-up versus other IFN has median
Δ +0.0055 (58 units, P = 0.22). T/NK-down versus other IFN has median
Δ −0.0099 (58 units, P = 0.18). NHEJ in those T/NK-down neighbourhoods
is lower (median Δ −0.011, 43/58 negative, P = 0.0016).

On the malignant-only graph at k = 30, d = 30, Q3+Q4, malignant-up
versus other IFN has median Δ **+0.023** (54 units, 32 positive,
P = 0.021). NHEJ on that contrast has P = 0.45 and STING has P = 0.28.

Patient-level deltas: `results/tables/nhood_gene_patient_deltas.tsv`.
Summary: `results/tables/nhood_gene_paired.tsv`.
Figure: `figures/nhood_ifn_delta_primary.png`.

## What this is for

Use it as neighbourhood-resolved abundance between CLDN4-high-rich and
low-rich patients: malignant states up, T/NK states down, on dissociated
cells. The patient-level T/NK association remains the locked Spearman.
Inside a sample, malignant-up neighbourhoods and other malignant
neighbourhoods have the same IFN score. Quote **n = 65** (31 vs 34)
for the primary contrast, **n = 35** (16 vs 19) for Q4 vs Q1, and
**n = 65** (16 vs 49) for Q4 versus the rest.

A wider search, k ∈ {5, 8, 10, 12, 15, 20, 25, 30, 40, 50, 60},
d ∈ {3, 5, 8, 10, 15, 20, 30}, and eight high-rich cuts that keep
all 65 units, is in `MAX_EFFECT.md`. The count of 65 above is the
maximum inside the 45-specification grid on this page. The primary
model is unchanged.

Figures: `figures/volcano_binary_tmm.png`, `figures/logfc_by_class_binary_tmm.png`,
`figures/grid_tnk_down_sfdr05.png`, `figures/nhood_ifn_delta_primary.png`.
Primary table: `results/tables/da_binary_tmm.tsv`.
Grid: `results/tables/grid_full_ranked.tsv`, `results/tables/grid_malignant_only.tsv`.
Gene scores: `results/tables/nhood_gene_paired.tsv`.
Design ranks: `results/tables/design_rank.tsv`.
