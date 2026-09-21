# Concordant-4 Milo: wider sweep of T/NK neighbourhoods down in CLDN4-high-rich

Additive. The locked patient-level result stays the DerSimonian–Laird
Spearman on malignant CLDN4 %pos vs T/NK fraction (ρ = −0.531,
P = 1.65×10⁻⁵, I² = 0%, N = 65). The primary Milo model stays
k = 30, d = 30, Q3+Q4 vs Q1+Q2. This note does not replace either.

The question for this grid was fixed before the new fits. Among
full-rank models that keep all 65 locked units, with at least 10
units in each arm, maximise both of:

- the number of T/NK-majority neighbourhoods with SpatialFDR < 0.05 and log2FC < 0
- the median |log2FC| of those neighbourhoods

There is no single specification that wins both. The result is the
Pareto front of that pair. Three labeled rows are also reported.
The count maximum breaks ties by a larger median |log2FC|, then by
a more negative median log2FC of all tested T/NK neighbourhoods.
The |log2FC| maximum is taken among specifications with at least 20
down-hits, so a handful of extreme neighbourhoods cannot win it.
The joint score is the product of the count and that median.

## What was fit

Same graph rules as the primary: Euclidean kNN on the Harmony
embedding, undirected, Milo median refinement, prop = 0.1, seed = 1.
A neighbourhood is tested when cells from at least 5 of the 65 units
fall in it. The GLM is edgeR `glmQLFit(robust=TRUE)` + `glmQLFTest`
(edgeR 4.0.16), TMM, library size = column sums, k-distance
SpatialFDR, design `~ dataset + high`. Positive logFC is higher
abundance in the CLDN4-high arm. The sample is the locked unit.
The cell count is not n.

k ∈ {5, 8, 10, 12, 15, 20, 25, 30, 40, 50, 60}.
d ∈ {3, 5, 8, 10, 15, 20, 30}, the first d of the 30 Harmony
dimensions. Harmony has no further dimensions, so d = 30 is the
edge of the embedding.

Every contrast labels all 65 units. None drops the middle of the
cohort. Within-cohort cuts rank `mal_CLDN4_pct` inside each dataset
with `method="first"`.

| contrast | high | low | rule |
|---|---:|---:|---|
| q4_vs_rest | 16 | 49 | locked quartile Q4 vs Q1+Q2+Q3 |
| q34 | 31 | 34 | locked Q3+Q4 vs Q1+Q2 |
| top15 | 12 | 53 | within-cohort top 15% vs the rest |
| top20 | 15 | 50 | within-cohort top 20% vs the rest |
| top25 | 19 | 46 | within-cohort top 25% vs the rest |
| top33 | 23 | 42 | within-cohort top third vs the rest |
| top40 | 28 | 37 | within-cohort top 40% vs the rest |
| half | 34 | 31 | within-cohort upper half vs lower half |

Per dataset the high arm is at least 2 and the low arm is at least 4.
Balance: `results/tables/max_effect_contrast_balance.tsv`.

77 graphs × 8 contrasts = 616 fits. All 616 were full rank (5/5)
with n_samples = 65. None was dropped by the arm-size rule.

## Checks against the earlier grid

The same code, on cells that the earlier grid already published,
reproduces those rows.

| k | d | contrast | T/NK down | T/NK up | median \|log2FC\| of the down hits |
|---:|---:|---|---:|---:|---:|
| 15 | 10 | q4_vs_rest | 65 | 0 | 1.991 |
| 30 | 30 | q4_vs_rest | 41 | 0 | 2.729 |
| 30 | 30 | q34 | 23 | 1 | 2.331 |
| 15 | 10 | half | 45 | 0 | 1.604 |

The last row is the earlier upper-half specification. Minimum
SpatialFDR on the k = 30, d = 30, Q3+Q4 refit is 1.70×10⁻⁸, the
same minimum as the primary table.

## Pareto front

18 of 616 specifications are undominated. A row is on the front
when no other eligible specification has at least as many T/NK-down
hits and at least as large a median |log2FC|, with one of the two
strictly larger. Full precision is
`results/tables/max_effect_pareto.tsv`.

| k | d | contrast | high / low | down | up | median \|log2FC\| of down hits | median log2FC, all tested T/NK | T/NK tested | down / tested |
|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 12 | 15 | top40 | 28 / 37 | 111 | 4 | 1.596 | −0.687 | 708 | 15.7% |
| 15 | 5 | top33 | 23 / 42 | 103 | 0 | 1.681 | −0.685 | 781 | 13.2% |
| 10 | 5 | q4_vs_rest | 16 / 49 | 99 | 1 | 1.708 | −0.670 | 780 | 12.7% |
| 12 | 5 | q4_vs_rest | 16 / 49 | 93 | 1 | 1.802 | −0.638 | 788 | 11.8% |
| 20 | 5 | q4_vs_rest | 16 / 49 | 92 | 0 | 1.917 | −0.720 | 773 | 11.9% |
| 10 | 30 | top20 | 15 / 50 | 83 | 2 | 1.961 | −0.759 | 613 | 13.5% |
| 25 | 3 | q4_vs_rest | 16 / 49 | 77 | 0 | 1.975 | −0.809 | 729 | 10.6% |
| 15 | 5 | q4_vs_rest | 16 / 49 | 74 | 0 | 2.013 | −0.701 | 781 | 9.5% |
| 15 | 5 | top20 | 15 / 50 | 73 | 0 | 2.015 | −0.674 | 781 | 9.3% |
| 20 | 30 | q4_vs_rest | 16 / 49 | 65 | 0 | 2.421 | −0.813 | 552 | 11.8% |
| 25 | 30 | q4_vs_rest | 16 / 49 | 43 | 0 | 2.539 | −0.777 | 534 | 8.1% |
| 30 | 30 | q4_vs_rest | 16 / 49 | 41 | 0 | 2.729 | −0.812 | 536 | 7.6% |
| 40 | 30 | q4_vs_rest | 16 / 49 | 32 | 0 | 2.908 | −0.760 | 476 | 6.7% |
| 50 | 30 | q4_vs_rest | 16 / 49 | 29 | 0 | 2.934 | −0.838 | 445 | 6.5% |
| 40 | 30 | top20 | 15 / 50 | 24 | 0 | 3.102 | −0.736 | 476 | 5.0% |
| 40 | 30 | top15 | 12 / 53 | 8 | 0 | 3.466 | −0.511 | 476 | 1.7% |
| 50 | 30 | top15 | 12 / 53 | 5 | 0 | 3.544 | −0.608 | 445 | 1.1% |
| 60 | 30 | top15 | 12 / 53 | 2 | 0 | 4.033 | −0.488 | 434 | 0.5% |

The earlier count maximum (k = 15, d = 10, Q4 vs the rest: 65 down,
median |log2FC| 1.991) is inside this grid and is not on the front.
k = 20, d = 30, Q4 vs the rest has the same 65 down-hits with median
|log2FC| 2.421 and 0 up.

The right-hand end of the front is two neighbourhoods at median
|log2FC| 4.033 (k = 60, d = 30, top 15%). That row fails the
pre-specified floor of 20 hits, so it is not the |log2FC| label.
It is on the front because nothing else with a larger median has
at least as many hits.

## Labeled rows

**Count maximum, which is also the joint-score maximum.**
k = 12, d = 15, within-cohort top 40% vs the rest (28 vs 37).
111 T/NK neighbourhoods down and 4 up, out of 708 tested T/NK
neighbourhoods (15.7%). Median |log2FC| of the 111 is 1.596
(range −2.592 to −1.199). Median log2FC of all tested T/NK
neighbourhoods is −0.687. Median size of the down-hits is 30 cells.
The smallest SpatialFDR among the 111 is 1.31×10⁻⁴. The smallest
SpatialFDR in the whole fit is 3.11×10⁻⁶. Joint score
111 × 1.596 = 177.1.
The next product on the front is k = 20, d = 5, Q4 vs the rest:
92 × 1.917 = 176.4. The count winner is not separated from that
neighbour by a large gap.

Those 111 neighbourhoods are shared. Patients per neighbourhood:
median 15, minimum 5, maximum 26; 96 of 111 contain cells from at
least 10 patients. Median fraction of cells from the 28 high units
is 0.13, and from the 37 low units is 0.87. Cells from GSE123902,
GSE131907, and GSE205335 appear in 108, 109, and 109 of the 111;
GSE189357 appears in 99. Malignant neighbourhoods in the same fit
are 85 up and 35 down (median log2FC +0.41).

**Median |log2FC| maximum among specifications with at least 20 down-hits.**
k = 40, d = 30, within-cohort top 20% vs the rest (15 vs 50).
24 T/NK neighbourhoods down and 0 up, out of 476 tested (5.0%).
Median |log2FC| 3.102 (range −3.968 to −2.375). Median log2FC of
all tested T/NK neighbourhoods is −0.736. Median size of the hits
is 98 cells. Every hit contains cells from at least 10 patients
(median 16, minimum 13, maximum 29). The smallest SpatialFDR
among the 24 is 2.05×10⁻³. The smallest SpatialFDR in the whole
fit is 7.01×10⁻¹¹. Malignant neighbourhoods in the same fit are
92 up and 35 down.

On the locked Q3+Q4 cut, the widest count in this grid is k = 15,
d = 5: 83 down and 2 up, median |log2FC| 1.506. The primary
specification remains 23 down and 1 up.

## SpatialFDR cuts on the same fits

Loosening the threshold adds hits with smaller |log2FC|. It was
recorded and was not the selection rule. On the count-maximum fit:

| SpatialFDR cut | T/NK down | median \|log2FC\| of those hits |
|---:|---:|---:|
| 0.01 | 38 | 2.015 |
| 0.05 | 111 | 1.596 |
| 0.10 | 166 | 1.476 |
| 0.20 | 267 | 1.251 |

## What to quote

Quote n = 65 for every row in this grid. The count maximum is
111 T/NK neighbourhoods down at SpatialFDR < 0.05 (k = 12, d = 15,
top 40% vs the rest, 28 vs 37), with median |log2FC| 1.596 and 4
neighbourhoods up. The largest median |log2FC| among specifications
with at least 20 down-hits is 3.102, on 24 neighbourhoods
(k = 40, d = 30, top 20% vs the rest, 15 vs 50), with 0 up.
The primary model remains k = 30, d = 30, Q3+Q4 vs Q1+Q2.

Figures: `figures/max_effect_pareto.png`,
`figures/max_effect_count_heatmaps.png`.
Grid: `results/tables/max_effect_grid_eligible.tsv`.
Down-hit tables: `results/tables/tnk_down_max_effect_count.tsv`,
`results/tables/tnk_down_max_effect_abslogfc.tsv`.
