# Methods

Concordant-4 only. The locked patient-level result (malignant CLDN4 % positive
versus T/NK fraction, ρ = −0.531, n = 65) is not re-fit here. This folder asks
a different question: which genes distinguish CLDN4-high from CLDN4-low
malignant cells, and which of those genes do so in every cohort.

## Cohorts and malignant cells

| cohort | unit | malignant definition |
|---|---|---|
| GSE123902 | donor, one tumour sample (primary if both exist; normal dropped) | (EPCAM or KRT8 or KRT18 or KRT19) > 0 and PTPRC = 0 |
| GSE189357 | patient | same marker gate |
| GSE131907 | locked tumour-bearing sample | author `Cell_subtype == Malignant cells` |
| GSE205335 | patient, non-normal tissues pooled | author `lineage.sub == Malignant cells` |

Counts are GEO UMIs. GSE205335 is the double-gzipped RDS. Symbols are
upper-cased; duplicate symbols are summed.

## CLDN4-high versus CLDN4-low

Inside each unit, malignant cells are ranked by log1p(CP10k CLDN4). Ties keep
the earlier cell at the lower rank, matching R `rank(..., ties.method="first")`.
Q1 is rank ≤ floor(0.25 n). Q4 is rank > ceiling(0.75 n).

A unit is used only when all of these hold:

- n malignant ≥ 40
- at least 10 malignant cells with CLDN4 count > 0
- each arm has at least 10 cells
- Q4 mean CLDN4 > Q1 mean CLDN4

That last gate drops units with no real CLDN4 range (a quartile split of all
zeros is an arbitrary label). P4001 (27 malignant cells) is out on the size
gate. The unit of the reported Δ is the patient, donor, or sample. Pooled cell
counts are only the cells that entered the ranking.

## COMET / XL-mHG

For each gene, arm cells (Q4 and Q1 only) are ranked by log1p(CP10k), high
expression first. Ties are broken against the positive class. The XL-mHG
statistic is the minimum hypergeometric tail over cutoffs, from xlmhg 2.5.4
(`X = 1`, `L = N`). Two tests are run:

- UP: positive class = CLDN4-high cells
- DOWN: positive class = CLDN4-low cells

The statistic is a ranking score. It is not used as a p-value. Exact XL-mHG
p-values are not required for the 4/4 call, because with this many cells they
are small whenever the effect is real. The patient-level Wilcoxon p on the
per-unit deltas is reported on the gene-list tables as a descriptive check.

COMET AND-pairs: for the top 20 UP genes, each gene is binarized at its XL-mHG
cutoff expression (threshold must be > 0). A pair passes in a cohort when the
AND set is hypergeometrically enriched for CLDN4-high cells (p < 0.01, fold >
1.5, at least 10 high cells in the draw). A 4/4 pair passes in every cohort.

## 4/4 rules

These floors were set before the lists were built.

**Strict patient-effect, UP.** In every cohort: patient-mean Δ ≥ 0.10, at least
60% of units have Δ > 0, pooled AUC ≥ 0.55, and the gene is detected in ≥ 5%
of CLDN4-high cells. DOWN mirrors this (Δ ≤ −0.10, ≥ 60% of units negative,
AUC ≤ 0.45, detected in ≥ 5% of CLDN4-low cells).

**COMET top decile.** Inside each cohort, among genes with the matching 5%
detection floor, rank by the XL-mHG statistic and then by cell Δ. The top 10%
is the decile. A gene is COMET 4/4 when it is in that decile in every cohort
and the patient-mean Δ has the matching sign in every cohort.

**Same-sign.** Patient-mean Δ has the same sign in all four cohorts, with no
effect floor. Ranked by the weakest cohort. This is the backup Enrichr input
when the strict list has fewer than 15 genes.

CLDN4 is a positive control (it has to be up in every used unit) and is not a
marker on these lists. A companion strict list drops MT- / RPL / RPS / MRPL /
MRPS symbols.

## Modules and Enrichr

Spearman correlations are computed on a stratified sample of Q4 and Q1 cells
(at most 80 cells per unit per arm, at most 2000 cells, seed 1). Correlations
are averaged on the Fisher z scale across the four cohorts. Average-linkage
clustering uses distance `1 − r` and is cut at 0.5. Clusters with fewer than
4 genes are not modules.

The module input is the strict list when it has at least 8 genes, otherwise
the COMET decile, otherwise the top 60 same-sign genes.

Enrichment uses gseapy `enrichr` (Enrichr API) on GO Biological Process 2023,
KEGG 2021 Human, MSigDB Hallmark 2020, and Reactome 2022. If the API is
unreachable, the same libraries are downloaded and tested with a one-sided
Fisher exact test. That fallback uses the genes observed in all four cohorts
as the universe and is labeled `fisher_tested_universe`. Enrichr's own
background is the API default.

Input to Enrichr: the strict list when it has at least 15 genes, otherwise the
top 100 same-sign genes. The COMET decile is enriched as well when it has at
least 15 genes and was not already the strict input. Module gene sets of at
least 8 genes are enriched separately.
