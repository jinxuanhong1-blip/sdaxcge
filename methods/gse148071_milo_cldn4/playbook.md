# kNN neighbourhood DA vs malignant CLDN4 (GSE148071)

**Scope.** Additive test of whether transcriptional neighbourhoods change with a
sample-level **malignant CLDN4** score on public GSE148071 counts (Wu et al.,
*Nat Commun* 2021; 42 advanced NSCLC diagnostic biopsies). Prior TACSTD2 / CLDN4
Milo on GSE207422 is taken as given.

**Not miloR.** This environment does not ship R / Bioconductor / miloR / edgeR.
`scripts/knn_nhood.py` reimplements neighbourhood construction and **SpatialFDR**
from Dann et al., *Nat Biotechnol* 2022 (miloR `graphSpatialFDR`, `k-distance`
weights; Lun et al. cydar). The DA *model* is a sample-level Spearman (primary)
or Welch t-test (median-split secondary) on neighbourhood proportions. Do not
cite these p-values as miloR output.

## Graph

1. HVG by raw-count dispersion (mean in [0.01, 50], top 2000) across all cells.
2. log1p(CP10k), gene-scale, PCA (`d=30`).
3. Subtract the per-sample PCA mean (one-step batch centering; **not** Harmony).
4. Euclidean kNN (`k=30`).
5. Milo refined index sampling (`prop=0.1`).
6. Neighbourhood = index ∪ its k nearest neighbours (size `k+1`).

## DA (fallback)

```
prop[s, i] = n_cells(s in i) / n_cells(s)
```

- **Malignant CLDN4 (primary):** Spearman of `prop` vs the sample’s mean
  log1p-CP10k CLDN4 in malignant-like cells (sample dropped if <10 such cells).
- **CLDN4 high vs low (secondary, same graph):** Welch t-test after a median
  split of the samples that have a malignant CLDN4 score.

A neighbourhood is testable only if enough samples contribute cells. Overlapping
neighbourhoods are **not** independent observations. On this atlas, cancer cells
cluster by patient (Wu 2021), so most neighbourhoods are patient-private and
fail the ≥5-sample rule. |ρ|=1 at the 5–6 sample floor is a scipy perfect-rank
artifact, not a cohort DA claim. Report the n_present ≥ 8 sensitivity.

## SpatialFDR (honest)

`w = 1 / (k-th NN distance of the index cell)`

```
adjp[order] = rev(cummin(rev(sum(w) * p / cumsum(w))))
```

Report **all** of: `n_nhoods`, `n_testable`, `n` samples, raw p<0.05, BH-FDR,
SpatialFDR. Do not quote ~90k cells as *n*.

## Composition

Dissociated kNN neighbourhoods are **transcriptional**, not spatial.

1. **Interface neighbourhoods** — ≥3 malignant and ≥3 T/NK cells. Spearman of
   malignant CLDN4 vs T/NK fraction.
2. **Disjoint subset** — greedy non-overlapping neighbourhoods (honest *n*).
3. **Sample-paired** — for each sample, T/NK fraction among that sample’s cells
   that fall in CLDN4-high vs CLDN4-low neighbourhoods (median split on
   neighbourhoods with ≥5 malignant cells). Wilcoxon signed-rank across samples.

(3) is the only composition test whose independent unit is the patient.

## Public dataset

| Dataset | Matrix | Labels | Response |
|---|---|---|---|
| GSE148071 (Wu 2021) | GEO `GSE148071_RAW.tar` per-sample count TSVs | marker-reconstructed (author barcodes / CNA IDs not public) | not an ICI-response cohort; advanced diagnostic biopsies |

GEO SOFT has age/sex only. Histology (LUAD/LUSC) is in the paper supplement, not
on the series matrix, and is not used as a DA covariate here.

Raw SRA FASTQ are not used.
