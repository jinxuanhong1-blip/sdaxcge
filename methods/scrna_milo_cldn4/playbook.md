# kNN neighbourhood DA vs malignant CLDN4 (Milo fallback)

**Scope.** Additive test of whether transcriptional neighbourhoods change with a
sample-level **malignant CLDN4** score on public GSE207422 UMI, and whether
CLDN4-high neighbourhoods are T/NK-poor. Prior TACSTD2 Milo is taken as given.

**Not miloR.** This environment does not ship R / Bioconductor / miloR / edgeR.
`scripts/knn_nhood.py` reimplements neighbourhood construction and **SpatialFDR**
from Dann et al., *Nat Biotechnol* 2022 (miloR `graphSpatialFDR`, `k-distance`
weights; Lun et al. cydar). The DA *model* is a sample-level Welch t-test or
Spearman on neighbourhood proportions. Do not cite these p-values as miloR output.

## Graph (same as prior TACSTD2 Milo)

1. HVG by raw-UMI dispersion (mean in [0.01, 50], top 2000).
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
  log1p-CP10k CLDN4 in malignant/epithelial cells (sample dropped if <10 such cells).
- **MPR vs NMPR (secondary, same graph):** Welch t-test. `logFC = log2(mean_NMPR) − log2(mean_MPR)`.

A neighbourhood is testable only if enough samples contribute cells. Overlapping
neighbourhoods are **not** independent observations.

## SpatialFDR (honest)

`w = 1 / (k-th NN distance of the index cell)`

```
adjp[order] = rev(cummin(rev(sum(w) * p / cumsum(w))))
```

Report **all** of: `n_nhoods`, `n_testable`, `n` samples, raw p<0.05, BH-FDR,
SpatialFDR. Do not quote 90k cells as *n*.

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
| GSE207422 (Hu 2023) | GEO UMI TSV, ~92k cells | marker-reconstructed (author barcodes not public) | post-treatment MPR (incl. pCR) vs NMPR |

Raw GSA-Human / SRA FASTQ are not used. GSE241934 is out of scope here.
