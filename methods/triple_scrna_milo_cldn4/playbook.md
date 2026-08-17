# Triple-merge kNN neighbourhood DA vs malignant CLDN4

**Scope.** One additive question: do transcriptional neighbourhoods change
with **sample-level malignant CLDN4** when GSE131907, GSE148071 and
GSE205335 are analysed together? Patient/sample is the unit.

**Not miloR.** This environment does not ship R / Bioconductor / miloR /
edgeR. `scripts/knn_nhood.py` reimplements neighbourhood construction and
**SpatialFDR** from Dann et al., *Nat Biotechnol* 2022 (miloR
`graphSpatialFDR`, k-distance weights). The DA model is a sample-level
Spearman (primary) or Welch t-test (median split). Do not cite these
p-values as miloR / edgeR QLF output.

## Graph (per dataset)

Joint Harmony on the concatenated public matrices is not run: ~370k cells
and a joint log-normalized HVG matrix do not fit in 15 GB. The user-allowed
fallback is used:

1. HVG by raw-UMI dispersion (mean in [0.01, 50], top 2000), **inside each dataset**.
2. log1p(CP10k), gene-scale, PCA (`d=30`).
3. Subtract the per-sample PCA mean (one-step batch centering; **not** Harmony).
4. Euclidean kNN (`k=30`).
5. Milo refined index sampling (`prop=0.1`).
6. Neighbourhood = index ∪ its k nearest neighbours (size `k+1`).

GSE131907 uses **tLung only** (primary lung tumor; epithelium + immune).
Normals and other sites are not mixed into that graph. GSE148071 uses every
diagnostic biopsy. GSE205335 drops author-labelled normal tissues and pools
biopsies from the same patient.

## DA

```
prop[s, i] = n_cells(s in i) / n_cells(s)
```

Spearman of `prop` vs mean log1p-CP10k CLDN4 in malignant cells (sample
dropped if <10 such cells). SpatialFDR on testable neighbourhoods
(`min_samples=5`). Report n_nhoods, n_testable, n samples, raw p<0.05,
BH-FDR and SpatialFDR. Do not quote cell count as *n*.

## What is not done

- Dual-high TACSTD2+CLDN4 gates.
- GSE207422 (already analysed elsewhere).
- miloR `testNhoods`.
- Treating overlapping neighbourhoods as independent patients.
