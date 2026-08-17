# Pairwise Milo vs malignant CLDN4 (GSE131907 + GSE148071)

**Scope.** Additive CLDN4-only neighbourhood DA on two public NSCLC scRNA
series. Prior single-cohort TACSTD2 / CLDN4 Milo folders are taken as given
and are not re-run.

**Not this folder:** GSE205335, the triple merge, dual-high TACSTD2+CLDN4,
Harmony integration of the two matrices, miloR / edgeR.

## Why two graphs, not one

The matrices are not interchangeable. GSE131907 has author `Cell_type` /
`Cell_subtype`. GSE148071 has per-sample count TSVs and no public barcode
CNA IDs. Site is already a confounder inside GSE131907 (tLung ≠ mBrain).
A joint kNN would mix protocol and label systems and would not create a
larger honest *n*.

## Graph (each cohort)

1. HVG by raw-UMI dispersion (mean in [0.01, 50], top 2000) on that graph.
2. log1p(CP10k), gene-scale, PCA (`d=30`).
3. Subtract the per-sample PCA mean (**not** Harmony).
4. Euclidean kNN (`k=30`).
5. Milo refined index sampling (`prop=0.1`).
6. Neighbourhood = index ∪ its k nearest neighbours (size `k+1`).

## Malignant CLDN4 (CLDN4-only)

| Dataset | Malignant definition | Score |
|---|---|---|
| GSE131907 | author `Cell_subtype` in `{Malignant cells, tS1, tS2, tS3}` | mean log1p-CP10k CLDN4 in those cells |
| GSE148071 | epithelial AND NOT (alveolar/club/ciliated log1p-CP10k ≥ 1) | same formula |

Drop a sample if it has fewer than 10 malignant(-like) cells. TACSTD2 is
recorded only as a marker presence check, never as a gate.

## DA (fallback)

```
prop[s, i] = n_cells(s in i) / n_cells(s)
```

- **Primary:** Spearman of `prop` vs sample malignant CLDN4 (`min_samples=5`).
- **Secondary:** Welch t-test after a median split of scored samples.

SpatialFDR = miloR `graphSpatialFDR` k-distance weights. Neighbourhoods
are overlapping; SpatialFDR is overlap-aware. Do not cite cell count as *n*.

## Pairwise (sample is the unit)

The pair *n* is GSE148071 biopsies with a CLDN4 score **plus** GSE131907
**tLung** samples with a CLDN4 score. That is the only combined test whose
independent unit is the patient/sample.

GSE131907 mBrain is reported separately. It is not added to the pair *n*.

## Composition

Dissociated kNN ≠ histology.

1. Interface neighbourhoods (≥3 malignant and ≥3 T/NK).
2. Greedy disjoint subset (honest neighbourhood *n*).
3. Sample-paired T/NK in CLDN4-high vs CLDN4-low neighbourhoods.

(3) is the only composition test whose unit is the patient.
