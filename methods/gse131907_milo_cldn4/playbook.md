# kNN neighbourhood DA vs malignant CLDN4 (GSE131907)

**Scope.** Additive test of whether transcriptional neighbourhoods change
with a sample-level **malignant CLDN4** score on the public Kim et al. 2020
LUAD atlas (GSE131907), and whether CLDN4-high neighbourhoods are T/NK-poor.

**Not miloR.** No R / Bioconductor / miloR / edgeR. `scripts/knn_nhood.py`
reimplements neighbourhood construction and **SpatialFDR** from Dann et al.,
*Nat Biotechnol* 2022 (miloR `graphSpatialFDR`, k-distance weights; Lun et
al. cydar). The DA model is a sample-level Spearman or Welch t-test on
neighbourhood proportions. Do not cite these p-values as miloR output.

GSE207422 Milo is a different agent. Do not run it from this folder.

## When this is the right tool

- Multiple samples, each with malignant and immune cells.
- DA that is not locked to discrete clusters.
- Sample, not cell, is the independent unit.

Wrong tool if you want spatial niches from dissociated scRNA, or the exact
edgeR QLF numbers from a published miloR run.

## Graph

1. Author cell labels (not marker reconstruction). Keep epithelium + immune.
2. HVG by raw-UMI dispersion on the cells in that graph (mean in [0.01, 50],
   top 2000).
3. log1p(CP10k), gene-scale, PCA (`d=30`).
4. Subtract the per-sample PCA mean (one-step batch centering; **not** Harmony).
5. Euclidean kNN (`k=30`).
6. Milo refined index sampling (`prop=0.1`).
7. Neighbourhood = index ∪ its k nearest neighbours (size `k+1`).

## Malignant CLDN4 (author labels)

Malignant = author `Cell_subtype` in `{Malignant cells, tS1, tS2, tS3}`.
tLung primaries are labeled tS1/tS2/tS3, not "Malignant cells". Metastases
and tL/B use "Malignant cells". nLung AT1/AT2/club/ciliated are not malignant.

Sample score = mean log1p-CP10k CLDN4 in that sample’s malignant cells.
Drop the sample if it has fewer than 10 malignant cells.

## DA (fallback)

```
prop[s, i] = n_cells(s in i) / n_cells(s)
```

- **Primary:** Spearman of `prop` vs sample malignant CLDN4.
- **Secondary:** Welch t-test after a median split of that score
  (`logFC = log2(mean_high) − log2(mean_low)`).

A neighbourhood is testable only if enough samples contribute cells.
Overlapping neighbourhoods are not independent observations.

## Cohorts (do not pool sites as if they were one n)

Site is a confounder: mBrain composition is not tLung composition.
Primary graph = **tLung** (same tissue). Secondary same-site graph = **mBrain**.
Do not quote a pooled SpatialFDR as the main result.

This atlas is treatment-naive. There is no ICI / MPR / RECIST label. Do not
invent an MPR contrast.

## SpatialFDR (honest)

`w = 1 / (k-th NN distance of the index cell)`

```
adjp[order] = rev(cummin(rev(sum(w) * p / cumsum(w))))
```

Report all of: `n_nhoods`, `n_testable`, **n samples**, raw p<0.05, BH-FDR,
SpatialFDR. Do not quote 200k cells as *n*.

## Composition

Dissociated kNN neighbourhoods are **transcriptional**, not spatial.

1. Interface neighbourhoods — ≥3 malignant and ≥3 T/NK cells. Spearman of
   neighbourhood malignant CLDN4 vs T/NK fraction.
2. Disjoint subset — greedy non-overlapping neighbourhoods (honest *n*).
3. Sample-paired — for each sample, T/NK fraction among that sample’s cells
   that fall in CLDN4-high vs CLDN4-low neighbourhoods. Wilcoxon signed-rank
   across samples.

(3) is the only composition test whose independent unit is the patient.

## Public files

| File | Role |
|---|---|
| `GSE131907_Lung_Cancer_cell_annotation.txt.gz` | Author Sample / Cell_type / Cell_subtype |
| `GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz` | Processed UMI (in budget) |
| `GSE131907_series_matrix.txt.gz` | Patient / stage / origin |

Skipped: 3 GB log2TPM text (same cells; UMI is the Milo input), RDS copies
(no R), EGA FASTQ (`EGAD00001005054`, controlled).
