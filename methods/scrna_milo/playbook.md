# kNN neighbourhood differential abundance (Milo fallback)

**Scope.** How to test whether transcriptional neighbourhoods change with pathologic response (MPR vs NMPR) or with a sample-level malignant TACSTD2 score, and whether TACSTD2-high neighbourhoods are T/NK-poor. Public processed UMI/MTX only.

**Not miloR.** This environment does not ship R / Bioconductor / miloR / edgeR. The scripts in `scripts/` reimplement neighbourhood construction and **SpatialFDR** from Dann et al., *Nat Biotechnol* 2022 (miloR `graphSpatialFDR`, `k-distance` weights; Lun et al. cydar). The DA *model* is **not** edgeR quasi-likelihood. It is a sample-level Welch t-test or Spearman on neighbourhood proportions. Do not cite these p-values as miloR output.

## When this is the right tool

- You have a processed count matrix with **multiple samples** in two conditions, or a sample-level continuous score.
- You want DA that is not locked to discrete clusters (Milo’s reason for existing).
- You are willing to treat **sample**, not cell, as the independent unit.

It is the wrong tool if you only have one sample, if you want spatial niches from dissociated scRNA, or if you need the exact edgeR QLF numbers from a published miloR run.

## Graph

1. HVG by raw-UMI dispersion (mean in [0.01, 50], top 2000).
2. log1p(CP10k), gene-scale, PCA (`d=30`).
3. Subtract the per-sample PCA mean (one-step batch centering; **not** Harmony).
4. Euclidean kNN (`k=30`).
5. Milo refined index sampling (`prop=0.1`): random seeds → centroid of each seed’s kNN → nearest member becomes the index.
6. Neighbourhood = index ∪ its k nearest neighbours (size `k+1`).

## DA (fallback)

For neighbourhood *i* and sample *s*:

```
prop[s, i] = n_cells(s in i) / n_cells(s)
```

- **MPR vs NMPR:** Welch t-test of `prop` across samples. `logFC = log2(mean_NMPR) − log2(mean_MPR)`.
- **Malignant TACSTD2:** Spearman of `prop` vs the sample’s mean log1p-CP10k TACSTD2 in malignant/epithelial cells (sample dropped if <10 such cells).

A neighbourhood is testable only if enough samples contribute cells (`min_samples`, `min_per_group`). Overlapping neighbourhoods are **not** independent observations.

## SpatialFDR (honest)

`w = 1 / (k-th NN distance of the index cell)`

Order p-values; carry weights; weighted BH as in miloR:

```
adjp[order] = rev(cummin(rev(sum(w) * p / cumsum(w))))
```

Report **all** of: `n_nhoods`, `n_testable`, `n` samples per arm, raw p<0.05, BH-FDR<0.1 / 0.05, SpatialFDR<0.1 / 0.05. Do not quote 90k cells as *n*.

Equal k-distances reduce SpatialFDR to ordinary BH (see `knn_nhood.self_test`).

## TACSTD2-high vs T/NK (taken as the user trend)

Dissociated kNN neighbourhoods are **transcriptional**, not spatial. Epithelium clusters with epithelium, so an unrestricted TACSTD2 vs T/NK correlation is partly geometry.

Primary composition tests in the runners:

1. **Interface neighbourhoods** — ≥3 malignant and ≥3 T/NK cells. Spearman of malignant TACSTD2 vs T/NK fraction.
2. **Disjoint subset** — greedy non-overlapping neighbourhoods (honest *n*).
3. **Sample-paired** — for each sample, T/NK fraction among that sample’s cells that fall in TACSTD2-high vs TACSTD2-low neighbourhoods (median split on neighbourhoods with ≥5 malignant cells). Wilcoxon signed-rank across samples.

(3) is the only composition test whose independent unit is the patient.

## Public datasets used here

| Dataset | Matrix | Labels | Response |
|---|---|---|---|
| GSE207422 (Hu 2023) | GEO UMI TSV, 92,330 cells | marker-reconstructed (author barcodes not public) | post-treatment MPR (incl. pCR) vs NMPR, n=4 vs 8 |
| GSE241934 IIT (Zhao 2024) | GEO MTX, ~79k cells | author `major_cell_type` | MPR vs non-MPR, n=4 vs 7 |
| GSE241934 RWC | GEO MTX, ~230k cells | author `major_cell_type` | MPR+pCR vs non-MPR |

Raw GSA-Human / SRA FASTQ are not used.

## Reproduce

```bash
pip install -r methods/scrna_milo/requirements.txt
python3 methods/scrna_milo/scripts/download.py --dataset GSE207422
python3 methods/scrna_milo/scripts/run_gse207422.py
python3 methods/scrna_milo/scripts/download.py --dataset GSE241934
python3 methods/scrna_milo/scripts/run_gse241934.py --cohorts IIT
# optional, ~1.2 GB MTX:
python3 methods/scrna_milo/scripts/run_gse241934.py --cohorts RWC
```

If miloR becomes available later, keep the same neighbourhoods and replace only the DA model with `testNhoods` (edgeR QLF), then recompute SpatialFDR. Do not mix the two p-value columns in one table without labelling the model.
