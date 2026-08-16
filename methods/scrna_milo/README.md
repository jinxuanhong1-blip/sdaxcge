# Methods · neighborhood DA (Milo-equivalent) on public lung ICI scRNA

Additive high-end scRNA. The user thesis (malignant TACSTD2 / TROP2 as an ICI-relevant epithelial state next to a T/NK-poor milieu) is taken as given. This folder tests that claim at **neighborhood** resolution, not as a new sample-level correlation.

## Dataset (public only)

**GSE207422** (Hu et al., *Genome Medicine* 2023, PMID 36869384). Neoadjuvant PD-1 + chemo, resectable IIIA NSCLC. GEO processed UMI matrix: 24,292 genes × **92,330** cells, 15 samples.

| Paper scRNA group | Samples | n |
|---|---|---|
| Treatment-naïve biopsy (TN) | P01, P05, P08 | 3 |
| Post-treatment MPR (includes pCR P06) | P03, P06, P11, P14 | 4 |
| Post-treatment NMPR | P02, P04, P07, P09, P10, P12, P13, P15 | 8 |

P05/P08 are pre-treatment biopsies whose later surgical pathology is NMPR. They are **TN** for scRNA grouping (Hu Fig. 1), not NMPR.

**GSE241934** (NEOTIDE / EGFR-mutant + real-world neoadjuvant IO) was the listed alternative. It was not used: the public matrices are ~1.6 GB (IIT + real-world MTX), and GSE207422 already has epithelium + immune + MPR on a single in-budget public UMI matrix.

Author CopyKAT malignant barcodes are **not public** (GEO sample sheet only; author GitHub is scripts; no barcode table in the paper ESM). Malignant-like here is a marker proxy, not CopyKAT.

## What is run

1. Stream the GEO TSV (never materialize the uncompressed text).
2. Keep 2,500 HVGs ∪ lineage/target markers.
3. Assign lineage = argmax of mean log1p canonical markers (Hu set). **Malignant-like** = epithelial AND zero UMI for `SFTPA2`, `AGER`, `SCGB1A1`, `SCGB3A1`, `TPPP3`.
4. Restrict the graph to **epithelium + immune** (drop fibroblast / endothelial / other).
5. `log1p(CP10k)` → scale → PCA (30 PCs). Harmony on those PCs with `sample` as batch (`harmonypy`).
6. Combinatorial Milo on **both** embeddings (`X_pca` and `X_harmony`).
7. Patient-level malignant TACSTD2 = mean `log1p(CP10k)` in malignant-like cells. Median split among samples with ≥10 malignant-like cells (ties dropped). DA: high vs low.
8. MPR vs NMPR on post-treatment samples only, run because both arms have n≥3 (4 vs 8). Underpowered; reported anyway.
9. Keep neighborhoods that are **interface** (≥3 epithelial and ≥3 T/NK cells), **TACSTD2-high epithelium** (epithelial TACSTD2 ≥ median of interface nhoods), and **T/NK-depleted** (T/NK fraction ≤ 25th percentile of interface nhoods).

“Next to” is **transcriptomic KNN co-membership** in the epithelium+immune embedding. It is not physical spatial proximity.

## Milo / SpatialFDR (honest)

Neighborhoods follow milopy / miloR (Dann et al. 2022):

- Sample `prop=0.1` of cells on the KNN graph (`n_neighbors=30`).
- Refine each sampled neighborhood to the cell nearest the median embedding of its neighbors.
- Neighborhood = binary KNN of the refined index (index included).
- Count cells per **sample** in each neighborhood. The replicate is the sample, not the cell.

DA is a **Python negative-binomial GLM + quasi-likelihood F test** with `log(sample cell count)` as offset. This is a stand-in for `edgeR::glmQLFTest` as used by miloR. It is **not** the Bioconductor binary.

**SpatialFDR** is the published cydar / miloR weighted BH:

`weight_i = 1 / kth_neighbor_distance(index_i)`

`q_(i) = min_{j≥i} ( (∑w) p_(j) / ∑_{k≤j} w_(k) )` after sorting p-values.

Ordinary BH (`FDR_BH`) is written next to it so the two corrections can be compared. Cell-level p-values are not reported.

n is the number of **samples in each arm of that contrast**, plus the number of testable neighborhoods. Those numbers are in `results/scrna_milo/summary.json` → `da_counts`.

## Reproduce

```bash
python3 -m pip install -r methods/scrna_milo/requirements.txt
python3 methods/scrna_milo/00_download.py
python3 methods/scrna_milo/01_stream_matrix.py
python3 methods/scrna_milo/test_spatial_fdr.py
python3 methods/scrna_milo/02_run_milo.py
```

Outputs: `results/scrna_milo/`. Large matrices stay in `data/GSE207422/` (gitignored).

## Result of this public run

See `results/scrna_milo/FINDING.md`. Short version: **0 neighborhoods at SpatialFDR < 0.1 or < 0.2** (PCA and Harmony; TACSTD2 high vs low n=5 vs 5; MPR vs NMPR n=4 vs 8). Nominal P<0.05 is not empty. 24 compositional kept nhoods per embedding are written anyway.

## What this does not claim

- It does not re-open the sample-level A3 tests (those stay NS in prior public recomputes; this run’s sample Spearman is ρ=−0.22, p=0.47, n=13).
- A SpatialFDR hit at n=4 vs 8 MPR, or at a median split of 5 vs 5 TACSTD2 classes, would still be a neighborhood-level composition shift among few patients. None was observed.
- Harmony vs PCA disagreement is reported; replication across embeddings is not assumed.
