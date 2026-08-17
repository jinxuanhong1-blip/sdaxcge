# kNN neighbourhood DA vs malignant CLDN4 (merged GSE131907 + GSE205335)

**Scope.** Additive CLDN4-only test of whether transcriptional neighbourhoods
change with **sample/patient malignant CLDN4** on the two public UMI matrices
that already appear in the PR #320 T/NK combo. That T/NK ρ is given and is
not re-audited here.

**Not miloR.** No R / Bioconductor / miloR / edgeR. `scripts/knn_nhood.py`
reimplements neighbourhoods + **SpatialFDR** (Dann et al. 2022; k-distance
weights). DA = sample/patient Spearman on neighbourhood proportions.

GSE207422 is out of scope (SpatialFDR-null at n=7). No dual-high score.

## Graph

**Primary: one graph per dataset** (allowed; 15 GB RAM cannot hold a dense
Harmony of ~76k + ~81k cells).

1. Author labels. Keep epithelium + immune (drop fibroblasts / endothelium
   / oligodendrocytes / undetermined).
2. HVG by raw-UMI dispersion on the cells in that graph (mean in [0.01, 50],
   top 2000).
3. log1p(CP10k), gene-scale, PCA (`d=30`).
4. Subtract the per-sample (GSE131907) or per-patient (GSE205335) PCA mean.
   This is one-step batch centering, **not** Harmony.
5. Euclidean kNN (`k=30`), Milo refined index sampling (`prop=0.1`).
6. Neighbourhood = index ∪ its k nearest neighbours (size `k+1`).

**Extra: Harmony-aligned joint graph** if both HVG matrices can be held
after a documented per-unit cell cap (default 1,200 cells / sample or
patient). Batch = dataset. Shared genes = intersection of the two HVG
sets that exist in both matrices. kNN is on Harmony-corrected PCs
(no second per-sample centering). If Harmony OOMs or the shared HVG set
is <500 genes, the joint graph is skipped and FINDING says so.

## Malignant CLDN4

| Dataset | Malignant | Unit | Score |
|---|---|---|---|
| GSE131907 | author `Cell_subtype` in {tS1, tS2, tS3, Malignant cells} | **sample** | mean log1p-CP10k CLDN4 in that sample’s malignant cells |
| GSE205335 | author `lineage.sub == Malignant cells` | **patient** (biopsies pooled) | same |

Drop the unit if it has fewer than 10 malignant cells.

GSE131907 graphs:

- `tLung` — same-tissue primaries (n samples ≈ 11; scored if ≥10 malignant).
- `tumor_no_brain` — tLung + tL/B + mLN. These are the author-malignant
  sites that enter the PR #320 GSE131907 n=21 slice. mBrain, nLung, nLN,
  and PE are excluded. Mixed site is documented; this is **not** one tissue.

GSE205335 graph: non-normal tissues only (same drop as the GSE205335-only
Milo). RECIST is not an MPR substitute and is not the primary contrast.

## DA

```
prop[s, i] = n_cells(s in i) / n_cells(s)
```

Primary: Spearman of `prop` vs unit malignant CLDN4 (`min_samples = 5`).
Secondary: Welch t-test after a median split of that score.

A neighbourhood is testable only if enough units contribute cells.
Overlapping neighbourhoods are not independent observations. Do not pool
the two per-dataset SpatialFDR tables into one *n*.

## SpatialFDR

`w = 1 / (k-th NN distance of the index cell)`

```
adjp[order] = rev(cummin(rev(sum(w) * p / cumsum(w))))
```

Report `n_nhoods`, `n_testable`, **n samples/patients**, raw p<0.05,
BH-FDR, SpatialFDR. Do not quote cell count as *n*.

## What is given

PR #320 author %pos vs T/NK on GSE131907+GSE205335 (k=2, N=43,
ρ=−0.479, p=0.00152) is **given**. This folder does not recompute or
re-rank that combinatorial T/NK pool. Composition plots here are extra
and are transcriptional kNN, not that audit.

## Public files (no >2 GB extras)

| File | Role |
|---|---|
| `GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz` | ~390 MB UMI text |
| `GSE131907_Lung_Cancer_cell_annotation.txt.gz` | author labels |
| `GSE205335_Lung_IO_UMI_matrix.rds.gz` | ~500 MB UMI dgCMatrix |
| `GSE205335_Lung_IO_CellIdentity.txt.gz` | author labels |

Skipped: 3 GB log2TPM, EGA FASTQ, GSE207422.
