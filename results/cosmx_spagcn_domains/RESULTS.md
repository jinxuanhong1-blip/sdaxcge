# CosMx spatial domains: CLDN4 enrichment vs immune mixing

Public CosMx SMI NSCLC, figshare 25976224 (`cosmx_human_nsclc_clustered.h5ad`): **765,771 cells, 960 genes, 8 sections, 5 donors** (He et al. 2022). CLDN4-only. No private cohort.

This layer does **not** replace the locked cell-level exclusion result (CLDN4-high tumor has fewer nearby cytotoxic cells at 50/100 µm; 8/8 sections, 5/5 donors, one-sided sign P = 0.031). It asks whether that exclusion is visible as spaGCN-style spatial domains whose CLDN4 enrichment tracks Shannon diversity and a normalized mixing index.

## Method

Per section, histology-free SpaGCN-style aggregation (Hu et al. 2021) on a sparse graph, because a dense SpaGCN adjacency is not tractable at ~10^5 cells:

- Coordinates: object pixels × 0.18 µm/px (calibrated to the public Giotto locs; median nearest centroid = 9.0 µm on Lung5_Rep1 FOV 1).
- Expression: raw `counts`, library-size normalize to 10^4, log1p. Up to 400 highest-variance genes detected in ≥1% of cells. **CLDN4 is removed from these features** before clustering.
- Spatial graph: 10 nearest neighbors, Gaussian weights with length scale l = median distance to the 10th neighbor, row-normalized, self included. One graph convolution, then PCA (30 components).
- Domains: KMeans k=8 (fixed so every section has the same domain count; SpaGCN's original Louvain step is the piece replaced here), seed 20260921.
- Giotto-style Potts HMRF refinement: Jacobi updates, emission = squared PCA distance / initial within-cluster median, Potts term = beta × fraction of the 10 neighbors sharing the label. Primary beta=2.0. If neighbor agreement stays below 0.55, beta is raised to 4.0 (contiguity QC, not an outcome search).
- CLDN4-high tumor cells, for the cell-level bridge only: within-section tumor CLDN4 count ≥ the tumor 75th percentile and ≥ 1. Domain split is separate: eligible domains (n≥100 and n_tumor≥30) are median-split on **mean CLDN4 of tumor cells in the domain**.
- Shannon: H = −Σ p log2 p on four compartments (tumor, immune, stromal, other). Evenness = H / 2.
- Normalized mixing: among tumor cells in the domain, the mean fraction of the 10 spatial neighbors that are immune, divided by the **section** immune fraction. Values < 1 mean tumor cells in that domain see fewer immune neighbors than the section base rate.
- Exclusion link: mean number of author cytotoxic cells (T CD8 naive, T CD8 memory, NK) inside 50 µm and 100 µm of tumor cells. Sample tallies use tumor-cell-weighted means of domains. Donor tallies average the sample values (Lung5 = 3 sections, Lung9 = 2).
- Tests: Wilcoxon signed-rank on 8 section deltas (high − low), and a sign test on 5 donor deltas. One-sided tests use the exclusion direction (high < low for Shannon, mixing, immune fraction, and cytotoxic counts). Two-sided p-values are reported next to them.
- Epithelial controls (EPCAM, KRT8, KRT18), same tumor-cell mean vs immune-neighbor association, are descriptive. They are not a claim that CLDN4 is the strongest surface gene.

Author `niche` labels already in the object (tumor interior, stroma, immune, …) are the published Giotto niches. Those niches were built from neighborhood composition, so immune fraction is partly by construction. Tumor-cell mean CLDN4 is also similar across tumor-interior and immune niches, because the immune niche still contains some CLDN4-high tumor cells. Niches are descriptive, not the primary test.

## Cell-level bridge to exclusion

Recomputed on this object (author cell types, cytotoxic = CD8 T + NK). Ratio = mean neighbor count around CLDN4-high tumor / mean around the remaining tumor. Q75 is the primary cut; the median cut is shown beside it.

| Section | Donor | 50 µm Q75 high | 50 µm Q75 low | ratio 50 | ratio 100 | ratio 50 median |
|---|---|---:|---:|---:|---:|---:|
| Lung12 | Lung12 | 0.263 | 0.569 | 0.462 | 0.636 | 0.462 |
| Lung13 | Lung13 | 1.544 | 1.680 | 0.919 | 0.977 | 0.928 |
| Lung5_Rep1 | Lung5 | 0.232 | 0.214 | 1.082 | 1.189 | 0.806 |
| Lung5_Rep2 | Lung5 | 0.244 | 0.237 | 1.033 | 1.214 | 0.696 |
| Lung5_Rep3 | Lung5 | 0.245 | 0.232 | 1.055 | 1.197 | 0.857 |
| Lung6 | Lung6 | 0.153 | 0.141 | 1.083 | 1.069 | 1.083 |
| Lung9_Rep1 | Lung9 | 0.376 | 0.632 | 0.596 | 0.753 | 0.525 |
| Lung9_Rep2 | Lung9 | 0.233 | 0.541 | 0.431 | 0.580 | 0.444 |

Q75 count-ratio < 1 in **4/8** sections at 50 µm and **4/8** at 100 µm; donor means < 1 in **3/5** (50 µm) and **3/5** (100 µm). Median section ratios: 50 µm 0.976, 100 µm 1.023. The locked 0.36 / 0.52 figures are the prior cell-level result and are not overwritten by these recomputed ratios.

## Domain result

Eligible domains: 61 of 64. Median Gaussian length scale 22.800 µm. Median post-HMRF neighbor agreement 0.849.

Section tallies (tumor-cell-weighted means of CLDN4-high vs CLDN4-low domains):

| Section | Donor | Shannon high | Shannon low | Δ Shannon | Mixing high | Mixing low | Δ mixing | Cyt50 high | Cyt50 low |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Lung12 | Lung12 | 1.010 | 1.495 | -0.484 | 0.294 | 0.973 | -0.679 | 0.221 | 0.875 |
| Lung13 | Lung13 | 0.570 | 1.054 | -0.484 | 0.353 | 0.551 | -0.199 | 1.578 | 1.795 |
| Lung5_Rep1 | Lung5 | 0.289 | 1.055 | -0.766 | 0.088 | 0.268 | -0.180 | 0.183 | 0.284 |
| Lung5_Rep2 | Lung5 | 0.560 | 1.074 | -0.514 | 0.111 | 0.968 | -0.857 | 0.201 | 0.828 |
| Lung5_Rep3 | Lung5 | 0.385 | 0.963 | -0.578 | 0.100 | 0.236 | -0.136 | 0.217 | 0.262 |
| Lung6 | Lung6 | 0.260 | 0.475 | -0.215 | 0.133 | 0.271 | -0.139 | 0.109 | 0.174 |
| Lung9_Rep1 | Lung9 | 0.315 | 1.330 | -1.015 | 0.100 | 1.197 | -1.097 | 0.426 | 1.784 |
| Lung9_Rep2 | Lung9 | 0.164 | 0.902 | -0.738 | 0.082 | 0.844 | -0.762 | 0.192 | 0.947 |

Donor tallies (unweighted mean of that donor's sections):

| Donor | Sections | Δ Shannon | Δ mixing | Δ immune frac | Δ cyt 50 µm | Δ cyt 100 µm |
|---|---:|---:|---:|---:|---:|---:|
| Lung12 | 1 | -0.484 | -0.679 | -0.379 | -0.654 | -1.853 |
| Lung13 | 1 | -0.484 | -0.199 | -0.179 | -0.217 | 0.048 |
| Lung5 | 3 | -0.620 | -0.391 | -0.308 | -0.258 | -0.774 |
| Lung6 | 1 | -0.215 | -0.139 | -0.030 | -0.065 | -0.372 |
| Lung9 | 2 | -0.877 | -0.930 | -0.376 | -1.056 | -3.124 |

### Tests (high − low; exclusion direction is negative)

| Metric | Sections <0 | Wilcoxon p two-sided | Wilcoxon p (high<low) | Donors <0 | Sign p one-sided | Sign p two-sided |
|---|---:|---:|---:|---:|---:|---:|
| Shannon | 8/8 | 0.0078 | 0.0039 | 5/5 | 0.0312 | 0.0625 |
| Normalized mixing | 8/8 | 0.0078 | 0.0039 | 5/5 | 0.0312 | 0.0625 |
| Immune fraction | 8/8 | 0.0078 | 0.0039 | 5/5 | 0.0312 | 0.0625 |
| Cytotoxic count 50 µm | 8/8 | 0.0078 | 0.0039 | 5/5 | 0.0312 | 0.0625 |
| Cytotoxic count 100 µm | 6/8 | 0.0391 | 0.0195 | 4/5 | 0.1875 | 0.3750 |

Within-section Spearman (eligible domains, tumor-mean CLDN4 vs metric):

| Metric | Sections with ρ<0 | Median ρ |
|---|---:|---:|
| Shannon | 6/8 | -0.452 |
| Normalized mixing | 5/8 | -0.310 |
| Immune fraction | 6/8 | -0.405 |
| Cytotoxic count 50 µm | 6/8 | -0.173 |

CLDN4 rank residualized on EPCAM within each section (epithelial-program control). A negative residual ρ would mean CLDN4 still tracks that metric after the shared epithelial axis is removed.

| Residual CLDN4 vs | Sections with ρ<0 | Median ρ |
|---|---:|---:|
| Shannon | 3/8 | 0.060 |
| Normalized mixing | 2/8 | 0.248 |
| Cytotoxic count 50 µm | 4/8 | 0.028 |

Epithelial control genes, within-section Spearman of the tumor-cell mean against normalized mixing:

| Gene | Sections with ρ<0 | Median ρ |
|---|---:|---:|
| EPCAM | 8/8 | -0.708 |
| KRT8 | 8/8 | -0.696 |
| KRT18 | 8/8 | -0.768 |

## Sensitivity to domain count

Same smoothed embedding, KMeans k = 6, 8, and 10, each with its own HMRF pass. Deltas are still high − low CLDN4.

| k | Metric | Sections <0 | Wilcoxon p (high<low) | Donors <0 | Sign p one-sided | Median section Δ |
|---:|---|---:|---:|---:|---:|---:|
| 6 | mean_cytotoxic_50um | 7/8 | 0.0078 | 5/5 | 0.0312 | -0.164 |
| 8 | mean_cytotoxic_50um | 8/8 | 0.0039 | 5/5 | 0.0312 | -0.422 |
| 10 | mean_cytotoxic_50um | 7/8 | 0.0078 | 4/5 | 0.1875 | -0.412 |
| 6 | normalized_mixing | 7/8 | 0.0078 | 5/5 | 0.0312 | -0.276 |
| 8 | normalized_mixing | 8/8 | 0.0039 | 5/5 | 0.0312 | -0.439 |
| 10 | normalized_mixing | 8/8 | 0.0039 | 5/5 | 0.0312 | -0.441 |
| 6 | shannon | 8/8 | 0.0039 | 5/5 | 0.0312 | -0.591 |
| 8 | shannon | 8/8 | 0.0039 | 5/5 | 0.0312 | -0.546 |
| 10 | shannon | 8/8 | 0.0039 | 5/5 | 0.0312 | -0.599 |

## How this sits with exclusion

Two different aggregations are in this file. Cell-level Q75 (each tumor cell scored by its own CLDN4 count) gives a cytotoxic-count ratio < 1 in 4/8 sections at 50 µm and 4/8 at 100 µm. That estimator is not the locked 0.36 / 0.52 result, and these recomputed ratios do not replace it. Domain-level aggregation puts tumor cells into spaGCN/HMRF domains and compares CLDN4-high vs CLDN4-low domains. Shannon is lower in the CLDN4-high domains in 8/8 sections and 5/5 donors; normalized mixing in 8/8 sections and 5/5 donors; 50 µm cytotoxic counts in 8/8 sections and 5/5 donors; 100 µm cytotoxic counts in 6/8 sections and 4/5 donors.

The domain contrast follows the epithelial program, not a CLDN4-only axis. EPCAM, KRT8, and KRT18 tumor-means correlate with lower normalized mixing at least as consistently as CLDN4, and the CLDN4 residual after EPCAM does not stay negative. Read the domain result as: CLDN4-enriched epithelial domains are the immune-poor domains (exclusion geography). Do not read it as public evidence that CLDN4 outranks other epithelial genes. Effector-cell state (GZMB/PRF1/NKG7/IFNG) is not retested here.

## Author niches (descriptive)

| Niche | Sections | Mean tumor CLDN4 | Shannon | Normalized mixing | Cytotoxic 50 µm | Immune fraction |
|---|---:|---:|---:|---:|---:|---:|
| immune | 8 | 1.672 | 1.450 | 1.182 | 2.452 | 0.605 |
| tumor interior | 8 | 1.615 | 0.882 | 0.131 | 0.293 | 0.098 |
| tumor-stroma boundary | 8 | 1.591 | 1.699 | 0.449 | 0.736 | 0.303 |
| lymphoid structure | 7 | 1.545 | 1.059 | 1.185 | 3.253 | 0.729 |
| myeloid-enriched stroma | 8 | 1.538 | 1.283 | 1.842 | 2.068 | 0.694 |
| plasmablast-enriched stroma | 6 | 1.423 | 1.248 | 0.909 | 0.889 | 0.686 |
| stroma | 8 | 1.386 | 1.396 | 1.115 | 1.830 | 0.526 |
| neutrophils | 7 | 1.277 | 1.345 | 1.325 | 1.023 | 0.645 |
| macrophages | 5 | 0.994 | 0.773 | 2.032 | 0.773 | 0.816 |

## Run

```bash
python3 scripts/download_cosmx_figshare_25976224.py
python3 scripts/cosmx_spagcn_domain_mixing.py
```

Seed 20260921. k=8. Cells used: 765771.

