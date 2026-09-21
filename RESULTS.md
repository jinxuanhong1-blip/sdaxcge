# CosMx spatial domains: CLDN4 enrichment vs immune mixing

Public CosMx SMI NSCLC, figshare 25976224 (`cosmx_human_nsclc_clustered.h5ad`): **765,771 cells, 960 genes, 8 sections, 5 donors** (He et al. 2022). CLDN4-only. No private cohort.

This layer does **not** replace the locked cell-level exclusion result (CLDN4-high tumor has fewer nearby cytotoxic cells at 50/100 µm; 8/8 sections, 5/5 donors, one-sided sign P = 0.031). Raw CLDN4 domain enrichment tracks Shannon diversity and normalized mixing. A separate battery asks whether any of that immune-poor signal is CLDN4-specific after epithelial covariates.

**FINAL: epithelial-domain exclusion.** Joint EPCAM+KRT residuals, residual-scored domains, TACSTD2∩CLDN4 versus TACSTD2-only, alternate k, and FOV holdouts do not produce a CLDN4-unique immune-poor signal at the bar (≥7/8 sections and ≥4/5 donors, ≥6 sections in the test, and ≥2 domains in each arm). Strongest adequate unique contrast: domain_resid_joint_k6 / shannon (6/8 sections, 4/5 donors, median Δ -0.474). The locked cell-level ratios 0.36 / 0.52 are unchanged.

Domains fit after dropping CLDN4, EPCAM, TACSTD2, KRT7, KRT8, KRT18, and KRT19 are still immune-poor when split on raw CLDN4 (8/8 sections, 5/5 donors at 50 µm). Those domains still track EPCAM (median within-section Spearman of domain-mean CLDN4 vs EPCAM 0.726), so the remaining panel still carries the epithelial program. That result is not counted as CLDN4-specific.

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

Same smoothed embedding, KMeans k = 4, 6, 8, 10, and 12, each with its own HMRF pass. Deltas are still high − low raw CLDN4.

| k | Metric | Sections <0 | Wilcoxon p (high<low) | Donors <0 | Sign p one-sided | Median section Δ |
|---:|---|---:|---:|---:|---:|---:|
| 4 | mean_cytotoxic_50um | 8/8 | 0.0039 | 5/5 | 0.0312 | -0.585 |
| 6 | mean_cytotoxic_50um | 7/8 | 0.0078 | 5/5 | 0.0312 | -0.164 |
| 8 | mean_cytotoxic_50um | 8/8 | 0.0039 | 5/5 | 0.0312 | -0.422 |
| 10 | mean_cytotoxic_50um | 7/8 | 0.0078 | 4/5 | 0.1875 | -0.412 |
| 12 | mean_cytotoxic_50um | 6/8 | 0.0742 | 3/5 | 0.5000 | -0.109 |
| 4 | normalized_mixing | 8/8 | 0.0039 | 5/5 | 0.0312 | -0.686 |
| 6 | normalized_mixing | 7/8 | 0.0078 | 5/5 | 0.0312 | -0.276 |
| 8 | normalized_mixing | 8/8 | 0.0039 | 5/5 | 0.0312 | -0.439 |
| 10 | normalized_mixing | 8/8 | 0.0039 | 5/5 | 0.0312 | -0.441 |
| 12 | normalized_mixing | 7/8 | 0.0742 | 4/5 | 0.1875 | -0.169 |
| 4 | shannon | 8/8 | 0.0039 | 5/5 | 0.0312 | -0.720 |
| 6 | shannon | 8/8 | 0.0039 | 5/5 | 0.0312 | -0.591 |
| 8 | shannon | 8/8 | 0.0039 | 5/5 | 0.0312 | -0.546 |
| 10 | shannon | 8/8 | 0.0039 | 5/5 | 0.0312 | -0.599 |
| 12 | shannon | 7/8 | 0.0391 | 5/5 | 0.0312 | -0.500 |

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

## CLDN4 after epithelium

Tumor cells only. Joint residual is log1p(CLDN4) after an intercept and log1p(EPCAM)+log1p(KRT8)+log1p(KRT18), fit inside each section. The broader residual adds KRT7 and KRT19. Residual-scored domains are the same spaGCN/HMRF domains, median-split on the mean of that cell-level residual. Dual-high is the top half of domain-mean CLDN4 and the top half of domain-mean TACSTD2; the contrast is dual-high minus TACSTD2-high/CLDN4-low. Epithelial-held-out domains drop CLDN4, EPCAM, TACSTD2, KRT7, KRT8, KRT18, and KRT19 from the features, then split on raw CLDN4. EPCAM tertiles split raw CLDN4 inside an EPCAM band and do **not** remove keratins. A negative delta is immune-poor in the high arm. Domain contrasts enter the sign test only when each arm has at least two domains. The beyond-epithelium bar is ≥7/8 sections and ≥4/5 donors on a test with at least 6 sections.

| Contrast | Outcome | Sections <0 | Wilcoxon p (high<low) | Donors <0 | Sign p one-sided | Median Δ |
|---|---|---:|---:|---:|---:|---:|
| cell_resid_EPCAM_KRT7_8_18_19 | cytotoxic_100um | 4/8 | 0.4219 | 3/5 | 0.5000 | 0.011 |
| cell_resid_EPCAM_KRT7_8_18_19 | cytotoxic_50um | 5/8 | 0.1914 | 2/5 | 0.8125 | -0.016 |
| cell_resid_EPCAM_KRT7_8_18_19 | immune_neighbor_frac | 5/8 | 0.1914 | 2/5 | 0.8125 | -0.009 |
| cell_resid_EPCAM_KRT8_KRT18 | cytotoxic_100um | 4/8 | 0.3203 | 3/5 | 0.5000 | -0.003 |
| cell_resid_EPCAM_KRT8_KRT18 | cytotoxic_50um | 6/8 | 0.0977 | 3/5 | 0.5000 | -0.020 |
| cell_resid_EPCAM_KRT8_KRT18 | immune_neighbor_frac | 6/8 | 0.1250 | 3/5 | 0.5000 | -0.010 |
| domain_resid_joint_k10 | mean_cytotoxic_50um | 4/8 | 0.2305 | 3/5 | 0.5000 | -0.085 |
| domain_resid_joint_k10 | normalized_mixing | 5/8 | 0.3711 | 4/5 | 0.1875 | -0.152 |
| domain_resid_joint_k10 | shannon | 5/8 | 0.1250 | 4/5 | 0.1875 | -0.348 |
| domain_resid_joint_k12 | mean_cytotoxic_50um | 5/8 | 0.1914 | 4/5 | 0.1875 | -0.256 |
| domain_resid_joint_k12 | normalized_mixing | 5/8 | 0.3203 | 4/5 | 0.1875 | -0.352 |
| domain_resid_joint_k12 | shannon | 5/8 | 0.1914 | 4/5 | 0.1875 | -0.303 |
| domain_resid_joint_k4 | mean_cytotoxic_50um | 5/5 | 0.0312 | 4/4 | 0.0625 | -0.597 |
| domain_resid_joint_k4 | normalized_mixing | 5/5 | 0.0312 | 4/4 | 0.0625 | -0.558 |
| domain_resid_joint_k4 | shannon | 5/5 | 0.0312 | 4/4 | 0.0625 | -0.699 |
| domain_resid_joint_k6 | mean_cytotoxic_50um | 5/8 | 0.1562 | 4/5 | 0.1875 | -0.316 |
| domain_resid_joint_k6 | normalized_mixing | 5/8 | 0.1562 | 4/5 | 0.1875 | -0.517 |
| domain_resid_joint_k6 | shannon | 6/8 | 0.1250 | 4/5 | 0.1875 | -0.474 |
| domain_resid_joint_k8 | mean_cytotoxic_50um | 5/8 | 0.3711 | 4/5 | 0.1875 | -0.290 |
| domain_resid_joint_k8 | normalized_mixing | 5/8 | 0.5273 | 4/5 | 0.1875 | -0.366 |
| domain_resid_joint_k8 | shannon | 5/8 | 0.1250 | 4/5 | 0.1875 | -0.321 |
| domains_epithelial_genes_held_out_k8 | mean_cytotoxic_50um | 8/8 | 0.0039 | 5/5 | 0.0312 | -0.232 |
| domains_epithelial_genes_held_out_k8 | normalized_mixing | 8/8 | 0.0039 | 5/5 | 0.0312 | -0.317 |
| domains_epithelial_genes_held_out_k8 | shannon | 8/8 | 0.0039 | 5/5 | 0.0312 | -0.607 |
| dual_high_vs_TACSTD2_only_k8 | mean_cytotoxic_50um | 0/0 | NA | 0/0 | NA | NA |
| dual_high_vs_TACSTD2_only_k8 | normalized_mixing | 0/0 | NA | 0/0 | NA | NA |
| dual_high_vs_TACSTD2_only_k8 | shannon | 0/0 | NA | 0/0 | NA | NA |
| epcam_tertile1_raw_CLDN4 | cytotoxic_100um | 3/4 | 0.1250 | 2/3 | 0.5000 | -0.305 |
| epcam_tertile1_raw_CLDN4 | cytotoxic_50um | 4/4 | 0.0625 | 3/3 | 0.1250 | -0.111 |
| epcam_tertile1_raw_CLDN4 | immune_neighbor_frac | 4/4 | 0.0625 | 3/3 | 0.1250 | -0.077 |
| epcam_tertile2_raw_CLDN4 | cytotoxic_100um | 3/6 | 0.4219 | 2/3 | 0.5000 | 8.34e-04 |
| epcam_tertile2_raw_CLDN4 | cytotoxic_50um | 4/6 | 0.1562 | 2/3 | 0.5000 | -0.045 |
| epcam_tertile2_raw_CLDN4 | immune_neighbor_frac | 3/6 | 0.2812 | 2/3 | 0.5000 | -0.003 |
| epcam_tertile3_raw_CLDN4 | cytotoxic_100um | 4/7 | 0.4062 | 3/4 | 0.3125 | -0.299 |
| epcam_tertile3_raw_CLDN4 | cytotoxic_50um | 4/7 | 0.1484 | 3/4 | 0.3125 | -0.075 |
| epcam_tertile3_raw_CLDN4 | immune_neighbor_frac | 4/7 | 0.5312 | 3/4 | 0.3125 | -0.007 |

TACSTD2∩CLDN4 versus TACSTD2-only: 3 sections had a non-empty TACSTD2-high/CLDN4-low arm, and 0 of those had at least two domains in each arm. Where the TACSTD2-only arm existed it was a single domain, so this contrast does not enter the sign test.

FOV holdout of the joint cell residual, 50 µm cytotoxic counts: 140/212 tested FOVs (at least 20 residual-high and 20 residual-low tumor cells) have Δ < 0. Leave-one-FOV-out sign flips of the section residual Δ: 1/232.

## Run

```bash
python3 scripts/download_cosmx_figshare_25976224.py
python3 scripts/cosmx_spagcn_domain_mixing.py
```

Seed 20260921. k=8. Cells used: 765771.

