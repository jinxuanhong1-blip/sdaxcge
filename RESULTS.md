# CosMx NSCLC: IFN module and STING-pathway neighborhood scores

ADDITIVE, **CLDN4-only**, He et al. 2022 CosMx 960-plex (figshare 25976224, `cosmx_human_nsclc_clustered.h5ad`; 8 sections / 5 patients). This layer scores **IFN-module** and **STING-pathway** genes in the spatial neighborhood of malignant cells. It does **not** replace the locked exclusion result (fewer cytotoxic cells at 50/100 µm around CLDN4-high tumor; nearby effectors not muzzled on GZMB/PRF1/NKG7/IFNG). No private 8-KL. No ICI labels.

## Design

- **Malignant index:** published `cell_type` matched to the section's tumor cluster (`tumor 5` in LUAD-5, `tumor 6` in LUSC-6, `tumor 9` in LUAD-9, `tumor 12` in LUAD-12, `tumor 13` in LUAD-13). Cells carrying another patient's tumor label are not used as the index. Generic `epithelial` is not called malignant.
- **CLDN4-high / low:** Q4 vs Q1 of log1p(CP10k) CLDN4 among those malignant cells, **per section**. Q1 is 0 where CLDN4 is zero-inflated, so the low arm is larger than a strict quartile.
- **Neighborhood:** other cells within **50 µm** and **100 µm** of the malignant centroid (global coordinates, 0.18 µm/pixel). The index cell is excluded. FOV edges are not censored.
- **Module score:** each gene is log1p(CP10k), then z-scored within the section across all cells. The cell score is the mean z of the module. The neighborhood score is the mean cell score of neighbors. A companion column is the mean log1p(CP10k) of the same genes (expression units).
- **Contexts:** all neighbors (primary); immune neighbors only (published immune types — program among immune cells that are present); the malignant cell itself (intrinsic, not a neighborhood).
- **Epithelial control:** EPCAM, KRT8, KRT18, KRT19, KRT7, same neighborhood machinery.
- **Inference:** two-sided exact Wilcoxon signed-rank on the 8 section means, and on the 5 patient means (a patient's sections are averaged with equal section weight). Honest n = 8 sections / 5 patients.

## Gene sets and panel limits

IFN module (26): STAT1, JAK1, JAK2, IFNGR1, IFNGR2, IFNAR1, IFNAR2, IFNG, IFNB1, IFNA1, IFIT1, IFITM1, IFITM3, MX1, OAS1, OAS2, OAS3, OASL, CXCL9, CXCL10, CCL5, CD274, IDO1, BST2, TAP1, TAP2.

STING-pathway neighborhood (10): IRF3, NFKB1, RELA, NFKBIA, IFNB1, IFNA1, CXCL10, CCL5, TNF, IL6.

Confirmed **absent** from this 960-plex object: STING1, TMEM173, CGAS, MB21D1, TBK1, MAVS, TREX1. The STING score is the on-panel downstream neighborhood (IRF3, NF-κB, type-I IFN, CXCL10/CCL5), not a measurement of STING1 or cGAS.

HLA-A/B/C and B2M are on the panel and were **not** put in the IFN module (too constitutive). GZMB, PRF1, and NKG7 were **not** put in it either (those belong to the locked effector readout).

CXCL10, CCL5, IFNB1, and IFNA1 sit in both modules because they are IFN-response genes and STING outputs. The two scores are therefore not independent.

## Inventory

| Section | Patient | QC cells | Malignant | CLDN4-high | CLDN4-low | CLDN4>0 | Immune | Other tumor label |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| LUAD-5 R1 | Lung5 | 98,002 | 17,837 | 4,460 | 5,226 | 0.707 | 54,540 | 763 |
| LUAD-5 R2 | Lung5 | 100,335 | 17,914 | 4,485 | 5,537 | 0.691 | 57,741 | 1,196 |
| LUAD-5 R3 | Lung5 | 97,809 | 15,893 | 3,974 | 5,908 | 0.628 | 57,173 | 662 |
| LUSC-6 | Lung6 | 89,975 | 66,196 | 16,323 | 49,873 | 0.247 | 13,039 | 785 |
| LUAD-9 R1 | Lung9 | 87,606 | 38,864 | 9,750 | 11,030 | 0.716 | 26,220 | 558 |
| LUAD-9 R2 | Lung9 | 139,504 | 94,876 | 23,771 | 36,954 | 0.611 | 29,338 | 1,946 |
| LUAD-12 | Lung12 | 71,304 | 18,267 | 4,567 | 11,007 | 0.397 | 30,969 | 245 |
| LUAD-13 | Lung13 | 81,236 | 26,030 | 6,519 | 6,853 | 0.737 | 39,735 | 281 |

## Primary neighborhood scores (all neighbors)

Δ is CLDN4-high minus CLDN4-low. A positive Δ means a **higher** module score around CLDN4-high malignant cells. “Higher” counts sections (or patients) with Δ > 0.

| Module | Radius | Median section z high | Median section z low | Median Δ z | Sections higher | Section p | Patients higher | Patient p | Median Δ log1p |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ifn | 50 | 0.053 | 0.033 | 0.036 | 8/8 | 0.0078 | 5/5 | 0.0625 | 0.055 |
| ifn | 100 | 0.039 | 0.021 | 0.029 | 6/8 | 0.0391 | 5/5 | 0.0625 | 0.044 |
| sting | 50 | 0.050 | 0.010 | 0.033 | 8/8 | 0.0078 | 5/5 | 0.0625 | 0.047 |
| sting | 100 | 0.035 | 0.001 | 0.027 | 8/8 | 0.0078 | 5/5 | 0.0625 | 0.038 |

Per-section Δ z (all neighbors):

| Section | Patient | IFN 50 | IFN 100 | STING 50 | STING 100 |
|---|---|---:|---:|---:|---:|
| LUAD-5 R1 | Lung5 | 0.004 | -0.003 | 0.007 | 0.003 |
| LUAD-5 R2 | Lung5 | 0.002 | -0.001 | 0.009 | 0.005 |
| LUAD-5 R3 | Lung5 | 0.012 | 0.006 | 0.015 | 0.011 |
| LUSC-6 | Lung6 | 0.032 | 0.028 | 0.027 | 0.024 |
| LUAD-9 R1 | Lung9 | 0.065 | 0.051 | 0.071 | 0.051 |
| LUAD-9 R2 | Lung9 | 0.040 | 0.030 | 0.040 | 0.029 |
| LUAD-12 | Lung12 | 0.057 | 0.040 | 0.057 | 0.040 |
| LUAD-13 | Lung13 | 0.060 | 0.049 | 0.043 | 0.033 |

## Immune-neighbor-only and intrinsic scores

Immune-neighbor scores use only malignant cells that have at least one immune neighbor. Intrinsic is the malignant cell's own module score and is not a spatial neighborhood.

| Module | Context | Radius | Median Δ z | Sections higher | Section p | Patients higher | Patient p |
|---|---|---:|---:|---:|---:|---:|---:|
| ifn | all_neighbors | 50 | 0.036 | 8/8 | 0.0078 | 5/5 | 0.0625 |
| ifn | all_neighbors | 100 | 0.029 | 6/8 | 0.0391 | 5/5 | 0.0625 |
| ifn | immune_neighbors | 50 | 0.029 | 8/8 | 0.0078 | 5/5 | 0.0625 |
| ifn | immune_neighbors | 100 | 0.021 | 8/8 | 0.0078 | 5/5 | 0.0625 |
| ifn | intrinsic | — | 0.090 | 8/8 | 0.0078 | 5/5 | 0.0625 |
| sting | all_neighbors | 50 | 0.033 | 8/8 | 0.0078 | 5/5 | 0.0625 |
| sting | all_neighbors | 100 | 0.027 | 8/8 | 0.0078 | 5/5 | 0.0625 |
| sting | immune_neighbors | 50 | 0.030 | 8/8 | 0.0078 | 5/5 | 0.0625 |
| sting | immune_neighbors | 100 | 0.021 | 8/8 | 0.0078 | 5/5 | 0.0625 |
| sting | intrinsic | — | 0.089 | 8/8 | 0.0078 | 5/5 | 0.0625 |
| epithelial_control | all_neighbors | 50 | 0.115 | 5/8 | 0.1094 | 4/5 | 0.1250 |
| epithelial_control | all_neighbors | 100 | 0.083 | 5/8 | 0.2500 | 4/5 | 0.1875 |
| epithelial_control | immune_neighbors | 50 | 0.044 | 7/8 | 0.0156 | 5/5 | 0.0625 |
| epithelial_control | immune_neighbors | 100 | 0.026 | 7/8 | 0.0156 | 5/5 | 0.0625 |
| epithelial_control | intrinsic | — | 0.296 | 8/8 | 0.0078 | 5/5 | 0.0625 |

Immune-neighbor **counts** are not this layer's endpoint (the locked cytotoxic ratio already covers exclusion). They are recorded so a higher module score is not misread as “more immune cells.” Median immune neighbors at 50 µm, CLDN4-high vs low:

| Section | Median immune neighbors, high | Median immune neighbors, low | Fraction of high cells with ≥1 | Fraction of low cells with ≥1 |
|---|---:|---:|---:|---:|
| LUAD-5 R1 | 3.0 | 2.0 | 0.741 | 0.701 |
| LUAD-5 R2 | 3.0 | 2.0 | 0.738 | 0.724 |
| LUAD-5 R3 | 3.0 | 2.0 | 0.743 | 0.697 |
| LUSC-6 | 0.0 | 0.0 | 0.406 | 0.396 |
| LUAD-9 R1 | 1.0 | 2.0 | 0.595 | 0.715 |
| LUAD-9 R2 | 1.0 | 2.0 | 0.522 | 0.672 |
| LUAD-12 | 4.0 | 8.0 | 0.796 | 0.877 |
| LUAD-13 | 12.0 | 14.0 | 0.948 | 0.975 |

## Per-gene all-neighbor expression (50 µm)

Median section Δ of mean neighbor log1p(CP10k). Full 50 and 100 µm tests are in `tables/gene_neighborhood_tests.csv`.

| Module | Gene | Median Δ log1p | Sections higher | Section p | Patients higher | Patient p |
|---|---|---:|---:|---:|---:|---:|
| ifn | STAT1 | 0.0973 | 8/8 | 0.0078 | 5/5 | 0.0625 |
| ifn | JAK1 | 0.0701 | 5/8 | 0.1953 | 4/5 | 0.1250 |
| ifn | JAK2 | 0.0138 | 5/8 | 0.5469 | 4/5 | 0.4375 |
| ifn | IFNGR1 | 0.0562 | 8/8 | 0.0078 | 5/5 | 0.0625 |
| ifn | IFNGR2 | 0.0671 | 6/8 | 0.1484 | 4/5 | 0.1250 |
| ifn | IFNAR1 | 0.0264 | 5/8 | 0.1484 | 4/5 | 0.1250 |
| ifn | IFNAR2 | 0.0424 | 7/8 | 0.0234 | 5/5 | 0.0625 |
| ifn | IFNG | 0.0097 | 8/8 | 0.0078 | 5/5 | 0.0625 |
| ifn | IFNB1 | 0.0111 | 8/8 | 0.0078 | 5/5 | 0.0625 |
| ifn | IFNA1 | 0.0105 | 8/8 | 0.0078 | 5/5 | 0.0625 |
| ifn | IFIT1 | 0.0162 | 8/8 | 0.0078 | 5/5 | 0.0625 |
| ifn | IFITM1 | 0.0503 | 7/8 | 0.0156 | 5/5 | 0.0625 |
| ifn | IFITM3 | 0.1096 | 6/8 | 0.0391 | 4/5 | 0.1250 |
| ifn | MX1 | 0.1071 | 8/8 | 0.0078 | 5/5 | 0.0625 |
| ifn | OAS1 | 0.0239 | 6/8 | 0.0391 | 5/5 | 0.0625 |
| ifn | OAS2 | 0.0388 | 6/8 | 0.0391 | 4/5 | 0.1250 |
| ifn | OAS3 | 0.0584 | 7/8 | 0.0234 | 5/5 | 0.0625 |
| ifn | OASL | 0.0099 | 8/8 | 0.0078 | 5/5 | 0.0625 |
| ifn | CXCL9 | 0.0157 | 8/8 | 0.0078 | 5/5 | 0.0625 |
| ifn | CXCL10 | 0.0188 | 8/8 | 0.0078 | 5/5 | 0.0625 |
| ifn | CCL5 | 0.0227 | 8/8 | 0.0078 | 5/5 | 0.0625 |
| ifn | CD274 | 0.0139 | 8/8 | 0.0078 | 5/5 | 0.0625 |
| ifn | IDO1 | 0.0357 | 8/8 | 0.0078 | 5/5 | 0.0625 |
| ifn | BST2 | 0.0202 | 5/8 | 0.2500 | 4/5 | 0.1875 |
| ifn | TAP1 | 0.0426 | 7/8 | 0.0156 | 5/5 | 0.0625 |
| ifn | TAP2 | 0.0531 | 8/8 | 0.0078 | 5/5 | 0.0625 |
| sting | IRF3 | 0.0743 | 5/8 | 0.2500 | 4/5 | 0.1875 |
| sting | NFKB1 | 0.0289 | 6/8 | 0.0547 | 4/5 | 0.1250 |
| sting | RELA | 0.0945 | 6/8 | 0.0547 | 4/5 | 0.1250 |
| sting | NFKBIA | 0.1102 | 8/8 | 0.0078 | 5/5 | 0.0625 |
| sting | IFNB1 | 0.0111 | 8/8 | 0.0078 | 5/5 | 0.0625 |
| sting | IFNA1 | 0.0105 | 8/8 | 0.0078 | 5/5 | 0.0625 |
| sting | CXCL10 | 0.0188 | 8/8 | 0.0078 | 5/5 | 0.0625 |
| sting | CCL5 | 0.0227 | 8/8 | 0.0078 | 5/5 | 0.0625 |
| sting | TNF | 0.0232 | 7/8 | 0.0156 | 4/5 | 0.1250 |
| sting | IL6 | 0.0231 | 8/8 | 0.0078 | 5/5 | 0.0625 |

## Reading this layer

CLDN4-high neighborhoods score **higher**, not lower, on both modules. All-neighbor IFN Δ z is 0.036 at 50 µm (8/8 sections higher, p=0.0078; 5/5 patients, p=0.0625) and 0.029 at 100 µm (6/8 higher, p=0.0391). All-neighbor STING Δ z is 0.033 at 50 µm (8/8 higher, p=0.0078; 5/5 patients, p=0.0625) and 0.027 at 100 µm (8/8 higher, p=0.0078). The shift is small (a few hundredths of a within-section z). Patient p = 0.0625 is the smallest two-sided exact Wilcoxon p at n = 5 when every patient has the same sign.

The same direction is present when the average is restricted to immune neighbors (50 µm IFN Δ z 0.029, 8/8, p=0.0078; STING Δ z 0.030, 8/8, p=0.0078) and in the malignant cell itself (IFN Δ z 0.090, 8/8, p=0.0078; STING Δ z 0.089, 8/8, p=0.0078). At 50 µm, 26/26 IFN genes and 10/10 STING-pathway genes have a positive median section Δ. Section-consistent genes (8/8 higher) include STAT1, IFNGR1, MX1, IFIT1, CXCL9, CXCL10, CCL5, CD274, IDO1, NFKBIA, and IL6. JAK1, JAK2, and IRF3 are positive at the median but not consistent across sections. Epithelial-marker all-neighbor Δ z at 50 µm is larger at the median (0.115) but only 5/8 sections are higher (p=0.1094). The IFN and STING shifts are smaller and more consistent across sections. The malignant cell's own epithelial score is higher (8/8, median Δ z 0.296).

Lung5 is the flat section: IFN all-neighbor Δ z is about +0.002 to +0.012 at 50 µm, and two of the three Lung5 replicates flip slightly negative at 100 µm. Lung6 is immune-poor (median immune-neighbor count 0 in both arms); its immune-neighbor score uses the ~40% of malignant cells that have at least one immune neighbor. Lung9 and Lung12 have fewer immune neighbors around CLDN4-high than around CLDN4-low, and the IFN score among those neighbors is still higher.

## What this does not claim

- Not a new estimate of the locked 50/100 µm cytotoxic-cell ratio, and not a GZMB/PRF1/NKG7/IFNG muzzling test.
- Not evidence that STING1 or cGAS was measured. Those genes are off-panel.
- Not ICI response, not private 8-KL, and not a ligand-receptor probability.
- Q4 vs Q1 is a within-section contrast. n = 8 sections (5 patients), not 8 independent patients. Lung5 and Lung9 contribute repeated sections.

## Reproduce

```bash
python3 scripts/download_cosmx_nsclc_h5ad.py
python3 scripts/cosmx_cldn4_ifn_sting_neighborhood.py
```

The h5ad is gitignored. Tables are under `results/cosmx_cldn4_ifn_sting/tables/`. Figures: `figures/paired_ifn_neighborhood.png`, `figures/paired_sting_neighborhood.png`, `figures/delta_context_50um.png`, `figures/gene_delta_50um.png`.

