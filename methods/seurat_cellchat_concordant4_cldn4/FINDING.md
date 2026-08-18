# FINDING — Seurat + CellChat concordant-four CLDN4-high vs low senders

ADDITIVE. **Thesis already correct. Ligands stay.**
CLDN4 only. No dual-high. Concordant four only
(GSE123902 + GSE131907 + GSE205335 + GSE189357).
Do **not** add GSE148071 / GSE127465 / CD45-only. This is **not** a full-pool.

Engine: **R + Seurat + CellChat** (not a Python reimplementation).
Seurat 5.5.1. CellChat 2.2.0.9001.
Outgoing communication probability is CellChat `computeCommunProb`
(truncatedMean, trim=0.1, population.size=TRUE).
Senders = malignant CLDN4-high vs CLDN4-low; receivers = T/NK.
Honest n = patient / locked sample.

Thesis:

- Barrier/inhibitory outgoing **UP from CLDN4-high** (F11R, NECTIN2–TIGIT, CDH1, LGALS9) is **ON-thesis**.
- IFN/T-recruit outgoing UP from CLDN4-low is the KD-like arm (often **not** detected; CXCL9/10 sparse).
  Do not bury barrier-up-in-high as a recruit-up skip.

Primary split is malignant **Q4 vs Q1**. Extra: median and %pos.

## Honest n

Locked four n=65 (13+21+22+9) is **not** the ligand n.

| gate | n | note |
|---|---:|---|
| Locked four | 65 | 13+21+22+9 |
| Inventory units loaded | 65 | GSE123902=13, GSE131907=21, GSE189357=9, GSE205335=22 |
| Both compartments (n_mal≥10, n_tnk≥20) | 65 | honest inventory |
| Q4 vs Q1 CellChat units | **61** | primary ligand n |

GSE123902 / GSE189357 malignant labels are thin → epithelium marker-malignant
(EPCAM|KRT8|KRT18|KRT19 > 0 and PTPRC == 0) → T/NK.
GSE131907 / GSE205335 use author malignant → T/NK.
TACSTD2 is never a gate. PVRL2 is aliased to NECTIN2 on GSE131907.

## Primary — barrier / inhibitory family ΔP (high − low)

ON-thesis. Primary table: `results/tables/ligand_table.tsv` (Q4 vs Q1).

| pair | axis | expect | n | 123902 | 131907 | 205335 | 189357 | mean ΔP | p_W | observed | agrees |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|---|
| JAM1_ITGAL_ITGB2 | F11R | high>low | 43 | 6 | 13 | 15 | 9 | +0.005 | 1.54e-08 | high>low | yes |
| NECTIN2_TIGIT | NECTIN2-TIGIT | high>low | 54 | 9 | 19 | 18 | 8 | +0.003 | 5.37e-10 | high>low | yes |
| CDH1_ITGAE_ITGB7 | CDH1 | high>low | 43 | 7 | 15 | 17 | 4 | +0.002 | 1.16e-08 | high>low | yes |
| CDH1_KLRG1 | CDH1 | high>low | 33 | 3 | 10 | 16 | 4 | +0.002 | 5.64e-07 | high>low | yes |
| LGALS9_HAVCR2 | LGALS9 | high>low | 23 | 6 | 8 | 8 | 1 | +0.001 | 0.00711 | high>low | yes |
| LGALS9_CD44 | LGALS9 | high>low | 43 | 10 | 13 | 16 | 4 | +0.004 | 8.60e-06 | high>low | yes |
| LGALS9_CD45 | LGALS9 | high>low | 44 | 10 | 14 | 16 | 4 | +0.005 | 1.70e-05 | high>low | yes |
| FAMILY_barrier_inhibitory | family | high>low | 59 | 12 | 19 | 19 | 9 | +0.003 | 1.01e-10 | high>low | yes |

Family aggregate: n=59 mean ΔP=+0.003 p_W=1.01e-10 agrees=yes.

Official CellChat `computeCommunProb` probabilities are on a smaller
absolute scale than the prior Python Hill reimplementation (that folder
reported family mean ΔP ≈ +0.099). Direction and Wilcoxon sign agree.
This run is the R CellChat primary.

## KD-like arm — IFN / T-recruit / MHC-I (expect low > high)

Do not file the barrier result as a recruit-up skip.

| pair | axis | expect | n | 123902 | 131907 | 205335 | 189357 | mean ΔP | p_W | observed | agrees |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|---|
| CXCL9_CXCR3 | CXCL9/10-CXCR3 | low>high | 1 | 1 | 0 | 0 | 0 | +0.000 | NA | high>low | opposite |
| CXCL10_CXCR3 | CXCL9/10-CXCR3 | low>high | 7 | 1 | 4 | 2 | 0 | -0.001 | 0.554 | low>high | yes |
| CCL5_CCR5 | CCL5 | low>high | 2 | 0 | 1 | 1 | 0 | -0.000 | 1 | low>high | yes |
| CCL5_CCR1 | CCL5 | low>high | 2 | 0 | 0 | 2 | 0 | -0.000 | 0.371 | low>high | yes |
| HLA-A_CD8A | HLA-CD8 | low>high | 46 | 8 | 17 | 14 | 7 | +0.002 | 0.0072 | high>low | opposite |
| HLA-B_CD8A | HLA-CD8 | low>high | 47 | 7 | 18 | 16 | 6 | +0.003 | 8.83e-05 | high>low | opposite |
| HLA-C_CD8A | HLA-CD8 | low>high | 43 | 8 | 16 | 14 | 5 | +0.003 | 3.43e-04 | high>low | opposite |
| FAMILY_ifn_recruit_mhci | family | low>high | 49 | 8 | 18 | 16 | 7 | +0.002 | 7.22e-05 | high>low | opposite |

## Extra: median-split barrier family

| pair | axis | expect | n | 123902 | 131907 | 205335 | 189357 | mean ΔP | p_W | observed | agrees |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|---|
| JAM1_ITGAL_ITGB2 | F11R | high>low | 42 | 5 | 12 | 16 | 9 | +0.004 | 1.71e-08 | high>low | yes |
| NECTIN2_TIGIT | NECTIN2-TIGIT | high>low | 52 | 8 | 17 | 19 | 8 | +0.002 | 8.87e-09 | high>low | yes |
| CDH1_ITGAE_ITGB7 | CDH1 | high>low | 42 | 6 | 15 | 17 | 4 | +0.002 | 4.99e-08 | high>low | yes |
| CDH1_KLRG1 | CDH1 | high>low | 33 | 2 | 10 | 17 | 4 | +0.001 | 5.64e-07 | high>low | yes |
| LGALS9_HAVCR2 | LGALS9 | high>low | 22 | 5 | 8 | 8 | 1 | +0.001 | 0.00705 | high>low | yes |
| LGALS9_CD44 | LGALS9 | high>low | 41 | 7 | 13 | 17 | 4 | +0.003 | 1.32e-04 | high>low | yes |
| LGALS9_CD45 | LGALS9 | high>low | 44 | 8 | 14 | 17 | 5 | +0.004 | 2.65e-04 | high>low | yes |
| FAMILY_barrier_inhibitory | family | high>low | 57 | 10 | 18 | 20 | 9 | +0.002 | 1.06e-08 | high>low | yes |

## What is not claimed

- TACSTD2 is not used to define high/low. This is not dual-high.
- GSE148071, GSE127465, and CD45-only libraries are not added.
- This is not a CellChat discovery screen and not a 7-pool.
- Cell-pooled tests are not reported. Honest n is the patient/sample.

## Reproduce

```bash
Rscript methods/seurat_cellchat_concordant4_cldn4/scripts/install_packages.R
bash methods/seurat_cellchat_concordant4_cldn4/scripts/download.sh /tmp/concordant4_raw
Rscript methods/seurat_cellchat_concordant4_cldn4/scripts/run_cellchat.R --raw=/tmp/concordant4_raw
```

CellChat: type=truncatedMean, trim=0.1, population.size=TRUE.
Arm floor ≥10; T/NK ≥20; Q4 vs Q1 also n_mal≥40.

