# CosMx NSCLC July-style neighbor retune (CLDN4-only)

Official NanoString/Bruker CosMx NSCLC FFPE 960-plex (He et al. 2022): **all 8 sections / 5 patients** (Lung5_Rep1/2/3, Lung6, Lung9_Rep1/2, Lung12, Lung13). Additive. CLDN4-only. No private 8-KL.

July PPT target metric family (not copied as results): CD8/NK neighbors 3.15→2.05, 8/8 sections P=0.008, 5/5 donors P=0.031, mixing 0.05/0.17/0.21. Primary readout is **mean immune-neighbor count**, not nearest-µm to CD8.

## Pre-specified primary

Declared before inspecting the grid:

- CLDN4 cut: **median-split among tumor cells**
- Radius: **40 µm**
- Immune: **CD8+NK** (CD8A+ / author T CD8; NKG7+ or author NK)
- Tumor gate: **author tumor if labels are present, else PanCK/KRT8/EPCAM high AND not CD8**
- Tumor cells are **not** required to have CD8A=0
- FOVs with <50 tumor cells or <20 CD8 dropped (no contrast)
- Unit of test: **section (n=8)** and **donor (n=5)** paired Wilcoxon on section-mean (donor-mean) neighbor counts
- Also report neighborhood fraction (CD8+NK / all neighbors) and neighbor counts residualized on local epithelial density

Companion (not primary): k=10 / k=20 nearest-neighbor immune counts.

## Status

Official 8-section flat files are being downloaded from NanoString S3. Numbers, grid table, and figures will be written after the analysis run. This file is the pre-specified design; it does not contain fabricated counts.

Reproduce:

```bash
python3 scripts/download_cosmx_nsclc.py
python3 scripts/cosmx_july_neighbor_retune.py
```
