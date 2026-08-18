# Seurat GSE189357 — CLDN4-only malignant vs T/NK and IFN/MHC/TJ

ADDITIVE public slice. Zhu et al. 2022, GEO **GSE189357**.
**CLDN4 only.** No TACSTD2∩CLDN4 dual-high. No GSE148071.
Honest unit = **patient** (TD1–TD9). **n may be 9 — say so.**

Primary is **R + Seurat** `CreateSeuratObject` on the public processed 10x MTX.
If Seurat cannot install, stop.

1. `scripts/download.sh` — GEO processed tar only.
2. `scripts/analyze.R` — Seurat object, marker-malignant + T/NK gates, patient Spearman, plots, `FINDING.md`.

Done when `FINDING.md`, `results/tables/patient_units.tsv`, and the figures exist.
