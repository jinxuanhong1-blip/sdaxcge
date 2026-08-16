# A3 GSE207422: public cell annotation hunt + malignant-only TACSTD2

**Verdict on author malignant / CopyKAT IDs: not found.**

A third-party public cell annotation **was** found (DRMref, not Hu/Zhang CopyKAT). Malignant-only TACSTD2 was recomputed on that label. Primary tests remain non-significant and in the same direction as the local-box epithelial result (NMPR higher; inverse vs T/NK, NS).

## What was missing

Hu et al., *Genome Medicine* 2023 (PMID 36869384) called malignant epithelium with CopyKAT on reclustered epithelial cells (clusters E0_DST, E3_PCNA, E4_TOP2A, E7_SERPINB9, and part of E1_KRT17). Those per-cell IDs are **not** on GEO, not in the paper supplements, not in the author GitHub repo, and not in figshare / Zenodo / TISCH2 / CELLxGENE / GSA processed files.

GEO `GSE207422_NSCLC_scRNAseq_metadata.xlsx` is **sample-level only** (15 rows: patient, MPR/NMPR, drug). The UMI matrix has 92,330 barcodes (`BD_immuneXX_<id>`) and no cell-type column.

Author GitHub `Junjie-Hu/NSCLC-immunotherapy` is scripts only. CopyKAT is run from local `epithelium.rds` / `stromal.rds` and writes `data_out/copykat_res.rds`, which was never deposited.

## What was found (not author CopyKAT)

DRMref (Liu et al., *NAR* 2024; https://ccsm.uth.edu/DRMref/) hosts annotated Seurat objects:

- `All_RData_after_Annotation/GSE207422_Tor_seurat_afterAnno.RDS` (8,690 cells; P02/P03/P04/P06/P07)
- `All_RData_after_Annotation/GSE207422_Sin_seurat_afterAnno.RDS` (22,187 cells; P09–P15)

These are **DRMref marker-based 16-type labels**, not CopyKAT aneuploid/diploid calls. There is no normal-epithelial class. All 12 **post-treatment** scRNA samples are covered (3 treatment-naive biopsies P01/P05/P08 are absent). After DRMref QC the union is 30,877 / 92,330 GEO cells. All 30,877 barcodes match the GEO UMI matrix.

Malignant cells: **2,051** (Tor 1,142 + Sin 909).

## Recompute (malignant-only)

TACSTD2 from GEO UMI matrix. Per-sample mean `log1p(CP10k)` in DRMref `Malignant cells`. pCR (P06) counted as MPR. Wilcoxon exact two-sided by enumerating C(12,4)=495 assignments. T/NK = CD8+ T + CD4+ T + NK.

| Test | Result |
|---|---|
| NMPR vs MPR, mean log1p(CP10k) TACSTD2 | NMPR 1.57 vs MPR 1.12; Δ=+0.45; U=25; **p=0.15** (n=8 vs 4) |
| NMPR vs MPR, mean log1p(UMI) TACSTD2 | NMPR 1.55 vs MPR 1.24; **p=0.57** |
| NMPR vs MPR, % TACSTD2+ malignant | 75% vs 66%; **p=0.68** |
| Spearman mean log1p(CP10k) vs T/NK fraction | **ρ=−0.49, p=0.21** (n=12) |
| Spearman % TACSTD2+ vs T/NK fraction | ρ=−0.64, p=0.048 (secondary metric) |

Sensitivity dropping P06 (only 15 malignant cells): NMPR vs MPR log1p(CP10k) p=0.048 (n=8 vs 3). That is a post-hoc filter on n=3 MPR and is **not** the primary result.

## Honest read vs the local-box epithelial mismatch

Local box (n=92,330, epithelial-not-CopyKAT): NMPR>MPR p=0.68; vs T/NK ρ=−0.28 NS.

DRMref malignant-only: same direction, still NS on the pre-specified tests (p=0.15 and ρ=−0.49 p=0.21). This does **not** convert the A3 mismatch into a significant malignant-only TACSTD2–MPR or TACSTD2–T/NK result. It also does **not** recover the authors’ CopyKAT malignant set. P06 has almost no DRMref-malignant cells (15), which is consistent with pCR residual disease being sparse, but that label is still not CopyKAT.

## Files

- `sample_level_malignant_TACSTD2.tsv` — per-patient malignant TACSTD2 and T/NK fractions
- `stats.tsv` — Wilcoxon / Spearman
- `celltype_counts_by_patient.tsv` — DRMref cell-type counts
- `drmref_cell_annotation.tsv.gz` — 30,877 cells with DRMref type + TACSTD2 UMI
- `hunt_log.md` — sources checked
- `geo_scRNAseq_sample_metadata.tsv` — GEO sample sheet (not cell-level)
