# Honest TACSTD2 / CLDN4 results (no size skip)

All numbers below were computed from downloaded open files with
`scipy.stats.spearmanr` (two-sided) or recorded as **absent**.
`p = 0.0` from scipy underflow is written as `<1e-300`.
No statistic was invented for a missing gene.

Canonical table: `master_stats.csv`.

## What was downloaded (open only)

| DOI | What | Size | TACSTD2 | CLDN4 |
|-----|------|------|---------|-------|
| 10.5281/zenodo.10731914 | human 8-LUAD count matrix (already on disk) | 58 MB | yes | yes |
| 10.5281/zenodo.10731914 | mouse LUAD count matrix | 251 MB | Tacstd2 yes | Cldn4 yes (almost never detected) |
| 10.5281/zenodo.8417887 | 7 Visium h5 + malignant txt + **Immune_cell_data.txt 2.67 GB** | ~3.1 GB | yes | yes |
| 10.5281/zenodo.11205626 | 24samples.h5.tar.gz (Cell protocol used) | 761 MB | yes | yes |
| 10.5281/zenodo.13947395 | Fig6_processed_scRNAseq.zip (T-cell h5ad) | 1.5 GB | yes in var | yes in var |
| 10.5281/zenodo.2635194 | Nanostring NSCLC anti-PD1 | 0.3 MB | **no** | **no** |
| 10.64898/2026.01.16.25342913 | I3LUNG_DATA.zip | 25 MB | **no** | **no** |
| 10.5281/zenodo.8041882 | IMC SCE (already on disk) | 160 MB | **no** | **no** |

Restricted Zenodo records were still skipped.

## Computed ρ (Spearman TACSTD2 vs CLDN4)

**Where both genes are actually expressed (epithelial / spatial / malignant):**

| subset | n | ρ | p |
|--------|---|-----|---|
| 8-LUAD all cells (CP10K) | 28457 | 0.759 | <1e-300 |
| 8-LUAD epithelial only (CP10K) | 9517 | 0.388 | <1e-300 |
| early LUAD Visium, 7 samples pooled | 26992 spots | 0.717 | <1e-300 |
| early LUAD Visium per sample | 3235–4576 | 0.50–0.77 | all p ≤ 1.4e-254 |
| early LUAD malignant cells | 3321 | 0.398 | 3.33e-126 |
| paired LUAD Cell-protocol, all | 74023 | 0.587 | <1e-300 |
| paired LUAD Cell-protocol, tumor | 34290 | 0.611 | <1e-300 |
| paired LUAD Cell-protocol, normal | 39733 | 0.563 | <1e-300 |

Mann-Whitney (8-LUAD, CP10K, epithelial n=9517 vs rest n=18940): TACSTD2 p <1e-300; CLDN4 p <1e-300. Medians 3.91 vs 0 and 4.22 vs 0.

**Where genes are present but essentially off (do not over-read ρ):**

| subset | n | % TACSTD2+ | % CLDN4+ | ρ | p |
|--------|---|------------|----------|-----|---|
| early LUAD immune cells | 34544 | 2.83 | 2.40 | 0.045 | 5.23e-17 |
| mouse LUAD all cells | 63987 | 0.22 | 0.0016 | −0.00018 | 0.963 |
| Bischoff T cells | 14808 | 1.74 | 2.24 | 0.022 | 0.0082 |
| HTAN-MSK T cells | 26877 | 1.02 | 1.35 | 0.010 | 0.086 |
| Laughney T cells | 7042 | 1.36 | 1.24 | 0.042 | 0.00046 |

The immune / T-cell p-values can be small because n is large; the effect size is near zero and detection is ~1–3%. Mouse Cldn4 is one-in-sixty-thousand; ρ is undefined in three treatment groups where Cldn4 is all zeros.

**Absent — ρ not computed:**

- Nanostring nCounter NSCLC anti-PD1 (10.5281/zenodo.2635194): 784-gene immune panel, n=55+36 patients. TACSTD2/CLDN4/TROP2 not on the panel.
- I3LUNG_DATA (10.64898/2026.01.16.25342913): clinical / radiomics / KRAS-P53-STK11 only.
- Anti-PD1 NSCLC IMC (10.5281/zenodo.8041882): 41-protein panel, 99659 cells.

## Not downloaded (honest reasons)

- 10.5281/zenodo.7227571 NSCLC neutrophil atlas: 14–66 GB tarballs (scanVI model / containers / full atlas). Not a single streamable gene table on a 15 GB RAM machine.
- 10.5281/zenodo.21877227 LuCaS+ `Cell Ranger Outputs.zip` (4.6 GB): 454-patient snRNA-seq. The accompanying `Metadata_cell_annotation.csv` has no gene columns. Not pulled this round after ~6 GB of other matrices.
- 10.5281/zenodo.14594702 OS*.rds: pediatric osteosarcoma pulmonary metastases, not lung carcinoma ICI.
- 10.5281/zenodo.20539275: melanoma brain IMC, not lung.
- Restricted Zenodo (18488822, 18729311, 10911472, 7990870).

## How to reproduce

```bash
pip install pandas scipy h5py matplotlib rdata
python3 scripts/fable_zenodo/download_noskip.py
python3 scripts/fable_zenodo/stats_luad_10731914.py
python3 scripts/fable_zenodo/analyze_noskip.py
python3 scripts/fable_zenodo/analyze_tcell_h5ad.py
python3 scripts/fable_zenodo/write_master_table.py
```

Downloads are gitignored; checksums are in `download_manifest.json`.
