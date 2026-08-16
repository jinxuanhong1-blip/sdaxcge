# W200-A3 · GSE207422 — malignant-only TACSTD2 (NMPR vs MPR, vs T/NK)

**Task:** Recompute Hu et al. 2023 malignant TACSTD2: NMPR vs MPR, and vs T/NK. Honest.

## Verdict

**Neither half of A3 is supported** on this public matrix.

- **NMPR > MPR malignant TACSTD2:** no. Sample-level MWU p=0.89 (n=8 vs 2). Two of four MPR samples (P11, P14) have **zero** malignant-like cells — all their epithelium expresses normal-lung markers, matching the paper’s MPR normal-epithelium expansion. All-epithelial TACSTD2 is also null (p=0.68, n=8 vs 4).
- **Malignant TACSTD2 vs T/NK anti-correlation (user ρ −0.40 to −0.50):** no. Post-eligible ρ=+0.18 (p=0.70, n=7); all-eligible ρ=−0.14 (p=0.70, n=10); all 15 samples epithelial ρ=−0.13 (p=0.65). Closest planned |ρ| is 0.30 (p=0.40).
- **TACSTD2 is epithelial-restricted:** median %pos 74.9% malignant-like vs 1.6% T/NK (paired Wilcoxon p=9.8e-4, n=10).
- A rejected cell-level Spearman on all 92,330 cells (TACSTD2 vs T/NK score) is ρ=−0.16, not −0.45. That test is compositional and is not the A3 claim.

## Data policy

- Processed UMI matrix **175.5 MB gzip < 2 GB** → recomputed (not catalog-only).
- 92,330 cells (header count matches the paper's post-QC 92,330).
- Raw GSA-Human HRA001033 was **not** downloaded.
- Author per-cell labels are **not public** (GEO sample sheet only; paper ESM has no barcode table; [Junjie-Hu/NSCLC-immunotherapy](https://github.com/Junjie-Hu/NSCLC-immunotherapy) has scripts, not RDS/metadata_scrublet.csv; TISCH2 NSCLC_GSE207422 CellMetainfo 404). See `label_search.json`.

## What was used instead of author labels

- Lineage = argmax of mean log1p marker scores (Hu canonical: EPCAM/KRTs, CD3D/E/CD2, NKG7/GNLY/FGFBP2, etc.).
- **Malignant-like** = epithelial AND zero UMI for normal-lung markers SFTPA2, AGER, SCGB1A1, SCGB3A1, TPPP3.
- This is **not** author CopyKAT aneuploid calls. Residual basal/unmarked normal epithelium can leak in.
- T/NK = T or NK lineage. Eligible sample: ≥10 malignant-like and ≥20 T/NK cells.
- Paper groups: TN = 3 pre-biopsies; post MPR includes pCR P06 (n=4); post NMPR n=8.
- Marker-based counts: epithelial=11019, malignant-like=6627, T/NK=36323.

## Results

- all eligible samples: malig-like TACSTD2 mean_log1p vs T/NK fraction: ρ=-0.139, p=0.701, n=10
- all eligible samples: malig-like TACSTD2 %pos vs T/NK fraction: ρ=-0.103, p=0.777, n=10
- all eligible samples: malig-like TACSTD2 mean log1p(CP10k) vs T/NK fraction: ρ=-0.042, p=0.907, n=10
- all eligible samples: all-epithelial TACSTD2 mean_log1p vs T/NK fraction: ρ=-0.115, p=0.751, n=10
- post-treatment eligible (paper MPR/NMPR): malig-like TACSTD2 mean_log1p vs T/NK fraction: ρ=0.179, p=0.702, n=7
- post-treatment eligible (paper MPR/NMPR): malig-like TACSTD2 %pos vs T/NK fraction: ρ=0.179, p=0.702, n=7
- post-treatment eligible (paper MPR/NMPR): malig-like TACSTD2 mean log1p(CP10k) vs T/NK fraction: ρ=0.179, p=0.702, n=7
- post-treatment eligible (paper MPR/NMPR): all-epithelial TACSTD2 mean_log1p vs T/NK fraction: ρ=0.179, p=0.702, n=7
- post-treatment all (no cell-count filter): malig-like TACSTD2 mean_log1p vs T/NK fraction: ρ=-0.164, p=0.651, n=10
- post-treatment all (no cell-count filter): malig-like TACSTD2 %pos vs T/NK fraction: ρ=0.024, p=0.947, n=10
- post-treatment all (no cell-count filter): malig-like TACSTD2 mean log1p(CP10k) vs T/NK fraction: ρ=-0.297, p=0.405, n=10
- post-treatment all (no cell-count filter): all-epithelial TACSTD2 mean_log1p vs T/NK fraction: ρ=-0.049, p=0.88, n=12
- all samples with ≥10 epithelial cells: malig-like TACSTD2 mean_log1p vs T/NK fraction: ρ=-0.231, p=0.448, n=13
- all samples with ≥10 epithelial cells: malig-like TACSTD2 %pos vs T/NK fraction: ρ=-0.069, p=0.823, n=13
- all samples with ≥10 epithelial cells: malig-like TACSTD2 mean log1p(CP10k) vs T/NK fraction: ρ=-0.231, p=0.448, n=13
- all samples with ≥10 epithelial cells: all-epithelial TACSTD2 mean_log1p vs T/NK fraction: ρ=-0.129, p=0.648, n=15
- post-treatment with ≥10 epithelial cells: malig-like TACSTD2 mean_log1p vs T/NK fraction: ρ=-0.164, p=0.651, n=10
- post-treatment with ≥10 epithelial cells: malig-like TACSTD2 %pos vs T/NK fraction: ρ=0.024, p=0.947, n=10
- post-treatment with ≥10 epithelial cells: malig-like TACSTD2 mean log1p(CP10k) vs T/NK fraction: ρ=-0.297, p=0.405, n=10
- post-treatment with ≥10 epithelial cells: all-epithelial TACSTD2 mean_log1p vs T/NK fraction: ρ=-0.049, p=0.88, n=12
- all 15 samples (no cell-count filter): malig-like TACSTD2 mean_log1p vs T/NK fraction: ρ=-0.231, p=0.448, n=13
- all 15 samples (no cell-count filter): malig-like TACSTD2 %pos vs T/NK fraction: ρ=-0.069, p=0.823, n=13
- all 15 samples (no cell-count filter): malig-like TACSTD2 mean log1p(CP10k) vs T/NK fraction: ρ=-0.231, p=0.448, n=13
- all 15 samples (no cell-count filter): all-epithelial TACSTD2 mean_log1p vs T/NK fraction: ρ=-0.129, p=0.648, n=15
- all samples: malig-like TACSTD2 mean_log1p vs residual tumor fraction: ρ=-0.042, p=0.896, n=12
- all samples: all-epithelial TACSTD2 mean_log1p vs residual tumor fraction: ρ=0.046, p=0.875, n=14

- post: NMPR vs MPR malig-like TACSTD2 mean_log1p (sample-level): median 1.624 vs 1.561; MWU p=0.889; n=8 vs 2
- post: NMPR vs MPR malig-like TACSTD2 %pos (sample-level): median 76.716 vs 85.657; MWU p=0.6; n=8 vs 2
- post: NMPR vs MPR malig-like TACSTD2 mean log1p(CP10k) (sample-level): median 1.865 vs 1.466; MWU p=1; n=8 vs 2
- post: NMPR vs MPR all-epithelial TACSTD2 mean_log1p (sample-level): median 1.637 vs 1.575; MWU p=0.683; n=8 vs 4
- post: NMPR vs MPR T/NK fraction (sanity vs paper Fig 1D): median 0.343 vs 0.507; MWU p=0.57; n=8 vs 4
- post: NMPR vs MPR malig-like CX3CL1 mean_log1p (paper: MPR higher): median 0.163 vs 0.299; MWU p=1; n=8 vs 2
- post: NMPR vs MPR malig-like HLA-DRA mean_log1p (paper: MPR higher): median 0.960 vs 2.990; MWU p=0.4; n=8 vs 2
- post: NMPR vs MPR malig-like AKR1C1 mean_log1p (paper: NMPR higher): median 0.717 vs 0.837; MWU p=1; n=8 vs 2
- paired samples: malig-like TACSTD2 > T/NK TACSTD2 (Wilcoxon signed-rank): n=10, p=0.000977; median %pos malig-like=74.9 vs T/NK=1.6.

## Honest limits

- The correlation unit is **samples (n=7–15)**, not ~90k cells. ~90k is the matrix size.
- n=12 post-treatment (4 MPR / 8 NMPR) is underpowered for a ρ≈−0.45 claim.
- MPR malignant compartment is nearly empty after normal-lung exclusion (P06=1, P11=0, P14=0 cells).
- Sample-level tests are the correct unit. Cell-level p-values are pseudoreplication (`rejected_cell_level.json`).
- Without public CopyKAT calls, “malignant-only” is a marker proxy.
- Expression is log1p(raw UMI) and log1p(CP10k), not author Seurat log-normalized values.

## Files

| File | Role |
|---|---|
| `feasibility.json` / `file_manifest.tsv` | Size budget |
| `label_search.json` | Author-label hunt |
| `sample_metadata.tsv` | GEO + paper TN/MPR/NMPR groups |
| `cell_type_composition.tsv` / `sample_composition.tsv` | Marker-based lineages |
| `per_sample_tacstd2.tsv` | Sample-level malignant / T/NK TACSTD2 |
| `association_statistics.tsv` | Spearman vs T/NK |
| `group_statistics.tsv` | NMPR vs MPR MWU |
| `summary.json` / `audit.json` / `sanity_checks.json` | Verdict |
| `rejected_cell_level.json` | Why a 92k-cell ρ is not the test |
| `fig_malig_tacstd2_vs_tnk_fraction.png` | A3-style correlation |
| `fig_tacstd2_nmpr_mpr_and_compartment.png` | Group + restriction |

## Reproduce

```bash
python3 scripts/w200/A3_GSE207422_malig/download.py
python3 scripts/w200/A3_GSE207422_malig/extract.py
python3 scripts/w200/A3_GSE207422_malig/analyze.py
```
