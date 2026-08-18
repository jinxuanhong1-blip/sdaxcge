# GSE207422 — Seurat CreateSeuratObject, Cldn4-only malignant vs T/NK and IFN/MHC

**Additive Seurat slice.** Thesis already correct. CLDN4 is the splitter. TACSTD2 is never a gate. Dual-high was not run. GSE148071 was not merged. Python was used only to stream the public UMI into 10x MTX and to write the GEO xlsx as TSV. Every scored value and every p-value comes from **R / Seurat 5.0.1** after `CreateSeuratObject` + `NormalizeData` (LogNormalize, scale.factor=1e4).

**Verdict (honest n=6 paired post patients, not 12):** in A3-malignant cells, CLDN4-high (within-tumor Q4) is higher than CLDN4-low (Q1) for the **IFN ISG** module in the same cells (6/6, exact Wilcoxon p=0.031). MHC-I APM goes the same way in 5/6 (p=0.063). Patient-mean CLDN4 vs IFN ISG is ρ=+0.94 (p=0.0048, n=6). Patient-mean CLDN4 vs T/NK fraction is **near-null** (n=6 min-20 ρ=+0.37 p=0.47; n=10 all post with any A3-malignant ρ=−0.13 p=0.73; epithelial complete-case n=12 ρ=−0.007 p=0.98). ICI labels are present (MPR / NMPR / pCR). After the occupancy floor, NMPR vs MPR malignant CLDN4 is 5 vs 1 and is **not tested**.

Do **not** cite n=12 for malignant tests. Six of twelve post patients have <20 A3-malignant cells (including 3 of 4 MPR). Cell-level p-values are not reported.

## Data and n

- Public GEO UMI only: **92,330** cells × **24,292** genes (`CreateSeuratObject` on the streamed MTX). Author CopyKAT / epithelium RDS barcodes are not on GEO. No processed Seurat RDS is deposited.
- **Unit of every primary test is the patient.** One post-treatment surgery sample per patient (Hu et al., *Genome Med* 2023, PMID 36869384). The three pre-treatment biopsies (P01 NE, P05/P08 NMPR-labeled) are excluded from tests.
- ICI labels from `GSE207422_NSCLC_scRNAseq_metadata.xlsx`: Pathologic Response. **pCR (P06) is folded into MPR** as in the paper (post MPR n=4, NMPR n=8). RECIST and PD-1 antibody (toripalimab / camrelizumab / sintilimab) are carried in `sample_metadata.tsv` and are not used as the primary split.
- Marker lineages (argmax of Seurat `data` layer marker means; same panels as the public CLDN4 slices): epithelial **10,987**; A3-malignant-like **6,608**; T **29,448**; NK **7,097**; T/NK **36,545**. Counts sit within ~2% of the prior Python UMI slices (those used raw log1p for lineage, not `NormalizeData`).
- A3-malignant-like = epithelial AND zero UMI for SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3 (**not** CopyKAT). CLDN4-only: TACSTD2 is scored as a companion and is **never a gate**.
- Primary floor: ≥20 A3-malignant cells. Paired tails additionally require ≥8 cells in each CLDN4 quartile. **Eligible post n=6:** P03 (MPR), P04, P07, P09, P10, P12 (NMPR).
- Dropped post (A3-malignant n): P02 (10), P06 (1), P11 (0), P13 (3), P14 (0), P15 (4). Three of four MPR are empty or one cell.
- IFN ISG 40/40 genes present. MHC-I APM 21/21 present. A3 normal-lung panel complete.

## Definitions

| Item | Rule |
|---|---|
| Object | `ReadMtx` → `CreateSeuratObject(min.cells=0, min.features=0)` → `NormalizeData` |
| CLDN4 / modules | mean of Seurat `data` layer = `log1p(CP10k)` |
| CLDN4 split | within-patient Q4 vs Q1 of malignant CLDN4 |
| IFN ISG | 40-gene type-I ISG core (ISG15, MX1, OAS*, IFIT*, STAT1/2, IRF7/9, …) |
| IFN core-6 | IFI27, OAS2, IFIT1, MX1, ISG15, HLA-A |
| MHC-I APM | HLA-A/B/C/E/F, B2M, TAP1/2, TAPBP, PSMB8/9/10, … |
| OXPHOS control | NDUFA1, COX5A, ATP5F1A/B, UQCRC1, SDHA, … |
| T/NK | lineage T or NK / all cells, per patient |
| AddModuleScore | Seurat-native sensitivity only (not the primary mean) |
| Primary n | post patients with ≥20 A3-malignant; paired tests also ≥8/tail |
| Dual-high | **not run** |
| GSE148071 | **not merged** |

## Primary table — patient unit

A3-malignant, 12 post patients attempted, **6 eligible** for min-20 / paired tails.

| Test | n | result |
|---|---:|---|
| malignant CLDN4 vs T/NK fraction | **6** | ρ=+0.37, p=0.47 |
| malignant CLDN4 %pos vs T/NK fraction | 6 | ρ=−0.03, p=0.96 |
| malignant CLDN4 vs IFN ISG | **6** | ρ=+0.94, p=0.0048 |
| malignant CLDN4 vs IFN core-6 | 6 | ρ=+0.94, p=0.0048 |
| malignant CLDN4 vs MHC-I APM | 6 | ρ=+0.26, p=0.62 |
| malignant CLDN4 vs OXPHOS | 6 | ρ=−0.71, p=0.11 |
| malignant CLDN4 vs AddModuleScore IFN | 6 | ρ=+0.83, p=0.042 |
| malignant CLDN4 vs AddModuleScore MHC-I | 6 | ρ=+0.26, p=0.62 |
| paired IFN ISG, CLDN4 Q4 vs Q1 | **6** | 0.441 vs 0.367; Δ median +0.053; 6/6; p=**0.031** |
| paired IFN core-6, CLDN4 Q4 vs Q1 | 6 | 1.267 vs 1.142; 6/6; p=**0.031** |
| paired MHC-I APM, CLDN4 Q4 vs Q1 | 6 | 1.412 vs 1.247; 5/6; p=0.063 |
| paired OXPHOS, CLDN4 Q4 vs Q1 | 6 | 0.947 vs 0.865; 4/6; p=0.16 |
| paired library size (nCount_RNA) | 6 | 14542 vs 6182; 6/6; p=0.031 |

Exact two-sided Wilcoxon p cannot go below 0.031 at n=6 when all signs agree. R reports V=21 (sum of positive ranks) for the 6/6 IFN test; that is the same test as a scipy W=0.

P07 IFN ISG delta is +0.0016 (visually flat). The other five IFN deltas are +0.046 to +0.182.

### Sensitivities that keep more patients (still patient unit)

| Test | n | result |
|---|---:|---|
| A3 CLDN4 vs T/NK, all post with any A3-malignant | 10 | ρ=−0.13, p=0.73 |
| epithelial CLDN4 vs T/NK (complete-case) | 12 | ρ=−0.007, p=0.98 |
| epithelial CLDN4 vs IFN ISG | 12 | ρ=+0.54, p=0.071 |
| epithelial CLDN4 vs MHC-I | 12 | ρ=+0.16, p=0.62 |
| epithelial paired IFN ISG Q4 vs Q1 | 11 | 10/11; p=0.014 (P13 tail<8 dropped) |
| epithelial paired MHC-I Q4 vs Q1 | 11 | 7/11; p=0.067 |

The n=10 T/NK ρ=−0.13 matches the prior public CLDN4/PD-L1 slice on this UMI. The min-20 n=6 set is 5 NMPR + P03; it is underpowered and is not a sign flip of the thesis.

## ICI labels (present)

Post-treatment surgery: MPR n=4 (P03, P06 pCR, P11, P14), NMPR n=8. Pre biopsies are not stacked into the contrast.

| Contrast | n | result |
|---|---|---|
| A3-malignant CLDN4, NMPR vs MPR, min-20 | 5 vs 1 | **too few — not tested** |
| A3-malignant CLDN4, NMPR vs MPR, any A3 | 8 vs 2 | 1.49 vs 0.51; p=0.18 |
| epithelial CLDN4, NMPR vs MPR | 8 vs 4 | 1.52 vs 1.50; p=0.93 |
| T/NK fraction, NMPR vs MPR | 8 vs 4 | 0.38 vs 0.52; p=0.68 |

Do not read an ICI response claim from 5 vs 1 or from two MPR samples with usable malignant CLDN4 (P03=734 cells, P06=1 cell). Epithelial NMPR vs MPR CLDN4 is null, matching prior public slices.

## Honest limits

1. **n=6, not 12.** Half the post cohort has almost no A3-malignant cells. MPR residual tumor is the main hole (P06=1, P11=0, P14=0). Header n stays 12 attempted / 6 tested.
2. Exact Wilcoxon p=0.031 is the smallest two-sided p at n=6. One IFN pair (P07) is a near-tie.
3. CLDN4-high cells have more UMI (6/6, p=0.031). Scores are CP10k-normalized. OXPHOS does not significantly rise (4/6, p=0.16), but depth is not fully ruled out.
4. This is not Hu et al. CopyKAT. A3-malignant-like can leak unmarked epithelium.
5. Dual-high (TACSTD2 AND CLDN4) was not run. TACSTD2 is reported only as a companion column on the per-patient table.
6. GSE148071 was not downloaded or merged.
7. Lineage is a marker argmax on the Seurat `data` layer, not author cell types (those barcodes are not public).
8. ICI contrast after the occupancy floor is 5 vs 1. Do not cite NMPR>MPR malignant CLDN4 from this object.

## Figures

- `figures/fig_cldn4_vs_tnk.png` — patient malignant CLDN4 vs T/NK; point size = A3-malignant n; color = ICI
- `figures/fig_cldn4_vs_ifn_mhc.png` — patient malignant CLDN4 vs IFN ISG and MHC-I (min-20)
- `figures/fig_paired_ifn_mhc.png` — same-cell paired Q4 vs Q1 (n=6)
- `figures/fig_honest_n.png` — A3-malignant occupancy; floor n=20 visible

## Files

- `tables/primary_table.tsv` — compact primary tests
- `tables/per_patient.tsv` — 15 samples; tests use post rows
- `tables/paired_high_low.tsv` / `tables/eligibility.tsv` / `tables/tests.tsv`
- `tables/lineage_counts.tsv`
- `sample_metadata.tsv` — GEO clinical + ICI labels
- `summary.json` / `sessionInfo.txt`
- Scripts: `scripts/download.py`, `scripts/00_write_metadata.py`, `scripts/01_umi_to_10x.py`, `scripts/02_analyze.R`

## Reproduce

```bash
python3 methods/seurat_gse207422_cldn4/scripts/download.py
python3 methods/seurat_gse207422_cldn4/scripts/00_write_metadata.py
python3 methods/seurat_gse207422_cldn4/scripts/01_umi_to_10x.py
Rscript methods/seurat_gse207422_cldn4/scripts/02_analyze.R \
  data/GSE207422/tenx \
  methods/seurat_gse207422_cldn4
```

Requires R ≥ 4.3 with Seurat ≥ 5.0 (`CreateSeuratObject`, `NormalizeData`, `AddModuleScore`). If Seurat or the public UMI matrix is missing, stop — do not substitute a Python-only primary.
