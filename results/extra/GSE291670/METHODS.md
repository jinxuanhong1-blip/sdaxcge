# Methods — EXTRA GSE291670 (not GSE207422)

## Why this series

Requested alternatives:

| Accession | Why not the extra malignant×MPR figure |
|---|---|
| GSE253013 | No matching live GEO record used here. Closest is **GSE243013** (234 NSCLC, post neoadjuvant chemo-IO): CD45+ immune atlas, **0 epithelial/malignant cells**, 6.6 GB MTX. |
| GSE266035 | Blood T cells from **one** NSCLC patient on PD-1/CTLA-4. No tumor epithelium. |
| GSE131907 tLung | Treatment-naive LUAD. No ICI / MPR. Already leftover-scored (tLung n=11, epithelial TACSTD2 vs T/NK ρ=+0.09). |

**GSE291670** (Xia et al., *J Transl Med* 2025; GEO 2025-03-31): 6 NSCLC surgical tumors after neoadjuvant **anlotinib + camrelizumab**. Sample titles are MPR-1/2/3 and Non-MPR-1/2/3. Open 10x MTX (122 MB RAW.tar). Contains epithelium and immune cells.

GSE207422 CopyKAT slide numbers are taken as given. This folder is an extra independent series.

## Data

- `GSE291670_RAW.tar` — six 10x MTX/features/barcodes.
- No author barcode-level cell types on GEO.
- Treatment-naive HRA001033 biopsies used in the paper are **not** in this GEO series and were not mixed in.

## Scoring

1. QC: ≥200 genes, <8000 genes, mitochondrial UMIs <20%.
2. Lineage = argmax of mean `log1p(CP10k)` marker scores (epithelial / T / NK / myeloid / fibroblast / endothelial / B / plasma / mast).
3. **Malignant** = epithelial AND below the sample’s epithelial 75th percentile of a normal-lung marker score (SFTPA/AGER/SCGB/TPPP3/FOXJ1). All-epithelial is a sensitivity.
4. T/NK = assigned T or NK. Sensitivity: CD3E≥1 or CD8A≥1 or NKG7≥1.
5. Per-patient malignant TACSTD2 and CLDN4: mean `log1p(CP10k)` and % UMI≥1.

## Tests (patient is the unit)

- NMPR vs MPR: exact two-sided Wilcoxon by enumerating C(6,3)=20 assignments. At 3 vs 3 the smallest two-sided exact p is **0.10**.
- vs T/NK: Spearman, n=6.

Cell-level p-values are not used.
