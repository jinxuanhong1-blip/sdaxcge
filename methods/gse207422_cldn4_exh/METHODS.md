# Methods — GSE207422 malignant CLDN4 vs T/NK cytotoxicity / exhaustion

Additive extra. Public GEO objects only. Patient is the unit. Not a TACSTD2 redo.

## Data

Hu et al., *Genome Medicine* 2023 (PMID 36869384). GEO [GSE207422](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207422).

- `GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz` (92,330 cells; ~176 MB gzip, under the 2 GB processed-file budget)
- `GSE207422_NSCLC_scRNAseq_metadata.xlsx`

Raw GSA-Human HRA001033 and author CopyKAT barcodes are not public and were not used.

Timepoint from `Resource` (Pre-treatment biopsy vs Post-treatment surgery). Pathologic response from `Pathologic Response` with pCR grouped as MPR. P01 is pathologic NE.

## Expression and scores

Per cell: `log1p(CP10k) = log1p(UMI / library_size × 10⁴)`. QC: library size ≥ 200 UMI.

| Score | Compartment | Genes |
|---|---|---|
| CLDN4 | malignant-like / zero-normal / all-epithelial | CLDN4 |
| Cytotoxicity | T/NK (or CD8) | GZMB, PRF1, GNLY, NKG7 |
| Exhaustion | T/NK (or CD8) | PDCD1, HAVCR2, LAG3, TIGIT, TOX |

Patient metric = mean of the cell-level gene (or mean of the genes in the set) inside the named compartment.

## Lineage and malignant proxies

Lineage = argmax of mean log1p marker modules:

- Epithelial: EPCAM, KRT8/18/19/5/7/17, ELF3, CDH1, MUC1
- T/NK: CD3D/E/G, TRAC, CD2, NKG7, GNLY, KLRD1
- plus B/plasma, myeloid, mast, endothelial, fibroblast

Normal-lung module: SFTPA1/A2/B/D, AGER, NAPSA, SCGB1A1, SCGB3A2, TPPP3, FOXJ1, CAPS.  
Tumor-epi module: EPCAM, KRT8/18/19/7, CEACAM5/6, MUC1, ELF3.

| Call | Rule |
|---|---|
| malignant-like (primary) | epithelial **and** tumor-epi > normal-lung **and** normal-lung < 0.4 |
| zero-normal-UMI | epithelial **and** zero UMI of SFTPA2, AGER, SCGB1A1, SCGB3A1, TPPP3 |
| all-epithelial | epithelial lineage, no normal-lung gate |
| CD8 | T/NK **and** CD8A or CD8B UMI > 0 |

## Patient filter and statistics

Keep a patient when ≥10 cells in the CLDN4 compartment and ≥10 T/NK (CD8 sensitivity: ≥10 CD8).  
Primary Spearman is **post-treatment only**. Pre biopsies are unpaired and are not stacked into that ρ.  
Tests are two-sided Spearman (or Mann–Whitney for MPR vs NMPR). Cells are not treated as n.

## How to run

```bash
python3 methods/gse207422_cldn4_exh/download.py
python3 methods/gse207422_cldn4_exh/extract.py
python3 methods/gse207422_cldn4_exh/analyze.py
```
