# Hunt log: GSE207422 cell-level annotation

Date: 2026-08-16. Goal: find **any** public per-cell annotation (author CopyKAT preferred). Full repo clones were not used.

## Author / GEO (no cell IDs)

| Source | What is there | Cell-level malignant IDs? |
|---|---|---|
| GEO GSE207422 supplementary | UMI matrix 92,330 cells; `GSE207422_NSCLC_scRNAseq_metadata.xlsx` (15 samples × clinical columns) | **No** |
| GEO GSM records GSM6287372–GSM6287410 | Sample phenotype only | **No** |
| GEO SOFT | Same 4 supplementary files. Raw FASTQ pointed to GSA **HRA001033** (controlled access) | **No** |
| Paper Additional file 1 (Table S1) | Clinical tables (scRNA / bulk / IHC / metabolomics) | **No** |
| Paper Additional file 2 | Supplementary figures PDF (~20 MB), including Fig. S3 CopyKAT UMAPs | Figures only, **no barcode table** |
| Paper Additional file 3 (Table S2) | Gene lists for modules (MHC-II, cytotoxicity, etc.) | **No** |
| Paper Additional file 4 (Table S3) | Steroid intensity matrix | **No** |
| Author GitHub `Junjie-Hu/NSCLC-immunotherapy` | R/Python scripts only. `03.1_Tumor_copykat_analysis.R` reads local `epithelium.rds` / `stromal.rds`, writes `data_out/copykat_res.rds` | **No deposited annotation file** |
| Author GitHub issues | #3 asks for Ensembl IDs; reply is `gencode.v22.annotation.gtf` | **No** |
| GSA HRA001033 | Controlled-access FASTQ via DAC HDAC000619 | Not a public annotation |

## Portals / reanalyses

| Source | Result |
|---|---|
| TISCH2 NSCLC gallery | 17 datasets; **GSE207422 not included** |
| CELLxGENE curation API | No GSE207422 / PMID 36869384 |
| Figshare search `GSE207422` | Unrelated NSCLC bulk tables; no Hu 2023 barcodes |
| Zenodo search `GSE207422` | No cell-annotation deposit. LOSTdb Part 6 has `Hu_2023.tar.gz` (12.5 GB expression+metadata bundle). LOSTdb GitHub `Other_data/Metadata in LOSTdb.csv` is **sample/dataset catalog**, not 92k barcodes. Not downloaded (full clone / huge tarball). |
| `yasmina-bioinfo/scRNA_LUAD_Immunotherapy` | Scripts; states matrices/annotations are not in the repo |
| J Gene Med 2024 reanalysis (doi 10.1002/jgm.3736) | Describes SingleR annotation of GSE207422; **no public barcode file** found |
| CancerSCEM / SCAR / Single Cell Portal / UCSC Cell Browser | No GSE207422 annotation file located |

## Found (used)

**DRMref** (https://ccsm.uth.edu/DRMref/ ; Liu et al. NAR 2024):

- Download index: `gene_search_result_0.cgi?type=all_dataset_info_download`
- `https://ccsm.uth.edu/DRMref/All_RData_after_Annotation/GSE207422_Tor_seurat_afterAnno.RDS` (306,232,023 bytes)
- `https://ccsm.uth.edu/DRMref/All_RData_after_Annotation/GSE207422_Sin_seurat_afterAnno.RDS` (799,394,628 bytes)

`meta.data$celltype` includes `Malignant cells`. This is DRMref’s own annotation, **not** author CopyKAT. No Camrelizumab/TN objects exist on that page.

TACSTD2 counts were taken from GEO `GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz` (all 30,877 annotated barcodes match).
