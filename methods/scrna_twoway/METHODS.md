# Methods — two-way lung scRNA extra

## Inclusion

A public series is two-way if it has (i) a timepoint axis (pre / post / treatment-naive vs treated), (ii) a response axis (MPR / NMPR / pCR or RECIST), and (iii) a processed count matrix that includes epithelial or malignant cells. Restricted archives (EGA, dbGaP, GSA HRA) without a public processed matrix are excluded. T-cell-only or CD45-sorted deposits are recorded but not used for malignant TACSTD2/CLDN4.

Hunt: NCBI eUtils GDS searches for lung/NSCLC + scRNA + immunotherapy/MPR/RECIST, plus a seed list of known ICI lung series. SOFT/series-matrix traits were read for the seed and for text-level two-axis hits.

## GSE207422

Hu et al., *Genome Medicine* 2023 (PMID 36869384). Public files:

- `GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz` (92,330 cells)
- `GSE207422_NSCLC_scRNAseq_metadata.xlsx`

Timepoint from `Resource` (Pre-treatment biopsy vs Post-treatment surgery). Pathologic response from `Pathologic Response` with pCR grouped as MPR. RECIST from `RECIST`. P01 is pathologic NE and is omitted from the MPR 2×2.

Lineages are assigned from the paper’s canonical markers (T, NK, B, plasma, myeloid, neutrophil, mast, epithelial, fibroblast, endothelial). Normal lung programs: alveolar (SFTPA1/2, SFTPC/B, AGER), club (SCGB1A1, SCGB3A2), ciliated (TPPP3, FOXJ1, CAPS). Malignant = epithelial and not a clear normal-lung program, with a chromosome-mean |z| dispersion score using fibroblasts/endothelia as the reference (author CopyKAT barcodes are not deposited). A broader epithelial-minus-normal call is also stored.

Expression: per-cell log1p(CP10k). Patient metric = mean over malignant cells. Patients with <10 malignant cells are dropped from tests (none dropped under the primary call). Tests are two-sided Mann–Whitney on patients. Cells are not treated as n.

## Stacking

If a second public two-way series existed, patient-level malignant means would be stacked with a series indicator. None was found.
