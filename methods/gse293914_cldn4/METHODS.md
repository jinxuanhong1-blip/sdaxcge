# Methods — ADDITIVE GSE293914 (H1975 xenograft)

## Why this series

Public EGFR-mutant **H1975 xenograft** snRNA after anti-PD-1 ± anti-CCL20 in a humanized NSG + PBMC model (Kwok et al., PMID 42549034). Requested as an additive CLDN4/TACSTD2 score in tumor vs T/NK, with honest n and an explicit xenograft label.

## Public files only

| File | Used |
|---|---|
| `GSE293914_series_matrix.txt.gz` (3.3 KB, 0 expression rows) | Sample title, characteristics, protocol |
| `GSM8893263_{barcodes,features,matrix}.mtx/tsv.gz` | Cell Ranger 7.1 raw Flex matrix (37,143 × 2,032,489; 111,813,751 nnz) |
| SRA / FASTQ / author RDS | **Not used** |

GEO lists **one** sample. Title: “H1975 xenograft, isotype, anti-pd-1, anti-ccl20, combination”. Characteristics: `treatment: isotype`. Extract protocol: 10x Chromium Fixed RNA Profiling multiplex (CG000565), 20,000 nuclei targeted per sample, GRCh38-2020-A, Cell Ranger 7.1.0. Author QC stated on GEO: genes <200 or >8,000 or MT% >10 removed.

## Probe barcodes

Cell barcodes are 24-nt Flex (16-nt GEM + 8-nt translated probe). Four IDs dominate the raw matrix and all QC cells:

| 8-mer | ID | Raw barcodes | QC cells |
|---|---|---:|---:|
| ACTTTAGG | BC001 | 263,379 | 4,154 |
| AACGGGAA | BC002 | 457,885 | 8,503 |
| AGTAGGCT | BC003 | 450,075 | 13,037 |
| ATGTTGAC | BC004 | 516,301 | 17,217 |

IDs follow the 10x Flex v1 translation table. The other 12 probe 8-mers are unused-channel raw barcodes and contribute **0** QC cells. GEO does not say which BC is isotype / aPD-1 / aCCL20 / combo.

## Scoring

1. Stream the MTX (do not load the full sparse matrix).
2. QC: ≥200 genes, <8,000 genes, mitochondrial UMIs <10% (author).
3. Lineage = argmax of mean `log1p(CP10k)` marker modules (epithelial / T / NK / B / plasma / myeloid / mast / fibroblast / endothelial). Same marker lists as the GSE291670 leftover scorer.
4. **Tumor** = epithelial. This is a cell-line xenograft; no normal-lung filter.
5. **T/NK** = assigned T or NK. Sensitivity: CD3E≥1 or CD8A≥1 or NKG7≥1; also PTPRC.
6. Per-cell TACSTD2 and CLDN4: `log1p(CP10k)` and UMI≥1.

## Tests

None. Honest n = 1 library. Probe groups are unlabeled multiplex channels, not biological replicates. Cell-level p-values are not reported.

## Reproduce

```bash
python3 methods/gse293914_cldn4/download.py --datadir data/GSE293914
python3 methods/gse293914_cldn4/analyze.py \
  --datadir data/GSE293914 \
  --outdir methods/gse293914_cldn4
```
