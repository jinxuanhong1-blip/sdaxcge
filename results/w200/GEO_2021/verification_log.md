# Verification log — GEO 2021 leftover TACSTD2/CLDN4

All accessions are real NCBI GSE IDs from E-utilities. None were invented.

## Search
- Script: `search_2021.py`
- Database: NCBI `gds`
- Filters: `"gse"[Entry Type]`, `"Homo sapiens"[Organism]`, PDAT 2021/01/01–2021/12/31
- Combined lung × ICI query: count=38
- Per-drug supplements (nivolumab, pembrolizumab, atezolizumab, durvalumab, ipilimumab, avelumab, cemiplimab): 0 new UIDs
- Outputs: `search_uids_2021.json`, `candidates_2021.json`

## Already done (not leftover)
- GSE111414: PBMC CD8, prior slice
- GSE182328: Akkermansia surrogate, prior slice

## Analyzed leftover
- GSE190265: `GSE190265_TPM_France3.csv.gz` + `GSE190265_samples_info_France3.csv.gz` (43/43 overlap). Series matrix has 34 GSM; 26 titles match TPM IDs. Histology available for 25 overlapping labeled samples.
- GSE190266: `GSE190266_TPM_France4.csv.gz` (70×16383, last gene MTMR14) + series-matrix 6-month PFS. TACSTD2/OPTN/TLR9 absent.

## Side audit
- GSE146100: `GSE146100_NormData.txt.gz` SHA-256 `b937b797bd701da06e1b616f72364ee9fdaecb6d3cf74be2b7ae667614c378b6`

## File checksums
See `input_manifest.csv`.
