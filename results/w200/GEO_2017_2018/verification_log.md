# Verification log — GEO 2017–2018 leftover

Date: 2026-08-16. All accessions from NCBI esearch/esummary. No invented IDs.

## Search
- Script: `scripts/w200/geo_2017_2018/01_search_geo.py`
- Combined lung × ICI query, PDAT 2017/01/01–2018/12/31, Homo sapiens, GSE: count=22
- Per-drug add-ons (nivolumab, pembrolizumab, atezolizumab, durvalumab,
  ipilimumab, avelumab, cemiplimab): 0 new UIDs
- Output: `search_uids.json` (22 UIDs), `candidates_metadata.json` (22 records)

## Deposited-file checks (this slice)
| Accession | File | What was verified |
|---|---|---|
| GSE93157 | `downloads/GSE93157/GSE93157_raw_data_values.txt.gz` | 770 gene symbols; TACSTD2/TROP2/CLDN4/CLAUDIN* = 0 hits |
| GSE93157 | `downloads/GSE93157/GSE93157_series_matrix.txt.gz` | 65 GSM; source counts MELANOMA 25 / LUNG NON-SQUAMOUS 22 / SQUAMOUS LUNG 13 / HEADNECK 5; `best.resp` + `pfs` present; GEO `response` = RC_RP_SD for all 35 lung including 14 PD |
| GSE110390 | `downloads/GSE110390/GSE110390_CP1108_NSCLC_21gene_expression.txt.gz` | 21 genes × 97 NSCLC samples; TACSTD2/CLDN4 not in list |
| GSE110390 | `downloads/GSE110390/GSE110390_series_matrix.txt.gz` | characteristics = tissue + needle biopsy only; `data_row_count=0`; SRA linked; no outcome fields |
| GSE91061 | series matrix | `tissue: melanoma` on all samples; response present; not lung |
| GSE99531, GSE100860, GSE90728, GSE124199 | series matrices | T-cell / DC designs; empty expression tables; no ICI endpoint on GSM |

GSE111360 and GSE99254 have no series-matrix on the standard FTP path
(404). Brief SOFT + suppl listings used instead (organoid scRNA; T-cell scRNA).

## What was not done
- No SRA / EGA reprocessing
- No association statistics
- No analysis of GSE91061 as a lung proxy
