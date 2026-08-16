# Verification log — GEO 2019 leftover

Date: 2026-08-16. All accessions below were returned by NCBI E-utilities or
read from NCBI SOFT / series matrix. None were invented.

## Search

- Script: `scripts/w200_geo_2019/01_search_geo.py`
- Database: `gds`, entry type GSE, organism Homo sapiens
- Date window: `("2019/01/01"[PDAT] : "2019/12/31"[PDAT])`
- Combined lung × ICI count = 24; union with response-wording query = 26
- Target-gene query (TACSTD2 / TROP2 / CLDN4 / sacituzumab / datopotamab / SKB264) = 0
- Named-drug queries (nivolumab, pembrolizumab, atezolizumab, durvalumab,
  ipilimumab, avelumab, cemiplimab, tremelimumab, camrelizumab, sintilimab,
  tislelizumab, toripalimab) each returned 0 in this 2019 × lung window
- Output: `results/w200/GEO_2019/search/search_uids.json`

Sanity checks (same Entrez API):

- `GSE135222[Accession]` exists and has PDAT 2019/08/02
- `GSE126044` PDAT 2020/02/03 (out of this slice; already in fable PR)
- `GSE111414` PDAT 2021/03/01 (out of this slice; already in fable PR)
- `GSE136961` PDAT 2020/02/03 (out of this slice; fable: no TACSTD2/CLDN4 on panel)
- `GSE182328` PDAT 2021/08/20 (out of this slice; already in fable PR)

## Metadata

- `02_fetch_metadata.py` pulled esummary + SOFT brief headers for all 26 UIDs
- Submission dates recovered from `!Series_submission_date` (not an Entrez field)

## Sample-level review

- `03_triage_and_soft.py` downloaded series matrices where NCBI hosts them
- Series with empty expression tables (header only): GSE117570, GSE118933,
  GSE120028, GSE127825, GSE133605, GSE135164, GSE142620, GSE139555
  (expression lives in SRA / supplementary files, not the matrix)
- GSM brief used when no series matrix: GSE111898, GSE117049, GSE120101,
  GSE127465, GSE129380, GSE129381, GSE138571
- GSE120028 sources inspected: gastric, GE, CRC (bronchus met), pancreatic
  liver met — no primary lung
- GSE135164 characteristics inspected: one of 15 PDCs annotated
  “Chemotherapy, immunotherapy”; no response field
- GSE139555 phenotypes inspected: 13 lung immune GSMs, no ICI field

## Skip list (fable 2019–2021 PR)

Do not re-analyze: GSE135222, GSE126044, GSE111414, GSE182328, GSE136961.
Only GSE135222 falls in the 2019 PDAT window.

## Outcome

`04_build_triage.py` → 26/26 reviewed, 0 leftover series marked
`analyzable_tacstd2_cldn4_ici = yes`.
