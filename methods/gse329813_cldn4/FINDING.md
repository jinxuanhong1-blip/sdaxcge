# GSE329813 leftover public matrix — CLDN4 vs T/NK or IFN

**Additive only. Verdict: EMPTY.**

Leftover from the durvalumab / PACIFIC public-matrix hunt. GEO [GSE329813](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE329813) is **not durvalumab**. It is post-neoadjuvant **pembrolizumab + platinum** GeoMx DSP (NCT05383716; single-arm phase II). TACSTD2 vs MPR / T/NK on this same matrix is already scored elsewhere and is treated as given. This slice is **CLDN4 vs T/NK or IFN**.

**CLDN4 is not on the public panel.** Pairwise n is **0**. No Spearman is computed. TACSTD2 and EPCAM are not used as stand-ins.

## Honest n

| Item | Public? | n | Note |
|---|---|---:|---|
| GSM / ROI titles | yes | **127** | GSM9711715–GSM9711841 |
| Titles parsed (ROI, patient, site, MPR/NMPR) | yes | **127** | all titles match |
| Patients | yes | **22** | P1–P22; design text agrees |
| Primary tumor-bed ROIs | yes | **72** | 20 patients × 3 ROI, 2 patients × 6 ROI |
| Lymph-node ROIs | yes | **55** | 18/22 patients have LN ROIs |
| Tumor-bed patients (MPR / NMPR) | yes | **22** | **11 / 11** |
| Processed matrix genes × ROIs | yes | **1812 × 127** | `GSE329813_processed_data_file_normalized_data.csv.gz` |
| Series-matrix expression rows | empty | **0** | header only; `!Sample_data_row_count` = 0 |
| Durvalumab / PACIFIC RNA | no | **0** | pembrolizumab + platinum leftover |
| CLDN4 finite values | no | **0** | no CLDN4; no claudin-family gene |
| **Primary pairwise n (CLDN4 + T/NK)** | empty | **0** | T/NK genes are on the panel; CLDN4 is not |
| **Primary pairwise n (CLDN4 + IFN)** | empty | **0** | IFN genes are on the panel; CLDN4 is not |

Do not write n=22 or n=127 for a CLDN4 correlation. Those are ROI / patient counts on a matrix that does not contain CLDN4.

## One-row table

| dataset | drug | assay | n ROI | n patients | CLDN4 | T/NK genes | IFN genes | CLDN4–T/NK n | CLDN4–IFN n | verdict |
|---|---|---|---:|---:|---|---|---|---:|---:|---|
| GSE329813 | pembrolizumab + platinum (NCT05383716); **not Durva** | GeoMx DSP, total-count norm (scale 10,210) | 127 | 22 | **ABSENT** | 13/13 | 12/12 | **0** | **0** | **EMPTY** |

Full row: `tables/one_row.tsv`.

## Panel check

Only public counts: the supplementary normalized CSV. The series matrix has no expression table.

| Query | On processed matrix | Aliases searched |
|---|---|---|
| **CLDN4** | **no** | CLDN4, Claudin4, Claudin-4, CLDN-4, CPE-R, CPER, CPETR, CPETR1, WBSCR8, NM_001305, ENSG00000189143 |
| Any CLDN / claudin family | **no** | substring `CLDN` or `CLAUDIN` on all 1812 IDs |
| TACSTD2 | yes | not a CLDN4 substitute |
| EPCAM | yes | not a CLDN4 substitute |

T/NK genes present (13/13): CD3D, CD3E, CD3G, CD2, CD8A, CD8B, NKG7, GNLY, KLRD1, GZMB, PRF1, GZMA, GZMK.

IFN genes present (12/12): IFNG, STAT1, CXCL9, CXCL10, CXCL11, IDO1, HLA-DRA, HLA-A, HLA-B, HLA-C, CD274, LAG3.

The immune axes are computable. The target is not. Nearby epithelial genes were not scored as CLDN4.

## Assigned tests (not run)

| pair | n | ρ | p | verdict |
|---|---:|---|---|---|
| CLDN4 vs T/NK mean-z | **0** | — | — | **EMPTY** |
| CLDN4 vs IFN / Ayers-like mean-z | **0** | — | — | **EMPTY** |
| CLDN4 vs CD8A | **0** | — | — | **EMPTY** |
| CLDN4 vs IFNG | **0** | — | — | **EMPTY** |

No Q4/Q1, no MPR contrast, no purity residual. Those require a CLDN4 column.

## What this is / is not

- **Leftover Durva hunt:** open post-treatment GeoMx after pembro + chemo. Not PACIFIC, not GSE248378, not GSE253564.
- **Not a CLDN4–immune law.** The public matrix cannot test it.
- **Not a re-score of TACSTD2.** Tumor-bed TACSTD2 vs T/NK (n=22) is already on this accession and is not this claim.
- **Not bulk RNA-seq or scRNA.** ROI-level GeoMx CTA-scale panel (~1.8k genes).

## Reproduce

```bash
python3 -m pip install -r methods/gse329813_cldn4/requirements.txt
python3 methods/gse329813_cldn4/analyze.py
```

Downloads (not committed) go to `$GSE329813_CLDN4_DATA` (default `/tmp/gse329813_cldn4`):

- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE329nnn/GSE329813/suppl/GSE329813_processed_data_file_normalized_data.csv.gz`
- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE329nnn/GSE329813/matrix/GSE329813_series_matrix.txt.gz`

## Files

- `analyze.py` — download, panel check, honest-n tables; stops without scoring CLDN4
- `tables/label_inventory.tsv` — public vs empty fields
- `tables/gene_coverage.tsv` — CLDN4 / T/NK / IFN / nearby presence
- `tables/cldn4_vs_tnk_ifn.tsv` — empty assigned tests
- `tables/one_row.tsv`
- `tables/roi_annotation.tsv` — 127 parsed titles
- `tables/summary.json`
