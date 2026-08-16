# GEO 2019 leftover: TACSTD2 / CLDN4 vs lung ICI

Honest result: **there is no leftover 2019 GEO cohort to analyze.**

After an independent 2019-only search and sample-level metadata review, the
only human lung ICI series with processed tumor expression and a deposited
per-sample ICI endpoint is **GSE135222**. That series was already tested for
TACSTD2 and CLDN4 in `cursor/fable-geo-2019-2021-c71e` (n=27, PFS; neither
gene associated). It is not repeated here.

## What was searched

NCBI `gds`, GSE, *Homo sapiens*, `PDAT` 2019-01-01 to 2019-12-31.

Queries (union): lung × ICI/drug terms (including camrelizumab / sintilimab /
tislelizumab / toripalimab); lung × TACSTD2/TROP2/CLDN4/ADC names; lung ×
responder / DCB / PFS wording.

| Query | Hits |
|---|---|
| lung × ICI | 24 |
| lung × TACSTD2/CLDN4/ADC | **0** |
| lung × response wording | 4 (2 new) |
| per-drug lung queries | 0 additional |
| **unique series** | **26** |

All 26 accessions are real NCBI GSE IDs. Submission dates from SOFT are
reported next to PDAT in `results/w200/GEO_2019/search/candidates_metadata.json`.

The four fable-analyzed 2019–2021 lung ICI series that are **not** in this
2019 PDAT slice (and were therefore out of scope, not “missed”) are
GSE126044 (2020), GSE136961 (2020; panel lacks both genes), GSE111414 (2021;
CD8 PBMC), and GSE182328 (2021).

## Closest leftovers (why they still fail)

| Series | Why it is not a TACSTD2/CLDN4 × lung ICI test |
|---|---|
| GSE135222 | Already analyzed in the fable 2019–2021 PR. Skip. |
| GSE119144 | Methylation companion of the same paper. GEO labels are only “tumor”. No ICI outcome to correlate with TACSTD2/CLDN4 CpGs. |
| GSE135164 | 15 NSCLC PDCs. One line (PDC6) says “Chemotherapy, immunotherapy”. No drug, no response, n=1. |
| GSE139555 | scRNA/TCR of CD3+/CD45+ cells; 13 GSMs are lung. No ICI field. Epithelial genes are the wrong compartment. |
| GSE120028 / GSE120101 | IT1208 anti-CD4 phase 1. No primary lung tumor. Not PD-1/PD-L1/CTLA-4. |
| GSE117570 / GSE127465 / GSE111894 | Lung immune scRNA, treatment-naive or no ICI. |
| GSE133605 | n=1 ALK+/PD-L1+ biopsy reanalysis. PD-L1 is a marker, not treatment. |
| GSE138571 | Mouse KP ATAC ± Asf1a; anti-PD-1 is a mouse experiment, not a human expression cohort. |
| Pirfenidone / MET / IL-1β / polycomb series | Cell lines. “Checkpoint” appears in the background text only. |
| GSE117049 | LUAD hypoxia lncRNA. Matched because “ICI” occurs as a substring in unrelated words. Not immunotherapy. |
| GSE71799 / GSE124574 / GSE127825 | Cystic fibrosis, melanoma, thymus. |

Full row-level reasons: `results/w200/GEO_2019/tables/triage_leftover_2019.csv`.

## What is not claimed

- No new TACSTD2 or CLDN4 association with ICI response, PFS, or OS.
- No claim that 2019 GEO “contains no lung ICI data” — GSE135222 exists and
  was already used.
- No methylation-vs-outcome analysis of GSE119144, because the outcome is not
  in that series.
- No n=1 PDC6 contrast dressed up as a cohort result.

If a later slice (2020+) is meant to pick up GSE126044 / GSE136961, that is
outside this 2019 leftover task.
