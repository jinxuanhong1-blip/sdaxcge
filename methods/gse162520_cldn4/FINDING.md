# GSE162520 TUMADOR — CLDN4-only vs CD8A / Immune / CD274 / IFN

**Verdict: panel missing. Stop.**

**CLDN4 is not on the HTG EdgeSeq Oncology Biomarker Panel (OBP)** used in GSE162520. No CLDN4 vs CD8A / Immune / CD274 / IFN scores. No LUAD vs LUSC split on CLDN4. No dual-high.

Foy et al., *Eur J Cancer* 2022 ([PMID 36038492](https://pubmed.ncbi.nlm.nih.gov/36038492/); GEO [GSE162520](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162520)). Public file `GSE162520_GEO_data_TUMADOR_log2cpm.csv.gz`: Centre Léon Bérard TUMADOR surgical NSCLC, HTG EdgeSeq OBP, log2 CPM. Series summary: 92 early-stage NSCLC treated with surgery; FFPE biopsy or resection collected before treatment.

GEO title says “advanced” NSCLC and PD-1/PD-L1 response. That is the paper’s ICI question. The deposited matrix is the **surgical TUMADOR** cohort used to test prognostic impact of the 27-gene HOT score, not an ICI-response series.

## Honest n

| Item | n | Source |
|---|---:|---|
| GEO HTG arrays / unique `CLB_*` columns | **92** | `GSE162520_GEO_data_TUMADOR_log2cpm.csv.gz` (CLB_1–CLB_96 with 27, 28, 35, 62 unused) |
| GEO GSM records | **92** | GSM4953658–GSM4953749 |
| Series summary / overall design | **92** | “92 patients with early stage NSCLC treated with surgery” |
| Histology labeled | ADK **47** / SCC **45** | SOFT `histology:` (`ADK` = LUAD, `SCC` = LUSC) |
| HOT phenotype (paper score, not CLDN4) | HOT 44 / COLD 48 | SOFT `hot phenotype:` |
| Sex | Male 58 / Female 34 | SOFT `Sex:` |

92 GEO arrays match the 92-patient series text. Age is `NA` for 9 arrays. Patient labels are unique `CLB_*` IDs; no extra hidden replicates in the matrix.

Do not treat this as an ICI-treated n=92. Public characteristics are diagnosis, tissue, age, sex, OS/PFS, histology, HOT score/phenotype. No ICI response labels on these GSM records.

## Panel check

HTG EdgeSeq OBP: **2560** unique gene symbols in the public matrix (matches the commercial 2560-gene OBP). Platform record is [GPL18573](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL18573) (Illumina NextSeq 500), not a gene table. The gene list is the deposited TUMADOR matrix.

Vendor OBP gene list (Omixys PDF of the HTG OBP list) also lists **CLDN3** and not **CLDN4**.

| Query | In TUMADOR OBP matrix (2560 genes) | Aliases searched |
|---|---|---|
| **CLDN4** | **no** | CLDN4, CLAUDIN, CLD4; regex `CLDN\|CLAUDIN\|CLD4\|TACSTD\|TROP` |
| CLDN3 | yes | only claudin on this panel |
| CD8A | yes | — |
| Immune (CXCL9/10, GZMB, PRF1, PDCD1) | yes | — |
| CD274 | yes | — |
| IFN (IFNG, STAT1) | yes | IFNGR1 also present |

No other claudin, TACSTD2, or TROP2 symbol is on the panel. CD8B, CD3E, HLA-C, and B2M are also absent (not required once CLDN4 is missing).

## Why stop

The assigned test is **CLDN4-only** vs CD8A / Immune / CD274 / IFN, with a LUAD vs LUSC split if labeled, honest n, no dual-high.

CLDN4 is absent, so the contrast cannot be scored from this series. **CLDN3 was not used as a stand-in.** Nearby immune genes (CD8A, CD274, IFNG) were not scored against a missing target. No extra figures.

Histology *is* labeled (ADK 47 / SCC 45), but that split is unused without CLDN4.

## Reproduce

```bash
python3 methods/gse162520_cldn4/check_panel.py
```

Writes `methods/gse162520_cldn4/panel_check.json` (downloads the GEO matrix and SOFT family file into `methods/gse162520_cldn4/_cache/`).
