# Noskip blood/PBMC ICI catalog — TACSTD2 / CLDN4

No blood/PBMC/plasma/platelet ICI series was dropped because it is blood.
If a gene is missing from the processed matrix or platform, the catalog says
**ABSENT**. Association tests were run only where the gene exists **and** a
response label exists. No p-values were invented.

Primary table: `catalog_blood_ici_series.csv`  
Tests: `association_tests.csv`  
Detectability: `detectability.csv`

## Series with a response contrast (gene present)

| Series | Gene | Contrast | n | statistic | p | Honest read |
|--------|------|----------|--:|-----------|--:|-------------|
| GSE285888 | TACSTD2 | CR vs PD (patient det. rate) | 7 vs 13 | U=46.5 | 0.96 | no association |
| GSE285888 | CLDN4 | CR vs PD (patient CP10K) | 7 vs 13 | U=37 | 0.54 | no association |
| GSE111414 | TACSTD2 | PR vs PD, baseline N1 counts | 5 vs 5 | U=11 | 0.80 | near-zero counts; no association |
| GSE111414 | CLDN4 | PR vs PD, baseline N1 counts | 5 vs 5 | — | NA | **all 10 N1 counts = 0**; test undefined |
| GSE202417 | TACSTD2 | R vs NR, pre-treatment | 6 vs 8 | U=18 | 0.49 | no association |
| GSE202417 | CLDN4 | R vs NR, pre-treatment | 6 vs 8 | U=35 | 0.18 | no association |

No public **survival** (OS/PFS time) table was present in any of these GEO
series matrices, so no Cox / log-rank test was run.

## Time contrasts (not response)

| Series | Gene | Contrast | n | p | Note |
|--------|------|----------|--:|--:|------|
| GSE202417 | TACSTD2 | paired pre vs post (nivo+bezafibrate) | 14 | **0.035** | time/treatment effect only; not a response test |
| GSE202417 | CLDN4 | paired pre vs post | 14 | 0.46 | NS |
| GSE141479 | TACSTD2 | paired pre vs post nivo | 33 | 0.49 | **no response labels in GEO** |
| GSE141479 | CLDN4 | paired pre vs post nivo | 33 | 0.76 | **no response labels in GEO** |

## Genes recorded as ABSENT (series kept)

| Series | Platform | TACSTD2 | CLDN4 | Why kept |
|--------|----------|---------|-------|----------|
| GSE216297 | TEP platelet RNA-seq, 3805-gene panel | ABSENT (ENSG00000184292) | ABSENT (ENSG00000189143) | 286 samples **do** have Responder/nonResponder labels; genes are not on the panel so no test |
| GSE235048 | PBMC bulk TPM | present (median TPM 0.10; 9/15 > 0) | **ABSENT** (other CLDNs are in the matrix) | CLDN4 not dropped from catalog |
| GSE310370 | plasma miRNA 4.0 (GPL21572) | ABSENT | ABSENT | miRNA platform; 4 samples all PD |
| GSE207715 | plasma/EV miRNA | ABSENT | ABSENT | miRNA platform; 282 samples |

## Detectability where the gene is on-platform (honest lows)

| Series | Compartment | TACSTD2 | CLDN4 | Immune control |
|--------|-------------|---------|-------|----------------|
| GSE285888 | PBMC scRNA | 0.0806% of 222,144 cells | 0.1819% | PTPRC 97.8% |
| GSE213902 | PBMC T-cell 10x pools | **0 / 3903 and 0 / 5013 cells** | 1 / 3903 and 0 / 5013 | PTPRC 75–76% |
| GSE111414 | PBMC CD8 bulk counts | 4 / 20 samples > 0 (max count 9) | 2 / 20 (max 3) | — |
| GSE100860 | blood CD8 FPKM | 4 / 14 files > 0 (median 0) | 4 / 14 (median 0) | PTPRC median FPKM ~608 |
| GSE235048 | PBMC TPM | 9 / 15 > 0 (median 0.10, max 1.51) | ABSENT | PTPRC median TPM 303 |
| GSE305086 | whole-blood array | background (20th–42nd percentile) | probe-discordant | PTPRC 99th percentile |

## Reproduce

```bash
export GEO_DIR=/tmp/geo
bash   scripts/fable_blood_ici/01_download.sh
bash   scripts/fable_blood_ici/05_download_noskip.sh
python scripts/fable_blood_ici/02_extract_gse285888.py
python scripts/fable_blood_ici/03_analyze_gse285888.py
python scripts/fable_blood_ici/04_analyze_gse305086.py
python scripts/fable_blood_ici/05_noskip_blood.py
```
