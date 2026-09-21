# GSE22493 NHEJ / STING / IFN / APM

**FINAL for this accession** (`accession_status` = `FINAL_DISCORDANT` in `key_stats.json`). The method sweep is in `method_sweep.tsv`. Do not add another analysis of GSE22493.

SKOV-3-IP-Luc ovarian arrays (GSM558700–702). Deposited series-matrix VALUE, read as log2(CLDN4 siRNA / CLDN4-overexpression control).

Primary panels are pre-specified. IFN and APM use the same gene lists as the earlier C4 slice of this accession. Set tests use genes measured on at least two arrays.

| Set | Role | Measured ≥2 / list | Median of gene-mean log2 | vs background | Mann–Whitney p |
|---|---|---:|---:|---|---:|
| NHEJ | primary | 6/8 | −0.298 | DOWN | 0.36 |
| STING | primary | 3/5 | +0.575 | UP | 0.36 |
| IFN | primary | 58/73 | −0.232 | DOWN | 0.29 |
| APM | primary | 12/16 | −0.454 | DOWN | 0.16 |
| NHEJ accessory | secondary | 5/7 | +0.584 | UP | 0.14 |
| STING regulators | secondary | 6/8 | −0.025 | FLAT | 0.50 |

Background median of gene-mean log2 = −0.027 (15,853 genes). No gene in these panels has within-panel BH q < 0.05. STING1 and TBK1 are absent from GPL10555. CLDN4 knockdown is not confirmed on the array. See `notes/w200/GSE22493_nhej_sting_ifn_apm/WRITEUP.md`.
