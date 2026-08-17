# FINDING — CLDN4-only high-end CellChat on GSE123902 + GSE189357

**Verdict:** On the given differing pair, outgoing CellChat scores from CLDN4-high
marker-malignant cells to T/NK are **higher**, not lower, for inhibitory and barrier
pairs (LGALS9–PTPRC/CD44, NECTIN2–TIGIT, PVR–TIGIT) and for CXCL16–CXCR6.
CXCL9/10 are not detected at scale. This is **not** a recruit-down copy of the given
%pos Spearman. Honest LR n is **21** at the median split (given n=22);
LX699 fails the median high floor (7 cells). Tertile / Q4 vs Q1 extras: 21 / 21.
Between-patient Q4 tails **7/5 are thin** and are not the CellChat n.

ADDITIVE. **CLDN4 only. No dual-high.** Patient/donor is the unit.

The pair that already **differs** is **taken as given** and is not re-audited:

- PR #459 `methods/scrna_cldn4_combo_enum`, pair **GSE123902+GSE189357 %pos**
- **n=22 · ρ=−0.638 · p=0.003 · I²=0%**
- Members: GSE123902 marker-malignant donors n=13 + GSE189357 marker-malignant patients n=9
- Between-patient Q4 vs Q1 on that same %pos vector is **r=−1.000 p=0.003 on tails 7/5 — thin**. That is not the CellChat n.

This slice asks a different question: which **outgoing ligand–receptor** pairs from
**CLDN4-high vs CLDN4-low marker-malignant cells** go to **same-unit T/NK**.

## Honest n (LR, not the given Spearman)

The given n=22 is the %pos Spearman unit. LR n is the number of patients/donors
who pass the within-unit CLDN4 split and T/NK floors. Cells are not n.
Tails may be thin — said here, not hidden.

| Item | n | Note |
| --- | ---: | --- |
| Given combo (do not re-audit) | 22 | 13 + 9; %pos ρ=−0.638 |
| Given between-patient Q4 vs Q1 | 7/5 | **thin tails**; r=−1 is not CellChat n |
| Eligible vs T/NK (median) | **21** | ≥10 high/low mal + ≥20 T/NK |
| Eligible vs T/NK (tertile extra) | 21 | middle third dropped |
| Eligible vs T/NK (Q4 vs Q1 extra) | 21 | ≥40 mal; high-end |
| Median units with thin arm (<20 high or low) | 1 | kept and flagged |
| Q4 vs Q1 units with thin arm | 1 | expected on GSE123902 mets |

GSE123902 LX699 (46 malignant; median high=7) **fails** the median, tertile, and Q4 floors.
LX701 (90 malignant; high=10) is **thin** and is flagged, not dropped.
Tertile / Q4 vs Q1 arms are mutually exclusive (CLDN4=0 ties stay in the low arm).

### Per-patient floors (median split)

| cohort | patient | n_mal | n_high | n_low | n_TNK | eligible_TNK | thin_tail |
|---|---|---|---|---|---|---|---|
| GSE123902 | LX255B | 228 | 80 | 148 | 1583 | True | False |
| GSE123902 | LX653 | 145 | 72 | 73 | 49 | True | False |
| GSE123902 | LX661 | 117 | 58 | 59 | 3726 | True | False |
| GSE123902 | LX666 | 786 | 393 | 393 | 177 | True | False |
| GSE123902 | LX675 | 618 | 194 | 424 | 1336 | True | False |
| GSE123902 | LX676 | 147 | 62 | 85 | 2606 | True | False |
| GSE123902 | LX679 | 497 | 215 | 282 | 625 | True | False |
| GSE123902 | LX680 | 258 | 129 | 129 | 344 | True | False |
| GSE123902 | LX681 | 1207 | 603 | 604 | 86 | True | False |
| GSE123902 | LX682 | 421 | 79 | 342 | 2558 | True | False |
| GSE123902 | LX684 | 173 | 86 | 87 | 92 | True | False |
| GSE123902 | LX699 | 46 | 7 | 39 | 2263 | False | False |
| GSE123902 | LX701 | 90 | 10 | 80 | 322 | True | True |
| GSE189357 | TD1 | 2186 | 1093 | 1093 | 9915 | True | False |
| GSE189357 | TD2 | 1207 | 430 | 777 | 11334 | True | False |
| GSE189357 | TD3 | 664 | 319 | 345 | 6525 | True | False |
| GSE189357 | TD4 | 491 | 149 | 342 | 7096 | True | False |
| GSE189357 | TD5 | 1736 | 640 | 1096 | 11055 | True | False |
| GSE189357 | TD6 | 1211 | 605 | 606 | 3392 | True | False |
| GSE189357 | TD7 | 849 | 311 | 538 | 1851 | True | False |
| GSE189357 | TD8 | 1019 | 462 | 557 | 7257 | True | False |
| GSE189357 | TD9 | 4889 | 2444 | 2445 | 4345 | True | False |

## CellChat-style outgoing CLDN4-high → T/NK (median, patient unit)

Honest n = **21 patients/donors**. Pairs scored = 115. median Δ<0: 17; Δ>0: 98; BH-FDR<0.05: 43.

FDR < 0.05 (sorted by |median Δ|):

| ligand | receptor | ligand_class | n_patients | median_delta | pval | padj |
|---|---|---|---|---|---|---|
| LAMB3 | CD44 | other | 19 | 0.253 | 3.81e-06 | 5.21e-05 |
| APP | CD74 | other | 21 | 0.22 | 1.91e-06 | 3.91e-05 |
| SFTPD | ADGRE5 | other | 15 | 0.209 | 0.000427 | 0.00222 |
| MDK | NCL | other | 21 | 0.2 | 9.54e-07 | 3.91e-05 |
| CD69 | KLRB1 | other | 16 | -0.193 | 0.000153 | 0.000962 |
| LAMA5 | CD44 | other | 20 | 0.171 | 1.91e-06 | 3.91e-05 |
| LAMC2 | CD44 | other | 15 | 0.143 | 0.000982 | 0.00375 |
| CD55 | ADGRE5 | other | 20 | 0.138 | 1.91e-06 | 3.91e-05 |
| COL4A3 | CD44 | other | 9 | 0.122 | 0.0117 | 0.0275 |
| COL1A1 | CD44 | other | 19 | 0.119 | 0.00117 | 0.00417 |
| LGALS9 | PTPRC | inhibitory | 16 | 0.109 | 0.00539 | 0.0142 |
| COL4A4 | CD44 | other | 13 | 0.0885 | 0.000244 | 0.00143 |
| LAMB2 | CD44 | other | 18 | 0.0843 | 3.81e-05 | 0.000284 |
| LGALS9 | CD44 | inhibitory | 16 | 0.0799 | 0.00539 | 0.0142 |
| LAMA3 | CD44 | other | 12 | 0.0775 | 0.000488 | 0.00222 |
| ICAM1 | SPN | other | 16 | 0.0757 | 3.05e-05 | 0.00025 |
| GDF15 | TGFBR2 | other | 18 | 0.0729 | 7.63e-06 | 8.94e-05 |
| LAMC1 | CD44 | other | 19 | 0.0659 | 3.81e-06 | 5.21e-05 |
| CLEC2B | KLRB1 | other | 17 | -0.0652 | 0.00209 | 0.00626 |
| ICAM1 | ITGAL | other | 14 | 0.0602 | 0.000122 | 0.000834 |

Pre-specified barrier / inhibitory / recruit / attack pairs (all, not just FDR):

| ligand | receptor | ligand_class | n_patients | median_delta | pval | padj |
|---|---|---|---|---|---|---|
| LGALS9 | PTPRC | inhibitory | 16 | 0.109 | 0.00539 | 0.0142 |
| LGALS9 | CD44 | inhibitory | 16 | 0.0799 | 0.00539 | 0.0142 |
| CEACAM6 | CEACAM6 | barrier | 3 | 0.0588 | nan | nan |
| NECTIN2 | TIGIT | barrier|inhibitory | 17 | 0.0336 | 1.53e-05 | 0.000156 |
| CXCL16 | CXCR6 | recruit | 9 | 0.033 | 0.00391 | 0.011 |
| CDH1 | KLRG1 | barrier|inhibitory | 6 | 0.0327 | 0.0312 | 0.0523 |
| LGALS9 | P4HB | inhibitory | 15 | 0.0193 | 0.011 | 0.0265 |
| HLA-F | CD8A | inhibitory | 18 | 0.0128 | 0.0237 | 0.0451 |
| CEACAM5 | CEACAM6 | barrier | 3 | 0.00973 | nan | nan |
| HLA-G | CD8A | inhibitory | 9 | 0.00972 | 0.0742 | 0.111 |
| HLA-E | CD8A | inhibitory | 18 | 0.0097 | 0.702 | 0.747 |
| NECTIN2 | CD226 | barrier|inhibitory | 6 | 0.00903 | 0.0312 | 0.0523 |
| TNF | TNFRSF1B | attack | 4 | -0.00892 | nan | nan |
| CXCL13 | CXCR3 | recruit | 3 | 0.00591 | nan | nan |
| NECTIN3 | TIGIT | barrier | 10 | 0.00372 | 0.00195 | 0.00616 |
| PVR | TIGIT | barrier|inhibitory | 11 | 0.00347 | 0.000977 | 0.00375 |

## Extra: LIANA/CellPhoneDB-style outgoing CLDN4-high → T/NK (median)

Honest n = **21 patients/donors**. Pairs scored = 102. median Δ<0: 22; Δ>0: 80; BH-FDR<0.05: 32.

FDR < 0.05 (sorted by |median Δ|):

| ligand | receptor | ligand_class | n_patients | median_delta | pval | padj |
|---|---|---|---|---|---|---|
| CXCL14 | CXCR4 | other | 7 | 1.535 | 0.0156 | 0.0338 |
| CD55 | ADGRE5 | other | 20 | 0.346 | 1.91e-06 | 6.39e-05 |
| DHCR7 | RORA | other | 15 | 0.34 | 6.1e-05 | 0.000584 |
| SFTPD | ADGRE5 | other | 15 | 0.337 | 0.00061 | 0.00241 |
| PTGES2 | PTGER4 | other | 17 | 0.31 | 0.000153 | 0.000929 |
| LPAR2 | ADGRE5 | other | 7 | 0.306 | 0.0156 | 0.0338 |
| PTGES | PTGER4 | other | 7 | 0.304 | 0.0156 | 0.0338 |
| CEACAM5 | CD8A | barrier | 10 | 0.267 | 0.00195 | 0.00595 |
| CDH1 | ITGAE+ITGB7 | barrier|inhibitory | 10 | 0.257 | 0.00195 | 0.00595 |
| DHCR24 | RORA | other | 20 | 0.246 | 1.91e-06 | 6.39e-05 |
| PVR | CD96 | barrier|inhibitory | 12 | 0.24 | 0.000488 | 0.00218 |
| ALCAM | CD6 | other | 17 | 0.221 | 1.53e-05 | 0.000256 |
| ICAM1 | SPN | other | 16 | 0.22 | 3.05e-05 | 0.000409 |
| DHCR7 | NR1H2 | other | 13 | 0.218 | 0.000244 | 0.00126 |
| F11R | ITGAL+ITGB2 | barrier | 14 | 0.213 | 0.000122 | 0.000818 |
| ICAM4 | ITGAL+ITGB2 | other | 12 | 0.207 | 0.000488 | 0.00218 |
| ICAM1 | ITGAL | other | 14 | 0.199 | 0.000122 | 0.000818 |
| ICAM1 | ITGAL+ITGB2 | other | 14 | 0.199 | 0.000122 | 0.000818 |
| NECTIN3 | TIGIT | barrier | 10 | 0.188 | 0.00195 | 0.00595 |
| PVR | TIGIT | barrier|inhibitory | 11 | 0.186 | 0.000977 | 0.00363 |

Pre-specified barrier / inhibitory / recruit / attack pairs (all, not just FDR):

| ligand | receptor | ligand_class | n_patients | median_delta | pval | padj |
|---|---|---|---|---|---|---|
| TNF | TNFRSF1B | attack | 4 | -0.332 | nan | nan |
| CDH1 | KLRG1 | barrier|inhibitory | 6 | 0.331 | 0.0312 | 0.0566 |
| CEACAM6 | CEACAM6 | barrier | 3 | 0.33 | nan | nan |
| CEACAM5 | CD8A | barrier | 10 | 0.267 | 0.00195 | 0.00595 |
| CDH1 | ITGAE+ITGB7 | barrier|inhibitory | 10 | 0.257 | 0.00195 | 0.00595 |
| PVR | CD96 | barrier|inhibitory | 12 | 0.24 | 0.000488 | 0.00218 |
| CD47 | SIRPG | inhibitory | 4 | 0.218 | nan | nan |
| F11R | ITGAL+ITGB2 | barrier | 14 | 0.213 | 0.000122 | 0.000818 |
| NECTIN3 | TIGIT | barrier | 10 | 0.188 | 0.00195 | 0.00595 |
| PVR | TIGIT | barrier|inhibitory | 11 | 0.186 | 0.000977 | 0.00363 |
| TGFB2 | TGFBR3 | inhibitory | 3 | 0.18 | nan | nan |
| PVR | CD226 | barrier|inhibitory | 5 | 0.138 | nan | nan |
| CXCL16 | CXCR6 | recruit | 9 | 0.108 | 0.00391 | 0.0109 |
| CXCL10 | CXCR3 | recruit | 3 | 0.101 | nan | nan |
| NECTIN2 | TIGIT | barrier|inhibitory | 17 | 0.0948 | 1.53e-05 | 0.000256 |
| CEACAM5 | CEACAM6 | barrier | 3 | 0.0761 | nan | nan |

## Extra: CellChat tertile (high-end vs low-end) → T/NK

Honest n = **21 patients/donors**. Pairs scored = 115. median Δ<0: 18; Δ>0: 97; BH-FDR<0.05: 43.

FDR < 0.05 (sorted by |median Δ|):

| ligand | receptor | ligand_class | n_patients | median_delta | pval | padj |
|---|---|---|---|---|---|---|
| SFTPD | ADGRE5 | other | 15 | 0.347 | 0.000305 | 0.00169 |
| LAMB3 | CD44 | other | 19 | 0.265 | 3.81e-06 | 6.41e-05 |
| MDK | NCL | other | 21 | 0.245 | 9.54e-07 | 4.01e-05 |
| COL4A3 | CD44 | other | 9 | 0.241 | 0.00391 | 0.0113 |
| APP | CD74 | other | 21 | 0.223 | 1.91e-06 | 4.01e-05 |
| LAMA5 | CD44 | other | 20 | 0.196 | 1.91e-06 | 4.01e-05 |
| CD69 | KLRB1 | other | 16 | -0.194 | 0.000153 | 0.00107 |
| LAMC2 | CD44 | other | 16 | 0.151 | 0.000655 | 0.00282 |
| CD55 | ADGRE5 | other | 20 | 0.148 | 1.91e-06 | 4.01e-05 |
| LGALS9 | PTPRC | inhibitory | 17 | 0.135 | 0.00446 | 0.0121 |
| LAMC1 | CD44 | other | 17 | 0.123 | 1.53e-05 | 0.00016 |
| COL1A1 | CD44 | other | 18 | 0.117 | 0.00233 | 0.00853 |
| CLEC2B | KLRB1 | other | 18 | -0.11 | 0.000671 | 0.00282 |
| LAMB2 | CD44 | other | 19 | 0.108 | 3.81e-05 | 0.00032 |
| LGALS9 | CD44 | inhibitory | 17 | 0.106 | 0.00446 | 0.0121 |
| COL4A4 | CD44 | other | 13 | 0.101 | 0.000244 | 0.00146 |
| LAMA3 | CD44 | other | 13 | 0.091 | 0.000488 | 0.00228 |
| GDF15 | TGFBR2 | other | 18 | 0.0863 | 7.63e-06 | 0.000107 |
| THBS1 | CD47 | other | 16 | 0.0837 | 0.000214 | 0.00138 |
| ICAM1 | SPN | other | 16 | 0.0828 | 3.05e-05 | 0.000285 |

Pre-specified barrier / inhibitory / recruit / attack pairs (all, not just FDR):

| ligand | receptor | ligand_class | n_patients | median_delta | pval | padj |
|---|---|---|---|---|---|---|
| LGALS9 | PTPRC | inhibitory | 17 | 0.135 | 0.00446 | 0.0121 |
| LGALS9 | CD44 | inhibitory | 17 | 0.106 | 0.00446 | 0.0121 |
| CEACAM6 | CEACAM6 | barrier | 3 | 0.0957 | nan | nan |
| CDH1 | KLRG1 | barrier|inhibitory | 6 | 0.042 | 0.0312 | 0.0536 |
| CXCL16 | CXCR6 | recruit | 9 | 0.0353 | 0.00391 | 0.0113 |
| NECTIN2 | TIGIT | barrier|inhibitory | 17 | 0.0333 | 1.53e-05 | 0.00016 |
| HLA-G | CD8A | inhibitory | 9 | 0.0234 | 0.129 | 0.177 |
| LGALS9 | P4HB | inhibitory | 16 | 0.0231 | 0.0146 | 0.032 |
| HLA-F | CD8A | inhibitory | 18 | 0.0186 | 0.0304 | 0.0536 |
| TNF | TNFRSF1B | attack | 5 | -0.0154 | nan | nan |
| CEACAM5 | CEACAM6 | barrier | 3 | 0.00944 | nan | nan |
| HLA-G | CD8B | inhibitory | 8 | 0.0084 | 0.195 | 0.241 |
| NECTIN2 | CD226 | barrier|inhibitory | 6 | 0.00802 | 0.0312 | 0.0536 |
| PVR | TIGIT | barrier|inhibitory | 10 | 0.00714 | 0.00195 | 0.00746 |
| HLA-E | CD8A | inhibitory | 18 | 0.00637 | 0.64 | 0.672 |
| CXCL13 | CXCR3 | recruit | 3 | 0.00591 | nan | nan |

## Extra: CellChat within-unit Q4 vs Q1 → T/NK (thin tails flagged)

Honest n = **21 patients/donors**. Pairs scored = 114. median Δ<0: 17; Δ>0: 97; BH-FDR<0.05: 42.

FDR < 0.05 (sorted by |median Δ|):

| ligand | receptor | ligand_class | n_patients | median_delta | pval | padj |
|---|---|---|---|---|---|---|
| SFTPD | ADGRE5 | other | 14 | 0.339 | 0.000244 | 0.00143 |
| LAMB3 | CD44 | other | 19 | 0.299 | 0.000196 | 0.00134 |
| MDK | NCL | other | 21 | 0.284 | 9.54e-07 | 5.21e-05 |
| COL4A3 | CD44 | other | 9 | 0.259 | 0.00391 | 0.0107 |
| APP | CD74 | other | 21 | 0.255 | 1.91e-06 | 5.21e-05 |
| LAMA5 | CD44 | other | 19 | 0.221 | 3.81e-06 | 7.82e-05 |
| CD69 | KLRB1 | other | 16 | -0.182 | 3.05e-05 | 0.000278 |
| CD55 | ADGRE5 | other | 20 | 0.156 | 1.91e-06 | 5.21e-05 |
| LAMC2 | CD44 | other | 16 | 0.148 | 0.000655 | 0.00298 |
| COL1A1 | CD44 | other | 17 | 0.142 | 0.00316 | 0.00996 |
| LGALS9 | PTPRC | inhibitory | 16 | 0.134 | 0.00314 | 0.00996 |
| CLEC2B | KLRB1 | other | 18 | -0.123 | 0.00105 | 0.00408 |
| LGALS9 | CD44 | inhibitory | 16 | 0.109 | 0.00314 | 0.00996 |
| LAMC1 | CD44 | other | 17 | 0.107 | 1.53e-05 | 0.000179 |
| LAMB2 | CD44 | other | 19 | 0.105 | 5.34e-05 | 0.000438 |
| GDF15 | TGFBR2 | other | 18 | 0.0996 | 7.63e-06 | 0.000125 |
| LAMA3 | CD44 | other | 12 | 0.0987 | 0.000977 | 0.00408 |
| ICAM1 | SPN | other | 16 | 0.0866 | 3.05e-05 | 0.000278 |
| ICAM1 | ITGAL | other | 14 | 0.0725 | 0.000122 | 0.00091 |
| THBS1 | CD47 | other | 16 | 0.0724 | 0.000427 | 0.00234 |

Pre-specified barrier / inhibitory / recruit / attack pairs (all, not just FDR):

| ligand | receptor | ligand_class | n_patients | median_delta | pval | padj |
|---|---|---|---|---|---|---|
| LGALS9 | PTPRC | inhibitory | 16 | 0.134 | 0.00314 | 0.00996 |
| LGALS9 | CD44 | inhibitory | 16 | 0.109 | 0.00314 | 0.00996 |
| CEACAM6 | CEACAM6 | barrier | 3 | 0.0959 | nan | nan |
| CDH1 | KLRG1 | barrier|inhibitory | 6 | 0.0485 | 0.0312 | 0.0523 |
| CXCL16 | CXCR6 | recruit | 9 | 0.0428 | 0.00391 | 0.0107 |
| NECTIN2 | TIGIT | barrier|inhibitory | 17 | 0.0395 | 1.53e-05 | 0.000179 |
| LGALS9 | P4HB | inhibitory | 15 | 0.0265 | 0.011 | 0.0258 |
| HLA-G | CD8A | inhibitory | 9 | 0.0258 | 0.129 | 0.173 |
| HLA-F | CD8A | inhibitory | 18 | 0.0187 | 0.0304 | 0.0523 |
| TNF | TNFRSF1B | attack | 5 | -0.0154 | nan | nan |
| HLA-E | CD8A | inhibitory | 18 | 0.0114 | 0.551 | 0.594 |
| NECTIN2 | CD226 | barrier|inhibitory | 6 | 0.00984 | 0.0312 | 0.0523 |
| HLA-G | CD8B | inhibitory | 8 | 0.00869 | 0.195 | 0.243 |
| PVR | TIGIT | barrier|inhibitory | 9 | 0.00825 | 0.00391 | 0.0107 |
| CEACAM5 | CEACAM6 | barrier | 3 | 0.00816 | nan | nan |
| CXCL13 | CXCR3 | recruit | 3 | 0.00591 | nan | nan |

## Methods (short)

- Public processed GEO only. Matrices are **not** concatenated across accessions.
- Marker-malignant and T/NK gates match PR #459 (not author labels).
- CLDN4-high / low is a **within-unit** split of marker-malignant `log1p(CP10k)` CLDN4.
  Median is primary. Tertile and Q4 vs Q1 are high-end extras.
- CellChat-style: 10% truncated mean, Hill Kh=0.5, CellChatDB v2 protein pairs, `expr_prop ≥ 0.10`.
- LIANA-style extra: CellPhoneDB mean-of-means on log1p(CP10k), CellPhoneDB v5 pairs, `expr_prop ≥ 0.10`.
  LIANA package status: `not_run (ModuleNotFoundError)`.
- Inference: paired Wilcoxon on patient scores (high vs low). BH-FDR within method × rule.
- TACSTD2 is a companion gene, never a gate.

## What is not claimed

- The given %pos ρ=−0.638 is not re-computed here.
- Between-patient Q4 vs Q1 r=−1 on 7/5 is thin and is not the LR n.
- This is not a dual-high TACSTD2∩CLDN4 gate and not the 7-cohort pool.
- Cell-pooled permutation p-values are not the inferential unit; the patient/donor is.

## Extra figures

- `figures/fig_n_per_patient.png` — per-unit high/low malignant and T/NK floors
- `figures/fig_honest_n.png` — given n vs CellChat n (thin tails called out)
- `figures/fig_top_cellchat_TNK_median.png` / `fig_focus_cellchat_TNK_median.png`
- tertile and Q4 vs Q1 extras under `figures/`

## Reproduce

```bash
python3 methods/pair_123902_189357_hiend_cldn4/scripts/download.py
python3 methods/pair_123902_189357_hiend_cldn4/scripts/analyze.py
```

Primary table: [`results/lr_table.tsv`](results/lr_table.tsv).

