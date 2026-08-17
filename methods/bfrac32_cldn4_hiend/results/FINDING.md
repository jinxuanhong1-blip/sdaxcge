# FINDING — CLDN4-only high-end LR on the given B-frac combo

**Verdict:** On the given B-frac combo, outgoing CellChat / LIANA scores from
CLDN4-high malignant cells to B, TLS-like, and T/NK are **higher**, not lower,
for inhibitory and barrier pairs (LGALS9–PTPRC/CD44, HLA-E–CD8, NECTIN2–TIGIT,
CDH1–ITGAE/ITGB7). CXCL16–CXCR6 is also up into T/NK (CellChat n=21, padj=0.003).
CXCL9/10 are not detected at scale. This is **not** a recruit-down copy of the
given B-fraction Spearman. Honest LR n is **25 (B) / 27 (TLS-like) / 28 (T/NK)**
at the median split; tertile high-end recovers **32/32** patients for T/NK.

Additive. **CLDN4 only** (no dual-high TACSTD2∩CLDN4). Patient is the unit.
The B-fraction combo that already differs is **taken as given** and is not re-audited:

- PR #290 `methods/scrna_cldn4_combo`, family `author/b/mean`
- **GSE131907 + GSE241934 IIT · k=2 · n=32 · ρ=−0.513 · p=0.00388 · I²=0%**
- Members: GSE131907 author-malignant n=21 (ρ=−0.465) + GSE241934 IIT author-epi n=11 (ρ=−0.609)
- GSE207422-only B-frac (PR #400) is NS and is **not** re-run.

This slice asks a different question: which **outgoing ligand–receptor** pairs from
**CLDN4-high vs CLDN4-low malignant cells** go to **B / TLS-like** and to **T/NK**
on those same two public cohorts.

## Honest n (LR, not the given Spearman)

The given n=32 is the B-fraction Spearman unit. LR n is the number of patients
who pass the within-patient CLDN4 split and destination floors. Cells are not n.

| Item | n | Note |
| --- | ---: | --- |
| Given combo patients (do not re-audit) | 32 | 21 + 11 |
| GSE131907 patients in matrix (author-malig samples) | 21 | sample≈patient on the 21 mets |
| GSE241934 IIT patients in matrix | 11 | author Epi |
| Patients with ≥10 high and ≥10 low malignant (median) | 28 | CLDN4-only median |
| Eligible vs B (median) | 25 | ≥10 B cells |
| Eligible vs TLS-like (median) | 27 | B + CXCL13+ T/NK; ≥10 |
| Eligible vs T/NK (median) | 28 | ≥20 T/NK |
| Patients with ≥10 high/low (tertile extra) | 32 | high-end extra |

TLS-like is a dissociated cellular proxy (B + CXCL13+ T/NK), not a histologic follicle.

Four GSE131907 patients (P1013, P1015, P1049, P3016) have **median malignant CLDN4 = 0**,
so the median rule puts every malignant cell in the high arm (low n=0). They are
dropped from the median LR, not hidden. The tertile high-end split recovers all four.
Three further patients fail the B floor (P1028 n_B=7; P3004 n_B=4; P3012 n_B=9).
IIT 11/11 pass median B and T/NK floors.

### Per-patient floors (median split)

| cohort | patient | n_mal | n_high | n_low | n_B | n_TLS | n_TNK | eligible_B | eligible_TLS | eligible_TNK |
|---|---|---|---|---|---|---|---|---|---|---|
| GSE131907 | P1006 | 1186 | 593 | 593 | 111 | 137 | 611 | True | True | True |
| GSE131907 | P1010 | 367 | 184 | 183 | 2366 | 2592 | 1530 | True | True | True |
| GSE131907 | P1011 | 203 | 102 | 101 | 1413 | 1424 | 809 | True | True | True |
| GSE131907 | P1012 | 257 | 129 | 128 | 176 | 191 | 328 | True | True | True |
| GSE131907 | P1013 | 376 | 376 | 0 | 1871 | 1976 | 1150 | False | False | False |
| GSE131907 | P1015 | 256 | 256 | 0 | 39 | 63 | 479 | False | False | False |
| GSE131907 | P1019 | 217 | 109 | 108 | 98 | 151 | 1427 | True | True | True |
| GSE131907 | P1028 | 4640 | 2320 | 2320 | 7 | 11 | 232 | False | True | True |
| GSE131907 | P1049 | 114 | 114 | 0 | 189 | 191 | 1126 | False | False | False |
| GSE131907 | P1051 | 1285 | 643 | 642 | 99 | 101 | 1143 | True | True | True |
| GSE131907 | P1058 | 460 | 230 | 230 | 162 | 320 | 1546 | True | True | True |
| GSE131907 | P3002 | 181 | 91 | 90 | 142 | 149 | 977 | True | True | True |
| GSE131907 | P3003 | 2713 | 1357 | 1356 | 16 | 24 | 145 | True | True | True |
| GSE131907 | P3004 | 1229 | 615 | 614 | 4 | 6 | 131 | False | False | True |
| GSE131907 | P3006 | 83 | 42 | 41 | 30 | 108 | 383 | True | True | True |
| GSE131907 | P3007 | 5108 | 2554 | 2554 | 54 | 59 | 166 | True | True | True |
| GSE131907 | P3012 | 2811 | 1406 | 1405 | 9 | 10 | 22 | False | True | True |
| GSE131907 | P3013 | 1352 | 676 | 676 | 656 | 867 | 1780 | True | True | True |
| GSE131907 | P3016 | 79 | 79 | 0 | 256 | 366 | 385 | False | False | False |
| GSE131907 | P3017 | 1112 | 556 | 556 | 70 | 84 | 410 | True | True | True |
| GSE131907 | P3019 | 755 | 378 | 377 | 74 | 76 | 370 | True | True | True |
| GSE241934_IIT | P343 | 188 | 94 | 94 | 648 | 1076 | 3307 | True | True | True |
| GSE241934_IIT | P438 | 63 | 32 | 31 | 935 | 1582 | 3994 | True | True | True |
| GSE241934_IIT | P519 | 41 | 21 | 20 | 2016 | 2565 | 3100 | True | True | True |
| GSE241934_IIT | P529 | 61 | 31 | 30 | 1592 | 2342 | 4761 | True | True | True |
| GSE241934_IIT | P531 | 32 | 16 | 16 | 1142 | 1544 | 4061 | True | True | True |
| GSE241934_IIT | P547 | 394 | 197 | 197 | 673 | 736 | 2680 | True | True | True |
| GSE241934_IIT | P549 | 298 | 149 | 149 | 1560 | 1665 | 2870 | True | True | True |
| GSE241934_IIT | P586 | 106 | 53 | 53 | 1937 | 2165 | 3250 | True | True | True |
| GSE241934_IIT | P589 | 368 | 184 | 184 | 137 | 163 | 3521 | True | True | True |
| GSE241934_IIT | P590 | 72 | 36 | 36 | 1561 | 1917 | 5663 | True | True | True |
| GSE241934_IIT | P591 | 76 | 38 | 38 | 855 | 945 | 3505 | True | True | True |

## CellChat-style outgoing CLDN4-high → B (median, patient unit)

Honest n = **25 patients**. Pairs scored = 75. median Δ<0: 7; Δ>0: 64; BH-FDR<0.05: 21.

FDR < 0.05 (sorted by median Δ):

| ligand | receptor | ligand_class | n_patients | median_delta | pval | padj |
|---|---|---|---|---|---|---|
| C3 | CR2 | other | 11 | -0.0531 | 0.0137 | 0.0362 |
| TNFSF13 | TNFRSF13B | other | 10 | -0.048 | 0.00977 | 0.0288 |
| LPAR2 | ADGRE5 | other | 14 | 0.00385 | 0.0107 | 0.03 |
| THBS3 | CD47 | other | 14 | 0.00876 | 0.00122 | 0.00625 |
| IGFBP3 | TMEM219 | other | 16 | 0.0242 | 0.00269 | 0.0109 |
| MDK | NCL | other | 25 | 0.0249 | 0.00203 | 0.00894 |
| LGALS9 | P4HB | inhibitory | 23 | 0.0253 | 0.000849 | 0.005 |
| APP | CD74 | other | 24 | 0.0258 | 0.000108 | 0.00108 |
| CD99 | CD99 | other | 24 | 0.0267 | 0.00481 | 0.015 |
| COL4A4 | CD44 | other | 16 | 0.0274 | 0.0182 | 0.0461 |
| LAMB2 | CD44 | other | 25 | 0.0297 | 0.0013 | 0.00625 |
| COL9A2 | CD44 | other | 14 | 0.0377 | 0.00305 | 0.0116 |
| CD55 | ADGRE5 | other | 22 | 0.0387 | 3.34e-05 | 0.000802 |
| COL6A1 | CD44 | other | 18 | 0.0391 | 0.000713 | 0.00472 |
| LAMC2 | CD44 | other | 24 | 0.0411 | 0.000108 | 0.00108 |
| LGALS9 | CD44 | inhibitory | 23 | 0.0515 | 0.00387 | 0.0137 |
| LAMB3 | CD44 | other | 25 | 0.0549 | 0.000162 | 0.00123 |
| LGALS9 | PTPRC | inhibitory | 23 | 0.0556 | 0.00434 | 0.0144 |
| LAMA5 | CD44 | other | 25 | 0.0596 | 4.54e-05 | 0.000802 |
| COL1A1 | CD44 | other | 20 | 0.0679 | 3.81e-06 | 0.000202 |

Pre-specified barrier / inhibitory / recruit / attack pairs (all, not just FDR):

| ligand | receptor | ligand_class | n_patients | median_delta | pval | padj |
|---|---|---|---|---|---|---|
| LGALS9 | PTPRC | inhibitory | 23 | 0.0556 | 0.00434 | 0.0144 |
| LGALS9 | CD44 | inhibitory | 23 | 0.0515 | 0.00387 | 0.0137 |
| LGALS9 | P4HB | inhibitory | 23 | 0.0253 | 0.000849 | 0.005 |
| TNFSF10 | TNFRSF10A | attack | 11 | 0.0154 | 0.0537 | 0.103 |
| HLA-F | LILRB1 | inhibitory | 11 | 0.00755 | 0.175 | 0.257 |
| CADM1 | CADM1 | barrier | 4 | 0.00185 | nan | nan |
| F11R | F11R | barrier | 4 | 0.00102 | nan | nan |
| HLA-G | LILRB1 | inhibitory | 5 | 0.000766 | nan | nan |
| CXCL9 | CXCR3 | recruit | 3 | 0 | nan | nan |

## CellChat-style outgoing CLDN4-high → TLS-like (median)

Honest n = **27 patients**. Pairs scored = 103. median Δ<0: 18; Δ>0: 77; BH-FDR<0.05: 26.

FDR < 0.05 (sorted by median Δ):

| ligand | receptor | ligand_class | n_patients | median_delta | pval | padj |
|---|---|---|---|---|---|---|
| TNFSF13 | TNFRSF13B | other | 10 | -0.0424 | 0.00586 | 0.0188 |
| LPAR2 | ADGRE5 | other | 16 | 0.00368 | 0.0076 | 0.0232 |
| THBS3 | CD47 | other | 16 | 0.00548 | 0.00131 | 0.00741 |
| HLA-F | CD8A | inhibitory | 12 | 0.00699 | 0.0186 | 0.0453 |
| HLA-C | CD8A | other | 13 | 0.0107 | 0.021 | 0.0493 |
| HLA-A | CD8A | other | 13 | 0.013 | 0.0122 | 0.0338 |
| HLA-B | CD8A | other | 13 | 0.013 | 0.00928 | 0.0269 |
| NECTIN2 | TIGIT | barrier|inhibitory | 7 | 0.0202 | 0.0156 | 0.0414 |
| HLA-E | CD8A | inhibitory | 13 | 0.0204 | 0.00146 | 0.00741 |
| LGALS9 | P4HB | inhibitory | 25 | 0.0214 | 0.00145 | 0.00741 |
| APP | CD74 | other | 26 | 0.0247 | 0.000591 | 0.00456 |
| IGFBP3 | TMEM219 | other | 18 | 0.0256 | 0.00158 | 0.00741 |
| COL9A2 | CD44 | other | 15 | 0.027 | 0.00336 | 0.0128 |
| CD99 | CD99 | other | 26 | 0.0272 | 0.00428 | 0.0145 |
| COL4A4 | CD44 | other | 16 | 0.0288 | 0.0182 | 0.0453 |
| MDK | NCL | other | 27 | 0.0289 | 0.00262 | 0.0106 |
| LAMB2 | CD44 | other | 27 | 0.0338 | 0.000598 | 0.00456 |
| COL6A1 | CD44 | other | 19 | 0.0376 | 0.000738 | 0.005 |
| CD55 | ADGRE5 | other | 24 | 0.039 | 9.9e-05 | 0.00101 |
| LAMC2 | CD44 | other | 26 | 0.0417 | 6.03e-05 | 0.000931 |

Pre-specified barrier / inhibitory / recruit / attack pairs (all, not just FDR):

| ligand | receptor | ligand_class | n_patients | median_delta | pval | padj |
|---|---|---|---|---|---|---|
| LGALS9 | PTPRC | inhibitory | 25 | 0.0588 | 0.00391 | 0.014 |
| LGALS9 | CD44 | inhibitory | 25 | 0.0534 | 0.00203 | 0.00882 |
| LGALS9 | P4HB | inhibitory | 25 | 0.0214 | 0.00145 | 0.00741 |
| HLA-E | CD8A | inhibitory | 13 | 0.0204 | 0.00146 | 0.00741 |
| NECTIN2 | TIGIT | barrier|inhibitory | 7 | 0.0202 | 0.0156 | 0.0414 |
| HLA-E | CD8B | inhibitory | 7 | 0.0188 | 0.0312 | 0.0561 |
| CXCL16 | CXCR6 | recruit | 6 | 0.0165 | 0.0312 | 0.0561 |
| TNFSF10 | TNFRSF10A | attack | 11 | 0.0153 | 0.042 | 0.0692 |
| CD274 | PDCD1 | inhibitory | 4 | 0.0115 | nan | nan |
| CLDN3 | CLDN3 | barrier | 3 | 0.00919 | nan | nan |
| HLA-F | CD8A | inhibitory | 12 | 0.00699 | 0.0186 | 0.0453 |
| NECTIN3 | TIGIT | barrier | 3 | 0.0042 | nan | nan |
| HLA-F | LILRB1 | inhibitory | 6 | 0.0029 | 0.562 | 0.635 |
| LGALS9 | HAVCR2 | inhibitory | 4 | 0.0017 | nan | nan |
| HLA-F | CD8B | inhibitory | 7 | 0.00128 | 0.0938 | 0.133 |
| HLA-G | CD8B | inhibitory | 3 | 0.00053 | nan | nan |

## CellChat-style outgoing CLDN4-high → T/NK (median)

Honest n = **28 patients**. Pairs scored = 105. median Δ<0: 12; Δ>0: 91; BH-FDR<0.05: 40.

FDR < 0.05 (sorted by median Δ):

| ligand | receptor | ligand_class | n_patients | median_delta | pval | padj |
|---|---|---|---|---|---|---|
| BAG6 | NCR3 | other | 12 | 0.00428 | 0.00488 | 0.0109 |
| THBS3 | CD47 | other | 16 | 0.00595 | 0.00131 | 0.00356 |
| NECTIN2 | CD226 | barrier|inhibitory | 11 | 0.00951 | 0.000977 | 0.00297 |
| LPAR2 | ADGRE5 | other | 19 | 0.0103 | 0.00141 | 0.0037 |
| CXCL16 | CXCR6 | recruit | 21 | 0.0109 | 0.00102 | 0.00298 |
| HLA-F | CD8B | inhibitory | 22 | 0.0112 | 0.0115 | 0.0236 |
| HLA-F | CD8A | inhibitory | 24 | 0.0127 | 0.00652 | 0.0141 |
| HLA-A | CD8A | other | 26 | 0.0159 | 1.33e-05 | 0.000112 |
| HLA-A | CD8B | other | 23 | 0.0187 | 0.000128 | 0.000607 |
| APP | SORL1 | other | 20 | 0.0195 | 3.81e-06 | 4.83e-05 |
| HLA-B | CD8A | other | 26 | 0.0195 | 1.64e-06 | 4.15e-05 |
| MDK | NCL | other | 28 | 0.02 | 0.00223 | 0.00531 |
| HLA-E | KLRC1 | inhibitory | 7 | 0.0207 | 0.0156 | 0.0312 |
| HLA-C | CD8A | other | 26 | 0.0217 | 2.62e-06 | 4.83e-05 |
| HLA-B | CD8B | other | 23 | 0.0237 | 4.94e-05 | 0.00025 |
| HLA-C | CD8B | other | 23 | 0.0238 | 4.03e-05 | 0.000236 |
| IGFBP3 | TMEM219 | other | 18 | 0.0258 | 0.00158 | 0.004 |
| CLEC2B | KLRB1 | other | 21 | 0.0264 | 0.0113 | 0.0236 |
| CD99 | CD99 | other | 27 | 0.0264 | 0.00325 | 0.00748 |
| ICAM1 | ITGAL | other | 21 | 0.0324 | 0.000197 | 0.000882 |

Pre-specified barrier / inhibitory / recruit / attack pairs (all, not just FDR):

| ligand | receptor | ligand_class | n_patients | median_delta | pval | padj |
|---|---|---|---|---|---|---|
| NECTIN2 | TIGIT | barrier|inhibitory | 11 | 0.0768 | 0.000977 | 0.00297 |
| LGALS9 | PTPRC | inhibitory | 26 | 0.0665 | 0.000935 | 0.00297 |
| LGALS9 | CD44 | inhibitory | 26 | 0.0651 | 0.00117 | 0.00328 |
| HLA-E | CD8A | inhibitory | 26 | 0.0443 | 9.83e-07 | 4.15e-05 |
| LGALS9 | P4HB | inhibitory | 26 | 0.0404 | 0.000835 | 0.00289 |
| HLA-E | CD8B | inhibitory | 23 | 0.0372 | 7.87e-06 | 8.54e-05 |
| CDH1 | KLRG1 | barrier|inhibitory | 19 | 0.0327 | 3.81e-06 | 4.83e-05 |
| HLA-E | KLRC2 | inhibitory | 4 | 0.0324 | nan | nan |
| NECTIN3 | TIGIT | barrier | 5 | 0.0235 | nan | nan |
| HLA-E | KLRC1 | inhibitory | 7 | 0.0207 | 0.0156 | 0.0312 |
| HLA-F | CD8A | inhibitory | 24 | 0.0127 | 0.00652 | 0.0141 |
| HLA-F | CD8B | inhibitory | 22 | 0.0112 | 0.0115 | 0.0236 |
| CXCL16 | CXCR6 | recruit | 21 | 0.0109 | 0.00102 | 0.00298 |
| NECTIN2 | CD226 | barrier|inhibitory | 11 | 0.00951 | 0.000977 | 0.00297 |
| HLA-G | CD8B | inhibitory | 6 | 0.00471 | 0.438 | 0.512 |
| TNF | TNFRSF1B | attack | 3 | 0.00443 | nan | nan |

## LIANA/CellPhoneDB-style outgoing CLDN4-high → B (median)

Honest n = **25 patients**. Pairs scored = 65. median Δ<0: 10; Δ>0: 55; BH-FDR<0.05: 8.

FDR < 0.05 (sorted by median Δ):

| ligand | receptor | ligand_class | n_patients | median_delta | pval | padj |
|---|---|---|---|---|---|---|
| ITGAV+ITGB1 | ADGRE5 | other | 22 | 0.0194 | 0.00593 | 0.0341 |
| LGALS9 | P4HB | inhibitory | 23 | 0.0395 | 0.00387 | 0.0257 |
| IGFBP3 | TMEM219 | other | 16 | 0.0609 | 0.00269 | 0.0257 |
| APP | CD74 | other | 24 | 0.0716 | 8.34e-06 | 0.000192 |
| CDH1 | ITGAE+ITGB7 | barrier|inhibitory | 20 | 0.104 | 1.91e-06 | 8.77e-05 |
| HLA-E | VSIR | inhibitory | 10 | 0.119 | 0.00391 | 0.0257 |
| PODXL2 | SELL | other | 19 | 0.14 | 0.00284 | 0.0257 |
| CD55 | ADGRE5 | other | 22 | 0.159 | 6.53e-05 | 0.001 |

Pre-specified barrier / inhibitory / recruit / attack pairs (all, not just FDR):

| ligand | receptor | ligand_class | n_patients | median_delta | pval | padj |
|---|---|---|---|---|---|---|
| HLA-E | VSIR | inhibitory | 10 | 0.119 | 0.00391 | 0.0257 |
| CXCL9 | CXCR3 | recruit | 3 | -0.119 | nan | nan |
| HLA-G | LILRB1 | inhibitory | 5 | 0.118 | nan | nan |
| CDH1 | ITGAE+ITGB7 | barrier|inhibitory | 20 | 0.104 | 1.91e-06 | 8.77e-05 |
| F11R | ITGAL+ITGB2 | barrier | 8 | 0.0923 | 0.0156 | 0.0653 |
| CADM1 | CADM1 | barrier | 4 | 0.0846 | nan | nan |
| TGFB2 | TGFBR1+TGFBR2 | inhibitory | 6 | 0.084 | 0.0938 | 0.216 |
| TNFSF10 | TNFRSF10A | attack | 11 | 0.0634 | 0.083 | 0.212 |
| TNFSF10 | TNFRSF10B | attack | 12 | 0.0581 | 0.0771 | 0.212 |
| LGALS9 | P4HB | inhibitory | 23 | 0.0395 | 0.00387 | 0.0257 |
| HLA-F | VSIR | inhibitory | 10 | 0.0296 | 0.232 | 0.371 |
| HLA-F | LILRB1 | inhibitory | 11 | 0.0259 | 0.24 | 0.371 |
| TGFB1 | TGFBR1+TGFBR2 | inhibitory | 10 | 0.00872 | 0.625 | 0.75 |

## LIANA/CellPhoneDB-style outgoing CLDN4-high → TLS-like (median)

Honest n = **27 patients**. Pairs scored = 87. median Δ<0: 15; Δ>0: 72; BH-FDR<0.05: 14.

FDR < 0.05 (sorted by median Δ):

| ligand | receptor | ligand_class | n_patients | median_delta | pval | padj |
|---|---|---|---|---|---|---|
| ITGAV+ITGB1 | ADGRE5 | other | 24 | 0.0194 | 0.00652 | 0.0316 |
| LPAR2 | ADGRE5 | other | 16 | 0.0195 | 0.00418 | 0.0243 |
| LGALS9 | P4HB | inhibitory | 25 | 0.0395 | 0.00203 | 0.0172 |
| IGFBP3 | TMEM219 | other | 18 | 0.0609 | 0.00105 | 0.0142 |
| APP | CD74 | other | 26 | 0.0644 | 1.33e-05 | 0.000453 |
| PODXL2 | SELL | other | 19 | 0.0815 | 0.00284 | 0.0214 |
| F11R | ITGAL+ITGB2 | barrier | 11 | 0.0928 | 0.00195 | 0.0172 |
| HLA-C | CD8A | other | 13 | 0.103 | 0.00464 | 0.0243 |
| CDH1 | ITGAE+ITGB7 | barrier|inhibitory | 20 | 0.104 | 1.91e-06 | 0.00013 |
| HLA-A | CD8A | other | 13 | 0.11 | 0.00464 | 0.0243 |
| HLA-B | CD8A | other | 13 | 0.117 | 0.00342 | 0.0232 |
| HLA-E | VSIR | inhibitory | 11 | 0.121 | 0.00195 | 0.0172 |
| ALCAM | CD6 | other | 11 | 0.123 | 0.000977 | 0.0142 |
| CD55 | ADGRE5 | other | 24 | 0.159 | 3.02e-05 | 0.000684 |

Pre-specified barrier / inhibitory / recruit / attack pairs (all, not just FDR):

| ligand | receptor | ligand_class | n_patients | median_delta | pval | padj |
|---|---|---|---|---|---|---|
| NECTIN3 | TIGIT | barrier | 3 | 0.15 | nan | nan |
| CXCL9 | CXCR3 | recruit | 4 | -0.144 | nan | nan |
| CD274 | PDCD1 | inhibitory | 4 | 0.128 | nan | nan |
| HLA-E | VSIR | inhibitory | 11 | 0.121 | 0.00195 | 0.0172 |
| CDH1 | ITGAE+ITGB7 | barrier|inhibitory | 20 | 0.104 | 1.91e-06 | 0.00013 |
| F11R | ITGAL+ITGB2 | barrier | 11 | 0.0928 | 0.00195 | 0.0172 |
| TGFB2 | TGFBR1+TGFBR2 | inhibitory | 6 | 0.0828 | 0.0938 | 0.178 |
| CXCL16 | CXCR6 | recruit | 6 | 0.0796 | 0.0312 | 0.085 |
| CEACAM5 | CD8A | barrier | 8 | 0.0668 | 0.0234 | 0.0722 |
| TNFSF10 | TNFRSF10B | attack | 11 | 0.0634 | 0.083 | 0.168 |
| TNFSF10 | TNFRSF10A | attack | 11 | 0.0634 | 0.083 | 0.168 |
| NECTIN2 | TIGIT | barrier|inhibitory | 7 | 0.0571 | 0.0156 | 0.0531 |
| CD47 | SIRPG | inhibitory | 10 | 0.052 | 0.131 | 0.234 |
| LGALS9 | P4HB | inhibitory | 25 | 0.0395 | 0.00203 | 0.0172 |
| LGALS9 | HAVCR2 | inhibitory | 4 | 0.031 | nan | nan |
| HLA-F | VSIR | inhibitory | 11 | 0.0259 | 0.206 | 0.334 |

## LIANA/CellPhoneDB-style outgoing CLDN4-high → T/NK (median)

Honest n = **28 patients**. Pairs scored = 111. median Δ<0: 14; Δ>0: 97; BH-FDR<0.05: 38.

FDR < 0.05 (sorted by median Δ):

| ligand | receptor | ligand_class | n_patients | median_delta | pval | padj |
|---|---|---|---|---|---|---|
| PTGES2 | PTGER4 | other | 24 | 0.00964 | 0.015 | 0.0356 |
| LPAR2 | ADGRE5 | other | 19 | 0.0123 | 0.00169 | 0.00554 |
| COL4A4 | ITGA1+ITGB1 | other | 12 | 0.0166 | 0.00684 | 0.0192 |
| ITGAV+ITGB1 | ADGRE5 | other | 27 | 0.0191 | 0.007 | 0.0192 |
| JAG1 | CD46 | other | 17 | 0.0282 | 0.0129 | 0.0331 |
| CXCL16 | CXCR6 | recruit | 21 | 0.0292 | 0.00118 | 0.004 |
| ICAM4 | ITGAL+ITGB2 | other | 19 | 0.0308 | 0.000523 | 0.00234 |
| COL1A1 | ITGA1+ITGB1 | other | 11 | 0.0374 | 0.00684 | 0.0192 |
| ALCAM | CD6 | other | 26 | 0.0381 | 4.17e-07 | 1.77e-05 |
| LGALS9 | P4HB | inhibitory | 26 | 0.0386 | 0.00104 | 0.0037 |
| DHCR24 | RORA | other | 27 | 0.0393 | 0.0151 | 0.0356 |
| LGALS9 | HAVCR2 | inhibitory | 14 | 0.0405 | 0.0134 | 0.0336 |
| PTGES | PTGER4 | other | 18 | 0.0485 | 0.000839 | 0.0034 |
| NECTIN2 | TIGIT | barrier|inhibitory | 11 | 0.0571 | 0.000977 | 0.0037 |
| TNC | ITGA4+ITGB1 | other | 16 | 0.0603 | 0.0214 | 0.0479 |
| IGFBP3 | TMEM219 | other | 18 | 0.0609 | 0.00105 | 0.0037 |
| APP | CD74 | other | 27 | 0.0629 | 7.99e-06 | 7.54e-05 |
| CD47 | SIRPG | inhibitory | 20 | 0.0685 | 0.000261 | 0.00139 |
| APP | SORL1 | other | 20 | 0.0776 | 3.81e-06 | 4.63e-05 |
| F11R | ITGAL+ITGB2 | barrier | 21 | 0.0819 | 1.91e-06 | 4.05e-05 |

Pre-specified barrier / inhibitory / recruit / attack pairs (all, not just FDR):

| ligand | receptor | ligand_class | n_patients | median_delta | pval | padj |
|---|---|---|---|---|---|---|
| NECTIN3 | TIGIT | barrier | 5 | 0.325 | nan | nan |
| CXCL9 | CXCR3 | recruit | 4 | -0.277 | nan | nan |
| CEACAM5 | CD8A | barrier | 15 | 0.135 | 0.00061 | 0.00259 |
| TGFB2 | TGFBR1+TGFBR2 | inhibitory | 5 | 0.122 | nan | nan |
| HLA-E | VSIR | inhibitory | 11 | 0.121 | 0.00195 | 0.00615 |
| CDH1 | KLRG1 | barrier|inhibitory | 19 | 0.119 | 3.81e-06 | 4.63e-05 |
| HLA-E | KLRD1 | inhibitory | 25 | 0.115 | 8.34e-07 | 2.36e-05 |
| HLA-E | KLRC2+KLRD1 | inhibitory | 4 | 0.11 | nan | nan |
| HLA-E | KLRC2 | inhibitory | 4 | 0.11 | nan | nan |
| HLA-E | KLRC1+KLRD1 | inhibitory | 6 | 0.108 | 0.0312 | 0.0627 |
| PVR | CD226 | barrier|inhibitory | 9 | -0.103 | 0.0547 | 0.0989 |
| CDH1 | ITGAE+ITGB7 | barrier|inhibitory | 24 | 0.103 | 1.19e-07 | 1.01e-05 |
| HLA-E | KLRC1 | inhibitory | 7 | 0.101 | 0.0156 | 0.0359 |
| TGFB2 | TGFBR3 | inhibitory | 6 | 0.0993 | 0.0938 | 0.15 |
| F11R | ITGAL+ITGB2 | barrier | 21 | 0.0819 | 1.91e-06 | 4.05e-05 |
| TNFSF10 | TNFRSF10A | attack | 8 | 0.0818 | 0.0781 | 0.13 |

## Extra: CellChat tertile (high-end vs low-end) → B

Honest n = **29 patients**. Pairs scored = 76. median Δ<0: 5; Δ>0: 70; BH-FDR<0.05: 22.

FDR < 0.05 (sorted by median Δ):

| ligand | receptor | ligand_class | n_patients | median_delta | pval | padj |
|---|---|---|---|---|---|---|
| C3 | CR2 | other | 11 | -0.0857 | 0.0137 | 0.0329 |
| PPIA | BSG | other | 29 | -0.00905 | 0.0129 | 0.0327 |
| LPAR2 | ADGRE5 | other | 19 | 0.00177 | 0.00247 | 0.0069 |
| THBS3 | CD47 | other | 16 | 0.0114 | 0.00629 | 0.0167 |
| LAMB2 | CD44 | other | 29 | 0.0244 | 0.000888 | 0.00262 |
| LAMC1 | CD44 | other | 27 | 0.0266 | 0.000744 | 0.00232 |
| LGALS9 | P4HB | inhibitory | 24 | 0.0289 | 7.63e-05 | 0.000505 |
| MDK | NCL | other | 29 | 0.0301 | 0.000332 | 0.00126 |
| APP | CD74 | other | 28 | 0.032 | 6.74e-06 | 8.92e-05 |
| IGFBP3 | TMEM219 | other | 19 | 0.0402 | 0.000126 | 0.000607 |
| CD55 | ADGRE5 | other | 25 | 0.0451 | 5.25e-06 | 8.92e-05 |
| LAMC2 | CD44 | other | 26 | 0.0454 | 7.02e-05 | 0.000505 |
| LAMA5 | CD44 | other | 29 | 0.0471 | 4.7e-05 | 0.000498 |
| COL9A2 | CD44 | other | 16 | 0.0499 | 6.1e-05 | 0.000505 |
| COL6A2 | CD44 | other | 12 | 0.0621 | 0.000488 | 0.00164 |
| COL6A1 | CD44 | other | 21 | 0.0641 | 1.91e-06 | 5.05e-05 |
| LAMB3 | CD44 | other | 27 | 0.0649 | 0.000162 | 0.000717 |
| CD99 | CD99 | other | 28 | 0.0692 | 0.00012 | 0.000607 |
| LGALS9 | PTPRC | inhibitory | 24 | 0.0705 | 0.000494 | 0.00164 |
| LGALS9 | CD44 | inhibitory | 24 | 0.0747 | 0.000278 | 0.00113 |

Pre-specified barrier / inhibitory / recruit / attack pairs (all, not just FDR):

| ligand | receptor | ligand_class | n_patients | median_delta | pval | padj |
|---|---|---|---|---|---|---|
| LGALS9 | CD44 | inhibitory | 24 | 0.0747 | 0.000278 | 0.00113 |
| LGALS9 | PTPRC | inhibitory | 24 | 0.0705 | 0.000494 | 0.00164 |
| LGALS9 | P4HB | inhibitory | 24 | 0.0289 | 7.63e-05 | 0.000505 |
| TNFSF10 | TNFRSF10A | attack | 11 | 0.0169 | 0.0674 | 0.108 |
| HLA-F | LILRB1 | inhibitory | 14 | 0.00841 | 0.0419 | 0.0765 |
| CADM1 | CADM1 | barrier | 4 | 0.00306 | nan | nan |
| F11R | F11R | barrier | 4 | 0.00141 | nan | nan |
| HLA-G | LILRB1 | inhibitory | 5 | 0.000394 | nan | nan |

## Extra: CellChat tertile → T/NK

Honest n = **32 patients**. Pairs scored = 103. median Δ<0: 7; Δ>0: 96; BH-FDR<0.05: 46.

FDR < 0.05 (sorted by median Δ):

| ligand | receptor | ligand_class | n_patients | median_delta | pval | padj |
|---|---|---|---|---|---|---|
| PPIA | BSG | other | 32 | -0.00804 | 0.0184 | 0.0331 |
| LGALS9 | HAVCR2 | inhibitory | 15 | 0.00185 | 0.0181 | 0.0331 |
| LPAR1 | ADGRE5 | other | 9 | 0.0041 | 0.00781 | 0.0149 |
| HLA-DRB1 | CD4 | other | 20 | 0.00684 | 0.00365 | 0.0081 |
| BAG6 | NCR3 | other | 14 | 0.0081 | 0.00061 | 0.00156 |
| CLEC2D | KLRB1 | other | 14 | 0.00811 | 0.00525 | 0.0102 |
| HLA-F | CD8B | inhibitory | 26 | 0.0116 | 0.000664 | 0.00165 |
| CXCL16 | CXCR6 | recruit | 22 | 0.0117 | 0.0021 | 0.00478 |
| TGM2 | ADGRG1 | other | 11 | 0.012 | 0.0186 | 0.0331 |
| THBS3 | CD47 | other | 17 | 0.0125 | 0.00385 | 0.00821 |
| NECTIN2 | CD226 | barrier|inhibitory | 11 | 0.013 | 0.00488 | 0.00977 |
| LPAR2 | ADGRE5 | other | 24 | 0.0135 | 0.000916 | 0.00221 |
| HLA-F | CD8A | inhibitory | 28 | 0.014 | 0.000381 | 0.00104 |
| HLA-E | KLRC1 | inhibitory | 9 | 0.0156 | 0.00391 | 0.00821 |
| HLA-A | CD8B | other | 27 | 0.025 | 1.35e-05 | 7.08e-05 |
| HLA-A | CD8A | other | 30 | 0.0253 | 1.68e-06 | 1.73e-05 |
| HLA-B | CD8B | other | 27 | 0.0268 | 8.2e-07 | 1.22e-05 |
| HLA-B | CD8A | other | 30 | 0.0277 | 6.15e-08 | 2.8e-06 |
| MDK | NCL | other | 32 | 0.0298 | 0.000306 | 0.000866 |
| HLA-C | CD8B | other | 27 | 0.031 | 1.04e-06 | 1.22e-05 |

Pre-specified barrier / inhibitory / recruit / attack pairs (all, not just FDR):

| ligand | receptor | ligand_class | n_patients | median_delta | pval | padj |
|---|---|---|---|---|---|---|
| LGALS9 | PTPRC | inhibitory | 27 | 0.114 | 4.09e-05 | 0.000168 |
| LGALS9 | CD44 | inhibitory | 27 | 0.0951 | 3.52e-05 | 0.000152 |
| NECTIN2 | TIGIT | barrier|inhibitory | 11 | 0.0949 | 0.00488 | 0.00977 |
| HLA-E | CD8A | inhibitory | 30 | 0.0628 | 1.02e-07 | 2.8e-06 |
| CDH1 | KLRG1 | barrier|inhibitory | 21 | 0.0624 | 0.000293 | 0.000857 |
| HLA-E | KLRC2 | inhibitory | 4 | 0.0559 | nan | nan |
| LGALS9 | P4HB | inhibitory | 27 | 0.0539 | 7.29e-05 | 0.00026 |
| HLA-E | CD8B | inhibitory | 27 | 0.0535 | 1.04e-06 | 1.22e-05 |
| NECTIN3 | TIGIT | barrier | 6 | 0.0285 | 0.0312 | 0.0545 |
| HLA-E | KLRC1 | inhibitory | 9 | 0.0156 | 0.00391 | 0.00821 |
| HLA-F | CD8A | inhibitory | 28 | 0.014 | 0.000381 | 0.00104 |
| NECTIN2 | CD226 | barrier|inhibitory | 11 | 0.013 | 0.00488 | 0.00977 |
| CXCL16 | CXCR6 | recruit | 22 | 0.0117 | 0.0021 | 0.00478 |
| HLA-F | CD8B | inhibitory | 26 | 0.0116 | 0.000664 | 0.00165 |
| TNF | TNFRSF1B | attack | 5 | 0.00984 | nan | nan |
| TNF | TNFRSF1A | attack | 3 | 0.0079 | nan | nan |

## Methods (short)

- Public processed GEO only. Matrices are **not** concatenated across accessions.
- CLDN4-high / low is a **within-patient** split of malignant (or IIT author-epi) `log1p(CP10k)` CLDN4.
  Median is primary. Tertile is the high-end extra (middle third dropped).
- CellChat-style: 10% truncated mean, Hill \(K_h=0.5\), CellChatDB v2 protein pairs, `expr_prop ≥ 0.10`.
- LIANA-style: CellPhoneDB mean-of-means on log1p(CP10k), CellPhoneDB v5 pairs, `expr_prop ≥ 0.10`.
  LIANA package status: `not_run (ModuleNotFoundError)`.
- Inference: paired Wilcoxon on patient scores (high vs low). BH-FDR within method × destination × rule.
- TACSTD2 is a companion gene, never a gate.

## What is not claimed

- The given B-fraction ρ=−0.513 is not re-computed here.
- This is not GSE207422 and not an ICI/MPR contrast on GSE131907 (treatment-naive atlas).
- TLS-like is not a follicle. CopyKAT was not re-run.
- Cell-pooled permutation p-values are not the inferential unit; the patient is.

## Reproduce

```bash
python3 methods/bfrac32_cldn4_hiend/scripts/download.py --outdir /tmp/bfrac32_cldn4_hiend
python3 methods/bfrac32_cldn4_hiend/scripts/analyze.py --data /tmp/bfrac32_cldn4_hiend --out methods/bfrac32_cldn4_hiend
```

Primary table: [`results/lr_table.tsv`](results/lr_table.tsv).

