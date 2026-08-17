# QUAD merge — LIANA/CellPhoneDB-style LR from CLDN4-high malignant to T/NK

**Verdict (paired n=75 honest):** On the patient-level CellPhoneDB-style score (paired n=75 across four public objects), CLDN4-high vs CLDN4-low outgoing T-recruit/IFN/MHC-I pairs: 8/27 have median Δ < 0; 13/27 reach FDR < 0.05. By axis: {'IFN': '1/1 Δ<0, median Δ=-0.063', 'MHC_I': '0/15 Δ<0, median Δ=+0.115', 'T_recruit': '7/11 Δ<0, median Δ=-0.005'}. This is a ranked co-expression list, not a recruitment mechanism.

ADDITIVE. **CLDN4 only.** TACSTD2 is not a gate. Dual-high is not run. GSE207422 is included because the user asked for it. This is **not** a re-audit of the GSE207422-only LIANA (PR #344).

## Question

Do CLDN4-high malignant cells show weaker outgoing T-recruit / MHC-I communication toward same-patient T/NK than CLDN4-low cells in the four-cohort public merge?

## What was run

| Method | Status |
|---|---|
| Documented CellPhoneDB-style score (mean of partner means; 2,920 pairs) | **primary** — patient-level paired Wilcoxon |
| LIANA `mt.cellphonedb` (per cohort, ≤1,500 cells/group, 20 permutations) | Not run. `liana` is not installed. CellPhoneDB-style score is primary. |
| CellChat | Not run. R unavailable. |
| Joint Harmony / concatenated expression LIANA | Not run. Platforms differ. |

## Honest n

The test unit is the **patient**. Cells are not n. CLDN4-high/low is the within-cohort malignant median. A patient enters if it has ≥10 malignant cells in both bins and ≥20 T/NK after pooling that patient’s tumor samples.

| Item | GSE207422 | GSE131907 | GSE148071 | GSE205335 | Merge |
|---|---:|---:|---:|---:|---:|
| Cells | 92,330 | 208,506 | 82,267 | 96,505 | 479,608 |
| Malignant | 6,627 | 31,136 | 48,118 | 28,512 | 114,393 |
| T / NK | 28,313 / 8,010 | 43,285 / 4,727 | 4,384 / 0 | 27,728 / 2,994 | — |
| CLDN4-high / low | 3,314 / 3,313 | 15,569 / 15,567 | 24,060 / 24,058 | 14,256 / 14,256 | — |
| Listed patients | 15 | 44 | 42 | 26 | 127 |
| **Paired patients (test n)** | **8** | **29** | **17** | **21** | **75** |
| Scale | log1p CP10k | log1p CP10k | TISCH2 log2(TPM/10+1) | log1p CP10k | deltas stacked |

**Paired patients (n=75):** GSE207422:P01 (BD_immune01/TN; high=386/low=182, T/NK=619); GSE207422:P03 (BD_immune03/MPR; high=227/low=519, T/NK=4905); GSE207422:P04 (BD_immune04/NMPR; high=28/low=23, T/NK=5345); GSE207422:P05 (BD_immune05/TN; high=397/low=962, T/NK=781); GSE207422:P07 (BD_immune07/NMPR; high=2036/low=1173, T/NK=984); GSE207422:P09 (BD_immune09/NMPR; high=99/low=65, T/NK=3133); GSE207422:P10 (BD_immune10/NMPR; high=70/low=119, T/NK=3876); GSE207422:P12 (BD_immune12/NMPR; high=50/low=256, T/NK=1272); GSE131907:P0006 (tLung; high=42/low=58, T/NK=1539); GSE131907:P0008 (tLung; high=20/low=22, T/NK=2128); GSE131907:P0018 (tLung; high=573/low=429, T/NK=1090); GSE131907:P0019 (tLung; high=67/low=155, T/NK=1712); GSE131907:P0020 (tLung; high=258/low=176, T/NK=2418); GSE131907:P0025 (tLung; high=82/low=141, T/NK=1599); GSE131907:P0028 (tLung; high=665/low=319, T/NK=1631); GSE131907:P0030 (tLung; high=169/low=615, T/NK=1414); GSE131907:P0031 (tLung; high=45/low=137, T/NK=3713); GSE131907:P0034 (tLung; high=999/low=1375, T/NK=366); GSE131907:P1006 (PE,tL/B; high=422/low=764, T/NK=3904); GSE131907:P1010 (mLN; high=94/low=273, T/NK=1530); GSE131907:P1011 (PE,mLN; high=73/low=130, T/NK=2757); GSE131907:P1012 (PE,mLN; high=64/low=193, T/NK=2413); GSE131907:P1015 (mLN; high=54/low=202, T/NK=479); GSE131907:P1019 (mLN; high=119/low=98, T/NK=1427); GSE131907:P1028 (tL/B; high=3586/low=1054, T/NK=232); GSE131907:P1049 (tL/B; high=17/low=97, T/NK=1126); GSE131907:P1051 (mLN; high=247/low=1038, T/NK=1143); GSE131907:P1058 (tL/B; high=51/low=409, T/NK=1546); GSE131907:P3002 (mBrain; high=19/low=162, T/NK=977); GSE131907:P3003 (mBrain; high=1904/low=809, T/NK=145); GSE131907:P3004 (mBrain; high=287/low=942, T/NK=131); GSE131907:P3006 (mBrain; high=31/low=52, T/NK=383); GSE131907:P3007 (mBrain; high=2782/low=2326, T/NK=166); GSE131907:P3012 (mBrain; high=1466/low=1345, T/NK=22); GSE131907:P3013 (mBrain; high=451/low=901, T/NK=1780); GSE131907:P3017 (mBrain; high=839/low=273, T/NK=410); GSE131907:P3019 (mBrain; high=141/low=614, T/NK=370); GSE148071:P1 (TISCH2; high=2234/low=1666, T/NK=59); GSE148071:P10 (TISCH2; high=1301/low=1309, T/NK=48); GSE148071:P13 (TISCH2; high=64/low=13, T/NK=241); GSE148071:P14 (TISCH2; high=230/low=1126, T/NK=31); GSE148071:P18 (TISCH2; high=940/low=623, T/NK=70); GSE148071:P19 (TISCH2; high=80/low=37, T/NK=30); GSE148071:P22 (TISCH2; high=123/low=16, T/NK=97); GSE148071:P23 (TISCH2; high=486/low=1291, T/NK=50); GSE148071:P24 (TISCH2; high=17/low=10, T/NK=73); GSE148071:P28 (TISCH2; high=510/low=162, T/NK=159); GSE148071:P38 (TISCH2; high=23/low=20, T/NK=177); GSE148071:P4 (TISCH2; high=98/low=691, T/NK=325); GSE148071:P40 (TISCH2; high=73/low=28, T/NK=283); GSE148071:P6 (TISCH2; high=1294/low=432, T/NK=58); GSE148071:P7 (TISCH2; high=31/low=118, T/NK=1197); GSE148071:P8 (TISCH2; high=231/low=14, T/NK=187); GSE148071:P9 (TISCH2; high=46/low=22, T/NK=67); GSE205335:P0031 (ADC/NE; high=117/low=202, T/NK=3824); GSE205335:P1006 (ADC/PR; high=735/low=547, T/NK=2998); GSE205335:P1015 (ADC/PD; high=88/low=203, T/NK=473); GSE205335:P1016 (SCLC/PR; high=2965/low=2371, T/NK=763); GSE205335:P1017 (SQ/SD; high=1141/low=1045, T/NK=2370); GSE205335:P1018 (ADC/NE; high=704/low=476, T/NK=2546); GSE205335:P1025 (SCLC/PD; high=717/low=263, T/NK=128); GSE205335:P1027 (ADC/PR; high=999/low=319, T/NK=2619); GSE205335:P1030 (ADC/PD; high=357/low=246, T/NK=140); GSE205335:P1037 (SQ/PR; high=885/low=2786, T/NK=611); GSE205335:P1056 (NUT/SD; high=733/low=1653, T/NK=1354); GSE205335:P1062 (ADC/PD; high=75/low=113, T/NK=1934); GSE205335:P1063 (ADC/NE; high=41/low=90, T/NK=3144); GSE205335:P1072 (SCLC/NE; high=135/low=271, T/NK=586); GSE205335:P1076 (ADC/PD; high=310/low=539, T/NK=2089); GSE205335:P1079 (ADC/NE; high=202/low=341, T/NK=501); GSE205335:P1084 (ADC/NE; high=157/low=39, T/NK=2350); GSE205335:P1089 (ADC/PD; high=899/low=164, T/NK=272); GSE205335:P1090 (SQ/PR; high=132/low=434, T/NK=247); GSE205335:P1115 (SCLC/PR; high=2597/low=1172, T/NK=271); GSE205335:P1119 (ADC/PD; high=261/low=961, T/NK=621).

**Dropped (not in the Wilcoxon header n):** GSE207422:P02 (malig=10, high=6, low=4, T/NK=2580); GSE207422:P06 (malig=1, high=0, low=1, T/NK=2809); GSE207422:P08 (malig=17, high=11, low=6, T/NK=2990); GSE207422:P11 (malig=0, high=0, low=0, T/NK=1637); GSE207422:P13 (malig=3, high=3, low=0, T/NK=1286); GSE207422:P14 (malig=0, high=0, low=0, T/NK=2686); GSE207422:P15 (malig=4, high=1, low=3, T/NK=1420); GSE131907:P0001 (malig=0, high=0, low=0, T/NK=0); GSE131907:P0009 (malig=5, high=2, low=3, T/NK=1981); GSE131907:P1013 (malig=376, high=0, low=376, T/NK=4318); GSE131907:P1064 (malig=0, high=0, low=0, T/NK=2777); GSE131907:P2001 (malig=0, high=0, low=0, T/NK=0); GSE131907:P2002 (malig=0, high=0, low=0, T/NK=0); GSE131907:P2003 (malig=0, high=0, low=0, T/NK=0); GSE131907:P2004 (malig=0, high=0, low=0, T/NK=0); GSE131907:P2005 (malig=0, high=0, low=0, T/NK=0); GSE131907:P2006 (malig=0, high=0, low=0, T/NK=0); GSE131907:P2007 (malig=0, high=0, low=0, T/NK=0); GSE131907:P2008 (malig=0, high=0, low=0, T/NK=0); GSE131907:P2011 (malig=0, high=0, low=0, T/NK=0); GSE131907:P2012 (malig=0, high=0, low=0, T/NK=0); GSE131907:P3016 (malig=79, high=0, low=79, T/NK=385); GSE148071:P11 (malig=62, high=2, low=60, T/NK=184); GSE148071:P12 (malig=75, high=67, low=8, T/NK=132); GSE148071:P15 (malig=2226, high=1289, low=937, T/NK=2); GSE148071:P16 (malig=950, high=755, low=195, T/NK=2); GSE148071:P17 (malig=7901, high=2699, low=5202, T/NK=0); GSE148071:P2 (malig=272, high=245, low=27, T/NK=11); GSE148071:P20 (malig=395, high=346, low=49, T/NK=4); GSE148071:P21 (malig=1811, high=1199, low=612, T/NK=11); GSE148071:P25 (malig=3138, high=1672, low=1466, T/NK=1); GSE148071:P26 (malig=24, high=19, low=5, T/NK=44); GSE148071:P27 (malig=10, high=9, low=1, T/NK=34); GSE148071:P29 (malig=815, high=234, low=581, T/NK=4); GSE148071:P3 (malig=7447, high=3057, low=4390, T/NK=13); GSE148071:P30 (malig=573, high=230, low=343, T/NK=4); GSE148071:P31 (malig=341, high=218, low=123, T/NK=4); GSE148071:P32 (malig=164, high=159, low=5, T/NK=55); GSE148071:P33 (malig=36, high=28, low=8, T/NK=3); GSE148071:P34 (malig=423, high=205, low=218, T/NK=16); GSE148071:P35 (malig=7, high=6, low=1, T/NK=62); GSE148071:P36 (malig=155, high=66, low=89, T/NK=1); GSE148071:P37 (malig=12, high=0, low=12, T/NK=48); GSE148071:P39 (malig=6, high=4, low=2, T/NK=278); GSE148071:P41 (malig=5672, high=3763, low=1909, T/NK=8); GSE148071:P42 (malig=243, high=7, low=236, T/NK=289); GSE148071:P5 (malig=1, high=0, low=1, T/NK=22); GSE205335:P2001 (malig=0, high=0, low=0, T/NK=0); GSE205335:P2009 (malig=0, high=0, low=0, T/NK=0); GSE205335:P2016 (malig=0, high=0, low=0, T/NK=0); GSE205335:P3032 (malig=0, high=0, low=0, T/NK=0); GSE205335:P4001 (malig=27, high=6, low=21, T/NK=1801).

Per-patient counts: `results/n_cells_patients.tsv`.

## Primary LR table — outgoing CLDN4-high malignant → T/NK

Median Δ = high − low. Negative = weaker from CLDN4-high. Full table: `results/lr_table_cldn4_outgoing_tnk.tsv`.

| Pathway | Pair | n patients | n cohorts | median Δ | Wilcoxon p | FDR |
|---|---|---:|---:|---:|---:|---:|
| IFN | IFNG–IFNGR1 | 3 | 2 | -0.063 | 0.25 | 0.34 |
| MHC_I | HLA-F–KIR3DL2 | 3 | 2 | +0.026 | 0.25 | 0.34 |
| MHC_I | HLA-B–CD8A | 59 | 4 | +0.082 | 5.06e-05 | 0.000229 |
| MHC_I | HLA-A–CD8A | 59 | 4 | +0.087 | 5.3e-06 | 3.43e-05 |
| MHC_I | HLA-C–CD8A | 59 | 4 | +0.091 | 2.56e-06 | 1.83e-05 |
| MHC_I | HLA-A–CD8B | 50 | 4 | +0.093 | 2.19e-05 | 0.00011 |
| MHC_I | HLA-B–CD8B | 50 | 4 | +0.096 | 0.000227 | 0.000857 |
| MHC_I | HLA-C–CD8B | 50 | 4 | +0.102 | 2.37e-05 | 0.000115 |
| MHC_I | HLA-E–KLRD1 | 57 | 4 | +0.115 | 6.92e-09 | 1.57e-07 |
| MHC_I | HLA-E–KLRC1 | 19 | 4 | +0.126 | 3.81e-06 | 2.59e-05 |
| MHC_I | HLA-E–KLRC1+KLRD1 | 17 | 4 | +0.126 | 1.53e-05 | 8.23e-05 |
| MHC_I | HLA-E–KLRC3+KLRD1 | 5 | 1 | +0.170 | 0.0625 | 0.106 |
| MHC_I | HLA-E–KLRC2 | 16 | 3 | +0.172 | 3.05e-05 | 0.000143 |
| MHC_I | HLA-E–KLRC2+KLRD1 | 15 | 3 | +0.174 | 6.1e-05 | 0.000268 |
| MHC_I | HLA-E–KLRK1 | 19 | 1 | +0.174 | 0.0141 | 0.0342 |
| MHC_I | HLA-B–KIR3DL2 | 3 | 2 | +0.184 | 0.25 | 0.34 |
| T_recruit | CXCL10–CXCR3 | 11 | 4 | -0.017 | 0.413 | 0.515 |
| T_recruit | CCL5–CCR1 | 8 | 4 | -0.016 | 0.547 | 0.63 |
| T_recruit | CCL5–CCR5 | 12 | 4 | -0.011 | 0.455 | 0.543 |
| T_recruit | CCL3–CCR5 | 5 | 3 | -0.008 | 0.438 | 0.527 |
| T_recruit | CXCL9–CXCR3 | 6 | 3 | -0.005 | 1 | 1 |
| T_recruit | CCL4–CCR5 | 8 | 4 | -0.005 | 0.383 | 0.487 |
| T_recruit | CCL3–CCR1 | 3 | 2 | -0.003 | 1 | 1 |
| T_recruit | CXCL11–CXCR3 | 5 | 3 | +0.000 | 1 | 1 |
| T_recruit | CXCL12–CXCR4 | 7 | 3 | +0.017 | 0.578 | 0.644 |
| T_recruit | CX3CL1–CX3CR1 | 7 | 3 | +0.032 | 0.0391 | 0.077 |
| T_recruit | CXCL16–CXCR6 | 39 | 4 | +0.040 | 3.06e-08 | 5.21e-07 |

8/27 focus pairs have median Δ < 0 (weaker from/to CLDN4-high). **13/27 reach FDR < 0.05**.

Incoming T/NK → malignant (focus):

| Pathway | Pair | n patients | n cohorts | median Δ | Wilcoxon p | FDR |
|---|---|---:|---:|---:|---:|---:|
| IFN | IFNG–IFNGR2 | 50 | 4 | +0.016 | 0.00151 | 0.00662 |
| IFN | IFNG–IFNGR1 | 50 | 4 | +0.017 | 0.00258 | 0.00839 |
| IFN | IFNG–IFNGR1+IFNGR2 | 50 | 4 | +0.017 | 0.00104 | 0.00508 |
| MHC_I | HLA-A–CD8B | 4 | 2 | +0.002 | 0.875 | 0.987 |
| MHC_I | HLA-B–CD8B | 4 | 2 | +0.002 | 0.875 | 0.987 |
| MHC_I | HLA-C–CD8B | 4 | 2 | +0.002 | 0.875 | 0.987 |
| MHC_I | HLA-A–CD8A | 7 | 3 | +0.004 | 1 | 1 |
| MHC_I | HLA-B–CD8A | 7 | 3 | +0.004 | 1 | 1 |
| MHC_I | HLA-C–CD8A | 7 | 3 | +0.004 | 1 | 1 |
| MHC_I | HLA-E–KLRD1 | 5 | 2 | +0.004 | 1 | 1 |

0/10 focus pairs have median Δ < 0 (weaker from/to CLDN4-high). **3/10 reach FDR < 0.05**.

## Method (short)

Partner expression = min of subunit means. Pair score = mean of the two partner means (Efremova 2020 / Garcia-Alonso 2022). `pass_expr_prop` = both partners in ≥10% of cells. Wilcoxon is paired across patients. FDR is BH within the contrast. This is not a CellChat probability and does not observe secretion or spatial contact.

LIANA was not run (`liana` package not installed). The primary method is the documented CellPhoneDB-style score.

## What this is not

- Not a re-audit of PR #344 (GSE207422-only LIANA).
- Not dual-high / TACSTD2-gated.
- Not CellChat. Not inferCNV/CopyKAT re-calls.
- Not ICI/MPR as the test (GSE131907 is treatment-naive; GSE205335 has RECIST, not MPR).
- Not a joint embedding. Expression was not concatenated.
- Cells are not the sample size.

## Extra figures

| File | Content |
|---|---|
| `results/figures/n_cells_by_patient.png` | Malignant vs T/NK per patient, by cohort |
| `results/figures/cldn4_outgoing.png` | Focus outgoing Δ |
| `results/figures/cldn4_incoming.png` | Focus incoming Δ |
| `results/figures/extra_n_by_cohort.png` | Listed vs paired patient n |
| `results/figures/extra_high_low_counts.png` | High vs low malignant n |
| `results/figures/extra_cohort_forest.png` | Cohort-stratified Δ |
| `results/figures/extra_patient_focus_heatmap.png` | Per-patient focus Δ |
| `results/figures/extra_paired_strips.png` | HLA-B–CD8A and CXCL16–CXCR6 |

## Reproduce

```bash
cd methods/quad_207422_liana_cldn4
python3 scripts/00_download.py
python3 scripts/01_build_pairs.py
python3 scripts/02_run_ccc.py
```

