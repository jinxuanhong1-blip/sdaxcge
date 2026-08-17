# Triple merge GSE131907+GSE148071+GSE205335 — CLDN4-only LIANA/CellPhoneDB

**CLDN4 only.** TACSTD2 is not a gate and is not used to define dual-high.
**Not GSE207422.** **Not** the 131907+205335-only pair merge.

**Verdict (paired n=72 honest):** MHC-I outgoing is higher from CLDN4-high malignant cells (0/15 pairs Δ<0; median Δ=+0.128; 12/15 FDR<0.05). T-recruit is mixed and dropout-limited (5/8 Δ<0; median Δ=−0.006; CXCL9/10 n≤8). CXCL16–CXCR6 is higher from CLDN4-high (n=37; FDR=1.4e-6). IFN outgoing is n=3. This is the observed rank in this public triple merge, not a general rule.

## What was run

| Method | Status |
|---|---|
| Documented CellPhoneDB-style score (mean of partner means on log1p CP10k; 2920 pairs) | **primary** — patient-level paired Wilcoxon |
| LIANA `mt.cellphonedb` (resource `cellphonedb`, 50 permutations, ≤2,000 cells/group) | **secondary** — ran_cellphonedb_method |
| CellChat | **not run** — R unavailable; no CellChat tables were written |

## Honest n

| Set | n |
|---|---:|
| Cells kept (malignant or T/NK) | 193294 |
| Malignant / malignant-like | 101801 |
| T / NK / T+NK | 81947 / 9546 / 91493 |
| CLDN4-high / low malignant | 50904 / 50897 |
| Patients in the extract | 101 |
| Patients in the paired LR test | **72** |

| Dataset | Deposited in extract | Paired | Malignant cells | T/NK |
|---|---:|---:|---:|---:|
| GSE131907 | 33 | 29 | 31136 | 48012 |
| GSE148071 | 42 | 22 | 42153 | 3606 |
| GSE205335 | 26 | 21 | 28512 | 39875 |

CLDN4 high = at or above the **global median** log1p(CP10k) among malignant cells in the merge (threshold = 1.464). Low = below. A patient enters the paired test if it has ≥10 malignant cells in **both** bins and ≥20 T/NK cells.

**Paired patients (n=72):** GSE131907|P0006 (GSE131907; high=68/low=32 malig, T/NK=1539), GSE131907|P0008 (GSE131907; high=29/low=13 malig, T/NK=2128), GSE131907|P0018 (GSE131907; high=773/low=229 malig, T/NK=1090), GSE131907|P0019 (GSE131907; high=119/low=103 malig, T/NK=1712), GSE131907|P0020 (GSE131907; high=329/low=105 malig, T/NK=2418), GSE131907|P0025 (GSE131907; high=119/low=104 malig, T/NK=1599), GSE131907|P0028 (GSE131907; high=819/low=165 malig, T/NK=1631), GSE131907|P0030 (GSE131907; high=281/low=503 malig, T/NK=1414), GSE131907|P0031 (GSE131907; high=79/low=103 malig, T/NK=3713), GSE131907|P0034 (GSE131907; high=1486/low=888 malig, T/NK=366), GSE131907|P1006 (GSE131907; high=609/low=577 malig, T/NK=3904), GSE131907|P1010 (GSE131907; high=147/low=220 malig, T/NK=1530), GSE131907|P1011 (GSE131907; high=109/low=94 malig, T/NK=2757), GSE131907|P1012 (GSE131907; high=125/low=132 malig, T/NK=2413), GSE131907|P1015 (GSE131907; high=75/low=181 malig, T/NK=479), GSE131907|P1019 (GSE131907; high=169/low=48 malig, T/NK=1427), GSE131907|P1028 (GSE131907; high=4117/low=523 malig, T/NK=232), GSE131907|P1049 (GSE131907; high=33/low=81 malig, T/NK=1126), GSE131907|P1051 (GSE131907; high=500/low=785 malig, T/NK=1143), GSE131907|P1058 (GSE131907; high=106/low=354 malig, T/NK=1546), GSE131907|P3002 (GSE131907; high=45/low=136 malig, T/NK=977), GSE131907|P3003 (GSE131907; high=2245/low=468 malig, T/NK=145), GSE131907|P3004 (GSE131907; high=614/low=615 malig, T/NK=131), GSE131907|P3006 (GSE131907; high=40/low=43 malig, T/NK=383), GSE131907|P3007 (GSE131907; high=3725/low=1383 malig, T/NK=166), GSE131907|P3012 (GSE131907; high=1968/low=843 malig, T/NK=22), GSE131907|P3013 (GSE131907; high=667/low=685 malig, T/NK=1780), GSE131907|P3017 (GSE131907; high=922/low=190 malig, T/NK=410), GSE131907|P3019 (GSE131907; high=264/low=491 malig, T/NK=370), GSE148071|P1 (GSE148071; high=1401/low=2035 malig, T/NK=140), GSE148071|P10 (GSE148071; high=856/low=1515 malig, T/NK=58), GSE148071|P13 (GSE148071; high=44/low=14 malig, T/NK=77), GSE148071|P14 (GSE148071; high=107/low=754 malig, T/NK=27), GSE148071|P18 (GSE148071; high=688/low=739 malig, T/NK=85), GSE148071|P19 (GSE148071; high=34/low=23 malig, T/NK=47), GSE148071|P22 (GSE148071; high=112/low=30 malig, T/NK=193), GSE148071|P23 (GSE148071; high=337/low=1188 malig, T/NK=25), GSE148071|P24 (GSE148071; high=34/low=48 malig, T/NK=54), GSE148071|P26 (GSE148071; high=37/low=22 malig, T/NK=65), GSE148071|P27 (GSE148071; high=17/low=16 malig, T/NK=79), GSE148071|P28 (GSE148071; high=235/low=141 malig, T/NK=62), GSE148071|P3 (GSE148071; high=2309/low=4701 malig, T/NK=24), GSE148071|P32 (GSE148071; high=163/low=29 malig, T/NK=74), GSE148071|P38 (GSE148071; high=27/low=46 malig, T/NK=146), GSE148071|P4 (GSE148071; high=19/low=212 malig, T/NK=281), GSE148071|P40 (GSE148071; high=31/low=24 malig, T/NK=330), GSE148071|P5 (GSE148071; high=40/low=77 malig, T/NK=45), GSE148071|P6 (GSE148071; high=902/low=447 malig, T/NK=75), GSE148071|P7 (GSE148071; high=28/low=107 malig, T/NK=772), GSE148071|P8 (GSE148071; high=67/low=15 malig, T/NK=55), GSE148071|P9 (GSE148071; high=985/low=882 malig, T/NK=65), GSE205335|P0031 (GSE205335; high=106/low=213 malig, T/NK=5182), GSE205335|P1006 (GSE205335; high=687/low=595 malig, T/NK=2998), GSE205335|P1015 (GSE205335; high=83/low=208 malig, T/NK=473), GSE205335|P1016 (GSE205335; high=2768/low=2568 malig, T/NK=763), GSE205335|P1017 (GSE205335; high=1103/low=1083 malig, T/NK=2370), GSE205335|P1018 (GSE205335; high=682/low=498 malig, T/NK=2546), GSE205335|P1025 (GSE205335; high=686/low=294 malig, T/NK=128), GSE205335|P1027 (GSE205335; high=994/low=324 malig, T/NK=2619), GSE205335|P1030 (GSE205335; high=352/low=251 malig, T/NK=140), GSE205335|P1037 (GSE205335; high=703/low=2968 malig, T/NK=611), GSE205335|P1056 (GSE205335; high=685/low=1701 malig, T/NK=1354), GSE205335|P1062 (GSE205335; high=72/low=116 malig, T/NK=1934), GSE205335|P1063 (GSE205335; high=39/low=92 malig, T/NK=3144), GSE205335|P1072 (GSE205335; high=122/low=284 malig, T/NK=586), GSE205335|P1076 (GSE205335; high=274/low=575 malig, T/NK=2089), GSE205335|P1079 (GSE205335; high=190/low=353 malig, T/NK=501), GSE205335|P1084 (GSE205335; high=156/low=40 malig, T/NK=2350), GSE205335|P1089 (GSE205335; high=881/low=182 malig, T/NK=272), GSE205335|P1090 (GSE205335; high=124/low=442 malig, T/NK=247), GSE205335|P1115 (GSE205335; high=2454/low=1315 malig, T/NK=271), GSE205335|P1119 (GSE205335; high=250/low=972 malig, T/NK=621).

**Dropped (empty or one-sided malignant bins, or thin T/NK):** GSE131907|P0009 (GSE131907; high=3/low=2 malig, T/NK=1981), GSE131907|P1013 (GSE131907; high=0/low=376 malig, T/NK=4318), GSE131907|P1064 (GSE131907; high=0/low=0 malig, T/NK=2777), GSE131907|P3016 (GSE131907; high=0/low=79 malig, T/NK=385), GSE148071|P11 (GSE148071; high=1/low=76 malig, T/NK=39), GSE148071|P12 (GSE148071; high=20/low=9 malig, T/NK=55), GSE148071|P15 (GSE148071; high=986/low=952 malig, T/NK=14), GSE148071|P16 (GSE148071; high=611/low=252 malig, T/NK=1), GSE148071|P17 (GSE148071; high=1288/low=4383 malig, T/NK=3), GSE148071|P2 (GSE148071; high=1/low=4 malig, T/NK=0), GSE148071|P20 (GSE148071; high=206/low=47 malig, T/NK=0), GSE148071|P21 (GSE148071; high=662/low=704 malig, T/NK=16), GSE148071|P25 (GSE148071; high=1415/low=1864 malig, T/NK=2), GSE148071|P29 (GSE148071; high=124/low=365 malig, T/NK=2), GSE148071|P30 (GSE148071; high=139/low=384 malig, T/NK=3), GSE148071|P31 (GSE148071; high=166/low=133 malig, T/NK=3), GSE148071|P33 (GSE148071; high=1/low=1 malig, T/NK=0), GSE148071|P34 (GSE148071; high=127/low=252 malig, T/NK=15), GSE148071|P35 (GSE148071; high=38/low=52 malig, T/NK=15), GSE148071|P36 (GSE148071; high=47/low=119 malig, T/NK=3), GSE148071|P37 (GSE148071; high=0/low=13 malig, T/NK=38), GSE148071|P39 (GSE148071; high=5/low=4 malig, T/NK=209), GSE148071|P41 (GSE148071; high=2591/low=2444 malig, T/NK=3), GSE148071|P42 (GSE148071; high=2/low=127 malig, T/NK=406), GSE205335|P2001 (GSE205335; high=0/low=0 malig, T/NK=2365), GSE205335|P2009 (GSE205335; high=0/low=0 malig, T/NK=606), GSE205335|P2016 (GSE205335; high=0/low=0 malig, T/NK=3652), GSE205335|P3032 (GSE205335; high=0/low=0 malig, T/NK=252), GSE205335|P4001 (GSE205335; high=5/low=22 malig, T/NK=1801).

The header n for the Wilcoxon is the paired count, not deposited n and not the number of GEO samples.
Per-sample counts: `results/n_cells_patients.tsv`.

## Score

On log1p(CP10k), each partner’s expression is the **minimum subunit mean** (CellPhoneDB complex rule). The pair score is the **mean of the two partner means** (Efremova et al. 2020 *Nat Protoc*; Garcia-Alonso et al. 2022 *Nat Protoc*). A pair is flagged `pass_expr_prop` when both partners are detected in ≥10% of cells in their group. Patient-level tests use that patient’s own T/NK and that patient’s CLDN4-high vs CLDN4-low malignant cells. Cells are not treated as replicates. FDR is Benjamini–Hochberg within each contrast.

This is **not** a CellChat communication probability.

## Primary LR table — outgoing CLDN4-high malignant → T/NK

Median patient Δ = high − low. Negative = weaker from the CLDN4-high state.

| Pathway | Pair | n patients | median Δ | Wilcoxon p | FDR |
|---|---|---:|---:|---:|---:|
| IFN | IFNG–IFNGR1 | 3 | -0.061 | 0.25 | 0.32 |
| MHC_I | HLA-F–KIR3DL2 | 3 | +0.030 | 0.25 | 0.32 |
| MHC_I | HLA-C–CD8A | 59 | +0.102 | 1.58e-06 | 1.07e-05 |
| MHC_I | HLA-B–CD8A | 59 | +0.109 | 3.65e-05 | 0.000141 |
| MHC_I | HLA-E–KLRD1 | 55 | +0.113 | 4.26e-08 | 5.97e-07 |
| MHC_I | HLA-C–CD8B | 48 | +0.114 | 4.26e-06 | 2.18e-05 |
| MHC_I | HLA-A–CD8A | 59 | +0.115 | 2.06e-06 | 1.2e-05 |
| MHC_I | HLA-A–CD8B | 48 | +0.118 | 6.15e-06 | 2.83e-05 |
| MHC_I | HLA-E–KLRC1 | 15 | +0.128 | 6.1e-05 | 0.000223 |
| MHC_I | HLA-B–CD8B | 48 | +0.137 | 0.000104 | 0.000352 |
| MHC_I | HLA-E–KLRC1+KLRD1 | 13 | +0.142 | 0.000244 | 0.000727 |
| MHC_I | HLA-E–KLRC3+KLRD1 | 6 | +0.156 | 0.0625 | 0.103 |
| MHC_I | HLA-E–KLRK1 | 19 | +0.176 | 0.0141 | 0.0316 |
| MHC_I | HLA-B–KIR3DL2 | 3 | +0.178 | 0.25 | 0.32 |
| MHC_I | HLA-E–KLRC2 | 14 | +0.182 | 0.000244 | 0.000727 |
| MHC_I | HLA-E–KLRC2+KLRD1 | 14 | +0.182 | 0.000244 | 0.000727 |
| T_recruit | CCL4–CCR5 | 4 | -0.028 | 0.25 | 0.32 |
| T_recruit | CCL5–CCR1 | 7 | -0.018 | 0.578 | 0.691 |
| T_recruit | CXCL10–CXCR3 | 8 | -0.009 | 0.547 | 0.66 |
| T_recruit | CCL3–CCR5 | 4 | -0.007 | 0.625 | 0.721 |
| T_recruit | CCL5–CCR5 | 11 | -0.006 | 0.765 | 0.849 |
| T_recruit | CXCL12–CXCR4 | 6 | +0.009 | 0.0312 | 0.0588 |
| T_recruit | CX3CL1–CX3CR1 | 7 | +0.031 | 0.0156 | 0.0328 |
| T_recruit | CXCL16–CXCR6 | 37 | +0.043 | 1.45e-07 | 1.43e-06 |

6/24 focus pairs have median Δ < 0 (weaker from/to CLDN4-high). **14/24 reach FDR < 0.05** in either direction.

Full ranked table (all pathways that passed filters): `results/lr_table_cldn4_outgoing_tnk.tsv`.

## Incoming T/NK → CLDN4-high vs CLDN4-low malignant (secondary)

| Pathway | Pair | n patients | median Δ | Wilcoxon p | FDR |
|---|---|---:|---:|---:|---:|
| IFN | IFNG–IFNGR1 | 47 | +0.020 | 6.39e-05 | 0.000371 |
| IFN | IFNG–IFNGR1+IFNGR2 | 44 | +0.024 | 0.000179 | 0.000917 |
| IFN | IFNG–IFNGR2 | 46 | +0.038 | 3.6e-06 | 2.85e-05 |
| MHC_I | HLA-A–CD8A | 4 | -0.045 | 0.625 | 0.745 |
| MHC_I | HLA-B–CD8A | 4 | -0.045 | 0.625 | 0.745 |
| MHC_I | HLA-C–CD8A | 4 | -0.045 | 0.625 | 0.745 |
| MHC_I | HLA-E–KLRD1 | 3 | -0.015 | 0.5 | 0.659 |
| MHC_I | HLA-E–KLRC2 | 3 | -0.010 | 0.75 | 0.87 |
| MHC_I | HLA-A–CD8B | 4 | +0.006 | 0.875 | 0.952 |
| MHC_I | HLA-B–CD8B | 4 | +0.006 | 0.875 | 0.952 |
| MHC_I | HLA-C–CD8B | 4 | +0.006 | 0.875 | 0.952 |

5/11 focus pairs have median Δ < 0 (weaker from/to CLDN4-high). **3/11 reach FDR < 0.05** in either direction.

## LIANA CellPhoneDB method (secondary, pooled / downsampled)

- CLDN4: LIANA_OK n_edges=525 file=liana_cellphonedb_cldn4.csv; downsampled groups={'Malig_CLDN4high': 2000, 'Malig_CLDN4low': 2000, 'NK': 2000, 'T': 2000}

Focus edges that cleared LIANA `expr_prop=0.10` from malignant senders to T/NK:

| Pair | target | high `lr_means` | low `lr_means` | higher in |
|---|---|---:|---:|---|
| APP–CD74 | NK | 1.238 | 1.135 | high |
| APP–CD74 | T | 1.254 | 1.151 | high |
| CD58–CD2 | NK | 0.525 | 0.519 | high |
| CD58–CD2 | T | 0.841 | 0.836 | high |
| CX3CL1–CX3CR1 | NK | 0.291 | nan | low |
| CXCL16–CXCR6 | T | 0.241 | 0.201 | high |
| HLA-B–KIR3DL2 | NK | 0.935 | 0.904 | high |
| HLA-C–KIR2DL3 | NK | 0.860 | 0.796 | high |
| HLA-E–KLRC1 | NK | 0.697 | 0.674 | high |
| HLA-E–KLRC1_KLRD1 | NK | 0.697 | 0.674 | high |
| HLA-E–KLRC2 | NK | 0.656 | 0.633 | high |
| HLA-E–KLRC2_KLRD1 | NK | 0.656 | 0.633 | high |

LIANA p-values are within-object specificity, not patient-level tests.

## Readout

On the patient-level CellPhoneDB-style score (paired n=72), MHC-I outgoing is higher from CLDN4-high (0/15 Δ<0; median Δ=+0.128; 12/15 FDR<0.05). T-recruit: 5/8 Δ<0, median Δ=−0.006 (CXCR3 ligands n=8; CXCL16–CXCR6 n=37 is higher from high). IFN outgoing n=3. 6/24 focus pairs have median Δ<0; 14/24 reach FDR<0.05, driven by MHC-I up rather than T-recruit down. This is the observed rank in this public triple merge, not a general rule.

## Honest limits

1. Paired n = 72. Deposited patients in the extract = 101. One-sided or empty malignant bins are reported, not patched.
2. CLDN4 only. Dual-high (TACSTD2 AND CLDN4) was not run. TACSTD2 is not a gate.
3. GSE207422 is excluded. This is not the 131907+205335-only pair.
4. Labels are mixed by design: author malignant/T/NK in GSE131907 and GSE205335; A3 marker-argmax malignant-like in GSE148071 (no GEO labels).
5. GSE131907 uses tumor-site cells only. Normal-lung epithelium is not counted as malignant.
6. Ambient RNA cannot be re-estimated from the processed matrices.
7. Chemokine dropout is high; read `n patients` and `pass_expr_prop` with the ranks.
8. CellChat was not run. LIANA p-values are not patient-level tests.
9. The three series are different platforms and clinical settings. The merge is additive, not batch-corrected Harmony/scVI.

## Files

| File | Role |
|---|---|
| `FINDING.md` | This note |
| `results/lr_table_cldn4_outgoing_tnk.tsv` | Primary LR table (patient-level ranks, outgoing) |
| `results/lr_table_cldn4_incoming_tnk.tsv` | Incoming ranks |
| `results/lr_table_focus_outgoing.tsv` | T-recruit / IFN / MHC-I outgoing subset |
| `results/n_cells_patients.tsv` | Per-patient cell counts and high/low bins |
| `results/patient_cldn4_outgoing.tsv.gz` | Per-patient pair scores |
| `results/summary.json` | Machine-readable n and method flags |
| `results/figures/` | Honest-n, density, focus Δ, extra ligand table |

## Reproduce

```bash
cd methods/triple_scrna_liana_cldn4
python3 scripts/00_download.py
python3 scripts/01_build_pairs.py
python3 scripts/02_extract.py
python3 scripts/03_run_ccc.py
```
