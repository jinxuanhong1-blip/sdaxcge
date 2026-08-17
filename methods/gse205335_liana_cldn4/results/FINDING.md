# GSE205335 — LIANA/CellPhoneDB-style LR from CLDN4-high malignant to T/NK

**CLDN4 only.** TACSTD2 is not a gate and is not used to define dual-high. Author `lineage.sub` / `lineage.total` labels are used as given. This is not GSE207422 (PR #344).

**Verdict (paired n=21 honest):** CLDN4-high author-malignant cells do not show a coordinated reduction of outgoing T-recruit / MHC-I communication toward same-patient T/NK. MHC-I outgoing median-of-pair Δ = +0.170 (0/13 Δ<0). T-recruit median-of-pair Δ = +0.066 (1/3 Δ<0). **6/16 focus pairs reach FDR < 0.05.**

## What was run

| Method | Status |
|---|---|
| Documented CellPhoneDB-style score (mean of partner means on log1p CP10k; 2,920 pairs) | **primary** — patient-level paired Wilcoxon |
| LIANA `mt.cellphonedb` (resource `cellphonedb`, 50 permutations, ≤2000 cells/group) | **secondary** — ran_cellphonedb_method |
| CellChat | **not run** — R unavailable; no CellChat tables were written |

## Honest n

| Set | n |
|---|---:|
| Cells in public UMI | 96,505 |
| GEO samples / patients | 33 / 26 |
| Tumor samples / patients (normals dropped) | 28 / 22 |
| Author malignant (tumor samples) | 28,512 |
| T / NK / T+NK | 27,728 / 2,994 / 31,642 |
| CLDN4-high / low malignant | 14,256 / 14,256 |
| Patients in the paired LR test | **21** |

CLDN4 high = at or above the **global median** log1p(CP10k) among author-malignant cells (threshold = 1.371). Low = below. A patient enters the paired test if it has ≥10 malignant cells in **both** bins and ≥20 T/NK cells after pooling that patient’s **tumor** samples.

**Paired patients (n=21):** P0031 (ADC/NE; high=117/low=202 malig, T/NK=3824), P1006 (ADC/PR; high=735/low=547 malig, T/NK=2998), P1015 (ADC/PD; high=88/low=203 malig, T/NK=473), P1016 (SCLC/PR; high=2965/low=2371 malig, T/NK=763), P1017 (SQ/SD; high=1141/low=1045 malig, T/NK=2370), P1018 (ADC/NE; high=704/low=476 malig, T/NK=2546), P1025 (SCLC/PD; high=717/low=263 malig, T/NK=128), P1027 (ADC/PR; high=999/low=319 malig, T/NK=2619), P1030 (ADC/PD; high=357/low=246 malig, T/NK=140), P1037 (SQ/PR; high=885/low=2786 malig, T/NK=611), P1056 (NUT/SD; high=733/low=1653 malig, T/NK=1354), P1062 (ADC/PD; high=75/low=113 malig, T/NK=1934), P1063 (ADC/NE; high=41/low=90 malig, T/NK=3144), P1072 (SCLC/NE; high=135/low=271 malig, T/NK=586), P1076 (ADC/PD; high=310/low=539 malig, T/NK=2089), P1079 (ADC/NE; high=202/low=341 malig, T/NK=501), P1084 (ADC/NE; high=157/low=39 malig, T/NK=2350), P1089 (ADC/PD; high=899/low=164 malig, T/NK=272), P1090 (SQ/PR; high=132/low=434 malig, T/NK=247), P1115 (SCLC/PR; high=2597/low=1172 malig, T/NK=271), P1119 (ADC/PD; high=261/low=961 malig, T/NK=621).

**Dropped from the paired test:** P4001 (ADC/SD; malig=27, high=6, low=21, T/NK=1801).

Normal-only GEO patients P2001 / P2009 / P2016 (Normal LN) and P3032 (Normal Brain) have 0 author-malignant cells and are out of the tumor extract. P0031 Normal Lung is excluded; P0031 tumor lung is kept. MPR/NMPR is unlabeled; RECIST is not used as MPR. The header n for the Wilcoxon is the paired count, not 26 and not 22.

Per-patient counts: `results/n_cells_patients.tsv`.

## Score

On log1p(CP10k), each partner’s expression is the **minimum subunit mean** (CellPhoneDB complex rule). The pair score is the **mean of the two partner means** (Efremova et al. 2020 *Nat Protoc*; Garcia-Alonso et al. 2022 *Nat Protoc*). A pair is flagged `pass_expr_prop` when both partners are detected in ≥10% of cells in their group. Patient-level tests use that patient’s own T/NK and that patient’s CLDN4-high vs CLDN4-low malignant cells. Cells are not treated as replicates. FDR is Benjamini–Hochberg within each contrast (outgoing or incoming).

This is **not** a CellChat communication probability.

## Primary LR table — outgoing CLDN4-high malignant → T/NK

Median patient Δ = high − low. Negative = weaker from the CLDN4-high state.

| Pathway | Pair | n patients | median Δ | Wilcoxon p | FDR |
|---|---|---:|---:|---:|---:|
| MHC_I | HLA-B–CD8A | 18 | +0.123 | 0.0539 | 0.0778 |
| MHC_I | HLA-A–CD8B | 17 | +0.132 | 0.0448 | 0.0672 |
| MHC_I | HLA-C–CD8A | 18 | +0.148 | 0.0385 | 0.0589 |
| MHC_I | HLA-A–CD8A | 18 | +0.149 | 0.0304 | 0.0483 |
| MHC_I | HLA-B–CD8B | 17 | +0.151 | 0.0714 | 0.0944 |
| MHC_I | HLA-C–CD8B | 17 | +0.158 | 0.0505 | 0.0744 |
| MHC_I | HLA-E–KLRC3+KLRD1 | 5 | +0.170 | 0.0625 | 0.0841 |
| MHC_I | HLA-E–KLRC1 | 5 | +0.174 | 0.0625 | 0.0841 |
| MHC_I | HLA-E–KLRC1+KLRD1 | 5 | +0.174 | 0.0625 | 0.0841 |
| MHC_I | HLA-E–KLRK1 | 19 | +0.174 | 0.0141 | 0.0274 |
| MHC_I | HLA-E–KLRD1 | 18 | +0.175 | 0.00193 | 0.00564 |
| MHC_I | HLA-E–KLRC2 | 10 | +0.190 | 0.00195 | 0.00564 |
| MHC_I | HLA-E–KLRC2+KLRD1 | 10 | +0.190 | 0.00195 | 0.00564 |
| T_recruit | CCL5–CCR5 | 3 | -0.106 | 0.25 | 0.287 |
| T_recruit | CX3CL1–CX3CR1 | 3 | +0.066 | 0.25 | 0.287 |
| T_recruit | CXCL16–CXCR6 | 11 | +0.075 | 0.00195 | 0.00564 |

1/16 focus pairs have median Δ < 0 (weaker from/to CLDN4-high). **6/16 reach FDR < 0.05** in either direction.

CXCL9/10/11–CXCR3 did not enter the paired table (fewer than 3 patients passed `expr_prop` on either arm). That dropout is the n, not a hidden negative.

Full ranked table (all pathways that passed filters): `results/lr_table_cldn4_outgoing_tnk.tsv`.

## Incoming T/NK → CLDN4-high vs CLDN4-low malignant (secondary)

| Pathway | Pair | n patients | median Δ | Wilcoxon p | FDR |
|---|---|---:|---:|---:|---:|
| IFN | IFNG–IFNGR1 | 13 | +0.060 | 0.000244 | 0.00127 |
| IFN | IFNG–IFNGR1+IFNGR2 | 13 | +0.060 | 0.000244 | 0.00127 |
| IFN | IFNG–IFNGR2 | 13 | +0.063 | 0.000488 | 0.00177 |

0/3 focus pairs have median Δ < 0 (weaker from/to CLDN4-high). **3/3 reach FDR < 0.05** in either direction.

## LIANA CellPhoneDB method (secondary, pooled / downsampled)

- import: ok
- CLDN4: LIANA_OK n_edges=938 file=liana_cellphonedb_cldn4.csv; downsampled groups={'Malig_CLDN4high': 2000, 'Malig_CLDN4low': 2000, 'NK': 2000, 'T': 2000}

Focus edges that cleared LIANA `expr_prop=0.10` from malignant senders to T or NK:

| Pair | target | high `lr_means` | low `lr_means` | higher in | cellphone_p high |
|---|---|---:|---:|---|---:|
| CXCL16–CXCR6 | T | 0.224 | 0.182 | high | 0 |
| CX3CL1–CX3CR1 | NK | 0.441 | NA | high only | 0 |
| HLA-E–KLRC1 | NK | 0.652 | 0.625 | high | 1 |
| MIF–CD74 | T | 2.253 | 2.287 | low | 0 |
| MIF–CD74 | NK | 2.201 | 2.235 | low | 0 |
| FAM3C–HLA-C | NK | 2.102 | 2.018 | high | 0 |
| FAM3C–HLA-C | T | 1.876 | 1.793 | high | 0 |
| MIF–TNFRSF14 | T | 1.530 | 1.564 | low | 0 |
| MIF–TNFRSF14 | NK | 1.502 | 1.537 | low | 0 |
| APP–CD74 | T | 1.485 | 1.286 | high | 0 |
| APP–CD74 | NK | 1.432 | 1.234 | high | 0 |
| COPA–CD74 | T | 1.277 | 1.175 | high | 0 |
| COPA–CD74 | NK | 1.224 | 1.123 | high | 0 |
| MDK–SORL1 | NK | 1.090 | 0.799 | high | 0 |
| HLA-E–KLRK1 | NK | 0.978 | 0.951 | high | 1 |

LIANA p-values are within-object specificity, not patient-level tests. CXCL9/10/11–CXCR3 typically fail the 10% expression filter in epithelium.

## Readout

On the patient-level CellPhoneDB-style score (paired n=21), CLDN4-high vs CLDN4-low outgoing T-recruit / MHC-I / IFN pairs: 1/16 have median Δ < 0; 6/16 reach FDR < 0.05. MHC_I: 0/13 Δ<0, median Δ=+0.170; T_recruit: 1/3 Δ<0, median Δ=+0.066. This is the observed rank in this public GSE205335 tumor extract, not a general rule.

## Honest limits

1. Paired n = 21. GEO has 26 patients / 33 samples; four normal-only patients have 0 malignant cells. That is reported, not patched.
2. CLDN4 only. Dual-high (TACSTD2 AND CLDN4) was not run. TACSTD2 is not a gate.
3. Author malignant labels, not public CopyKAT / inferCNV calls.
4. Histology is mixed (ADC / SQ / SCLC / NUT). SCLC is not hidden.
5. MPR is unlabeled. RECIST is not MPR.
6. Ambient RNA cannot be re-estimated from the processed matrix.
7. Chemokine dropout is high; read `n patients` and `pass_expr_prop` with the ranks.
8. CellChat was not run. This is not the Q4 vs Q1 between-patient table in PR #362.
9. This is not GSE207422 (PR #344).

## Extra figures

| File | Content |
|---|---|
| `results/figures/n_cells_by_patient.png` | Author malignant vs T/NK per patient |
| `results/figures/cldn4_outgoing.png` | Focus-axis outgoing Δ |
| `results/figures/cldn4_incoming.png` | Focus-axis incoming Δ |
| `results/figures/extra_cldn4_distribution.png` | Malignant CLDN4 high/low histogram |
| `results/figures/extra_patient_focus_heatmap.png` | Per-patient MHC-I / T-recruit Δ |
| `results/figures/extra_paired_strips.png` | Paired high vs low scores for key pairs |
| `results/figures/extra_paired_by_histology.png` | Paired vs dropped by histology |
| `results/figures/extra_high_low_counts.png` | High vs low malignant n per patient |

## Files

| File | Role |
|---|---|
| `FINDING.md` | This note |
| `results/lr_table_cldn4_outgoing_tnk.tsv` | Primary LR table (patient-level ranks, outgoing) |
| `results/lr_table_cldn4_incoming_tnk.tsv` | Incoming ranks |
| `results/lr_table_focus_outgoing.tsv` | T-recruit / IFN / MHC-I outgoing subset |
| `results/n_cells_patients.tsv` | Per-patient cell counts and high/low bins |
| `results/patient_cldn4_outgoing.tsv.gz` | Per-patient pair scores (full list) |
| `results/pooled_outgoing_cldn4.tsv` | All-cell descriptive scores |
| `results/liana_cellphonedb_cldn4.csv` | Full LIANA CellPhoneDB output (if run) |
| `results/summary.json` | Machine-readable n and method flags |
| `results/figures/` | n-cell bars, pathway Δ, and extra figures |

## Reproduce

```bash
cd methods/gse205335_liana_cldn4
python3 scripts/00_download.py
python3 scripts/01_build_pairs.py
python3 scripts/02_run_ccc.py
```
