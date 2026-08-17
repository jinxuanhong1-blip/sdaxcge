# GSE207422 — LIANA/CellPhoneDB-style LR from CLDN4-high malignant to T/NK

**CLDN4 only.** TACSTD2 is not a gate and is not used to define dual-high. User A3 lineage and malignant-like rules are taken as given.

**Verdict (paired n=8 honest):** CLDN4-high malignant-like cells do not show a coordinated reduction of outgoing T-recruit / MHC-I communication toward T/NK. MHC-I outgoing is higher from CLDN4-high (7/7 pairs; none FDR < 0.05). T-recruit is mixed and dropout-limited (CXCR3 ligands n=3). Incoming IFNG–IFNGR is flat.

## What was run

| Method | Status |
|---|---|
| Documented CellPhoneDB-style score (mean of partner means on log1p CP10k; 2,920 pairs) | **primary** — patient-level paired Wilcoxon |
| LIANA `mt.cellphonedb` (resource `cellphonedb`, 50 permutations, ≤2,000 cells/group) | **secondary** — ran (904 edges) |
| CellChat | **not run** — R unavailable; no CellChat tables were written |

## Honest n

| Set | n |
|---|---:|
| Cells in public UMI | 92,330 |
| Samples in matrix | 15 |
| Epithelial (A3) | 11,019 |
| Malignant-like (A3) | 6,627 |
| T / NK / T+NK | 28,313 / 8,010 / 36,323 |
| CLDN4-high / low malignant-like | 3,314 / 3,313 |
| Patients in the paired LR test | **8** |

CLDN4 high = at or above the **global median** log1p(CP10k) among malignant-like cells (threshold = 1.457). Low = below. A patient enters the paired test if it has ≥10 malignant-like cells in **both** bins and ≥20 T/NK cells.

**Paired patients (n=8):** P01 (TN; high=386/low=182 malig, T/NK=619), P03 (MPR; high=227/low=519 malig, T/NK=4905), P04 (NMPR; high=28/low=23 malig, T/NK=5345), P05 (TN; high=397/low=962 malig, T/NK=781), P07 (NMPR; high=2036/low=1173 malig, T/NK=984), P09 (NMPR; high=99/low=65 malig, T/NK=3133), P10 (NMPR; high=70/low=119 malig, T/NK=3876), P12 (NMPR; high=50/low=256 malig, T/NK=1272).

**Dropped (empty or one-sided malignant bins):** P02 (NMPR; malig=10, high=6, low=4), P06 (MPR; malig=1, high=0, low=1), P08 (TN; malig=17, high=11, low=6), P11 (MPR; malig=0, high=0, low=0), P13 (NMPR; malig=3, high=3, low=0), P14 (MPR; malig=0, high=0, low=0), P15 (NMPR; malig=4, high=1, low=3). P06/P11/P14 are MPR samples with 1 / 0 / 0 malignant-like cells. The header n for the Wilcoxon is 8, not 15 and not 12.

Per-sample counts: `results/n_cells_patients.tsv`.

## Score

On log1p(CP10k), each partner’s expression is the **minimum subunit mean** (CellPhoneDB complex rule). The pair score is the **mean of the two partner means** (Efremova et al. 2020 *Nat Protoc*; Garcia-Alonso et al. 2022 *Nat Protoc*). A pair is flagged `pass_expr_prop` when both partners are detected in ≥10% of cells in their group. Patient-level tests use that patient’s own T/NK and that patient’s CLDN4-high vs CLDN4-low malignant-like cells. Cells are not treated as replicates. FDR is Benjamini–Hochberg within each contrast (outgoing or incoming).

This is **not** a CellChat communication probability.

## Primary LR table — outgoing CLDN4-high malignant → T/NK

Median patient Δ = high − low. Negative = weaker from the CLDN4-high state.

| Pathway | Pair | n patients | median Δ | Wilcoxon p | FDR |
|---|---|---:|---:|---:|---:|
| MHC_I | HLA-B–CD8A | 8 | +0.101 | 0.0156 | 0.145 |
| MHC_I | HLA-B–CD8B | 7 | +0.102 | 0.0312 | 0.169 |
| MHC_I | HLA-A–CD8A | 8 | +0.105 | 0.109 | 0.263 |
| MHC_I | HLA-E–KLRD1 | 8 | +0.109 | 0.00781 | 0.145 |
| MHC_I | HLA-C–CD8A | 8 | +0.112 | 0.0156 | 0.145 |
| MHC_I | HLA-C–CD8B | 7 | +0.129 | 0.0312 | 0.169 |
| MHC_I | HLA-A–CD8B | 7 | +0.146 | 0.0469 | 0.203 |
| T_recruit | CXCL10–CXCR3 | 3 | -0.092 | 0.5 | 0.677 |
| T_recruit | CXCL9–CXCR3 | 3 | -0.004 | 1 | 1 |
| T_recruit | CCL5–CCR5 | 4 | +0.001 | 0.875 | 0.981 |
| T_recruit | CCL4–CCR5 | 4 | +0.006 | 0.875 | 0.981 |
| T_recruit | CXCL16–CXCR6 | 7 | +0.050 | 0.0156 | 0.145 |

2/12 focus pairs have median Δ < 0 (weaker from/to CLDN4-high). **0/12 reach FDR < 0.05** in either direction.

Full ranked table (all pathways that passed filters): `results/lr_table_cldn4_outgoing_tnk.tsv`.

## Incoming T/NK → CLDN4-high vs CLDN4-low malignant (secondary)

| Pathway | Pair | n patients | median Δ | Wilcoxon p | FDR |
|---|---|---:|---:|---:|---:|
| IFN | IFNG–IFNGR2 | 7 | +0.011 | 0.578 | 0.85 |
| IFN | IFNG–IFNGR1 | 7 | +0.018 | 0.297 | 0.558 |
| IFN | IFNG–IFNGR1+IFNGR2 | 7 | +0.018 | 0.219 | 0.558 |

0/3 focus pairs have median Δ < 0 (weaker from/to CLDN4-high). **0/3 reach FDR < 0.05** in either direction.

## LIANA CellPhoneDB method (secondary, pooled / downsampled)

- import: ok
- CLDN4: LIANA_OK n_edges=904 file=liana_cellphonedb_cldn4.csv; downsampled groups={'Malig_CLDN4high': 2000, 'Malig_CLDN4low': 2000, 'NK': 2000, 'T': 2000}

Focus edges that cleared LIANA `expr_prop=0.10` from malignant-like senders to T or NK:

| Pair | target | high `lr_means` | low `lr_means` | higher in | cellphone_p high |
|---|---|---:|---:|---|---:|
| CXCL16–CXCR6 | T | 0.416 | 0.411 | high | 0 |
| CXCL16–CXCR6 | NK | 0.385 | 0.380 | high | 0 |
| CX3CL1–CX3CR1 | NK | 0.401 | 0.396 | high | 0 |
| HLA-E–KLRC1 | NK | 0.944 | 0.836 | high | 1 |
| FAM3C–HLA-C | T | 2.307 | 2.214 | high | 0 |
| FAM3C–HLA-C | NK | 2.414 | 2.321 | high | 0 |
| APP–CD74 | T | 1.864 | 1.848 | high | 0 |
| APP–CD74 | NK | 1.976 | 1.960 | high | 0 |
| CD58–CD2 | T | 1.285 | 1.290 | low | 0 |
| CD58–CD2 | NK | 0.935 | 0.940 | low | 0 |

LIANA p-values are within-object specificity, not patient-level tests. CXCL9/10/11–CXCR3 typically fail the 10% expression filter in epithelium.

## Readout

On the patient-level CellPhoneDB-style score (paired n=8), CLDN4-high vs CLDN4-low outgoing T-recruit / MHC-I pairs: **2/12 have median Δ < 0; 0/12 reach FDR < 0.05**. MHC-I: 0/7 Δ<0, median Δ=+0.109. T-recruit: 2/5 Δ<0, median Δ=+0.001 (CXCR3 ligands n=3; CXCL16–CXCR6 n=7 is higher from CLDN4-high). This is the observed rank in this public 15-sample BD Rhapsody object, not a general rule.

## Honest limits

1. Paired n = 8. The matrix has 15 samples; MPR residual tumors P06/P11/P14 are empty or one-cell under the A3 malignant-like rule. That is reported, not patched.
2. CLDN4 only. Dual-high (TACSTD2 AND CLDN4) was not run. TACSTD2 is not a gate.
3. BD Rhapsody, not 10x. GEO deposited no per-cell labels; A3 marker gates are used as given.
4. Malignant-like is a normal-lung-marker exclusion, not public CopyKAT calls.
5. Ambient RNA cannot be re-estimated from the processed matrix.
6. Chemokine dropout is high; read `n patients` and `pass_expr_prop` with the ranks.
7. CellChat was not run.

## Files

| File | Role |
|---|---|
| `FINDING.md` | This note |
| `results/lr_table_cldn4_outgoing_tnk.tsv` | Primary LR table (patient-level ranks, outgoing) |
| `results/lr_table_cldn4_incoming_tnk.tsv` | Incoming ranks |
| `results/lr_table_focus_outgoing.tsv` | T-recruit / IFN / MHC-I outgoing subset |
| `results/n_cells_patients.tsv` | Per-sample cell counts and high/low bins |
| `results/patient_cldn4_outgoing.tsv` | Per-patient pair scores (full list) |
| `results/pooled_outgoing_cldn4.tsv` | All-cell descriptive scores |
| `results/liana_cellphonedb_cldn4.csv` | Full LIANA CellPhoneDB output (if run) |
| `results/liana_cldn4_outgoing_tnk.csv` | LIANA edges malignant → T/NK |
| `results/summary.json` | Machine-readable n and method flags |
| `results/figures/` | n-cell bars and pathway Δ plots |

## Reproduce

```bash
cd methods/scrna_liana_cldn4
python3 scripts/00_download.py
python3 scripts/01_build_pairs.py
python3 scripts/02_run_ccc.py
```
