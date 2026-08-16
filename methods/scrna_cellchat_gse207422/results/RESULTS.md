# GSE207422 malignant TACSTD2 / CLDN4 ligand–receptor communication

Public processed UMI from Hu et al., *Genome Med* 2023 (GSE207422; PMID 36869384).
User A3 lineage and malignant-like rules are used as given. This note reports
outgoing signals from malignant-like cells to T/NK and incoming signals from T/NK
to malignant-like cells, split on TACSTD2 and (separately) CLDN4.

## What was run

| Method | Status |
|---|---|
| Documented CellPhoneDB-style score (mean of partner means on log1p CP10k; 2,920 pairs) | **primary** — patient-level paired Wilcoxon |
| LIANA `mt.cellphonedb` (resource `cellphonedb`, 50 permutations, ≤2,000 cells/group) | **secondary** — ran (923 TACSTD2 edges; 904 CLDN4 edges) |
| CellChat | **not run** — R unavailable; no CellChat tables were written |

## n cells / patients

| Set | n |
|---|---:|
| Cells in matrix | 92,330 |
| Patients | 15 |
| Epithelial (A3) | 11,019 |
| Malignant-like (A3) | 6,627 |
| T / NK / T+NK | 28,313 / 8,010 / 36,323 |
| TACSTD2-high / low malignant-like | 3,314 / 3,313 (median log1p CP10k = 2.027) |
| CLDN4-high / low malignant-like | 3,314 / 3,313 (median log1p CP10k = 1.457) |
| Patients in the paired test | **8** (same 8 for TACSTD2 and CLDN4) |

Paired patients (≥10 malignant-like cells in both bins and ≥20 T/NK):
**P01 (TN), P03 (MPR), P04 (NMPR), P05 (TN), P07 (NMPR), P09 (NMPR), P10 (NMPR), P12 (NMPR).**

Dropped for empty or one-sided malignant bins: P02, P06, P08, P11, P13, P14, P15
(P06/P11/P14 are MPR samples with 1 / 0 / 0 malignant-like cells). Per-sample
counts: `n_cells_patients.tsv`.

## Score

On log1p(CP10k), each partner’s expression is the **minimum subunit mean**
(CellPhoneDB complex rule). The pair score is the **mean of the two partner means**
(Efremova et al. 2020 *Nat Protoc*; Garcia-Alonso et al. 2022 *Nat Protoc*).
A pair is flagged `pass_expr_prop` when both partners are detected in ≥10% of
cells in their group. Patient-level tests use that patient’s own T/NK and that
patient’s high vs low malignant-like cells. Cells are not treated as replicates.
FDR is Benjamini–Hochberg within each contrast.

## TACSTD2-high vs TACSTD2-low

Median patient Δ = high − low. Negative = weaker from/to the high state.

### Outgoing malignant → T/NK

| Pathway | Pair | n patients | median Δ | Wilcoxon p | FDR |
|---|---|---:|---:|---:|---:|
| MHC_I | HLA-A–CD8A | 8 | +0.061 | 0.25 | 0.414 |
| MHC_I | HLA-A–CD8B | 7 | +0.072 | 0.109 | 0.313 |
| MHC_I | HLA-B–CD8A | 8 | +0.093 | 0.078 | 0.259 |
| MHC_I | HLA-C–CD8A | 8 | +0.122 | 0.016 | 0.123 |
| MHC_I | HLA-B–CD8B | 7 | +0.136 | 0.109 | 0.313 |
| MHC_I | HLA-C–CD8B | 7 | +0.137 | 0.031 | 0.176 |
| MHC_I | HLA-E–KLRD1 | 8 | +0.144 | 0.0078 | 0.098 |
| T_recruit | CCL5–CCR5 | 4 | −0.020 | 0.25 | 0.414 |
| T_recruit | CCL4–CCR5 | 4 | −0.001 | 0.875 | 1 |
| T_recruit | CXCL10–CXCR3 | 3 | +0.010 | 1 | 1 |
| T_recruit | CXCL9–CXCR3 | 3 | +0.028 | 0.50 | 0.643 |
| T_recruit | CXCL16–CXCR6 | 7 | +0.111 | 0.016 | 0.123 |

2/12 focus pairs have median Δ < 0. **0/12 reach FDR < 0.05.**
MHC-I is uniformly higher from TACSTD2-high (7/7; median of pair Δ = +0.122).
CXCR3-axis ligands are sparse (n=3 patients pass `expr_prop` on either side).
CXCL16–CXCR6 is the strongest T-recruit edge and is **higher** from TACSTD2-high.

IFN ligands are T/NK products, so they do not appear as malignant outgoing pairs.

### Incoming T/NK → malignant

| Pathway | Pair | n patients | median Δ | Wilcoxon p | FDR |
|---|---|---:|---:|---:|---:|
| IFN | IFNG–IFNGR1+IFNGR2 | 7 | +0.002 | 1 | 1 |
| IFN | IFNG–IFNGR1 | 7 | +0.021 | 0.109 | 0.456 |
| IFN | IFNG–IFNGR2 | 7 | +0.023 | 0.688 | 0.849 |

Incoming IFNG–IFNGR is flat (0/3 Δ < 0; 0 FDR < 0.05).

## CLDN4-high vs CLDN4-low

### Outgoing malignant → T/NK

| Pathway | Pair | n patients | median Δ | Wilcoxon p | FDR |
|---|---|---:|---:|---:|---:|
| MHC_I | HLA-B–CD8A | 8 | +0.101 | 0.016 | 0.145 |
| MHC_I | HLA-B–CD8B | 7 | +0.102 | 0.031 | 0.169 |
| MHC_I | HLA-A–CD8A | 8 | +0.105 | 0.109 | 0.263 |
| MHC_I | HLA-E–KLRD1 | 8 | +0.109 | 0.0078 | 0.145 |
| MHC_I | HLA-C–CD8A | 8 | +0.112 | 0.016 | 0.145 |
| MHC_I | HLA-C–CD8B | 7 | +0.129 | 0.031 | 0.169 |
| MHC_I | HLA-A–CD8B | 7 | +0.146 | 0.047 | 0.203 |
| T_recruit | CXCL10–CXCR3 | 3 | −0.092 | 0.50 | 0.677 |
| T_recruit | CXCL9–CXCR3 | 3 | −0.004 | 1 | 1 |
| T_recruit | CCL5–CCR5 | 4 | +0.001 | 0.875 | 0.981 |
| T_recruit | CCL4–CCR5 | 4 | +0.006 | 0.875 | 0.981 |
| T_recruit | CXCL16–CXCR6 | 7 | +0.050 | 0.016 | 0.145 |

Same pattern as TACSTD2: MHC-I higher from the high state (7/7); CXCR3 sparse
(n=3); CXCL16–CXCR6 higher from CLDN4-high. **0/12 FDR < 0.05.**

### Incoming T/NK → malignant

| Pathway | Pair | n patients | median Δ | Wilcoxon p | FDR |
|---|---|---:|---:|---:|---:|
| IFN | IFNG–IFNGR2 | 7 | +0.011 | 0.578 | 0.850 |
| IFN | IFNG–IFNGR1 | 7 | +0.018 | 0.297 | 0.558 |
| IFN | IFNG–IFNGR1+IFNGR2 | 7 | +0.018 | 0.219 | 0.558 |

Incoming IFNG–IFNGR is again flat.

## LIANA CellPhoneDB method (secondary, pooled / downsampled)

Focus edges that cleared LIANA `expr_prop=0.10` from malignant-like senders to T or NK:

| Split | Pair | high `lr_means` | low `lr_means` | higher in |
|---|---|---:|---:|---|
| TACSTD2 | CXCL16–CXCR6 → T | 0.440 | 0.393 | high |
| TACSTD2 | CXCL16–CXCR6 → NK | 0.409 | 0.361 | high |
| TACSTD2 | HLA-E–KLRC1 → NK | 0.946 | 0.852 | high (p=1, not specific) |
| TACSTD2 | CX3CL1–CX3CR1 → NK | — | 0.414 | low only |
| CLDN4 | CXCL16–CXCR6 → T | 0.416 | 0.411 | high (near-tie) |
| CLDN4 | HLA-E–KLRC1 → NK | 0.944 | 0.836 | high (p=1, not specific) |
| CLDN4 | CX3CL1–CX3CR1 → NK | 0.401 | 0.396 | high (near-tie) |

CXCL9/10/11–CXCR3 did **not** clear LIANA’s 10% expression filter in epithelium.
LIANA p-values are within-object specificity, not patient-level tests.

## Readout on a TROP2-high reduced T-recruit / IFN / MHC-I pattern

In these 8 paired patients, TACSTD2-high malignant-like cells do **not** show a
coordinated reduction on those three axes:

- **MHC-I outgoing** is higher from TACSTD2-high (7/7 pairs; median Δ +0.12; min FDR 0.098).
- **T-recruit outgoing** is mixed and dropout-limited. The only well-powered chemokine
  edge (CXCL16–CXCR6, n=7) is stronger from TACSTD2-high. CXCR3 ligands are
  detected in 3 patients and are not lower after FDR.
- **IFN** is an incoming T/NK → malignant edge (IFNG–IFNGR) and is flat.

CLDN4-high tracks the same direction. This is the rank in this public
15-patient BD Rhapsody object.

## Limits

- BD Rhapsody, not 10x. GEO deposited no per-cell labels; A3 marker gates are used as given.
- Malignant-like is a normal-lung-marker exclusion, not public CopyKAT calls.
- Ambient RNA cannot be re-estimated from the processed matrix.
- Chemokine dropout is high; read `n patients` and `pass_expr_prop` with the ranks.
- Paired n = 8. MPR malignant-like compartments are nearly empty (P06/P11/P14).
- CellChat was not run.

## Files

| File | Role |
|---|---|
| `n_cells_patients.tsv` | Per-sample cell counts and high/low bins |
| `ranks_*_{outgoing,incoming}.tsv` | Patient-level pair ranks |
| `pathway_summary.tsv` | T-recruit / IFN / MHC-I axis rollup |
| `patient_*.tsv` | Per-patient pair scores (full CellPhoneDB list) |
| `pooled_*.tsv` | All-cell descriptive scores |
| `liana_cellphonedb_*.csv` | Full LIANA CellPhoneDB output |
| `liana_*_{outgoing,incoming}_tnk.csv` | LIANA edges involving T/NK |
| `figures/` | n-cell bars and pathway Δ plots |
| `summary.json` | Machine-readable n and method flags |
