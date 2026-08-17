# ICI trio CLDN4-only: GSE207422 + GSE205335 + GSE179994

Additive **CLDN4-only** (no dual-high, no TACSTD2 gate). Patient is the unit.

Intended merge: **GSE207422** (Hu *Genome Med* 2023, neoadjuvant PD-1 + chemo) +
**GSE205335** (palliative ICI biopsies) + **GSE179994** (pembro + carboplatin +
pemetrexed) if GSE179994 has a public processed all-cell matrix <2 GB.

## GSE179994 is unusable — trio runs as 207422 + 205335

The only public processed expression file under 2 GB is
`GSE179994_all.Tcell.rawCounts.rds.gz` (**421 MB**, under the cap) but it is
**T cells only**. Author metadata (`GSE179994_Tcell.metadata.tsv.gz`): 150,849
barcodes, 36 patients, 47 samples; `celltype` is CD8 / CD4 / NA. No epithelial,
malignant, or EPCAM compartment. No all-cell UMI, MTX, h5, or TISCH extract is
on GEO. Series text: raw data not provided.

Without a malignant compartment, CLDN4 mean/%pos and a same-patient T/NK
**fraction** are undefined (a T-cell-only matrix would give T/NK ≈ 1 by
construction).

| Quantity | Value |
|---|---:|
| Analysis n (malignant CLDN4 + all-cell T/NK denominator) | **0** |
| T-cell barcodes in metadata | 150,849 (not n) |
| T-cell patients in metadata | 36 (not n) |

The T-cell RDS was **not** downloaded. 207422 is included. The merge below is
the **ICI duo GSE207422 + GSE205335**.

## Honest n — patient-level malignant CLDN4 vs T/NK

Gates: ≥10 malignant cells and ≥20 T/NK. Score = malignant CLDN4 **mean
log1p(CP10k)** vs same-patient T/NK fraction. Cells are not replicates.

| Cohort | Malignant def | Immune def | Eligible n | Used in primary combo | Note |
|---|---|---|---:|---:|---|
| GSE207422 | A3 marker malignant-like (epithelial AND zero UMI of SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3) | marker T or NK | 10 (all) / **7** (post-ICI) | 7 | Author per-cell labels are not on GEO. MPR P06/P11/P14 have 1 / 0 / 0 malignant-like cells and drop. TN P01/P05/P08 are eligible but excluded from the ICI primary row. |
| GSE205335 | author `Malignant cells` | author `T/NK cells` | **22** | 22 | Tumor tissues only (Normal LN/lung/brain dropped). RECIST only; no MPR. NSCLC ADC+SQ = 17; SCLC+NUT = 5. |
| GSE179994 | none | none | **0** | 0 | T-cell-only public matrix. |
| **Primary combo** | mixed (as above) | mixed | **29** | **29** | k=2. 179994 not in the pool. |

Q4 vs Q1 is not done (would need n≥16 per cohort or a pre-specified pooled
quartile; not claimed).

## Combo table — malignant CLDN4 vs T/NK

DerSimonian–Laird random-effects on Fisher-z of Spearman ρ. p-values are
descriptive.

| Subset | k | N | ρ | p | I² | Note |
|---|---:|---:|---:|---:|---:|---|
| GSE207422 post-ICI | 1 | 7 | +0.286 | 0.535 | — | A3 malignant-like; small n |
| GSE207422 all eligible | 1 | 10 | −0.018 | 0.960 | — | includes 3 TN biopsies |
| GSE205335 all | 1 | 22 | −0.200 | 0.371 | — | author malignant |
| GSE205335 NSCLC | 1 | 17 | −0.022 | 0.933 | — | drop SCLC/NUT |
| GSE179994 | 0 | 0 | NA | NA | — | unusable |
| **combo_207422post_205335all** | **2** | **29** | **−0.116** | **0.576** | **0%** | **primary ICI duo; 207422 included** |
| combo_207422post_205335NSCLC | 2 | 24 | +0.048 | 0.838 | 0% | drop SCLC/NUT |
| combo_207422all_205335all | 2 | 32 | −0.152 | 0.434 | 0% | 207422 includes TN |

Primary row: **N=29, ρ=−0.116 (95% CI −0.48 to +0.28), p=0.576, I²=0%**.
GSE207422 post-ICI is a **positive** point estimate on 7 patients; GSE205335 is
weakly negative on 22. The duo does **not** support a patient-level inverse
CLDN4–T/NK claim.

Secondary score (malignant CLDN4 **%pos** vs T/NK), same patients: GSE205335
ρ=−0.435 p=0.043 n=22; duo DL RE ρ=−0.357 p=0.074 N=29. That is a different
score (detection, not mean) and is not the primary row.

Full numbers: `tables/combo_cldn4_tnk.tsv`. Patient rows:
`tables/patient_cldn4_tnk.tsv`.

## CellChat-style outgoing CLDN4-high malignant → T/NK

Not CellChat (R / CellChat were not run). Documented CellPhoneDB-style score:
partner expression = min(subunit means) on log1p(CP10k); pair score = mean of
the two partner means (Efremova 2020; Garcia-Alonso 2022). High/low =
**within-patient median** of malignant CLDN4. A patient enters if both bins
have ≥10 malignant cells and ≥20 T/NK. Negative Δ = weaker outgoing score from
the CLDN4-high state. Receiver T/NK is the same within a patient, so Δ is the
ligand difference on the malignant side.

Paired n: **GSE207422 = 6** (P02 has only 10 malignant-like cells and cannot
split), **GSE205335 = 22**, **combo = 28**. GSE179994 paired n = 0.

### LR table — pairs with `expr_prop` ≥10% in ≥3 patients (combo)

Same-ligand MHC pairs share Δ (receiver T/NK is fixed). One canonical receptor
is shown; KIR partners of the same HLA are omitted (they do not pass
`expr_prop`). IFN ligands are T/NK products and are **not** treated as
malignant outgoing edges (IFNG is detected in ≥10% of malignant cells in only
3/28 patients).

| Pathway | Pair | n | n pass `expr_prop` | median Δ | n Δ<0 | Wilcoxon p | FDR |
|---|---|---:|---:|---:|---:|---:|---:|
| T_recruit | CCL5–CCR5 | 28 | 6 | −0.028 | 20 | 0.0014 | 0.006 |
| T_recruit | CCL4–CCR5 | 28 | 6 | −0.003 | 19 | 0.032 | 0.043 |
| T_recruit | CCL3–CCR5 | 28 | 3 | −0.003 | 17 | 0.041 | 0.053 |
| T_recruit | CXCL9–CXCR3 | 28 | 3 | 0.000 | 13 | 0.761 | 0.843 |
| T_recruit | CXCL10–CXCR3 | 28 | 6 | +0.000 | 10 | 0.475 | 0.545 |
| T_recruit | CXCL11–CXCR3 | 28 | 3 | +0.001 | 5 | 0.011 | 0.016 |
| T_recruit | CX3CL1–CX3CR1 | 28 | 3 | +0.025 | 3 | 2.6e-05 | 4.1e-04 |
| T_recruit | CXCL16–CXCR6 | 28 | 16 | +0.064 | 2 | 1.4e-07 | 4.4e-06 |
| MHC_I | HLA-B–CD8A | 28 | 25 | +0.083 | 7 | 0.010 | 0.016 |
| MHC_I | HLA-A–CD8A | 28 | 25 | +0.101 | 7 | 0.010 | 0.016 |
| MHC_I | HLA-C–CD8A | 28 | 25 | +0.111 | 7 | 0.008 | 0.016 |
| MHC_I | HLA-E–KLRD1 | 28 | 25 | +0.130 | 4 | 0.004 | 0.013 |

Readout: MHC-I classical/non-classical ligands are **higher** from CLDN4-high
malignant cells (7/7 displayed pairs; 21–24/28 patients Δ>0). The only
well-powered chemokine edge is **CXCL16–CXCR6** (16/28 pass `expr_prop`), also
**higher** from CLDN4-high. CXCR3 ligands are sparse. CCL5–CCR5 is the one
T-recruit pair with median Δ<0 and FDR<0.05, but it passes `expr_prop` in only
6/28 patients.

This is **not** a coordinated reduction of T-recruit / MHC-I outgoing from
CLDN4-high cells on this duo.

Full per-cohort and combo pairs (including dropout-only edges):
`tables/lr_outgoing_cldn4high_tnk.tsv`. Patient-level scores:
`tables/lr_patient_scores.tsv`.

## What was not done

- No dual-high / TACSTD2 gate.
- No CellChat communication probability.
- GSE179994 T-cell RDS not downloaded (filename + author metadata already
  establish T-cell-only content).
- Author CopyKAT / DRMref labels for GSE207422 are not public; A3 marker
  malignant-like is used. That empties several MPR samples and is why 207422
  combo n is 7, not the 12-patient author-label table used elsewhere.
- GSE205335 has no MPR field; RECIST is not substituted for MPR.

## Files

| File | Role |
|---|---|
| `tables/combo_cldn4_tnk.tsv` | Combo + single-cohort Spearman / DL RE |
| `tables/patient_cldn4_tnk.tsv` | Patient-level malignant CLDN4 and T/NK |
| `tables/lr_outgoing_cldn4high_tnk.tsv` | LR table (per-cohort and combo) |
| `tables/lr_patient_scores.tsv` | Per-patient outgoing scores |
| `tables/gse179994_feasibility.json` | Why 179994 is out |
| `figures/scatter_cldn4_tnk.png` | Patient scatter |
| `figures/lr_outgoing_combo.png` | Combined LR Δ |

## Reproduce

```bash
python3 -m pip install -r methods/ici_trio_207422_205335_179994_cldn4/requirements.txt
python3 methods/ici_trio_207422_205335_179994_cldn4/download.py
python3 methods/ici_trio_207422_205335_179994_cldn4/analyze.py
```
