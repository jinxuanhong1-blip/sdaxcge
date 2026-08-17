# FINDING — Harmony n=73 CLDN4-only high-end (CellChat-style + LIANA-style)

**ADDITIVE. CLDN4 only. No dual-high.**

The Harmony / multi-cohort patient table is taken as given and was **not** re-audited:
malignant-like *CLDN4* vs T/NK **n=73, ρ=−0.27** (GSE131907 + GSE253013 + GSE148071 + GSE127465;
`data/given/per_donor_metrics.tsv`, `association_statistics.tsv`).

On the Harmony n=73 table (given ρ=−0.27, not re-audited), outgoing LR was scored in **32** paired patients from GSE131907+GSE148071+GSE127465 (GSE253013 n=9 skipped: 9 GB RDS not downloaded). Focus pairs: 0/11 have median Δ<0; **4/11 reach FDR<0.05**. Recruit median Δ=0.015; barrier/inhibitory median Δ=0.041.

## What was run

| Method | Status |
|---|---|
| Given Harmony malignant CLDN4 vs T/NK (n=73, ρ=−0.27) | **taken as given — not re-audited** |
| CellChat-style Hill P (10% truncated means, K_h=0.5, CellChatDB v2 protein) | **run** — patient-level Δ, then meta |
| LIANA / CellPhoneDB-style score (mean of min-subunit means on log1p) | **primary** — patient-level Wilcoxon + cohort DL |
| LIANA R/Python `liana.mt.cellphonedb` | **not run** — no joint object; documented score used instead |
| CellChat R | **not run** — R package not used |
| Milo-style nhood on the Harmony joint embedding | **not run** — embedding is not stored; GSE253013 RDS not downloaded so it cannot be rebuilt |

## Honest n

| Set | n |
|---|---:|
| Given Harmony patients (malignant CLDN4 vs T/NK) | **73** |
| GSE253013 in that table (no public processed UMI used here) | 9 |
| Patients with a public processed matrix (131907+148071+127465) | 64 |
| Patients scored for LR (seen in a matrix) | 64 |
| **Paired LR patients** (≥8 CLDN4-high **and** ≥8 CLDN4-low malignant-like, ≥20 T/NK) | **32** |
| Paired GSE131907 / GSE148071 / GSE127465 | 14 / 17 / 1 |
| Dropped (one-sided / thin malignant or T/NK) | 32 |
| CellChatDB v2 protein pairs scored | 2239 |
| CellPhoneDB v5 pairs scored | 2920 |

The header n for the Wilcoxon is the paired n, **not 73**. GSE253013's 9 patients stay in the given ρ=−0.27 row and are absent from the LR table.

CLDN4-high / low is the **patient-specific median** of malignant-like *CLDN4* (same Harmony malignant-like rule: epithelial marker-argmax and near-zero normal-lung score). TACSTD2 is not a gate.

### Paired patients

| Dataset | Donor | mal high | mal low | T/NK |
|---|---|---:|---:|---:|
| GSE131907 | BRONCHO_11 | 57 | 56 | 789 |
| GSE131907 | EBUS_06 | 306 | 306 | 581 |
| GSE131907 | EBUS_10 | 83 | 82 | 1487 |
| GSE131907 | EBUS_12 | 60 | 60 | 335 |
| GSE131907 | EBUS_28 | 1400 | 1399 | 229 |
| GSE131907 | EFFUSION_11 | 14 | 14 | 1927 |
| GSE131907 | LUNG_T18 | 111 | 110 | 1054 |
| GSE131907 | LUNG_T28 | 29 | 29 | 1535 |
| GSE131907 | NS_02 | 66 | 65 | 986 |
| GSE131907 | NS_03 | 1073 | 1072 | 138 |
| GSE131907 | NS_04 | 171 | 171 | 123 |
| GSE131907 | NS_06 | 39 | 38 | 369 |
| GSE131907 | NS_07 | 1961 | 1961 | 140 |
| GSE131907 | NS_19 | 202 | 201 | 354 |
| GSE148071 | P1 | 1722 | 1722 | 132 |
| GSE148071 | P10 | 1180 | 1180 | 54 |
| GSE148071 | P13 | 27 | 27 | 99 |
| GSE148071 | P18 | 744 | 744 | 62 |
| GSE148071 | P19 | 25 | 25 | 45 |
| GSE148071 | P22 | 71 | 71 | 198 |
| GSE148071 | P26 | 33 | 32 | 76 |
| GSE148071 | P27 | 13 | 12 | 81 |
| GSE148071 | P28 | 188 | 188 | 67 |
| GSE148071 | P3 | 3437 | 3437 | 20 |
| GSE148071 | P32 | 90 | 89 | 70 |
| GSE148071 | P38 | 32 | 32 | 147 |
| GSE148071 | P40 | 25 | 24 | 346 |
| GSE148071 | P5 | 54 | 54 | 51 |
| GSE148071 | P6 | 634 | 633 | 67 |
| GSE148071 | P8 | 23 | 23 | 60 |
| GSE148071 | P9 | 934 | 933 | 61 |
| GSE127465 | p3 | 247 | 247 | 2001 |

Per-patient counts: `results/n_patients.tsv`.

## Score

- **LIANA-style (primary):** on log1p(CP10k) for GSE131907 / GSE148071, and log1p of the deposited normalized counts for GSE127465. Partner expression = minimum subunit mean. Pair score = mean of the two partner means (Efremova 2020; Garcia-Alonso 2022). `pass_expr_prop` = both partners in ≥10% of cells.
- **CellChat-style (companion):** 10% truncated means, geometric mean of subunits, Hill P = LR / (0.5 + LR). Not a CellChat R communication probability.
- Unit = **patient**. Receiver T/NK are that patient's own T/NK (so outgoing Δ is ligand-driven). Cells are not replicates.
- Meta: Wilcoxon signed-rank on patient Δ (high−low); BH-FDR within method. Cohort DerSimonian–Laird on cohort mean Δ when ≥2 cohorts have ≥4 patients.

## Primary ligand table — outgoing CLDN4-high malignant → T/NK (LIANA-style)

Median patient Δ = high − low. Negative = weaker from the CLDN4-high state.

| Class | Pair | n patients | n cohorts | median Δ | Wilcoxon p | FDR | n+ / n− |
|---|---|---:|---:|---:|---:|---:|---:|
| MHC_I | HLA-C–CD8A | 24 | 3 | +0.057 | 0.00958 | 0.060 | 19/5 |
| MHC_I | HLA-A–CD8A | 24 | 3 | +0.064 | 0.00252 | 0.035 | 18/6 |
| MHC_I | HLA-B–CD8A | 24 | 3 | +0.089 | 0.00434 | 0.036 | 20/4 |
| barrier | PVR–TIGIT | 11 | 2 | +0.012 | 0.320 | 0.480 | 8/3 |
| barrier | NECTIN2–TIGIT | 12 | 2 | +0.015 | 0.077 | 0.157 | 8/4 |
| barrier | F11R–ITGAL+ITGB2 | 17 | 3 | +0.073 | 0.031 | 0.103 | 14/3 |
| barrier | CDH1–ITGAE+ITGB7 | 13 | 2 | +0.082 | 0.000244 | 0.00421 | 13/0 |
| inhibitory | LGALS9–HAVCR2 | 6 | 1 | +0.011 | 0.562 | 0.723 | 3/3 |
| inhibitory | HLA-E–KLRD1 | 22 | 3 | +0.068 | 0.085 | 0.164 | 16/6 |
| recruit | CCL5–CCR5 | 6 | 2 | +0.007 | 0.844 | 0.872 | 4/2 |
| recruit | CXCL16–CXCR6 | 16 | 2 | +0.022 | 0.00516 | 0.036 | 13/3 |

0/11 focus pairs have median Δ < 0. **4/11 reach FDR < 0.05** (CDH1–ITGAE/ITGB7; HLA-A/B–CD8A; CXCL16–CXCR6).

**Not detected** at `expr_prop ≥ 0.10` in ≥4 paired patients: CD274–PDCD1, CXCL9–CXCR3, CXCL10–CXCR3, CXCL11–CXCR3, CX3CL1–CX3CR1, CEACAM1–CD8A, LGALS9–PTPRC (LIANA-style; LGALS9–PTPRC is present in the CellChat-style companion).

The given abundance association (n=73, ρ=−0.27) is **not** accompanied by weaker outgoing T-recruit / MHC-I scores from CLDN4-high malignant cells in the 32 patients who can be paired. CXCL16–CXCR6 is higher from the high arm (n=16, FDR=0.036), same direction as the GSE207422 / GSE148071 CellChat extras.

32 of 64 matrix patients were dropped because the patient-specific *CLDN4* median was 0 (high dropout → empty low arm) or a bin had <8 malignant-like cells. That is reported, not patched.

Full ranked table: `results/ligand_table_liana_outgoing.tsv` (primary) and `results/ligand_table_cellchat_outgoing.tsv`.
Patient-level deltas (`pass_expr_prop` on either arm): `results/patient_deltas.tsv`.

## CellChat-style companion (same patients)

| Class | Pair | n | median ΔP | Wilcoxon p | FDR |
|---|---|---:|---:|---:|---:|
| MHC_I | HLA-A–CD8A | 24 | +0.016 | 0.016 | 0.075 |
| MHC_I | HLA-C–CD8A | 24 | +0.021 | 0.014 | 0.075 |
| MHC_I | HLA-B–CD8A | 24 | +0.022 | 0.00652 | 0.064 |
| barrier | PVR–TIGIT | 11 | +0.004 | 0.206 | 0.334 |
| barrier | NECTIN2–TIGIT | 12 | +0.008 | 0.110 | 0.223 |
| inhibitory | LGALS9–HAVCR2 | 6 | +0.003 | 0.062 | 0.143 |
| inhibitory | HLA-E–CD8A | 24 | +0.029 | 0.00792 | 0.064 |
| inhibitory | LGALS9–PTPRC | 22 | +0.036 | 0.021 | 0.080 |
| recruit | CCL5–CCR5 | 6 | +0.001 | 0.688 | 0.793 |
| recruit | CXCL16–CXCR6 | 16 | +0.004 | 0.00214 | 0.039 |

## Milo / joint embedding

The Harmony joint UMAP/PCA is **not** in the repository (only donor scores and PNG figures).
Rebuilding it requires the GSE253013 Garnett RDS (~9 GB), which was not downloaded.
No Milo-style neighborhood table is claimed. File: `results/nhood_skip.json`.

## What is not claimed

- The given n=73 ρ=−0.27 was not recomputed.
- Dual-high (TACSTD2 **and** CLDN4) was not run.
- Cell-level p-values are not reported.
- GSE253013 LR is absent (honest skip).
- A rebuilt Harmony embedding / Milo DA test.
- ICI response (these four series are treatment-naïve / diagnostic atlases).

## Files

| File | Role |
|---|---|
| `FINDING.md` | This note |
| `results/ligand_table_liana_outgoing.tsv` | Primary ligand table (honest n) |
| `results/ligand_table_cellchat_outgoing.tsv` | CellChat-style companion |
| `results/ligand_table_focus.tsv` | Curated outgoing pairs |
| `results/n_patients.tsv` | Per-patient cell counts |
| `results/patient_deltas.tsv` | Per-patient pair Δ |
| `results/nhood_skip.json` | Why Milo was not run |
| `results/summary.json` | Machine-readable n |
| `results/figures/` | Extra figures |

## Reproduce

```bash
cd methods/harmony73_cldn4_hiend
python3 scripts/00_download.py   # does not fetch GSE253013
python3 scripts/01_run_lr.py
```

