# FINDING — CellChat-style outgoing CLDN4-high malignant on the CXCL13+ trio

**Verdict (patient-level, CLDN4-only):** on the given CXCL13+ trio, outgoing
**NECTIN2–TIGIT**, **LGALS9–CD45**, **HLA-E–CD8A**, **JAM1–ITGAL/ITGB2**, and
**MDK–NCL** are higher from CLDN4-high epithelium toward CXCL13+ T and toward
T/NK. **CXCL16–CXCR6 is also higher from high** (recruit-up; not a chemokine-low
story). **CD274–PDCD1 is not differential** (median ΔP = 0; thin detection).
Communication n is **50** (Mal→T/NK) / **34** (Mal→CXCL13+ T), not the given
Spearman n=60.

**CLDN4 only. No dual-high. Not GSE207422-only.** Patient is the unit.

The CXCL13+ trio that already differs is taken as given and is **not re-audited**:
GSE148071 + GSE207422 + GSE253013, **n=60, ρ=−0.425, p=0.00121**
(`methods/scrna_cldn4_combo`, PR #290, family `tls/cxcl13pos/mean`).

This folder adds **outgoing CellChat-style probability** from **CLDN4-high vs CLDN4-low**
putative malignant epithelium toward **CXCL13+ T** and toward **T/NK**, scored
**within each patient** (median CLDN4 gate), then tested across patients
(Wilcoxon on ΔP). CellChat R and LIANA were not run.

---

## Honest n (communication, not the given Spearman n=60)

Given combo n=60 is mean-CLDN4 vs CXCL13+ T *fraction*. Communication n is the
number of patients with enough cells on **both** CLDN4 arms **and** the receiver.

| Cohort | Deposited (this extract) | Eligible Mal→T/NK | Eligible Mal→CXCL13+ T |
| --- | ---: | ---: | ---: |
| GSE148071 | 42 | 27 | 12 |
| GSE207422 | 15 | 15 | 14 |
| GSE253013 | 9 | 8 | 8 |
| **Trio** | **66** | **50** | **34** |

Floors: ≥10 CLDN4-high and ≥10 CLDN4-low epithelial cells
(within-patient median of `log1p(CP10k)` CLDN4); ≥15 T/NK; ≥8 CXCL13+ T
(T lineage and CXCL13 UMI > 0). Missing this run: **none**.
Cohorts scored: GSE148071, GSE207422, GSE253013.

Do **not** read the LR table as n=60. The given 60 is the Spearman set
(36+15+9). This extract lists all 42+15+9 GEO patients; 15 GSE148071 biopsies
fail the communication floor (thin epithelium and/or T/NK). GSE253013 **MRC003**
has 19 author tumor epithelial cells (9/10 arms) and 2 CXCL13+ T cells and is
not scored. Patients below the floor stay in `results/patient_n.tsv`.

## Gate

- Marker: **CLDN4 only**. TACSTD2 is recorded in the stream and is not used.
- No dual-high TACSTD2×CLDN4 quadrant.
- Split: **within-patient median** among epithelial cells (putative malignant;
  marker-argmax; GSE253013 uses author epithelial/malignant when present).
- Probability: Jin et al. 2021 CellChat Hill / mass-action on 10% truncated
  means, CellChatDB v2 protein pairs, `expr_prop ≥ 0.10`, Kh=0.5.
- Test: Wilcoxon signed-rank on per-patient ΔP = P(high→receiver) − P(low→receiver),
  among patients where the pair is detected on at least one arm (n≥6).
  BH-FDR within receiver. p-values are descriptive.

## Ligand table — outgoing CLDN4-high → CXCL13+ T

Full table: `results/ligand_table.tsv`.

**Higher from CLDN4-high (median ΔP > 0):**

| interaction_name | pathway_name | ligand_class | n_wilcoxon | n_cohorts | median_delta | wilcoxon_p | n_delta_pos |
| --- | --- | --- | --- | --- | --- | --- | --- |
| NECTIN2_TIGIT | NECTIN | barrier|inhibitory | 31 | 3 | +0.056 | 8.42e-07 | 28 |
| CD55_ADGRE5 | ADGRE | other | 32 | 3 | +0.022 | 1.29e-06 | 28 |
| CXCL16_CXCR6 | CXCL | recruit | 23 | 3 | +0.007 | 6.03e-05 | 20 |
| GDF15_TGFBR2 | GDF | other | 18 | 3 | +0.001 | 7.63e-05 | 17 |
| NECTIN2_CD226 | NECTIN | barrier|inhibitory | 22 | 3 | +0.002 | 9.87e-05 | 20 |
| HLA-E_CD8B | MHC-I | inhibitory | 24 | 3 | +0.008 | 0.000205 | 20 |
| ITGAV_ITGB1_ADGRE5 | ADGRE | other | 28 | 3 | +0.013 | 0.00187 | 18 |
| LGALS9_CD45 | GALECTIN | inhibitory | 28 | 3 | +0.041 | 0.00204 | 21 |
| HLA-E_CD8A | MHC-I | inhibitory | 29 | 3 | +0.007 | 0.00299 | 22 |
| CD99_CD99 | CD99 | other | 34 | 3 | +0.025 | 0.00339 | 25 |
| ICAM1_ITGAL | ICAM | other | 28 | 3 | +0.008 | 0.00344 | 20 |
| JAM1_ITGAL_ITGB2 | JAM | barrier | 28 | 3 | +0.015 | 0.00478 | 22 |

**Higher from CLDN4-low (median ΔP < 0):**

No pair with n_wilcoxon≥6 has a meaningfully negative median ΔP. PPIA–BSG is
−0.000 (p=0.78). CD274–PDCD1 is detected in only 8 / 6 patients (CXCL13+ T /
T/NK) with median ΔP = 0.

## Ligand table — outgoing CLDN4-high → T/NK

**Higher from CLDN4-high (median ΔP > 0):**

| interaction_name | pathway_name | ligand_class | n_wilcoxon | n_cohorts | median_delta | wilcoxon_p | n_delta_pos |
| --- | --- | --- | --- | --- | --- | --- | --- |
| NECTIN2_TIGIT | NECTIN | barrier|inhibitory | 42 | 3 | +0.022 | 1.27e-08 | 37 |
| CD55_ADGRE5 | ADGRE | other | 37 | 3 | +0.007 | 1.61e-06 | 31 |
| MDK_NCL | MK | other | 50 | 3 | +0.056 | 2.43e-05 | 36 |
| JAM1_ITGAL_ITGB2 | JAM | barrier | 33 | 3 | +0.005 | 0.000145 | 29 |
| ICAM1_SPN | ICAM | other | 39 | 3 | +0.004 | 0.000192 | 31 |
| HLA-E_CD8A | MHC-I | inhibitory | 39 | 3 | +0.005 | 0.000388 | 29 |
| LGALS9_CD45 | GALECTIN | inhibitory | 38 | 3 | +0.018 | 0.000442 | 29 |
| MDK_ITGA4_ITGB1 | MK | other | 38 | 3 | +0.003 | 0.000643 | 26 |
| LGALS9_CD44 | GALECTIN | inhibitory | 38 | 3 | +0.006 | 0.000726 | 29 |
| HLA-C_CD8A | MHC-I | other | 39 | 3 | +0.001 | 0.000943 | 30 |
| APP_CD74 | APP | other | 50 | 3 | +0.015 | 0.00432 | 34 |
| CD99_CD99 | CD99 | other | 49 | 3 | +0.020 | 0.00466 | 35 |

**Higher from CLDN4-low (median ΔP < 0):**

_no pair with n_wilcoxon≥6 and median ΔP<0_ (same one-sided pattern as CXCL13+ T).

## What is not supported

- Re-using **n=60** as the communication n.
- A GSE207422-only CellChat (already in `methods/scrna_cellchat_cldn4`).
- Dual-high TACSTD2×CLDN4 gates.
- Cell-pooled permutation as a patient claim.
- ICI response or histology as the unit of this table.

## Files

`results/ligand_table.tsv` (pooled patient-level LR),
`results/lr_pairs_patient.tsv`, `results/patient_n.tsv`,
`results/fig_n_per_patient.png`, `results/fig_n_eligible.png`,
`results/fig_extra_ligand_table.png`, `results/fig_top_outgoing_CXCL13T.png`,
`results/fig_top_outgoing_TNK.png`, `results/summary.json`.
