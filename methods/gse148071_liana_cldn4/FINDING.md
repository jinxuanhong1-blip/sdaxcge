# FINDING — GSE148071 LIANA/LR from CLDN4-high malignant to T/NK

**Object:** TISCH2 `NSCLC_GSE148071` (Wu et al. 2021 *Nat Commun*, PMID 33953163; GEO GSE148071).
**Sender:** TISCH2 `Malignant` cells with CLDN4 ≥ global malignant median (1.121 on TISCH2 log2(TPM/10+1)).
**Receiver:** TISCH2 T/NK major-lineage. In this object that is **CD8T + Tprolif only** — **0 NK / CD4T / Treg / NKT cells** are labeled.
**Unit of inference:** patient (n=42; one Sample each). Cells are not replicates.

## Honest n

| Item | n | Note |
| --- | ---: | --- |
| Patients in TISCH2 object | **42** | Wu et al. advanced NSCLC biopsies |
| Cells aligned to h5 | 82,267 | all Tissue=Tumor |
| Malignant | 48,118 | TISCH2 major-lineage |
| T/NK total | **4,384** | CD8T=3,750; Tprolif=634; NK=0 |
| CLDN4-high / low malignant | 24,060 / 24,058 | median split among malignant |
| Patients with ≥10 CLDN4-high malignant **and** ≥20 T/NK | **20** | used for per-patient scores |
| Patients also having ≥10 CLDN4-low malignant | **17** | used for high−low Wilcoxon |
| Patients dropped (too few T/NK or senders) | 22 | see `n_cells_patients.tsv` |
| CPDB pairs scored / passing 10% expr_prop | 1,885 / **25** | pooled CLDN4-high → T/NK |
| LIANA | LIANA_OK | LIANA_OK n_edges=125 downsampled={'TNK': 2000, 'Malig_CLDN4high': 2000, 'Malig_CLDN4low': 2000} file=liana_cellphonedb.csv |
| CellChat | not run | R unavailable |

Per-patient counts: `results/n_cells_patients.tsv`. Do not treat cell counts as the sample size.

Usable outgoing patients (n=20): P1, P10, P12, P13, P14, P18, P19, P22, P23, P24, P26, P28, P32, P38, P4, P40, P6, P7, P8, P9.
Paired high vs low (n=17): P1, P10, P13, P14, P18, P19, P22, P23, P24, P28, P38, P4, P40, P6, P7, P8, P9.

## What was run

Primary score is the documented CellPhoneDB mean (Efremova 2020; Garcia-Alonso 2022) on the public TISCH2/MAESTRO matrix **log2(TPM/10+1)**: partner expression = min(subunit means); pair score = mean of the two partner means. A pair passes `expr_prop` when both partners are detected (>0) in ≥10% of cells in their group. This is **not** a CellChat probability and **not** a raw-UMI reprocess.

LIANA `mt.cellphonedb` status: `LIANA_OK`. LIANA p-values (if present) are within-object specificity after downsampling, not patient-level tests.

## LR table — pooled CLDN4-high malignant → T/NK

Top 25 pairs by CPDB mean score among those passing 10% expression in **both** partners. Sender n=24,060 cells; receiver n=4,384 cells. Pooled ranks are descriptive.

| pathway | ligand | receptor | ligand_frac | receptor_frac | cpdb_mean_score |
|---|---|---|---:|---:|---:|
| other | APP | CD74 | 0.306 | 0.665 | 1.074 |
| MHC_I | HLA-A | CD8A | 0.780 | 0.165 | 0.976 |
| other | CD58 | CD2 | 0.187 | 0.597 | 0.923 |
| MHC_I | HLA-B | CD8A | 0.729 | 0.165 | 0.915 |
| other | PPIA | BSG | 0.784 | 0.164 | 0.853 |
| T_recruit_partial | CXCL14 | CXCR4 | 0.389 | 0.236 | 0.804 |
| MHC_I | HLA-C | CD8A | 0.689 | 0.165 | 0.781 |
| other | SPP1 | ITGA4+ITGB1 | 0.407 | 0.160 | 0.625 |
| other | CD44 | TYROBP | 0.476 | 0.256 | 0.618 |
| other | CD55 | ADGRE5 | 0.375 | 0.138 | 0.422 |
| other | F11R | ITGAL+ITGB2 | 0.409 | 0.136 | 0.364 |
| MHC_I_partial | CEACAM5 | CD8A | 0.165 | 0.165 | 0.363 |
| other | IGFBP3 | TMEM219 | 0.250 | 0.147 | 0.347 |
| checkpoint | NECTIN2 | TIGIT | 0.175 | 0.223 | 0.324 |
| other | DHCR24 | RORA | 0.272 | 0.161 | 0.306 |
| other | JAG1 | CD46 | 0.238 | 0.159 | 0.276 |
| other | CD47 | SIRPG | 0.302 | 0.108 | 0.273 |
| other | ICAM1 | SPN | 0.118 | 0.148 | 0.213 |
| other | DHCR7 | RORA | 0.122 | 0.161 | 0.211 |
| other | PLAUR | ITGA4+ITGB1 | 0.106 | 0.160 | 0.205 |
| other | LIPA | RORA | 0.107 | 0.161 | 0.203 |
| other | ICAM1 | ITGAL | 0.118 | 0.136 | 0.200 |
| other | ICAM1 | ITGAL+ITGB2 | 0.118 | 0.136 | 0.200 |
| other | ITGAV+ITGB1 | ADGRE5 | 0.142 | 0.138 | 0.199 |
| other | LPAR1 | ADGRE5 | 0.101 | 0.138 | 0.175 |

### Focus axes that cleared 10% (T-recruit / checkpoint / MHC-I)

| pathway | ligand | receptor | ligand_frac | receptor_frac | cpdb_mean_score |
|---|---|---|---:|---:|---:|
| MHC_I | HLA-A | CD8A | 0.780 | 0.165 | 0.976 |
| MHC_I | HLA-B | CD8A | 0.729 | 0.165 | 0.915 |
| MHC_I | HLA-C | CD8A | 0.689 | 0.165 | 0.781 |
| checkpoint | NECTIN2 | TIGIT | 0.175 | 0.223 | 0.324 |

## Patient-level high vs low (honest paired n)

Paired Wilcoxon across **17** patients with ≥10 CLDN4-high, ≥10 CLDN4-low malignant, and ≥20 T/NK. Δ = high − low CPDB score. FDR is BH within this contrast. Negative Δ = weaker from the CLDN4-high state.

| pathway | ligand | receptor | n_patients | median_delta | pval | padj |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| MHC_I | HLA-B | CD8B | 6 | -0.045 | 0.844 | 1 |
| MHC_I | HLA-C | CD8B | 6 | -0.022 | 0.844 | 1 |
| MHC_I | HLA-A | CD8B | 6 | -0.005 | 1 | 1 |
| MHC_I | HLA-B | CD8A | 10 | +0.015 | 1 | 1 |
| MHC_I | HLA-C | CD8A | 10 | +0.019 | 0.77 | 1 |
| MHC_I | HLA-A | CD8A | 10 | +0.032 | 0.492 | 0.803 |
| MHC_I | HLA-E | KLRD1 | 10 | +0.041 | 0.275 | 0.69 |
| MHC_I | HLA-E | KLRC1+KLRD1 | 3 | +0.075 | 0.25 | 0.662 |
| MHC_I | HLA-E | KLRC2 | 3 | +0.075 | 0.25 | 0.662 |
| MHC_I | HLA-E | KLRC1 | 4 | +0.092 | 0.125 | 0.602 |
| MHC_I_partial | HLA-F | VSIR | 5 | +0.058 | 0.0625 | 0.473 |
| MHC_I_partial | HLA-E | VSIR | 8 | +0.060 | 0.195 | 0.662 |
| MHC_I_partial | CEACAM5 | CD8A | 5 | +0.173 | 0.125 | 0.602 |
| T_recruit | CCL5 | CCR1 | 3 | -0.011 | 1 | 1 |
| T_recruit | CCL5 | CCR5 | 5 | +0.015 | 0.438 | 0.803 |
| T_recruit | CXCL16 | CXCR6 | 4 | +0.057 | 0.25 | 0.662 |
| T_recruit_partial | CXCL14 | CXCR4 | 12 | -0.008 | 0.91 | 1 |
| checkpoint | PVR | TIGIT | 5 | -0.025 | 0.625 | 0.896 |
| checkpoint | NECTIN2 | TIGIT | 11 | +0.015 | 0.00977 | 0.166 |
| checkpoint_partial | PVR | CD96 | 6 | -0.025 | 0.312 | 0.69 |
| checkpoint_partial | LGALS9 | P4HB | 8 | +0.007 | 0.945 | 1 |

Among the 21 focus/partial pairs with ≥3 paired patients: 7 have median Δ < 0; **0 reach FDR < 0.05**.

## LIANA CellPhoneDB (secondary, pooled / downsampled)

Edges with `source=Malig_CLDN4high` and `target=TNK` after ≤2,000 cells/group and 50 permutations.
These p-values are within-object specificity, not patient-level tests.

| ligand | receptor | lr_means | cellphone_pvals |
|---|---|---:|---:|
| APP | CD74 | 1.103 | 0 |
| CD58 | CD2 | 0.927 | 0 |
| SPP1 | CD44 | 0.810 | 0.02 |
| SPP1 | ITGA4_ITGB1 | 0.629 | 0 |
| HBEGF | CD44 | 0.498 | 0 |
| NECTIN1 | CD96 | 0.462 | 0 |
| AREG | ICAM1 | 0.278 | 0 |
| JAG1 | CD46 | 0.274 | 1 |

## Readout

Pooled, 25 pairs pass 10% expression from CLDN4-high malignant to T/NK (top: APP–CD74, HLA-A–CD8A, CD58–CD2, HLA-B–CD8A, PPIA–BSG). Patient-level high vs low uses n=17 patients; 0 focus-axis pairs reach FDR < 0.05. T/NK in this object is CD8T+Tprolif (NK=0). This is the rank in public TISCH2 GSE148071, not a CellChat run.

## Limits

- TISCH2 labels and log-normalized matrix, not a GEO raw-UMI / CellRanger reprocess.
- **No NK cells** are present under TISCH2 major-lineage in GSE148071; T/NK = CD8T + Tprolif.
- Tprolif is a cycling T-cell bin, not a separate lineage proof.
- Malignant is the TISCH2 call, not a re-run of inferCNV/CopyKAT.
- Pooled LR ranks mix patients; cite the paired-n table for inference.
- Chemokine dropout is high; read `ligand_frac` / `n_patients` with every rank.
- CellChat was not run.

## Files

| File | Role |
| --- | --- |
| `results/n_cells_patients.tsv` | Honest per-patient cell counts |
| `results/lr_cldn4high_to_tnk.tsv` | Full pooled LR table |
| `results/lr_cldn4high_to_tnk_pass.tsv` | Pairs passing 10% expr_prop |
| `results/patient_cldn4_outgoing.tsv` | Per-patient high vs low scores |
| `results/ranks_cldn4_outgoing.tsv` | Patient-level Wilcoxon ranks |
| `results/liana_cellphonedb.csv` | Full LIANA CellPhoneDB output (125 edges) |
| `results/liana_cldn4high_to_tnk.csv` | LIANA edges with source=Malig_CLDN4high, target=TNK |
| `results/summary.json` | Machine-readable n and method flags |
| `results/figures/` | n-cell bars and top-pair plot |

Generated 2026-08-17T17:16:39.169229+00:00.

