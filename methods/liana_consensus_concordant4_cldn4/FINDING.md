# FINDING — LIANA+ / CellPhoneDB / Connectome consensus, concordant-4

ADDITIVE. **CLDN4 only. No dual-high.** Concordant four only
(GSE123902 + GSE131907 + GSE205335 + GSE189357).
Do **not** add GSE148071 / GSE127465 / GSE154826 / GSE200563 / E-MTAB-13526.
This is an expression ligand–receptor contrast. It is **not** a spatial exclusion test.

Engine: Python liana 1.10.0 `rank_aggregate`
(CellPhoneDB + Connectome + log2FC + NATMI + SingleCellSignalR; RobustRankAggregate).
The three reported scores, from that one call, are:

- **CellPhoneDB** magnitude `lr_means` (Δ = high − low).
- **Connectome** magnitude `expr_prod` (Δ = high − low).
- **LIANA+** `magnitude_rank` (RobustRankAggregate ρ; Δ = ρ_low − ρ_high, so positive means CLDN4-high ranks stronger).

Senders = malignant CLDN4 Q4 vs Q1. Receivers = T/NK and myeloid, run separately.
Honest n = patient / locked sample. Primary test = Wilcoxon signed-rank of Δ across units.
BH is within method × receiver across the pre-specified edge panel, and separately
across the family tests, and within each unit for CellPhoneDB permutation p-values.

Pre-specified classes:

- **Barrier / exclusion** (expect high > low): junction and inhibitory edges
  (F11R, NECTIN2/PVR–TIGIT/CD96, CDH1, LGALS9, PD-1 ligands, HLA-E–NKG2A,
  CD47–SIRPA, CD24–SIGLEC10, TGFB1, MIF–CD74).
- **Effector recruitment** (expect low > high; the KD-like arm): CXCL9/10/11–CXCR3,
  CXCL16–CXCR6, CCL5–CCR1/CCR5, CX3CL1–CX3CR1.
- **Myeloid recruitment** (no directional thesis): CCL2/CCL7–CCR2, CXCL1/2/8–CXCR1/2,
  CCL3–CCR1, CSF1–CSF1R. Direction is reported, not scored as pass/fail.
- HLA–CD8 and CXCL12–CXCR4 are scored in the table and are not in the support ranks.

## Honest n

| gate | n | note |
|---|---:|---|
| Locked four | 65 | 13+21+22+9 |
| Inventory units | 65 | GSE123902=13, GSE131907=21, GSE189357=9, GSE205335=22 |
| Q4 eligible (n_mal≥40 and a receiver ≥20) | 64 | not the ligand n |

GSE123902 / GSE189357: epithelium marker-malignant (EPCAM\|KRT8\|KRT18\|KRT19 > 0 and PTPRC = 0);
T/NK = CD3D\|CD3E\|CD8A\|NKG7\|GNLY\|KLRD1 > 0 and not malignant;
myeloid = LYZ\|CD68\|CD14\|FCGR3A\|CSF1R\|AIF1 > 0 and not malignant or T/NK.
GSE131907 / GSE205335: author malignant, T/NK, and myeloid labels. Normal tissue in GSE205335 is out.
TACSTD2 is not a gate. PVRL2→NECTIN2, JAM1→F11R, IL8→CXCL8 when the official symbol is absent.
Label QC against the CellChat inventory: n_mal and n_tnk match on every locked unit.

Each sender or receiver group is capped at 400 cells (seed 1337) inside LIANA.
CellPhoneDB permutations: n_perms=1000. Expression is log1p CP10k using the full-library size.

## Family tests (pre-specified; BH across this 18-test family)

Δ is the within-unit mean of detected edge deltas. Positive = higher from CLDN4-high.

| receiver | class | method | n | 123902 | 131907 | 205335 | 189357 | mean Δ | p_W | q_BH | observed | agrees |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|---|---|
| TNK | barrier_exclusion | cpdb | 64 | 13 | 21 | 21 | 9 | +0.060 | 1.77e-10 | 1.60e-09 | high>low | yes |
| TNK | barrier_exclusion | conn | 64 | 13 | 21 | 21 | 9 | +0.049 | 3.25e-10 | 1.81e-09 | high>low | yes |
| TNK | barrier_exclusion | liana | 64 | 13 | 21 | 21 | 9 | +0.068 | 4.03e-10 | 1.81e-09 | high>low | yes |
| TNK | recruit_effector | cpdb | 64 | 13 | 21 | 21 | 9 | +0.009 | 1.02e-05 | 2.30e-05 | high>low | opposite |
| TNK | recruit_effector | conn | 64 | 13 | 21 | 21 | 9 | +0.005 | 7.18e-06 | 1.85e-05 | high>low | opposite |
| TNK | recruit_effector | liana | 64 | 13 | 21 | 21 | 9 | +0.008 | 0.0558 | 0.0836 | high>low | opposite |
| TNK | recruit_myeloid | cpdb | 64 | 13 | 21 | 21 | 9 | +0.000 | 0.893 | 0.926 | high>low | observed |
| TNK | recruit_myeloid | conn | 64 | 13 | 21 | 21 | 9 | +0.000 | 0.5 | 0.693 | high>low | observed |
| TNK | recruit_myeloid | liana | 64 | 13 | 21 | 21 | 9 | -0.001 | 0.593 | 0.762 | low>high | observed |
| MYE | barrier_exclusion | cpdb | 64 | 13 | 21 | 21 | 9 | +0.044 | 7.02e-11 | 1.26e-09 | high>low | yes |
| MYE | barrier_exclusion | conn | 64 | 13 | 21 | 21 | 9 | +0.040 | 3.25e-06 | 9.74e-06 | high>low | yes |
| MYE | barrier_exclusion | liana | 64 | 13 | 21 | 21 | 9 | +0.074 | 1.37e-08 | 4.91e-08 | high>low | yes |
| MYE | recruit_effector | cpdb | 64 | 13 | 21 | 21 | 9 | -0.002 | 0.926 | 0.926 | low>high | yes |
| MYE | recruit_effector | conn | 64 | 13 | 21 | 21 | 9 | -0.001 | 0.75 | 0.844 | low>high | yes |
| MYE | recruit_effector | liana | 64 | 13 | 21 | 21 | 9 | +0.002 | 0.677 | 0.812 | high>low | opposite |
| MYE | recruit_myeloid | cpdb | 64 | 13 | 21 | 21 | 9 | +0.009 | 0.00206 | 0.00413 | high>low | observed |
| MYE | recruit_myeloid | conn | 64 | 13 | 21 | 21 | 9 | +0.003 | 0.014 | 0.0252 | high>low | observed |
| MYE | recruit_myeloid | liana | 64 | 13 | 21 | 21 | 9 | +0.017 | 0.0241 | 0.0395 | high>low | observed |

## How to read the ranks

Barrier/exclusion is the class that replicates across methods and receivers.
Effector chemokines CXCL9/10/11–CXCR3 and CCL5 do not.
The effector-family mean on T/NK is pulled by CXCL16–CXCR6, which is higher from CLDN4-high.
Myeloid-recruitment edges do not clear the two-method bar one by one; the family mean toward myeloid is higher from CLDN4-high.

- T/NK barrier hits (high>low, ≥2 methods): 11 (NECTIN2-CD96, NECTIN2-TIGIT, CDH1-integrin, F11R-LFA1, CDH1-KLRG1, PVR-CD96, PVR-TIGIT, LGALS9-CD44, LGALS9-CD45, HLA-E-NKG2A, MIF-CD74).
- myeloid barrier hits (high>low, ≥2 methods): 9 (CD47-SIRPA, F11R-F11R, CD24-SIGLEC10, F11R-LFA1, CDH1-integrin, LGALS9-TIM3, LGALS9-CD44, LGALS9-CD45, MIF-CD74).
- T/NK CXCL9-CXCR3: n=49, Δ CPDB +0.002 (q=0.375), Connectome +0.000 (q=0.393), LIANA+ +0.007 (q=0.492), call=none.
- T/NK CXCL10-CXCR3: n=58, Δ CPDB +0.004 (q=0.375), Connectome -0.002 (q=0.657), LIANA+ +0.001 (q=0.917), call=none.
- T/NK CXCL11-CXCR3: n=34, Δ CPDB +0.001 (q=0.677), Connectome +0.001 (q=0.677), LIANA+ -0.001 (q=0.725), call=none.
- T/NK CCL5-CCR5: n=62, Δ CPDB -0.003 (q=0.303), Connectome -0.003 (q=0.303), LIANA+ -0.004 (q=0.601), call=none.
- T/NK CXCL16-CXCR6: n=62, Δ CPDB +0.047 (q=9.35e-06), Connectome +0.031 (q=6.61e-06), LIANA+ +0.037 (q=0.0543), call=high>low.
- Myeloid CSF1–CSF1R: n=64, Δ CPDB +0.045 (q=0.0121), Connectome +0.011 (q=0.0668), LIANA+ +0.059 (q=0.261), call=none.
- T/NK PD-L1–PD-1 is not a barrier hit: n=62, Δ CPDB +0.008 (q=0.213), Connectome +0.003 (q=0.123), LIANA+ +0.002 (q=0.725), call=none.
- T/NK TGFB1–TGFBR is not counted as barrier support. n=64, Δ CPDB -0.028 (q=0.0187), Connectome -0.008 (q=0.0117), LIANA+ -0.019 (q=0.3), call=low>high.
- T/NK HLA-A-CD8A is outside the support ranks (MHC continuity, not recruitment): n=64, Δ CPDB +0.147 (q=2.69e-05), Connectome +0.187 (q=1.47e-05), LIANA+ +0.058 (q=4.43e-04), call=high>low.
- T/NK HLA-B-CD8A is outside the support ranks (MHC continuity, not recruitment): n=64, Δ CPDB +0.123 (q=2.70e-04), Connectome +0.170 (q=8.35e-05), LIANA+ +0.056 (q=8.83e-04), call=high>low.
- T/NK HLA-C-CD8A is outside the support ranks (MHC continuity, not recruitment): n=64, Δ CPDB +0.113 (q=5.27e-04), Connectome +0.155 (q=4.11e-05), LIANA+ +0.053 (q=0.00788), call=high>low.
- CellPhoneDB is not the same sign in every cohort with ≥3 units for: MYE MIF-CD74 (GSE123902:+;GSE131907:-;GSE189357:+;GSE205335:+); TNK HLA-E-NKG2A (GSE123902:+;GSE131907:+;GSE189357:-;GSE205335:+); TNK MIF-CD74 (GSE123902:+;GSE131907:-;GSE189357:+;GSE205335:+).

## Consensus hits

A hit needs ≥2 of {CellPhoneDB, Connectome, LIANA+} with BH q<0.05, the same sign,
and no method significant in the opposite direction. LIANA+ ρ already uses CellPhoneDB
and Connectome magnitudes, so the three columns are not independent votes.
The call is corroboration of the rank aggregate by the two component magnitudes.

| receiver | edge | class | supports | n | Δ CPDB | q | Δ Connectome | q | Δ LIANA+ ρ | q | z | sig |
|---|---|---|---|---:|---:|---|---:|---|---:|---|---:|---|
| TNK | NECTIN2-CD96 | barrier_exclusion | barrier_exclusion | 64 | +0.142 | 2.03e-09 | +0.126 | 2.25e-09 | +0.242 | 3.61e-08 | +11.05 | cpdb+conn+liana |
| MYE | CD47-SIRPA | barrier_exclusion | barrier_exclusion | 64 | +0.141 | 4.23e-09 | +0.087 | 2.30e-08 | +0.160 | 3.44e-08 | +10.67 | cpdb+conn+liana |
| TNK | NECTIN2-TIGIT | barrier_exclusion | barrier_exclusion | 63 | +0.126 | 8.65e-09 | +0.115 | 7.27e-09 | +0.198 | 6.28e-07 | +10.30 | cpdb+conn+liana |
| TNK | CDH1-integrin | barrier_exclusion | barrier_exclusion | 64 | +0.150 | 3.64e-08 | +0.085 | 3.64e-08 | +0.182 | 6.28e-07 | +9.95 | cpdb+conn+liana |
| MYE | F11R-F11R | barrier_exclusion | barrier_exclusion | 64 | +0.097 | 4.04e-07 | +0.031 | 5.97e-07 | +0.200 | 1.89e-06 | +9.32 | cpdb+conn+liana |
| TNK | F11R-LFA1 | barrier_exclusion | barrier_exclusion | 64 | +0.118 | 2.35e-07 | +0.094 | 1.59e-07 | +0.165 | 7.36e-06 | +9.27 | cpdb+conn+liana |
| MYE | CD24-SIGLEC10 | barrier_exclusion | barrier_exclusion | 64 | +0.133 | 3.21e-06 | +0.043 | 3.69e-06 | +0.147 | 2.45e-06 | +8.75 | cpdb+conn+liana |
| MYE | F11R-LFA1 | barrier_exclusion | barrier_exclusion | 63 | +0.096 | 5.56e-06 | +0.034 | 5.56e-06 | +0.220 | 1.07e-05 | +8.38 | cpdb+conn+liana |
| TNK | CDH1-KLRG1 | barrier_exclusion | barrier_exclusion | 64 | +0.120 | 3.34e-06 | +0.064 | 3.34e-06 | +0.146 | 5.87e-05 | +8.34 | cpdb+conn+liana |
| MYE | CDH1-integrin | barrier_exclusion | barrier_exclusion | 63 | +0.082 | 3.49e-05 | +0.029 | 3.49e-05 | +0.162 | 3.49e-05 | +7.72 | cpdb+conn+liana |
| TNK | PVR-CD96 | barrier_exclusion | barrier_exclusion | 63 | +0.063 | 9.59e-06 | +0.023 | 2.06e-05 | +0.025 | 0.0543 | +6.65 | cpdb+conn |
| MYE | LGALS9-TIM3 | barrier_exclusion | barrier_exclusion | 64 | +0.044 | 4.65e-04 | +0.026 | 0.00104 | +0.161 | 2.19e-04 | +6.61 | cpdb+conn+liana |
| TNK | PVR-TIGIT | barrier_exclusion | barrier_exclusion | 62 | +0.054 | 9.35e-06 | +0.021 | 2.06e-05 | +0.014 | 0.233 | +6.28 | cpdb+conn |
| MYE | LGALS9-CD44 | barrier_exclusion | barrier_exclusion | 64 | +0.122 | 0.00118 | +0.077 | 0.0013 | +0.170 | 0.00234 | +6.00 | cpdb+conn+liana |
| MYE | LGALS9-CD45 | barrier_exclusion | barrier_exclusion | 64 | +0.073 | 0.00117 | +0.045 | 0.00216 | +0.163 | 0.00266 | +5.88 | cpdb+conn+liana |
| TNK | LGALS9-CD44 | barrier_exclusion | barrier_exclusion | 64 | +0.092 | 9.72e-04 | +0.061 | 0.00226 | +0.108 | 0.0052 | +5.83 | cpdb+conn+liana |
| TNK | LGALS9-CD45 | barrier_exclusion | barrier_exclusion | 64 | +0.124 | 0.00217 | +0.089 | 0.00341 | +0.113 | 0.00758 | +5.52 | cpdb+conn+liana |
| TNK | HLA-E-NKG2A | barrier_exclusion | barrier_exclusion | 64 | +0.040 | 0.00337 | +0.031 | 0.00418 | +0.035 | 0.0543 | +4.95 | cpdb+conn |
| TNK | MIF-CD74 | barrier_exclusion | barrier_exclusion | 64 | +0.051 | 0.0129 | +0.204 | 0.0101 | +0.013 | 0.218 | +4.16 | cpdb+conn |
| MYE | MIF-CD74 | barrier_exclusion | barrier_exclusion | 64 | +0.049 | 0.0173 | +0.363 | 0.0268 | +0.017 | 0.284 | +3.79 | cpdb+conn |
| MYE | CX3CL1-CX3CR1 | recruit_effector | recruitment_from_high | 56 | +0.009 | 0.0343 | +0.003 | 0.0377 | +0.024 | 0.254 | +3.66 | cpdb+conn |
| TNK | CXCL16-CXCR6 | recruit_effector | recruitment_from_high | 62 | +0.047 | 9.35e-06 | +0.031 | 6.61e-06 | +0.037 | 0.0543 | +6.88 | cpdb+conn |

## Ranked within class (top 5 by support score, including non-hits)

Barrier support score = consensus z (high>low). Effector-recruitment support score = −z
(low>high). Myeloid-recruitment support score = |z|.

| receiver | edge | class | supports | n | Δ CPDB | q | Δ Connectome | q | Δ LIANA+ ρ | q | z | sig |
|---|---|---|---|---:|---:|---|---:|---|---:|---|---:|---|
| TNK | NECTIN2-CD96 | barrier_exclusion | barrier_exclusion | 64 | +0.142 | 2.03e-09 | +0.126 | 2.25e-09 | +0.242 | 3.61e-08 | +11.05 | cpdb+conn+liana |
| TNK | NECTIN2-TIGIT | barrier_exclusion | barrier_exclusion | 63 | +0.126 | 8.65e-09 | +0.115 | 7.27e-09 | +0.198 | 6.28e-07 | +10.30 | cpdb+conn+liana |
| TNK | CDH1-integrin | barrier_exclusion | barrier_exclusion | 64 | +0.150 | 3.64e-08 | +0.085 | 3.64e-08 | +0.182 | 6.28e-07 | +9.95 | cpdb+conn+liana |
| TNK | F11R-LFA1 | barrier_exclusion | barrier_exclusion | 64 | +0.118 | 2.35e-07 | +0.094 | 1.59e-07 | +0.165 | 7.36e-06 | +9.27 | cpdb+conn+liana |
| TNK | CDH1-KLRG1 | barrier_exclusion | barrier_exclusion | 64 | +0.120 | 3.34e-06 | +0.064 | 3.34e-06 | +0.146 | 5.87e-05 | +8.34 | cpdb+conn+liana |
| TNK | CCL5-CCR5 | recruit_effector | none | 62 | -0.003 | 0.303 | -0.003 | 0.303 | -0.004 | 0.601 | -1.82 | nan |
| TNK | CCL5-CCR1 | recruit_effector | none | 63 | -0.000 | 0.677 | -0.001 | 0.677 | -0.000 | 0.725 | -0.77 | nan |
| TNK | CXCL11-CXCR3 | recruit_effector | none | 34 | +0.001 | 0.677 | +0.001 | 0.677 | -0.001 | 0.725 | +0.26 | nan |
| TNK | CXCL10-CXCR3 | recruit_effector | none | 58 | +0.004 | 0.375 | -0.002 | 0.657 | +0.001 | 0.917 | +0.30 | nan |
| TNK | CXCL9-CXCR3 | recruit_effector | none | 49 | +0.002 | 0.375 | +0.000 | 0.393 | +0.007 | 0.492 | +1.73 | nan |
| TNK | CSF1-CSF1R | recruit_myeloid | none | 56 | -0.002 | 0.368 | -0.000 | 0.657 | -0.004 | 0.725 | -1.18 | nan |
| TNK | CCL3-CCR1 | recruit_myeloid | none | 63 | +0.002 | 0.253 | +0.001 | 0.253 | -0.002 | 0.492 | +0.97 | nan |
| TNK | CCL2-CCR2 | recruit_myeloid | none | 52 | +0.000 | NA | +0.000 | NA | +0.000 | NA | NA | nan |
| TNK | CXCL8-CXCR2 | recruit_myeloid | none | 50 | +0.000 | NA | +0.000 | NA | +0.000 | NA | NA | nan |
| TNK | CXCL8-CXCR1 | recruit_myeloid | none | 37 | +0.000 | NA | +0.000 | NA | +0.000 | NA | NA | nan |
| MYE | CD47-SIRPA | barrier_exclusion | barrier_exclusion | 64 | +0.141 | 4.23e-09 | +0.087 | 2.30e-08 | +0.160 | 3.44e-08 | +10.67 | cpdb+conn+liana |
| MYE | F11R-F11R | barrier_exclusion | barrier_exclusion | 64 | +0.097 | 4.04e-07 | +0.031 | 5.97e-07 | +0.200 | 1.89e-06 | +9.32 | cpdb+conn+liana |
| MYE | CD24-SIGLEC10 | barrier_exclusion | barrier_exclusion | 64 | +0.133 | 3.21e-06 | +0.043 | 3.69e-06 | +0.147 | 2.45e-06 | +8.75 | cpdb+conn+liana |
| MYE | F11R-LFA1 | barrier_exclusion | barrier_exclusion | 63 | +0.096 | 5.56e-06 | +0.034 | 5.56e-06 | +0.220 | 1.07e-05 | +8.38 | cpdb+conn+liana |
| MYE | CDH1-integrin | barrier_exclusion | barrier_exclusion | 63 | +0.082 | 3.49e-05 | +0.029 | 3.49e-05 | +0.162 | 3.49e-05 | +7.72 | cpdb+conn+liana |
| MYE | CCL5-CCR5 | recruit_effector | none | 64 | -0.009 | 0.163 | -0.001 | 0.216 | +0.004 | 0.676 | -1.57 | nan |
| MYE | CCL5-CCR1 | recruit_effector | none | 64 | -0.012 | 0.663 | -0.007 | 0.375 | -0.015 | 0.745 | -1.04 | nan |
| MYE | CX3CL1-CX3CR1 | recruit_effector | recruitment_from_high | 56 | +0.009 | 0.0343 | +0.003 | 0.0377 | +0.024 | 0.254 | +3.66 | cpdb+conn |
| MYE | CXCL9-CXCR3 | recruit_effector | none | 56 | +0.000 | NA | +0.000 | NA | +0.000 | NA | NA | nan |
| MYE | CXCL10-CXCR3 | recruit_effector | none | 58 | +0.000 | NA | +0.000 | NA | +0.000 | NA | NA | nan |
| MYE | CSF1-CSF1R | recruit_myeloid | none | 64 | +0.045 | 0.0121 | +0.011 | 0.0668 | +0.059 | 0.261 | +3.72 | cpdb |
| MYE | CCL3-CCR1 | recruit_myeloid | none | 64 | +0.012 | 0.421 | +0.009 | 0.37 | +0.041 | 0.291 | +1.89 | nan |
| MYE | CXCL8-CXCR2 | recruit_myeloid | none | 59 | +0.002 | 0.392 | +0.000 | 0.37 | +0.009 | 0.37 | +1.73 | nan |
| MYE | CXCL2-CXCR2 | recruit_myeloid | none | 59 | +0.003 | 0.392 | +0.001 | 0.37 | +0.013 | 0.37 | +1.73 | nan |
| MYE | CCL7-CCR2 | recruit_myeloid | none | 54 | +0.000 | 0.392 | +0.000 | 0.37 | +0.005 | 0.37 | +1.73 | nan |

Full rank: `results/tables/consensus_ranked.tsv`.
Per-unit scores: `results/tables/per_patient_edges.tsv`.
CellPhoneDB permutation BH (within unit × sender × receiver): `results/tables/cpdb_specificity_bh.tsv`.
Figure: `results/figures/consensus_barrier_vs_recruitment.png`.

## What is not claimed

- This does not measure spatial exclusion, contact, or muzzling.
- TACSTD2 is not used. This is not dual-high.
- GSE148071, GSE127465, and the other non-concordant sets are not added.
- Stouffer z ranks edges. It is not a calibrated meta-analytic p-value:
  the three scores share the same cells, and LIANA+ already aggregates CellPhoneDB and Connectome.
- Myeloid chemokine direction was not given a pass/fail expectation.
- Cell-pooled tests are not the honest n.

## Reproduce

```bash
bash methods/liana_consensus_concordant4_cldn4/scripts/download.sh /tmp/concordant4_raw
python3 methods/liana_consensus_concordant4_cldn4/scripts/run_consensus.py
```

liana 1.10.0; anndata 0.13.4; scipy 1.18.1.
GSE205335 RDS is read in R (Matrix) because it is a double-gzipped dgCMatrix.
