# FINDING — Multi-sample NicheNet, CLDN4-high vs low malignant → same-patient T/NK

ADDITIVE **CLDN4-only**. No dual-high TACSTD2×CLDN4. Patient is the unit.
Sender = author-malignant cells split by CLDN4 `log1p(CP10k)` at the **cohort-wide
malignant median**. Receiver = **same-patient** T/NK. Run on the PR #320 winning
public merge **GSE131907+GSE205335** plus **GSE207422** as a third cohort.

**Method:** Documented NicheNet-v2 ligand–target prior (`ligand_target_matrix_nsga2r_final`, Zenodo 10.5281/zenodo.7074291) scored in Python (`rdata`); MultiNicheNet-style equal-weight prioritization (Bonte/Browaeys vignette). R packages `nichenetr` and `multinichenetr` were **not** installed.

## Verdict

Malignant CLDN4 %pos is negatively associated with same-patient T/NK fraction in the three-cohort patient pool (ρ=-0.420, p=0.000886, N=64). Ligand ranks below are prior + paired sender DE, not proof of causation.

## Honest n

| Item | n | Note |
| --- | ---: | --- |
| Program patients (all 3) | **64** | ≥20 malignant + ≥20 T/NK; unit of Spearman / Q4Q1 |
| Paired patients (all 3) | **60** | ≥10 CLDN4-high + ≥10 low + ≥20 T/NK; unit of ligand Δ |
| GSE131907 program / paired | 31 / 29 | tumor-origin samples pooled per `patient_id` |
| GSE205335 program / paired | 22 / 21 | tumor samples pooled; normals dropped |
| GSE207422 program / paired | 11 / 10 | post-tx DRMref; pCR counted as MPR in the note column |
| Winning merge program N | **53** | GSE131907+GSE205335 (PR #320 given) |
| Three-cohort program N | **64** | +GSE207422 |

Cells are counts, not n. GSE131907 is **patient-pooled** (mBrain/mLN/PE included); that is not the
PR #320 sample-level n=21 lock. GEO patient-ID strings can collide across studies — n is the
sum of per-cohort patients, not `nunique` across the merge.
Q4 vs Q1 uses quartile **tails only**; thin tails (n<8 or a tail <3) are flagged and not treated as a pool.

## Patient-level T/NK vs malignant CLDN4 %pos

| cohort | n | ρ T/NK frac (p) | Q4 vs Q1 r (p; n_Q1/n_Q4) | ρ cyto (p) | ρ IFN (p) |
| --- | ---: | --- | --- | --- | --- |
| GSE131907 | 31 | -0.343 (0.0588) | -0.406 (0.195; 8/8) | -0.103 (0.582) | -0.213 (0.251) |
| GSE205335 | 22 | -0.435 (0.0429) | -0.778 (0.026; 6/6) | -0.347 (0.113) | -0.239 (0.284) |
| GSE207422 | 11 | -0.618 (0.0426) | -0.778 (0.2; 3/3) | 0.009 (0.979) | 0.118 (0.729) |
| **GSE131907+GSE205335** | 53 | **-0.381 (0.0059, I²=0%)** | — | -0.205 (0.154) | -0.223 (0.12) |
| **+GSE207422 (3 cohorts)** | 64 | **-0.42 (0.000886, I²=0%)** | — | -0.175 (0.191) | -0.175 (0.19) |

p-values are descriptive. Fisher-z pool uses patient n per cohort (k=2 or 3).

## Ligand-activity table (primary)

Full table: `results/ligand_activity.tsv`. Rank = MultiNicheNet-style equal-weight
mean of min-max scaled (1) paired sender Δ, (2) NicheNet-v2 activity vs empirical
T/NK DE (CLDN4-high vs low **patients**), (3) receptor expression in T/NK,
(4) fraction of paired patients with ligand detected in ≥10% of CLDN4-high senders.

| rank | ligand | n patients | n cohorts | median Δ | paired p | activity r | prio |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | COPA | 60 | 3 | 0.077 | 1.81e-09 | 0.295 | 0.807 |
| 2 | APP | 60 | 3 | 0.138 | 3.24e-09 | 0.146 | 0.769 |
| 3 | CD59 | 60 | 3 | 0.133 | 3.78e-07 | 0.158 | 0.727 |
| 4 | HMGB1 | 60 | 3 | -0.001 | 0.503 | 0.132 | 0.694 |
| 5 | CD9 | 60 | 3 | 0.324 | 4.02e-11 | 0.219 | 0.684 |
| 6 | NECTIN2 | 31 | 2 | 0.239 | 1.86e-09 | 0.184 | 0.680 |
| 7 | LAMB2 | 60 | 3 | 0.033 | 3.82e-08 | 0.129 | 0.679 |
| 8 | CD58 | 60 | 3 | 0.019 | 3.99e-05 | 0.308 | 0.675 |
| 9 | CD47 | 60 | 3 | 0.143 | 1.73e-07 | 0.363 | 0.670 |
| 10 | F11R | 60 | 3 | 0.184 | 2.21e-11 | 0.114 | 0.664 |
| 11 | HLA-E | 60 | 3 | 0.233 | 2.12e-08 | 0.090 | 0.663 |
| 12 | CD55 | 60 | 3 | 0.222 | 7.20e-09 | 0.187 | 0.659 |
| 13 | HLA-DMA | 60 | 3 | 0.011 | 0.239 | 0.119 | 0.650 |
| 14 | CLDN3 | 60 | 3 | 0.366 | 1.63e-11 | 0.145 | 0.647 |
| 15 | LGALS9 | 60 | 3 | 0.025 | 0.0026 | 0.175 | 0.644 |

Activity Pearson is prior-structure (unsigned). A high IFN/MHC ligand score does
**not** mean CLDN4-high cells induce that program — test direction on the patient rows above.

## What this is not

- Not R `nichenetr` or `multinichenetr` (packages absent; Python prior scoring is used).
- Not dual-high TACSTD2×CLDN4. TACSTD2 is never a gate.
- Not CellChat / LIANA. Those live in other folders.
- Not a full-transcriptome NicheNet run. Background = panel ∩ prior ∩ T/NK-expressed genes.
- GSE131907 has no ICI/MPR labels. GSE205335 RECIST is not substituted for MPR.
- Cells are not replicates.

## Extra figures

- `figures/n_patients_inclusion.png`
- `figures/cldn4_vs_tnk_frac.png`
- `figures/cldn4_vs_tnk_programs.png`
- `figures/forest_spearman_tnk_frac.png`
- `figures/q4q1_tnk_frac_box.png`
- `figures/ligand_prioritization_top.png`
- `figures/ligand_activity_empirical.png`
- `figures/ligand_delta_top.png`
- `figures/extra_n_dropped.png`
- `figures/extra_ligand_activity_ifn.png`

## Reproduce

```bash
cd methods/multinichenet_cldn4
python3 scripts/00_download.py
python3 scripts/01_convert_prior.py
python3 scripts/02_extract_panels.py
python3 scripts/03_analyze.py
```

