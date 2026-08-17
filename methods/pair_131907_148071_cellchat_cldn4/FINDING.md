# FINDING — pair GSE131907 + GSE148071: malignant CLDN4 vs T/NK + CellChat

ADDITIVE. **CLDN4 only.** No dual-high TACSTD2×CLDN4. **GSE205335 is not
in this pair.** Patient is the unit. Combo ρ is DerSimonian–Laird on
Fisher-z of the two within-cohort Spearmans. p-values are descriptive.
Prior single-dataset CellChat (PR #347, #348) and the Q4 meta that
pairs GSE131907 with GSE205335 (PR #320) are given and are not re-ranked.

## Honest n

- GSE131907 (Kim et al. 2020): locked author-`Malignant cells` extract
  (PR #320) with ≥20 malignant **and** ≥20 T/NK. After the floor this
  is **21 patients** (tL/B + mLN + mBrain; one sample each). tLung is
  out of this extract because primary tumor epithelium is labeled tS*
  not `Malignant cells`. PE / nLung / nLN out. No ICI / MPR labels.
- GSE148071 (Wu et al. 2021): **25 / 42** biopsies with
  ≥25 marker-argmax epithelial (putative malignant) **and**
  ≥25 T/NK. No histology / ICI labels. Epithelium is
  putative; CopyKAT IDs are not on GEO.
- Combo N = **46**. Q4 vs Q1 uses **within-cohort** tails:
  **n=13 vs 11** (n_compared=24), not the full N.
- Stacked 10x + Singleron Spearman is a companion only. It is not the combo ρ.

## Combo ρ (primary deliverable)

| analysis | combo | k | N | effect (p, I²) |
|---|---|---:|---:|---|
| Spearman %pos vs T/NK | GSE131907+GSE148071 | 2 | 46 | ρ=-0.242 (0.446, 76%) |
| Spearman mean vs T/NK | GSE131907+GSE148071 | 2 | 46 | ρ=-0.105 (0.73, 73%) |
| Q4 vs Q1 r %pos | GSE131907+GSE148071 | 2 | 24 compared | r=-0.317 (0.342, 53%) |

Stacked companion (not the combo): ρ=-0.172 p=0.253 n=46.

Primary cut is **%pos**. Mean is the same patients, secondary.

## Per-cohort malignant CLDN4 vs same-patient T/NK

| cohort | score | n | Spearman ρ (p) | Q4 vs Q1 r (p; n_Q1/n_Q4) | median T/NK Q1 | median T/NK Q4 | Δ |
|---|---|---:|---|---|---:|---:|---:|
| GSE131907 | pct_pos | 21 | -0.522 (0.0152) | -0.600 (0.126; 6/5) | 0.390 | 0.048 | -0.342 |
| GSE148071 | pct_pos | 25 | +0.069 (0.742) | +0.000 (1; 7/6) | 0.080 | 0.074 | -0.006 |
| GSE131907 | mean | 21 | -0.396 (0.0755) | -0.533 (0.177; 6/5) | 0.346 | 0.105 | -0.240 |
| GSE148071 | mean | 25 | +0.189 (0.365) | -0.048 (0.945; 7/6) | 0.080 | 0.092 | +0.012 |

### Quartile tails (within-cohort CLDN4 %pos)

| tail | cohort | unit | origin | n_mal | n_TNK | CLDN4 %pos | T/NK frac |
|---|---|---|---|---:|---:|---:|---:|
| Q1 | GSE131907 | P3006 | mBrain | 83 | 383 | 60.2 | 0.358 |
| Q1 | GSE131907 | P3002 | mBrain | 181 | 977 | 58.6 | 0.422 |
| Q1 | GSE131907 | P1049 | tL/B | 114 | 1126 | 49.1 | 0.634 |
| Q1 | GSE131907 | P1015 | mLN | 256 | 479 | 47.3 | 0.433 |
| Q1 | GSE131907 | P3016 | mBrain | 79 | 385 | 11.4 | 0.292 |
| Q1 | GSE131907 | P1013 | mLN | 376 | 1150 | 0.3 | 0.324 |
| Q1 | GSE148071 | P27 | biopsy | 34 | 87 | 50.0 | 0.123 |
| Q1 | GSE148071 | P24 | biopsy | 88 | 64 | 47.7 | 0.166 |
| Q1 | GSE148071 | P23 | biopsy | 1466 | 32 | 40.8 | 0.013 |
| Q1 | GSE148071 | P14 | biopsy | 690 | 31 | 34.9 | 0.017 |
| Q1 | GSE148071 | P4 | biopsy | 112 | 304 | 25.9 | 0.067 |
| Q1 | GSE148071 | P42 | biopsy | 151 | 454 | 9.3 | 0.454 |
| Q1 | GSE148071 | P11 | biopsy | 106 | 46 | 2.8 | 0.080 |
| Q4 | GSE131907 | P1028 | tL/B | 4640 | 232 | 96.3 | 0.045 |
| Q4 | GSE131907 | P1019 | mLN | 217 | 1427 | 94.0 | 0.676 |
| Q4 | GSE131907 | P3003 | mBrain | 2713 | 145 | 93.1 | 0.048 |
| Q4 | GSE131907 | P3007 | mBrain | 5108 | 166 | 90.4 | 0.029 |
| Q4 | GSE131907 | P3004 | mBrain | 1229 | 131 | 89.3 | 0.069 |
| Q4 | GSE148071 | P8 | biopsy | 422 | 66 | 90.3 | 0.049 |
| Q4 | GSE148071 | P32 | biopsy | 381 | 79 | 90.0 | 0.133 |
| Q4 | GSE148071 | P22 | biopsy | 153 | 231 | 86.3 | 0.359 |
| Q4 | GSE148071 | P13 | biopsy | 165 | 90 | 83.6 | 0.096 |
| Q4 | GSE148071 | P12 | biopsy | 288 | 46 | 81.9 | 0.052 |
| Q4 | GSE148071 | P28 | biopsy | 693 | 64 | 81.8 | 0.027 |

Q4 vs Q1 is within-cohort so 10x and Singleron are not ranked on one scale.

## CellChat-style ligands

Not scored in this write-up (matrix step skipped). Re-run without
`--skip-cellchat` to fill `results/ligand_table.tsv`.

## Files

- `results/combo_rho.tsv` — combo ρ / Q4 vs Q1 r (this pair only)
- `results/q4q1_tnk.tsv` — per-cohort Spearman + Q4 vs Q1
- `results/patients_with_quartiles.tsv` — honest patient table
- `results/ligand_table.tsv` — CellChat-style outgoing Mal → T/NK
- `figures/fig_combo_rho.png` — extra combo-ρ forest
- `figures/q4q1_tnk_pct.png` — Q4 vs Q1 T/NK box
- `figures/scatter_cldn4_tnk.png` — extra scatter
- `figures/fig_extra_ligand_table.png` — extra ligand-table figure
- `METHODS.md` — floors, quartiles, Hill probability

