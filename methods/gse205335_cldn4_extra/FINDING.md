# FINDING — GSE205335 malignant CLDN4 vs same-patient T/NK Q4 + CellChat-style ligands

ADDITIVE. **CLDN4 only.** Patient is the unit. Prior TACSTD2 A3 and the
multi-cohort Q4 meta (PR #320) are given and are not re-ranked here.
p-values are descriptive.

## Honest n

Locked extract: **22 patients** with ≥20 author-malignant and ≥20 T/NK
cells (PR #279 / #320). Four GEO patients have (near-)zero captured
malignant cells and are already out (3 PR + 1 PD). Q4 vs Q1 uses the
quartile **tails only**: **n=6 vs 6** (n_compared=12), not 22.
MPR/NMPR is unlabeled; RECIST is not used as MPR.

## Malignant CLDN4 vs same-patient T/NK

| score | n | Spearman ρ (p) | Q4 vs Q1 r (p; n_Q1/n_Q4) | median T/NK Q1 | median T/NK Q4 | Δ |
|---|---:|---|---|---:|---:|---:|
| CLDN4 %pos | 22 | -0.435 (0.0429) | -0.778 (0.026; 6/6) | 0.514 | 0.116 | -0.398 |
| CLDN4 mean log1p(CP10k) | 22 | -0.200 (0.371) | -0.556 (0.132; 6/6) | 0.331 | 0.166 | -0.165 |

Primary cut is **%pos** (same row as PR #320 GSE205335 single:
r=−0.778, p=0.026, 6/6). Mean is the same 22 patients, weaker.
NSCLC-only (ADC+SQ) is not the verdict; it was already thinner in the
given extract.

### Quartile tails (CLDN4 %pos)

| tail | patients (histology, RECIST) | n_mal | n_TNK | CLDN4 %pos | T/NK frac |
|---|---|---:|---:|---:|---:|
| Q4 | P1025 (SCLC, PD) | 980 | 128 | 92.9 | 0.105 |
| Q4 | P1115 (SCLC, PR) | 3769 | 271 | 88.9 | 0.062 |
| Q4 | P1089 (ADC, PD) | 1063 | 272 | 87.8 | 0.169 |
| Q4 | P1084 (ADC, NE) | 196 | 2350 | 86.7 | 0.665 |
| Q4 | P1037 (SQ, PR) | 3671 | 611 | 86.1 | 0.123 |
| Q4 | P1016 (SCLC, PR) | 5336 | 763 | 82.6 | 0.110 |
| Q1 | P1063 (ADC, NE) | 131 | 3144 | 35.1 | 0.732 |
| Q1 | P1119 (ADC, PD) | 1222 | 621 | 36.1 | 0.260 |
| Q1 | P1090 (SQ, PR) | 566 | 247 | 38.2 | 0.253 |
| Q1 | P1015 (ADC, PD) | 291 | 473 | 41.6 | 0.396 |
| Q1 | P1062 (ADC, PD) | 188 | 1934 | 45.2 | 0.633 |
| Q1 | P4001 (ADC, SD) | 27 | 1801 | 48.1 | 0.780 |

Q4 mixes SCLC (P1025, P1115, P1016) with ADC. That histology mix is
part of the honest n, not hidden.

## CellChat-style ligands

Not scored in this write-up (matrix step skipped or failed).
Re-run without `--skip-cellchat` to fill `results/ligand_table.tsv`.

## Files

- `results/q4q1_tnk.tsv` — Spearman + Q4 vs Q1 T/NK
- `results/patients_with_quartiles.tsv` — 22-patient table + quartile labels
- `results/ligand_table.tsv` — CellChat-style differential pairs
- `figures/q4q1_tnk_pct.png` — Q4 vs Q1 T/NK box
- `figures/fig_extra_ligand_table.png` — extra ligand-table figure
- `METHODS.md` — quartiles, author labels, Hill probability

