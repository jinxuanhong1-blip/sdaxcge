# FINDING — leftover TISCH NSCLC merge, CLDN4-only combo + high-end

ADDITIVE. **CLDN4-only.** No dual-high TACSTD2×CLDN4. Public TISCH2
`expression.h5` + `CellMetainfo` only. Files >2 GB skipped. Series with
neither a CLDN4 compartment (epithelial/malignant) nor T/NK skipped.

Hunt list (not in the 131907/148071/205335/207422 set):
GSE117570, GSE139555, GSE146100, GSE149655, GSE143423, GSE127471.

Patient is the unit. Sample/nodule rows are a sensitivity (GSE146100 is
1 patient / 3 nodules). Eligible unit: ≥20 scored
malignant cells (else epithelial-like) and ≥20 T/NK.
Tumor-like tissue only. TISCH values are MAESTRO `log2(TPM/10+1)`.

p-values are descriptive. Combo table:
[`tables/highlighted_combos.tsv`](tables/highlighted_combos.tsv)
(empty header-only is an honest empty hunt).

## Verdict

Highlighted leftover combos: **none**. Empty hunt. No leftover series or merge reached ρ<0 with p<0.05, and no clear (non-thin) Q4 vs Q1 drop.

Patient-level leftover merge does **not** reach n≥8
(malignant-only n=4; mixed malignant+epi-like n=7).
Sample/nodule merge n=9 is labeled as samples, not patients.

- **malignant leftover patients:** n=4 (malignant 4, epi-like 0) · ρ not computed (n<5) · Q4 vs Q1 not computed · n≥8 merge=no · n=4<5; Spearman not a claim
- **mixed leftover patients:** n=7 (malignant 4, epi-like 3) · ρ=-0.250 p=0.589 · Q4 vs Q1 r=-1.000 p=0.333 (n_Q1=2, n_Q4=2; thin) · n≥8 merge=no · n=7<8; merge threshold not met; Q4 vs Q1 thin
- **mixed leftover samples/nodules:** n=9 (malignant 4, epi-like 5) · ρ=-0.333 p=0.381 · Q4 vs Q1 r=-0.333 p=0.800 (n_Q1=3, n_Q4=2; thin) · n≥8 merge=yes · Q4 vs Q1 thin

## Series audit

| Series | cells | Malignant | Epi-like | T/NK | eligible patients | eligible samples | CLDN4 in h5 | skip |
|---|---:|---:|---:|---:|---:|---:|---|---|
| GSE117570 | 11453 | 2721 | 501 | 3721 | 2 | 2 | yes | usable |
| GSE139555 | 78829 | 0 | 0 | 67655 | 0 | 0 | n/a | no epithelial/malignant cells; CLDN4 cannot be scored |
| GSE146100 | 10996 | 0 | 975 | 6854 | 1 | 3 | yes | usable |
| GSE149655 | 9591 | 0 | 4119 | 1776 | 2 | 2 | yes | usable |
| GSE143423 | 12193 | 9237 | 0 | 99 | 2 | 2 | yes | usable |
| GSE127471 | 1108 | 0 | 0 | 521 | 0 | 0 | n/a | no epithelial/malignant cells; CLDN4 cannot be scored |

GSE139555 is T/immune-sorted (T/NK present, no epithelium → CLDN4 not
scorable). GSE127471 is PBMC-only. Both fail the “CLDN4 present” gate;
h5 was not downloaded. GSE146100 and GSE149655 have **no TISCH
Malignant call** — they enter the mixed merge as epithelial-like only.

## Eligible leftover patients

| series | patient | def | n_scored | n_TNK | CLDN4 mean | CLDN4 %pos | frac T/NK |
|---|---|---|---:|---:|---:|---:|---:|
| GSE117570 | P1 | Malignant | 128 | 577 | 0.000 | 0.0 | 0.316 |
| GSE117570 | P2 | Malignant | 199 | 52 | 1.440 | 87.4 | 0.040 |
| GSE143423 | lbm2 | Malignant | 854 | 44 | 1.552 | 98.7 | 0.026 |
| GSE143423 | lbm3 | Malignant | 4217 | 42 | 1.674 | 84.6 | 0.007 |
| GSE146100 | Patient 1 | epithelial_like | 975 | 6854 | 1.634 | 86.7 | 0.623 |
| GSE149655 | Patient 01 | epithelial_like | 1277 | 52 | 1.755 | 81.2 | 0.029 |
| GSE149655 | Patient 02 | epithelial_like | 407 | 905 | 1.579 | 90.7 | 0.340 |

GSE117570 **P1** malignant CLDN4 mean and %pos are 0 in the TISCH h5
(128 scored malignant cells). That is a data fact, not a filter.
P3/P4 have malignant CLDN4 but fewer than 20 T/NK and are not eligible.

## Per-series Spearman (patient-level, CLDN4 mean vs T/NK)

Spearman only if n≥5. Smaller n is listed, not tested.

| series | def | n | ρ (p) | Q4 vs Q1 | note |
|---|---|---:|---|---|---|
| GSE117570 | Malignant | 2 | — (—) | — | n=2<5; Spearman not a claim |
| GSE143423 | Malignant | 2 | — (—) | — | n=2<5; Spearman not a claim |
| GSE146100 | epithelial_like | 1 | — (—) | — | n=1<5; Spearman not a claim |
| GSE149655 | epithelial_like | 2 | — (—) | — | n=2<5; Spearman not a claim |

### Malignant-only leftover combos (patient)

| k | N | cohorts | ρ mean (p) | ρ %pos (p) | Q4 vs Q1 r (p; n_Q1/n_Q4) | trigger | note |
|---:|---:|---|---|---|---|---|---|
| 1 | 2 | GSE117570 | — (—) | — (—) | — | no | n=2<5; Spearman not a claim |
| 1 | 2 | GSE143423 | — (—) | — (—) | — | no | n=2<5; Spearman not a claim |
| 2 | 4 | GSE117570+GSE143423 | -1.000 (0) | -0.400 (0.600) | — | no | n=4<5; Spearman not a claim |

### Mixed leftover combos (patient; malignant + epi-like)

| k | N | cohorts | ρ mean (p) | ρ %pos (p) | Q4 vs Q1 r (p; n_Q1/n_Q4) | trigger | note |
|---:|---:|---|---|---|---|---|---|
| 1 | 2 | GSE117570 | — (—) | — (—) | — | no | n=2<5; Spearman not a claim |
| 1 | 2 | GSE143423 | — (—) | — (—) | — | no | n=2<5; Spearman not a claim |
| 1 | 2 | GSE149655 | — (—) | — (—) | — | no | n=2<5; Spearman not a claim |
| 1 | 1 | GSE146100 | — (—) | — (—) | — | no | n=1<5; Spearman not a claim |
| 2 | 4 | GSE117570+GSE143423 | -1.000 (0) | -0.400 (0.600) | — | no | n=4<5; Spearman not a claim |
| 2 | 4 | GSE117570+GSE149655 | -0.400 (0.600) | +0.400 (0.600) | — | no | n=4<5; Spearman not a claim |
| 2 | 4 | GSE143423+GSE149655 | +0.000 (1.000) | +0.000 (1.000) | — | no | n=4<5; Spearman not a claim |
| 2 | 3 | GSE117570+GSE146100 | +0.500 (0.667) | -0.500 (0.667) | — | no | n=3<5; Spearman not a claim |
| 2 | 3 | GSE143423+GSE146100 | -0.500 (0.667) | +0.500 (0.667) | — | no | n=3<5; Spearman not a claim |
| 2 | 3 | GSE146100+GSE149655 | -0.500 (0.667) | +0.500 (0.667) | — | no | n=3<5; Spearman not a claim |
| 3 | 6 | GSE117570+GSE143423+GSE149655 | -0.429 (0.397) | -0.086 (0.872) | -1.000 (0.333; 2/2) thin | no | n=6<8; merge threshold not met; Q4 vs Q1 thin |
| 3 | 5 | GSE117570+GSE143423+GSE146100 | -0.400 (0.505) | -0.200 (0.747) | — | no | n=5<8; merge threshold not met |
| 3 | 5 | GSE117570+GSE146100+GSE149655 | -0.100 (0.873) | +0.300 (0.624) | — | no | n=5<8; merge threshold not met |
| 3 | 5 | GSE143423+GSE146100+GSE149655 | -0.100 (0.873) | +0.100 (0.873) | — | no | n=5<8; merge threshold not met |
| 4 | 7 | GSE117570+GSE143423+GSE146100+GSE149655 | -0.250 (0.589) | +0.000 (1.000) | -1.000 (0.333; 2/2) thin | no | n=7<8; merge threshold not met; Q4 vs Q1 thin |

### Malignant-only leftover combos (sample/nodule)

| k | N | cohorts | ρ mean (p) | ρ %pos (p) | Q4 vs Q1 r (p; n_Q1/n_Q4) | trigger | note |
|---:|---:|---|---|---|---|---|---|
| 1 | 2 | GSE117570 | — (—) | — (—) | — | no | n=2<5; Spearman not a claim |
| 1 | 2 | GSE143423 | — (—) | — (—) | — | no | n=2<5; Spearman not a claim |
| 2 | 4 | GSE117570+GSE143423 | -1.000 (0) | -0.400 (0.600) | — | no | n=4<5; Spearman not a claim |

### Mixed leftover combos (sample/nodule; not independent patients)

| k | N | cohorts | ρ mean (p) | ρ %pos (p) | Q4 vs Q1 r (p; n_Q1/n_Q4) | trigger | note |
|---:|---:|---|---|---|---|---|---|
| 1 | 3 | GSE146100 | -0.500 (0.667) | -0.500 (0.667) | — | no | n=3<5; Spearman not a claim |
| 1 | 2 | GSE117570 | — (—) | — (—) | — | no | n=2<5; Spearman not a claim |
| 1 | 2 | GSE143423 | — (—) | — (—) | — | no | n=2<5; Spearman not a claim |
| 1 | 2 | GSE149655 | — (—) | — (—) | — | no | n=2<5; Spearman not a claim |
| 2 | 5 | GSE117570+GSE146100 | +0.300 (0.624) | +0.000 (1.000) | — | no | n=5<8; merge threshold not met |
| 2 | 5 | GSE143423+GSE146100 | -0.600 (0.285) | +0.100 (0.873) | — | no | n=5<8; merge threshold not met |
| 2 | 5 | GSE146100+GSE149655 | -0.800 (0.104) | +0.100 (0.873) | — | no | n=5<8; merge threshold not met |
| 2 | 4 | GSE117570+GSE143423 | -1.000 (0) | -0.400 (0.600) | — | no | n=4<5; Spearman not a claim |
| 2 | 4 | GSE117570+GSE149655 | -0.400 (0.600) | +0.400 (0.600) | — | no | n=4<5; Spearman not a claim |
| 2 | 4 | GSE143423+GSE149655 | +0.000 (1.000) | +0.000 (1.000) | — | no | n=4<5; Spearman not a claim |
| 3 | 7 | GSE117570+GSE143423+GSE146100 | -0.286 (0.535) | +0.000 (1.000) | -0.500 (0.667; 2/2) thin | no | n=7<8; merge threshold not met; Q4 vs Q1 thin |
| 3 | 7 | GSE117570+GSE146100+GSE149655 | -0.214 (0.645) | +0.286 (0.535) | -0.500 (0.667; 2/2) thin | no | n=7<8; merge threshold not met; Q4 vs Q1 thin |
| 3 | 7 | GSE143423+GSE146100+GSE149655 | -0.500 (0.253) | +0.071 (0.879) | -1.000 (0.333; 2/2) thin | no | n=7<8; merge threshold not met; Q4 vs Q1 thin |
| 3 | 6 | GSE117570+GSE143423+GSE149655 | -0.429 (0.397) | -0.086 (0.872) | -1.000 (0.333; 2/2) thin | no | n=6<8; merge threshold not met; Q4 vs Q1 thin |
| 4 | 9 | GSE117570+GSE143423+GSE146100+GSE149655 | -0.333 (0.381) | +0.100 (0.798) | -0.333 (0.800; 3/2) thin | no | Q4 vs Q1 thin |

## CellChat-style

CellChat-style **not run**. Trigger is ρ<0 with p<0.05 or a clear Q4 vs Q1 drop on a leftover series/merge. That trigger did not fire.

## What was not done

- No dual-high TACSTD2×CLDN4 score.
- GSE131907 / GSE148071 / GSE205335 / GSE207422 were not re-scored.
- Other TISCH NSCLC objects already in `methods/tisch_nsclc_pool/`
  (GSE127465, GSE153935, GSE162498, GSE150660, EMTAB6149) were not
  added to this hunt list.
- Files >2 GB skipped. GSE139555 / GSE127471 h5 skipped (no CLDN4 compartment).
- Cell-level p-values are not primary evidence.

## Extra figures

- `fig_inventory_eligible_patients.png`
- `fig_lineage_bars.png`
- `fig_scatter_leftover_patients.png`
- `fig_scatter_leftover_samples.png`
- `fig_combo_patient_mixed_bars.png`
- `fig_combo_sample_mixed_bars.png`
- `fig_q4q1_largest_merge.png`
- `fig_patient_table.png`

Reproduce:

```bash
python3 -m pip install -r methods/tisch_leftover_merge_cldn4/requirements.txt
python3 methods/tisch_leftover_merge_cldn4/download.py
python3 methods/tisch_leftover_merge_cldn4/analyze.py
```

Generated: 2026-08-17T17:30:48.519825+00:00
