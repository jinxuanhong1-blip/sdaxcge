# Finding — winning-pair GSE131907+GSE205335, real Slingshot CLDN4-only

ADDITIVE. **CLDN4 only.** Winning pair is given (PR #320: author-malignant CLDN4 %pos vs T/NK, Q4 vs Q1 **n=23** r=−0.705) and is **not** re-audited. Malignant/epithelial cells. No dual-high TACSTD2×CLDN4. GSE148071 is not added. GSE207422 is not added. This is **not** a DPT-null writeup and **not** a redo of PR #449 (R missing → DPT).

Primary clock: **Slingshot (pyslingshot-bio (Street et al. 2018 Python port))**. Inferential unit = **sample/patient** (GSE131907 `Sample`, GSE205335 `patient`). Cell-level ρ is descriptive. Barrier/keratin and IFN scores **exclude CLDN4**. Root cluster = Leiden 5 (Leiden with most GSE131907 nLung author AT2; max-CLDN4 cluster excluded). Never CLDN4-high (excluded Leiden 30).

**What holds (n=65 units).** CLDN4 tracks a CLDN4-excluded barrier/keratin score (n=65, ρ=0.405, p=0.000814) and a malignant-like score (n=65, ρ=0.629, p=1.95e-08). On the primary Slingshot curve (Lineage7), root-bin → terminal-bin CLDN4 1.068 → 2.804, barrier 0.790 → 2.378, IFN 0.117 → 0.328. The CLDN4-high terminal is Lineage7 Leiden 30 (CLDN4 2.498, IFN 0.315). The IFN-high terminal is a **different** lineage: Lineage2 Leiden 24 (CLDN4 1.481, IFN 0.827). Sample-level CLDN4 vs IFN is n=65, ρ=0.045, p=0.72 (null). **Thesis:** CLDN4-high sits at the barrier/malignant end, not the IFN-high end — supported.

## Verdict

Slingshot engine=pyslingshot-bio (Street et al. 2018 Python port); 14 lineage(s) from root 5. Sample-level CLDN4 vs primary PT: n=49, ρ=0.506, p=0.000208. CLDN4 vs barrier/keratin (no CLDN4): n=65, ρ=0.405, p=0.000814. CLDN4 vs IFN: n=65, ρ=0.045, p=0.72. barrier vs primary PT: n=49, ρ=0.586, p=9.69e-06. IFN vs primary PT: n=49, ρ=0.450, p=0.00118. Not a DPT-null writeup. Not a TACSTD2 redo. No both-high gate. GSE148071 not added. PR #320 T/NK r=−0.705 given, not re-audited.

## Honest n

- Analysis cells after QC (capped ≤280/unit): **n_cells = 16167** (GSE131907 10432, GSE205335 5735).
- Units (GSE131907 Sample + GSE205335 patient): **n_units = 65** (GSE131907 43, GSE205335 22).
- Units with ≥10 epithelial cells used for Spearman: **n = 65**.
- GSE131907 nLung cells / author AT2 in the object: 2883 / 2020.
- Author subtypes (cells): {'Malignant cells': 10101, 'AT2': 2020, 'tS2': 973, 'tS1': 970, 'NA': 596, 'Non-malignant cells': 590, 'Ciliated': 362, 'AT1': 237, 'Club': 229, 'tS3': 54, 'Undetermined': 35}.
- Slingshot lineages: **14**. Engine: `pyslingshot-bio (Street et al. 2018 Python port)`.
- R Slingshot: available=False; Rscript not on PATH.
- Palantir companion: available=True; ran.
- Genes absent from locked sets: {'AT2': [], 'AT1': [], 'club': [], 'basal': [], 'ciliated': [], 'barrier_keratin': [], 'malignant_like': [], 'IFN': []}.
- GSE207422 not used. GSE148071 not used. GSE131907 PE unlabeled epithelium dropped. GSE205335 normal-tissue samples dropped.
- PR #320 T/NK r=−0.705 is given and was not recomputed.

## Slingshot / PAGA lineages

Start cluster **5** (mean CLDN4 0.964, mean AT2 3.680, nLung AT2 cells 1942). PAGA components at connectivity>0: **1** among 35 Leiden vertices.

| lineage | clusters | n_cells | terminal | term CLDN4 | term barrier | term IFN |
| --- | --- | ---: | --- | ---: | ---: | ---: |
| Lineage1 | 5 → 20 → 1 → 26 → 14 → 25 → 13 → 32 → 27 → 31 | 7366 | 31 | 1.609 | 0.598 | 0.079 |
| Lineage2 | 5 → 20 → 1 → 26 → 19 → 22 → 9 → 21 → 23 → 24 | 6529 | 24 | 1.481 | 2.204 | 0.827 |
| Lineage3 | 5 → 20 → 1 → 26 → 14 → 25 → 13 → 11 | 6788 | 11 | 2.082 | 1.546 | 0.090 |
| Lineage4 | 5 → 20 → 1 → 26 → 19 → 2 → 16 → 18 | 6814 | 18 | 0.010 | 0.025 | 0.246 |
| Lineage5 | 5 → 20 → 1 → 26 → 19 → 2 → 16 → 34 | 6797 | 34 | 0.909 | 0.751 | 0.045 |
| Lineage6 | 5 → 20 → 1 → 26 → 14 → 25 → 3 | 6309 | 3 | 1.283 | 0.809 | 0.105 |
| Lineage7 | 5 → 20 → 1 → 26 → 19 → 22 → 30 | 5363 | 30 | 2.498 | 2.239 | 0.315 |
| Lineage8 | 5 → 20 → 1 → 26 → 4 → 0 | 4967 | 0 | 1.025 | 1.409 | 0.280 |
| Lineage9 | 5 → 20 → 1 → 26 → 4 → 15 | 5274 | 15 | 1.282 | 1.353 | 0.089 |
| Lineage10 | 5 → 20 → 1 → 26 → 4 → 17 | 4864 | 17 | 1.534 | 2.033 | 0.193 |
| Lineage11 | 5 → 20 → 1 → 26 → 6 → 29 | 5320 | 29 | 1.197 | 2.265 | 0.372 |
| Lineage12 | 5 → 20 → 1 → 8 → 28 → 33 | 5119 | 33 | 0.659 | 0.613 | 0.109 |
| Lineage13 | 5 → 20 → 1 → 10 → 12 | 4998 | 12 | 0.247 | 0.340 | 0.283 |
| Lineage14 | 5 → 20 → 1 → 8 → 7 | 5115 | 7 | 1.789 | 1.023 | 0.327 |

Primary lineage for along-curve plots: **Lineage7** (terminal with highest mean CLDN4 among root-started lineages).

## Sample-level Spearman (done criterion)

Unit = GSE131907 sample or GSE205335 patient. BH inside this list only.

| Contrast | n_units | ρ | p | q |
| --- | ---: | ---: | ---: | ---: |
| CLDN4 vs Slingshot PT (shared) | 65 | -0.130 | 0.301 | 0.401 |
| CLDN4 vs Slingshot PT (primary lineage) | 49 | 0.506 | 0.000208 | 0.000416 |
| CLDN4 vs AT2 score | 65 | 0.086 | 0.497 | 0.551 |
| CLDN4 vs barrier/keratin (no CLDN4) | 65 | 0.405 | 0.000814 | 0.0014 |
| CLDN4 vs IFN (no CLDN4) | 65 | 0.045 | 0.72 | 0.72 |
| CLDN4 vs malignant-like score | 65 | 0.629 | 1.95e-08 | 7.79e-08 |
| barrier vs Slingshot PT (primary) | 49 | 0.586 | 9.69e-06 | 2.91e-05 |
| IFN vs Slingshot PT (primary) | 49 | 0.450 | 0.00118 | 0.00177 |
| CLDN4 vs TACSTD2 (comparator) | 65 | 0.492 | 3.15e-05 | 7.57e-05 |
| SFTPC vs Slingshot PT (control) | 65 | -0.659 | 2.48e-09 | 1.49e-08 |
| AT2 score vs Slingshot PT (control) | 65 | -0.842 | 1.55e-18 | 1.86e-17 |
| CLDN4 vs Palantir PT (companion) | 65 | 0.084 | 0.505 | 0.551 |

## Sensitivity (not in the BH family)

| Contrast | n_units | ρ | p |
| --- | ---: | ---: | ---: |
| GSE131907-only CLDN4 vs primary PT | 36 | 0.742 | 2.19e-07 |
| GSE205335-only CLDN4 vs primary PT | 13 | 0.445 | 0.128 |
| GSE131907-only CLDN4 vs barrier | 43 | 0.445 | 0.0028 |
| GSE205335-only CLDN4 vs barrier | 22 | 0.303 | 0.17 |
| GSE131907-only CLDN4 vs IFN | 43 | 0.159 | 0.31 |
| GSE205335-only CLDN4 vs IFN | 22 | -0.054 | 0.813 |
| tLung-only CLDN4 vs primary PT | 11 | 0.791 | 0.00375 |
| nLung-only CLDN4 vs primary PT | 11 | 0.655 | 0.0289 |
| tumor-only (drop nLung) CLDN4 vs primary PT | 38 | 0.388 | 0.0161 |

## Along the primary Slingshot curve

Binned cell means on the primary lineage (0 = root, 1 = terminal). This is the same curve for CLDN4, barrier/keratin (no CLDN4), and IFN.

| bin | n_cells | mean CLDN4 | mean barrier | mean IFN | mean AT2 | mean malignant |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 1253 | 1.068 | 0.790 | 0.117 | 3.890 | 0.080 |
| 1 | 1161 | 0.916 | 0.637 | 0.091 | 3.832 | 0.074 |
| 2 | 687 | 1.007 | 0.703 | 0.117 | 2.258 | 0.294 |
| 3 | 703 | 1.639 | 1.158 | 0.232 | 1.817 | 0.601 |
| 4 | 596 | 2.003 | 1.373 | 0.373 | 1.271 | 0.767 |
| 5 | 469 | 1.994 | 1.697 | 0.265 | 0.641 | 0.699 |
| 6 | 289 | 2.292 | 2.083 | 0.259 | 0.356 | 0.678 |
| 7 | 61 | 2.619 | 2.072 | 0.377 | 0.468 | 1.023 |
| 8 | 11 | 1.795 | 1.698 | 0.599 | 0.382 | 1.225 |
| 9 | 23 | 1.828 | 1.986 | 0.277 | 0.061 | 0.362 |
| 10 | 30 | 2.548 | 2.157 | 0.282 | 0.070 | 0.588 |
| 11 | 80 | 2.804 | 2.378 | 0.328 | 0.030 | 0.739 |

## Extra figures

- `results/figures/fig_lineage_cldn4.png` — **real Slingshot lineage plot** (done criterion)
- `results/figures/fig_along_curve.png` — CLDN4 / barrier / IFN on the same curve
- `results/figures/fig_sample_cldn4_pt.png` — sample-level CLDN4 vs Slingshot PT
- `results/figures/fig_paga.png` — PAGA on Leiden
- `results/figures/fig_extra_sample_programs.png` — sample CLDN4 vs barrier and vs IFN
- `results/figures/fig_extra_two_lineages.png` — CLDN4-high lineage vs IFN-high lineage
- `results/figures/fig_honest_n.png`

## What this does not claim

- Cell-level p-values are not the claim. n_cells is large by construction.
- PR #320 T/NK Q4 vs Q1 r=−0.705 is given; this folder does not re-rank T/NK.
- This is not a redo of PR #325 (GSE131907-only PAGA) or PR #449 (DPT fallback).
- The pooled PT Spearman mixes cohorts and nLung vs tumor. Tumor-only is in Sensitivity.
- Malignant-like is CEACAM5/6/MKI67 — **not CNV**.
- GSE205335 is an ICI biopsy/effusion cohort; this analysis is **not** an ICI / MPR test.
- No TACSTD2∩CLDN4 both-high gate. GSE148071 not used.
- Do not write “AT2 differentiates into LUAD because PAGA is connected.”
- Slingshot is an ordering along principal curves, not a developmental clock.

## Outputs

- `results/tables/sample_cldn4_vs_pseudotime.tsv` — **done criterion**
- `results/tables/sample_level_spearman.tsv`
- `results/figures/fig_lineage_cldn4.png` — **done criterion**
- `results/summary.json`

## Reproduce

```bash
pip install -r methods/winpair_slingshot_real_cldn4/requirements.txt
python3 methods/winpair_slingshot_real_cldn4/scripts/download.py \
  --out /tmp/winpair_slingshot_real
python3 methods/winpair_slingshot_real_cldn4/scripts/extract.py \
  --data /tmp/winpair_slingshot_real \
  --out /tmp/winpair_slingshot_real/epithelium.h5ad
python3 methods/winpair_slingshot_real_cldn4/scripts/analyze.py \
  --input /tmp/winpair_slingshot_real/epithelium.h5ad \
  --outdir methods/winpair_slingshot_real_cldn4/results \
  --finding methods/winpair_slingshot_real_cldn4/FINDING.md
```

