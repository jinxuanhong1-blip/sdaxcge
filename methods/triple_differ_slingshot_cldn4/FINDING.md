# Finding — triple that differs (GSE123902+GSE131907+GSE205335), CLDN4-only trajectory

ADDITIVE. **CLDN4 only.** The triple that **differs** is taken as given from PR #459 (author %pos GSE123902+GSE131907+GSE205335, n=56, ρ=−0.522 vs T/NK). This folder does **not** redo that T/NK Spearman. Do **not** add GSE148071. This is **not** the 7-pool. No TACSTD2∩CLDN4 dual-high gate.

Primary clock: **Slingshot (Street 2018)**. Graph mode: **harmony_joint**. Inferential unit = **sample/patient/donor-sample**. Cell-level ρ is descriptive. Barrier/keratin score **excludes CLDN4**. Root is never CLDN4-high (prefer GSE131907 nLung author AT2).

**What holds (n=82 units).** CLDN4 vs CLDN4-excluded barrier/keratin: n=82, ρ=0.518, p=6.41e-07. CLDN4 vs malignant-like: n=82, ρ=0.601, p=2.32e-09. Within-unit CLDN4-high vs low barrier/keratin: n=77, W=1.0, Δmed=0.423, p=2.56e-14. **What does not automatically hold.** Stacked CLDN4 vs slingshot: n=82, ρ=-0.125, p=0.261. CLDN4 vs AT2: n=82, ρ=0.159, p=0.154.

## Verdict

Stacked sample-level CLDN4 vs AT2-rooted slingshot: n=82, ρ=-0.125, p=0.261. CLDN4 vs AT2 score: n=82, ρ=0.159, p=0.154. CLDN4 vs barrier/keratin (CLDN4 excluded from the score): n=82, ρ=0.518, p=6.41e-07. CLDN4 vs malignant-like: n=82, ρ=0.601, p=2.32e-09. CLDN4 vs TACSTD2 (comparator only): n=82, ρ=0.633, p=1.82e-10. GSE123902-only CLDN4 vs PT: n=17, ρ=-0.275, p=0.286. GSE131907-only CLDN4 vs PT: n=43, ρ=0.138, p=0.378. GSE205335-only CLDN4 vs PT: n=22, ρ=0.091, p=0.687. Paired CLDN4-high vs low barrier/keratin: n=77, W=1.0, Δmed=0.423, p=2.56e-14. Paired CLDN4-high vs low AT2: n=77, W=1097.0, Δmed=0.063, p=0.0581. PAGA has 1 component(s) at connectivity>0 among 32 Leiden vertices. Graph mode=harmony_joint. The stacked pseudotime correlation mixes cohorts and normal vs tumor and is not a within-tumor progression test. Not a TACSTD2 redo. No both-high gate. GSE148071 not added. Not the 7-pool.

## Honest n

- Analysis cells after QC (capped ≤350/unit): **n_cells = 23724** (GSE123902 4759, GSE131907 11986, GSE205335 6979).
- Units (GSE123902 donor-sample + GSE131907 Sample + GSE205335 patient): **n_units = 82** (GSE123902 17, GSE131907 43, GSE205335 22).
- Units with ≥10 epithelial cells used for Spearman: **n = 82**.
- GSE131907 nLung cells / author AT2 in the object: 3141 / 2020.
- GSE123902 NORMAL cells: 1343.
- Author/marker subtypes (cells): {'Malignant cells': 12118, 'marker_epithelial': 3613, 'AT2': 2020, 'tS2': 1171, 'AT2_like': 1146, 'tS1': 1103, 'Non-malignant cells': 711, 'NA': 665, 'Ciliated': 463, 'AT1': 310, 'Club': 304, 'tS3': 56, 'Undetermined': 44}.
- CLDN4 tertile cells: low 8471, mid 7345, high 7908.
- Units with ≥8 cells in both CLDN4-high and CLDN4-low arms: **n = 77**.
- Genes absent from locked sets: {'AT2': [], 'AT1': [], 'club': [], 'basal': [], 'ciliated': [], 'barrier_keratin': [], 'malignant_like': []}.
- GSE148071 not used. 7-pool extras not used. GSE131907 PE unlabeled epithelium dropped. GSE205335 normal-tissue samples dropped. GSE123902 36.5 GB H5 skipped; marker epithelium used.
- Graph: harmony_joint. Harmony: harmonypy on PCA, batch=dataset.
- Slingshot: available=True; [1] ‘2.10.0’. ran=True.
- Palantir: available=False; python package palantir not installed.

## Locked choices

- Leiden resolution 0.6; HVG 3000; neighbors 30; PCs 30.
- Batch: harmonypy on PCA, batch=dataset.
- Root: GSE131907 nLung author AT2 (median AT2; not CLDN4-high) (root cell index 5031, unit GSE131907:LUNG_N18, CLDN4 tertile mid).
- Slingshot fit: Street 2018 on a stratified subsample (2560 / 23724 cells, ≤80/Leiden, `approx_points=150`), then 5-NN projection to all cells. Start cluster = 7.
- PAGA components at connectivity>0: **1** among 32 Leiden vertices.
- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).

## Primary (stacked sample-level Spearman, BH inside this list)

| Contrast | n_units | ρ | p | q |
| --- | ---: | ---: | ---: | ---: |
| CLDN4 vs pseudotime | 82 | -0.125 | 0.261 | 0.261 |
| CLDN4 vs AT2 score | 82 | 0.159 | 0.154 | 0.174 |
| CLDN4 vs club score | 82 | 0.176 | 0.114 | 0.147 |
| CLDN4 vs basal score | 82 | 0.435 | 4.42e-05 | 6.63e-05 |
| CLDN4 vs barrier/keratin (no CLDN4) | 82 | 0.518 | 6.41e-07 | 1.44e-06 |
| CLDN4 vs malignant-like score | 82 | 0.601 | 2.32e-09 | 6.96e-09 |
| CLDN4 vs TACSTD2 (comparator) | 82 | 0.633 | 1.82e-10 | 8.17e-10 |
| SFTPC vs pseudotime (control) | 82 | -0.472 | 7.74e-06 | 1.39e-05 |
| AT2 score vs pseudotime (control) | 82 | -0.682 | 1.69e-12 | 1.52e-11 |

## Sensitivity (not in the BH family)

| Contrast | n_units | ρ | p |
| --- | ---: | ---: | ---: |
| GSE123902-only CLDN4 vs pseudotime | 17 | -0.275 | 0.286 |
| GSE131907-only CLDN4 vs pseudotime | 43 | 0.138 | 0.378 |
| GSE205335-only CLDN4 vs pseudotime | 22 | 0.091 | 0.687 |
| GSE123902-only CLDN4 vs AT2 | 17 | 0.559 | 0.0197 |
| GSE131907-only CLDN4 vs AT2 | 43 | -0.182 | 0.243 |
| GSE205335-only CLDN4 vs AT2 | 22 | 0.102 | 0.651 |
| GSE123902-only CLDN4 vs barrier/keratin (no CLDN4) | 17 | 0.789 | 0.000166 |
| GSE131907-only CLDN4 vs barrier/keratin (no CLDN4) | 43 | 0.440 | 0.00314 |
| GSE205335-only CLDN4 vs barrier/keratin (no CLDN4) | 22 | 0.273 | 0.219 |
| tLung/PRIMARY-only CLDN4 vs pseudotime | 19 | -0.325 | 0.175 |
| nLung/NORMAL-only CLDN4 vs pseudotime | 15 | -0.125 | 0.657 |
| tumor-only (drop nLung/NORMAL) CLDN4 vs pseudotime | 67 | -0.163 | 0.188 |
| GSE205335 ADC+SQ CLDN4 vs pseudotime | 17 | -0.069 | 0.794 |

## Extra figure — CLDN4-high vs CLDN4-low (unit-paired)

Emitted: **True**. Rule: sample-level Spearman(CLDN4, barrier_keratin_no_CLDN4) ρ>0 and p<0.05, or any paired tertile Wilcoxon p<0.05, or n_paired≥4. Observed Spearman n=82, ρ=0.518, p=6.41e-07.

| Paired contrast (high − low) | n_units | Δ median | p |
| --- | ---: | ---: | ---: |
| barrier/keratin (no CLDN4) high vs low | 77 | 0.423 | 2.56e-14 |
| AT2 high vs low | 77 | 0.063 | 0.0581 |
| malignant-like high vs low | 77 | 0.204 | 1.55e-12 |
| pseudotime high vs low | 77 | 6.222 | 2.96e-06 |
| club high vs low | 77 | 0.086 | 1.17e-07 |

## What this does not claim

- Cell-level p-values are not the claim. n_cells is large by construction.
- The PR #459 n=56 T/NK Spearman is given and is **not** re-audited here.
- The stacked pseudotime Spearman mixes cohorts and normal vs tumor and is **not** a within-tumor progression test.
- Malignant-like is author tS1/tS2/tS3 / Malignant cells and/or CEACAM5/6/MKI67 — **not CNV**.
- GSE123902 epithelium is marker-based (36.5 GB author H5 skipped), not author cell types.
- GSE205335 is an ICI biopsy/effusion cohort, but this analysis is **not** an ICI / MPR / RECIST test.
- No TACSTD2∩CLDN4 both-high gate.
- GSE148071 was not added. This is not the 7-pool.
- Do not write “AT2 differentiates into LUAD because PAGA is connected.”
- Do not write “CLDN4 marks the malignant terminal.”
- Slingshot ran (Street 2018, subsample + kNN projection). It is still an ordering along inferred lineages, not a developmental clock.

## Outputs

- `results/tables/stacked_sample_cldn4_pseudotime.tsv` — **done criterion**
- `results/tables/sample_level_spearman.tsv`
- `results/tables/sample_means.tsv`
- `results/figures/fig_trajectory_cldn4.png`
- `results/figures/fig_stacked_cldn4_pseudotime.png`
- `results/figures/fig_extra_cldn4_tertile.png`
- `results/figures/fig_honest_n.png`
- `results/summary.json`

## Reproduce

```bash
pip install -r methods/triple_differ_slingshot_cldn4/requirements.txt
python3 methods/triple_differ_slingshot_cldn4/scripts/download.py \
  --outdir /tmp/triple_differ_raw
python3 methods/triple_differ_slingshot_cldn4/scripts/extract_epithelium.py \
  --data /tmp/triple_differ_raw \
  --out /tmp/triple_differ_raw/epithelium.h5ad
python3 methods/triple_differ_slingshot_cldn4/scripts/analyze.py \
  --input /tmp/triple_differ_raw/epithelium.h5ad \
  --outdir methods/triple_differ_slingshot_cldn4/results \
  --finding methods/triple_differ_slingshot_cldn4/FINDING.md
```

