# Finding — pair GSE123902+GSE131907 epithelium, CLDN4-only REAL Slingshot/PAGA

ADDITIVE. **CLDN4 only.** Pair that **differs** in PR #459 (malignant CLDN4 %pos GSE123902+GSE131907, n=34, ρ=−0.575 vs T/NK). That T/NK rho is **not re-audited**. This folder asks a different question: where do **CLDN4**, a CLDN4-excluded **barrier/keratin** score, and a compact **IFN** score sit on a real Slingshot lineage. Does **not** redo GSE131907-only PAGA (PR #325) or the winning-pair GSE131907+GSE205335 Slingshot/DPT (PR #449). GSE148071 is not added. No TACSTD2∩CLDN4 dual-high gate.

Primary clock: **Slingshot 2.10.0** (16 lineage(s); primary `Lineage1`). PAGA is geometry only. DPT is a companion ordering, not the claim clock. Inferential unit = **GSE123902 donor + GSE131907 sample**. Cell-level ρ is descriptive. Barrier/keratin score **excludes CLDN4**. Root is GSE131907 nLung author AT2, never CLDN4-high.

**What holds (n=55 units).** CLDN4 vs Slingshot PT is n=55, ρ=0.089, p=0.517. Barrier vs Slingshot PT is n=55, ρ=0.520, p=4.79e-05. IFN vs Slingshot PT is n=55, ρ=0.314, p=0.0194. CLDN4 tracks barrier/keratin (n=56, ρ=0.471, p=0.000246) and IFN (n=56, ρ=-0.069, p=0.611). Within-unit CLDN4-high vs low barrier: n=50, W=1.0, Δmed=0.308, p=3.55e-15; IFN: n=50, W=24.0, Δmed=0.066, p=1.35e-12.

## Verdict

Sample-level CLDN4 vs AT2-rooted Slingshot PT: n=55, ρ=0.089, p=0.517. Barrier/keratin (no CLDN4) vs Slingshot PT: n=55, ρ=0.520, p=4.79e-05. IFN vs Slingshot PT: n=55, ρ=0.314, p=0.0194. CLDN4 vs barrier/keratin (CLDN4 excluded): n=56, ρ=0.471, p=0.000246. CLDN4 vs IFN: n=56, ρ=-0.069, p=0.611. CLDN4 vs malignant-like: n=56, ρ=0.108, p=0.428. GSE131907-only CLDN4 vs Slingshot: n=42, ρ=0.241, p=0.125. GSE123902-only CLDN4 vs Slingshot: n=13, ρ=-0.011, p=0.972. Paired CLDN4-high vs low barrier/keratin: n=50, W=1.0, Δmed=0.308, p=3.55e-15. Paired CLDN4-high vs low IFN: n=50, W=24.0, Δmed=0.066, p=1.35e-12. PAGA has 1 component(s) at connectivity>0 among 31 Leiden vertices. Slingshot produced 16 lineage(s); primary=Lineage1. The pooled Slingshot correlation mixes cohorts and nLung vs tumor and is not a within-tumor progression test. REAL Slingshot was run. Not a TACSTD2 redo. No both-high gate. GSE148071 not added.

## Honest n

- Analysis cells after QC (capped ≤350/unit): **n_cells = 14651** (GSE131907 11983, GSE123902 2668).
- Units (GSE131907 Sample + GSE123902 donor): **n_units = 56** (GSE131907 43, GSE123902 13).
- Units with ≥10 epithelial cells used for Spearman: **n = 56**.
- GSE131907 nLung cells / author AT2 in the object: 3138 / 2019.
- Author / marker subtypes (cells): {'Malignant cells': 5850, 'marker_epithelial': 2668, 'AT2': 2019, 'tS2': 1171, 'tS1': 1103, 'NA': 665, 'Ciliated': 463, 'AT1': 310, 'Club': 302, 'tS3': 56, 'Undetermined': 44}.
- CLDN4 tertile cells: low 4884, mid 4883, high 4884.
- Units with ≥8 cells in both CLDN4-high and CLDN4-low arms: **n = 50**.
- Cells on primary Slingshot lineage: **6928**.
- Genes absent from locked sets: {'AT2': ['SFTPC', 'SFTPA1'], 'AT1': [], 'club': ['SCGB1A1', 'SCGB3A2', 'SCGB3A1'], 'basal': ['KRT5', 'TP63', 'NGFR'], 'ciliated': [], 'barrier_keratin': [], 'malignant_like': ['CEACAM5', 'CEACAM6'], 'IFN': ['CXCL9', 'CXCL11', 'IDO1'], 'single_genes': ['SFTPC', 'SCGB1A1', 'KRT5']}.
- GSE148071 not used. GSE123902 NORMAL samples dropped. GSE131907 PE unlabeled epithelium dropped. Author 36.5 GB H5 skipped.
- Slingshot: available=True; version=2.10.0.
- Start cluster: 2 (mean CLDN4=1.0836125568474073); CLDN4-high clusters refused as root: ['16', '18', '21', '4', '9'].

## Locked choices

- Leiden resolution 0.6; HVG 3000; neighbors 30; PCs 30. Slingshot curves use the first 10 Harmony PCs and `approx_points=150`.
- Batch: harmonypy on PCA, batch=dataset.
- Shared-gene intersection dropped SFTPC / SCGB1A1 / CEACAM5/6 / CXCL9 / CXCL11 / IDO1 (GSE123902 dense CSV). AT2 score still uses SFTPB/NAPSA/LAMP3/ABCA3. SFTPC-vs-PT control is therefore NA.
- Slingshot / DPT root: GSE131907 nLung author AT2 (median AT2, not CLDN4-high Leiden) start_cluster=2 (root cell index 4718, unit GSE131907:LUNG_N09).
- PAGA components at connectivity>0: **1** among 31 Leiden vertices.
- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).
- IFN: compact ISG / IFNG-response core (STAT1, IRF1/7/9, ISG15, IFIT1-3, MX1/2, OAS*, CXCL9/10/11, IDO1, GBP*, IFI*, RSAD2, B2M, TAP1, PSMB8/9). **CLDN4 out**.

## Lineages (done criterion)

Primary lineage `Lineage1` starts at Leiden 2 (nLung AT2, not CLDN4-high). 16 Slingshot lineage(s).

| lineage | start | end | n_clusters | n_cells | path |
| --- | --- | --- | ---: | ---: | --- |
| Lineage1 | 2 | 4 | 9 | 6928 | `2>3>6>7>16>15>27>29>4` |
| Lineage2 | 2 | 9 | 9 | 6957 | `2>3>6>7>16>15>27>29>9` |
| Lineage3 | 2 | 30 | 9 | 6940 | `2>3>6>7>16>15>27>29>30` |
| Lineage4 | 2 | 14 | 7 | 7281 | `2>3>6>7>17>13>14` |
| Lineage5 | 2 | 19 | 7 | 6796 | `2>3>6>7>21>8>19` |
| Lineage6 | 2 | 20 | 7 | 6295 | `2>3>6>7>16>15>20` |
| Lineage7 | 2 | 22 | 7 | 7282 | `2>3>6>7>17>13>22` |
| Lineage8 | 2 | 28 | 7 | 6101 | `2>3>6>7>16>15>28` |
| Lineage9 | 2 | 12 | 6 | 5235 | `2>3>6>0>24>12` |
| Lineage10 | 2 | 25 | 6 | 5048 | `2>3>6>0>24>25` |
| Lineage11 | 2 | 26 | 6 | 5772 | `2>3>6>7>16>26` |
| Lineage12 | 2 | 5 | 5 | 5558 | `2>3>6>7>5` |
| Lineage13 | 2 | 10 | 5 | 5332 | `2>3>6>11>10` |
| Lineage14 | 2 | 18 | 5 | 5720 | `2>3>6>7>18` |
| Lineage15 | 2 | 23 | 4 | 4971 | `2>3>6>23` |
| Lineage16 | 2 | 1 | 3 | 3827 | `2>3>1` |

Machine table: `results/tables/slingshot_lineages.tsv`.

## Primary (sample-level Spearman, BH inside this list)

| Contrast | n_units | ρ | p | q |
| --- | ---: | ---: | ---: | ---: |
| CLDN4 vs Slingshot PT | 55 | 0.089 | 0.517 | 0.632 |
| barrier/keratin (no CLDN4) vs Slingshot PT | 55 | 0.520 | 4.79e-05 | 0.000176 |
| IFN vs Slingshot PT | 55 | 0.314 | 0.0194 | 0.0427 |
| CLDN4 vs barrier/keratin (no CLDN4) | 56 | 0.471 | 0.000246 | 0.000675 |
| CLDN4 vs IFN | 56 | -0.069 | 0.611 | 0.673 |
| CLDN4 vs AT2 score | 56 | 0.133 | 0.328 | 0.602 |
| CLDN4 vs malignant-like score | 56 | 0.108 | 0.428 | 0.632 |
| CLDN4 vs TACSTD2 (comparator) | 56 | 0.619 | 3.64e-07 | 2e-06 |
| SFTPC vs Slingshot PT (control) | 0 | NA | NA | NA |
| AT2 score vs Slingshot PT (control) | 55 | -0.633 | 2.19e-07 | 2e-06 |
| CLDN4 vs DPT (companion) | 56 | 0.089 | 0.514 | 0.632 |

## Sensitivity (not in the BH family)

| Contrast | n_units | ρ | p |
| --- | ---: | ---: | ---: |
| GSE131907-only CLDN4 vs Slingshot PT | 42 | 0.241 | 0.125 |
| GSE123902-only CLDN4 vs Slingshot PT | 13 | -0.011 | 0.972 |
| GSE131907-only barrier vs Slingshot PT | 42 | 0.623 | 1.07e-05 |
| GSE123902-only barrier vs Slingshot PT | 13 | 0.242 | 0.426 |
| GSE131907-only IFN vs Slingshot PT | 42 | 0.338 | 0.0285 |
| GSE123902-only IFN vs Slingshot PT | 13 | -0.165 | 0.59 |
| GSE131907-only CLDN4 vs IFN | 43 | 0.083 | 0.596 |
| GSE123902-only CLDN4 vs IFN | 13 | -0.352 | 0.239 |
| tLung-only CLDN4 vs Slingshot PT | 11 | 0.745 | 0.00845 |
| nLung-only CLDN4 vs Slingshot PT | 11 | 0.464 | 0.151 |
| tumor-only (drop nLung) CLDN4 vs Slingshot PT | 44 | 0.067 | 0.665 |
| GSE131907-only CLDN4 vs DPT (companion) | 43 | 0.172 | 0.271 |
| GSE123902-only CLDN4 vs DPT (companion) | 13 | 0.005 | 0.986 |

## Extra figure — CLDN4-high vs CLDN4-low (unit-paired)

Emitted: **True**. Rule: sample-level Spearman(CLDN4, barrier_keratin_no_CLDN4) ρ>0 and p<0.05, or any paired tertile Wilcoxon p<0.05, or n_paired≥4. Observed Spearman n=56, ρ=0.471, p=0.000246.

| Paired contrast (high − low) | n_units | Δ median | p |
| --- | ---: | ---: | ---: |
| barrier/keratin (no CLDN4) high vs low | 50 | 0.308 | 3.55e-15 |
| IFN high vs low | 50 | 0.066 | 1.35e-12 |
| AT2 high vs low | 50 | 0.061 | 0.605 |
| malignant-like high vs low | 50 | 0.006 | 0.00683 |
| Slingshot PT high vs low | 47 | -0.597 | 0.686 |
| DPT high vs low | 50 | 0.022 | 1.05e-05 |

## What this does not claim

- Cell-level p-values are not the claim. n_cells is large by construction.
- Do not write n_cells as the inferential n. GSE123902 is **donor**; GSE131907 is **sample**.
- This is not a redo of PR #325 (GSE131907-only PAGA) or PR #449 (GSE131907+GSE205335 DPT).
- The pooled Slingshot Spearman mixes cohorts and nLung vs tumor and is **not** a within-tumor progression test.
- GSE123902 epithelium is a marker gate (EPCAM/KRT+, PTPRC−), not the skipped 36.5 GB author H5.
- Malignant-like is author tS1/tS2/tS3 / Malignant cells and/or CEACAM5/6/MKI67 — **not CNV**.
- PR #459 T/NK ρ is given and was not re-audited.
- No TACSTD2∩CLDN4 both-high gate. GSE148071 not added.
- Do not write “AT2 differentiates into LUAD because PAGA is connected.”
- Do not write “CLDN4 marks the malignant terminal.”
- Slingshot is an ordering of Harmony space, not a developmental clock.

## Outputs

- `results/tables/slingshot_lineages.tsv` — **done criterion (lineage)**
- `results/tables/sample_level_spearman.tsv` — **done criterion (sample-level)**
- `results/tables/sample_means.tsv`
- `results/figures/fig_trajectory_cldn4.png`
- `results/figures/fig_extra_along_pseudotime.png`
- `results/figures/fig_extra_cldn4_tertile.png`
- `results/figures/fig_honest_n.png`
- `results/summary.json`

## Reproduce

```bash
bash methods/pair_123902_131907_slingshot_cldn4/scripts/install_r_slingshot.sh
pip install -r methods/pair_123902_131907_slingshot_cldn4/requirements.txt
python3 methods/pair_123902_131907_slingshot_cldn4/scripts/download.py \
  --out /tmp/pair_123902_131907
python3 methods/pair_123902_131907_slingshot_cldn4/scripts/extract_epithelium.py \
  --data /tmp/pair_123902_131907 \
  --out /tmp/pair_123902_131907/epithelium.h5ad
python3 methods/pair_123902_131907_slingshot_cldn4/scripts/analyze.py \
  --input /tmp/pair_123902_131907/epithelium.h5ad \
  --outdir methods/pair_123902_131907_slingshot_cldn4/results \
  --finding methods/pair_123902_131907_slingshot_cldn4/FINDING.md
```

