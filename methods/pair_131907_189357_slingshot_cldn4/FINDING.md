# Finding — pair GSE131907+GSE189357 epithelium, CLDN4-only REAL Slingshot/PAGA

ADDITIVE. **CLDN4 only.** Pair that already differs in PR #459 (malignant CLDN4 %pos vs T/NK, n=30, ρ=−0.542; Q4 r=−0.619). That T/NK rho is **not re-audited**. Does **not** redo GSE131907-only PAGA (PR #325) or the winning-pair Slingshot/DPT (PR #449). GSE148071 is not added. No TACSTD2∩CLDN4 dual-high gate.

Primary clock: **real Slingshot** (Street et al. 2018) on Harmony PCA + Leiden, start cluster = GSE131907 nLung author AT2, **never CLDN4-high**. PAGA is the graph geometry. Inferential unit = **GSE131907 sample + GSE189357 patient**. Cell-level ρ is descriptive. Barrier/keratin **excludes CLDN4**. IFN = Hallmark IFNα ∪ IFNγ mean (CLDN4 not a member).

**What holds (n=52 units).** CLDN4 vs Slingshot PT n=45, ρ=0.330, p=0.0268. Barrier vs Slingshot PT n=45, ρ=0.426, p=0.00351. IFN vs Slingshot PT n=45, ρ=0.118, p=0.441. CLDN4 vs barrier (CLDN4 excluded) n=52, ρ=0.518, p=8.5e-05. Within-unit CLDN4-high vs low barrier n=49, W=1.0, Δmed=0.323, p=7.11e-15; IFN n=49, W=25.0, Δmed=0.056, p=3.21e-12. **Honest split.** GSE131907-only n=36, ρ=0.360, p=0.031; GSE189357-only n=9, ρ=0.550, p=0.125 (GSE189357 n is the 9-patient arm — do not overclaim).

## Verdict

Sample-level CLDN4 vs AT2-rooted Slingshot PT: n=45, ρ=0.330, p=0.0268. Barrier/keratin (no CLDN4) vs Slingshot PT: n=45, ρ=0.426, p=0.00351. IFN (Hallmark α∪γ) vs Slingshot PT: n=45, ρ=0.118, p=0.441. CLDN4 vs barrier/keratin: n=52, ρ=0.518, p=8.5e-05. CLDN4 vs IFN: n=52, ρ=0.147, p=0.3. CLDN4 vs AT2: n=52, ρ=-0.149, p=0.291. GSE131907-only CLDN4 vs Slingshot PT: n=36, ρ=0.360, p=0.031. GSE189357-only CLDN4 vs Slingshot PT: n=9, ρ=0.550, p=0.125. Paired CLDN4-high vs low barrier: n=49, W=1.0, Δmed=0.323, p=7.11e-15. Paired CLDN4-high vs low IFN: n=49, W=25.0, Δmed=0.056, p=3.21e-12. Slingshot n_lineages=25 start=5. PAGA has 1 component(s) at connectivity>0 among 33 Leiden vertices. The pooled Slingshot correlation mixes cohorts and nLung vs tumor and is not a within-tumor progression test. Real Slingshot was run. Not a TACSTD2 redo. No both-high gate. GSE148071 not added.

## Honest n

- Analysis cells after QC (capped ≤350/unit): **n_cells = 15109** (GSE131907 11986, GSE189357 3123).
- Units (GSE131907 Sample + GSE189357 patient): **n_units = 52** (GSE131907 43, GSE189357 9). Do not write these as one patient n.
- Units with ≥10 epithelial cells used for Spearman: **n = 52**.
- Units with finite primary-lineage Slingshot PT (used for PT Spearman): **n = 45**. The other 7 units have cells, but none assigned to Lineage1.
- GSE131907 nLung cells / author AT2 in the object: 3141 / 2020.
- GSE189357 stages (cells): {'tumor': 8845, 'nLung': 3141, 'AIS': 1046, 'MIA': 1044, 'IAC': 1033}.
- Author / marker subtypes (cells): {'Malignant cells': 5850, 'marker_epithelium': 3123, 'AT2': 2020, 'tS2': 1171, 'tS1': 1103, 'NA': 665, 'Ciliated': 463, 'AT1': 310, 'Club': 304, 'tS3': 56, 'Undetermined': 44}.
- CLDN4 tertile cells: low 5037, mid 5036, high 5036.
- Units with ≥8 cells in both CLDN4-high and CLDN4-low arms: **n = 49**.
- Genes absent from locked sets: {'AT2': [], 'AT1': [], 'club': [], 'basal': [], 'ciliated': [], 'barrier_keratin': [], 'malignant_like': [], 'IFN_core': []}.
- IFN Hallmark union genes present: 220 / 224.
- GSE148071 not used. GSE131907 PE unlabeled epithelium dropped.
- Slingshot: available=True; version=2.10.0; n_lineages=25.
- Root: GSE131907 nLung author AT2 (median AT2 score); start ≠ CLDN4-high cluster (cell index 4011, unit GSE131907:LUNG_N09, tertile mid, start Leiden 5).

## Locked choices

- Leiden resolution 0.6; HVG 3000; neighbors 30; PCs 30.
- Batch: harmonypy on PCA, batch=dataset.
- Slingshot start: Leiden 5 (rejected CLDN4-high start=False).
- PAGA components at connectivity>0: **1** among 33 Leiden vertices.
- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).
- IFN: Hallmark IFNα ∪ IFNγ mean (**CLDN4 out**).

Slingshot returned **25 lineages** because it draws one path from start Leiden 5 to every MST leaf among 33 clusters. That is the real Street fit, not a curated two-path developmental model. Lineage1 is the primary clock (longest early path `5->9->2->12->6->8`).

## Lineage table (done criterion)

| lineage | start | end | n_clusters | n_cells | path | ρ CLDN4 | ρ barrier | ρ IFN |
| --- | --- | --- | ---: | ---: | --- | ---: | ---: | ---: |
| Lineage1 | 5 | 8 | 6 | 5970 | `5->9->2->12->6->8` | 0.305 | 0.307 | 0.103 |
| Lineage2 | 5 | 20 | 6 | 5392 | `5->9->31->3->10->20` | 0.157 | 0.479 | 0.196 |
| Lineage3 | 5 | 22 | 6 | 5161 | `5->9->31->3->10->22` | 0.072 | 0.410 | 0.123 |
| Lineage4 | 5 | 1 | 5 | 5166 | `5->9->2->12->1` | 0.066 | 0.359 | 0.064 |
| Lineage5 | 5 | 7 | 5 | 4809 | `5->9->31->3->7` | -0.087 | 0.280 | -0.114 |
| Lineage6 | 5 | 4 | 5 | 4534 | `5->9->31->3->4` | -0.276 | 0.217 | -0.078 |
| Lineage7 | 5 | 11 | 5 | 4887 | `5->9->2->12->11` | 0.184 | 0.311 | 0.180 |
| Lineage8 | 5 | 13 | 5 | 4709 | `5->9->31->3->13` | -0.103 | 0.302 | 0.024 |
| Lineage9 | 5 | 14 | 5 | 4795 | `5->9->31->3->14` | -0.387 | -0.081 | -0.226 |
| Lineage10 | 5 | 15 | 5 | 4826 | `5->9->31->3->15` | -0.180 | 0.289 | 0.047 |
| Lineage11 | 5 | 16 | 5 | 4822 | `5->9->31->3->16` | -0.123 | 0.281 | -0.102 |
| Lineage12 | 5 | 17 | 5 | 4805 | `5->9->31->3->17` | -0.048 | 0.342 | -0.095 |
| Lineage13 | 5 | 18 | 5 | 4690 | `5->9->31->3->18` | -0.189 | 0.290 | 0.044 |
| Lineage14 | 5 | 19 | 5 | 4628 | `5->9->31->3->19` | -0.274 | 0.137 | -0.138 |
| Lineage15 | 5 | 21 | 5 | 4786 | `5->9->31->3->21` | -0.151 | 0.314 | 0.020 |
| Lineage16 | 5 | 23 | 5 | 4796 | `5->9->31->3->23` | -0.049 | 0.321 | -0.048 |
| Lineage17 | 5 | 24 | 5 | 4785 | `5->9->31->3->24` | -0.156 | 0.261 | -0.299 |
| Lineage18 | 5 | 25 | 5 | 4798 | `5->9->31->3->25` | -0.178 | 0.301 | 0.033 |
| Lineage19 | 5 | 28 | 5 | 4636 | `5->9->31->3->28` | -0.180 | 0.275 | 0.021 |
| Lineage20 | 5 | 29 | 5 | 4814 | `5->9->31->3->29` | -0.111 | 0.201 | -0.269 |
| Lineage21 | 5 | 30 | 5 | 4652 | `5->9->31->3->30` | -0.178 | 0.252 | -0.041 |
| Lineage22 | 5 | 32 | 5 | 4794 | `5->9->31->3->32` | -0.156 | 0.260 | -0.081 |
| Lineage23 | 5 | 0 | 4 | 3965 | `5->9->2->0` | -0.060 | 0.042 | 0.009 |
| Lineage24 | 5 | 26 | 4 | 3994 | `5->9->2->26` | -0.215 | -0.053 | 0.158 |
| Lineage25 | 5 | 27 | 3 | 3655 | `5->9->27` | -0.387 | -0.241 | -0.021 |

## Primary (sample-level Spearman vs Slingshot PT, BH inside this list)

| Contrast | n_units | ρ | p | q |
| --- | ---: | ---: | ---: | ---: |
| CLDN4 vs Slingshot PT | 45 | 0.330 | 0.0268 | 0.0382 |
| barrier/keratin (no CLDN4) vs Slingshot PT | 45 | 0.426 | 0.00351 | 0.00585 |
| IFN (Hallmark α∪γ) vs Slingshot PT | 45 | 0.118 | 0.441 | 0.441 |
| CLDN4 vs barrier/keratin (no CLDN4) | 52 | 0.518 | 8.5e-05 | 0.00017 |
| CLDN4 vs IFN (Hallmark α∪γ) | 52 | 0.147 | 0.3 | 0.333 |
| CLDN4 vs AT2 score | 52 | -0.149 | 0.291 | 0.333 |
| CLDN4 vs malignant-like score | 52 | 0.643 | 2.7e-07 | 2.7e-06 |
| CLDN4 vs TACSTD2 (comparator) | 52 | 0.604 | 2.13e-06 | 1.06e-05 |
| SFTPC vs Slingshot PT (control) | 45 | -0.611 | 8.23e-06 | 2.74e-05 |
| AT2 score vs Slingshot PT (control) | 45 | -0.577 | 3.28e-05 | 8.2e-05 |

## Sensitivity (not in the BH family)

| Contrast | n_units | ρ | p |
| --- | ---: | ---: | ---: |
| GSE131907-only CLDN4 vs Slingshot PT | 36 | 0.360 | 0.031 |
| GSE189357-only CLDN4 vs Slingshot PT | 9 | 0.550 | 0.125 |
| GSE131907-only barrier vs Slingshot PT | 36 | 0.462 | 0.00456 |
| GSE189357-only barrier vs Slingshot PT | 9 | -0.217 | 0.576 |
| GSE131907-only IFN vs Slingshot PT | 36 | 0.089 | 0.607 |
| GSE189357-only IFN vs Slingshot PT | 9 | 0.517 | 0.154 |
| tLung-only CLDN4 vs Slingshot PT | 11 | 0.364 | 0.272 |
| nLung-only CLDN4 vs Slingshot PT | 11 | 0.564 | 0.071 |
| tumor-only (drop nLung) CLDN4 vs Slingshot PT | 34 | 0.286 | 0.102 |
| GSE189357 AIS/MIA/IAC CLDN4 vs Slingshot PT | 9 | 0.550 | 0.125 |
| CLDN4 vs DPT (companion, not the clock) | 52 | 0.118 | 0.405 |
| IFN core vs Slingshot PT | 45 | 0.111 | 0.469 |

## Extra figure — CLDN4-high vs CLDN4-low (unit-paired)

Emitted: **True**. Observed barrier Spearman n=52, ρ=0.5176299837787074, p=8.500507989316485e-05.

| Paired contrast (high − low) | n_units | Δ median | p |
| --- | ---: | ---: | ---: |
| barrier/keratin (no CLDN4) high vs low | 49 | 0.323 | 7.11e-15 |
| IFN (Hallmark a+g) high vs low | 49 | 0.056 | 3.21e-12 |
| AT2 high vs low | 49 | 0.034 | 0.767 |
| Slingshot PT high vs low | 38 | 10.699 | 1e-06 |

## What this does not claim

- Cell-level p-values are not the claim. n_cells is large by construction.
- This is not a redo of PR #325 (GSE131907-only PAGA) or PR #449 (winning-pair DPT).
- The pooled Slingshot Spearman mixes cohorts and nLung vs tumor and is **not** a within-tumor progression test.
- GSE189357 epithelium is **marker-gated**, not author-labeled. Malignant-like is CEACAM5/6/MKI67 — **not CNV**.
- PR #459 T/NK ρ is given and was not re-audited.
- No TACSTD2∩CLDN4 both-high gate. GSE148071 not added.
- Do not write “AT2 differentiates into LUAD because PAGA is connected.”
- Do not write “CLDN4 marks the malignant terminal.”
- Slingshot is an ordering, not a developmental clock.

## Outputs

- `results/tables/slingshot_lineages.tsv` — **done criterion** (also copied to `lineage_table.tsv`)
- `results/tables/sample_level_spearman.tsv`
- `results/figures/fig_trajectory_cldn4.png`
- `results/figures/fig_along_pseudotime.png`
- `results/figures/fig_extra_cldn4_tertile.png`
- `results/figures/fig_honest_n.png`
- `results/summary.json`

## Reproduce

```bash
pip install -r methods/pair_131907_189357_slingshot_cldn4/requirements.txt
# R 4.3 + Bioconductor slingshot (user library ~/R/library)
python3 methods/pair_131907_189357_slingshot_cldn4/scripts/run_all.py
```

