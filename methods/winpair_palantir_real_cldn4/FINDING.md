# Finding — winning-pair REAL Palantir, CLDN4 only (GSE131907+GSE205335)

ADDITIVE. **CLDN4 only.** Winning pair from the CLDN4-first combinatorial search (PR #290: author %pos GSE131907+GSE205335 vs T/NK). This folder does **not** redo GSE131907-only PAGA (PR #325) or the Slingshot/DPT fallback (PR #449). GSE148071 is not added. GSE207422 is not added. No TACSTD2∩CLDN4 dual-high gate.

Primary method: **Palantir** (Setty et al. 2019), installed and run. Inferential unit = **patient**. Cell-level ρ is descriptive. Barrier/keratin score **excludes CLDN4**. IFN core **excludes CLDN4**. Early cell is GSE131907 nLung author AT2 and is **never CLDN4-high**. DPT is a companion only; empty DPT is not a stop.

**Thesis tested.** CLDN4-high = barrier end; CLDN4-low = IFN-higher end of the same malignant trajectory. Palantir auto-terminals put the highest-barrier and highest-IFN means on the **same** destiny (destiny_1); branching into two destinies is not required for a two-end continuum. Thesis is **partial**: barrier end holds; IFN-higher CLDN4-low end is not supported at the pre-specified tests.

**What holds (n=48 patients).** CLDN4 vs barrier (CLDN4 excluded): n=48, ρ=0.307, p=0.0341. CLDN4 vs IFN: n=48, ρ=-0.083, p=0.574. Within-patient CLDN4-high vs low barrier n=47, W=1.0, Δmed=0.421, p=2.84e-14; IFN n=47, W=96.0, Δmed=0.054, p=4.98e-08. **Palantir product.** 3 destinies, 18972 cells with finite pseudotime. **What is mixed / null.** CLDN4 vs Palantir PT n=48, ρ=0.264, p=0.0697; IFN vs Palantir PT n=48, ρ=-0.148, p=0.314.

## Verdict

Patient-level CLDN4 vs Palantir PT: n=48, ρ=0.264, p=0.0697. Barrier (no CLDN4) vs Palantir PT: n=48, ρ=-0.196, p=0.182. IFN vs Palantir PT: n=48, ρ=-0.148, p=0.314. CLDN4 vs barrier (no CLDN4): n=48, ρ=0.307, p=0.0341. CLDN4 vs IFN: n=48, ρ=-0.083, p=0.574. IFN vs barrier: n=48, ρ=0.558, p=3.85e-05. Malignant-restricted CLDN4 vs IFN: n=46, ρ=-0.108, p=0.475. Malignant-restricted CLDN4 vs barrier: n=46, ρ=0.315, p=0.0332. Malignant-restricted IFN vs barrier: n=46, ρ=0.478, p=0.000791. Paired CLDN4-high vs low barrier: n=47, W=1.0, Δmed=0.421, p=2.84e-14. Paired CLDN4-high vs low IFN: n=47, W=96.0, Δmed=0.054, p=4.98e-08. Palantir destinies=3; early cell not CLDN4-high (mid). DPT companion empty=False. Not a TACSTD2 redo. No both-high gate. GSE148071 not added. GSE207422 not added.

## Honest n

- Analysis cells after QC (capped ≤350/sample): **n_cells = 18972** (GSE131907 11986, GSE205335 6986).
- Patients: **n_patients = 48** (GSE131907 26, GSE205335 22).
- Patients with ≥10 epithelial cells used for Spearman: **n = 48**.
- Patients with ≥10 malignant cells: **n_malignant_patients = 46**.
- GSE131907 nLung cells / author AT2 in the object: 3141 / 2020.
- Author subtypes (cells): {'Malignant cells': 12125, 'AT2': 2020, 'tS2': 1171, 'tS1': 1103, 'Non-malignant cells': 711, 'NA': 665, 'Ciliated': 463, 'AT1': 310, 'Club': 304, 'tS3': 56, 'Undetermined': 44}.
- CLDN4 tertile cells: low 6324, mid 6324, high 6324.
- Patients with ≥8 cells in both CLDN4-high and CLDN4-low arms: **n = 47**.
- Palantir destinies: **3** (destiny_1, destiny_2, destiny_3). Finite Palantir pseudotime cells: 18972.
- DPT companion: ran=True empty=False n_finite=18972; ok.
- Genes absent from locked sets: {'AT2': [], 'AT1': [], 'club': [], 'basal': [], 'ciliated': [], 'barrier_keratin': [], 'malignant_like': [], 'IFN': []}.
- GSE148071 not used. GSE207422 not used. GSE131907 PE unlabeled epithelium dropped. GSE205335 normal-tissue samples dropped.
- Palantir: available=True; version=1.4.5.

## Locked choices

- Leiden resolution 0.6; HVG 3000; neighbors 30; PCs 30.
- Batch: harmonypy on PCA, batch=dataset.
- Palantir early cell: GSE131907 nLung author AT2, not CLDN4-high, median AT2 score (cell `GSE131907:CTCGTCAAGCGATAGC_LUNG_N28`, patient GSE131907:28, CLDN4 tertile mid).
- Palantir waypoints: 500; diffusion components 10.
- Terminals: DM-boundary fallback after auto ARPACK failure (not defined by CLDN4 / barrier / IFN).
- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).
- IFN core: STAT1, IRF1/7/9, ISG15, IFIT1/2/3, OAS1/2, MX1/2, CXCL9/10/11, IDO1, TAP1, GBP1, IFI44L, IFI27, IFI44, RSAD2, USP18, EPSTI1, SAMD9 (**CLDN4 out**; LAMP3 and CDKN1A out).

## Destinies (cell assignment = argmax fate probability)

| Destiny | n_cells | n_patients | % malignant | mean CLDN4 | mean barrier (no CLDN4) | mean IFN | mean PT | top subtype |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| destiny_1 | 18283 | 48 | 0.79 | 1.294 | 1.238 | 0.222 | 0.087 | Malignant cells |
| destiny_2 | 342 | 2 | 1.00 | 0.966 | 0.773 | 0.050 | 0.811 | Malignant cells |
| destiny_3 | 347 | 1 | 1.00 | 1.681 | 0.612 | 0.080 | 0.913 | Malignant cells |

Auto terminal detection failed (`ArpackNoConvergence`). Destinies are Palantir absorption probabilities to three diffusion-map boundary cells (not chosen by CLDN4 / barrier / IFN). **destiny_1** holds 18,283/18,972 cells (all 48 patients). **destiny_2** (n=342, 2 patients) and **destiny_3** (n=347, 1 patient) are GSE205335-only terminals. They are not a second IFN arm. Highest-barrier and highest-IFN *means* both sit on destiny_1.

## Primary (patient-level Spearman, BH inside this list)

| Contrast | n_patients | ρ | p | q |
| --- | ---: | ---: | ---: | ---: |
| CLDN4 vs Palantir PT | 48 | 0.264 | 0.0697 | 0.296 |
| barrier (no CLDN4) vs Palantir PT | 48 | -0.196 | 0.182 | 0.344 |
| IFN vs Palantir PT | 48 | -0.148 | 0.314 | 0.411 |
| CLDN4 vs IFN | 48 | -0.083 | 0.574 | 0.61 |
| CLDN4 vs barrier (no CLDN4) | 48 | 0.307 | 0.0341 | 0.193 |
| IFN vs barrier (no CLDN4) | 48 | 0.558 | 3.85e-05 | 0.000654 |
| CLDN4 vs AT2 score | 48 | 0.086 | 0.563 | 0.61 |
| SFTPC vs Palantir PT (control) | 48 | -0.215 | 0.142 | 0.343 |
| CLDN4 vs fate destiny_1 | 48 | 0.235 | 0.108 | 0.343 |
| barrier vs fate destiny_1 | 48 | 0.181 | 0.217 | 0.369 |
| IFN vs fate destiny_1 | 48 | 0.017 | 0.909 | 0.909 |
| CLDN4 vs fate destiny_2 | 48 | -0.097 | 0.51 | 0.61 |
| barrier vs fate destiny_2 | 48 | -0.211 | 0.15 | 0.343 |
| IFN vs fate destiny_2 | 48 | -0.331 | 0.0217 | 0.185 |
| CLDN4 vs fate destiny_3 | 48 | 0.153 | 0.3 | 0.411 |
| barrier vs fate destiny_3 | 48 | -0.205 | 0.162 | 0.343 |
| IFN vs fate destiny_3 | 48 | -0.163 | 0.268 | 0.411 |

## Sensitivity (not in the BH family)

| Contrast | n_patients | ρ | p |
| --- | ---: | ---: | ---: |
| GSE131907-only CLDN4 vs Palantir PT | 26 | 0.086 | 0.674 |
| GSE205335-only CLDN4 vs Palantir PT | 22 | 0.448 | 0.0366 |
| GSE131907-only CLDN4 vs IFN | 26 | -0.020 | 0.922 |
| GSE205335-only CLDN4 vs IFN | 22 | -0.109 | 0.629 |
| GSE131907-only CLDN4 vs barrier | 26 | 0.328 | 0.102 |
| GSE205335-only CLDN4 vs barrier | 22 | 0.273 | 0.219 |
| malignant-cells CLDN4 vs Palantir PT | 46 | 0.139 | 0.358 |
| malignant-cells CLDN4 vs IFN | 46 | -0.108 | 0.475 |
| malignant-cells CLDN4 vs barrier | 46 | 0.315 | 0.0332 |
| malignant-cells IFN vs barrier | 46 | 0.478 | 0.000791 |
| malignant-cells IFN vs Palantir PT | 46 | -0.121 | 0.421 |
| malignant-cells barrier vs Palantir PT | 46 | -0.306 | 0.0387 |
| malignant-cells CLDN4 vs fate destiny_1 (barrier-high destiny) | 46 | 0.203 | 0.176 |
| malignant-cells IFN vs fate destiny_1 (IFN-high destiny) | 46 | 0.010 | 0.948 |

## Extra figure — CLDN4-high vs CLDN4-low (patient-paired)

Emitted: **True**. Observed patient Spearman(CLDN4, barrier) n=48, ρ=0.307, p=0.0341.

| Paired contrast (high − low) | n_patients | Δ median | p |
| --- | ---: | ---: | ---: |
| barrier/keratin (no CLDN4) high vs low | 47 | 0.421 | 2.84e-14 |
| IFN high vs low | 47 | 0.054 | 4.98e-08 |
| AT2 high vs low | 47 | 0.030 | 0.596 |
| Palantir PT high vs low | 47 | 0.009 | 1.98e-06 |
| fate destiny_1 high vs low | 47 | 0.000 | 0.18 |
| fate destiny_2 high vs low | 47 | 0.000 | 0.18 |
| fate destiny_3 high vs low | 47 | 0.000 | 0.317 |

## What this does not claim

- Cell-level p-values are not the claim. n_cells is large by construction.
- This is not a redo of PR #325 or PR #449.
- Palantir destinies are phenotypic terminals, not proven lineages.
- Barrier score excludes CLDN4; IFN core excludes CLDN4. Overlap of programs is still possible.
- Malignant-like is author tS1/tS2/tS3 / Malignant cells and/or CEACAM5/6/MKI67 — **not CNV**.
- GSE205335 is an ICI biopsy/effusion cohort, but this analysis is **not** an ICI / MPR / RECIST test.
- No TACSTD2∩CLDN4 both-high gate. GSE148071 not used.
- Do not write “AT2 differentiates into LUAD because Palantir ran.”
- DPT was not the product. Empty DPT would not have stopped this run.

## Outputs

- `results/tables/palantir_pseudotime.tsv` — **done criterion**
- `results/tables/palantir_destinies.tsv` — **done criterion**
- `results/tables/palantir_fate_probabilities.tsv`
- `results/tables/patient_means.tsv`
- `results/tables/patient_level_spearman.tsv`
- `results/figures/fig_palantir_destinies_programs.png`
- `results/figures/fig_extra_fate_vs_programs.png`
- `results/summary.json`

## Reproduce

```bash
pip install -r methods/winpair_palantir_real_cldn4/requirements.txt
python3 methods/winpair_palantir_real_cldn4/scripts/download.py \
  --out /tmp/winpair_131907_205335
python3 methods/winpair_palantir_real_cldn4/scripts/extract_epithelium.py \
  --data /tmp/winpair_131907_205335 \
  --out /tmp/winpair_131907_205335/epithelium.h5ad
python3 methods/winpair_palantir_real_cldn4/scripts/analyze.py \
  --input /tmp/winpair_131907_205335/epithelium.h5ad \
  --outdir methods/winpair_palantir_real_cldn4/results \
  --finding methods/winpair_palantir_real_cldn4/FINDING.md
```

