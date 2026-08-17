# Finding — winning-pair GSE131907+GSE205335 epithelium, CLDN4-only trajectory

ADDITIVE. **CLDN4 only.** Winning pair from the CLDN4-first combinatorial search (PR #290: author %pos GSE131907+GSE205335, n=43, ρ=−0.479 vs T/NK). This folder does **not** redo GSE131907-only PAGA (PR #325). GSE207422 is not added. No TACSTD2∩CLDN4 dual-high gate.

Primary clock: **documented AT2-rooted diffusion pseudotime (scanpy DPT; Slingshot R missing)**. Inferential unit = **sample/patient**. Cell-level ρ is descriptive. Barrier/keratin score **excludes CLDN4**. Root is GSE131907 nLung author AT2, never CLDN4-high.

**What holds (n=65 units).** CLDN4 tracks a CLDN4-excluded barrier/keratin score (ρ=0.376, p=0.00204, q=0.00368) and a malignant-like score (ρ=0.584, p=3.3e-07). Within the same unit, CLDN4-high cells are more barrier/keratin (paired n=62, Δmed +0.340, p=7.96e-12). The AT2-rooted DPT control works: AT2 vs DPT ρ=−0.786; SFTPC vs DPT ρ=−0.673.

**What does not hold.** Pooled sample-level CLDN4 vs DPT is null (n=65, ρ=0.078, p=0.538). CLDN4 vs AT2 is null (ρ=0.086, p=0.498). Per-cohort DPT tests are also null (GSE131907 n=43 p=0.289; GSE205335 n=22 p=0.294). Tumor-only (drop nLung) CLDN4 vs DPT is null (n=54, ρ=−0.032). The PR #325 n=22 tLung+nLung DPT association is **not** recovered on this mixed winning-pair object. tLung-only remains small-n (n=11, ρ=0.700, p=0.0165).

## Verdict

Sample-level CLDN4 vs AT2-rooted DPT: n=65, ρ=0.078, p=0.538. CLDN4 vs AT2 score: n=65, ρ=0.086, p=0.498. CLDN4 vs barrier/keratin (CLDN4 excluded from the score): n=65, ρ=0.376, p=0.00204. CLDN4 vs malignant-like: n=65, ρ=0.584, p=3.3e-07. CLDN4 vs TACSTD2 (comparator only): n=65, ρ=0.507, p=1.62e-05. GSE131907-only CLDN4 vs DPT: n=43, ρ=0.165, p=0.289. GSE205335-only CLDN4 vs DPT: n=22, ρ=0.234, p=0.294. Paired CLDN4-high vs low barrier/keratin: n=62, W=1.0, Δmed=0.340, p=7.96e-12. Paired CLDN4-high vs low AT2: n=62, W=946.0, Δmed=0.017, p=0.831. PAGA has 1 component(s) at connectivity>0 among 39 Leiden vertices. The pooled DPT correlation mixes cohorts and nLung vs tumor and is not a within-tumor progression test. Slingshot R was not run; DPT is the documented clock. Not a TACSTD2 redo. No both-high gate. GSE207422 not added.

## Honest n

- Analysis cells after QC (capped ≤350/unit): **n_cells = 18972** (GSE131907 11986, GSE205335 6986).
- Units (GSE131907 Sample + GSE205335 patient): **n_units = 65** (GSE131907 43, GSE205335 22).
- Units with ≥10 epithelial cells used for Spearman: **n = 65**.
- GSE131907 nLung cells / author AT2 in the object: 3141 / 2020.
- Author subtypes (cells): {'Malignant cells': 12125, 'AT2': 2020, 'tS2': 1171, 'tS1': 1103, 'Non-malignant cells': 711, 'NA': 665, 'Ciliated': 463, 'AT1': 310, 'Club': 304, 'tS3': 56, 'Undetermined': 44}.
- CLDN4 tertile cells: low 6324, mid 6324, high 6324.
- Units with ≥8 cells in both CLDN4-high and CLDN4-low arms: **n = 62**.
- Genes absent from locked sets: {'AT2': [], 'AT1': [], 'club': [], 'basal': [], 'ciliated': [], 'barrier_keratin': [], 'malignant_like': []}.
- GSE207422 not used. GSE131907 PE unlabeled epithelium dropped. GSE205335 normal-tissue samples dropped.
- Slingshot: available=False; Rscript not on PATH.
- Palantir: available=False; python package palantir not installed.

## Locked choices

- Leiden resolution 0.6; HVG 3000; neighbors 30; PCs 30.
- Batch: harmonypy on PCA, batch=dataset.
- DPT root: GSE131907 nLung author AT2 (median AT2 score) (root cell index 3038, unit GSE131907:LUNG_N34).
- PAGA components at connectivity>0: **1** among 39 Leiden vertices.
- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).

## Primary (sample-level Spearman, BH inside this list)

| Contrast | n_units | ρ | p | q |
| --- | ---: | ---: | ---: | ---: |
| CLDN4 vs DPT | 65 | 0.078 | 0.538 | 0.538 |
| CLDN4 vs AT2 score | 65 | 0.086 | 0.498 | 0.538 |
| CLDN4 vs club score | 65 | 0.260 | 0.0366 | 0.0471 |
| CLDN4 vs basal score | 65 | 0.266 | 0.0325 | 0.0471 |
| CLDN4 vs barrier/keratin (no CLDN4) | 65 | 0.376 | 0.00204 | 0.00368 |
| CLDN4 vs malignant-like score | 65 | 0.584 | 3.3e-07 | 9.91e-07 |
| CLDN4 vs TACSTD2 (comparator) | 65 | 0.507 | 1.62e-05 | 3.65e-05 |
| SFTPC vs DPT (control) | 65 | -0.673 | 7.88e-10 | 3.55e-09 |
| AT2 score vs DPT (control) | 65 | -0.786 | 8.71e-15 | 7.84e-14 |

## Sensitivity (not in the BH family)

| Contrast | n_units | ρ | p |
| --- | ---: | ---: | ---: |
| GSE131907-only CLDN4 vs DPT | 43 | 0.165 | 0.289 |
| GSE205335-only CLDN4 vs DPT | 22 | 0.234 | 0.294 |
| GSE131907-only CLDN4 vs AT2 | 43 | -0.179 | 0.252 |
| GSE205335-only CLDN4 vs AT2 | 22 | 0.102 | 0.651 |
| GSE131907-only CLDN4 vs barrier/keratin (no CLDN4) | 43 | 0.440 | 0.00314 |
| GSE205335-only CLDN4 vs barrier/keratin (no CLDN4) | 22 | 0.273 | 0.219 |
| tLung-only CLDN4 vs DPT | 11 | 0.700 | 0.0165 |
| nLung-only CLDN4 vs DPT | 11 | 0.509 | 0.11 |
| tumor-only (drop nLung) CLDN4 vs DPT | 54 | -0.032 | 0.817 |
| GSE205335 ADC+SQ CLDN4 vs DPT | 17 | 0.248 | 0.338 |

## Extra figure — CLDN4-high vs CLDN4-low (unit-paired)

Emitted: **True**. Rule: sample-level Spearman(CLDN4, barrier_keratin_no_CLDN4) ρ>0 and p<0.05, or any paired tertile Wilcoxon p<0.05, or n_paired≥4. Observed Spearman n=65, ρ=0.376, p=0.00204.

| Paired contrast (high − low) | n_units | Δ median | p |
| --- | ---: | ---: | ---: |
| barrier/keratin (no CLDN4) high vs low | 62 | 0.340 | 7.96e-12 |
| AT2 high vs low | 62 | 0.017 | 0.831 |
| malignant-like high vs low | 62 | 0.156 | 2e-09 |
| DPT high vs low | 62 | 0.009 | 5.86e-08 |
| club high vs low | 62 | 0.040 | 8.95e-05 |

## What this does not claim

- Cell-level p-values are not the claim. n_cells is large by construction.
- This is not a redo of PR #325 (GSE131907-only PAGA).
- The pooled DPT Spearman mixes cohorts, nLung and tumor, and is **not** a within-tumor progression test.
- Malignant-like is author tS1/tS2/tS3 / Malignant cells and/or CEACAM5/6/MKI67 — **not CNV**.
- GSE205335 is an ICI biopsy/effusion cohort, but this analysis is **not** an ICI / MPR / RECIST test.
- No TACSTD2∩CLDN4 both-high gate.
- Do not write “AT2 differentiates into LUAD because PAGA is connected.”
- Do not write “CLDN4 marks the malignant terminal.”
- Slingshot R was not run. DPT is an ordering, not a clock.

## Outputs

- `results/tables/sample_level_spearman.tsv` — **done criterion**
- `results/tables/sample_means.tsv`
- `results/figures/fig_trajectory_cldn4.png`
- `results/figures/fig_extra_cldn4_tertile.png`
- `results/figures/fig_honest_n.png`
- `results/summary.json`

## Reproduce

```bash
pip install -r methods/winpair_131907_205335_slingshot_cldn4/requirements.txt
python3 methods/winpair_131907_205335_slingshot_cldn4/scripts/download.py \
  --out /tmp/winpair_131907_205335
python3 methods/winpair_131907_205335_slingshot_cldn4/scripts/extract_epithelium.py \
  --data /tmp/winpair_131907_205335 \
  --out /tmp/winpair_131907_205335/epithelium.h5ad
python3 methods/winpair_131907_205335_slingshot_cldn4/scripts/analyze.py \
  --input /tmp/winpair_131907_205335/epithelium.h5ad \
  --outdir methods/winpair_131907_205335_slingshot_cldn4/results \
  --finding methods/winpair_131907_205335_slingshot_cldn4/FINDING.md
```

