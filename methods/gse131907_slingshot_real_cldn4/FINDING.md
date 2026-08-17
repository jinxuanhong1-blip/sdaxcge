# Finding — REAL Slingshot on GSE131907 epithelium/malignant, CLDN4-only

ADDITIVE. **CLDN4 only.** Kim et al., *Nat Commun* 2020, PMID 32385277 (GSE131907). Epithelium / malignant cells only. This folder does **not** redo PR #325 (DPT/PAGA on nLung+tLung). It does **not** redo PR #449 (winning-pair DPT fallback; Slingshot R missing). GSE207422 and GSE205335 are not added. No TACSTD2∩CLDN4 dual-high gate.

Primary clock: **Slingshot (Street 2018 Python equivalent: MST + principal curves)**. Inferential unit = **sample**. Cell-level ρ is descriptive. Barrier/keratin and IFN scores **exclude CLDN4**. Root is nLung author AT2, never CLDN4-high.

**What holds (sample n=43).** CLDN4 vs Slingshot n=43, ρ=0.419, p=0.00519. CLDN4 vs barrier/keratin (CLDN4 excluded) n=43, ρ=0.435, p=0.00355. CLDN4 vs IFN n=43, ρ=0.179, p=0.251. Within-sample CLDN4-high vs low barrier n=40, W=0.0, Δmed=0.257, p=1.82e-12; IFN n=40, W=3.0, Δmed=0.065, p=9.09e-12.

## Verdict

Sample-level CLDN4 vs AT2-rooted Slingshot: n=43, ρ=0.419, p=0.00519. CLDN4 vs AT2 score: n=43, ρ=-0.185, p=0.235. CLDN4 vs barrier/keratin (CLDN4 excluded from the score): n=43, ρ=0.435, p=0.00355. CLDN4 vs IFN (CLDN4 excluded): n=43, ρ=0.179, p=0.251. CLDN4 vs malignant-like: n=43, ρ=0.637, p=4.44e-06. tLung-only CLDN4 vs Slingshot: n=11, ρ=0.791, p=0.00375. Paired CLDN4-high vs low barrier/keratin: n=40, W=0.0, Δmed=0.257, p=1.82e-12. Paired CLDN4-high vs low IFN: n=40, W=3.0, Δmed=0.065, p=9.09e-12. Paired CLDN4-high vs low AT2: n=40, W=218.0, Δmed=-0.042, p=0.009. Slingshot has 9 lineage(s) from Leiden 2. PAGA has 1 component(s) at connectivity>0 among 31 Leiden vertices. The mixed nLung+tLung+met Slingshot correlation is not a within-tumor progression test. Not ICI. Not a TACSTD2 redo. No both-high gate. Not a redo of PR #325 DPT.

## Honest n

- Catalog epithelium/malignant (before cap): **n_cells_catalog = 36071** across **n_samples_catalog = 43** (origins mBrain 15463, tLung 7270, tL/B 6582, nLung 3703, mLN 3053).
- Analysis cells after QC (cap ≤400/sample, nLung AT2 protected): **n_cells = 13031** (nLung 3305, tLung 3120, other malignant sites 6606).
- Samples: **n_samples = 43** (nLung 11, tLung 11, other 21).
- Patients: **n_patients = 26**.
- Samples with ≥10 cells used for Spearman: **n = 43**.
- Author AT2 in the object: **2020**. Author basal cells: **0**.
- Author subtypes (cells): {'Malignant cells': 6466, 'AT2': 2020, 'tS2': 1309, 'tS1': 1188, 'NA': 704, 'Ciliated': 531, 'AT1': 365, 'Club': 338, 'tS3': 59, 'Undetermined': 51}.
- CLDN4 tertile cells: low 4344, mid 4343, high 4344.
- Samples with ≥8 cells in both CLDN4-high and CLDN4-low arms: **n = 40**.
- Genes absent from locked sets: {'AT2': [], 'AT1': [], 'club': [], 'basal': [], 'ciliated': [], 'barrier_keratin': [], 'malignant_like': [], 'IFN': []}.
- Slingshot lineages: **9** from start cluster 2 (root sample LUNG_N20).
- Slingshot engine: python:slingshot_py.
- PAGA components at connectivity>0: **1** among 31 Leiden vertices.
- GSE207422 / GSE205335 not used. PE unlabeled epithelium dropped.

## Locked choices

- Leiden resolution 0.6; HVG 3000; neighbors 30; PCs 30.
- Slingshot embedding: first 5 PCs (UMAP is visualization only).
- Root: nLung author AT2 (median AT2 score); never CLDN4-high (root cell index 2871, sample LUNG_N20, subtype AT2, CLDN4 tertile mid).
- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).
- IFN genes: compact ISG panel in `scripts/gene_sets.py` (**CLDN4 out**).

## Primary (sample-level Spearman, BH inside this list)

| Contrast | n_samples | ρ | p | q |
| --- | ---: | ---: | ---: | ---: |
| CLDN4 vs Slingshot | 43 | 0.419 | 0.00519 | 0.00693 |
| CLDN4 vs AT2 score | 43 | -0.185 | 0.235 | 0.274 |
| CLDN4 vs club score | 43 | 0.111 | 0.478 | 0.478 |
| CLDN4 vs basal score | 43 | 0.484 | 0.00102 | 0.00175 |
| CLDN4 vs barrier/keratin (no CLDN4) | 43 | 0.435 | 0.00355 | 0.00532 |
| CLDN4 vs IFN (no CLDN4) | 43 | 0.179 | 0.251 | 0.274 |
| CLDN4 vs malignant-like score | 43 | 0.637 | 4.44e-06 | 1.78e-05 |
| CLDN4 vs TACSTD2 (comparator) | 43 | 0.507 | 0.00052 | 0.00104 |
| SFTPC vs Slingshot (control) | 43 | -0.653 | 2.07e-06 | 1.24e-05 |
| AT2 score vs Slingshot (control) | 43 | -0.680 | 5.29e-07 | 6.35e-06 |
| IFN vs Slingshot | 43 | 0.542 | 0.000172 | 0.000473 |
| barrier/keratin vs Slingshot | 43 | 0.538 | 0.000197 | 0.000473 |

## Sensitivity (not in the BH family)

| Contrast | n_samples | ρ | p |
| --- | ---: | ---: | ---: |
| tLung-only CLDN4 vs Slingshot | 11 | 0.791 | 0.00375 |
| nLung-only CLDN4 vs Slingshot | 11 | 0.573 | 0.0655 |
| tumor/met (drop nLung) CLDN4 vs Slingshot | 32 | 0.315 | 0.0788 |
| nLung+tLung only CLDN4 vs Slingshot | 22 | 0.660 | 0.000829 |
| tLung-only CLDN4 vs barrier/keratin (no CLDN4) | 11 | 0.645 | 0.032 |
| tLung-only CLDN4 vs IFN | 11 | 0.818 | 0.00208 |
| tLung-only CLDN4 vs AT2 | 11 | -0.882 | 0.00033 |
| CLDN4 vs DPT (sensitivity clock) | 43 | 0.258 | 0.0949 |
| nLung+tLung CLDN4 vs DPT (PR #325-like slice) | 22 | 0.287 | 0.195 |

## Extra figure — CLDN4-high vs CLDN4-low (sample-paired)

Emitted: **True**. Rule: sample-level Spearman(CLDN4, barrier_keratin_no_CLDN4) ρ>0 and p<0.05, or any paired tertile Wilcoxon p<0.05, or n_paired≥4. Observed Spearman n=43, ρ=0.435, p=0.00355.

| Paired contrast (high − low) | n_samples | Δ median | p |
| --- | ---: | ---: | ---: |
| barrier/keratin (no CLDN4) high vs low | 40 | 0.257 | 1.82e-12 |
| IFN (no CLDN4) high vs low | 40 | 0.065 | 9.09e-12 |
| AT2 high vs low | 40 | -0.042 | 0.009 |
| malignant-like high vs low | 40 | 0.129 | 1.64e-09 |
| Slingshot high vs low | 40 | 0.160 | 1.82e-12 |
| DPT high vs low | 40 | 0.013 | 5.03e-09 |

## What this does not claim

- Cell-level p-values are not the claim. n_cells is large by construction.
- This is not a redo of PR #325 (GSE131907 DPT/PAGA).
- The pooled Slingshot Spearman mixes nLung, tLung, and metastatic sites and is **not** a within-tumor progression test.
- Malignant-like is author tS1/tS2/tS3 / Malignant cells and/or CEACAM5/6/MKI67 — **not CNV**.
- GSE131907 is treatment-naive. Do not write ICI / MPR / RECIST language.
- No TACSTD2∩CLDN4 both-high gate.
- Do not write “AT2 differentiates into LUAD because Slingshot/PAGA is connected.”
- Do not write “CLDN4 marks the malignant terminal.”
- Slingshot is an ordering on an embedding, not a developmental clock.

## Outputs

- `results/tables/sample_level_spearman.tsv` — **done criterion**
- `results/tables/sample_means.tsv`
- `results/tables/slingshot_lineages.tsv`
- `results/tables/lineage_bin_means.tsv` — CLDN4 / barrier / IFN along lineages
- `results/tables/leiden_paga_vertices.tsv`
- `results/figures/fig_trajectory_cldn4.png`
- `results/figures/fig_extra_cldn4_tertile.png`
- `results/figures/fig_extra_lineage_programs.png`
- `results/figures/fig_extra_sample_cldn4_programs.png`
- `results/figures/fig_honest_n.png`
- `results/summary.json`

## Reproduce

```bash
pip install -r methods/gse131907_slingshot_real_cldn4/requirements.txt
python3 methods/gse131907_slingshot_real_cldn4/scripts/download.py \
  --out /tmp/gse131907_slingshot_data
python3 methods/gse131907_slingshot_real_cldn4/scripts/extract_epithelium.py \
  --data /tmp/gse131907_slingshot_data \
  --out /tmp/gse131907_slingshot_data/epithelium_malignant.h5ad
python3 methods/gse131907_slingshot_real_cldn4/scripts/analyze.py \
  --input /tmp/gse131907_slingshot_data/epithelium_malignant.h5ad \
  --outdir methods/gse131907_slingshot_real_cldn4/results \
  --finding methods/gse131907_slingshot_real_cldn4/FINDING.md
```

