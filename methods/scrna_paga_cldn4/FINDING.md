# Finding — PAGA / diffusion on LUAD epithelium scored by CLDN4

Additive public slice. **Not a TACSTD2 redo.** Primary readout is **CLDN4** on author-labeled GSE131907 LUAD epithelium (Kim et al., *Nat Commun* 2020, PMID 32385277). The public raw UMI matrix is the same file already referenced by `methods/scrna_paga/`. GSE207422 was not pooled (NSCLC mixed histology; no author epithelial labels on GEO).

PAGA + AT2-rooted diffusion pseudotime. Inferential unit = **sample**. Cell-level ρ is descriptive. Barrier/keratin score **excludes CLDN4** (no circularity). Root is nLung AT2, never CLDN4-high.

**What holds (sample n=22, patients n=12).** CLDN4 rises with AT2-rooted DPT (ρ=0.618, p=0.00216) and with a CLDN4-excluded barrier/keratin score (ρ=0.775, p=2.26e-05). It is anti-correlated with AT2 (ρ=−0.833). Within the same sample, CLDN4-high cells are more barrier/keratin (Δmed +0.304, p=4.77e-07), later on DPT (Δmed +0.148, p=9.54e-07), and less AT2 (Δmed −1.003, p=2.38e-06). Club is a null (ρ=−0.187, p=0.405; paired p=0.483). Author basal = 0 cells.

**What does not hold.** This is not a within-tumor progression proof: the n=22 DPT Spearman mixes nLung and tLung. tLung-only CLDN4 vs DPT is n=11, ρ=0.645, p=0.032 — small-n, not a lineage claim. PAGA is one connected component; that does not mean AT2 becomes LUAD. Not ICI.

## Verdict

Sample-level CLDN4 vs AT2-rooted DPT: n=22, ρ=0.618, p=0.00216. CLDN4 vs AT2 score: n=22, ρ=-0.833, p=1.51e-06. CLDN4 vs barrier/keratin (CLDN4 excluded from the score): n=22, ρ=0.775, p=2.26e-05. CLDN4 vs malignant-like: n=22, ρ=0.721, p=0.000153. CLDN4 vs TACSTD2 (comparator only): n=22, ρ=0.863, p=2.3e-07. tLung-only CLDN4 vs DPT: n=11, ρ=0.645, p=0.032 (p<0.05). Paired CLDN4-high vs low barrier/keratin: n=22, W=0.0, Δmed=0.304, p=4.77e-07. Paired CLDN4-high vs low AT2: n=22, W=3.0, Δmed=-1.003, p=2.38e-06. PAGA has 1 component(s) at connectivity>0 among 21 Leiden vertices. Author basal n_cells=0; basal score is a KRT5/KRT15/TP63/NGFR proxy, not an author state. The mixed nLung+tLung DPT correlation is not a within-tumor progression test. Not ICI. Not a TACSTD2 redo. No both-high gate.

## Honest n

- Cells after QC: **n_cells = 10973** (nLung 3703, tLung 7270).
- Samples: **n_samples = 22** (nLung 11, tLung 11).
- Patients: **n_patients = 12** (LUNG_Nxx / LUNG_Txx numeric id; unpaired N01 and T25).
- Samples with ≥10 epithelial cells used for Spearman: **n = 22**.
- Author subtypes (cells): AT2 2020, Club 439, AT1 530, Ciliated 654, tS1 3270, tS2 3018, tS3 64, NA/other 978. **Author basal = 0.**
- CLDN4 tertile cells: low 3658, mid 3657, high 3658.
- Samples with ≥8 cells in both CLDN4-high and CLDN4-low arms (paired extra test): **n = 22**.
- Genes absent from locked sets: **none**.
- GSE207422 not used (NSCLC, not LUAD-only; author cell labels not on GEO). GSE253013 skipped (processed RDS ~9.3 GB).

## Locked choices

- Leiden resolution 0.6; HVG 3000; neighbors 30; PCs 30.
- DPT root: nLung author AT2 (median AT2 score) (root cell index 1217, sample LUNG_N34).
- PAGA components at connectivity>0: **1** among 21 Leiden vertices.
- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).

## Primary (sample-level Spearman, BH inside this list)

| Contrast | n_samples | ρ | p | q |
| --- | ---: | ---: | ---: | ---: |
| CLDN4 vs DPT | 22 | 0.618 | 0.00216 | 0.00278 |
| CLDN4 vs AT2 score | 22 | -0.833 | 1.51e-06 | 4.53e-06 |
| CLDN4 vs club score | 22 | -0.187 | 0.405 | 0.405 |
| CLDN4 vs basal score | 22 | 0.610 | 0.00255 | 0.00287 |
| CLDN4 vs barrier/keratin (no CLDN4) | 22 | 0.775 | 2.26e-05 | 5.08e-05 |
| CLDN4 vs malignant-like score | 22 | 0.721 | 0.000153 | 0.000275 |
| CLDN4 vs TACSTD2 (comparator) | 22 | 0.863 | 2.3e-07 | 2.07e-06 |
| SFTPC vs DPT (control) | 22 | -0.674 | 0.000589 | 0.000883 |
| AT2 score vs DPT (control) | 22 | -0.850 | 5.59e-07 | 2.51e-06 |

## Sensitivity (origin-stratified; not in the BH family)

| Contrast | n_samples | ρ | p |
| --- | ---: | ---: | ---: |
| tLung-only CLDN4 vs DPT | 11 | 0.645 | 0.032 |
| nLung-only CLDN4 vs DPT | 11 | 0.682 | 0.0208 |
| tLung-only CLDN4 vs barrier/keratin (no CLDN4) | 11 | 0.673 | 0.0233 |
| tLung-only CLDN4 vs AT2 | 11 | -0.864 | 0.000612 |
| tLung-only CLDN4 vs TACSTD2 | 11 | 0.791 | 0.00375 |

## Extra figure — CLDN4-high vs CLDN4-low (sample-paired)

Emitted: **True**. Rule: sample-level Spearman(CLDN4, barrier_keratin_no_CLDN4) ρ>0 and p<0.05, or any paired tertile Wilcoxon p<0.05. Observed Spearman n=22, ρ=0.775, p=2.26e-05.

| Paired contrast (high − low) | n_samples | Δ median | p |
| --- | ---: | ---: | ---: |
| barrier/keratin (no CLDN4) high vs low | 22 | 0.304 | 4.77e-07 |
| AT2 high vs low | 22 | -1.003 | 2.38e-06 |
| malignant-like high vs low | 22 | 0.121 | 2.62e-05 |
| DPT high vs low | 22 | 0.148 | 9.54e-07 |
| club high vs low | 22 | 0.039 | 0.483 |

## What this does not claim

- Cell-level p-values are not the claim. n_cells is large by construction.
- Malignant-like is author tS1/tS2/tS3 and/or CEACAM5/6/MKI67 — **not CNV**.
- Author basal n_cells=0; basal score is a KRT5/KRT15/TP63/NGFR proxy.
- GSE131907 is treatment-naive. Do not write ICI / MPR / RECIST language.
- No TACSTD2∩CLDN4 both-high gate.
- Do not write “AT2 differentiates into LUAD because PAGA is connected.”
- Do not write “CLDN4 marks the malignant terminal.”
- The mixed nLung+tLung DPT correlation is not a within-tumor progression test.

## Outputs

- `results/figures/fig_trajectory_cldn4.png` — PAGA + UMAP CLDN4 + DPT + sample CLDN4 vs DPT
- `results/figures/fig_extra_cldn4_tertile.png` — extra CLDN4-high vs low programs
- `results/figures/fig_honest_n.png` — author subtype × origin counts
- `results/tables/sample_means.tsv`, `sample_level_spearman.tsv`, `leiden_paga_vertices.tsv`
- `results/summary.json`

## Reproduce

```bash
pip install -r methods/scrna_paga_cldn4/requirements.txt
bash methods/scrna_paga_cldn4/scripts/download.sh /tmp/scrna_paga_cldn4_data
python3 methods/scrna_paga_cldn4/scripts/extract_epithelium.py \
  --data /tmp/scrna_paga_cldn4_data \
  --out /tmp/scrna_paga_cldn4_data/epithelium.h5ad
python3 methods/scrna_paga_cldn4/scripts/analyze_paga_cldn4.py \
  --input /tmp/scrna_paga_cldn4_data/epithelium.h5ad \
  --outdir methods/scrna_paga_cldn4/results \
  --finding methods/scrna_paga_cldn4/FINDING.md
```

