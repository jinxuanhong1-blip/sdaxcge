# GSE131907 PAGA / diffusion: TACSTD2 and CLDN4 on LUAD epithelium

**Slice:** `methods/scrna_paga/` + `scripts/scrna_paga/` + `results/scrna_paga/` only.
**Additive.** Does not rewrite `methods/scrna/` or `methods/trajectory/`.
**Cohort:** GSE131907 (Kim et al. 2020), LUAD, public UMI. tLung epithelium + nLung epithelium (AT2/club prior).
**ICI labels:** none. This slice cannot test ICI / MPR / RECIST.

## What was tested

PAGA connectivities and diffusion pseudotime on author-labeled LUAD epithelium, then sample-level Spearman of TACSTD2 and CLDN4 (separately) versus DPT and versus pre-specified AT2 / club / basal / barrier-keratin / malignant-like scores.

## Coverage (honest missingness)

- Cells in object: **n_cells = 10973** (nLung 3703, tLung 7270).
- Samples: **n_samples = 22** (nLung 11, tLung 11).
- Patients: **n_patients = 12** (LUNG_Nxx / LUNG_Txx numeric id; unpaired N01 and T25).
- Author subtypes (cells): AT2 2020, Club 439, AT1 530, Ciliated 654, tS1 3270, tS2 3018, tS3 64, NA 918. **Author basal = 0.**
- Genes absent from locked sets: {'AT2': [], 'AT1': [], 'club': [], 'basal': [], 'ciliated': [], 'barrier_keratin': [], 'malignant_like': []}.
- Samples with ≥10 epithelial cells used for sample-level Spearman: **n = 22**.
- GSE253013 skipped (processed RDS ~9.3 GB, over public-file cap). GSE207422 not pooled (NSCLC, not LUAD-only).

## Locked choices

- Leiden resolution 0.6; HVG 3000; neighbors 30; PCs 30.
- DPT root: nLung author AT2 (median AT2 score) (root cell index 1217).
- PAGA components at connectivity>0: 1.

## Primary (sample-level Spearman)

| Contrast | n_samples | ρ | p | q |
| --- | ---: | ---: | ---: | ---: |
| TACSTD2 vs DPT | 22 | 0.617 | 0.00221 | 0.00289 |
| CLDN4 vs DPT | 22 | 0.618 | 0.00216 | 0.00289 |
| TACSTD2 vs CLDN4 | 22 | 0.863 | 2.3e-07 | 2.94e-06 |
| TACSTD2 vs AT2 score | 22 | -0.843 | 8.43e-07 | 3.94e-06 |
| TACSTD2 vs club score | 22 | 0.045 | 0.844 | 0.844 |
| TACSTD2 vs basal score | 22 | 0.697 | 0.00031 | 0.000619 |
| TACSTD2 vs barrier/keratin score | 22 | 0.854 | 4.2e-07 | 2.94e-06 |
| TACSTD2 vs malignant-like score | 22 | 0.616 | 0.00227 | 0.00289 |
| CLDN4 vs AT2 score | 22 | -0.833 | 1.51e-06 | 5.29e-06 |
| CLDN4 vs club score | 22 | -0.187 | 0.405 | 0.436 |
| CLDN4 vs basal score | 22 | 0.610 | 0.00255 | 0.00298 |
| CLDN4 vs barrier/keratin score | 22 | 0.817 | 3.48e-06 | 9.74e-06 |
| CLDN4 vs malignant-like score | 22 | 0.721 | 0.000153 | 0.000357 |
| SFTPC vs DPT (control) | 22 | -0.674 | 0.000589 | 0.00103 |

## Sensitivity (origin-stratified; not in the BH family)

| Contrast | n_samples | ρ | p |
| --- | ---: | ---: | ---: |
| tLung-only TACSTD2 vs DPT | 11 | 0.345 | 0.298 |
| tLung-only CLDN4 vs DPT | 11 | 0.645 | 0.032 |
| nLung-only TACSTD2 vs DPT | 11 | 0.836 | 0.00133 |
| nLung-only CLDN4 vs DPT | 11 | 0.682 | 0.0208 |
| tLung-only TACSTD2 vs barrier/keratin | 11 | 0.627 | 0.0388 |

## Extra barrier/keratin figure

Emitted: **True**. Rule: sample-level Spearman(TACSTD2, barrier_keratin) ρ>0 and p<0.05. Observed n=22, ρ=0.854, p=4.2e-07.

## Verdict

Sample-level TACSTD2 vs DPT: n=22, ρ=0.617, p=0.00221. CLDN4 vs DPT: n=22, ρ=0.618, p=0.00216. TACSTD2 vs CLDN4 co-expression: n=22, ρ=0.863, p=2.3e-07. TACSTD2 vs barrier/keratin: n=22, ρ=0.854, p=4.2e-07; extra figure yes. PAGA has 1 component among 21 Leiden vertices at connectivity>0. tLung-only TACSTD2 vs DPT: n=11, ρ=0.345, p=0.298 (inconclusive). Author basal n_cells=0; basal score is a KRT5/KRT15/TP63/NGFR proxy, not an author state. Club vs TACSTD2 is a sample-level null. Do not write that TACSTD2 marks club on this object. The n=22 DPT correlations mix nLung and tLung; they are not a within-tumor progression test.

## Caveats

- Cell-level p-values are not the claim. n_cells is large by construction.
- Malignant-like is author tS1/tS2/tS3 and/or a weak CEACAM5/6/MKI67 score — **not CNV**.
- Club/basal are airway programs. If PAGA disconnects them from AT2, that is a discrete-state result, not a failed download.
- GSE131907 is treatment-naive. Do not write ICI language.
- No TACSTD2∩CLDN4 both-high gate.
- Do not write “AT2 differentiates into LUAD because PAGA is connected.”
- Do not write “TROP2 marks the malignant terminal.”

## Reproduce

```bash
pip install -r requirements.txt
bash scripts/scrna_paga/download.sh /tmp/scrna_paga_data
python3 scripts/scrna_paga/extract_epithelium.py --data /tmp/scrna_paga_data --out /tmp/scrna_paga_data/epithelium.h5ad
python3 scripts/scrna_paga/analyze_paga.py --input /tmp/scrna_paga_data/epithelium.h5ad --outdir results/scrna_paga
```

