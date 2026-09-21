# scVI/scANVI concordant-4: continuous CLDN4 vs T/NK

ADDITIVE. **CLDN4-only.** The locked patient-level result is not replaced: malignant CLDN4 % positive vs T/NK on GSE123902 + GSE131907 + GSE205335 + GSE189357 is still the n=65 association (published ρ=−0.531). This folder asks whether that association holds for **continuous scVI latent CLDN4** after a **patient-batch scVI / scANVI** integration, with a **mixed model whose random effect is the locked unit**.

Not GSE148071, GSE127465, GSE154826, GSE207422, or GSE200563. No dual-high. No TACSTD2 gate. Do not quote the integration cell count as n.

## Honest n

- **n_units = 65** (13 donors + 21 samples + 22 patients + 9 patients).
- Integration subsample (QC, then ≤220 malignant + ≤140 T/NK + ≤40 other per unit): **24218 cells**, 2001 genes, 65 patient batches.
- T/NK fraction and % positive are computed on the **full unit**, before the cap.

## Integration

- scVI 1.3.3 negative binomial, n_latent=20, 2 layers.
- **Batch key = patient** (`dataset|unit_id`). HVGs (2,000, seurat_v3) are selected within dataset; CLDN4 is forced into the model.
- scANVI is initialized from that scVI model. GSE131907 and GSE205335 keep author labels (malignant / T/NK / other). GSE123902 and GSE189357 are **Unknown** during training and are annotated by `predict`.
- Continuous CLDN4 = mean log1p of scVI normalized expression at library size 10,000, inside scANVI-malignant cells of that unit (≥10 cells, otherwise the marker/author seed call).
- Epochs: scVI 174, scANVI 13.

## 1. Replication of ρ=−0.53 (full unit, locked malignant definition)

% positive is 100 × mean(CLDN4 UMI > 0) in author-malignant (GSE131907, GSE205335) or marker-malignant (EPCAM|KRT8|KRT18|KRT19 > 0 and PTPRC == 0) cells. T/NK uses the same gates as the Seurat concordant-4 table. Pooling is DerSimonian–Laird on Fisher-z.

| score | meta |
|---|---|
| malignant CLDN4 %pos | ρ=-0.531 (p=1.65e-05, I²=0.0%, -0.697 to -0.312, N=65) |
| malignant mean log1p(UMI) | ρ=-0.394 (p=0.00245, I²=0.0%, -0.595 to -0.146, N=65) |
| stacked within-cohort %pos Q4 vs Q1 rank-biserial | r=-0.724 (n_Q1=19, n_Q4=16, p=0.000288) |

The published %pos pool was ρ=−0.531, p=1.65×10⁻⁵, I²=0%, N=65. The number above is recomputed from the public matrices, not copied from that table.

### Cohort Spearmans

| cohort | n | %pos ρ | %pos p | seed-latent ρ | seed p | scANVI-latent ρ | scANVI p |
|---|---:|---:|---:|---:|---:|---:|---:|
| GSE123902 | 13 | -0.659 | 0.0142 | -0.544 | 0.0546 | 0.137 | 0.655 |
| GSE131907 | 21 | -0.522 | 0.0152 | -0.653 | 0.00132 | -0.653 | 0.00132 |
| GSE205335 | 22 | -0.435 | 0.0429 | -0.409 | 0.0585 | -0.399 | 0.0657 |
| GSE189357 | 9 | -0.600 | 0.0876 | -0.300 | 0.433 | -0.133 | 0.732 |

## 2. Continuous scVI latent CLDN4

Seed-malignant cells are the locked definition: author malignant where the authors labeled cells, and the marker gate on GSE123902 and GSE189357. scANVI-malignant cells are `predict()` after scANVI saw author labels and treated the two marker cohorts as Unknown. Both averages use the scVI normalized CLDN4 on the integration subsample.

| score | meta |
|---|---|
| seed-malignant scVI CLDN4 | ρ=-0.516 (p=3.2e-05, I²=0.0%, -0.686 to -0.293, N=65) |
| seed Q4 vs Q1 rank-biserial | r=-0.671 (n_Q1=19, n_Q4=16, p=0.000777) |
| scANVI-malignant scVI CLDN4 | ρ=-0.342 (p=0.0762, I²=49.5%, -0.635 to 0.037, N=65) |
| scANVI Q4 vs Q1 rank-biserial | r=-0.375 (n_Q1=19, n_Q4=16, p=0.0614) |

The continuous latent score on the locked malignant cells is the replication of ρ=−0.53 (same sign, I²=0, n=65). Restricting to scANVI-malignant cells moves the pooled Spearman because scANVI, trained on author malignant labels, calls a fraction of marker-gate epithelial cells "other" in the two unlabeled cohorts. That is reported, not folded back into the locked number. Latent column `latent_CLDN4` uses the scANVI call when a unit has at least 10 such cells ({'scanvi': 65}).

## 3. Mixed models (patient random effect)

Each locked unit is one binomial observation: `cbind(n_T/NK, n_other) ~ scale(CLDN4) + dataset + (1 | unit_id)`, logit link. The random intercept is the patient / donor / sample. It is an observation-level random effect (logit-normal extra-binomial variance). The likelihood ratio compares that model to the same model without CLDN4. A plain binomial GLM would treat every cell as independent; the random intercept is what keeps the test at the unit level. Coefficients are per 1 SD of the named CLDN4 score.

- **%pos, patient RE, Wald:** β=-0.6188 per SD, SE=0.1431, stat=-4.32, p=1.53e-05, n_units=65. patient RE sd=1.0357 singular=FALSE family=binomial Wald
- **%pos, patient RE, LRT:** β=-0.6188 per SD, SE=0.1431, stat=16.4, p=5.02e-05, n_units=65. patient RE sd=1.0357 singular=FALSE family=binomial Wald LRT vs dataset + (1|unit)
- **Seed-latent CLDN4, patient RE, Wald:** β=-0.6353 per SD, SE=0.1527, stat=-4.16, p=3.17e-05, n_units=65. patient RE sd=1.0446 singular=FALSE family=binomial Wald
- **Seed-latent CLDN4, patient RE, LRT:** β=-0.6353 per SD, SE=0.1527, stat=15.4, p=8.86e-05, n_units=65. patient RE sd=1.0446 singular=FALSE family=binomial Wald LRT vs dataset + (1|unit)
- **scANVI-latent CLDN4, patient RE, Wald:** β=-0.4574 per SD, SE=0.1397, stat=-3.28, p=0.00106, n_units=65. patient RE sd=1.0893 singular=FALSE family=binomial Wald
- **scANVI-latent CLDN4, patient RE, LRT:** β=-0.4574 per SD, SE=0.1397, stat=9.93, p=0.00162, n_units=65. patient RE sd=1.0893 singular=FALSE family=binomial Wald LRT vs dataset + (1|unit)
- **Gaussian LMM** `frac_tnk ~ seed-latent + (1 | dataset)`: β=-0.1033 per SD, SE=0.0247, stat=-4.19, p=8.92e-05, n_units=65. dataset RE sd=0.0294 singular=FALSE gaussian t df=63
- **Equal-weight OLS** `frac_tnk ~ seed-latent + dataset`: β=-0.1013 per SD, SE=0.0281, stat=-3.61, p=0.000632, n_units=65. OLS equal-unit weight, dataset fixed effect

Cell n is not the sample size. `n_units` above is the number of random-effect levels.

## 4. scANVI annotation

Author/marker seed vs scANVI prediction (cells in the subsample):

| dataset | seed | scANVI | n |
|---|---|---|---:|
| GSE123902 | T/NK | T/NK | 1330 |
| GSE123902 | T/NK | malignant | 42 |
| GSE123902 | T/NK | other | 244 |
| GSE123902 | malignant | T/NK | 155 |
| GSE123902 | malignant | malignant | 1595 |
| GSE123902 | malignant | other | 452 |
| GSE123902 | other | T/NK | 117 |
| GSE123902 | other | malignant | 29 |
| GSE123902 | other | other | 374 |
| GSE131907 | T/NK | T/NK | 2808 |
| GSE131907 | T/NK | malignant | 4 |
| GSE131907 | T/NK | other | 1 |
| GSE131907 | malignant | malignant | 4176 |
| GSE131907 | malignant | other | 1 |
| GSE131907 | other | T/NK | 4 |
| GSE131907 | other | malignant | 10 |
| GSE131907 | other | other | 826 |
| GSE189357 | T/NK | T/NK | 916 |
| GSE189357 | T/NK | malignant | 64 |
| GSE189357 | T/NK | other | 280 |
| GSE189357 | malignant | T/NK | 205 |
| GSE189357 | malignant | malignant | 1258 |
| GSE189357 | malignant | other | 517 |
| GSE189357 | other | T/NK | 52 |
| GSE189357 | other | malignant | 13 |
| GSE189357 | other | other | 295 |
| GSE205335 | T/NK | T/NK | 3040 |
| GSE205335 | T/NK | malignant | 27 |
| GSE205335 | T/NK | other | 1 |
| GSE205335 | malignant | T/NK | 21 |
| GSE205335 | malignant | malignant | 4477 |
| GSE205335 | malignant | other | 4 |
| GSE205335 | other | T/NK | 25 |
| GSE205335 | other | malignant | 49 |
| GSE205335 | other | other | 806 |

Leiden 0.6 on the scANVI latent called **13 / 23** clusters malignant. A cluster is malignant when ≥50% of its author-labeled cells are author-malignant, else when ≥50% of its cells are scANVI-malignant. Table: `results/tables/cluster_annotation.tsv`.

## What this does not say

- The latent neighborhood of a CLDN4-high cell is not a tissue neighborhood. This is not spatial exclusion.
- GSE131907's locked unit remains the tumor-bearing sample. GSE205335 pools a patient's libraries before the test.
- I² and the DL interval describe the four cohort Spearmans. The GLMM asks a different question (log-odds of T/NK per SD of CLDN4, with a unit-level random intercept) and is not a second copy of ρ.

## Reproduce

```bash
python3 methods/scanvi_concordant4_cldn4/download.py --out /tmp/geo_c4
python3 methods/scanvi_concordant4_cldn4/prepare.py --raw /tmp/geo_c4
python3 methods/scanvi_concordant4_cldn4/integrate_model.py
```

Raw matrices are not committed. `results/tables/locked_replication_check.tsv` is the cell-count / %pos diff against the Seurat concordant-4 patient table.
