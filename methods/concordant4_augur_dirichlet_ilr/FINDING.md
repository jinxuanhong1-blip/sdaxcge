# Concordant-4: Dirichlet / ILR composition and Augur-style priority

ADDITIVE. **CLDN4-only.** The locked malignant CLDN4 % positive and the within-cohort quartiles are copied, not recomputed. Cohorts stay GSE123902 + GSE131907 + GSE205335 + GSE189357. Not GSE148071, GSE127465, GSE154826, GSE200563, or E-MTAB-13526.

The locked patient-level result is reproduced and then split. It is not replaced. Unit = donor / sample / patient. Do not quote cell counts as n.

## Honest n

- **n = 65** units (13 donors + 21 samples + 22 patients + 9 patients).
- Author T + NK equals the locked T/NK count on all 43 GSE131907 and GSE205335 units. Malignant counts and cell totals match the locked table on all 65.
- Marker cohorts (GSE123902, GSE189357): the locked T/NK gate also matches on all 22 units. The T and NK counts in the models add CD4/CD3G T cells and NCAM1/NCR1 NK cells that the locked gate left out (12,095 cells across those 22 units).
- Quartile contrast uses the locked labels: **Q1 = 19, Q4 = 16**.
- Augur uses those Q1/Q4 units only (**35** units, 4,095 subsampled cells). That cell count is not the test n.

## 1. Fractions (DerSimonian–Laird on Spearman ρ)

Primary covariate = locked malignant CLDN4 % positive. I² = 0 on every row below.

| fraction | N | ρ (p, 95% CI) | stacked Q4 vs Q1 r (19/16, p) |
|---|---:|---|---|
| locked T/NK gate | 65 | −0.531 (1.65×10⁻⁵, −0.697 to −0.312) | −0.724 (2.88×10⁻⁴) |
| T | 65 | −0.516 (3.2×10⁻⁵, −0.686 to −0.293) | −0.625 (0.0018) |
| NK | 65 | −0.351 (0.0076, −0.562 to −0.097) | −0.520 (0.0093) |
| myeloid | 65 | −0.220 (0.10, −0.456 to +0.046) | −0.395 (0.049) |
| malignant | 65 | +0.556 (5.1×10⁻⁶, +0.343 to +0.714) |  |

The locked gate row matches the published concordant-4 number (ρ = −0.531, p = 1.65×10⁻⁵, I² = 0). T carries that association. NK moves the same way and is weaker. The myeloid fraction is not a clear continuous association. GSE189357 myeloid ρ is +0.32 (n = 9); the other three cohorts are negative. Malignant fraction goes up with CLDN4 % positive.

## 2. ILR balances

4-part simplex (T, NK, myeloid, Rest), half-cell pseudocount, cohort indicator (reference GSE131907) plus within-cohort CLDN4 z-score. Primary p-values are 1999 within-cohort permutations of that z-score. Rest contains the malignant cells.

| balance | β per SD | perm p | Q4 vs Q1 rank-biserial (perm p) |
|---|---:|---:|---|
| T/NK/myeloid vs Rest | −0.375 | 0.002 | −0.572 (0.001) |
| T/NK vs myeloid | −0.181 | 0.045 | −0.349 (0.037) |
| T vs NK | −0.120 | 0.20 | −0.026 (0.83) |

Joint Wilks test on the three balances: perm p = 0.004.

The T-versus-NK balance does not move (Spearman meta ρ = −0.054, p = 0.69). The lymphoid-versus-myeloid balance is the modest reshape: cohort-adjusted ILR perm p = 0.045, while the rank meta-analysis of the same balance is ρ = −0.226, p = 0.093. Those two internal balances are algebraically the same on the renormalized T/NK/myeloid simplex, so they are not a second test. The within-immune Dirichlet below is the non-redundant model, because its total is T+NK+myeloid rather than all cells. Joint Wilks on that 2-balance block: perm p = 0.040.

## 3. Dirichlet-multinomial and CLR

Dirichlet-multinomial, shared precision, 299 within-cohort permutations. Fraction change is the mean predicted proportion at +0.5 SD minus −0.5 SD.

4-part, reference = Rest, φ = 9.63:

| part | log-ratio vs Rest, per SD | perm p | Δ fraction |
|---|---:|---:|---:|
| T | −0.526 | 0.0033 | −0.085 |
| NK | −0.331 | 0.013 | −0.006 |
| myeloid | −0.270 | 0.020 | −0.006 |
| Rest | reference |  | +0.096 |

Within T/NK/myeloid, reference = myeloid, φ = 11.5: T versus myeloid β = −0.286, perm p = 0.023 (T's share of that pool −0.055; myeloid's share +0.053). NK versus myeloid β = −0.165, perm p = 0.12.

Centered log-ratio (each part versus the geometric mean of all four), Benjamini–Hochberg across the four parts:

| part | CLR β per SD | perm p | q |
|---|---:|---:|---:|
| T | −0.267 | 0.001 | 0.004 |
| Rest | +0.325 | 0.005 | 0.010 |
| NK | −0.098 | 0.36 | 0.48 |
| myeloid | +0.040 | 0.53 | 0.53 |

Myeloid's centered log-ratio is flat. Its negative Dirichlet coefficient versus Rest is the expansion of Rest, and the absolute fraction change is about half a percentage point. T is the part that falls relative to the geometric mean.

## 4. Augur-style cell-type priority

Q1 versus Q4 only. Up to 25 cells per type per unit, then 12 per unit inside each repeat. Features: shared 343-gene activity panel, log library-size normalized, cohort-residualized inside the training fold. CLDN4 is held out everywhere. Malignant cells also drop EPCAM, KRT7/8/18/19, CLDN3/7, and TACSTD2. Classifier: L2 logistic (C = 0.2), StratifiedGroupKFold grouped by unit. Primary score = AUC of the unit-mean out-of-fold probability (12 repeats). Null = 99 within-cohort shuffles of the Q label. A random forest was fit first; its null scores sat below 0.5, so it is not the reported scale (`augur_priority_forest_sensitivity.tsv`). In that sensitivity only malignant beat the forest null.

| cell type | units Q1/Q4 | unit AUC | cell AUC | perm p | BH q |
|---|---|---:|---:|---:|---:|
| malignant | 19/16 | 0.813 | 0.732 | 0.030 | 0.15 |
| B | 19/12 | 0.660 | 0.591 | 0.13 | 0.33 |
| myeloid | 19/16 | 0.601 | 0.554 | 0.27 | 0.36 |
| T | 19/16 | 0.579 | 0.547 | 0.29 | 0.36 |
| NK | 18/14 | 0.288 | 0.383 | 0.95 | 0.95 |

Label-shuffle means sit at 0.47–0.51, so 0.5 is the right null center for this logistic. Malignant cells are the top rank after the CLDN4/epithelial genes are removed. Across five cell types the BH q is 0.15. T, NK, and myeloid are not prioritized. NK's point estimate is below 0.5 (two-sided perm p = 0.09) and is not read as a signed program.

## What this adds

CLDN4-high units in this locked set have a higher malignant fraction and a lower T fraction. NK declines in the same direction without a T-versus-NK rebalancing. Myeloid is not selectively lost. A patient-blocked classifier does not rank T, NK, or myeloid transcriptomes as the compartment that tracks the quartile. The malignant transcriptome still separates Q4 from Q1 after CLDN4 and the gate genes are removed, at a nominal permutation p that does not clear a five-test BH threshold of 0.10.

## Reproduce

```bash
bash methods/concordant4_augur_dirichlet_ilr/scripts/download.sh /tmp/geo_c4
python3 methods/concordant4_augur_dirichlet_ilr/scripts/test_compositional.py
python3 methods/concordant4_augur_dirichlet_ilr/scripts/build_counts.py
python3 methods/concordant4_augur_dirichlet_ilr/scripts/fit_compositional.py
Rscript methods/concordant4_augur_dirichlet_ilr/scripts/extract_gse205335.R
python3 methods/concordant4_augur_dirichlet_ilr/scripts/extract_augur.py
python3 methods/concordant4_augur_dirichlet_ilr/scripts/run_augur.py
python3 methods/concordant4_augur_dirichlet_ilr/scripts/make_figures.py
```
