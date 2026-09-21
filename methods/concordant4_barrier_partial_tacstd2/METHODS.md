# Methods — barrier outgoing, TACSTD2 partial dependence on CLDN4

ADDITIVE. Concordant four only: GSE123902, GSE131907, GSE205335, GSE189357.
Does not replace the locked CLDN4 CellChat table or the locked +25.6 percentage-point
result. GSE148071, GSE127465, GSE154826, GSE200563, and E-MTAB-13526 are not added.
This is an expression ligand–receptor contrast. It is not a spatial exclusion test.
The patient is the unit. Cell-pooled p-values are not reported.

Senders are malignant cells. The receiver is T/NK.
Barrier ligands are F11R, NECTIN2, CDH1, and LGALS9, the same edges as the locked
concordant-4 communication runs:

- F11R–LFA-1 (ITGAL+ITGB2), F11R–F11R
- NECTIN2–TIGIT, NECTIN2–CD96
- CDH1–integrin (ITGAE+ITGB7), CDH1–KLRG1
- LGALS9–TIM3 (HAVCR2), LGALS9–CD44, LGALS9–CD45 (PTPRC)

PVRL2 and JAM1 are renamed to NECTIN2 and F11R when the official symbol is absent.
An edge is counted only when the receptor complex is in ≥10% of T/NK cells and the
ligand is in ≥10% of at least one sender arm. A complex receptor uses the geometric
mean of subunit means and the minimum subunit proportion. A ligand score is the mean
of its detected edges. An undetected ligand contributes 0. The family score is the
mean of the four ligands.

## Scores (fixed; not ranked)

| score | definition |
|---|---|
| expr_prop_pp | 100 × (fraction ligand-positive in the high arm − fraction in the low arm) |
| cellchat_hill | Hill probability on 10% trimmed log1p(CP10k), Kh = 0.5, no population-size weight. Δ = Hill(Lhigh)·Hill(R) − Hill(Llow)·Hill(R) |
| cpdb_means | CellPhoneDB `lr_means` on arithmetic log1p means. The receptor term cancels, so Δ = (Lhigh − Llow) / 2 |
| liana_logfc | (mean log1p high − mean log1p low) / ln(2) |

CP10k uses the full-transcriptome library size. Expression is not capped. There is
no cell downsample. Co-primary readouts for the shrink call are `expr_prop_pp` and
`cellchat_hill`. The other two scores are reported with the same rule.

## Cells

Labels match the locked CellChat inventory (65 units).

- GSE123902 and GSE189357: epithelium marker-malignant (EPCAM, KRT8, KRT18, or KRT19 count > 0, and PTPRC = 0). T/NK = CD3D, CD3E, CD8A, NKG7, GNLY, or KLRD1, and not malignant. GSE123902 keeps the primary tumor when a donor also has a metastasis. Normal is out.
- GSE131907: author subtype in {Malignant cells, tS1, tS2, tS3}. T/NK = author cell type in {T lymphocytes, NK cells}.
- GSE205335: author `lineage.sub` malignant. T/NK = author `lineage.total`. Normal tissue is out.

A unit enters a gate when malignant cells ≥ 40, T/NK ≥ 20, and the gate gene separates:
ordinal Q4 vs Q1 (low = rank ≤ floor(n/4), high = rank > ceil(0.75 n)), each arm ≥ 10,
and the high-arm mean of the gate gene is greater than the low-arm mean. Decile
(outer 10%, same arm floor, so n_mal ≥ 100) is a sensitivity split. It is not used
for the shrink call.

## Contrasts

The gate is TACSTD2 or CLDN4. The covariate is the other gene.

1. **crude** — Q4 vs Q1 of the gate. Same cells and same rank rule as the locked CLDN4 run.
2. **strata_resid** (primary residual) — inside the patient, malignant cells are split into up to four quantile strata of the covariate. Tied covariate values stay in one stratum. Each barrier ligand's positivity and log1p(CP10k) is demeaned inside those strata. Log1p is then shifted back by the malignant grand mean so the Hill input stays on the original scale. The Q4 vs Q1 contrast is recomputed on those residualized ligand values. Detection still uses the raw arms. Same cells as crude.
3. **linear_resid** — the same Q4 vs Q1 contrast after an ordinary-least-squares residual of the ligand on the covariate log1p(CP10k). If the covariate has no variance, the residual contrast equals the crude contrast.
4. **matched** — 1:1 match, without replacement, of gate-Q4 cells to gate-Q1 cells that sit in the same covariate stratum. The larger arm inside a stratum is subsampled with seed 3979. A patient counts when ≥10 pairs are formed. The match rate is pairs / n_Q4. This is a different cell set from crude.
5. **strata_local** — inside each covariate stratum, median-split the gate gene (ordinal rank, ≥5 cells per side). The communication score is computed in the stratum and averaged across strata. A patient counts when ≥2 strata qualify. This is a local gate contrast, not the Q4 vs Q1 contrast.

The crude CLDN4 expression-proportion family mean is checked against the locked
Q4 vs Q1 result (+25.6 percentage points, n=64). A mismatch larger than 1 percentage
point stops the run.

## Shrink rule (fixed before the cohorts are scored)

On the Q4 vs Q1 split, same patients, for each score:

- shrink fraction = 1 − (mean adjusted / mean crude), defined when the crude mean is positive
- paired two-sided Wilcoxon on (crude − adjusted)

**Shrinks** means all three of: crude mean > 0, shrink fraction ≥ 0.20, and paired p < 0.05.
**Does not shrink** means the adjusted mean is at least the crude mean.
Otherwise the call is **does not materially shrink**.
The adjusted contrast is also tested against zero. A remainder can stay positive after a shrink.

The headline partial-dependence test is TACSTD2 crude vs TACSTD2 strata-residualized
on CLDN4, on both co-primary scores. Linear residual and within-stratum matching are
the pre-specified checks. The reverse (CLDN4 residualized or matched on TACSTD2) is
reported so a shared epithelial program is not read as a TACSTD2-only effect.

Co-expression is the within-patient Spearman of malignant TACSTD2 and CLDN4 log1p,
plus the between-stratum R² of each gene on the other gene's strata. The covariate
gap is the difference in covariate log1p between the gate arms. Matching is balanced
only if that gap is near zero in the matched cells while the gate gap remains.

Sign-flip null: 10,000 flips, seed 3979, one-sided for a positive mean, on the
Q4 vs Q1 crude, strata-residual, and linear-residual contrasts for the two
co-primary scores. Cohort heterogeneity is DerSimonian–Laird I² on the four cohort
means, weighted by the square of the patient-level standard error.

## Not done

Dual-high gating. Extra cohorts. Spatial distance. A population-scaled CellChat
probability. A claim that a percentage-point gap is a CellChat probability.
Cell-pooled tests. Choosing the score or the contrast after seeing which one shrinks.
