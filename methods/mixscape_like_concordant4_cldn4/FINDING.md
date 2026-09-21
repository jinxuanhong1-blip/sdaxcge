# Mixscape-like CLDN4 signature vs Hallmark IFN (concordant-4)

ADDITIVE. **CLDN4-only.** Not a CRISPR screen and not Seurat `RunMixscape`.
There are no gRNAs, no non-targeting guides, and no escapee mixture model.
CLDN4-low malignant cells are an observational **KD-like** proxy;
CLDN4-high malignant cells are the **NT-like** neighbor pool.
The four locked cohorts only: **GSE123902 + GSE131907 + GSE205335 + GSE189357**.
Not GSE148071, GSE127465, GSE207422, GSE154826, or GSE200563.

This does not replace the locked patient-level result
(malignant CLDN4 %pos vs T/NK ρ = −0.531, n = 65) or the malignant
pseudobulk IFN logFC of −0.609. Those are different estimands.
The unit remains the patient / donor / sample. Do not quote cell counts as n.

## Method, in one paragraph

Inside each unit, malignant cells are split by CLDN4 count into value quartiles
(ties stay in the lower bin; zero-inflated genes can leave Q2 or Q3 empty).
Q4 is the NT-like pool. For every malignant cell the local signature is
log1p(CP10k) minus the mean of its k nearest Q4 neighbors (k = 20, or n_Q4 − 1)
in a within-unit PCA (15 PCs). PCA genes are highly variable genes after removing
CLDN4, Hallmark IFNα, Hallmark IFNγ, the custom MHC-I/APM list, Hallmark
spermatogenesis, and the detection-matched control genes, so the neighbor graph
is not built on the readout. The IFN score is the mean residual of
Hallmark IFNα ∪ IFNγ. Similarity to that set is the KD-like residual of IFN
minus a detection-matched gene set of the same size. Spermatogenesis, matched
the same way, is the negative-control hallmark. p-values are descriptive.

## Who entered

Locked units = 65. Status counts: ok=60, too_few_nt=5.
Units with a within-unit Spearman (non-Q4 cells, ≥2 CLDN4 values) = **59**.
Out of the signature test: GSE123902:LX699 (too_few_nt); GSE123902:LX701 (too_few_nt); GSE131907:EBUS_13 (too_few_nt); GSE131907:NS_16 (too_few_nt); GSE205335:P4001 (too_few_nt).
Those units do not have 15 or more CLDN4-high (Q4) malignant cells, so k nearest
NT-like neighbors are not defined. One additional scored unit has a Spearman of NA
because every non-Q4 cell shares the same CLDN4 count (ties at zero).

## 1. Dose across CLDN4 quartiles

Equal-unit mean of the per-cell Hallmark IFN local residual.
Q4 is a self-neighborhood null (each Q4 cell minus other Q4 neighbors), so it is
expected to sit near zero even when Q1–Q3 are real contrasts against that same pool.
The dose among Q1–Q3 uses one shared NT-like pool and is the fairer slope.

| quartile | role | n units | mean IFN residual | mean IFN − matched | mean sperm residual |
|---|---|---:|---:|---:|---:|
| Q1 | lowest CLDN4, KD-like | 60 | -0.069 | 0.004 | -0.023 |
| Q2 | lower-mid | 41 | -0.049 | -0.004 | -0.013 |
| Q3 | upper-mid | 59 | -0.026 | -3.259e-04 | -0.007 |
| Q4 | highest CLDN4, NT-like self-null | 60 | -0.002 | 4.931e-04 | 6.696e-05 |

Paired Wilcoxon on unit means, IFN residual Q1 − Q3: median diff -0.047, n = 59, p = 4.443e-08.
IFN residual Q1 − Q4 (KD-like vs self-null): median diff -0.060, n = 60, p = 7.541e-10.

## 2. Within-unit Spearman (primary dose test)

One Spearman per unit, CLDN4 count vs IFN local residual, **non-Q4 cells only**
(Q4 residuals are the self-null and are not in this correlation).
Units are equally weighted. The cell-level Spearman p inside a unit is not used.

CLDN4-higher malignant cells carry a higher Hallmark IFN local residual than CLDN4-lower cells in the same unit.

Equal-unit mean ρ = **0.331** (median 0.371; 5 negative / 54 positive; n = 59). One-sample t on Fisher z: p = 1.389e-12. Wilcoxon signed-rank: p = 2.366e-09.

| cohort | unit | n | mean ρ | median ρ | n negative |
|---|---|---:|---:|---:|---:|
| GSE123902 | donor | 10 | 0.167 | 0.192 | 1 |
| GSE131907 | sample | 19 | 0.508 | 0.467 | 0 |
| GSE205335 | patient | 21 | 0.326 | 0.358 | 2 |
| GSE189357 | patient | 9 | 0.093 | 0.033 | 2 |

DerSimonian–Laird across the 4 cohort means (weight = 1 / SE² of the within-cohort Fisher z): ρ = 0.284 (p = 0.01682, I² = 95.0%, 95% CI 0.053 to 0.486, N = 59).
I² is high: the sign is shared more than the size. GSE131907 is the steep cohort;
GSE189357 is shallow. The equal-unit mean is not a single common effect size.

This within-unit local residual is a different estimand from the locked
between-patient pseudobulk (IFN logFC −0.609 in CLDN4-high patients).
Here, CLDN4-low cells are lower, not higher, on the Hallmark IFN residual
than the CLDN4-high neighbors they were matched to.

## 3. Similarity of the KD-like signature to Hallmark IFN

KD-like cells are the lowest occupied quartile with at least 20 cells (Q1 when it qualifies).
The competitive score is mean IFN residual minus mean residual of genes matched one-to-one
on detection rate. Spermatogenesis is scored the same way after IFN genes are removed from it.

Wilcoxon on the KD-like IFN − matched delta: median -0.001, n = 60, p = 0.9179.
Wilcoxon on the KD-like spermatogenesis − matched delta: median 0.004, n = 60, p = 0.001329.
Cosine of the KD-like gene residual with the Hallmark IFN indicator: median -0.080 (n = 60). The cosine is at or below zero.
A gene-label shuffle is not a centered null here, because KD-like residuals are
broadly negative, so that permutation p-value is not used as an IFN-enrichment call.

The similarity call is the detection-matched delta. It sits on zero.
Hallmark IFN does not rise above genes with the same detection rate in the KD-like
local signature. The quartile slope of the raw IFN residual is therefore not an
IFN-specific transfer.

## 4. Boundaries

- Not true CRISPR. Do not call these cells knockouts, escapees, or gRNA-perturbed.
- Not the locked T/NK exclusion result and not a spatial analysis.
- Quartiles are within-unit CLDN4 values. Many malignant cells have CLDN4 = 0, so
  Q2/Q3 are often empty. Empty bins are omitted, not filled by random tie breaks.
- Q4’s residual is mechanically shrunk by matching Q4 cells to other Q4 cells.
- Neighbor matching removes shared state only along the PCA axes that were kept.
  It is not a causal estimate of CLDN4 deletion.

## Reproduce

```
bash methods/mixscape_like_concordant4_cldn4/scripts/download.sh /tmp/geo_mixscape
python3 methods/mixscape_like_concordant4_cldn4/scripts/mixscape_like.py
python3 methods/mixscape_like_concordant4_cldn4/scripts/summarize.py
```

