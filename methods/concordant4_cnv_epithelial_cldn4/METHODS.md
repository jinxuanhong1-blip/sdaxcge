# CNV malignant call on concordant-4

Additive to the locked annotation test (malignant CLDN4 % positive vs T/NK
fraction, DerSimonian–Laird ρ = −0.531, 65 units). The four cohorts stay
GSE123902, GSE131907, GSE205335, and GSE189357. No other accession is merged.

## Units

Same patient / donor / sample list as `methods/seurat_concordant4_cldn4`:

- GSE123902: 13 tumor donors (primary or metastasis; normal lung held out)
- GSE131907: 21 tumor-bearing samples with author malignant n ≥ 20
- GSE205335: 22 patients, normal-tissue cells removed from the denominator
- GSE189357: 9 tumor patients

T/NK fraction is `n_T/NK / n_cells` on that full unit. It is not recomputed
inside the CNV draw.

## Epithelium and the annotation malignant set

- GSE123902 and GSE189357: `(EPCAM|KRT8|KRT18|KRT19) > 0` and `PTPRC == 0`.
  That marker gate is the locked malignant definition, so CNV here is a
  filter on those cells. Matched normal-lung donors in GSE123902 are the
  diploid epithelial controls and are not test units.
- GSE131907: author `Epithelial cells`. Annotation malignant is
  `Cell_subtype == Malignant cells`. T/NK is T lymphocytes plus NK cells.
  Normal-lung (`nLung`) epithelium is the control.
- GSE205335: author `Epithelial cells` in non-normal tissue. Annotation
  malignant is `lineage.sub == Malignant cells`. T/NK is `lineage.total ==
  T/NK cells`. Normal Lung / LN / Brain epithelium is the control.

## Profiles

Autosomes only, hg38 refGene coordinates, longest transcript per symbol.
Chromosome 6 is dropped (MHC / IFN expression artifact). Remaining genes
detected in at least 5% of the unit's drawn cells are averaged in 5 Mb bins.

Drawn cells, fixed seed per unit: up to 800 epithelial and 250 T/NK
(GSE131907: 500 and 150, because the genes-by-cells text is materialized
while it is streamed). Annotation CLDN4 % uses every malignant cell, not
the draw. CNV CLDN4 % uses the CNV+ cells inside the draw.

**InferCNV preliminary (HMM off).** log2(CP10k+1) with the cell's full
library size, subtract the unit's T/NK mean, moving mean of 101 genes
within each chromosome, subtract the cell mean.

**CopyKAT transform (Gao et al. 2021).** `log(sqrt(x)+sqrt(x+1))`, subtract
the cell mean, local-level Kalman smooth (`dV=0.16`, `dW=0.001`, the
`dlmModPoly` setting in the R package), subtract the T/NK median, subtract
the cell mean.

## Normal-epithelium null

Immune-referenced profiles still carry the epithelial expression program.
Each cohort's null is the median 5 Mb profile of its normal epithelial
cells. A cell's CNV residual is its profile minus that null. The positive
threshold is the 95th percentile of residual mean squared error among those
normal epithelial cells, computed once and then applied to test cells.
GSE189357 has no normal epithelium in the GEO matrix; it uses the
GSE123902 normal-lung null and threshold.

## Labels

- **InferCNV+**: residual MSE above that normal-epithelium percentile.
- **CopyKAT aneuploid**: Ward clustering (`k` = 4, 3, or 2, minimum cluster
  size > 10) on the CopyKAT residual. The cluster that holds the largest
  share of the below-median residual cells is diploid. Remaining clusters
  are aneuploid when their median profile is closer, by 1-Wasserstein
  distance, to the aneuploid consensus than to the diploid consensus.
  If the two consensus profiles correlate ≥ 0.6, every cell is diploid
  (CopyKAT v1.2). Correlation in [0.4, 0.6) is labeled low confidence and
  still used. A cluster call is kept only if the cell also exceeds the
  normal-epithelium CopyKAT threshold.
- **Consensus**: InferCNV+ and CopyKAT aneuploid.

CLDN4 % positive = `100 * mean(UMI > 0)` in the called cells. A unit enters
a Spearman test only when that set has at least 20 cells.

## Association

Within each cohort, Spearman of CLDN4 % positive vs the full-unit T/NK
fraction. Cohorts are combined by DerSimonian–Laird on Fisher z, the same
estimator as the locked annotation meta-analysis. Within-cohort CLDN4
quartiles are stacked and compared by Mann–Whitney (rank-biserial r).
The annotation test is repeated on the units that have a CNV call so the
comparison is paired.

## What was not run

The R `copykat` MCMC breakpoint sampler and inferCNV's HMM / Numbat allele
model were not executed. Genomic bins are fixed 5 Mb windows rather than
CopyKAT's posterior segments. Calls are expression inferences, not DNA copy
number.
