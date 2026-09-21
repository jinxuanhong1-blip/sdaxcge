# Methods

Observational test of the claim **“TROP2–immune depends partly on CLDN4.”**

Exposure is malignant TACSTD2 (TROP2). The candidate mediator is malignant CLDN4. The immune outcome is pre-specified and different in each platform. No coefficient was chosen after a search over genes, radii, or cut points.

## Concordant-4

Patient table `data/concordant4_patient_scores.tsv` (source in `data/SOURCE.txt`). Four cohorts only: GSE123902 (n=13), GSE131907 (n=21), GSE205335 (n=22), GSE189357 (n=9). Unit = patient / donor / sample. N = 65. Cell counts are not n.

Two readouts of the same malignant cells:

- percent detected (raw UMI > 0)
- mean log1p(raw UMI)

Outcome: `frac_tnk`.

Within each cohort, Spearman ρ. Partial Spearman holding the other gene fixed:

ρ_{XY·Z} = (ρ_XY − ρ_XZ ρ_YZ) / √[(1−ρ_XZ²)(1−ρ_YZ²)]

The same partial from Pearson correlation of rank residuals matched the formula to numerical noise. Cohort coefficients are pooled with DerSimonian–Laird on Fisher z. Variance is 1/(n−3) for a raw Spearman and 1/(n−4) for a partial Spearman with one control.

ACME is the product-method estimator on variables standardized inside the cohort (sample SD). a is the slope of CLDN4 on TACSTD2. b is the slope of T/NK on CLDN4 in the regression that also contains TACSTD2. ACME = a×b. The 95% interval is a patient bootstrap (4,000 resamples, seed 4639). Cohort ACMEs are pooled with DerSimonian–Laird using the bootstrap standard errors.

A Spearman of all 65 rows, ignoring cohort, is written as a sensitivity. It is not the primary pool.

Pipeline check: DerSimonian–Laird ρ for CLDN4 percent detected versus T/NK must equal −0.5311678045689989. The script stops otherwise.

## He 2022 CosMx

figshare 25976224, `cosmx_human_nsclc_clustered.h5ad` (2,755,776,882 bytes). Malignant cells are the author `cell_type` matched to the section (`tumor 5/6/9/12/13`): 295,877 cells, 8 sections, 5 donors. Expression is log1p(count / n_counts × 10,000). Neighbors are CD8 or NK cells inside 50 µm and inside 100 µm. Coordinates are section pixels × 0.18 µm/px. Each section has its own tree. Median nearest-neighbor distance was required to fall between 3 and 40 µm.

The coefficient is computed inside each donor (cells from that donor’s sections pooled) and, separately, inside each section. The matrix reports the unweighted mean of the five donor coefficients. The interval resamples the five donors. Cell-level p-values are stored and are not used as the claim p-value, because cells inside a section are not independent patients.

The continuous CLDN4 Spearman versus CD8+NK count at 50 µm, averaged across sections, is checked against the prior sweep (−0.036980967). That check is the continuous association. The locked high-versus-low cytotoxic ratios (0.36 at 50 µm, 0.52 at 100 µm) are a different contrast and are not recomputed.

A between-section Spearman of section-mean expression versus section-mean neighbor count (n=8) is a sensitivity. Leave-one-section results are in `results/tables/cosmx_between_section_loo.tsv`.

## Attenuation and the call

Percent attenuation = 100 × (ρ_TACSTD2 − ρ_TACSTD2|CLDN4) / ρ_TACSTD2. It is left blank when |ρ_TACSTD2| < 0.05. A value above 100% means the partial crossed or passed zero. A negative value means the partial moved away from zero. Those values are shown and are not read as a mediated proportion.

The call uses the 95% interval of the raw TACSTD2–immune Spearman and the ACME interval.

- **support:** the raw interval lies entirely below 0, attenuation is between 20% and 100%, the ACME interval lies entirely below 0, and the CLDN4 partial stays negative.
- **partial:** raw ρ ≤ −0.20, the partial is less negative than the raw coefficient, and the ACME point estimate is negative, while at least one support interval fails.
- **opposite:** the raw interval lies entirely above 0, or the TACSTD2-given-CLDN4 partial is positive with p < 0.05 while the raw association is not significantly negative.
- **null:** every remaining row.

## What this is

The ACME is an observational product. It is not an experiment that changes TACSTD2 or CLDN4. Donor and patient are the independent units.
