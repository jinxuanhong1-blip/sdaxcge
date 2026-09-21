# Methods

Patient-level CD8 and NK state versus malignant CLDN4 in the locked concordant-4.
Not a spatial test. Not a re-estimate of the T/NK fraction.

## Question

CosMx NSCLC (He 2022) found fewer cytotoxic neighbors around CLDN4-high tumor cells, and no drop in GZMB, PRF1, NKG7, or IFNG in the neighbors that remain (hi/lo about 1.11–1.22, 0/8 sections down). scRNA has no neighbor graph. The composition analog is: among the same concordant-4 units, do the CD8 cells and the NK cells that are present in high-malignant-CLDN4 patients carry a lower effector score, or a higher exhaustion score, than those in low-malignant-CLDN4 patients?

A lower T/NK fraction is the locked exclusion result. It is not an intensity result.

## Units and exposure

Locked table `data/locked_patient_units.tsv` (n=65):

| cohort | unit | n |
|---|---|---:|
| GSE123902 | donor, tumor or metastasis file | 13 |
| GSE131907 | sample, tumor-bearing | 21 |
| GSE205335 | patient, normal tissue removed | 22 |
| GSE189357 | patient, tumor | 9 |

Exposure is the locked malignant CLDN4 percent positive and the locked within-cohort quartile. Recomputed percent positive is a gate check only. Tests do not switch to the recomputed value.

## Compartments

Author cohorts

- GSE131907 CD8: `Cell_subtype` in Cytotoxic CD8+ T, Exhausted CD8+ T, Naive CD8+ T, CD8 low T. Mixed CD4/CD8 is not included.
- GSE131907 NK: `Cell_type` NK cells, excluding cells whose subtype is one of the four CD8 labels (those stay CD8).
- GSE205335 CD8: `lineage.sub` CD8+ T cells, dropping celltype AMB cells.
- GSE205335 NK: `lineage.sub` NK cells, dropping AMB cells.
- GSE205335 malignant: `lineage.sub` Malignant cells. Normal-tissue samples are out before this cut.

Marker cohorts (GSE123902, GSE189357), same malignant gate as the locked analysis:

- Malignant: (EPCAM or KRT8 or KRT18 or KRT19) > 0 and PTPRC = 0.
- CD8: not malignant, (CD8A or CD8B) > 0, (CD3D or CD3E) > 0, and not CD4-dominant (CD4 greater than both CD8A and CD8B).
- NK: not malignant, not CD8, CD3D = 0, CD3E = 0, CD8A = 0, CD8B = 0, and (KLRD1 or NCR1 or KLRF1 or GNLY or (NCAM1 and PTPRC)) > 0.

GZMB, PRF1, NKG7, and IFNG are not gate genes. GNLY is an NK identity gene and is not part of the four-gene effector score. Gating NK on GNLY does not put GNLY into that score. Two GSE123902 libraries (LX653, LX684) lack NCAM1 and NCR1; every effector and exhaustion gene used in the scores is present in all four cohorts. LX653 and LX684 NK compartments are below the 10-cell cutoff anyway.

## Scores

Per cell, log1p(count / library size × 10,000), library size = total counts in that cell. Per unit, the mean over cells, then the mean over genes present in that matrix. A missing score gene is dropped, not filled with zero.

| score | genes |
|---|---|
| effector_cosmx | GZMB, PRF1, NKG7, IFNG |
| exhaustion | PDCD1, HAVCR2, LAG3, TIGIT, TOX |
| naive_memory | IL7R, TCF7, CCR7, SELL, LEF1 |
| effector_no_nkg7 | GZMB, PRF1, IFNG |

Primary slice uses every cell in the compartment with library size > 0. A unit is eligible when that slice has at least 10 cells. Sensitivity: library size ≥ 500; at least 30 cells; CD8 restricted to GNLY > 0.

## Primary tests

Four tests, pre-specified:

1. CD8 effector_cosmx
2. NK effector_cosmx
3. CD8 exhaustion
4. NK exhaustion

Within each cohort, Spearman ρ of the score versus locked malignant CLDN4 percent positive. Meta-analysis is DerSimonian–Laird on Fisher z, variance 1/(n−3), cohorts with n≥5. Quartile contrast: ordinary least squares of score on cohort indicators plus a Q4 indicator, fit only on locked Q1 and Q4 units that still pass the cell filter. A stacked Mann–Whitney rank-biserial is reported beside it. p-values are two-sided and descriptive.

Under a patient-level muzzling hypothesis the effector ρ would be negative and the exhaustion ρ positive. The CosMx neighbor result did not show an effector decrease.

## Sensitivities (not primary)

- Effector score without NKG7.
- Exhaustion plus CTLA4 and ENTPD1.
- Naive/memory score in CD8.
- Within-cohort residual of CD8 effector on CD8 naive/memory, then the same Spearman meta. This asks whether an effector association is a naive-versus-effector mix.
- Cell-weighted CD8+NK effector score, and that score residualized on the NK fraction of CD8+NK. The pool confounds lineage mix with intensity.
- NK fraction of CD8+NK versus CLDN4 (composition, not intensity).
- Gene-level means for each CosMx gene and each exhaustion gene.
- Author cell-state fractions and within-state effector means where the paper published those labels (GSE131907 subtypes; GSE205335 celltypes). Those labels are not meta-analyzed together.
- Leave-one-cohort-out of the four primaries.

## Statistics check

The same DerSimonian–Laird and stacked rank-biserial code is run on the locked malignant CLDN4 percent positive versus T/NK fraction, which is already published for these 65 units (ρ about −0.53). That run is a calibration of the code, not a new result.
