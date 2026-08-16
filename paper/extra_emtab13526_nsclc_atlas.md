# Extra scRNA n · E-MTAB-13526 Cvejic NSCLC atlas (not ICI)

**Placement.** Extra NSCLC TME scRNA panel. Do **not** replace or re-analyze
user A3 (GSE207422 malignant TACSTD2 vs T/NK).

## Honest cohort identity

E-MTAB-13526 (De Zuani, Xue, Park et al., *Nat Commun* 2024, PMID 38782901;
Cvejic lab) is a public treatment-naive NSCLC scRNA atlas: tumor and matched
non-involved lung from 24 deposited patients (paper n=25) plus 2 healthy
donors. ArrayExpress SDRF has disease, FACS fraction, sampling site, and TNM.
**No ICI / MPR / pCR / R / NR labels.** This is extra n for the malignant
TACSTD2/CLDN4 vs T/NK figure, not a response cohort.

The author-annotated h5ads (`10X_Lung_Tumour_Annotated_v2.h5ad` 58.7 GB;
background/healthy 45.5 GB) were not used. Analysis used the deposited
processed 10x Cell Ranger matrices for **CD235a− tumor** lanes only (RBC-depleted,
not CD45-sorted): 15 lanes, 12 patients, 276,018 cells after author-like QC
(UMI 400–100,000, genes 180–6,000, mito ≤ 20%). TACSTD2 and CLDN4 are both
present in the features table. CD45+/MDSC-sorted tumor lanes were omitted
because they lack malignant epithelium and inflate T/NK fractions.

## Results (patient unit)

Malignant-like = marker epithelial minus normal-lung markers (*SFTPA2*,
*AGER*, *SCGB1A1*, *SCGB3A1*, *TPPP3*). Eligible: ≥10 malignant-like and
≥20 T/NK cells. P16 had 2 epithelial cells; P22 had 969 epithelial but 1
malignant-like cell.

| Contrast | n | ρ | p |
|---|---:|---:|---:|
| Malignant-like TACSTD2 mean log1p UMI vs T/NK fraction | 10 | +0.45 | 0.19 |
| Malignant-like TACSTD2 mean log1p(CP10k) vs T/NK | 10 | +0.16 | 0.65 |
| Malignant-like TACSTD2 %pos vs T/NK | 10 | +0.27 | 0.45 |
| All-epithelial TACSTD2 mean log1p vs T/NK | 11 | +0.091 | 0.79 |
| All 12 CD235a− patients, malig-like TACSTD2 vs T/NK | 12 | −0.025 | 0.94 |
| Malignant-like TACSTD2 vs T/NK, n_malig ≥ 50 | 9 | +0.50 | 0.17 |
| Malignant-like TACSTD2 vs T/NK, LUSC eligible | 5 | +0.70 | 0.19 |
| Malignant-like TACSTD2 vs T/NK, LUAD eligible | 3 | — | too few |
| Malignant-like CLDN4 mean log1p UMI vs T/NK | 10 | +0.42 | 0.23 |
| Malignant-like CLDN4 mean log1p(CP10k) vs T/NK | 10 | +0.14 | 0.70 |
| All-epithelial CLDN4 mean log1p vs T/NK | 11 | +0.055 | 0.87 |
| Response / MPR / R vs TACSTD2 or CLDN4 | 0 | — | no public labels |

Both genes are epithelial-restricted (eligible n=10, paired Wilcoxon):
median %pos malignant-like TACSTD2 45.5 vs T/NK 2.71 (p=0.00098); CLDN4
32.6 vs 1.79 (p=0.00098).

## Suggested results text

In an independent public treatment-naive NSCLC scRNA atlas (E-MTAB-13526;
Cvejic / De Zuani *Nat Commun* 2024; CD235a− tumor lanes, 276,018 QC cells,
12 patients; **not** an ICI-response cohort), malignant-like *TACSTD2* was
not inversely associated with the per-patient tumor T/NK fraction (Spearman
ρ=+0.45, p=0.19, n=10 eligible). *CLDN4* was likewise non-significant
(ρ=+0.42, p=0.23). Including the two patients with almost no malignant-like
cells (n=12) moved both correlations to approximately zero. No major-pathologic-
response or radiographic-response labels are released with this accession.
Lineages are marker-based; author cell-type labels live only in the unused
58.7 GB annotated tumour h5ad.

## Methods (paper methods)

Processed 10x matrices were downloaded from ArrayExpress/BioStudies
E-MTAB-13526 (public FTP). Unfiltered Cell Ranger 3.1.0 mtx files for
CD235a− tumor lanes were streamed to a lineage/target gene panel
(`methods/emtab13526_tacstd2/`). Author QC gates were applied. Lineages
were assigned by argmax of mean log1p marker scores. Malignant-like
epithelium was epithelial cells with near-zero normal-lung markers.
Per-patient T/NK fraction used T (*CD3D*/*CD3E*/*CD2*) plus NK
(*NKG7*/*GNLY*/*FGFBP2*). Associations used Spearman rank correlation on
patients. GSE207422 was not re-analyzed.

## Figure

`methods/emtab13526_tacstd2/fig_emtab13526_tacstd2_cldn4_vs_tnk.png`

**Extra scRNA n.** Malignant-like *TACSTD2* and *CLDN4* versus per-patient
T/NK fraction in E-MTAB-13526 treatment-naive NSCLC (CD235a− tumor; 12
patients, 10 eligible). Not ICI. Spearman, patient unit: *TACSTD2* n=10,
ρ=+0.45, p=0.19; *CLDN4* n=10, ρ=+0.42, p=0.23.
