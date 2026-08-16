# A3 — GSE205335 malignant TACSTD2 versus T/NK; MPR not labeled

## Result

The requested malignant-versus-T/NK contrast is feasible from open processed
GEO data. TACSTD2 was strongly enriched in author-labelled malignant cells
relative to T/NK cells in every one of 22 evaluable patients:

- mean log1p(CP10K), median: **0.929 malignant vs 0.021 T/NK**
- TACSTD2-positive cells, median: **58.4% malignant vs 1.2% T/NK**
- paired two-sided Wilcoxon signed-rank: **W = 0, p = 4.77e-7** for both
  measures

This supports malignant-compartment restriction of TACSTD2 in this cohort. It
does not establish tumor specificity outside the sampled compartments, and low
T/NK signal may include ambient RNA.

**MPR vs NMPR cannot be run.** GEO sample characteristics and the author cell
table contain RECIST 1.1 only (`PR` / `SD` / `PD` / `NE`). There is no
pathologic-response or residual-viable-tumor field. This is a palliative ICI
biopsy/effusion cohort, not a neoadjuvant resection series. RECIST is not
substituted for MPR. The scan is recorded in `mpr_gate.json`.

The secondary RECIST comparison was null: six responders versus ten
non-responders had p = 0.96 for mean normalized expression, p = 0.71 for
percent positive, and p = 0.64 for pseudobulk CPM. An ADC+SQ-only sensitivity
is underpowered and is not the verdict. These small, heterogeneous groups
cannot exclude a moderate response association.

An exploratory Spearman of malignant TACSTD2 versus T/NK fraction among
(malignant + T/NK) cells was also null (n = 22, ρ ≈ 0.18, p ≈ 0.41). That is
not the GSE207422 MPR claim and is not used as a substitute MPR test.

## Data and methods

- Dataset: [GEO GSE205335](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205335)
- Matrix: author-processed 33,714-gene x 96,505-cell UMI `dgCMatrix`
- Compartments: the authors' `lineage.sub == "Malignant cells"` and
  `lineage.total == "T/NK cells"` annotations
- Unit: patient; samples from the same patient were pooled
- Eligibility: at least 20 cells in each compartment
- Expression: per-cell `log1p(TACSTD2 UMI / total UMI x 10,000)`
- Primary test: paired, two-sided Wilcoxon signed-rank across patients
- MPR gate: string search of GEO SOFT fields, sample characteristics, and
  identity-table columns for MPR/NMPR/pathologic-response tokens

The analysis used 28,512 malignant and 39,875 T/NK cells. Marker checks agreed
with the labels: EPCAM was detected in 69.3% of malignant versus 2.5% of T/NK
cells; PTPRC was detected in 4.2% versus 80.9%, respectively.

## Honest data gate

Only the open GEO processed matrix, cell identities, and GEO SOFT metadata were
downloaded. The controlled EGA raw data (`EGAD00001008703`) were not accessed.
The 524 MB processed matrix was cached under `/tmp` and is not committed.
Checksums of all analyzed inputs are in `provenance.tsv`.

Four patients lacked enough captured malignant cells. Dropout was asymmetric
(three responders and one non-responder), so the RECIST comparison is
selection-prone. Tissues, histologies, and 10x 3'/5' chemistries are mixed, and
malignancy was not independently re-called by CNV.

## Files

- `stats_results.txt`: readable statistical report, including the MPR skip
- `summary.json`: structured results and data-policy statement
- `mpr_gate.json`: explicit MPR/NMPR unavailability scan
- `per_patient_tacstd2.csv`: patient-by-compartment metrics
- `per_sample_tacstd2.csv`: sample-by-compartment metrics
- `tnk_fraction_spearman.tsv`: exploratory malignant TACSTD2 vs T/NK fraction
- `gsm_sample_metadata.csv`: GEO sample mapping and clinical fields
- `tacstd2_summary.png`: RECIST and paired-compartment panels
- `tacstd2_per_patient_bars.png`: malignant TACSTD2 by patient
- `tacstd2_vs_tnk_fraction.png`: exploratory T/NK-fraction scatter
- `provenance.tsv`: input sizes and SHA-256 checksums

Reproduction code is in `scripts/w200/A3_GSE205335/`.
