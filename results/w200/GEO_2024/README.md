# 2024 GEO lung ICI screen: TACSTD2 / CLDN4

## Bottom line

The pair is convincingly co-expressed as an epithelial/tumor program, but these
2024-public GEO data do **not** establish TACSTD2 or CLDN4 as an ICI-response
biomarker.

The cleanest dataset is [GSE205335](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205335).
Across the 11 independent core patients, malignant-cell pseudobulk expression
of TACSTD2 and CLDN4 was strongly correlated (Spearman rho = 0.827,
p = 0.00168). Neither gene separated the four responders from the seven
non-responders:

| Patient-level metric | Responder median | Non-responder median | Exact permutation p |
|---|---:|---:|---:|
| TACSTD2 log2(CPM + 1) | 9.976 | 9.220 | 0.830 |
| CLDN4 log2(CPM + 1) | 8.904 | 9.433 | 0.733 |
| TACSTD2 cell detection fraction | 0.732 | 0.559 | 0.103 |
| CLDN4 cell detection fraction | 0.735 | 0.587 | 0.385 |

These are small, exploratory comparisons—not negative validation. There are
only 11 patients, tissue sites and histologies vary, and one included sample
has only 27 annotated malignant cells.

## What was actually screened

“2024” was defined as **GEO public-release date from 2024-01-01 through
2024-12-31**, not publication year. “Leftover series” is not GEO terminology,
so the screen retained human lung ICI series that measure gene expression and
explicitly documented near-misses. The decisions are in
[`data/series_screen.tsv`](data/series_screen.tsv).

Four records were examined numerically:

1. **GSE205335 (direct, scRNA-seq):** 33 samples from 26 ICI-treated patients.
   The outcome analysis follows the paper's prespecified 14 pretreatment core
   samples (12,975 malignant cells, 11 patients). Counts were summed within
   each sample, converted to malignant-cell log2(CPM + 1), and multiple samples
   from one patient were averaged before inference. Response p-values enumerate
   all 330 assignments of four responder labels among 11 patients.
2. **[GSE248249](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248249)
   (direct, tumor array):** TACSTD2 and CLDN4 were moderately correlated at
   baseline (rho = 0.665, n = 13) and after acquired resistance (rho = 0.552,
   n = 29). Among 13 paired patients, median post-minus-pre RMA signal was
   -0.161 for TACSTD2 (6/13 increased) and -0.885 for CLDN4 (2/13 increased).
   This is descriptive: all post samples represent acquired resistance, sites
   can differ, and GEO notes that the deposited SST-RMA values are not the
   authors' final ComBat/LOESS-normalized analysis matrix.
3. **GSE225620 and GSE260770 (direct treatment, wrong compartment):** both are
   blood datasets. In GSE225620, median raw counts were 4 for TACSTD2 and 0 for
   CLDN4. In GSE260770, median FPKM was 0 for both genes in both response
   groups. These are failed tumor-marker readouts, not biological evidence
   against the pair.
4. **[GSE265899](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE265899)
   (context only):** the discovery tumors were stratified by PD-L1 and immune
   context, not by ICI outcome. Tumor AOIs had higher median signal
   (TACSTD2 177.7; CLDN4 106.3) and greater pair correlation (rho = 0.455,
   n = 48) than immune AOIs (24.3; 17.9; rho = 0.132, n = 47). This confirms
   compartment specificity, not ICI prediction.

The two 2024 SCLC chemoimmunotherapy GeoMx series (GSE261345 and GSE261348)
cannot test the pair: their targeted panels contain TACSTD2 but not CLDN4.
PBMC, miRNA, cfDNA-hydroxymethylation, cell-line, and mouse series were not
treated as tumor-expression validation.

## The 2024 CLDN4 paper is not an ICI GEO series

The paper that explicitly reports a positive CLDN4–TACSTD2 relationship is
[Zhang et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11022306/), a
malignant-pleural-effusion study. Its 16 baseline scRNA-seq samples were from
patients subsequently receiving mutation-directed therapy or platinum
chemotherapy—not ICI. The data are controlled access at
[HRA006761](https://ngdc.cncb.ac.cn/gsa-human/browse/HRA006761), not GEO.

The paper found CLDN4 correlated with TACSTD2 in metastatic cancer cells, but
in its separate 64-patient qPCR cohort only CLDN4 differed between recurrent
and non-recurrent effusions (p = 0.027); TACSTD2 did not. That result supports
CLDN4 as a candidate recurrence marker in MPE. It should not be relabeled as
evidence for ICI response.

## Reproducibility and files

- `analyze_gse205335.R` downloads the 500 MB GEO RDS and annotation, verifies
  the published 12,975-cell core subset, and writes patient/sample summaries.
  It requires R and the Matrix package.
- `analyze_other_series.py` uses only the Python standard library and
  downloads the smaller GEO matrices into a temporary directory.
- `data/GSE205335_core_sample_targets.tsv` — 14 sample-level pseudobulks.
- `data/GSE205335_core_patient_targets.tsv` — 11 patient-level values.
- `data/GSE205335_statistics.tsv` — correlation and exact response tests.
- `data/other_series_target_values.tsv` and
  `data/other_series_summary.tsv` — audits of the other measurable series.

Run from the repository root:

```bash
Rscript results/w200/GEO_2024/analyze_gse205335.R results/w200/GEO_2024/data
python3 results/w200/GEO_2024/analyze_other_series.py
```

No cross-series pooled effect was computed because scRNA UMI pseudobulks, RMA
array signal, blood FPKM/counts, and GeoMx Q3-normalized AOIs are not a common
measurement scale.
