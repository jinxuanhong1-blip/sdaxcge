# 2024 GEO lung ICI screen: TACSTD2 / CLDN4

## Bottom line

The pair is convincingly co-expressed as an epithelial/tumor program, but these
2024-public GEO data do **not** establish TACSTD2 or CLDN4 as an ICI-response
biomarker. The leftover residual-tumor series GSE241934 repeats the same
pattern: strong co-expression in leftover epithelium, no significant
MPR/pCR association.

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

The first pass examined four records numerically:

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
cannot test the pair: their targeted panels contain TACSTD2 but not CLDN4,
and the leftover xlsx files also lack a response column. PBMC, miRNA,
cfDNA-hydroxymethylation, and mouse series were not treated as
tumor-expression validation.

## Leftover 2024 series

The leftover pass finished human lung IO series that the first pass only
classified. It does not change the bottom line. The pair remains an
epithelial/tumor program. These leftover files still do not establish either
gene as an ICI-response biomarker.

### GSE241934 — residual tumor after neoadjuvant PD-1 plus chemotherapy

Public 2024-05-17. This is the highest-value leftover: post-resection 10x
scRNA-seq from [Zhang et al., Cell Rep Med 2024](https://doi.org/10.1016/j.xcrm.2024.101615)
(NEOTIDE/CTONG2104). Two cohorts were kept separate.

- **IIT:** 11 EGFR-mutant patients, all sintilimab × 3 cycles, 1,699
  epithelial cells. All 11 patients have ≥32 epithelial cells (4 MPR / 7
  non-MPR).
- **Real/RWC:** 34 EGFR-WT patients, mixed PD-1 agents, 12,785 epithelial
  cells in 33 patients (one patient has zero epithelial cells). Many
  MPR/pCR residuals have 1–17 leftover epithelial cells. The primary filter
  keeps patients with ≥20 epithelial cells (24 patients: 6 MPR/pCR, 18
  non-MPR). A stricter ≥50-cell filter leaves 20 patients (5 / 15).

This is a **residual-tumor** readout, not a pretreatment biomarker test.
MPR/pCR patients have few remaining epithelial cells by definition. IIT and
Real were not pooled.

Patient-level epithelial-cell log2(CPM + 1) of TACSTD2 versus CLDN4:

| Cohort | Filter | n | Spearman rho |
|---|---|---:|---:|
| IIT EGFR-mutant | ≥20 epithelial cells | 11 | 0.936 |
| Real EGFR-WT | ≥20 epithelial cells | 24 | 0.342 |
| Real EGFR-WT | ≥50 epithelial cells | 20 | 0.223 |

Exact permutation tests of MPR/pCR versus non-MPR (means shown; p enumerates
all assignments):

| Cohort / filter | Metric | MPR/pCR mean | Non-MPR mean | Exact p |
|---|---|---:|---:|---:|
| IIT ≥20 | TACSTD2 log2(CPM + 1) | 9.840 | 9.385 | 0.261 |
| IIT ≥20 | CLDN4 log2(CPM + 1) | 9.861 | 9.319 | 0.194 |
| IIT ≥20 | TACSTD2 detection | 0.806 | 0.815 | 0.900 |
| IIT ≥20 | CLDN4 detection | 0.861 | 0.837 | 0.488 |
| Real ≥20 | TACSTD2 log2(CPM + 1) | 9.602 | 9.286 | 0.399 |
| Real ≥20 | CLDN4 log2(CPM + 1) | 9.274 | 9.641 | 0.316 |
| Real ≥50 | TACSTD2 log2(CPM + 1) | 9.657 | 9.185 | 0.264 |
| Real ≥50 | CLDN4 log2(CPM + 1) | 9.552 | 9.561 | 0.980 |

IIT co-expression is stronger than the GSE205335 pretreatment core. The
response contrasts remain small-n and non-significant. Real-world
correlation is weaker, which is expected from mixed agents and from
filtering leftover epithelium after a good pathologic response.

### Other leftover numerical audits

- **GSE253564** (public 2024-01-18): pretreatment tumor FPKM, durvalumab ±
  SBRT, n = 32. Median TACSTD2 79.0, CLDN4 52.8, Spearman rho = 0.687. GEO
  has Arm1/Arm2 only (16 vs 16); no MPR or radiographic labels. Arm
  medians differ (TACSTD2 69.4 vs 93.9; CLDN4 50.7 vs 59.4) but that is
  not a response test.
- **GSE255144** (public 2024-12-31): one-patient hyperprogression culture
  model, 3 baseline (ADK17) versus 3 HPD (ADK18) replicates. Median
  TACSTD2 falls from 564.7 to 5.3; CLDN4 from 133.4 to 51.3. This is a
  cell-plasticity experiment, not a clinical cohort. n = 3 Spearman values
  are not reported as biology.
- **GSE260598** leftover lung AOIs (cases 43728 and 43881 only): tumor-nest
  CD68-negative AOIs have higher signal and pair correlation (median
  TACSTD2 17.8, CLDN4 3.52, rho = 0.771, n = 6) than stroma CD68-negative
  AOIs (2.24, 1.30, rho = 0.086, n = 6). Compartment leftover only; not an
  ICI-outcome cohort.
- **GSE270711** (public 2024-09-30): 3 surgical LUAD tumors and 3 adjacent
  samples, not ICI-treated. Tumor medians TACSTD2 47.3 / CLDN4 2.25.
  Spearman rho = 1.0 at n = 3 is not interpretable.

### Leftover series that remain unquantifiable or the wrong compartment

- **GSE271689:** first-line IO GeoMx WTA, 586 AOIs. GEO deposits raw DCC
  only. No processed counts and no response labels. Honest status:
  unquantifiable from the GEO deposit.
- **GSE280232:** neoadjuvant ICB scRNA-seq. Most samples are sorted T
  cells; mixed CD45-positive/negative samples have no cell annotation.
  Not an epithelial pair test.
- **GSE266219:** leftover PBMC/CD8 scRNA-seq. No epithelial tumor cells.

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
  downloads the smaller first-pass GEO matrices into a temporary directory.
- `analyze_gse241934.py` streams only the TACSTD2 and CLDN4 rows from the
  leftover 10x MTX files (IIT 445 MB; Real 1.2 GB) and writes patient-level
  residual-epithelium summaries. Standard library only.
- `analyze_leftover_series.py` audits GSE253564, GSE255144, GSE260598 lung
  AOIs, and GSE270711.
- `data/GSE205335_core_sample_targets.tsv` — 14 sample-level pseudobulks.
- `data/GSE205335_core_patient_targets.tsv` — 11 patient-level values.
- `data/GSE205335_statistics.tsv` — correlation and exact response tests.
- `data/other_series_target_values.tsv` and
  `data/other_series_summary.tsv` — first-pass audits of the other
  measurable series.
- `data/GSE241934_IIT_patient_targets.tsv` — 11 IIT residual-epithelium
  patients.
- `data/GSE241934_Real_patient_targets.tsv` — 33 Real patients with any
  epithelial cells.
- `data/GSE241934_statistics.tsv` — leftover IIT/Real correlation and
  exact MPR tests.
- `data/leftover_series_target_values.tsv` and
  `data/leftover_series_summary.tsv` — leftover bulk/GeoMx audits.

Run from the repository root:

```bash
Rscript results/w200/GEO_2024/analyze_gse205335.R results/w200/GEO_2024/data
python3 results/w200/GEO_2024/analyze_other_series.py
python3 results/w200/GEO_2024/analyze_leftover_series.py
python3 results/w200/GEO_2024/analyze_gse241934.py
```

No cross-series pooled effect was computed because pretreatment malignant
pseudobulks, residual-tumor epithelial pseudobulks, RMA array signal, blood
FPKM/counts, tumor FPKM, culture counts, and GeoMx AOIs are not a common
measurement scale.
