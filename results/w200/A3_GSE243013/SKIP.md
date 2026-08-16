# W200-A3 · GSE243013 — SKIPPED (honest skip)

**Task:** A3-analog analysis — TACSTD2 expression in malignant cells, MPR vs non-MPR — on
GSE243013 ("A single-cell atlas of immune heterogeneity in anti-PD1-treated non-small cell
lung cancer"; 234 NSCLC patients, 243 GEO samples, post-neoadjuvant chemo-immunotherapy),
conditional on the processed data being < 2 GB.

**Decision: SKIP.** Two independent blockers, each sufficient on its own:

1. **Size budget exceeded.** The only gene-expression file on GEO is
   `GSE243013_NSCLC_immune_scRNA_counts.mtx.gz` at **6.63 GB gzipped** (7,123,039,063 bytes),
   more than 3× the 2 GB budget. There is no per-sample or gene-subset expression alternative:
   `GSE243013_RAW.tar` (516 MB) contains only per-patient **TCR** tarballs (no expression), and
   even streaming just the TACSTD2 row would still require transferring the full 6.6 GB matrix.

2. **No malignant cells in the dataset.** This is a CD45+/immune-focused atlas. The full cell
   metadata (39 MB, 1,254,749 cells) contains exactly three major cell types — T/NK (766,574),
   B (297,076), Myeloid (191,099) — and none of its 51 sub-cell types is epithelial, malignant,
   tumor, or cancer. Malignant-cell TACSTD2 therefore cannot be measured from this dataset at
   any download size, so relaxing the budget would not rescue the endpoint.

## What was verified (all within budget; ~40 MB downloaded total)

- HTTP HEAD sizes of all 8 supplementary files → `file_manifest.tsv`
- Full cell metadata (`GSE243013_NSCLC_immune_scRNA_metadata.csv.gz`, 39 MB) → cell-type census
  in `cell_type_composition.tsv`; zero epithelial/malignant cells.
- Gene index (`GSE243013_genes.csv.gz`, 112 KB): TACSTD2 **is** in the matrix gene index, but its
  expression is locked inside the 6.6 GB matrix (and would be immune-cell TACSTD2, not malignant).
- Outcome labels: `pathological_response` per patient is available (pCR 85, MPR 45, non-MPR 112,
  unknown 1) → `patient_response_summary.tsv`. So the cohort would have supported an MPR contrast
  had a malignant compartment and an in-budget matrix existed.
- Machine-readable verdict: `feasibility.json`.

## Reproduce

```bash
python3 scripts/w200_a3_gse243013/feasibility_check.py
```

## Note for future tasks

GSE243013 remains valuable for **immune-compartment** questions vs MPR (e.g. myeloid/T-cell
programs; NMF patient groups in `GSE243013_NMF_all_group_5.csv.gz`), but any expression-level
analysis requires accepting the 6.6 GB matrix download. It is not usable for tumor-cell-intrinsic
TACSTD2/TROP2 questions.
