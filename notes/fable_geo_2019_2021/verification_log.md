# Verification log — GEO 2019–2021 human lung ICI series

This log documents the exhaustive search, verification, and provenance for the
parallel slice. All GEO accessions are **real** identifiers returned by NCBI
E-utilities; none were invented. Outputs are confined to
`notes/fable_geo_2019_2021/`, `scripts/fable_geo_2019_2021/`, and
`results/fable_geo_2019_2021/`.

## 1. Search strategy (`scripts/.../01_search_geo.py`)

Database: NCBI `gds` (GEO). Filters applied to every query:
- Entry type: `"gse"[Entry Type]` (series only)
- Organism: `"Homo sapiens"[Organism]`
- Publication date: `("2019/01/01"[PDAT] : "2021/12/31"[PDAT])`

Term = (lung / NSCLC / LUAD / LUSC / SCLC / adenocarcinoma / squamous …) AND
(immune checkpoint / checkpoint inhibitor / immunotherapy / PD-1 / PD-L1 /
CTLA-4 / nivolumab / pembrolizumab / atezolizumab / durvalumab / avelumab /
ipilimumab / cemiplimab …).

Per-drug supplementary queries were also run to maximise recall. Result:
**90 unique GEO series** (`results/.../search_uids.json`,
`results/.../candidates_metadata.json`).

## 2. Triage of all 90 candidates (`scripts/.../06_triage.py`)

Full auditable table: `results/.../tables/triage_all_candidates.csv`.

| category | n | meaning |
|---|---|---|
| ANALYZED | 4 | human lung + ICI outcome + processed matrix + genes measurable |
| VERIFIED_NO_TARGET_GENES | 1 | human lung + ICI, but assay panel lacks TACSTD2/CLDN4 |
| LUNG_ICI_other | 20 | lung+ICI but no per-sample bulk-tumor outcome (blood/panel/descriptive/anti-CD4/chemo) |
| LUNG_ICI_single_cell | 16 | lung+ICI single-cell (not patient-level bulk outcome) |
| LUNG_ICI_methylation | 6 | lung+ICI methylation/epigenetic assay (no mRNA of target genes) |
| LUNG_ICI_cellline_invitro | 11 | lung+ICI cell-line / in-vitro |
| EXCLUDED_not_ICI | 11 | lung but no ICI treatment context |
| EXCLUDED_not_lung | 21 | other tumor types (melanoma/urothelial/breast/prostate/glioma), COVID lung, etc. |

## 3. Verified datasets used (downloaded, open, processed, < 2 GB)

Downloaded from the official GEO FTP `suppl/` and `matrix/` directories
(`scripts/.../03_download.py`, 2 GB per-file hard cap enforced via HTTP HEAD).

| GSE | PMID | tissue / assay | n | processed file (size) | ICI outcome in GEO | TACSTD2/CLDN4 |
|---|---|---|---|---|---|---|
| GSE126044 | 32879421 | NSCLC tumor biopsy, bulk RNA-seq | 16 | `GSE126044_counts.txt.gz` (0.56 MB) | responder / non-responder (anti-PD-1) | present (symbols) |
| GSE135222 | 31537801 | NSCLC tumor, bulk RNA-seq (TPM) | 27 | `..._omicslab_exp.tsv.gz` (1.66 MB) | PFS event + time (anti-PD-1/PD-L1) | present (Ensembl) |
| GSE182328 | — | advanced/limited lung tumor, bulk RNA-seq | 44 | `GSE182328_Gene_counts_matrix.txt.gz` (1.10 MB) | **no per-sample response/PFS in GEO**; only Akkermansia (`akk_metaominer`) + stage | present (symbols) |
| GSE111414 | 30765392 | PBMC CD8+ T cells, bulk RNA-seq | 20 | `GSE111414_gene_counts.csv.gz` (1.02 MB) | responder / non-responder (nivolumab) | present but ~0 (epithelial genes) |
| GSE136961 | 31959763 | NSCLC tumor, Oncomine Immune Response panel | 21 | `GSE136961_TPM.tsv.gz` (0.07 MB) | DCB / NDB + survival | **absent** (395-gene immune panel) |

Per-sample clinical tables parsed from `*_series_matrix.txt.gz` are in
`results/.../clinical/` (`scripts/.../04_build_clinical.py`).

## 4. Caveats

- **Small samples**: n = 16–44; every test is underpowered.
- **No multiple-testing correction** was applied to the exploratory p-values.
- **GSE182328** has no per-patient ICI response/survival deposited in GEO; the
  Akkermansia-detectable stratifier is a *surrogate* prognostic biomarker (from
  the study's endpoint), not a direct outcome. Treated as exploratory only.
- **GSE111414** measures peripheral CD8+ T cells; TACSTD2/CLDN4 are epithelial
  genes and are ~0 here — included only as a negative/QC control.
- **GSE136961** is verified as a genuine NSCLC anti-PD-1 cohort but its targeted
  panel does not measure the two genes, so it cannot contribute to the analysis.
