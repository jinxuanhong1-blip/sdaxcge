# scRNA+TCR playbook — public lung ICI (additive)

This folder **runs** the analysis. The methods-only TCR protocol lives in
`methods/tcr/` (PR 213). Here the same estimands are applied to public files
and the numbers are written under `results/`.

Public only. **Do not download EGA** (Caushi raw FASTQ:
EGAS00001005343 / EGAD00001007728).

---

## Estimands (frozen before EDA)

| ID | Estimand | Unit | Contrast |
| --- | --- | --- | --- |
| **E1** | Abundance of **expanded CXCL13+ T cells** | **patient / sample** | MPR-any (GEO `pCR`+`MPR`) vs `non-MPR` |
| **C** | Combinatorial: CXCL13+ (or Tex) rate **among expanded vs non-expanded** cells | patient (paired) | within-sample, then vs MPR |
| **R** | Same E1/C endpoints vs **residual** TACSTD2/CLDN4 detection | patient | Spearman; **not** a tumor score |
| **E2** | E1 vs **tumor-epithelial** TACSTD2/CLDN4 | patient | only if a paired epithelial/bulk/IHC table exists |

E2 is **not estimable** on any cohort in this folder.

---

## Cohorts (what is actually public)

### GSE243013 (primary)

Liu et al., *Cell* 2025 (PMID 40147443). Post-neoadjuvant chemo-IO NSCLC,
**CD45+ immune-only** scRNA + TCR. MPR labels on GEO.

| File | Use |
| --- | --- |
| `GSE243013_T_with_TCR_annotation.csv.gz` (15 MB) | E1 + C (author clusters + CDR3) |
| `GSE243013_NSCLC_immune_scRNA_metadata.csv.gz` (39 MB) | MPR, histology |
| Residual table from the public MTX (this repo) | R only |
| `GSE243013_NSCLC_immune_scRNA_counts.mtx.gz` (6.6 GB) | already used for residual TACSTD2/CLDN4; not re-downloaded |
| NGDC / GSA-Human FASTQ | **list only** |

Verified TCR clusters: `CD4T_Tfh_CXCL13`, `CD4T_Th1-like_CXCL13`.
**There is no `CD8T_Tex_CXCL13`** in the public table. CD8 Tex is
`CD8T_Tex_HAVCR2` + `CD8T_terminal_Tex_LAYN` (secondary).

### Caushi / GSE176022 (if GEO processed exists)

Caushi et al., *Nature* 2021 (PMID 34290408). SuperSeries GSE173351.

| File | What it is | Use |
| --- | --- | --- |
| **GSE176022_RAW.tar** (10.9 MB) | **GEO processed exists.** MiXCR **bulk** TCR from MANAFEST / virus cultures. GEO title is bulk TCR-seq, not scRNA. | leftover TCR catalog; no CXCL13, no TACSTD2 |
| `GSE176021_CD3/CD8_annotations.rds.gz` | CellType + UMAP only | Tfh **proxy**, not CXCL13 RNA |
| per-sample `*.vdj.tar.gz` on GEO | processed 10x VDJ (double-gzipped tar) | clone expansion vs MPR |
| `GSE176021_RAW.tar` (4.1 GB GEX) | processed GEX; not required for E1 | not downloaded |
| EGA EGAS00001005343 / EGAD00001007728 | raw FASTQ | **skip** |

### Leftover open TCR

| Accession | Use |
| --- | --- |
| **GSE179994** (Liu *Nat Cancer* 2022) | public `all.scTCR.tsv.gz` + T-cell metadata; Tex cluster; **no response/MPR on GEO** |
| GSE185204 | n=3 post-ICB lung regions; ~900 MB; **listed only** |
| GSE236581, GSE162025 | not lung ICI; excluded |

---

## Definitions

**Clone key (primary):** amino-acid `TRA_cdr3|TRB_cdr3` within `sampleID`.
Never merge CDR3s across patients.

**Expanded (primary):** clone size ≥ 2 in that sample. Sensitivity: size ≥ 3;
author `expansion` flag on GSE243013.

**CXCL13+ (GSE243013 primary):** author `sub_cell_type` ∈
`{CD4T_Tfh_CXCL13, CD4T_Th1-like_CXCL13}`.
RNA gate is not used (6.6 GB MTX not re-streamed for CXCL13).

**CXCL13+ (Caushi):** not estimable from the public RDS (no expression).
`CellType` containing `Tfh` is a **proxy secondary**.

**Minimum depth:** drop samples with `< 50` TCR-qualified cells. Report n.

**Unit:** patient/sample. Wilcoxon / Spearman on patient rows only.
Do not Wilcoxon cells labeled MPR.

**Residual TACSTD2/CLDN4:** immune-compartment detection fraction / mean CPM
from the CD45+ MTX. Ambient / residual-tumor RNA is a live explanation.
**Not E2.**

---

## Primary endpoint

`exp_cxcl13_per_k_cd8` =
`n_expanded_CXCL13+_clones / n_CD8_TCR_cells * 1000`.

Combinatorial endpoints (pre-registered):

- `frac_cxcl13_among_expanded` vs `frac_cxcl13_among_nonexpanded` (paired Wilcoxon)
- same for CD8 Tex
- patient-level logOR of CXCL13+ given expanded vs not, then vs MPR

---

## What this playbook will not do

- Download EGA or NGDC FASTQ.
- Treat residual / ambient TACSTD2 as tumor-epithelial E2.
- Use cells as MPR replicates.
- Invent GSE179994 response labels that are not on GEO.
- Equate RECIST with MPR.
- Claim a `CD8T_Tex_CXCL13` cluster that is absent from the public table.
