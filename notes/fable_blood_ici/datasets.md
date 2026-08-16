# fable_blood_ici — dataset provenance & methods

Scope: detectability and ICI-response/irAE association of **TACSTD2** (TROP2) and
**CLDN4** in the *blood compartment* of human NSCLC patients. These are epithelial /
tumour-cell surface genes (antibody–drug-conjugate targets), so the biological prior
is that they should be low or absent in circulating immune cells; the point of this
slice is to quantify that honestly on real GEO data.

All raw downloads go to `$GEO_DIR` (default `/tmp/geo`) and are **not** committed.
Every processed file downloaded is < 2 GB.

## Datasets used

### GSE285888 — baseline PBMC scRNA-seq (primary, named in task)
- Title: *Single-Cell RNA Sequencing of Baseline PBMCs Predicts ICI efficacy and
  irAE Severity in NSCLC Patients*. PMID 40404203 (2025).
- Platform: 10x Genomics, Illumina NovaSeq (GPL16791). Homo sapiens.
- Content: one processed dense matrix `GSM8712033_matrix.txt.gz` (~486 MB,
  28,911 genes × 222,144 cells) + `GSM8712033_metadata.csv.gz` (per-cell
  annotation). Values are raw UMI counts.
- Labels (per patient, from `orig.ident` / `SubType`): **CR** complete response
  (7 pts), **DR** durable response (9 pts), **PD** progressive disease (13 pts),
  and 5 patients labelled only by irAE. irAE severity column: non / mild / severe.
- Why: only processed data is offered; matrix < 2 GB. Contains genuine ICI
  response and irAE labels — exactly the axis the task asks about.

### GSE305086 — whole-blood bulk microarray (complementary)
- Title: *Blood gene expression under immunotherapy as potential biomarker for
  immune checkpoint blockade in non-small-cell lung cancer*.
- Platform: Affymetrix HG-U133 Plus 2.0 (GPL570). PAXgene whole blood.
- Content: `GSE305086_Expression_matrix_final.csv.gz` (~25 MB, 26,452 probes ×
  173 samples; `;`-separated; log2 RMA-normalised intensities). 76 baseline +
  76 follow-up (after 4 IO cycles) + 21 age-matched controls (10 non-smoker,
  11 smoker).
- Labels available in GEO: treatment regimen only (1LIO / 2LIO / CHTIO) plus
  sample type (baseline / follow-up / control). No per-patient response/PFS
  label is public, so this dataset is used for **detectability in bulk whole
  blood** and disease/treatment contrasts, not for response prediction.
- Why: a bulk platform is the fairer place to ask "is the target detectable at
  all in blood", because array probes register low-level transcripts that a
  droplet scRNA assay drops to zero.

## Gene panel
- Targets: TACSTD2, CLDN4.
- Epithelial reference (expected near-background in blood): EPCAM, KRT8, KRT18.
- Immune lineage reference (expected high in blood): PTPRC/CD45, CD3D, CD8A,
  MS4A1, CD14, NKG7.
- Effector / study-highlighted: GZMB, PRF1, IL1B, CXCL8 (GSE285888);
  GATA3, PDCD1 as T-cell positive controls (GSE305086; FOXP3 probes are absent
  from the filtered array matrix).

Probe→symbol mapping for GPL570 was taken from the official GEO platform table
(`acc.cgi?acc=GPL570&targ=self&form=text&view=full`), restricted to the panel.

## Statistics
- **GSE285888.** Detection rate = fraction of cells with UMI count > 0
  (normalisation-independent). Magnitude = CP10K (count / nCount_RNA × 1e4).
  Association tests use **patient-level pseudobulk** (one value per patient) so
  the statistical unit is the patient, not the cell; groups compared with the
  two-sided **Mann-Whitney U** test. Contrasts: CR vs PD, responder(CR+DR) vs PD,
  and severe-irAE vs non.
- **GSE305086.** Detectability = percentile rank of each probe's mean intensity
  among all 26,452 probes (left tail = array background = effectively absent).
  Disease effect: baseline vs control (**Mann-Whitney U**). Treatment effect:
  paired baseline vs follow-up per patient (**Wilcoxon signed-rank**).

## Pipeline validity (positive controls)
The pipeline recovers the source studies' own real signals, confirming the
null result for the targets is not a processing artefact:
- GSE285888: IL1B and CXCL8 are significantly higher in severe-irAE vs non-irAE
  patients (p = 0.026 each) — the study's reported irAE predictors.
- GSE305086: GATA3 is strongly reduced in patient baseline blood vs controls
  (p = 1.3e-10) — consistent with the study's reported T-cell gene depletion.

## Reproduce
```bash
export GEO_DIR=/tmp/geo
bash   scripts/fable_blood_ici/01_download.sh
python scripts/fable_blood_ici/02_extract_gse285888.py
python scripts/fable_blood_ici/03_analyze_gse285888.py
python scripts/fable_blood_ici/04_analyze_gse305086.py
```
