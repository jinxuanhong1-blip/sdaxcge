# Methods — CXCL13+ trio, CLDN4-only CellChat-style outgoing LR

Additive folder. The CXCL13+ trio combo in `methods/scrna_cldn4_combo`
(GSE148071 + GSE207422 + GSE253013, n=60, ρ=−0.425) is **taken as given**
and is not re-audited. Prior single-cohort CellChat folders
(`methods/gse148071_cellchat_cldn4`, `methods/scrna_cellchat_cldn4`) are
not re-run as the answer.

**CLDN4 only.** TACSTD2 is not used to define groups. No dual-high
TACSTD2×CLDN4 gate. This is not a GSE207422-only analysis.

## Cohorts

Identified from the prior combo page (`tls/cxcl13pos/mean`):

| Cohort | Citation | Public files |
| --- | --- | --- |
| GSE148071 | Wu et al., *Nat Commun* 2021, PMID 33953163 | `GSE148071_RAW.tar` (42 Singleron UMI matrices) |
| GSE207422 | Hu et al. NSCLC neoadjuvant PD-1 | `GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz` + sample sheet |
| GSE253013 | Dysregulated NOTCH3 LUAD atlas, PMID 38335304 | `GSE253013_all_luad_garnett_temp.rds.gz` (streamed gene panel; full Seurat load does not fit 16 GB) |

GSE131907 / GSE205335 are **not** this trio. GSE207422 is one member, not the only cohort.

## Lineage and receivers

Canonical marker scores on `log1p(CP10k)`; assigned lineage = argmax
(same panel as the GSE148071 / GSE207422 CellChat folders). T vs NK
broken by CD3E when scores are close.

- **Malignant / sender** = epithelial lineage (putative; marker-argmax) on
  GSE148071 and GSE207422. GSE253013 uses author **Epithelial** (tumor),
  not Airway Epithelium. Adjacent lung (`ANT`/`NAT`) is dropped.
- **T/NK** = T or NK, merged.
- **CXCL13+ T** = T lineage **and** CXCL13 UMI > 0.

## CLDN4-only gate (patient-level)

Among epithelial cells **of that patient**, split `log1p(CP10k)` CLDN4
at the **within-patient median** (high = above median, low = at/below).
TACSTD2 is ignored.

A patient is eligible for Mal→T/NK if it has ≥10 CLDN4-high, ≥10
CLDN4-low, and ≥15 T/NK cells. Eligible for Mal→CXCL13+ T if the same
malignant arms plus ≥8 CXCL13+ T cells. Patients below the floor are
reported in `results/patient_n.tsv` and are not scored. That n is **not**
the given Spearman n=60.

## CellChat-style probability

CellChat R was not installed. Probability follows Jin et al. 2021 on
CellChatDB v2 protein pairs:

1. 10% truncated mean of `log1p(CP10k)` per gene; complexes = geometric
   mean of subunits (0 if any subunit mean is 0).
2. Detected if each complex has expressing-cell fraction ≥ 0.10 (AND).
3. \(P = (L \cdot R) / (K_h + L \cdot R)\) with \(K_h = 0.5\).

Outgoing only: Mal_high → receiver and Mal_low → receiver.
ΔP = P_high − P_low per patient.

## Patient-level test

For each LR pair and receiver, Wilcoxon signed-rank on per-patient ΔP
among patients where the pair is detected on at least one arm (n≥6).
BH-FDR within receiver. p-values are descriptive. Cell-pooled label
permutation is **not** the primary claim.

## Software

Python (numpy / pandas / scipy / matplotlib). CellChatDB v2 parsed from
the public `CellChatDB.human.rda`. Raw FASTQ was not used. CellChat R
and LIANA were not run.
