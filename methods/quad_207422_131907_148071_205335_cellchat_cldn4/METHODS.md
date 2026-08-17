# Methods — QUAD GSE207422+GSE131907+GSE148071+GSE205335 CLDN4 CellChat

ADDITIVE. **CLDN4 only.** TACSTD2 does not define groups. Patient is the unit.
GSE207422 is **in** the merge. GSE207422-only CLDN4 vs T/NK and PR #320
(GSE131907+GSE205335 Q4) are given and are not re-audited.

## Datasets (public processed UMIs)

| cohort | accession | processed file | skipped |
|---|---|---|---|
| GSE207422 | Hu et al., *Genome Med* 2023 | GEO UMI matrix + sample sheet | author Seurat/CopyKAT (not public) |
| GSE131907 | Kim et al., *Nat Commun* 2020 | GEO raw UMI + author annotation | 2.86 GB log2TPM; EGA FASTQ |
| GSE148071 | Wu et al., *Nat Commun* 2021 | GEO per-biopsy `*_exp.txt.gz` | no author cell labels on GEO |
| GSE205335 | Hu et al. lung IO scRNA | GEO `dgCMatrix` RDS + identity | EGA raw `EGAD00001008703` |

Any file >2 GB is skipped.

## Malignant and T/NK

Eligibility: **≥20 malignant and ≥20 T/NK cells in the same patient**.

- **GSE207422:** post-treatment resections only. Epithelial = marker-argmax
  (EPCAM/KRTs/CDH1); T/NK = T or NK marker-argmax. CopyKAT IDs are not public.
- **GSE131907:** author `Cell_type == Epithelial cells` and subtype in
  `{Malignant cells, tS1, tS2, tS3}`; T/NK = author T lymphocytes or NK.
  Tumor-origin only (`tLung`, `tL/B`, `mLN`, `mBrain`). Cells from one
  patient are pooled (patient is the unit, not sample).
- **GSE148071:** same marker-argmax rule as GSE207422; one biopsy per patient.
- **GSE205335:** author `lineage.sub == Malignant cells` and
  `lineage.total == T/NK cells`.

## Combo rho and Q4 vs Q1

Primary CLDN4 score = % of malignant cells with UMI > 0. Secondary = mean
`log1p(CP10k)` (written on the patient table, not used to rank the merge).

- Spearman ρ on the **pooled eligible patients** (honest n).
- Leave-one-cohort-out: drop each cohort and recompute the pooled Spearman.
- Fisher-z meta of the four within-cohort Spearmans (I² reported).
- Q4 vs Q1: `pd.qcut` on average ranks **within each cohort**, then
  Mann–Whitney on the pooled tails. Rank-biserial \(r = 2U/(n_4 n_1)-1\);
  \(r<0\) means Q4 has lower T/NK. Thin if n<8 or a tail <3.
- Global quartiles are written on the patient table as a companion column
  only. p-values are descriptive.

## CellChat-style ligands

CellChat R and LIANA are not run. Probability follows Jin et al. 2021 on
CellChatDB v2 protein pairs:

1. Keep a pair only if every ligand and receptor subunit is in that cohort's UMI matrix.
2. Per patient / compartment: 10% truncated mean of `log1p(CP10k)`; complexes =
   geometric mean of subunits (0 if any subunit mean is 0).
3. Detected if each complex has expressing-cell fraction ≥ 0.10 (AND rule).
4. \(P = (L \cdot R) / (K_h + L \cdot R)\) with \(K_h = 0.5\).
5. Outgoing = malignant → that patient's T/NK. Incoming = T/NK → malignant.
6. Q4 vs Q1 is Mann–Whitney on per-patient *P* using the within-cohort
   CLDN4 %pos tails. A pair is listed if detected in ≥3 Q4 and ≥3 Q1 patients.

Cell-pooled (stacked) truncated means are not the test.

## Software

Python (numpy / pandas / scipy / matplotlib / rdata). CellChatDB v2 from
`jinworks/CellChat`. Raw GEO matrices are not stored in git.
