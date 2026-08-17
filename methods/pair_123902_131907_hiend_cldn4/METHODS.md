# Methods — GSE123902 + GSE131907 CLDN4-only high-end CellChat / LIANA

ADDITIVE. **CLDN4 only.** No dual-high TACSTD2×CLDN4. Patient/donor is
the unit (GSE131907 sample). Do not add GSE148071. Do not pile the
seven-cohort pool.

## Locked combo (not re-audited)

PR #459 pair that **differs**:

| combo | score | N | ρ (p, I²) | Q4 vs Q1 r (p) |
|---|---|---:|---|---|
| GSE123902+GSE131907 | %pos | 34 | −0.575 (0.001, 0%) | −0.700 (0.012; 10 vs 8) |

| unit | malignant def | unit of analysis | n |
|---|---|---|---:|
| GSE123902 | marker-malignant (EPCAM\|KRT8\|KRT18\|KRT19 > 0 and PTPRC == 0) | donor | 13 |
| GSE131907 | author `Malignant cells` (not tS*) | sample | 21 |

T/NK: GSE123902 = (CD3D\|CD3E\|CD8A\|NKG7\|GNLY\|KLRD1) > 0 and not
malignant; GSE131907 = author `T lymphocytes` + `NK cells`. This folder
does not re-rank or re-pool that Spearman.

## What is new

Outgoing communication **CLDN4-high malignant → same-unit T/NK**, then
a meta-analysis of **patient Δ** across the two units.

High-end split (primary): within each unit, malignant cells are cut on
`log1p(CP10k)` **CLDN4** into Q4 vs Q1 (middle half dropped). Sensitivity:
median and tertile. TACSTD2 is not used.

Δ = score_high − score_low. One delta per unit per ligand–receptor pair.
Unit means are DerSimonian–Laird random-effects pooled. A patient-level
Wilcoxon on the stacked deltas is also reported. p-values are descriptive.

Extra between-unit cut: quartiles of the **given** locked CLDN4 %pos
scores (same vectors as the n=34 row). Per-unit score uses all malignant
cells → that unit’s T/NK. This is not a Spearman re-audit.

## CellChat-style probability

Jin et al. 2021 Hill / mass-action on CellChatDB v2 protein pairs:

1. Keep a pair only if every ligand and receptor subunit is in the matrix.
2. Per cell group: 10% truncated mean of `log1p(CP10k)`; complexes =
   geometric mean (0 if any subunit mean is 0).
3. Detected if each complex has expressing-cell fraction ≥ 0.10.
4. \(P = (L \cdot R) / (0.5 + L \cdot R)\).
5. Arms with fewer than 10 cells are not scored. Within-unit Q4 vs Q1
   also requires ≥40 malignant cells.

Outgoing only is the ligand table (Mal → T/NK). CellChat R was not run.

## LIANA-style score

CellPhoneDB-style **mean of partner means** on the same CellChat
complexes (Efremova 2020 / Garcia-Alonso 2022): \(S = 0.5 \times (L + R)\).
Same detect gate and patient-Δ meta as CellChat. LIANA
`mt.cellphonedb` permutations were not run.

## Honest n and skipped files

Given combo n=34 is not the LR n. Units that fail the malignant or T/NK
cell floor are dropped from the LR table and are reported in
`results/patient_coverage.tsv`.

GSE148071 is not downloaded. The other five PR #459 cohorts are not
piled. No FASTQ / EGA / SRA. No file ≥ 2 GB. The 36.5 GB Laughney H5
and the 2.9 GB GSE131907 log2TPM text are skipped.

## Software

Python (numpy / pandas / scipy / matplotlib). Public GEO processed
files only: `GSE123902_RAW.tar` (90 MB dense CSV) and
`GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz` (390 MB) + author
annotation (1.8 MB).
