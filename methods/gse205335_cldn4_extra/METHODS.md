# Methods — GSE205335 malignant CLDN4 Q4 vs same-patient T/NK + CellChat-style ligands

ADDITIVE. **CLDN4 only.** TACSTD2 is recorded but does not define groups.
Patient is the unit. Prior TACSTD2 A3 and the multi-cohort Q4 meta are given.

## Dataset

GEO [GSE205335](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205335)
(Hu et al., palliative ICI biopsy/effusion scRNA). Open processed UMI
`dgCMatrix` (33,714 genes × 96,505 cells), author cell-identity table, and
GEO SOFT. Controlled EGA raw (`EGAD00001008703`) was not accessed.

MPR/NMPR is not labeled. RECIST is not substituted for MPR.

## Given patient extract

`data/GSE205335_patients.tsv` is the locked 22-patient table from the
malignant T/NK extract (PR #279 / PR #320). Author `lineage.sub ==
"Malignant cells"` and `lineage.total == "T/NK cells"`. Eligibility was
≥20 cells in each compartment. Four patients with (near-)zero captured
malignant cells are already dropped (asymmetric: 3 PR + 1 PD).

## Q4 vs Q1 vs same-patient T/NK

Within the 22 patients, rank malignant CLDN4 (`pct_pos` primary; `mean`
log1p(CP10k) secondary) and cut quartiles with `pd.qcut` on average ranks
(same rule as `methods/cldn4_malig_q4_tnk`). Compare same-patient T/NK
fraction (`n_tnk / n_cells`) in Q4 vs Q1 by two-sided Mann–Whitney.
Rank-biserial \(r = 2U/(n_4 n_1) - 1\); \(r < 0\) means Q4 has lower T/NK.
Spearman ρ on the full n=22 is reported beside the quartile tails.
p-values are descriptive. n_compared = n_Q1 + n_Q4 (6+6), not 22.

## CellChat-style ligands

CellChat R and LIANA are not run. Probability follows Jin et al. 2021 on
**CellChatDB v2** protein pairs (same parse as `methods/scrna_cellchat_cldn4`):

1. Keep a pair only if every ligand and receptor subunit is in the UMI matrix.
2. Per cell group: 10% truncated mean of `log1p(CP10k)`; complexes =
   geometric mean of subunits (0 if any subunit mean is 0).
3. Detected if each complex has expressing-cell fraction ≥ 0.10 (AND rule).
4. \(P = (L \cdot R) / (K_h + L \cdot R)\) with \(K_h = 0.5\).
5. A patient is scored only if it has ≥20 malignant and ≥20 T/NK cells.

**Same-patient pairing.** Q4 vs Q1 is the **patient** split on malignant
CLDN4 %pos (primary). For each of those 12 patients, outgoing = Mal → that
patient's T/NK; incoming = that patient's T/NK → Mal. Q4 vs Q1 is then
Mann–Whitney on the per-patient probabilities (honest n = 6 vs 6). A pair
is listed if it is detected in ≥3 Q4 and ≥3 Q1 patients and
\(P_{\mathrm{Q4}} \ne P_{\mathrm{Q1}}\) on the median.

A pooled (cell-stacked) Q4 vs Q1 score is written as a companion column
only. Cell-pooled truncated means are not a patient-level mixed model.

## Honest n

- Unit = patient.
- Q4 vs Q1 tails are 6 vs 6. That is thin; it is reported, not pooled away.
- One effusion/LN/liver mix and mixed 3'/5' chemistry can dominate a
  cell-pooled mean. Per-patient cell counts are in `results/per_patient_cells.tsv`.
- Malignancy is the authors' label, not an independent CNV re-call.

## Software

Python (numpy / pandas / scipy). CellChatDB v2 from `jinworks/CellChat`.
The 500 MB GEO matrix is not stored in git.
