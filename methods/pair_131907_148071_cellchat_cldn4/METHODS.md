# Methods — pair GSE131907 + GSE148071, CLDN4-only vs T/NK + CellChat

ADDITIVE. **CLDN4 only.** TACSTD2 is recorded but does not define groups.
**GSE205335 is not merged.** Patient is the unit. Prior single-dataset
CellChat (GSE131907 PR #347, GSE148071 PR #348) and the Q4 meta that
pairs GSE131907 with GSE205335 (PR #320) are given.

## Datasets

| GEO | Citation | Malignant | T/NK | Floor |
|---|---|---|---|---|
| [GSE131907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907) | Kim et al., *Nat Commun* 2020 | author `Epithelial cells` ∩ `{Malignant cells, tS1, tS2, tS3}` | author T lymphocytes + NK | ≥20 malignant **and** ≥20 T/NK; tumor origins tLung / tL/B / mLN / mBrain |
| [GSE148071](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148071) | Wu et al., *Nat Commun* 2021 | marker-argmax epithelium (putative) | marker-argmax T+NK | ≥25 epithelial **and** ≥25 T/NK |

GSE131907 PE / nLung / nLN are out (PE epithelium is unlabeled; normals have no author-malignant cells). After the floor, each kept GSE131907 sample is one patient. GSE148071 has no author labels or CopyKAT on GEO.

No ICI / MPR labels on either series. EGA raw was not accessed.

## Patient-level CLDN4 vs T/NK

Within each cohort, rank malignant CLDN4 (`pct_pos` primary; `mean` log1p(CP10k) secondary) and cut quartiles with `pd.qcut` on average ranks (same rule as `methods/cldn4_malig_q4_tnk`). Compare same-patient T/NK fraction in Q4 vs Q1 by two-sided Mann–Whitney. Rank-biserial \(r = 2U/(n_4 n_1) - 1\); \(r < 0\) means Q4 has lower T/NK. Spearman ρ on the full eligible n is reported beside the tails. n_compared = n_Q1 + n_Q4, not the full n.

**Combo ρ** is DerSimonian–Laird random-effects on Fisher-z of the two within-cohort Spearmans, back-transformed to ρ. Combo Q4 vs Q1 is the same pool on the rank-biserial *r* values (n = n_compared). A stacked 10x+Singleron Spearman is written as a companion only.

## CellChat-style ligands

CellChat R and LIANA are not run. Probability follows Jin et al. 2021 on **CellChatDB v2** protein pairs:

1. Keep a pair only if every ligand and receptor subunit is present for that patient.
2. Per cell group: 10% truncated mean of `log1p(CP10k)`; complexes = geometric mean of subunits (0 if any subunit mean is 0).
3. Detected if each complex has expressing-cell fraction ≥ 0.10 (AND rule).
4. \(P = (L \cdot R) / (K_h + L \cdot R)\) with \(K_h = 0.5\).
5. A patient is scored only if it is in a within-cohort Q4 or Q1 tail **and** has ≥20 malignant and ≥20 T/NK cells.

Outgoing = Mal → that patient's T/NK. Incoming = T/NK → Mal. Q4 vs Q1 is Mann–Whitney on the per-patient probabilities. A pair is listed if it is detected in ≥3 Q4 and ≥3 Q1 patients and \(P_{\mathrm{Q4}} \ne P_{\mathrm{Q1}}\) on the median.

Cell-pooled (stacked) truncated means are not the test.

## Honest n

- Unit = patient.
- GSE131907 malignant is the authors' label, not a new CNV call. tS1–tS3 are the tLung tumor-epithelial labels.
- GSE148071 epithelium is marker-argmax (putative malignant).
- Quartiles are within-cohort so the two chemistries are not ranked on one scale.
- Thin tails stay flagged. p-values are descriptive.

## Software

Python (numpy / pandas / scipy / matplotlib). CellChatDB v2 from `jinworks/CellChat`. GEO matrices are not stored in git.
