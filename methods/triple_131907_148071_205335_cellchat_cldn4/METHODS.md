# Methods — triple GSE131907 + GSE148071 + GSE205335, CLDN4-only

ADDITIVE. **CLDN4 only.** No dual-high TACSTD2×CLDN4. No GSE207422. The
131907+205335-only CellChat (separate agent) is not re-run.

## Unit and floors

Patient is the unit. GSE131907 T/NK extract is **sample-level** (PR #320);
that is stated on those rows. A unit is scored for the patient table if it
has malignant (or putative malignant) cells and same-unit T/NK cells.

| cohort | malignant | T/NK | floor | locked n |
|---|---|---|---|---:|
| GSE131907 | author `Malignant cells` (not tS1–tS3) | author T + NK | PR #320 extract | 21 samples |
| GSE148071 | putative epithelium (marker-argmax; no GEO labels) | T + NK merged | ≥25 epi and ≥25 T/NK | 25 patients |
| GSE205335 | author `lineage.sub` = Malignant cells | author `lineage.total` = T/NK cells | ≥20 / ≥20 | 22 patients |

GSE148071 n=25 partial ρ=−0.49 and PR #320 131907+205335 Q4 *r*=−0.705 are
**given** and are not re-audited. The new row is the **triple** Fisher-z pool.

## Combo rho

Within-cohort Spearman of malignant CLDN4 vs same-unit T/NK fraction.
Primary cut is CLDN4 **% positive** (PR #320 author/tnk/pct family) plus the
given GSE148071 partial. Pool = DerSimonian–Laird on Fisher-z. Q4 vs Q1 is
Mann–Whitney rank-biserial *r* on within-cohort quartile tails (n≥8 and both
tails ≥3 to pool). p-values are descriptive.

## CellChat-style ligands

CellChat R was not required. Communication probability follows Jin et al.
2021 on CellChatDB v2 protein pairs:

1. Keep a pair only if every ligand and receptor subunit is in that cohort’s gene list.
2. Per patient, per compartment: 10% truncated mean of `log1p(CP10k)`; complexes = geometric mean (0 if any subunit mean is 0).
3. Detected if each complex has expressing-cell fraction ≥ 0.10 (AND).
4. \(P = (L \cdot R) / (K_h + L \cdot R)\) with \(K_h = 0.5\).
5. Outgoing = malignant → **same-patient** T/NK. Incoming is recorded but is not the ligand table.

CLDN4-high vs low uses **within-cohort** quartiles of malignant CLDN4 %pos
(do not mix cohort scales). Test = Mann–Whitney on per-patient *P* among
Q4 vs Q1 patients with the pair detected (arm n≥3). Combined Q4 vs Q1
stacks the three cohorts’ tails. Cell-pooled truncated means are not the test.

## Data

Public GEO processed UMI only (GSE131907 raw UMI text; GSE148071 per-biopsy
count matrices; GSE205335 UMI `dgCMatrix` RDS). The 3.07 GB GSE131907
log2TPM matrix is skipped. TISCH metadata are used only as a lineage check
for GSE148071. No EGA FASTQ. No dual-high score.
