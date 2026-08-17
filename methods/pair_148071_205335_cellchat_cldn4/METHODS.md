# Methods — pairwise GSE148071 + GSE205335, CLDN4-only CellChat

ADDITIVE. **CLDN4 only.** No dual-high TACSTD2×CLDN4 score. **GSE131907 is
not merged.** Patient is the unit. Prior single-cohort folders (GSE148071
CellChat, GSE205335 Q4 extra) and the multi-cohort Q4 meta are taken as
given.

## Pair

| Cohort | Malignant | T/NK | Eligible n | Source |
| --- | --- | --- | ---: | --- |
| GSE148071 | marker-argmax epithelium (putative) | marker T+NK | 25 / 42 | PR #348 per-sample table |
| GSE205335 | author `Malignant cells` | author `T/NK cells` | 22 | PR #279 / #320 / #362 |

Eligibility: GSE148071 ≥25 epithelial **and** ≥25 T/NK; GSE205335 ≥20
author-malignant **and** ≥20 T/NK (locked extract; four GEO patients with
near-zero malignant cells are already out).

GSE148071 has no histology / ICI / MPR labels on GEO. GSE205335 RECIST is
not substituted for MPR.

## Patient-level CLDN4 vs T/NK

T/NK fraction = `n_TNK / n_cells` (all captured cells). CLDN4 is scored
on the malignant / putative-malignant compartment only (`%pos` primary;
mean `log1p(CP10k)` secondary).

Quartiles are cut **within cohort** with `pd.qcut` on average ranks (same
rule as `methods/cldn4_malig_q4_tnk`). Q4 vs Q1 is two-sided Mann–Whitney
on T/NK fraction; rank-biserial \(r = 2U/(n_4 n_1) - 1\). \(r < 0\) means
Q4 has lower T/NK. n_compared = n_Q1 + n_Q4, not the full eligible n.

**Combo ρ** is DerSimonian–Laird random-effects on Fisher-z of the two
cohort Spearmans (not a stack of raw %pos). A stacked Spearman on
within-cohort ranks and a stacked Q4 vs Q1 (within-cohort labels) are
reported as companions. p-values are descriptive.

## CellChat-style probability

CellChat R and LIANA are not run. Probability follows Jin et al. 2021 on
**CellChatDB v2** protein pairs:

1. Keep a pair only if every ligand and receptor subunit is in that
   cohort's gene universe.
2. Per cell group: 10% truncated mean of `log1p(CP10k)`; complexes =
   geometric mean of subunits (0 if any subunit mean is 0).
3. Detected if each complex has expressing-cell fraction ≥ 0.10 (AND).
4. \(P = (L \cdot R) / (K_h + L \cdot R)\) with \(K_h = 0.5\).

**Same-patient pairing.** For each eligible patient, outgoing = Mal →
that patient's T/NK; incoming = T/NK → Mal. Q4 vs Q1 uses the
within-cohort CLDN4 %pos tails, then stacks the two cohorts. A pair is
listed if it is detected in ≥3 Q4 and ≥3 Q1 patients. The ligand table
is **outgoing** Mal → T/NK.

GSE205335 per-patient *P* is reused from the locked PR #362 table. GSE148071
is re-scored from the public `GSE148071_RAW.tar` matrices with the same
lineage rule as PR #348.

Cell-pooled (stacked cells) truncated means are not the test.

## Honest n

- Unit = patient.
- Eligible n is 25 + 22, not 42 + 26.
- Q4 vs Q1 tails are thinner than 47.
- One biopsy can dominate a cell-pooled mean; that is why the pair test
  is per-patient.
- Malignancy defs differ and are stated on every row.

## Software

Python (numpy / pandas / scipy / matplotlib). CellChatDB v2 from
`jinworks/CellChat`. Raw FASTQ and EGA were not used.
