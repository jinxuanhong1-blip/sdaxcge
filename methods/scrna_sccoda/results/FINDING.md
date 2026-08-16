# Combinatorial finding (public only)

**Question.** Which **cohort × annotation** recovers **T/NK or TLS down** in
patients with high malignant TACSTD2, using a Dirichlet-multinomial / ALR
composition model and honest FDR / n?

**Answer.** **None**, after BH on the pre-specified primary family
(52 patient-level permutation tests). Several cells of the grid point in the
hypothesized direction; none survive q < 0.10.

## Primary contract

- Exposure: within-cohort median split of malignant / author-Epi TACSTD2 mean log1p.
- Test: ALR slope, **patient permutation p** (exact when `C(n, n_high) ≤ 2000`).
- Family: TNK / TLS (or T, NK, B, plasma when not collapsed) × cohort × annotation.
- Recovery: effect < 0 and BH q < 0.10. Bonferroni on the same 52 tests is 1 for every row.
- scCODA HMC: `not_available` (no anndata/tensorflow). Not faked.
- DM-GLM LRT: diagnostic only. It scales with n_cells and is anti-conservative.

## Honest n

| Cohort | Patients | TACSTD2-eligible | MPR-eligible | Notes |
|---|---:|---:|---:|---|
| GSE207422 | 15 | 10 | 12 post | Hu markers recomputed; DRMref n=7–12 |
| GSE241934 IIT | 11 | 11 | 11 | EGFR-mut NEOTIDE; author labels |
| GSE241934 RWC | 34 | 29 | 34 | EGFR-WT real-world; author labels |
| GSE241934 pooled | 45 | 40 | 45 | same patients as IIT+RWC (not independent) |
| GSE253013 | 9 tumor | 9 | **0** | no public MPR; 9.3 GB RDS not downloaded |

Unit = patient. Cells are library size only.

## Closest unadjusted cells (not recoveries)

| Cohort × annotation | Compartment | ALR effect | perm p | BH q | n |
|---|---|---:|---:|---:|---:|
| GSE207422 × `hu_markers` | B (TLS) | −2.56 | 0.032 | 0.52 | 10 |
| GSE207422 × `hu_markers` | T (TNK) | −1.59 | 0.040 | 0.52 | 10 |
| GSE253013 × `marker_coarse` | TNK | −1.07 | 0.063 | 0.52 | 9 |
| GSE241934_RWC × `author_collapsed` | TNK | −0.98 | 0.068 | 0.52 | 29 |

GSE241934 author-major / pooled n=40 is near null (TNK p≈0.35). IIT n=11 is
underpowered and often opposite in sign for TNK.

## MPR (separate family)

No primary TACSTD2-style recovery is claimed for MPR. Companion ALR / MWU
tests of TNK/TLS ~ MPR are mostly null (smallest fraction-MWU p≈0.048 for
GSE241934 author-major NK, unadjusted, n=45). Hu et al. already reported
higher T/B in MPR without a significant p at n=4 vs 8.

## What this does *not* say

- It does not refute a biological TACSTD2–exclusion hypothesis; it says the
  public scRNA grid does not recover it after honest FDR / n.
- It does not use author CopyKAT IDs (not public on GSE207422).
- GSE253013 cannot test TLS (bundled table has no B counts) or MPR.
- Pooled GSE241934 is not an independent replication of IIT or RWC.
