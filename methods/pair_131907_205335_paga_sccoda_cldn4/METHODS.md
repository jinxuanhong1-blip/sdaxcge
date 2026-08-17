# Methods — pair PAGA + scCODA, CLDN4 only

## Object

Public processed UMI only.

| Cohort | Malignant | T/NK | Unit |
|---|---|---|---|
| GSE131907 Kim 2020 | author `Cell_subtype == Malignant cells` | `Cell_type ∈ {T lymphocytes, NK cells}` | **sample** with n_malignant ≥ 20 on tLung / tL/B / mLN / PE / mBrain |
| GSE205335 Hu 2023 | author `lineage.sub == Malignant cells` | `lineage.total == T/NK cells` | **patient**; drop `Normal*` libraries |

GSE131907 tLung author-malignant is empty in the given extract; the 21
winning samples are metastases. That is stated, not hidden.

Skipped: 2.86 GB GSE131907 log2TPM, EGA FASTQ, GSE148071, GSE207422,
7-cohort pool.

## CLDN4 exposure

Author-malignant **% positive**, the same vector as PR #320.

- Primary: within-cohort quartiles; Q4 vs Q1 tails (merged n=23).
- Sensitivity: within-cohort median split.
- No TACSTD2 gate. No dual-high score.

## scCODA

scCODA HMC is not available (no TensorFlow). Not faked.

Primary test = ALR of T/NK or B vs Other, **unit-level permutation p**
(exact enumeration when `C(n, n_high) ≤ 3000`). BH on the 6 Q4 vs Q1
tests. Fraction MWU and continuous Spearman are companions. Dirichlet-
multinomial two-group LRT is diagnostic only (library-size inflated).

GSE205335 B = B+plasma, matching PR #320 / #339.

## PAGA

A prior run Harmony/UMAP'd **all 53,296** malignant cells and died on
~15 GB RAM. This retry caps **150 malignant and 80 T/NK cells per unit**
before matrix subset (stated in FINDING).

- HVG 2000, 20 PCs, 20 neighbors, Leiden 0.6
- Harmony on `dataset` when it imports; PCA fallback otherwise
- DPT root = lowest-CLDN4 cell in the top AT2-score tercile of
  **malignant** cells. Never a CLDN4-high cell. Not comparable to
  PR #325 nLung author-AT2 root (those cells are not in this object).
- Barrier/keratin excludes CLDN4.
- Inferential tests = Spearman on unit means (n = samples/patients
  with ≥8 cells in the subsample). Cell-level ρ is not reported as n.

T/NK PAGA is an extra embedding (no AT2-rooted DPT).
