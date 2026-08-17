# Methods — merged GSE131907 + GSE205335 CLDN4-only LIANA / CellPhoneDB

ADDITIVE. **CLDN4 only.** No dual-high TACSTD2×CLDN4 gate. No GSE207422.

## Slice (given; not re-audited)

PR #320 author-malignant %pos vs T/NK on **GSE131907 + GSE205335**:
Spearman N=43 ρ=−0.479; **Q4 vs Q1 n=23 (12/11) r=−0.705**. That
between-patient quartile cut is taken as given and is not re-run.

Locked extracts reused here:

| Cohort | Unit in PR #320 | Gate | n |
|---|---|---|---:|
| GSE131907 | sample | tumor-origin ∩ author malignant ≥20 | 21 |
| GSE205335 | patient | author malignant ≥20 and T/NK ≥20 | 22 |

This folder asks a **different** question on the same slice: outgoing
ligand–receptor scores from **CLDN4-high vs CLDN4-low malignant cells
to same-patient T/NK**, tested with a **patient-level paired Wilcoxon**.

## Labels (author, not recomputed)

**GSE131907** (Kim et al., *Nat Commun* 2020, PMID 32385277). 10x.
Treatment-naive LUAD. No ICI / MPR labels.

- Tumor-origin: `Sample_Origin` ∈ {tLung, tL/B, mLN, mBrain, PE}
- Malignant: tumor-origin **and** `Cell_subtype` ∈ {Malignant cells, tS1, tS2, tS3}
- T/NK: `Cell_type` ∈ {T lymphocytes, NK cells}
- PE unlabeled epithelium and nLung AT1/AT2/Club/Ciliated are not malignant.

**GSE205335** (Hu et al. Lung IO). Author `lineage.sub == Malignant cells`
and `lineage.total == T/NK cells`. MPR/NMPR is unlabeled.

## Patient unit

- GSE205335: GEO `patient` (already the PR #320 unit).
- GSE131907: GEO `patient_id`. Cells from locked tumor-origin samples
  belonging to the same patient are **pooled**. Unique-patient n is the
  Wilcoxon n. Sample n (21) is reported separately because the PR #320
  T/NK extract is sample-level.

A patient enters the paired test if it has ≥10 CLDN4-high malignant,
≥10 CLDN4-low malignant, and ≥20 T/NK cells. Cells are not n.

## CLDN4 split

Among malignant cells **within each cohort**, CLDN4 is `log1p(CP10k)`
using the cell’s full UMI total. **High** = at or above that cohort’s
global malignant median; **low** = below. Thresholds are not shared
across platforms. TACSTD2 is stored and is not a gate.

## Score (primary)

CellPhoneDB v5 interactors (`resources/cellphonedb_v5_lr_pairs.tsv`) plus
a curated T-recruit / IFN / MHC-I overlay. For sender S and receiver R,
each partner is the **minimum of subunit means** on log1p(CP10k). The
pair score is the **arithmetic mean of the two partner means**
(Efremova et al., *Nat Protoc* 2020; Garcia-Alonso et al., *Nat Protoc*
2022). `pass_expr_prop` requires both partners in ≥10% of cells in their
group.

Outgoing = CLDN4-high (or low) malignant → **same-patient** T/NK.
Median Δ = high − low. Negative = weaker from CLDN4-high.

Paired Wilcoxon (two-sided) on patients who pass the detect gate for
that pair (`pass_either` for focus axes; `pass_both` for other pairs,
≥3 patients). FDR is Benjamini–Hochberg within the outgoing contrast.

This is **not** a CellChat communication probability.

## LIANA (secondary)

`liana.mt.cellphonedb` if the package imports, `expr_prop=0.10`,
`n_perms=50`, ≤2,000 cells/group, run **separately per cohort**.
LIANA p-values are within-object specificity, not patient tests.
If import fails, the documented CellPhoneDB mean is the result.

## CellChat

Not run. R / CellChat is not required for the done criterion.
