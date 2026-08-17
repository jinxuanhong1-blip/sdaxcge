# Methods — QUAD merge CLDN4-only LIANA / CellPhoneDB

ADDITIVE. **CLDN4 only.** TACSTD2 is not a gate and does not define dual-high.
This is **not** a re-audit of the GSE207422-only LIANA (PR #344). GSE207422 is
included in the merge because the user asked for it.

## Question

Do CLDN4-high malignant cells show weaker outgoing T-recruit / MHC-I
communication toward same-patient T/NK than CLDN4-low cells, when the four
public lung scRNA objects are scored with one documented CellPhoneDB-style
rule and one patient-level test?

## Cohorts (public processed files only)

| Accession | Paper | Platform / object | Malignant | T/NK |
|---|---|---|---|---|
| GSE207422 | Hu et al., *Genome Med* 2023 | BD Rhapsody UMI | A3: epithelial AND zero UMI of SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3 | A3 T or NK |
| GSE131907 | Kim et al., *Nat Commun* 2020 | 10x raw UMI | tumor-origin AND `Cell_subtype` ∈ {Malignant cells, tS1, tS2, tS3} | T lymphocytes + NK cells |
| GSE148071 | Wu et al., *Nat Commun* 2021 | TISCH2 `log2(TPM/10+1)` | TISCH2 major-lineage `Malignant` | TISCH2 T/NK labels (NK is often 0) |
| GSE205335 | Hu / Ahn / Lee, GEO 2024 | 10x UMI `dgCMatrix` | tumor sample AND `lineage.sub == Malignant cells` | `lineage.total == T/NK cells` |

Controlled raw (GSA-Human HRA001033; EGA EGAD00001005054 / EGAD00001008703)
is not accessed. Labels are taken as given from each object’s public
annotation or the A3 rule. They are not re-called with inferCNV/CopyKAT.

## What “merge” means

Expression is **not** concatenated across platforms. BD Rhapsody UMI, 10x UMI,
and TISCH2 TPM are different scales. Harmony/scVI latents are not used for
ligand–receptor scores.

Each cohort is scored on its documented scale:

- GSE207422 / GSE131907 / GSE205335: `log1p(CP10k)` from the cell’s full UMI total.
- GSE148071: TISCH2 `log2(TPM/10+1)` as stored.

CLDN4-high / low is the **within-cohort** global median among that cohort’s
malignant cells. Patients are then stacked for a paired Wilcoxon. That is the
merge. Cells are not n.

## Patient unit

A **patient** enters the paired test if, after pooling that patient’s tumor
samples in that cohort, it has ≥10 malignant cells in both CLDN4 bins and
≥20 T/NK cells.

- GSE207422: one BD sample = one patient (P01–P15).
- GSE131907: GEO `patient_id`; multiple tumor sites are one unit. Unique
  sample n is reported separately.
- GSE148071: TISCH2 `Patient`.
- GSE205335: GEO `patient`; normal-only samples are excluded from that
  patient’s malignant and T/NK pool.

GSE131907 is treatment-naive (no ICI/MPR). GSE205335 has RECIST, not MPR.
Response is recorded and is **not** the test.

## Pair list

CellPhoneDB v5 (`ventolab/cellphonedb-data`, `interaction_input.csv`).
The `interactors` string is split only when both sides match the gene list,
so hyphenated symbols such as HLA-A are preserved. Complexes use `+`
subunits. A small curated overlay adds canonical T-recruit / IFN / MHC-I
pairs if a parse miss would drop them.

Pathway labels organize the ranks. They do not subset the score.

## Documented score (primary)

For sender group S and receiver group R, each partner’s expression is the
**minimum of subunit means** (CellPhoneDB complex rule). The pair score is
the **arithmetic mean of the two partner means** (Efremova et al., *Nat
Protoc* 2020; Garcia-Alonso et al., *Nat Protoc* 2022). A pair **passes
`expr_prop`** when both partners are detected in ≥10% of cells in their group.

The test is a **paired Wilcoxon** of the high vs low score across patients
in the merge. FDR is Benjamini–Hochberg within each contrast (outgoing or
incoming). Per-cohort ranks are written as extra tables.

This is **not** a CellChat communication probability.

## LIANA (secondary, if installed)

`liana.mt.cellphonedb` on a downsampled **within-cohort** object
(`Malig_CLDN4high` / `Malig_CLDN4low` / `T` / `NK`; ≤1,500 cells/group;
20 permutations). Joint LIANA on concatenated expression is not run.
LIANA p-values are within-object specificity, not patient-level tests.

## CellChat

Not run. This environment has no R.

## Software

Python 3; numpy, pandas, scipy, statsmodels, matplotlib, seaborn, pyyaml,
h5py, rdata. Optional scanpy + liana. See `requirements.txt`.
Matrices are not stored in git.
