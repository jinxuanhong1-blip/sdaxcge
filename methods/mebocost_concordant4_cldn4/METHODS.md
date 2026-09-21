# Methods — MEBOCOST concordant-four CLDN4-high vs low metabolite senders

ADDITIVE. This is a metabolite communication layer on the same concordant-four
units as the CellChat protein ligand–receptor run. It does not replace that
result. CLDN4 only. No dual-high. Not a full-pool.
Patient / locked sample is the unit.

Concordant four: GSE123902 + GSE131907 + GSE205335 + GSE189357.
Do not add GSE148071 / GSE127465 / CD45-only.

## Tool choice

Requested order was CINE or MEBOCOST (metabolite-mediated CCC), else SpatialDM
on a spatial subset, else MultiNicheNet if installation failed.

CINE was not installed. Web and GitHub searches did not return a
metabolite-mediated CCC package named CINE. The name collides with unrelated
tools (CineMA cardiac MRI; comet infrared `cine`). Those were not used.

MEBOCOST (kaifuchenlab, cloned from GitHub, human sensor table
`human_met_sensor_update_Oct21_2025.tsv`) imported on Python 3.12:

```python
from mebocost import mebocost
```

SpatialDM and MultiNicheNet were not run. Concordant-four is scRNA-seq, and
the metabolite tool installed.

## Engine

`mebocost.create_obj(..., species='human', sensor_type='All', thread=1)`
then `infer_commu`. Metabolite levels use the package default
`met_est='mebocost'`: arithmetic mean of product-direction enzyme expression
minus substrate-direction enzyme expression. COMPASS flux and scFEA were not
applied. The communication score is the package `Commu_Score` (product of
sender metabolite level and receiver sensor expression). This script does not
reimplement that product.

Settings that differ from the package defaults, and why:

- `n_shuffle=10`, `seed=12345`. The package default is 1000. The estimand is
  the within-unit `Commu_Score` difference, which is computed before the
  shuffle. Shuffle count changes the within-sample permutation p-value and
  `Norm_Commu_Score`, not `Commu_Score`. Ten draws are stored so the official
  `infer_commu` path runs; they are not the primary test.
- `min_cell_number=10`. The package default is 50. Quartile arms on these
  units are often below 50. The cell floor matches the CellChat arm floor
  (≥10 per arm, T/NK ≥20, Q4 vs Q1 also n_mal ≥40).
- `pval_cutoff=1.01` so the call returns every pair. The table used for Δ is
  `original_result` (scores before low-proportion p-values are set to 1).
- `cutoff_exp='auto'`, `cutoff_met='auto'`, `cutoff_prop=0.15` stay at the
  package defaults and affect presence proportions, not the decision to read
  `Commu_Score`.

Input expression is log1p of counts per 10,000 using each cell's
full-transcriptome UMI sum (GSE131907 streams that sum; GSE205335 uses
`Matrix::colSums` of the RDS before gene subsetting; CSV/MTX use the loaded
matrix sum). The object passed to MEBOCOST contains only enzyme, sensor, and
lineage-marker genes. That is intentional: metabolite estimation only reads
enzyme genes, and the sensor table only reads sensor genes.

## Matrices

Public processed GEO only. Same files as the CellChat concordant-four run.

- GSE123902 `GSE123902_RAW.tar` (Laughney dense UMI CSVs). Donor-level.
  PRIMARY is kept when a donor also has METASTASIS. NORMAL is dropped.
- GSE131907 raw UMI TXT + author cell annotation. The log2TPM text and RDS
  files are not used.
- GSE205335 UMI RDS is double-gzipped. R 4.3.3 + Matrix reads the dgCMatrix
  and writes a gene-subset Matrix Market plus library size. Python does not
  reimplement MEBOCOST; R is only the RDS reader. GEO does not ship MTX for
  this series.
- GSE189357 `GSE189357_RAW.tar` (10x MTX, TD1–TD9).

## Malignant / T/NK

- GSE123902 and GSE189357: marker-malignant = EPCAM|KRT8|KRT18|KRT19 count > 0
  and PTPRC == 0. T/NK = CD3D|CD3E|CD8A|NKG7|GNLY|KLRD1 count > 0 and not
  malignant.
- GSE131907: author `Cell_subtype` in {Malignant cells, tS1, tS2, tS3};
  T/NK = `Cell_type` in {T lymphocytes, NK cells}. Samples with
  n_malignant > 0 in the locked sample table.
- GSE205335: author `lineage.sub` == Malignant cells; T/NK =
  `lineage.total` == T/NK cells. Rows whose gsm-map tissue starts with
  "Normal " are excluded. Locked patients with n_malignant > 0.

TACSTD2 is never a gate.

## Gates and splits

CLDN4-only. Within-unit malignant cells are ranked by log1p(CLDN4 CP10k),
`rankdata(..., method='ordinal')`, which matches R `rank(..., ties.method='first')`.

- Primary: Q4 vs Q1. Requires n_mal ≥ 40. High = rank > ceil(0.75 n).
  Low = rank ≤ floor(0.25 n).
- Extra: median split (high > median, low ≤ median) and %pos (CLDN4 > 0 vs = 0).

Honest n requires T/NK ≥ 20 and each arm ≥ 10. Middle malignant cells are
not a sender. Per unit, MEBOCOST sees three groups: `CLDN4_high`,
`CLDN4_low`, `TNK`.

## Pre-specified family

Not a discovery screen. Pairs are rows in the installed human sensor table.
Expectation was fixed as high > low before looking at Δ: immunosuppressive
metabolite signal from CLDN4-high malignant cells toward T/NK.

| pair | metabolite (DB name) | sensor | HMDB |
|---|---|---|---|
| Adenosine_ADORA2A | Adenosine | ADORA2A | HMDB0000050 |
| Adenosine_ADORA2B | Adenosine | ADORA2B | HMDB0000050 |
| PGE2_PTGER2 | Prostaglandin E2 | PTGER2 | HMDB0001220 |
| PGE2_PTGER4 | Prostaglandin E2 | PTGER4 | HMDB0001220 |
| DLactate_HCAR1 | D-Lactic acid | HCAR1 | HMDB0001311 |
| LLactate_SLC16A1 | L-Lactic acid | SLC16A1 | HMDB0000190 |

HCAR1 (GPR81) is curated on D-lactic acid in this database. The product-direction
enzymes for that HMDB id are HAGH and HAGHL, not LDHA. L-lactic acid
(HMDB0000190) has only substrate-direction reactions in
`metabolite_associated_gene_reaction_HMDB_summary.tsv`. MEBOCOST emits a
metabolite only when a product-direction enzyme is present, so
L-lactate–SLC16A1 is pre-specified and then not scored. That is a database
limit. Kynurenine–AHR is also absent (closest row is kynurenic acid–GPR35)
and is not scored.

PGE2 product-direction genes in this database include PTGES, PTGES2, PTGES3,
CBR1, and CBR3. Adenosine product-direction genes include NT5E and the other
5'-nucleotidases.

A pair is detected in a unit when both arm scores are finite and at least
one is > 0. Δ = score_high − score_low. The family value in a unit is the
mean of detected pair Δs. The primary test is a two-sided Wilcoxon
signed-rank of those unit values against 0 (`scipy.stats.wilcoxon`,
`zero_method='wilcox'`, zeros discarded, `method='auto'`).

`descriptive_all_pairs_q4q1.tsv` aggregates every MEBOCOST metabolite–sensor
with a detected high→TNK or low→TNK score. It is not a claim.

## Not done

Dual-high, GSE148071, GSE127465, CD45-only, 7-pool, cell-pooled tests,
COMPASS, scFEA, SpatialDM, MultiNicheNet, kynurenine–AHR, discovery screen.
