# Methods — max-effect CellPhoneDB / Connectome barrier family

ADDITIVE. CLDN4 only. No dual-high. Not a spatial test.
Concordant four: GSE123902 + GSE131907 + GSE205335 + GSE189357.
Do not add GSE148071, GSE127465, GSE154826, GSE200563, or E-MTAB-13526.
The unit is the locked patient or sample.

## Question

Which pre-specified CLDN4 contrast maximizes the CellPhoneDB magnitude Δ of the
barrier ligands F11R, NECTIN2, CDH1, and LGALS9 from malignant cells to T/NK,
while both CellPhoneDB and Connectome stay positive in every concordant-4 cohort?

## Scores

liana 1.10.0 formulas, checked in `scripts/run_max_effect.py` against
`liana.method.cellphonedb` before any cohort is scored.

Expression is log1p counts-per-10k using the full-matrix library size.
For a complex, the side mean is the minimum subunit mean, and the side is
detected only when every subunit has nonzero proportion at least `expr_prop`
(liana `complex_policy='min'`). CellPhoneDB `lr_means` is the average of the
ligand and receptor means, and is 0 if either side is 0. Connectome `expr_prod`
is the product of those two means.

An edge enters a patient's ligand score only when the receptor complex is
detected on that patient's receiver. The ligand score is the mean of its
detected edges. The family score is the unweighted mean of the ligands that
have a detected edge. A patient is in the family test only when at least 3 of
the 4 ligands are detected. Undetected ligands are omitted, not filled with 0.

## Grid (locked)

Gates, within malignant cells of one patient:

| name | high | low |
|---|---|---|
| median | top 50% by ordinal CLDN4 rank | bottom 50% |
| tertile | top 1/3 | bottom 1/3 |
| q4q1 | top 25% (`rank > ceil(0.75 n)`) | bottom 25% (`rank <= floor(0.25 n)`) |
| quintile | top 20% | bottom 20% |
| d15 | top 15% | bottom 15% |
| decile | top 10% | bottom 10% |
| ventile | top 5% | bottom 5% |
| pos_vs_neg | CLDN4 count > 0 | CLDN4 count = 0 |

`expr_prop` is 0.05, 0.10, or 0.20. Ties in CLDN4 use ordinal ranks, the same
rule as the concordant-4 CellChat quartile.

A patient is scored for a gate when n malignant ≥ 40, both arms ≥ 10, and the
receiver has ≥ 20 cells. T/NK is the selection receiver. Myeloid is scored with
the same gate and is not used to pick the winner.

## Selection

A T/NK grid row is eligible when n ≥ 40, every cohort has n ≥ 3, and every
cohort mean is > 0 for both CellPhoneDB and Connectome. The winner is the
eligible row with the largest mean CellPhoneDB family Δ. Ties break to a higher
fraction of patients with Δ>0, then larger n, then `expr_prop` closer to 0.10,
then the milder tail.

Wilcoxon signed-rank is two-sided across patients. BH is across the T/NK grid,
including rows that are not eligible. That q is the multiplicity figure for the
search. It is not a single pre-specified test.

The both-arms sensitivity keeps edges that are nonzero on both malignant arms.
It is reported for the winner and is not the selection score.

## Labels

Same malignant and T/NK rules as the CellChat concordant-4 inventory.
The script stops if n malignant or n T/NK disagrees with that inventory.
TACSTD2 is not a gate. Aliases apply only when the official symbol is absent:
PVRL2→NECTIN2, JAM1→F11R.

## Reproduce

```bash
bash methods/liana_max_barrier_family_cldn4/scripts/download.sh /tmp/concordant4_raw
python3 methods/liana_max_barrier_family_cldn4/scripts/run_max_effect.py
```

GSE205335 RDS is double-gzipped. `export_gse205335.R` peels it and writes a
gene-subset Matrix Market.
