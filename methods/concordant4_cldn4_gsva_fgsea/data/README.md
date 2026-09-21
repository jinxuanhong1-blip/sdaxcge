# Inputs (locked concordant-4 extracts)

Patient / donor / sample pseudobulk UMI sums of malignant cells, plus the
quartile labels already assigned on malignant CLDN4 % positive.

Copied from the concordant-4 CLDN4 analysis so this script does not re-gate
cells and does not re-cut quartiles.

- `GSE*_malignant_counts.tsv.gz` — genes × units, UMI sums.
- `tnk_units.tsv` — one row per unit: cohort, CLDN4 %pos, within-cohort
  quartile, and `in_count_matrix` (P4001 is Q1 on the T/NK vector and is
  absent from the GSE205335 malignant sum).
- `a8_sets.json` — frozen Hallmark / KEGG / GO / custom MHC-I lists.

Not a mega-merge. Not GSE148071, GSE127465, GSE207422, GSE154826, or CD45+
extracts.
