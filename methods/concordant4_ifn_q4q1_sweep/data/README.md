# Inputs

Patient / donor / sample malignant UMI sums and the locked within-cohort
CLDN4 %pos quartiles. Same extracts as the concordant-4 GSVA/fgsea run.
P4001 is absent from the GSE205335 matrix.

- `GSE*_malignant_counts.tsv.gz` — genes × units.
- `tnk_units.tsv` — cohort, CLDN4 %pos, quartile, `in_count_matrix`.
- `locked_sets.json` — Hallmark IFN-alpha, Hallmark IFN-gamma, and the
  frozen 21-gene MHC-I panel. CLDN4 is not in these lists.
- `public_sets.gmt` — interferon, MHC, and chemokine sets copied from
  Enrichr libraries Reactome_2022, KEGG_2021_Human,
  WikiPathway_2021_Human, and GO_Biological_Process_2023. The second
  column is the library and the original set name.

Derived sets (Hallmark intersection, lineage-depleted Hallmark, epithelial
ISG/APM, inflammatory chemokine ligands) are built in `analyze.py` from
these files. They are not chosen from the differential-expression result.
