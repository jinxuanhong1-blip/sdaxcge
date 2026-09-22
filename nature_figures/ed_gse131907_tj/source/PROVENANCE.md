# Provenance

Patient-level tables are copied from PR #773 (`cursor/gse131907-maxeffect-tj-t-744c`), directory `methods/gse131907_maxeffect_tj_t/results/tables/`.

`tj_contrasts_for_figure.tsv` is a subset of `tj_sweep.tsv`, `tj_gene_sweep.tsv` and `tj_residual_on_epi.tsv` from that pull request. Rows with a `patient_tag` are recomputed from `tj_patient_deltas.tsv` when the figure is drawn, and the script stops if the recomputed mean does not match the published value.

`published_constants.tsv` records sweep-level counts that cannot be recomputed from the patient tables alone (Spearman family size, Bonferroni threshold, permutation P).

No cell-by-gene matrix is included. Inferential n is patients.
