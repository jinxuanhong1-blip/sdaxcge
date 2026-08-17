# Results (this run)

| File | What |
| --- | --- |
| `n_table.tsv` | Every n used in FINDING |
| `sample_metrics.tsv` | 12 post-tx patients |
| `patient_level_tests.tsv` | E1/E2 Wilcoxon and Spearman |
| `tnk_patient_gene_means_NMPR_vs_MPR.tsv` | Gene-wise patient means |
| `ligand_activity_all.tsv` | All settings × gene sets |
| `top_ligands.tsv` | Top 15 per setting × set |
| `top_ligands_primary.tsv` | Pooled CLDN4-high → all T/NK, IFN and cytotoxicity only |
| `cldn4_high_ligand_patient_means.tsv` | Ligand `log1p(CP10k)` in CLDN4-high cells |
| `cldn4_high_ligand_NMPR_vs_MPR.tsv` | Patient Wilcoxon on those means |
| `ligand_vs_tnk_program_spearman.tsv` | Ligand vs T/NK scores, n=12 |
| `combinatorial_ifn_activity_MPR_vs_NMPR.tsv` | Shared-ligand Pearson (identical if in both sets) |
| `summary.json` | Machine-readable n + top lists |
| `fig1_patient_programs.png` | CLDN4 / IFN / cytotoxicity by MPR |
| `fig2_cldn4_vs_tnk_ifn.png` | Patient scatter |
| `fig3_ligand_activity_top.png` | Prior Pearson bars |
| `fig4_combinatorial_MPR_vs_NMPR.png` | MPR vs NMPR IFN activity (diagonal) |
