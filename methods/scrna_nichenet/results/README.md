# Results (this run)

| File | What |
| --- | --- |
| `n_table.tsv` | Every n used in FINDING |
| `sample_metrics.tsv` | 12 post-tx patients |
| `patient_level_tests.tsv` | E1/E2 Wilcoxon and Spearman |
| `tnk_patient_gene_means_NMPR_vs_MPR.tsv` | Gene-wise patient means |
| `ligand_activity_all.tsv` | All settings × gene sets |
| `top_ligands.tsv` | Top 15 per setting × set |
| `top_ligands_primary.tsv` | Pooled high-barrier → all T/NK, a priori sets only |
| `high_barrier_ligand_patient_means.tsv` | Ligand `log1p(CP10k)` in high-barrier cells |
| `high_barrier_ligand_NMPR_vs_MPR.tsv` | Patient Wilcoxon on those means |
| `ligand_vs_tnk_program_spearman.tsv` | Ligand vs T/NK scores, n=12 |
| `combinatorial_exhaustion_activity_MPR_vs_NMPR.tsv` | Shared-ligand Pearson (identical if in both sets) |
| `summary.json` | Machine-readable n + top lists |
| `fig1_patient_programs.png` | Barrier / cytotoxicity / exhaustion by MPR |
| `fig2_barrier_vs_exh_minus_cyto.png` | Patient scatter |
| `fig3_ligand_activity_top.png` | Prior Pearson bars |
| `fig4_combinatorial_MPR_vs_NMPR.png` | MPR vs NMPR activity (diagonal) |
