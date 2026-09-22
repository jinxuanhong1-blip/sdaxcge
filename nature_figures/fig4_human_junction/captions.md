# Fig. 4 — Human public TACSTD2-high malignant cells and apical junction

**Fig. 4 | Hallmark apical junction is the leading public ORA term in TACSTD2-high malignant cells.**

Concordant-4 malignant pseudobulk only (GSE123902, GSE131907, GSE205335, GSE189357). The split is within-cohort TACSTD2 expression quartile, not CLDN4. Differential expression is cohort-adjusted OLS on log2(TMM-CPM+1). Positive log2 fold change is higher in TACSTD2 Q4. Expression n = 64 units; the Q4 versus Q1 contrast is 15 versus 19. These panels are the public result in PR #741. They are not a private 8KL KEGG tight-junction rank 1.

**a**, Volcano of 24,083 genes. Dashed lines are the ORA query rule (nominal P < 0.01 and |logFC| > 0.25). TACSTD2 logFC = +3.90, P = 7.7 × 10⁻⁷, FDR = 0.019. It is the only gene with FDR < 0.05. Blue points are the 13 Hallmark apical-junction genes that overlap the up-gene query. Teal points are tight-junction ORA overlaps and the focal junction genes. CLDN4 (logFC +1.63, nominal P = 0.013) is labeled and is not in the apical-junction overlap.

**b**, Over-representation of 274 up genes (TACSTD2 held out) on a background of 24,083 genes. Hallmark apical junction is public ORA rank 1 (13/194 genes, enrichment 5.89, FDR 1.28 × 10⁻⁵). GO tight-junction organization is rank 11 (enrichment 5.94, FDR 8.94 × 10⁻³). KEGG tight junction is rank 31 (enrichment 2.28, FDR 0.20) and is not rank 1. Bubble color is −log10 FDR; size is overlap count. The same contrast, preranked GSEA (OLS t, TACSTD2 removed from the ranking, 1,000 permutations, seed 42), places Hallmark apical junction at NES +2.55, rank 8 of positive NES (FDR 1.46 × 10⁻³), and KEGG tight junction at NES +2.04, rank 22 (FDR 1.46 × 10⁻³). GSEA is not drawn as a separate panel.

**c**, Within-cohort median log2(TMM-CPM+1) for the same gene groups, row-centered within each gene for display (color truncated at ±2.5). Column n is the published Q1 or Q4 count. Numbers at the right are the published cohort-adjusted OLS logFC; bold marks nominal P < 0.01. Medians are descriptive. They are not a second differential-expression test. GSE189357 Q4 is two patients.

Source tables are copies from PR #741 (`methods/concordant4_tacstd2_malignant_deg/tables/`). The heatmap matrix is log2(TMM-CPM+1) from that page’s malignant UMI sums; TACSTD2 matches `units_tacstd2.tsv`.
