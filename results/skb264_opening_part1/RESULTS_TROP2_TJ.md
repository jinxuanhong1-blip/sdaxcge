# TROP2-high malignant cells enrich tight junction programs

TROP2-high malignant cells up-regulate a junction program. In the locked concordant-4 pool, over-representation of genes higher in the TACSTD2 upper quartile places Hallmark apical junction first among the tested sets. That is the top-pathway statement in this subsection. It is an over-representation rank on the malignant-cell contrast, not a claim that every GSEA universe ranks KEGG tight junction first, and it does not name a single junction gene.

## Design

The split is malignant TACSTD2, not a claudin. Cohorts are the locked four only: GSE123902, GSE131907, GSE205335, and GSE189357. The unit is the patient, donor, or sample. Malignant pseudobulk expression is log2(TMM-CPM+1). Within each cohort the upper quartile was compared with the lower quartile and the contrasts were stacked, giving 19 versus 15. The expression matrix has 64 units. The quartile contrast is not an n=65 differential-expression test. TACSTD2 was held out of the gene-set query. Private KL matrices were not used.

## Junction programs are higher, and apical junction is the top over-represented term

Family scores on the same model are higher in TACSTD2-high malignant cells: the 206-gene tight-junction list (TACSTD2 removed) has logFC +0.223 (FDR 0.004187), epithelial adhesion logFC +0.734 (FDR 0.004187), and keratin logFC +0.929 (FDR 0.00545). Interferon and MHC-I/APM family scores are not lower on this split.

Genes with p<0.01 and absolute logFC above 0.25 (274 up) were tested by over-representation against the pre-specified collection. Hallmark apical junction is rank 1 (overlap 13, enrichment 5.89, p=3.79×10⁻⁷, FDR 1.28×10⁻⁵). Hallmark epithelial–mesenchymal transition is rank 2 (enrichment 5.74) and GO keratinization is rank 3 (enrichment 13.98), at the same FDR. The top upregulated pathway in this test is therefore a junction term, with an EMT term and a keratinization term immediately behind it.

Preranked GSEA on the OLS t statistic (TACSTD2 dropped, 1,000 permutations) enriches the same programs and does not put them at NES rank 1. Among positive NES values, GO keratinization is rank 2 (NES +2.825), the epithelial-adhesion set is rank 7 (NES +2.585), Hallmark apical junction is rank 8 (NES +2.555), GO tight-junction organization is rank 20 (NES +2.118), and KEGG tight junction is rank 22 (NES +2.040). Those junction FDRs are 0.00146. The leading positive NES term is Hallmark TNFα signaling via NF-κB (NES +3.200). Enrichment of the junction sets is the GSEA result; rank 1 is not.

## Mouse epithelium, same direction

On public GEMM epithelium (GSE154977, GSE180963, GSE165641; private KL mice not merged), a frozen 7-gene tight-junction module was higher in Tacstd2-detected than in Tacstd2-undetected epithelial cells in 7 of 7 scored mice (mean Δ +0.2608). The one-sided binomial p on direction is 0.007812, and the across-mouse rank test against background is p=3.601×10⁻⁹. Tacstd2 was not used to call epithelium. This is a within-mouse module contrast, not a second copy of the human quartile test.

## What this subsection does not use

GSE137244 is bulk cell-line RNA-seq (five KL libraries versus five KP libraries). It is the locked KL-versus-KP measurement for Tacstd2 (Δ +3.238 on log2(FPKM+1), exact Mann–Whitney p=0.00794, complete separation), not a malignant-cell high-versus-low differential analysis. A median split on Tacstd2 in that matrix is the genotype contrast. Bulk Hallmark GSEA on that KL-versus-KP ranking does not place apical junction first.

The broad tight-junction score in concordant-4 does not track fewer T/NK cells (DerSimonian–Laird ρ=+0.137, p=0.60, N=64). Immune-negative associations of TACSTD2 itself are a separate result and are not inferred here from the junction score. No claudin is prioritized in this subsection.
