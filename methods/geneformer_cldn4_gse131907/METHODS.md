# Methods — CLDN4-high vs low malignant embeddings (GSE131907)

Pre-specified before the test statistics were computed. Numbers live in `FINDING.md`, which `analyze.py` writes from the tables.

## Cohort

GSE131907 (Kim et al.), author label `Cell_subtype == "Malignant cells"`, sample origin in {tLung, tL/B, mLN, PE, mBrain}, samples with at least 20 malignant cells. In this annotation the malignant subtype is present in mLN, tL/B, and mBrain only. tLung cells labeled tS1/tS2 are not reclassified as malignant.

The unit is the **sample** (the Kim sample id). It is not a patient id, and it is not the concordant-4 n=65.

## Cell draw

This is a high-vs-low contrast, not a prevalence sample.

1. On all author-malignant cells, CLDN4-high means raw UMI > 0 and CLDN4-low means raw UMI = 0. Full-sample %pos is reported from this census.
2. A sample enters the draw only if both arms have at least 5 cells.
3. Seed 1: up to 70 cells per arm, then QC (nUMI ≥ 500, nGenes ≥ 200, mitochondrial fraction < 0.20), then cap at 50 cells per arm.
4. The paired test uses samples that still have at least 5 embedded cells in each arm. Wilcoxon needs at least 6 such samples before a p-value is reported.

## Embeddings

**Geneformer V1-10M** (`ctheodoris/Geneformer`, gc30M dictionaries). Per cell: total UMI (all genes) scales counts to 10,000, divide by the pretrained nonzero median, rank detected genes in the vocabulary, truncate at 2048. No CLS/EOS (V1). Cell vector = mean of the second-to-last hidden state over gene tokens (official `emb_layer=-1`).

Primary tokenization **removes ENSG00000189143 (CLDN4) before ranking**. A second run keeps it.

**scGPT whole-human** (`wanglab/scGPT-human`). Within-cell 51-bin quantiles of nonzero raw counts (monotone, so log1p library-size normalization does not change bins). `<cls>` is prepended with value −2. Maximum length 1200. If more than 1199 genes fall in the vocabulary, a seeded subset is kept. Cell vector = L2-normalized final hidden state at `<cls>`. The checkpoint was trained with flash-attn; this run uses the same Wqkv / output / FFN weights with eager attention and post-norm ReLU blocks, dropout off. Primary run removes the CLDN4 token.

**Geneformer V2-104M_CLcancer** is benchmarked on 8 cells (CLS+EOS, max 4096, second-to-last layer). The full cohort is embedded only if that probe is at or under 1.5 s/cell. **V2-316M** is not loaded.

## Separation test

For each embedding, the direction for sample *s* is the mean of (centroid high − centroid low) over the other samples that pass the arm gate, L2-normalized. Cells in *s* are projected onto that direction. The paired difference is the mean projection of high cells minus low cells. Wilcoxon signed-rank, two-sided, across samples.

BH q-values for the foundation models are computed only across the holdout tests that were actually fit.

Sensitivity: residualize the V1 holdout projection on log1p(nUMI) within each sample and repeat the paired test.

## Programs

Mean log1p(CP10k) of the genes in the set that are present. CLDN4 is removed from every multi-gene set, including tight junction (KEGG ∪ GOBP organization ∪ focal genes that are not in the keratin list) and Hallmark apical junction. IFN is Hallmark interferon alpha ∪ gamma. MHC-I/APM is the custom list in `a8_sets.json`. Barrier/keratin is KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR. Epithelial states are AT1, AT2, club, basal, ciliated. The other sets are the Hallmark EMT, hypoxia, G2M, E2F, TNFα/NF-κB, and inflammatory response lists, plus the chemokine panel used in the concordant-4 script.

The 17 programs are one BH family. CLDN4 and TACSTD2 are controls and are not in that family.

A second test asks whether the program score tracks the Geneformer holdout axis: within-sample Spearman, then Wilcoxon of those coefficients across samples, BH across the 17 programs.

Nearest genes: mean within-sample Spearman of log1p(CP10k) vs the holdout projection, genes detected in at least 5% of embedded cells. A separate table lists Geneformer token states (same layer, seen in ≥30 cells) by cosine to the high−low vector.

## PCA and diffusion map

CLDN4 is removed from the feature matrix. 2000 genes with the largest variance of log1p(CP10k) among genes detected in ≥10 cells, z-scored, 30 principal components (full SVD). Diffusion map: 15-neighbor adaptive Gaussian kernel on those PCs, symmetrized, degree-normalized, 10 nontrivial components. Each component is sign-flipped so the median paired difference is non-negative; the two-sided p-value does not change. BH is separate for the 10 PCs and the 10 DCs.
