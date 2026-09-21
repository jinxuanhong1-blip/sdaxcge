# GSE274940 pathway sweep: multi-claudin null, not a CLDN4-only knockout

Mouse mammary EpH4, day 7, three wild-type libraries versus three Cldn-null libraries (GSM8462258–GSM8462263; PMID 41171911). The genotype is multiplex CRISPR of Cldn3, Cldn4, Cldn7, Cldn8, Cldn9, Cldn12, Cldn23, and Cldn25. It is not a CLDN4-only knockout. Cldn4 mRNA is not depleted: mean raw counts 15,074 versus 16,733 (111% of WT), 128% of WT CPM after library-size normalization (log2 fold-change +0.32, 95% CI −0.65 to +1.30, Welch p = 0.37). Cldn3 mRNA falls to 14% of WT CPM and Cldn7 to 37%. The RNA contrast below is that multi-claudin, Cldn3/Cldn7-led pattern.

No proteomics matrix is deposited. GEO for this series is the RNA-seq count file only. A PRIDE search for the paper returns no project. Zenodo [10.5281/zenodo.16886088](https://doi.org/10.5281/zenodo.16886088) is the authors' source data for figures (immunoblot images, immunofluorescence, genome traces, electron microscopy, physiology), not a protein table. The physiological barrier measurement in that deposit (Figure 4A, conductance, representative clone, n = 3) is 0.6, 0.6, and 0.6 mS/cm² in wild-type EpH4 and 40.6, 32.7, and 25.4 mS/cm² in the Cldn-null cells. Protein loss remains the immunoblot in the paper. These tables are RNA.

The sweep asked whether that barrier-null RNA moves NHEJ down and cGAS–STING / IFN / antigen presentation up. Every fit uses Cldn-null minus wild type. Numbers are from the deposited expected counts. Nothing here was chosen by flipping the sign after the fact.

## What is actually up, and survives FDR

Two transcripts go up in limma-voom, edgeR quasi-likelihood, and DESeq2 Wald, and pass Benjamini–Hochberg in at least one of those fits (13,905 genes with count ≥ 10 in at least two libraries).

| Gene | Role in the thesis | limma log2FC (q) | edgeR log2FC (q) | DESeq2 log2FC (q) |
|---|---|---:|---:|---:|
| Cgas | cGAS sensor, should rise | +1.68 (0.11) | +1.75 (0.089) | +1.76 (0.031) |
| Ifitm3 | ISG, should rise | +0.92 (0.018) | +0.93 (0.016) | +0.94 (2.0×10⁻⁵) |

Cgas is the FDR hit inside the cGAS–STING core. Its wild-type CPM is about 0.9 and the null CPM is about 3.1, so the fold-change is large on a low baseline. Ifitm3 is expressed (CPM about 20 to 38) and is the IFN-side hit that all three fits agree on.

## The same genes that block a pathway call

The 4-gene core (Cgas, Sting1, Tbk1, Irf3) has mean limma log2 fold-change **−0.23**. Sting1 falls: limma −2.33 (q = 0.026), edgeR −2.12 (q = 0.041), DESeq2 −2.12 (q = 3.4×10⁻⁴). Tbk1 and Irf3 are flat. Dropping Sting1 and retesting the other three, which is the leave-one-gene search, still does not clear a set test (best roast one-sided p = 0.33; z-mean exact permutation p = 0.45).

Hallmark interferon sets are enriched **down**, not up. Preranked GSEA, 2,000 gene-label permutations, all six libraries:

| Set | Ranking | NES | Permutation p | Thesis (up) |
|---|---|---:|---:|---|
| Hallmark IFN-γ | limma t | −1.60 | 0.001 | no |
| Hallmark IFN-α | limma logFC | −1.54 | 0.0065 | no |
| Compact ISG list | limma logFC | −1.73 | 0.003 | no |
| Hallmark IFN-γ | DESeq2 Wald | −1.45 | 0.0125 | no |
| Hallmark IFN-α | limma t | −1.35 | 0.044 | no |

Stat1 log2 fold-change is −0.09 and Isg15 is −0.15. Ifng, Cxcl9, Cxcl10, and Cxcl11 are essentially uncounted, as in the first pass.

Antigen presentation does not rise. The 16-gene MHC-I/APM core has preranked NES −1.61 on the DESeq2 Wald statistic (p = 0.028). Psmb8 is the extreme gene: log2 fold-change about −4.8 to −5.1, q < 0.005 in all three fits. B2m (+0.17 to +0.23) and H2-K1 (+0.10 to +0.13) do not move.

## NHEJ

The six-gene core (Prkdc, Lig4, Xrcc4, Xrcc5, Xrcc6, Nhej1) has mean limma log2 fold-change **−0.07**. Prkdc, Lig4, Xrcc5, and Nhej1 are slightly negative, Xrcc4 is about zero, Xrcc6 is +0.29 (p = 0.08). No core gene has q < 0.05. z-mean Welch p = 0.89. Exact reassignment of the six libraries gives permutation p = 0.40. Leaving Xrcc6 out, the gene that goes up, still leaves the best z-mean permutation p at 0.25.

The broad Reactome NHEJ set (67 genes, including histone genes) is enriched **up**. Limma logFC ranking: NES +2.00, 0 of 2,000 permutations (p < 1/2,000). Leading edge: H2bc8, H2bc6, Bard1, Brca1, H2bc4, and Polm. Polm itself is up in every fit (log2 fold-change +1.07, DESeq2 q = 1.1×10⁻⁵). That is an NHEJ polymerase going up, not the core ligase complex going down. Hallmark DNA repair mean log2 fold-change is −0.02.

## Least-null score that points the thesis way

On all six libraries, the smallest thesis-directed p in the pre-specified primary family (190 tests: camera, roast, fry, z-mean, GSVA, ssGSEA, and preranked GSEA on limma t, DESeq2 Wald, and signal-to-noise) is GSVA of GO cGAS–STING: score difference +0.12, Welch p = 0.14, one-sided p = 0.069, exact 20-assignment permutation p = 0.10. The minimum Benjamini–Hochberg q in that family is 1.00. The same GO set has mean limma log2 fold-change −0.05, and fry calls the direction Down (two-sided p = 0.054). GSVA and the mean log fold-change do not agree, and neither is a significant up.

Sensitivity that looks sharper does not survive the sample permutation. Dropping the noisiest quartile of genes in that GO set gives roast one-sided p = 0.038, but the z-mean exact permutation p stays 0.10. Dropping library KO2 gives a Reactome IFN-α/β z-mean Welch p = 0.028; with five libraries the exact permutation floor is 0.10, and that is the permutation p.

## Sweep

| Block | What was run |
|---|---|
| Fits | limma-voom (robust eBayes), edgeR quasi-likelihood, DESeq2 Wald with poscounts size factors. RSEM expected counts rounded to integers. |
| Rankings | limma t, limma logFC, edgeR signed √F, edgeR logFC, DESeq2 Wald, DESeq2 logFC, signal-to-noise. |
| Set tests | camera at estimated correlation, correlation 0, and 0.01; roast mean (4,999 rotations) in the thesis direction; fry; GSVA; ssGSEA; z-mean of log expression with a Welch test and an exact sample permutation. |
| Contrasts | all null vs WT; drop each library once. |
| Filters | count ≥ 10 in ≥ 2 libraries (main), ≥ 5 in ≥ 2, ≥ 10 in ≥ 3. |
| Set edits | as published; drop genes with mean CPM < 1; drop the highest-CV quartile inside the set; leave-one-gene-out for sets of 8 or fewer. |
| Gene sets | MSigDB mouse v2026.1 Hallmark IFN-α, IFN-γ, DNA repair; Reactome NHEJ, cytosolic DNA sensors, STING, IFN-α/β, IFN-γ, MHC-I; GO NHEJ, cGAS–STING, MHC-I peptide presentation, type-I IFN production; plus the six-gene NHEJ core, four-gene cGAS–STING core, compact ISG list, and MHC-I/APM core. |

fgsea 1.28 did not compile against the installed BH headers (a boost constexpr error). The preranked statistics are the weighted enrichment score with 2,000 gene-label permutations, the same null fgsea's simple procedure uses. That p can be small when a sample permutation cannot: with n = 3 vs 3 the exact reassignment p cannot go below 0.05 one-sided. Sample-level rows carry that permutation p in the `extra` column.

The merged table is `sweep_merged.tsv` (1,549 limma/edgeR/DESeq2/GSVA rows and 1,197 preranked rows). `sweep_headline.tsv` is one row per key set for the full contrast. Per-gene fits are `sweep_focus_genes.tsv`.

## What this does and does not say

Cgas and Ifitm3 are higher in the Cldn-null RNA, in three independent fits, and those two calls are the part of the barrier-loss → STING/IFN direction that the counts support. They sit next to Sting1 down, Psmb8 down, Hallmark IFN enriched down, and Reactome NHEJ enriched up. The six-gene NHEJ core is only a small negative point estimate. None of the pathway tests that were aimed at the thesis survive the primary-family correction, and the strongest pathway enrichments in the sweep go the other way.

Reproduce with `Rscript scripts/gse274940/sweep_pathways.R` and then `python3 scripts/gse274940/gsea_prerank.py` and `python3 scripts/gse274940/plot_sweep.py`. The first-pass Cldn4 count check remains `python3 scripts/gse274940/analyze_gse274940.py`.
