# Sweep — GSE207704 c-NHEJ / STING / IFN

Status: **FINAL discordant cancer-line KD**

Positive log2FC means higher after CLDN4 knockout. A thesis slice, fixed before reading these set scores, is a c-NHEJ set and a STING/IFN set in the **same** quantification and the **same** contrast, each with at least 3 measured genes, mean log2FC at least 0.10 in the thesis direction, and a majority of genes on that side. Sign-only results inside ±0.10 are near misses.

## What was tried

- Cufflinks group-mean FPKM (the GEO supplementary table), T47D, MCF7, and the mean of the two lines.
- Open SRA runs SRR20029118–SRR20029125 (public S3 `sra-pub-run-odp`). ENA HTTPS FASTQ failed an SSL handshake. `fasterq-dump` segfaulted; `fastq-dump` in 8 million-read chunks did not. Spots are 51 bp single-end, not the 100 bp length named in the paper. Kallisto 0.51.1 on Ensembl 90 cDNA, single-end fragment prior `-l 200 -s 20`, transcript counts summed to gene symbols. Set tests use size-factor log2FC and drop genes with mean estimated count < 10 in that contrast. PyDESeq2 Wald is per line.
- HGNC aliases for STING1 (approved symbol, previous symbol TMEM173, aliases STING, ERIS, MPYS, MITA, NET23, FLJ38577) and CGAS (MB21D1, C6orf150), matched to symbols actually in each matrix. The cufflinks locus window chr5:139475534–139482790 (Ensembl 90 TMEM173) has no row.
- Hallmark IFN-α, IFN-γ, and DNA repair from MSigDB 2023.1, 2023.2, 2024.1, and 2025.1, plus the repo’s Hallmark 2020 lists. IFN-α membership is identical across 2023.1–2025.1. IFN-γ differs by METTL7B (2023.1) versus TMT1B (2024.1/2025.1).
- Reactome 2022 via Enrichr: NHEJ, alt-NHEJ/MMEJ, STING, cytosolic DNA, IFN-α/β, IFN-γ, IFN signaling, DSB response. Reactome STING includes PRKDC, XRCC5, and XRCC6; those three were also dropped in a separate set.
- QuickGO GO:0006303 (NHEJ), GO:0060337 (type I IFN signaling), GO:0035456 (response to IFNβ).
- QuickGO GO:0031040 (micronucleus) and GO:0032125 (micronucleus organization): the terms exist and returned **zero** annotations, so no micronucleus gene set was scored.
- Methods other than GSEA: mean, median, sign test, Mann–Whitney against the rest of the rank, and Fisher exact overlap with genes at log2FC < −0.25 or > +0.25. Prerank GSEA (p=1, 1000 perms, seed 42) only when at least 8 genes are in the rank.
- Line emphasis is the T47D column versus the MCF7 column, not a pooled test.

## STING1 aliases found

| source | alias | found_as |
|---|---|---|
| kallisto_ensembl90_norm_log2FC | TMEM173 | TMEM173 |

## Focused set scores

| source | set | contrast | n_measured | mean_log2FC | n_pos | n_neg | direction | nes |
|---|---|---|---|---|---|---|---|---|
| cufflinks_group_mean_FPKM | APM_21_repo | MCF7 | 9 | 0.088 | 4 | 5 | mixed_or_flat | 0.862 |
| cufflinks_group_mean_FPKM | APM_21_repo | T47D | 9 | -0.206 | 2 | 7 | down | -1.319 |
| cufflinks_group_mean_FPKM | GO_NHEJ_0006303 | MCF7 | 68 | 0.052 | 38 | 30 | sign_up_below_bar | 0.801 |
| cufflinks_group_mean_FPKM | GO_NHEJ_0006303 | T47D | 68 | 0.006 | 38 | 30 | sign_up_below_bar | 0.490 |
| cufflinks_group_mean_FPKM | STING_core_5 | MCF7 | 4 | 0.080 | 3 | 1 | sign_up_below_bar | NA |
| cufflinks_group_mean_FPKM | STING_core_5 | T47D | 4 | 0.036 | 3 | 1 | sign_up_below_bar | NA |
| cufflinks_group_mean_FPKM | enzymatic_cNHEJ_5 | MCF7 | 5 | 0.044 | 3 | 2 | sign_up_below_bar | NA |
| cufflinks_group_mean_FPKM | enzymatic_cNHEJ_5 | T47D | 5 | 0.262 | 3 | 2 | up | NA |
| cufflinks_group_mean_FPKM | hallmark_DNA_REPAIR_2025.1 | MCF7 | 135 | 0.038 | 82 | 53 | sign_up_below_bar | 0.657 |
| cufflinks_group_mean_FPKM | hallmark_DNA_REPAIR_2025.1 | T47D | 135 | -0.043 | 69 | 66 | mixed_or_flat | -0.839 |
| cufflinks_group_mean_FPKM | hallmark_IFNa_2023.1 | MCF7 | 63 | 0.013 | 35 | 28 | sign_up_below_bar | -1.031 |
| cufflinks_group_mean_FPKM | hallmark_IFNa_2023.1 | T47D | 63 | -0.319 | 19 | 44 | down | -1.807 |
| cufflinks_group_mean_FPKM | hallmark_IFNa_2025.1 | MCF7 | 63 | 0.013 | 35 | 28 | sign_up_below_bar | -1.030 |
| cufflinks_group_mean_FPKM | hallmark_IFNa_2025.1 | T47D | 63 | -0.319 | 19 | 44 | down | -1.819 |
| cufflinks_group_mean_FPKM | hallmark_IFNg_2025.1 | MCF7 | 118 | -0.035 | 61 | 57 | mixed_or_flat | -1.204 |
| cufflinks_group_mean_FPKM | hallmark_IFNg_2025.1 | T47D | 118 | -0.182 | 46 | 72 | down | -1.523 |
| cufflinks_group_mean_FPKM | reactome_NHEJ | MCF7 | 43 | 0.026 | 20 | 23 | mixed_or_flat | 0.595 |
| cufflinks_group_mean_FPKM | reactome_NHEJ | T47D | 43 | 0.043 | 26 | 17 | sign_up_below_bar | 0.831 |
| cufflinks_group_mean_FPKM | reactome_STING | MCF7 | 9 | 0.083 | 7 | 2 | sign_up_below_bar | 0.599 |
| cufflinks_group_mean_FPKM | reactome_STING | T47D | 9 | -0.006 | 6 | 3 | mixed_or_flat | -0.753 |
| cufflinks_group_mean_FPKM | reactome_STING_drop_NHEJ_genes | MCF7 | 6 | 0.120 | 5 | 1 | up | NA |
| cufflinks_group_mean_FPKM | reactome_STING_drop_NHEJ_genes | T47D | 6 | -0.061 | 4 | 2 | mixed_or_flat | NA |
| cufflinks_group_mean_FPKM | reactome_altNHEJ_MMEJ | MCF7 | 12 | 0.096 | 8 | 4 | sign_up_below_bar | 0.726 |
| cufflinks_group_mean_FPKM | reactome_altNHEJ_MMEJ | T47D | 12 | 0.049 | 10 | 2 | sign_up_below_bar | 0.675 |
| cufflinks_group_mean_FPKM | reactome_cytosolic_DNA | MCF7 | 9 | 0.077 | 5 | 4 | sign_up_below_bar | 0.669 |
| cufflinks_group_mean_FPKM | reactome_cytosolic_DNA | T47D | 9 | -0.167 | 4 | 5 | down | -1.193 |
| cufflinks_group_mean_FPKM | user_cNHEJ_7 | MCF7 | 7 | 0.088 | 5 | 2 | sign_up_below_bar | NA |
| cufflinks_group_mean_FPKM | user_cNHEJ_7 | T47D | 7 | 0.240 | 5 | 2 | up | NA |
| kallisto_ensembl90_norm_log2FC | APM_21_repo | MCF7 | 18 | 0.269 | 11 | 7 | up | 1.427 |
| kallisto_ensembl90_norm_log2FC | APM_21_repo | T47D | 18 | -0.098 | 3 | 15 | sign_down_below_bar | -0.841 |
| kallisto_ensembl90_norm_log2FC | GO_NHEJ_0006303 | MCF7 | 72 | 0.081 | 48 | 24 | sign_up_below_bar | 0.800 |
| kallisto_ensembl90_norm_log2FC | GO_NHEJ_0006303 | T47D | 73 | 0.003 | 40 | 33 | sign_up_below_bar | 0.557 |
| kallisto_ensembl90_norm_log2FC | STING_core_5 | MCF7 | 4 | 0.147 | 3 | 1 | up | NA |
| kallisto_ensembl90_norm_log2FC | STING_core_5 | T47D | 4 | 0.087 | 3 | 1 | sign_up_below_bar | NA |
| kallisto_ensembl90_norm_log2FC | enzymatic_cNHEJ_5 | MCF7 | 5 | 0.135 | 4 | 1 | up | NA |
| kallisto_ensembl90_norm_log2FC | enzymatic_cNHEJ_5 | T47D | 5 | 0.076 | 3 | 2 | sign_up_below_bar | NA |
| kallisto_ensembl90_norm_log2FC | hallmark_DNA_REPAIR_2025.1 | MCF7 | 144 | 0.086 | 107 | 37 | sign_up_below_bar | 0.802 |
| kallisto_ensembl90_norm_log2FC | hallmark_DNA_REPAIR_2025.1 | T47D | 144 | -0.057 | 73 | 71 | mixed_or_flat | -0.922 |
| kallisto_ensembl90_norm_log2FC | hallmark_IFNa_2023.1 | MCF7 | 77 | 0.076 | 44 | 33 | sign_up_below_bar | 0.965 |
| kallisto_ensembl90_norm_log2FC | hallmark_IFNa_2023.1 | T47D | 82 | -0.297 | 25 | 57 | down | -1.757 |
| kallisto_ensembl90_norm_log2FC | hallmark_IFNa_2025.1 | MCF7 | 77 | 0.076 | 44 | 33 | sign_up_below_bar | 0.990 |
| kallisto_ensembl90_norm_log2FC | hallmark_IFNa_2025.1 | T47D | 82 | -0.297 | 25 | 57 | down | -1.745 |
| kallisto_ensembl90_norm_log2FC | hallmark_IFNg_2025.1 | MCF7 | 153 | 0.051 | 81 | 72 | sign_up_below_bar | 1.020 |
| kallisto_ensembl90_norm_log2FC | hallmark_IFNg_2025.1 | T47D | 149 | -0.202 | 55 | 94 | down | -1.597 |
| kallisto_ensembl90_norm_log2FC | reactome_NHEJ | MCF7 | 48 | 0.039 | 27 | 21 | sign_up_below_bar | 0.614 |
| kallisto_ensembl90_norm_log2FC | reactome_NHEJ | T47D | 48 | -0.006 | 27 | 21 | mixed_or_flat | -0.427 |
| kallisto_ensembl90_norm_log2FC | reactome_STING | MCF7 | 13 | 0.108 | 9 | 4 | up | 0.741 |
| kallisto_ensembl90_norm_log2FC | reactome_STING | T47D | 13 | -0.069 | 5 | 8 | sign_down_below_bar | -0.756 |
| kallisto_ensembl90_norm_log2FC | reactome_STING_drop_NHEJ_genes | MCF7 | 10 | 0.114 | 7 | 3 | up | 0.795 |
| kallisto_ensembl90_norm_log2FC | reactome_STING_drop_NHEJ_genes | T47D | 10 | -0.113 | 3 | 7 | down | -0.905 |
| kallisto_ensembl90_norm_log2FC | reactome_altNHEJ_MMEJ | MCF7 | 12 | 0.149 | 9 | 3 | up | 0.861 |
| kallisto_ensembl90_norm_log2FC | reactome_altNHEJ_MMEJ | T47D | 12 | 0.009 | 9 | 3 | sign_up_below_bar | -0.646 |
| kallisto_ensembl90_norm_log2FC | reactome_cytosolic_DNA | MCF7 | 13 | 0.117 | 9 | 4 | up | 0.782 |
| kallisto_ensembl90_norm_log2FC | reactome_cytosolic_DNA | T47D | 12 | -0.159 | 4 | 8 | down | -1.127 |
| kallisto_ensembl90_norm_log2FC | user_cNHEJ_7 | MCF7 | 7 | 0.161 | 6 | 1 | up | NA |
| kallisto_ensembl90_norm_log2FC | user_cNHEJ_7 | T47D | 7 | 0.099 | 5 | 2 | sign_up_below_bar | NA |

## Thesis slices (mean |log2FC| ≥ 0.10, majority of genes, ≥3 genes, same source and contrast)

_No pair met the bar._

## Sign-only near misses (right signs, |mean| < 0.10)

_none_

Full rows: `tables/sweep_set_scores.tsv`, `tables/sweep_joint.tsv`, `tables/sweep_aliases.tsv`, `tables/kallisto_panel_counts.tsv`, `tables/kallisto_alignment.tsv`, `tables/deseq2_T47D.tsv`, `tables/deseq2_MCF7.tsv`.

## Closest calls

No c-NHEJ set was direction-down in either quantification. The joint table and the near-miss table are empty.

Kallisto pseudoalignment was 42.3–47.8% (51 bp spots, k = 31, transcriptome index, no genome decoy). CLDN4 itself is down with both knockout replicates below both wild-type replicates in both lines (PyDESeq2 log2FC T47D −1.038, MCF7 −0.762). PyDESeq2 warned that residual degrees of freedom are under 3, so the dispersion prior is weak. Shared replicate order is the result that matters. `padj` is secondary.

**T47D, after reprocessing, keeps the wrong-sign interferon result.** Hallmark IFN-α mean log2FC −0.297, NES −1.745, BH-FDR 0.0077 (82 genes in the rank). IFN-γ mean −0.202, NES −1.597, FDR 0.0077. Both knockout replicates are below both wild-type replicates for ISG15, BST2, OAS1, IFIT1, and MX1. PyDESeq2: ISG15 −1.963 (padj 4.7×10⁻²¹), BST2 −5.455 (padj 6.7×10⁻¹⁰). The user c-NHEJ panel mean is +0.099 (5 of 7 genes positive). PRKDC, LIG4, RIF1, TP53BP1, and XRCC5 are higher in both knockout replicates.

**MCF7 is the only place an interferon NES is positive, and it is not a thesis slice.** Hallmark IFN-α/γ means are +0.076 / +0.051 (below the 0.10 bar). NES +0.99 / +1.02, nominal p about 0.26–0.31, BH-FDR 0.56. STING-core mean +0.147 and Reactome STING mean +0.108 (NES +0.74, FDR 0.59) sit next to a c-NHEJ panel that is up (+0.161, 6 of 7 genes positive). PRKDC alone is lower in MCF7 (size-factor log2FC −0.106; PyDESeq2 −0.164, padj 3.6×10⁻⁵; both knockout replicates below both wild-type), while XRCC4, XRCC5, XRCC6, RIF1, and TP53BP1 are higher in both knockout replicates. One gene is not the set.

**STING1 was measured.** In Ensembl 90 it is TMEM173. T47D estimated counts are about 1–5 per sample (PyDESeq2 baseMean 2.5); that ratio is excluded by the mean-count floor of 10. MCF7 baseMean 9.9, log2FC +0.16, p = 0.85. CGAS (MB21D1) is modestly higher in both lines (T47D +0.36, MCF7 +0.26) with hundreds of counts. TBK1 and STAT1 stay flat. Hallmark versions 2023.1 through 2025.1 do not flip the T47D interferon sign (IFN-α membership is identical; IFN-γ swaps METTL7B for TMT1B).

QuickGO has GO:0031040 (micronucleus) and GO:0032125 (micronucleus organization) and returned zero gene annotations, so those sets were not scored. The scored DNA-sensing sets are Reactome cytosolic DNA and Reactome STING. Cytosolic DNA is down in T47D (mean −0.159) and only weakly up in MCF7 (mean +0.117), where NHEJ is also up.

This accession is a **FINAL discordant cancer-line KD** for KD → NHEJ down, STING/IFN up.
