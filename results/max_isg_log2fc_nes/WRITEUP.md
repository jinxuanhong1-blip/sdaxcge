# Largest Hallmark ISG log2FC and NES, by cell line

DESeq2 and pre-ranked GSEA were rerun on the public count matrices. The numbers below are from that rerun. Each maximum belongs to one cell line and one contrast.

ISGs are the union of MSigDB Hallmark Interferon Alpha Response (97 genes) and Interferon Gamma Response (200 genes), Enrichr `MSigDB_Hallmark_2020`. The reported gene is the largest DESeq2 log2 fold change among those genes with baseMean ≥ 10 and Benjamini–Hochberg padj < 0.05. NES is gseapy prerank on the Wald statistic (1,000 permutations, seed 123). Ku/DNA-PKcs/53BP1 symbols are not members of these two Hallmark sets, so they were not stripped by hand.

## Primary maxima

| Cell line | Contrast | Largest ISG | log2FC | padj | Largest NES |
|---|---|---|---:|---:|---|
| HCT116 | Ku80-AID, CDKi day 4, GSE294709 | XAF1 | +6.56 | 2.5×10⁻²⁵ | IFN-α 3.36 (FDR 0) |
| HCT116 | DNA-PKcs KO, normoxia, GSE285698 | IFI44 | +3.93 | 0.020 | IFN-α 1.72 (FDR 0.0016) |
| MCF-7 | 53BP1 KO, untreated, GSE84986 | IFIT1 | +2.60 | 4.7×10⁻⁴ | IFN-α 2.82 (FDR 0) |

These three rows are different cell lines or different genes. They are not one NHEJ effect.

## HCT116 Ku80, GSE294709

Proliferation-matched degron: IAA+dox+CDKi versus CDKi alone, day 4, AID clones 1, 7 and 9. Design is clone plus condition.

XAF1 log2FC +6.56 (SE 0.60, baseMean 110). Raw counts are 1, 6, 2 versus 74, 436, 115. All three clones rise. XAF1 is in the interferon-gamma Hallmark only. The next genes are IFI44L +5.75, OAS1 +5.22 (baseMean 434, padj 1.0×10⁻⁷⁶), IFI44 +4.84, IFIT1 +4.22 (padj 1.5×10⁻¹⁷⁵), and ISG15 +3.35 (padj 3.4×10⁻⁶⁹).

Hallmark IFN-α NES 3.36 (FDR 0, 92 genes in the ranked list). Hallmark IFN-γ NES 3.26 (FDR 0, 170 genes). The larger NES is interferon-alpha. XRCC5 mRNA log2FC is −0.22 (padj 0.24). The AID tag removes Ku80 protein; the mRNA is not the on-target readout.

Same HCT116 line, Ku86-flox interaction (4OHT effect in CreERT2 beyond 4OHT without Cre): largest ISG is IFI27 +2.74 (padj 2.0×10⁻⁹, baseMean 155). IFN-α NES 2.63 (FDR 0). XRCC5 log2FC −3.19. The unadjusted Cre-versus-ethanol contrast is larger (SAMD9L +3.44, IFN-α NES 2.73) and is not the Ku-specific estimate: 4OHT without Cre still has IFN-α NES 1.91, with SECTM1 as its largest significant ISG at +0.86.

Same series, HEK293 Ku70 dox withdrawal (clones Sa11, SB, TI): IFN-α NES −1.32 (FDR 0.080) and IFN-γ NES −1.00 (FDR 0.46). XRCC6 log2FC −2.74. SOCS1 is +2.71 (padj 6.7×10⁻⁸) inside a set whose NES is negative. That row is HEK293, not HCT116.

## HCT116 DNA-PKcs knockout, GSE285698

Normoxia, H2O vehicle, knockout versus wild type, n=3. Library ids A1/B1/C1 are wild type and A2/B2/C2 are knockout (GEO sample descriptions). PRKDC log2FC −1.96 (padj 7.9×10⁻⁵⁷).

HLA-G has a larger point estimate, log2FC +6.84, and it fails the padj cutoff (padj 0.10, SE 3.32). Its counts are 0, 0, 0, 0, 0, 68: one knockout library. It is not the reported maximum.

The largest gene that passes padj < 0.05 is IFI44, log2FC +3.93 (SE 1.42, baseMean 13, padj 0.020). Counts are 0, 5, 0 versus 12, 11, 56. The direction is shared; the fold change is pulled by the third knockout library. CMPK2 +2.19 (baseMean 121, padj 8.3×10⁻⁸; counts 50, 54, 27 versus 171, 213, 194) and IFIT1 +1.70 (baseMean 207, padj 2.9×10⁻¹⁰; counts 80, 125, 93 versus 263, 319, 349) are the better-measured rises. ISG15 moves the other way: log2FC −0.69 (padj 0.0045).

IFN-α NES 1.72 (FDR 0.0016). IFN-γ NES 1.64 (FDR 0.0025). The interferon-gamma leading edge includes MT2A and CDKN1A along with IFIT1.

CoCl2 (200 µM, 12 h), same knockout: IFI44 log2FC +7.42 (padj 0.0017; counts 0, 0, 0 versus 44, 42, 22) and IFN-α NES 1.68 (FDR 0.0039). ISG15 log2FC −1.19. That arm is hypoxia, not the unstressed knockout.

## MCF-7 53BP1 knockout, GSE84986

Untreated wild type (n=3) versus two null clones (n=3 each, pooled). TP53BP1 log2FC −0.96 (padj 2.9×10⁻¹²). The six knockout columns are two clones of one parental line.

Largest ISG: IFIT1 log2FC +2.60 (SE 0.53, baseMean 45, padj 4.7×10⁻⁴). Next are IFIT2 +2.43, MX1 +2.02, and OAS1 +1.93. ISG15 log2FC +0.58 (padj 0.68). Clone 1 alone: IFIT1 +2.64, IFN-α NES 2.58. Clone 2 alone: IFIT2 +2.60, IFN-α NES 2.66. Both clones move.

IFN-α NES 2.82 (FDR 0, 90 genes). IFN-γ NES 2.42 (FDR 0, 160 genes). The 5 Gy, 4 h arm in the same line: IFIT1 +2.56 and IFN-α NES 2.70. That arm is irradiated.

## Methods

- GSE294709 HCT116 counts: `GSE294709_RNAseq-counts-on-gene-AMPSEQ.txt.gz`. HEK293 arm of the same series: kallisto estimated counts, rounded to integers.
- GSE285698: `GSE285698_raw_counts.txt.gz`. Symbols are the text after the Ensembl id.
- GSE84986: per-sample gene counts in `GSE84986_RAW.tar`, summed to HGNC symbols with the HGNC Ensembl map. Genes absent from that map are not in the MCF-7 test.
- Genes with fewer than 10 counts across the contrast samples were dropped before DESeq2. Size factors use `poscounts`. The log2 fold change is the Wald estimate, not an apeglm shrink.
- Ranking for GSEA adds 10⁻⁴ times the log2 fold change so tied Wald statistics have a stable order. gseapy still reported a small fraction of exact ties (0.03% to 5%, higher when n=2). For the HCT116 Ku80 day-4 contrast the IFN-α NES matches an independent prerank of the same contrast at 3.36.
- Residual degrees of freedom are below 3 in the clone-adjusted fits, so dispersion priors are coarse. Sample sizes are 2–3 per arm except the pooled MCF-7 contrast.

Reproduce: `python3 scripts/max_isg_log2fc_nes/analyze.py`

Tables: `tables/contrast_maxima.tsv`, `tables/top_isg_padj05.tsv`, `tables/gsea_nes.tsv`, `tables/isg_log2fc.tsv.gz`.
