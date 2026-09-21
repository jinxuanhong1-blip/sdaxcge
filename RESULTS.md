# GSE207704 CLDN4 knockout: IFN, MHC-I/APM, and tight junctions

## Result

In GSE207704, CLDN4 knockout does not upregulate interferon or MHC-I/APM genes. Interferon Alpha Response is NES rank 43 of 50 (NES -1.49, mean log2FC -0.12, call weak). Interferon Gamma Response is NES rank 37 of 50 (NES -1.40, mean log2FC -0.08, call weak). The MHC-I/APM panel is null (mean log2FC 0.03). The tight-junction panel, with CLDN4 removed, is down (unshrunk mean log2FC -0.32, shrunk mean -0.05). Negative NES means the gene set is enriched among genes that fall after knockout.

This is a breast-cancer cell-line knockout (T47D and MCF7), not a lung-cancer knockdown. It does not re-estimate the private CLDN4-knockdown IFN/APM result.

## Question

A private CLDN4 knockdown was described as upregulating IFN/APM (347 genes up, 736 down; 49/52 IFN/APM genes shared with an SKB264 signature, r = 0.52). No public lung-cancer CLDN4 knockdown or knockout transcriptome is deposited. The nearest public CLDN4 loss-of-function RNA-seq series is GSE207704 (Kage et al., Breast Cancer Research 2023, PMID 37059993): CRISPR knockout of CLDN4 in T47D and MCF7, two replicates each, compared with wild type.

## Data

Eight single-end NovaSeq libraries, SRA SRR20029118–SRR20029125 (BioProject PRJNA856719). SRA runinfo reports a mean read length of 51 bp and about 44–47 million reads per library. The paper text describes index-trimmed 100 bp reads; the public runs are 51 bp.

Reads were extracted with sratoolkit 3.4.1 `fastq-dump` in chunks of at most 5 million spots (1 million when a chunk aborted), then concatenated. Full-file `fasterq-dump` and `fastq-dump` runs segfaulted in this environment; the merged files contain the same spot counts as the SRA runinfo. Reads were quantified with kallisto 0.52.0 against Ensembl release 116 GRCh38 cDNA (`Homo_sapiens.GRCh38.cdna.all.fa.gz`), k = 31, single-end fragment prior `-l 200 -s 30`. Transcript estimated counts were summed to gene symbols. A low transcript-level unique-mapping rate is expected for 51 bp reads across isoforms; the differential-expression input is the gene-level sum of kallisto estimated counts. When two Ensembl genes shared a symbol, the gene with more total reads was kept (2494 symbols had a dropped gene id).

Pseudoalignment rate:

| Sample | Genotype | Reads | Pseudoaligned |
|---|---|---:|---:|
| T47D_WT_1 | T47D WT | 44,062,012 | 39.0% |
| T47D_WT_2 | T47D WT | 45,214,824 | 38.9% |
| T47D_KO_1 | T47D KO | 44,424,212 | 41.0% |
| T47D_KO_2 | T47D KO | 44,387,540 | 38.9% |
| MCF7_WT_1 | MCF7 WT | 46,699,152 | 39.3% |
| MCF7_WT_2 | MCF7 WT | 45,458,325 | 40.3% |
| MCF7_KO_1 | MCF7 KO | 47,026,612 | 43.9% |
| MCF7_KO_2 | MCF7 KO | 44,087,640 | 44.1% |

The GEO supplementary file `GSE207704_CLDN4_RNAseq.txt.gz` is a Cufflinks-style FPKM matrix with four columns (the two replicates are already averaged) and it does not contain HLA-A, HLA-B, B2M, TAP1, TAP2, or MX1. It is used only as a sensitivity check on the genes it does contain.

## Statistics

Counts were analyzed with PyDESeq2 0.5.4. The primary model is `~ cell_line + condition` on all eight samples, log2 fold change = KO / WT, wild type as the condition reference. T47D and MCF7 were also fit separately (`~ condition`, two versus two). Cook's filtering was off, matching DESeq2 when a group has fewer than three replicates. Wald statistics used maximum-likelihood log2 fold changes. apeglm shrinkage was applied only to the reported shrunk log2 fold change.

Hallmark enrichment is preranked GSEA (gseapy 1.3.1, fgsea multilevel, 1000 permutations, seed 207704) on the Wald statistic. The gene sets are Enrichr `MSigDB_Hallmark_2020`, the Hallmark collection of Liberzon et al. (Cell Systems 2015), with the large sets capped at 200 genes in that export. FDR is across the 50 Hallmark sets. A second rank is the hypergeometric overlap of each set with the top 10% of genes by the same Wald statistic.

A set is called **up** or **down** only when three things agree: the combined-model mean |log2FC| is at least 0.25, the Hallmark (or panel) GSEA FDR is below 0.05 with a matching NES sign, and the median log2FC is in that same direction in both cell lines. Opposite medians of at least 0.25 log2FC are **discordant**. A significant GSEA with a smaller mean is **weak**. A mean past the gate without GSEA FDR support is **suggestive**. Otherwise the call is **null**. These gates were written into the script before the counts were scored.

## Knockout check

CLDN4 combined log2FC = -0.89 (shrunk -0.87, Wald p = 3.67e-27, FDR = 3.50e-24). T47D log2FC = -1.06; MCF7 log2FC = -0.76. Author FPKM log2FC (pseudocount 0.1) was -1.07 in T47D and -0.73 in MCF7. CLDN4 RNA is lower in the knockout in both lines, so the sample labels and the quantification agree with the genotype.

## Global differential expression

Combined model, 18,812 genes with total count at least 10. Uncorrected p < 0.05: 1,864 up and 1,906 down. FDR < 0.05: 1,010 up and 1,060 down. FDR < 0.10: 2,752 genes. The combined model has four knockout and four wild-type libraries, so single-gene FDRs are not empty. The per-line models have two replicates and much less power.

Spearman correlation of unshrunk log2FC, T47D versus MCF7, genes tested in both: r = -0.005, p = 0.548, n = 16,563. Among genes with combined FDR < 0.05 that were tested in both lines, r = 0.775 (p <1e-300, n = 2,056, 99.6% same sign). The genome-wide correlation is near zero because most genes do not move. The genes the combined model calls significant do move in the same direction in both lines.

## Hallmark rank

NES rank 1 is the Hallmark most enriched among genes that rise in the knockout. Interferon Alpha Response is NES rank 43 of 50 (NES -1.49, FDR 0.035; overlap-up rank 37 of 50; 7/92 genes in the top 10%; 17/92 in the bottom 10%; bottom-decile FDR 0.029). Interferon Gamma Response is NES rank 37 of 50 (NES -1.40, FDR 0.039; overlap-up rank 34 of 50; 16/182 genes in the top 10%; 28/182 in the bottom 10%; bottom-decile FDR 0.037).

| NES rank | Hallmark | NES | FDR | Overlap-up rank | Top-decile overlap | Call |
|---:|---|---:|---:|---:|---:|---|
| 1 | Myc Targets V1 | 2.07 | 1.48e-07 | 1 | 54/199 | weak |
| 2 | mTORC1 Signaling | 1.62 | 0.007 | 2 | 48/199 | weak |
| 3 | Protein Secretion | 1.56 | 0.024 | 4 | 22/95 | weak |
| 4 | E2F Targets | 1.34 | 0.073 | 7 | 34/200 | null |
| 5 | Estrogen Response Early | 1.31 | 0.073 | 3 | 40/198 | null |
| 11 | Estrogen Response Late | 1.20 | 0.157 | 8 | 33/195 | null |
| 19 | Inflammatory Response | -0.80 | 0.931 | 44 | 12/175 | null |
| 20 | Allograft Rejection | -0.84 | 0.880 | 25 | 18/157 | null |
| 25 | Complement | -1.04 | 0.413 | 36 | 15/180 | null |
| 37 | Interferon Gamma Response | -1.40 | 0.039 | 34 | 16/182 | weak |
| 39 | Apical Junction | -1.42 | 0.035 | 49 | 11/184 | weak |
| 40 | TNF-alpha Signaling via NF-kB | -1.47 | 0.026 | 46 | 12/182 | weak |
| 41 | Bile Acid Metabolism | -1.47 | 0.035 | 42 | 7/102 | weak |
| 43 | Interferon Alpha Response | -1.49 | 0.035 | 37 | 7/92 | weak |
| 44 | Fatty Acid Metabolism | -1.50 | 0.035 | 28 | 16/149 | weak |
| 46 | Cholesterol Homeostasis | -1.54 | 0.033 | 24 | 9/73 | weak |
| 47 | UV Response Dn | -1.55 | 0.024 | 39 | 11/139 | weak |
| 48 | Apoptosis | -1.57 | 0.014 | 26 | 17/148 | weak |
| 49 | IL-6/JAK/STAT3 Signaling | -1.70 | 0.014 | 45 | 4/74 | weak |
| 50 | TGF-beta Signaling | -1.85 | 0.012 | 48 | 2/52 | weak |

Full 50-set tables are in `results/gse207704/`.

## Interferon and MHC-I/APM

| Set | Call | Mean log2FC | Median T47D | Median MCF7 | NES | FDR | Genes log2FC > 0 |
|---|---|---:|---:|---:|---:|---:|---:|
| Interferon Alpha Response | weak | -0.12 | -0.16 | 0.06 | -1.49 | 0.035 | 43% of 92 |
| Interferon Gamma Response | weak | -0.08 | -0.12 | -0.00 | -1.40 | 0.039 | 43% of 182 |
| MHC-I/APM | null | 0.03 | -0.22 | 0.02 | -1.20 | 0.230 | 35% of 20 |

Interferon-alpha is a weak coordinated shift (FDR < 0.05, mean |log2FC| < 0.25), and the enrichment is negative (mean log2FC -0.12). Interferon-gamma is a weak coordinated shift (FDR < 0.05, mean |log2FC| < 0.25), and the enrichment is negative (mean log2FC -0.08). The MHC-I/APM panel is null (mean log2FC 0.03). Negative NES means the set is enriched among genes that fall in the knockout.

MHC-I/APM genes on the panel: 20 of 22 were tested (total count ≥ 10 in the combined matrix). Mean shrunk log2FC = 0.03. Mann–Whitney p versus other genes = 0.488 (median minus background -0.10). Panel genes absent from the tested matrix: HLA-G, PSMB8.

Per-gene log2 fold changes for the interferon Hallmarks, the MHC-I/APM panel, and the tight-junction panel are in `results/gse207704/genes_ifn_apm_tj.tsv`.

Largest absolute combined Wald statistics inside the two interferon Hallmarks:

| Gene | Set membership | log2FC | Shrunk log2FC | FDR |
|---|---|---:|---:|---:|
| ARL4A | gamma | -0.64 | -0.63 | 8.95e-15 |
| CDKN1A | gamma | 0.65 | 0.63 | 6.13e-13 |
| NFKBIA | gamma | -0.44 | -0.42 | 1.56e-08 |
| PARP9 | alpha | -0.57 | -0.54 | 2.31e-07 |
| PNPT1 | alpha+gamma | 0.33 | 0.32 | 1.07e-06 |
| PSME1 | alpha+gamma | -0.38 | -0.36 | 2.99e-06 |
| OAS3 | gamma | -0.67 | -0.62 | 1.61e-05 |
| IFI44 | alpha+gamma | -2.08 | -1.90 | 2.02e-05 |
| ARID5B | gamma | -0.40 | -0.38 | 2.67e-05 |
| BST2 | alpha+gamma | -3.15 | -2.78 | 8.60e-05 |
| FAS | gamma | 0.93 | 0.83 | 1.04e-04 |
| HLA-A | gamma | 1.20 | 1.06 | 3.88e-04 |
| SRI | gamma | -0.40 | -0.36 | 6.21e-04 |
| RNF213 | gamma | -0.29 | -0.27 | 8.79e-04 |
| USP18 | alpha+gamma | -0.44 | -0.39 | 9.96e-04 |

## Tight junctions

The tight-junction panel excludes CLDN4. Call: **down**. Unshrunk mean log2FC = -0.32; apeglm-shrunk mean = -0.05 (median T47D -0.29, median MCF7 -0.17). The shrunk mean is closer to zero because several low-count claudins have large unshrunk fold changes. NES = -1.76, panel GSEA FDR = 0.008 (this FDR is within the custom panels, not the Hallmark family). 28 of 31 panel genes were tested. Absent: CLDN6, CLDN8, CLDN17.

Hallmark Apical Junction, which is broader than the claudin/TJ list, is NES rank 39 (NES -1.42, FDR 0.035, call weak).

Epithelial context genes (not part of the TJ call):

| Gene | Combined log2FC | T47D log2FC | MCF7 log2FC |
|---|---:|---:|---:|
| EPCAM | 0.21 | 0.32 | 0.08 |
| TACSTD2 | -0.84 | -0.65 | -1.04 |
| CDH1 | -0.06 | 0.12 | -0.26 |
| KRT8 | -0.01 | 0.19 | -0.24 |
| KRT18 | -0.14 | 0.42 | -0.71 |
| KRT19 | 0.15 | -0.15 | 0.43 |
| CLDN4 | -0.89 | -1.06 | -0.76 |

## Author FPKM sensitivity

On genes present in both the kallisto model and the deposited FPKM table, Spearman correlation of the two-line mean log2FC is r = 0.845 (p <1e-300, n = 10,752). Author-matrix GSEA on that mean FPKM log2FC places Interferon Alpha Response at NES rank 47 (NES -1.59) and Interferon Gamma Response at NES rank 46 (NES -1.49).

Of the interferon-alpha Hallmark, 37 genes are missing from the deposited table. Of interferon-gamma, 85 are missing. Of the MHC-I/APM panel, 13 are missing (HLA-A, HLA-B, HLA-E, HLA-F, HLA-G, B2M, TAP1, TAP2, TAPBP, PSMB8, PSMB9, ERAP2, NLRC5). A conclusion drawn only from that table would not have measured those genes.

## Limits

Two biological replicates per genotype. The cultures have no immune cells, so this is tumor-cell-intrinsic transcription only. The perturbation is a coding-sequence knockout in two ER-positive breast lines, not a knockdown in lung cancer. Ensembl 116 cDNA does not include every noncoding transcript. Single-end 51 bp reads limit isoform assignment; gene-level sums are the endpoint. The private knockdown counts and the SKB264 correlation were not recomputed here.

## Files

- `scripts/quantify_gse207704.sh` — kallisto index and quantification
- `scripts/analyze_gse207704.py` — DESeq2, GSEA, tables, this document
- `resources/MSigDB_Hallmark_2020.gmt` — Hallmark sets used for the ranks
- `results/gse207704/` — count matrix, DE tables, GSEA tables, figures

