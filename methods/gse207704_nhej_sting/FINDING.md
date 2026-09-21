# GSE207704 breast CLDN4 KO — c-NHEJ, STING, Hallmark IFN / APM

Public cufflinks table only (T47D and MCF7, CLDN4 CRISPR KO vs parental WT). The question is the direction against **KD → NHEJ down, STING/IFN up**.

c-NHEJ transcripts stay flat to slightly higher (0/7 consensus down). The four measured STING-core genes stay flat, and STING1 is absent from the deposit. Hallmark IFN-α and IFN-γ are lower after knockout in T47D (both BH-FDR < 0.05) and are not higher in MCF7. MHC-I/APM is mostly missing from the file; the eight genes that remain are not higher together.

## Design

Murakami et al., Breast Cancer Research 2023 (GEO [GSE207704](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207704)). Two breast lines, CLDN4−/− vs WT. GEO lists 2 biological replicates per genotype (GSM6310640–GSM6310647). The only expression file collapses each genotype to one FPKM. Ranks below are those group means. There is no gene-level sample FDR and no FASTQ re-quantification.

log2FC = log2((KO + 0.5) / (WT + 0.5)). One symbol can have more than one cufflinks locus; the reported value is the locus with the higher mean FPKM. A second column in `gene_panel.tsv` sums loci. UP / DOWN requires |log2FC| > 0.25 and max FPKM ≥ 1; otherwise FLAT. Consensus UP requires one line UP and the other not DOWN (same rule in the other direction).

CLDN4 QC, higher-mean locus: MCF7 -0.750 (FPKM 85.9 → 50.8), T47D -1.039 (FPKM 42.5 → 20.4). Both lines are DOWN. Residual mRNA is about 59% (MCF7) and 48% (T47D).

## Scorecard

| Axis | After CLDN4 KO, thesis says | What this file shows |
|---|---|---|
| c-NHEJ (PRKDC, XRCC4, LIG4, RIF1, TP53BP1, XRCC5, XRCC6) | down | **Not decreased.** 7/7 genes present. Both-line mean log2FC 0.164; 6/7 positive. Consensus DOWN: 0. Consensus UP: LIG4, PRKDC (T47D only). |
| STING core (CGAS, STING1, TBK1, IRF3, STAT1) | up | **Not increased.** STING1 / TMEM173 absent. The other four are FLAT (\|log2FC\| ≤ 0.25). Their mean log2FC is 0.058. |
| Hallmark IFN-α / IFN-γ | up | **Not increased.** T47D NES is negative (IFN-α -1.829, IFN-γ -1.567; both BH-FDR < 0.05). MCF7 NES is negative and not FDR < 0.05. |
| MHC-I / APM (21-gene list; Hallmark has no APM set) | up | **Not increased.** 8/21 genes are in the file. NES T47D -1.232, MCF7 0.868, neither FDR < 0.05. |

## 1. c-NHEJ

| Gene | Mapped | MCF7 log2FC | T47D log2FC | Mean | MCF7 | T47D | Both lines | vs thesis |
|---|---|---:|---:|---:|---|---|---|---|
| PRKDC | PRKDC | -0.249 | 0.324 | 0.038 | FLAT | UP | UP | not decreased (up in one line) |
| XRCC4 | XRCC4 | 0.238 | -0.038 | 0.100 | FLAT | FLAT | FLAT | not decreased |
| LIG4 | LIG4 | -0.046 | 1.038 | 0.496 | FLAT | UP | UP | not decreased (up in one line) |
| RIF1 | RIF1 | 0.241 | 0.168 | 0.205 | FLAT | FLAT | FLAT | not decreased |
| TP53BP1 | TP53BP1 | 0.155 | 0.200 | 0.178 | FLAT | FLAT | FLAT | not decreased |
| XRCC5 | XRCC5 | 0.168 | 0.164 | 0.166 | FLAT | FLAT | FLAT | not decreased |
| XRCC6 | XRCC6 | 0.108 | -0.181 | -0.037 | FLAT | FLAT | FLAT | not decreased |

Call rule on the both-line mean: 6 positive, 1 negative. Two-sided sign test p = 0.125. Wilcoxon signed-rank p = 0.031; that smaller p is the single negative gene (XRCC6, mean log2FC -0.037) having the smallest magnitude. These p-values are gene concordance on one collapsed rank, not replicate tests. They do not establish an NHEJ increase. They do show the panel is not the predicted decrease.

The only |log2FC| > 0.25 calls are T47D-limited and point up: LIG4 1.038 (MCF7 -0.046) and PRKDC 0.324 (MCF7 -0.249, just short of the −0.25 line). XRCC4, RIF1, TP53BP1, XRCC5, and XRCC6 stay FLAT in both lines.

PRKDC and TBK1 each have a second, lower-FPKM cufflinks fragment of the same locus. Summing fragments does not create a both-line NHEJ decrease: PRKDC sum-of-loci log2FC is MCF7 -0.273, T47D 0.308. Primary numbers stay on the higher-mean locus so they match the earlier rank.

c-NHEJ has 7 genes, under the prerank minimum of 8, so there is no NES for this panel.

## 2. STING pathway

| Gene | Mapped | MCF7 log2FC | T47D log2FC | Mean | MCF7 | T47D | Both lines | vs thesis |
|---|---|---:|---:|---:|---|---|---|---|
| CGAS | MB21D1 | 0.197 | 0.138 | 0.167 | FLAT | FLAT | FLAT | not increased |
| STING1 | absent | NA | NA | NA | ABSENT | ABSENT | ABSENT | unmeasured |
| TBK1 | TBK1 | 0.062 | 0.006 | 0.034 | FLAT | FLAT | FLAT | not increased |
| IRF3 | IRF3 | 0.116 | -0.103 | 0.006 | FLAT | FLAT | FLAT | not increased |
| STAT1 | STAT1 | -0.053 | 0.104 | 0.025 | FLAT | FLAT | FLAT | not increased |

CGAS is in the file as **MB21D1** (Entrez 115004): MCF7 0.197, T47D 0.138. STING1 and TMEM173 (Entrez 340061) are not in the table, so the receptor is unmeasured rather than a zero. TBK1, IRF3, and STAT1 are expressed (FPKM tens) and FLAT. Of the four measured genes, 4/4 have a positive both-line mean, sign-test p = 0.125, and every |log2FC| is ≤ 0.25. That is not STING up.

## 3. Hallmark IFN and APM

Prerank GSEA, weighted KS p = 1, 1000 gene-set permutations, seed 42, same engine as the CLDN4-loss prerank. Positive NES = the set sits at the CLDN4-KO end of the rank. BH-FDR is within these three sets for that contrast. Enrichment scores match the earlier GSE207704 rows (same rank and the same set order). The earlier FDR was BH across five headline sets, so the FDR column here is not the same number even though NES is.

| Set | Contrast | NES | nominal p | BH-FDR | Genes in rank | Mean log2FC | Call vs IFN/APM up |
|---|---|---:|---:|---:|---:|---:|---|
| Hallmark IFN-γ | T47D | -1.567 | 0.0080 | 0.0120 | 115 | -0.192 | decreased |
| Hallmark IFN-γ | MCF7 | -1.110 | 0.1109 | 0.3012 | 115 | -0.012 | not increased |
| Hallmark IFN-γ | mean of the two lines | -1.542 | 0.0070 | 0.0105 | 115 | -0.102 | decreased |
| Hallmark IFN-α | T47D | -1.829 | 0.0010 | 0.0030 | 60 | -0.339 | decreased |
| Hallmark IFN-α | MCF7 | -1.008 | 0.2008 | 0.3012 | 60 | 0.017 | not increased |
| Hallmark IFN-α | mean of the two lines | -1.637 | 0.0020 | 0.0060 | 60 | -0.161 | decreased |
| MHC-I / APM | T47D | -1.232 | 0.1439 | 0.1439 | 8 | -0.187 | not increased |
| MHC-I / APM | MCF7 | 0.868 | 0.3347 | 0.3347 | 8 | 0.062 | not increased |
| MHC-I / APM | mean of the two lines | -0.819 | 0.3357 | 0.3357 | 8 | -0.063 | not increased |

Coverage in this cufflinks annotation: Hallmark IFN-α 60/97, Hallmark IFN-γ 115/200, MHC-I/APM 8/21. Genes missing from the deposit cannot move the rank. The T47D IFN decrease is carried by genes that are present. MCF7 IFN-α mean log2FC is 0.017 while its NES is negative: the set is not shifted up. Leading edge, T47D IFN-α (KO-low end): TRIM5,UBE2L6,USP18,PSME2,IFIH1,NUB1,PARP12,HELZ2,IRF1,SLC25A28,HLA-C,OAS1,ISG20,TMEM140,IL4R,PARP9,TXNIP,TRIM21,LGALS3BP,IL7,IFI35,ISG15,IFI44,BST2. T47D IFN-γ: HLA-DRB1,DDX58,RAPGEF6,BPGM,FCGR1A,NMI,RNF213,PSMA2,NFKBIA,PML,ARID5B,UBE2L6,VAMP5,USP18,CASP7,PSME2,IFIH1,PARP12,HELZ2,SSPN,IRF5,IRF1,SLC25A28,ITGB7,ARL4A.

APM genes absent here include HLA-A, HLA-B, HLA-E, HLA-F, HLA-G, B2M, TAP1, TAP2, TAPBP, NLRC5, PSMB8, PSMB9, and ERAP2. Of the eight that are present, ERAP1 is down in both lines (MCF7 −0.729, T47D −0.386). IRF1 is down in T47D (−0.477) and flat in MCF7 (−0.154). HLA-C is up in MCF7 (+1.287) and down in T47D (−0.499). PSMB10 is up in MCF7 only (+0.357). No APM gene is up in both lines.

The “mean of the two lines” rank averages the two group-mean log2FC vectors. It is a summary, not a third cohort.

## Reading

Relative to KD → NHEJ down and STING/IFN up, GSE207704’s breast CLDN4 KO is not a supporting public example. NHEJ is flat to slightly higher (LIG4 and PRKDC up in T47D only). STING core is flat where it is measured. Hallmark interferon moves down in T47D. APM cannot be scored for classical MHC-I genes because they are absent, and the genes that remain are not up.

Numbers are descriptive directions on pooled FPKM. They are not a replicate-level DESeq2/edgeR result.

## Reproduce

```bash
python3 scripts/gse207704_nhej_sting/analyze.py
```

Inputs: `methods/gse207704_nhej_sting/data/GSE207704_CLDN4_RNAseq.txt.gz` (downloaded from GEO if missing) and `ifn_apm_sets.json`. Figure: `methods/gse207704_nhej_sting/figures/fig_nhej_sting_ifn.png`.
