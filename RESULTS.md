# Signed NHEJ set from Yamamoto 2022, scored in CLDN4-loss transcriptomes

Paper: Yamamoto, Webb, et al., *Loss of Claudin-4 Reduces DNA Damage Repair and Increases Sensitivity to PARP Inhibitors*, Molecular Cancer Therapeutics 2022;21:647–657. PMID 35373300, PMC8988515, DOI 10.1158/1535-7163.MCT-21-0827.

The paper reports that CLDN4 knockdown lowers NHEJ activity in an I-SceI GFP reporter and lowers 53BP1 foci, and that RPPA protein for 53BP1 and XRCC1 falls in OVCAR3 shCLDN4. This note asks whether a signed NHEJ gene set taken from the paper’s own tables moves in the same direction in two public CLDN4-loss mRNA datasets.

## Deposited omics

The manuscript data-availability statement says data are available on request from the corresponding author. Entrez searches done for this analysis returned 0 records in GEO (`gds`), SRA, and BioProject for `Yamamoto AND Bitler AND CLDN4`, for PMID 35373300, and for the manuscript id MCT-21-0827. Europe PMC marks the article as having supplementary files; the text-mined accession types are Cellosaurus and DOI. There is no GEO or SRA accession for the OVCAR3 shCLDN4 RPPA or for a knockdown transcriptome.

PMC hosts the author-manuscript supplements (`NIHMS1777898-supplement-1` through `-9`):

| File | What it is |
|---|---|
| supplement-1.docx | Legends for supplementary figures 1–4 |
| supplement-2 to -5.png | Supplementary figures (CanSAR/PrognoScan, cell-line blots, mass-spec correlations, ex vivo Ki67) |
| supplement-6.xlsx | **Table S1.** TCGA ovarian Firehose Legacy, CLDN4-high vs CLDN4-low. 19,834 genes. 1,582 with q < 0.05 (614 higher in CLDN4-high, 968 higher in CLDN4-low) |
| supplement-7.xlsx | Table S2. Cell-line BRCA1/2 status and claudin-4 protein level |
| supplement-8.xlsx | Table S3. OVCAR3 shCLDN4 drug-screen viabilities and DepMap/GDSC2 correlations |
| supplement-9.xlsx | Ex vivo tumor table (GTFB id, histology, claudin-4, BRCA, HRD) |

Table S1 is a reanalysis of public TCGA tumors. It is the gene list used below. Tables S2–S4 do not list knockdown genes. Figure 3A is an RPPA heatmap (P < 0.05, FDR < 15%). The supplement does not include that protein matrix. The text names two DNA-repair proteins as reduced after CLDN4 shRNA: 53BP1 (`TP53BP1`) and XRCC1.

The extracted Table S1 is `nhej_cldn4/input/yamamoto2022_table_s1.tsv`.

## Signed NHEJ gene set

Membership is classical NHEJ plus the 53BP1 axis, written in `nhej_cldn4/analysis.py`:

- KEGG 2021 Human “Non-homologous end-joining”: `DCLRE1C`, `DNTT`, `FEN1`, `LIG4`, `MRE11`, `NHEJ1`, `POLL`, `POLM`, `PRKDC`, `RAD50`, `XRCC4`, `XRCC5`, `XRCC6`
- Reactome R-HSA-5693571 factors on the 53BP1/DSB-response side of NHEJ: `TP53BP1`, `RIF1`, `RNF8`, `RNF168`, `H2AX`, `MDC1`, `PAXIP1`, `ATM`, `NBN`, `TDP1`, `TDP2`
- Shieldin: `MAD2L2`, `SHLD1`, `SHLD2`, `SHLD3`

Histone peptides and the BRCA1-A complex are left out of this list. They sit inside the broad Reactome NHEJ reaction record and are a different repair choice. `XRCC1` is kept beside the set because the paper reports it, and it is scored on its own. XRCC1 is single-strand break repair, so it is not given an NHEJ sign.

A gene enters the signed set when it is on that NHEJ list and Table S1 gives q < 0.05. The sign is the Table S1 column “Higher expression in”: +1 if higher in CLDN4-high, −1 if higher in CLDN4-low. Seven genes pass. All seven are higher in CLDN4-low tumors (sign −1).

| Gene | Membership | log2(CLDN4 high / low) | q | Sign with CLDN4 |
|---|---|---:|---:|---:|
| ATM | Reactome NHEJ DSB response | −0.29 | 0.0060 | −1 |
| MDC1 | Reactome NHEJ DSB response | −0.26 | 0.0185 | −1 |
| MRE11 | KEGG NHEJ | −0.22 | 0.0287 | −1 |
| PRKDC | KEGG NHEJ | −0.27 | 0.0044 | −1 |
| RNF168 | Reactome NHEJ DSB response | −0.26 | 0.0217 | −1 |
| SHLD2 | shieldin | −0.23 | 0.0298 | −1 |
| TDP1 | Reactome NHEJ | −0.24 | 0.0052 | −1 |

`TP53BP1` is on the NHEJ list and is the protein the paper highlights. In Table S1 its log2 ratio is −0.15 and q = 0.13, so it stays out of the q < 0.05 set. `XRCC1` log2 ratio is +0.04, q = 0.76. Both are reported below as the RPPA pair, with the paper’s protein direction (down after CLDN4 loss).

The paper also names `POLR2J`, `BRIP1`, and `NUP160` as examples from the Reactome cell-cycle / DNA-repair / DNA-replication enrichments of the 1,582 genes. Those three are in Table S1 at q < 0.05 (`POLR2J` +0.29, q = 0.0045, higher in CLDN4-high; `BRIP1` −0.50, q = 0.0031; `NUP160` −0.24, q = 0.0012). They are not NHEJ genes and are not in the score.

`STING1` (`TMEM173` in Table S1) has log2 ratio +0.37, q = 0.053, higher in CLDN4-high. It sits just outside the paper’s q < 0.05 call and is not an NHEJ gene. Neither knockdown matrix contains `STING1` or `TMEM173`, so it has no KD log2 fold change here.

Full membership with Table S1 statistics: `nhej_cldn4/output/nhej_signed_geneset.tsv`.

## How the public contrasts were scored

**GSE207704** (Murakami et al., Breast Cancer Research 2023). CRISPR CLDN4 knockout in T47D and MCF7, two replicates each. The deposited file `GSE207704_CLDN4_RNAseq.txt.gz` has one FPKM per genotype; the replicates are already collapsed. FPKM was summed across cufflinks rows that share a gene symbol. log2 fold change is log2((KO + 0.1) / (WT + 0.1)). A gene is scored when the higher of the two FPKM values is at least 1.

**GSE22493** (SKOV-3-IP-Luc, three two-color arrays). GEO channel 2 (Cy5) is CLDN4 knockdown and channel 1 (Cy3) is the CLDN4-overexpressing control. The series-matrix value matches the ScanArray normalized Ch2 log ratio on GSM558700 (Pearson r = 0.9998 on 30,477 paired probes), so a negative value is lower in the knockdown channel. `CLDN4` itself is −1.74 and −0.71 on the two arrays that report it. Probe symbols come from the Operon V3 name prefix. Within a replicate, probes for one gene are summarized by the median; the gene log2 ratio is the mean of replicates that have a value, and a gene needs two such replicates. Old symbols were mapped before scoring: `MRE11A`→`MRE11`, `NBS1`→`NBN`, `G22P1`→`XRCC6`, `TTRAP`→`TDP2`, `FAM35A`→`SHLD2`, `H2AFX`→`H2AX`, `TMEM173`→`STING1`.

**Signed concordance.** For each measured set gene, contribution = sign_with_CLDN4 × (−log2FC). The set score is the mean of those contributions. Because every gene in this set has sign −1, the observed score equals the mean KD log2FC. A positive score means these transcripts are higher after CLDN4 loss, which is the direction Table S1 records for them.

The null is 10,000 draws of an equal-sized gene set from the other Table S1 q < 0.05 genes measured in that contrast (seed 20220401; gene lists sorted; two-sided p = (n_as_or_more_extreme + 1) / 10001). A second null ignores sign and tests the mean log2FC. A third, separate readout is the mean log2FC of every detected NHEJ-list gene, not only the seven with q < 0.05, permuted against all genes that pass the expression filter. That third number is the direct mRNA reading of “NHEJ transcripts down.”

`MDC1` has no row in the GSE207704 FPKM table. `RNF168` has no Operon V3 probe, so it is absent from GSE22493. Each contrast therefore scores 6 of the 7 genes.

## Results

CLDN4 moves down in every contrast, which is the expected direction for these knockdowns and knockouts.

| Contrast | CLDN4 log2FC | Signed-set n | Concordant genes | Mean signed concordance | Permutation p | Set genes with log2FC < 0 | Mean log2FC of detected NHEJ list | NHEJ-list permutation p |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| T47D CLDN4−/− | −1.073 | 6 | 4/6 | +0.070 | 0.772 | 2/6 | +0.018 (25 genes, 44% down) | 0.873 |
| MCF7 CLDN4−/− | −0.726 | 6 | 3/6 | +0.027 | 0.872 | 3/6 | +0.020 (25 genes, 48% down) | 0.856 |
| Mean of T47D and MCF7 | −0.900 | 6 | 3/6 | +0.048 | 0.747 | 3/6 | +0.019 (25 genes, 48% down) | 0.781 |
| SKOV-3 CLDN4 knockdown | −1.225 | 6 | 4/6 | +0.170 | 0.636 | 2/6 | −0.009 (21 genes, 52% down) | 0.968 |

The mean of T47D and MCF7 is the average of the two breast log2FC values. It is not a third experiment. Background sizes for the signed-set null are 1,102 (T47D), 1,125 (MCF7), 1,091 (mean), and 843 (SKOV-3) Table S1 q < 0.05 genes.

Per-gene KD log2FC for the signed set:

| Gene | T47D | MCF7 | SKOV-3 |
|---|---:|---:|---:|
| ATM | +0.261 | +0.509 | +0.611 |
| MDC1 | absent | absent | +1.212 |
| MRE11 | −0.041 | +0.105 | +0.047 |
| PRKDC | +0.311 | −0.276 | +0.359 |
| RNF168 | +0.029 | −0.067 | absent |
| SHLD2 | −0.174 | +0.130 | −0.762 |
| TDP1 | +0.032 | −0.237 | −0.449 |

On the SKOV-3 arrays, one-sample t-tests across replicates are descriptive (n = 2 or 3). `MDC1` mean +1.21, p = 0.13. `SHLD2` mean −0.76, p = 0.65 (replicates −2.00 and +0.48). `TDP1` mean −0.45, p = 0.65. `ATM` mean +0.61, p = 0.36. `PRKDC` mean +0.36, p = 0.31. `MRE11` mean +0.05, p = 0.94. `CLDN4` is negative on both reporting arrays (−1.74, −0.71; mean −1.23; t-test p = 0.25 with n = 2).

### 53BP1 and XRCC1 mRNA

The paper’s OVCAR3 RPPA result is a protein decrease for both. In these mRNA contrasts:

| Gene | T47D log2FC | MCF7 log2FC | SKOV-3 log2FC (replicates) |
|---|---:|---:|---|
| TP53BP1 | +0.210 | +0.161 | −0.333 ( +0.444, −0.415, −1.029 ); t-test p = 0.52 |
| XRCC1 | +0.113 | −0.182 | +0.687 ( +0.687, +0.651, +0.723 ); t-test p = 9.1×10⁻⁴ |

SKOV-3 `XRCC1` mRNA is higher on all three arrays. T47D and MCF7 `TP53BP1` mRNA are higher. SKOV-3 `TP53BP1` mRNA is lower, with replicates on both sides of zero.

## Reading the scores

Table S1 puts the seven q < 0.05 NHEJ genes on the CLDN4-low side of ovarian tumors. In T47D, MCF7, and SKOV-3, the mean log2FC of the genes that could be measured is positive and small (+0.03 to +0.17). That is the Table S1 direction. The same genes are not shifted relative to the other 1,582-gene calls (permutation p 0.64–0.87).

The wider detected NHEJ list (25 genes in each breast line, 21 in SKOV-3) has a mean log2FC of +0.018, +0.020, and −0.009. About half the genes are lower after CLDN4 loss. Those means sit inside the permutation null (p 0.78–0.97).

The OVCAR3 measurements in the paper are protein (RPPA) and a functional NHEJ reporter. The two public datasets are mRNA after CLDN4 loss in T47D, MCF7, and SKOV-3. CLDN4 itself is lower in all three. The signed NHEJ mRNA set and the 53BP1/XRCC1 transcripts do not add a coordinated decrease on top of that.

Numbers in this file match `nhej_cldn4/output/set_scores.tsv`, `gene_scores.tsv`, and `gse22493_gene_ttests.tsv` from `python3 nhej_cldn4/analysis.py`.
