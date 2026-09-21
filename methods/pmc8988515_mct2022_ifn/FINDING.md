# PMC8988515 IFN side

Yamamoto, Webb, Davis, Baumgartner, Woodruff, Guntupalli, Neville, Behbakht, Bitler. Loss of Claudin-4 Reduces DNA Damage Repair and Increases Sensitivity to PARP Inhibitors. Mol Cancer Ther 2022;21:647–657. PMC8988515. DOI 10.1158/1535-7163.MCT-21-0827.

The paper's data-availability line is "Data are available upon request from the corresponding author." The NHEJ / PARP-inhibitor experiments are left as published. This page scores the interferon side of the expression tables that could be downloaded.

## Repository hunt

| Place | Author or paper query | Hits |
|---|---|---|
| GEO `gds` | Bitler + CLDN4; Yamamoto + CLDN4 + ovarian; PMC8988515 / MCT-21-0827 | 0 |
| ArrayExpress | Bitler or Yamamoto + CLDN4, `collection:arrayexpress` | 0 |
| PRIDE | Bitler claudin; Yamamoto OVCAR3 CLDN4; PMC8988515 | 0 |
| NODE browse | CLDN4, claudin-4, Bitler, PMC8988515, MCT-21-0827, OVCAR3 CLDN4 | 0 |

NODE is answering queries: the control string `lung cancer` returns 4,747 records.

Cited public deposits, not new author submissions:

- Qian 2020 blueprint (ref 17): ArrayExpress **E-MTAB-8107**, **E-MTAB-6149**, **E-MTAB-6653**. BioStudies lists 48 processed count CSVs on E-MTAB-8107. `ftp.ebi.ac.uk` returned HTTP 502, TLS failure, and FTP 421 timeout from this environment, so those counts were not scored.
- Coscia 2016 (ref 22): PRIDE **PXD003668**. The archive files are RAW. The scored matrix is the Nature processed log2 MaxLFQ table (Supplementary Data 1, MOESM1558).
- Fig 1 transcriptome: TCGA ovarian serous cystadenocarcinoma, Firehose Legacy, median CLDN4 split. Author table is figshare file 39984958 (Table S1). Independent matrix is UCSC Xena `TCGA.OV.sampleMap/HiSeqV2` (log2 RSEM norm_count+1).

AACR supplements S2–S4 are BRCA status, drug viability, and ex vivo histology. They do not contain an RPPA matrix or an shCLDN4 RNA-seq matrix.

## Table S1 (their published rank)

19,834 genes. **1,582** with q < 0.05, the same count as the paper. CLDN4 log2(high/low) = **1.2**.

Prerank GSEA on that log2 ratio (weighted KS, p = 1, 1,000 gene-set permutations, seed 42). Positive NES means the set sits toward CLDN4-high. BH-FDR is inside these four sets.

| Set | NES | FDR | nom p | genes in rank | mean log2(high/low) |
|---|---:|---:|---:|---:|---:|
| Hallmark IFN-α | +2.089 | 0.0013 | 0.0010 | 97 | +0.066 |
| Hallmark IFN-γ | +1.575 | 0.0013 | 0.0010 | 200 | +0.029 |
| MHC-I / APM | +1.445 | 0.021 | 0.021 | 21 | +0.063 |
| KEGG tight junction | +1.535 | 0.0013 | 0.0010 | 167 | +0.023 |

Genes inside each set that themselves pass q < 0.05:

| Set | q < 0.05 | higher in CLDN4-high | higher in CLDN4-low | median log2 |
|---|---:|---:|---:|---:|
| IFN-γ | 7 / 200 | 4 | 3 | +0.04 |
| IFN-α | 1 / 97 | 0 | 1 | +0.08 |
| MHC-I / APM | 1 / 21 | 1 | 0 | +0.09 |
| Tight junction | 24 / 167 | 18 | 6 | −0.01 |

Tight junction after dropping CLDN4 from the rank: NES **+1.439**, nominal p = 0.001 (166 genes). The leading edge is other claudins (CLDN3, CLDN7, CLDN9, TJP3).

The NES is a coherent nudge of a large set. The per-gene log2 ratios inside IFN-γ are mostly a few hundredths, and almost none are q < 0.05.

## Xena TCGA-OV primary tumors

308 columns in HiSeqV2. Four recurrent tumors (`-02`) are excluded: TCGA-29-2414-02, TCGA-13-1489-02, TCGA-29-1770-02, TCGA-61-2008-02. **n = 304** primary. Median split drops ties: **152 vs 152**. The paper's cBioPortal freeze (20 Dec 2020) was 154 low and 153 high.

Welch-t prerank GSEA, same engine, positive = CLDN4-high:

| Set | NES | FDR | nom p | genes | mean t |
|---|---:|---:|---:|---:|---:|
| Hallmark IFN-α | +2.694 | 0.0013 | 0.0010 | 92 | +0.79 |
| Hallmark IFN-γ | +2.153 | 0.0013 | 0.0010 | 196 | +0.45 |
| MHC-I / APM | +2.026 | 0.0020 | 0.0020 | 21 | +0.94 |
| KEGG tight junction | +1.685 | 0.0013 | 0.0010 | 163 | +0.16 |

Mean gene-level Spearman ρ versus continuous CLDN4 (same 304 tumors): IFN-γ **+0.035**, IFN-α **+0.050**, MHC-I **+0.079**, tight junction **+0.017**.

Tumor-level mean z-score versus CLDN4 (n = 304):

| Set | ρ | p | partial ρ given KRT8/18/19 | partial p |
|---|---:|---:|---:|---:|
| IFN-γ | +0.067 | 0.25 | +0.060 | 0.30 |
| IFN-α | +0.085 | 0.14 | +0.068 | 0.24 |
| MHC-I / APM | +0.108 | 0.060 | +0.099 | 0.085 |
| Tight junction | +0.108 | 0.060 | +0.035 | 0.54 |

CLDN4 versus the KRT8/KRT18/KRT19 mean is ρ = **+0.266**, p = **2.6×10⁻⁶**, n = 304. Keratin adjustment leaves the IFN correlations small.

## Coscia LFQ (PXD003668 processed table)

CLDN4 protein is quantified in **12 of 30** lines, the same count the MCT paper cites. Eleven are ovarian-cancer lines (IGROV1, IGROV1CP, JHOS4, KURAMOCHI, OVCA433, OVCAR5, OVSAHO, PEO1, PEO4, SKOV3IP1, SNU119). The twelfth is ME180C13 (cervical). OVCAR3, the line used for the paper's shCLDN4 RPPA, has no CLDN4 value in this table.

Spearman of the set score versus CLDN4 protein:

| Cohort | Set | n | ρ | p |
|---|---|---:|---:|---:|
| 12 lines with CLDN4 | IFN-γ | 12 | +0.266 | 0.40 |
| 12 lines with CLDN4 | IFN-α | 12 | +0.490 | 0.11 |
| 12 lines with CLDN4 | MHC-I | 12 | +0.538 | 0.071 |
| 11 ovarian lines | IFN-γ | 11 | +0.364 | 0.27 |
| 11 ovarian lines | IFN-α | 11 | +0.482 | 0.13 |
| 11 ovarian lines | MHC-I | 11 | +0.500 | 0.12 |

On the same 12 lines, TP53BP1 ρ = **+0.50** (p = 0.10) and XRCC1 ρ = **−0.48** (p = 0.11). The signs match the paper's supplement (53BP1 positive, XRCC1 not). The p-values on this processed table are above 0.05.

## What this adds

In the TCGA ovarian split the paper published, IFN-α, IFN-γ, and MHC-I sit slightly toward CLDN4-high. The tumor-level IFN score does not track CLDN4 (IFN-γ ρ = +0.067, p = 0.25, n = 304). This is bulk HGSOC RNA, with CLDN4 also tracking KRT8/18/19. It is not an shCLDN4 transcriptome, and the Qian blueprint counts were not retrieved.
