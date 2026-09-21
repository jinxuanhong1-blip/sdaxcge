# ArrayExpress / BioStudies CLDN4 knockdown–knockout wave

Search date: 2026-09-21. Live BioStudies, ArrayExpress-collection, PRIDE, and GEO queries. No accessions were invented.

## Verdict

No open CLDN4 knockdown, knockout, CRISPR, or siRNA matrix supports the four-axis pattern **NHEJ down, STING up, IFN up, and APM up**.

The only ArrayExpress transcriptomes that are actual CLDN4-loss experiments are **E-GEOD-50927** (mouse lung Cldn4 knockout) and **E-GEOD-22493** (SKOV-3 CLDN4 siRNA versus CLDN4 overexpression). The human CRISPR matrix **GSE207704** (MCF7 and T47D) is GEO-only. Its paper is BioStudies **S-EPMC10105442**, which does not hold the expression file and has no E-GEOD accession. It is scored so the ArrayExpress result is not mistaken for the cancer-cell test.

| Contrast | NHEJ (expect down) | STING (expect up) | IFN (expect up) | APM (expect up) |
|---|---|---|---|---|
| E-GEOD-50927 uninjured KO vs WT | null, median −0.011 | supports, median +0.129 | supports, median +0.459 | supports, median +0.368 |
| E-GEOD-50927 VILI-low KO vs WT | null | null | null | null |
| E-GEOD-50927 VILI-high KO vs WT | null | null | supports, median +0.245 | supports, median +0.261 |
| E-GEOD-22493 siRNA vs overexpression | null, median +0.359 | directional only, median +0.646, Wilcoxon p=0.109 | null, median −0.232 | null, median −0.212 |
| GSE207704 T47D CRISPR | null, median +0.071 | insufficient (3 genes) | opposite, median −0.404 | insufficient (5 genes, 4 down) |
| GSE207704 MCF7 CRISPR | null, median +0.108 | insufficient (3 genes) | opposite-directional, median −0.445 | insufficient (5 genes) |
| GSE207704 mean of the two lines | null, median +0.025 | insufficient | opposite, median −0.454 | insufficient |

`supports` means the set median log2FC is at least 0.10 in the expected direction and the one-sided Wilcoxon signed-rank p across genes is < 0.05. `directional` is the same magnitude with at least two-thirds of genes in that direction, without p < 0.05. These are gene-set tests. They are not sample-level replicate tests.

## What was found

ArrayExpress collection query `CLDN4` on 2026-09-21 returned four studies. Only E-GEOD-50927 and E-GEOD-22493 are CLDN4 genetic perturbations. E-GEOD-60885 is DNA methylation. E-GEOD-84742 is a colonic differentiation series, not a Cldn4 knockout. Broader ArrayExpress queries returned 0 E-MTAB, 0 E-MEXP, and 0 E-PROT records. BioStudies `S-BSST*` returned one CLDN4 study, S-BSST1967. PRIDE keyword search for CLDN4, claudin-4, and Cldn4 returned 0 projects. Full reject reasons are in `inventory.tsv`.

BioStudies HTTP file links returned 404, and `ftp.ebi.ac.uk` TLS failed from this environment. Processed matrices for the two ArrayExpress RNA experiments were therefore taken from the GEO deposits of the same accessions (GSE50927, GSE22493). GSE207704 was taken from its GEO FPKM table.

## E-GEOD-50927, uninjured lung

Cldn4 logFC −6.061, author FDR 4.07×10⁻²⁶. Design is 1 knockout lung versus 1 wild-type lung. The deposited edgeR table is used as published. VILIwtGenes is injury in wild-type lung (Cldn4 logFC +3.959) and is not a knockout contrast.

NHEJ core (11/12 genes; Paxx absent): median −0.011, 6 down / 5 up, Wilcoxon p=0.38, 0 genes FDR<0.05. Trp53bp1 logFC −0.012, FDR 1. NHEJ transcription does not fall.

STING set (9/9): median +0.129, 7 up / 1 down, Wilcoxon p=0.039, Mann–Whitney versus background p=0.029. The shift is the chemokine outputs: Ccl5 +2.710 (FDR 0.0011) and Cxcl10 +2.543 (FDR 0.094). The sensor/adapter subset (Cgas/Mb21d1, Tmem173, Tbk1, Ikbke, Irf3, Irf7) is only directional: median +0.123, 5 up / 1 down, Wilcoxon p=0.16, 0 genes FDR<0.05. Tmem173 +0.129 (FDR 1), Irf3 +0.045, Ifnb1 0. This is not a cGAS–STING sensor induction.

IFN (24 measured): median +0.459, 22 up / 2 down, Wilcoxon p=3.3×10⁻⁶, 7 genes FDR<0.05, including Isg15 +1.075 (FDR 0.0031), Oas2 +1.489 (FDR 0.0015), Oas3 +1.816 (FDR 0.0016), Ifit1 +0.576 (FDR 0.028).

APM (15/15): median +0.368, 13 up / 2 down, Wilcoxon p=1.5×10⁻⁴, 5 genes FDR<0.05 (B2m +0.655 FDR 0.034, Psmb9 +0.902 FDR 0.037, Calr +0.605 FDR 0.0014, H2-T23 +1.045 FDR 3.7×10⁻⁴, H2-M3 +1.191 FDR 8.0×10⁻⁶). Classical H2-K1 is not up (−0.071, FDR 1).

VILI-low is null on all four axes. VILI-high keeps a weaker IFN up (median +0.245, Wilcoxon p=0.0075, 2 genes FDR<0.05) and APM up (median +0.261, Wilcoxon p=0.0013, 0 genes FDR<0.05). NHEJ and STING stay null. Cldn4 remains strongly down (−12.69 and −10.16).

## E-GEOD-22493, ovarian siRNA

Three two-color arrays. Cy5 is CLDN4 siRNA and Cy3 is CLDN4 overexpression, so the value is log2(siRNA/overexpression), not log2(siRNA/scramble). CLDN4 itself is −1.225 (one-sample t p=0.25, 1 probe). Knockdown is directionally present and not statistically clean. TACSTD2 is not on GPL10555.

NHEJ median +0.359 (5 up / 4 down, Wilcoxon p=0.71). TP53BP1 −0.334 (p=0.52). STING median +0.646 is carried by noisy point estimates (CXCL10 +3.025, t p=0.54; IRF3 +1.065, p=0.45); signed-rank p=0.109 and no gene is FDR<0.05, so the call stays directional. IFN median −0.232. APM median −0.212, with TAP1/TAP2/B2M/TAPBP down. No axis meets the support rule.

## GSE207704, breast CRISPR

Collapsed FPKM, replicates already pooled, so there is no gene-level FDR. CLDN4 log2FC is −1.053 in T47D and −0.755 in MCF7 (pseudocount 0.1).

NHEJ does not fall (T47D median +0.071, MCF7 median +0.108). TP53BP1 is up in both lines (+0.210 T47D, +0.161 MCF7).

STING1/TMEM173, IFNB1, CXCL10, CCL5, IRF7, and IKBKE are absent from the deposit. The three measured genes (MB21D1, TBK1, IRF3) are too few to call.

IFN moves down. T47D median −0.404, 12 down / 6 up, Wilcoxon p=0.010, Mann–Whitney versus background p=0.0075. The two-line mean is −0.454 (Wilcoxon p=0.012). Large drops include T47D BST2 −4.68, IFI44 −3.39, ISG15 −1.94. MCF7 is the same direction (median −0.445, 13 down / 5 up) but the signed-rank p is 0.071 because a few genes go the other way (HERC5 +2.02, IFIH1 +1.89).

APM is mostly missing (HLA-A/B, B2M, TAP1/2, PSMB8/9, NLRC5 are absent). Of five measured genes, T47D has four down (HLA-C −0.521, ERAP1 −0.398, CALR −0.205, PSMB10 −0.043). That is not an APM increase, and it is below the six-gene floor for a set call.

## Not scored

- **S-BSST1967** (Xu et al., Cell Reports Medicine 2025, PMC12281411). BioStudies holds a DESeq2 RNA-seq table and a proteomics workbook. The paper’s bulk transcriptome is MHCC97H treated with the CLDN4-palmitoylation peptide CPP-S4 versus mock, and the proteomics are lenvatinib-versus-mock plus a palmitoyl-proteome. Raw sequence is CNGB CNP0006650. CLDN4 knockdown is in the paper and is not this deposited contrast. The ftp files were not opened because `ftp.ebi.ac.uk` TLS failed.
- **S-EPMC8988515** (Yamamoto et al., 2022). CLDN4 shRNA was used for NHEJ reporter assays. Table S1 is TCGA CLDN4-high versus low (1,582 transcripts), not a knockdown transcriptome. No GEO accession. The xlsx supplements were not opened; PMC returned a recaptcha page after the article HTML identified the table.
- **S-EPMC12603150** (2025). CLDN4 CRISPRi in OVCA429 and OVCAR3 with functional STING/IFN assays. The data-availability statement says raw data are available on request. No public matrix.

## Gene sets

NHEJ core: XRCC6, XRCC5, PRKDC, LIG4, XRCC4, NHEJ1, PAXX, DCLRE1C, POLL, POLM, PNKP, APTX. TP53BP1 is reported alone and is not inside the NHEJ call.

STING: CGAS/MB21D1, STING1/TMEM173, TBK1, IKBKE, IRF3, IRF7, IFNB1, CXCL10, CCL5. STING_sensor drops IFNB1, CXCL10, and CCL5.

IFN: ISG15, MX1, MX2, OAS1, OAS2, OAS3, OASL, IFIT1, IFIT2, IFIT3, IFI6, IFI27, IFI44, IFI44L, RSAD2, BST2, STAT1, STAT2, IRF1, IRF9, IFITM1, IFITM3, USP18, HERC5, IFIH1, DDX58, CXCL9, CXCL11, ISG20.

Mouse APM: H2-K1, H2-D1, B2m, Tap1, Tap2, Tapbp, Psmb8, Psmb9, Psmb10, Nlrc5, Erap1, Calr, Pdia3, H2-T23, H2-M3. Human APM uses HLA-A/B/C, B2M, TAP1/2, TAPBP, PSMB8/9/10, NLRC5, ERAP1/2, CALR, PDIA3, HLA-E.

Duplicate symbols keep the highest-logCPM row. Array probes collapse by median. GSE207704 duplicate symbols keep the highest total FPKM row.

## Files

- `inventory.tsv` — every ArrayExpress/BioStudies hit and why it was or was not scored
- `signature_scores.tsv` — one row per contrast and axis
- `de_genes.tsv` — gene-level log2FC for the four sets plus TP53BP1
- `qc_perturbation.tsv` — CLDN4 (and Cldn3/Cldn7) in each contrast
- `fig_signature_medians.png` — median log2FC for the four primary contrasts
- `key_stats.json`

Reproduce: `python3 scripts/ae_cldn4_kd_wave/analyze.py`
