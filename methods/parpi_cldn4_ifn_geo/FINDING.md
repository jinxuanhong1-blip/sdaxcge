# PARP inhibitor RNA-seq: DNA repair to IFN, scored on CLDN4 and tight junctions

Public only. No private 8-KL matrices. Locked lung results (CosMx exclusion, concordant-4 ρ, GSE137244 KL>KP, TCGA keratin correction, TISMO Tacstd2 49/64) are unchanged.

GEO has **no series indexed as CLDN4 and a PARP inhibitor**. The CLDN4–NHEJ paper deposited no RNA-seq. The bridge is therefore two published mechanisms plus PARPi RNA-seq matrices scored here for CLDN4, a core tight-junction set, and a core IFN set in the same samples.

## Literature bridge

**CLDN4 to DNA repair.** Yamamoto et al., Mol Cancer Ther 2022 (PMID 35373300). In BRCA-wild-type HGSOC cells, CLDN4 knockdown lowers 53BP1 and XRCC1 protein, cuts NHEJ reporter activity about twofold, and leaves homology-directed repair statistically unchanged. OVCAR3 and OVCA429 become more sensitive to olaparib and rucaparib. In primary ovarian tumor slices, CLDN4-low tumors showed a stronger drop in Ki67 after olaparib than CLDN4-high tumors. TCGA ovarian tumors split at median CLDN4 were enriched for DNA repair, replication, and cell-cycle genes. DepMap: higher CLDN4 mRNA tracked with higher AUC (less sensitivity) for olaparib, niraparib, rucaparib, and cisplatin. Data availability in the paper is “upon request.”

**DNA repair deficit to IFN.** Three papers, two without a GEO series we could recompute:

- Ding et al., Cell Reports 2018 (PMID 30540933). Olaparib in Brca1-deficient mouse ovarian tumors induces STING-dependent type I IFN signaling and CD8 responses. Rescored below as GSE120500.
- Pantelidou et al., Cancer Discovery 2019 (PMID 31015319). Olaparib recruits CD8 T cells through tumor-cell cGAS/STING in BRCA-deficient triple-negative breast cancer, more strongly than in HR-proficient models. No GEO record for this PMID.
- Sen et al., Cancer Discovery 2019 (PMID 30777870). DDR inhibition in small-cell lung cancer activates STING and T cells. No GEO record for this PMID. This is SCLC, not LUAD and not KL.

## How the matrices were scored

Pre-specified genes only. Counts were converted to log2(CPM+1) using the full library size. FPKM and TPM used log2(x+1). GSE120792 was left on its deposited scale (ACTB mean 9.94, the log2 band). Family log2FC is the mean of gene-level log2FC. Exact Mann–Whitney p is reported for n≥3 per group. For n=4 versus 4 the smallest two-sided exact p is 0.029, and for n=6 versus 6 it is 0.0022; those values mean the group scores do not overlap. For n=3 the smallest exact p is 0.10, so those contrasts are reported as direction and effect size. n<3 has no p.

Core IFN is 28 ISGs and chemokines (ISG15, IFIT1–3, MX1/2, OAS1–3, OASL, CXCL9/10/11, CCL5, CD274, STAT1/2, and related sensors). Core TJ excludes CLDN4 (CLDN1/3/7/18, OCLN, TJP1–3, F11R, CGN, JAM2/3, and related barrier genes). The broad thesis IFN and TJ lists from the concordant-4 work are in the tables; they dilute the ISG signal because they contain many non-ISG genes. Repair is TP53BP1, XRCC1, PRKDC, LIG4, XRCC4, NHEJ1, PARP1, BRCA1/2, RAD51, RAD51C, BRIP1.

## Acute PARP inhibition

| Contrast | core IFN | CLDN4 | core TJ | recruit |
|---|---:|---:|---:|---:|
| UWB1.289 BRCA1-mutant, olaparib 3.5 µM, 96 h, n=4 (GSE237361) | +2.14 (26/26 up), p=0.029 | −0.005, p=0.89 | +0.09, p=0.49 | +1.86 (5/5), p=0.029 |
| OVCAR3 BRCA1-WT, olaparib 7.5 µM, 96 h, n=4 (GSE237361) | +0.40 (21/26), p=0.029 | −0.43, p=0.029 | −0.02, p=0.69 | +0.37, p=0.029 |
| UWB1.289, 24 h, n=3 (GSE243208) | +0.33 (22/26) | +0.13 | −0.01 | +0.16 |
| Brca1-deficient bulk tumors, olaparib, n=6 (GSE120500) | +2.32 (25/25), p=0.0022 | absent from the panel | +0.48 (4 genes on panel), p=0.13 | +3.09, p=0.0022 |
| OVCAR3 talazoparib, n=3 (GSE285827) | +0.40 (22/27) | −0.53 | +0.25 | +0.29 |
| CAOV3 talazoparib, n=3 (GSE285827) | +0.76 (25/27) | −0.33 | −0.17 | +0.31 |
| OVCAR3 veliparib, n=3 (GSE285827) | +0.08 | −0.22 | +0.07 | +0.08 |
| SBC5 SCLC, olaparib 1 µM, 3 d, n=3 (GSE233820) | +1.33 (22/27) | +0.74, p=0.40 | median +0.35 | +2.53 (5/5) |

Gene-level anchors, same GSE237361 96 h contrast, all with non-overlapping groups (p=0.029) unless noted:

- UWB1.289: CXCL10 +3.81, IFIT1 +3.74, MX1 +3.72, OAS1 +2.50, ISG15 +2.31, CCL5 +2.10, CD274 +0.71, STAT1 +1.42. TP53BP1 −0.05 (p=0.69). XRCC1 −0.29 (p=0.20). IFNB1 is undetectable in both arms.
- OVCAR3: CXCL10 +0.85, ISG15 +0.83, MX1 +0.78, STAT1 +0.26. CLDN4 moves from log2(CPM+1) 7.80 to 7.37. F11R +0.32 and TJP1 +0.35, so the flat core-TJ mean is a mix of small moves.

Both lines already express CLDN4 at 96 h baseline (UWB1.289 DMSO 7.47, OVCAR3 DMSO 7.80). This RNA-seq does not show a large CLDN4 gap between them.

GSE120500 is an Ion AmpliSeq panel of 4,604 genes, not a full transcriptome. Cldn4 is absent. Cxcl10 +3.33, Ccl5 +3.22, Ifit1 +3.13, Mx1 +2.99, Ifng +1.73, Cd274 +1.33, all p=0.0022. F11r +1.32 (p=0.041). Columns were assigned to vehicle versus olaparib by matching library sums to SRA spot counts (within 2% of that group and more than 3% from the other group). Control-3/4/5 are too close to name individual replicates; the six-versus-six split is the contrast. The tumor is bulk, so IFN includes infiltrating cells. That matches the Ding paper’s in vivo claim.

OVCAR3 received the higher concentration (7.5 µM versus 3.5 µM) and still shows the smaller IFN shift. Talazoparib, the stronger trapper, moves IFN in OVCAR3 and CAOV3; veliparib in the same series stays near zero (CAOV3 veliparib core IFN −0.15). The 24 h UWB1.289 profile is a smaller IFN shift than 96 h.

SBC5 (SCLC, not LUAD) raises CCL5 +3.04 and CD274 +2.73 with olaparib alone. Adding 4 Gy reverses the direction (core IFN −0.51, CLDN4 −0.86, n=3).

Repair mRNAs stay nearly flat under acute olaparib (UWB1.289 family −0.03, p=0.49; tumor panel +0.13, p=0.18). OVCAR3 repair transcripts edge up (+0.15, p=0.029). Yamamoto’s NHEJ result is protein and foci, and these RNA-seq matrices do not show a transcriptional shutdown of 53BP1 or XRCC1.

## Acquired resistance is a different state

A2780 olaparib-resistant versus untreated parental cells (GSE153867, FPKM, n=8). CLDN4 log2(FPKM+1) goes from 0.55 to 6.72 (log2FC +6.17, p=1.6×10⁻⁴). TACSTD2 goes from 0.014 to 8.24 (log2FC +8.22). CLDN7 +4.96. The parental line is near the floor for these genes, so this is an off-to-on epithelial gain in the resistant line. Core IFN mean is only +0.17, and the genes split: ISG15 −2.79, MX1 +1.50, STAT1 +0.91. CXCL10 is flat (−0.004). This is not the acute BRCA1-mutant IFN program.

Drug-free UWB1.289 olaparib-resistant cells versus parental cells (GSE235980, n=2, descriptive): BRCA1-deficient resistant cells, core IFN −1.02 (author DESeq2: MX1 −2.56, IFIT1 −2.24, CXCL10 −2.42). BRCA1-restored resistant cells, core IFN +0.84 (author MX1 +3.10, CLDN4 −0.45, padj=0.0057). n=2, no test.

## Genotype without drug

GSE120792, BRCA1-mutant versus BRCA1-wild-type ovarian cells, n=3, no olaparib column. Core IFN +0.46 (22/27 up). CLDN4 −0.19 (p=0.40). ISG15 −0.85. The genotype shift is smaller than olaparib added on top of BRCA1 loss in GSE237361.

## What this adds next to the CLDN4 thesis

CLDN4 protein can hold NHEJ together in BRCA-wild-type ovarian cells (Yamamoto). PARP inhibition, when homologous recombination is already broken, raises the core IFN set by log2FC +2.14 in UWB1.289 cells and +2.32 in Brca1-deficient tumors (about 4-fold and 5-fold; GSE237361 and GSE120500). CXCL10 moves further (+3.81 and +3.33). In that BRCA1-mutant cell contrast, CLDN4 mRNA stays put. The IFN-low side of the locked lung CLDN4 result has a public mechanistic neighbor: repair competence, including the NHEJ activity CLDN4 supports, sits on the quiet side of the STING–IFN axis, and PARP inhibition is a documented way to open the IFN side.

These series are ovarian and SCLC models. They leave the locked CosMx exclusion ratios, the concordant-4 patient correlation, and the KL-versus-KP cell-line result as they are. HeLa n=1 (GSE298546, core IFN +0.76, CLDN4 +0.02) and sarcoma n=2 (GSE239639, core IFN within ±0.12) are descriptive only.

## Not scored

- GSE295677 and GSE293216 (PMID 41572432): LIN28B plus BMN673, malignant ascites and vascular permeability. The expression table columns are not keyed to sample titles.
- GSE247622: CDC7 plus PARPi, sample codes not decoded.
- GSE117765: PEO1 resistant clones versus adherent cells.
- GSE165548: pancreatic radiation codes not decoded.

## Files

- `scripts/score_parpi_cldn4_ifn.py` — downloads GEO supplements and writes the tables
- `results/tables/family_contrasts.tsv`, `gene_contrasts.tsv`, `family_log2fc_wide.tsv`
- `results/tables/literature_catalog.tsv`
- `results/figures/parpi_ifn_cldn4_tj.png`
