# FINDING — GSE207704 CLDN4 knockout: NHEJ, STING, IFN/APM scores vs WT

Additive public evidence. This file does not replace the CosMx, concordant-4, GSE137244, TCGA keratin, or TISMO results. It is not a lung dataset and it is not SKB264.

## Design

GSE207704 (PMID 37059993) is CRISPR **CLDN4−/− vs parental WT** in **T47D and MCF7** breast cancer lines. GEO lists 2 biological replicates per genotype (GSM6310640–GSM6310647). The only processed matrix, `GSE207704_CLDN4_RNAseq.txt.gz`, is **group-collapsed cufflinks FPKM** (one column per genotype per line). Sample-level *t* tests and gene-level FDR are not available from this file.

The parental WT arm is described as genotype WT, “No treatment”. That is the negative-control comparator used here. The series has **no siRNA non-targeting arm**. Calling the experiment “KD vs NC” overstates the design: it is a knockout versus parental WT.

Rank and scores use log2((FPKM + 0.5) / (FPKM_WT + 0.5)) after keeping the highest-mean-FPKM locus per symbol. Positive = higher after CLDN4 loss. Per-gene UP/DOWN uses |log2FC| > 0.25 and max FPKM ≥ 1. The gene score is the unweighted mean of log2(FPKM + 0.5) on genes that clear that expression floor. Δ = score(KO) − score(WT).

## QC

| Line | CLDN4 log2FC (KO vs WT) | TACSTD2 log2FC |
|---|---:|---:|
| T47D | -1.039 | -0.623 |
| MCF7 | -0.750 | -1.008 |

CLDN4 mRNA falls in both lines but is not abolished (residual FPKM remains). TACSTD2 falls with it. Symbols in the rank: 12635.

## NHEJ core (PRKDC, LIG4, XRCC4, XRCC5, XRCC6, NHEJ1)

All six genes are in the deposit. A consensus of UP means one line clears |log2FC| > 0.25 and the other line is not DOWN. No core NHEJ gene is UP in both lines.

| Gene | T47D log2FC | MCF7 log2FC | mean | consensus |
|---|---|---|---:|---|
| PRKDC | 0.324 UP | -0.249 FLAT | 0.038 | UP |
| LIG4 | 1.038 UP | -0.046 FLAT | 0.496 | UP |
| XRCC4 | -0.038 FLAT | 0.238 FLAT | 0.100 | FLAT |
| XRCC5 | 0.164 FLAT | 0.168 FLAT | 0.166 | FLAT |
| XRCC6 | -0.181 FLAT | 0.108 FLAT | -0.037 | FLAT |
| NHEJ1 | 0.083 FLAT | -0.112 FLAT | -0.014 | FLAT |

NHEJ_CORE: near_zero. T47D Δ 0.232 (n=6/6 expressed, up 2, down 0, flat 4); MCF7 Δ 0.018 (n=6/6 expressed, up 0, down 0, flat 6)

The 6-gene panel is under the 8-gene GSEA floor. Its NES is reported with a nominal permutation p and **no BH-FDR**.

## STING axis (CGAS, STING1, TBK1, IRF3)

CGAS is deposited as **MB21D1** (GRCh38 v90 symbol). STING1 (TMEM173, Entrez 340061) is **not in the file**, so the 4-gene axis cannot be scored complete. The three measured genes sit near zero.

| Gene | T47D log2FC | MCF7 log2FC | mean | consensus |
|---|---|---|---:|---|
| CGAS (MB21D1) | 0.138 FLAT | 0.197 FLAT | 0.167 | FLAT |
| STING1 | absent | absent | absent | ABSENT |
| TBK1 | 0.006 FLAT | 0.062 FLAT | 0.034 | FLAT |
| IRF3 | -0.103 FLAT | 0.116 FLAT | 0.006 | FLAT |

STING_AXIS: near_zero. T47D Δ 0.014 (n=3/4 expressed, up 0, down 0, flat 3); MCF7 Δ 0.125 (n=3/4 expressed, up 0, down 0, flat 3)

## IFN and APM gene scores

Same gene lists as the earlier GSE207704 IFN/MHC-I panel. HLA-A is in the APM score, not the IFN score. Genes missing from the cufflinks table are excluded from the mean, not entered as zero. Many classical APM genes (HLA-A/B, B2M, TAP1/2, PSMB8/9, NLRC5) are absent, so the APM score is only the genes the file actually contains.

IFN: both_lines_down. T47D Δ -0.469 (n=28/60 expressed, up 2, down 12, flat 14); MCF7 Δ -0.167 (n=25/60 expressed, up 5, down 12, flat 8)

APM: near_zero. T47D Δ -0.155 (n=9/22 expressed, up 0, down 3, flat 6); MCF7 Δ 0.000 (n=9/22 expressed, up 2, down 2, flat 5)

## GSEA

Engine: weighted KS *p* = 1, 1000 gene-set permutations, seed 42. Positive NES = enriched among genes that rise after CLDN4 loss. BH-FDR is computed **inside the pathway sets** (Reactome plus the IFN and APM panels) for each contrast. Named NHEJ and STING panels are not in that FDR family.

### Pathway sets (size floor 8)

| Contrast | Set | NES | nominal p | FDR | n in rank | mean log2FC |
|---|---|---:|---:|---:|---:|---:|
| T47D | APM | -1.192 | 0.149 | 0.347 | 9 | -0.155 |
| T47D | IFN | -1.842 | 0.002 | 0.010 | 29 | -0.451 |
| T47D | REACTOME_CLASS_I_MHC_PEPTIDE_LOADING | 0.721 | 0.341 | 0.397 | 19 | 0.036 |
| T47D | REACTOME_INTERFERON_ALPHA_BETA_SIGNALING | -1.794 | 0.003 | 0.010 | 38 | -0.294 |
| T47D | REACTOME_INTERFERON_GAMMA_SIGNALING | -1.000 | 0.282 | 0.397 | 47 | -0.096 |
| T47D | REACTOME_NONHOMOLOGOUS_END_JOINING | 0.773 | 0.323 | 0.397 | 29 | 0.052 |
| T47D | REACTOME_STING_MEDIATED_INDUCTION | -0.729 | 0.448 | 0.448 | 8 | -0.007 |
| MCF7 | APM | 0.784 | 0.400 | 0.526 | 9 | 0.000 |
| MCF7 | IFN | -1.524 | 0.015 | 0.052 | 29 | -0.235 |
| MCF7 | REACTOME_CLASS_I_MHC_PEPTIDE_LOADING | 1.126 | 0.150 | 0.262 | 19 | 0.127 |
| MCF7 | REACTOME_INTERFERON_ALPHA_BETA_SIGNALING | -1.540 | 0.005 | 0.035 | 38 | -0.190 |
| MCF7 | REACTOME_INTERFERON_GAMMA_SIGNALING | -1.248 | 0.063 | 0.147 | 47 | -0.074 |
| MCF7 | REACTOME_NONHOMOLOGOUS_END_JOINING | 0.592 | 0.526 | 0.526 | 29 | 0.043 |
| MCF7 | REACTOME_STING_MEDIATED_INDUCTION | 0.598 | 0.485 | 0.526 | 8 | 0.091 |
| mean_both_lines | APM | -0.874 | 0.324 | 0.402 | 9 | -0.077 |
| mean_both_lines | IFN | -1.893 | 0.002 | 0.007 | 29 | -0.343 |
| mean_both_lines | REACTOME_CLASS_I_MHC_PEPTIDE_LOADING | 1.244 | 0.083 | 0.145 | 19 | 0.082 |
| mean_both_lines | REACTOME_INTERFERON_ALPHA_BETA_SIGNALING | -1.865 | 0.001 | 0.007 | 38 | -0.242 |
| mean_both_lines | REACTOME_INTERFERON_GAMMA_SIGNALING | -1.264 | 0.072 | 0.145 | 47 | -0.085 |
| mean_both_lines | REACTOME_NONHOMOLOGOUS_END_JOINING | 0.800 | 0.345 | 0.402 | 29 | 0.048 |
| mean_both_lines | REACTOME_STING_MEDIATED_INDUCTION | 0.668 | 0.423 | 0.423 | 8 | 0.042 |


### Named panels (size floor 3, FDR not applied)

| Contrast | Set | NES | nominal p | FDR | n in rank | mean log2FC |
|---|---|---:|---:|---:|---:|---:|
| T47D | NHEJ_CORE | 1.253 | 0.108 | NA | 6 | 0.232 |
| T47D | STING_AXIS | 0.589 | 0.424 | NA | 3 | 0.014 |
| MCF7 | NHEJ_CORE | -0.451 | 0.445 | NA | 6 | 0.018 |
| MCF7 | STING_AXIS | 0.891 | 0.326 | NA | 3 | 0.125 |
| mean_both_lines | NHEJ_CORE | 1.151 | 0.165 | NA | 6 | 0.125 |
| mean_both_lines | STING_AXIS | 0.955 | 0.270 | NA | 3 | 0.069 |


Reactome STING (R-HSA-1834941) also contains PRKDC, XRCC5, and XRCC6, so that pathway is not a pure cGAS–STING1–TBK1–IRF3 test. The STING_AXIS score above is the four-gene readout.

## What this accession supports

- **NHEJ core is not a shared program after CLDN4 loss.** T47D Δ = +0.232 is carried by LIG4 (+1.038) and PRKDC (+0.324); the other four T47D genes are FLAT. MCF7 Δ = +0.018 and every MCF7 core gene is FLAT (PRKDC −0.249 sits just under the 0.25 cutoff). XRCC5 is slightly higher in both lines (+0.164, +0.168) and still FLAT. Named-panel NES is +1.253 (nominal p 0.108) in T47D and −0.451 (nominal p 0.445) in MCF7. Reactome NHEJ (29 genes in rank) is not enriched (FDR 0.397 and 0.526).
- **STING axis does not move, and STING1 is missing.** CGAS is MB21D1 in this table and is FLAT (+0.138 T47D, +0.197 MCF7). TBK1 and IRF3 are FLAT. STING1/TMEM173 is absent, so 3 of 4 genes are scored (T47D Δ +0.014, MCF7 Δ +0.125). The 3-gene NES is nominal only. Reactome STING has 8 of 16 genes in the rank, mixes in PRKDC/XRCC5/XRCC6, and is not enriched (FDR 0.448 and 0.526).
- **IFN gene score is lower in both lines.** T47D Δ −0.469 (12 down, 2 up among 28 expressed genes). MCF7 Δ −0.167 (12 down, 5 up among 25). Genes down in both lines include OAS3, USP18, IFI6, and IFI44. HERC5 is the IFN gene up in both lines (+1.108 T47D, +1.924 MCF7). GSEA on the mapped IFN panel: T47D NES −1.842 FDR 0.010; MCF7 NES −1.524 FDR 0.052; mean of lines NES −1.893 FDR 0.007. Reactome interferon alpha/beta signaling is down in both lines (T47D NES −1.794 FDR 0.010; MCF7 NES −1.540 FDR 0.035). The GSEA mean log2FC is a little more negative than the expressed-gene score because low-FPKM genes stay in the rank.
- **APM is not opened.** Nine of 22 requested genes are present. HLA-A/B, B2M, TAP1/2, PSMB8/9, and NLRC5 are absent. Among genes that are present, ERAP1 is down in both lines, HLA-C is discordant (T47D down, MCF7 up), PSMB10 is up only in MCF7, and CALR/CANX/PDIA3/SEC61A1 are flat. APM GSEA FDR is 0.347 (T47D) and 0.526 (MCF7).

Permutation p-values test where the gene set sits on this collapsed rank. They are not tests of the n=2 biological replicates. Replicate-level uncertainty is not estimable from the deposited file.
