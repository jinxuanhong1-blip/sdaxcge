# GSE135274 XRCC4 KO ± mirin — IFN / STING / cGAS / APM / NHEJ

Public counts only. This does not change the locked CosMx, concordant-4, GSE137244, TCGA keratin, or TISMO numbers.

## Result

XRCC4 knockout versus scrambled gRNA lowers Hallmark IFN. IFN-α NES -2.70 and IFN-γ NES -2.58 (nominal p=0.001 and 0.001, BH FDR=0.002 and 0.002 across 50 Hallmarks). Of 41 measured ISGs, 1 concordant UP and 28 concordant DOWN (median log2FC -0.79). STING/TMEM173 is DOWN (-0.65). cGAS/MB21D1 is FLAT (-0.23). APM: 0 UP / 7 DOWN of 19 (HLA-A, HLA-B, and B2M are DOWN). IFNB1 max CPM is 0.08 (LOW).

Mirin versus scrambled raises Hallmark IFN. IFN-α NES +2.46 and IFN-γ NES +2.35 (nominal p=0.001, BH FDR=0.005). 30 of 41 measured ISGs are concordant UP (median +0.64; ISG15 +0.86, IFIT1 +1.25, RSAD2 +2.81, STAT1 +0.76). cGAS/MB21D1 is FLAT (-0.29) and STING/TMEM173 is FLAT (+0.05). APM stays mostly flat (1 UP / 0 DOWN of 19). The mirin signal is the ISG cassette, not a rise in cGAS or STING mRNA.

XRCC4 KO + mirin versus scrambled stays IFN-low: IFN-α NES -2.07, IFN-γ NES -1.97 (FDR=0.005 and 0.005); ISGs 4 UP / 19 DOWN (median -0.31). Mirin on the KO background raises IFN relative to KO alone (IFN-α NES +2.10, IFN-γ NES +1.94, FDR=0.003; ISGs 17 UP / 3 DOWN) and does not put the double-block arm above scrambled.

NHEJ mRNA does not collapse and is not a large backup program. In the KO, Ku80/XRCC5 +0.40, DNA-PKcs/PRKDC +0.45, and XLF/NHEJ1 +0.35 are concordant UP; LIG4, Ku70, and alt-NHEJ are flat. Hallmark DNA repair is shifted toward the KO (see the GSEA table) with a small mean log2FC inside the set. PAXX is absent from this CLC gene table. IFNA1, IFNA2, IFNG, CXCL9, and HLA-G are absent too.

The same rank recovers the paper's non-IFN biology on KO+mirin versus scrambled: 8/8 listed glycolysis genes are concordant DOWN, CA9 -4.94, CDKN1A +1.24, ZMAT3 +1.24. Experiment 1 vs Experiment 2 log2FC Spearman on that contrast is 0.69. Hallmark hypoxia is negative in every contrast; Hallmark p53 is positive where mirin is in the numerator. That is the check that the IFN signs are not a processing artifact. The authors' Metascape term "defense response to virus" is a different gene set and is not re-tested here.

With 1000 permutations the smallest nominal p is 0.001. These are gene-set permutations of a two-experiment rank.

## Dataset

- **GSE135274**, PMID 35054780 (Benjamin et al., *Int J Mol Sci* 2022). HeLa, not lung.
- 2×2: scrambled gRNA vs XRCC4(−/−) clone 2G3, each ± mirin. Authors used **100 µM mirin**. Protein loss of XRCC4 is their Western blot; this table is mRNA.
- RNA at **48 h after TALEN + NHEJ-reporter transfection** in every arm. The contrast is NHEJ blockade on a TALEN-DSB background, not resting HeLa.
- Two biological experiments. GEO columns `*_1` and `*_2` are read 1 and read 2 quantified separately (technical). They were summed; CPM is unchanged if they are averaged instead.
- CLC mismatch bins `mis_0/1/2` were summed per gene (a partition of assigned reads).
- **Biological n = 2 per arm.** Gene-level Welch p-values are descriptive. Calls use sign concordance, not FDR.

## Method

- log2FC = difference of log2(CPM+1), paired inside Experiment 1 and Experiment 2, then averaged.
- **UP / DOWN**: both experiments |log2FC| > 0.25 and the same sign, and max CPM ≥ 1. Opposite signs beyond that cut = DISCORDANT. Otherwise FLAT. Max CPM < 1 = LOW.
- GSEA: MSigDB Hallmark 2023.2.Hs, preranked by mean paired log2FC, weighted KS (p=1), 1000 gene-set permutations, seed 42. Positive NES = enriched in the numerator (treated / KO) arm. BH FDR is across the 50 Hallmark sets inside that contrast.
- Genes in the rank: symbol with no spaces, and ≥10 counts in at least 2 of 8 libraries.

## QC

- Read1 vs read2 Spearman within each library: min 0.972, median 0.974.
- Genes after symbol filter: 21969. Genes in the GSEA rank: 14154.
- Genome-wide Spearman of Experiment 1 vs Experiment 2 log2FC (KO+mirin vs WT, the authors' strongest contrast): 0.692.

XRCC4 mRNA (mean CPM). The knockout is a small indel (authors: −2 bp / −10 bp) with **no protein**; mRNA is only partly lower.

| arm | Exp1 CPM | Exp2 CPM |
|---|---:|---:|
| WT scrambled | 5.79 | 5.78 |
| WT + mirin | 6.38 | 7.07 |
| XRCC4 KO | 5.39 | 2.63 |
| XRCC4 KO + mirin | 4.00 | 5.80 |

XRCC4 paired log2FC, KO vs WT: Exp1 -0.088, Exp2 -0.901, mean **-0.495** (FLAT). Exp1 is nearly unchanged, so the both-experiment rule does not call it DOWN.

Author-reported genes on KO+mirin vs scrambled (glycolysis genes were described as down; CA9, CDKN1A, ENO2, DUSP5, ZMAT3 were qPCR-checked):

| gene | Exp1 | Exp2 | mean log2FC | call |
|---|---:|---:|---:|---|
| CA9 | -5.015 | -4.870 | -4.943 | DOWN |
| CDKN1A | +1.199 | +1.282 | +1.241 | UP |
| ENO2 | -2.850 | -2.920 | -2.885 | DOWN |
| DUSP5 | -1.823 | -1.388 | -1.605 | DOWN |
| ZMAT3 | +1.178 | +1.295 | +1.237 | UP |
| PGK1 | -1.674 | -1.480 | -1.577 | DOWN |
| ALDOC | -2.263 | -1.711 | -1.987 | DOWN |
| PFKFB4 | -2.080 | -1.753 | -1.917 | DOWN |
| TPI1 | -1.508 | -1.159 | -1.334 | DOWN |
| ENO3 | -1.883 | -1.517 | -1.700 | DOWN |
| PFKP | -1.561 | -1.122 | -1.342 | DOWN |
| HK2 | -1.195 | -1.068 | -1.131 | DOWN |

## Hallmark GSEA

Positive NES means the set is higher in the first arm of the contrast.

| contrast | set | NES | nominal p | BH FDR (50 Hallmarks) | mean log2FC in set |
|---|---|---:|---:|---:|---:|
| XRCC4 KO vs scrambled | IFN-α | -2.705 | 0.001 | 0.002 | -0.484 |
| XRCC4 KO vs scrambled | IFN-γ | -2.578 | 0.001 | 0.002 | -0.350 |
| XRCC4 KO vs scrambled | p53 | -1.075 | 0.291 | 0.338 | -0.044 |
| XRCC4 KO vs scrambled | hypoxia | -2.298 | 0.001 | 0.002 | -0.263 |
| XRCC4 KO vs scrambled | DNA repair | +1.720 | 0.001 | 0.002 | +0.087 |
| mirin vs scrambled | IFN-α | +2.464 | 0.001 | 0.005 | +0.430 |
| mirin vs scrambled | IFN-γ | +2.348 | 0.001 | 0.005 | +0.304 |
| mirin vs scrambled | p53 | +1.656 | 0.001 | 0.005 | +0.180 |
| mirin vs scrambled | hypoxia | -2.858 | 0.001 | 0.005 | -0.152 |
| mirin vs scrambled | DNA repair | +0.975 | 0.497 | 0.591 | +0.115 |
| XRCC4 KO + mirin vs scrambled | IFN-α | -2.074 | 0.001 | 0.005 | -0.294 |
| XRCC4 KO + mirin vs scrambled | IFN-γ | -1.969 | 0.001 | 0.005 | -0.231 |
| XRCC4 KO + mirin vs scrambled | p53 | +1.469 | 0.001 | 0.005 | +0.053 |
| XRCC4 KO + mirin vs scrambled | hypoxia | -2.563 | 0.001 | 0.005 | -0.466 |
| XRCC4 KO + mirin vs scrambled | DNA repair | +1.356 | 0.003 | 0.012 | +0.074 |
| mirin on XRCC4 KO vs KO | IFN-α | +2.103 | 0.001 | 0.003 | +0.190 |
| mirin on XRCC4 KO vs KO | IFN-γ | +1.938 | 0.001 | 0.003 | +0.120 |
| mirin on XRCC4 KO vs KO | p53 | +1.836 | 0.001 | 0.003 | +0.097 |
| mirin on XRCC4 KO vs KO | hypoxia | -2.243 | 0.001 | 0.003 | -0.203 |
| mirin on XRCC4 KO vs KO | DNA repair | -0.718 | 0.573 | 0.578 | -0.013 |

Full 50-set tables: `tables/gsea_hallmark_all.tsv`.

## Panel calls

Measured = not ABSENT and not LOW. UP/DOWN are both-experiment calls.

| contrast | category | measured | UP | DOWN | FLAT | DISCORDANT | LOW | ABSENT | median log2FC |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| XRCC4 KO vs scrambled | qc | 2 | 0 | 0 | 2 | 0 | 0 | 0 | -0.202 |
| XRCC4 KO vs scrambled | cgas_sting | 14 | 0 | 2 | 12 | 0 | 0 | 0 | -0.065 |
| XRCC4 KO vs scrambled | ifn_ligand | 1 | 0 | 1 | 0 | 0 | 1 | 3 | -2.121 |
| XRCC4 KO vs scrambled | ifn_isg | 41 | 1 | 28 | 12 | 0 | 5 | 1 | -0.793 |
| XRCC4 KO vs scrambled | ifn_signaling | 13 | 0 | 4 | 9 | 0 | 0 | 0 | -0.141 |
| XRCC4 KO vs scrambled | apm | 19 | 0 | 7 | 11 | 1 | 1 | 1 | -0.260 |
| XRCC4 KO vs scrambled | nhej | 8 | 3 | 0 | 5 | 0 | 1 | 1 | +0.178 |
| XRCC4 KO vs scrambled | alt_nhej | 7 | 0 | 0 | 7 | 0 | 0 | 0 | +0.132 |
| mirin vs scrambled | qc | 2 | 1 | 0 | 1 | 0 | 0 | 0 | +0.277 |
| mirin vs scrambled | cgas_sting | 14 | 0 | 1 | 13 | 0 | 0 | 0 | +0.107 |
| mirin vs scrambled | ifn_ligand | 1 | 0 | 1 | 0 | 0 | 1 | 3 | -0.606 |
| mirin vs scrambled | ifn_isg | 41 | 30 | 1 | 10 | 0 | 5 | 1 | +0.637 |
| mirin vs scrambled | ifn_signaling | 13 | 4 | 1 | 8 | 0 | 0 | 0 | +0.239 |
| mirin vs scrambled | apm | 19 | 1 | 0 | 18 | 0 | 1 | 1 | +0.103 |
| mirin vs scrambled | nhej | 8 | 0 | 0 | 8 | 0 | 1 | 1 | +0.160 |
| mirin vs scrambled | alt_nhej | 7 | 2 | 0 | 5 | 0 | 0 | 0 | +0.207 |
| XRCC4 KO + mirin vs scrambled | qc | 2 | 0 | 0 | 2 | 0 | 0 | 0 | -0.038 |
| XRCC4 KO + mirin vs scrambled | cgas_sting | 14 | 1 | 4 | 9 | 0 | 0 | 0 | -0.137 |
| XRCC4 KO + mirin vs scrambled | ifn_ligand | 1 | 0 | 1 | 0 | 0 | 1 | 3 | -2.212 |
| XRCC4 KO + mirin vs scrambled | ifn_isg | 41 | 4 | 19 | 16 | 2 | 5 | 1 | -0.310 |
| XRCC4 KO + mirin vs scrambled | ifn_signaling | 13 | 1 | 3 | 8 | 1 | 0 | 0 | +0.005 |
| XRCC4 KO + mirin vs scrambled | apm | 19 | 0 | 10 | 9 | 0 | 1 | 1 | -0.377 |
| XRCC4 KO + mirin vs scrambled | nhej | 8 | 1 | 0 | 7 | 0 | 1 | 1 | +0.157 |
| XRCC4 KO + mirin vs scrambled | alt_nhej | 7 | 1 | 0 | 6 | 0 | 0 | 0 | +0.183 |
| mirin on XRCC4 KO vs KO | qc | 2 | 0 | 0 | 1 | 1 | 0 | 0 | +0.164 |
| mirin on XRCC4 KO vs KO | cgas_sting | 14 | 1 | 1 | 11 | 1 | 0 | 0 | -0.034 |
| mirin on XRCC4 KO vs KO | ifn_ligand | 1 | 0 | 0 | 1 | 0 | 1 | 3 | -0.090 |
| mirin on XRCC4 KO vs KO | ifn_isg | 41 | 17 | 3 | 19 | 2 | 5 | 1 | +0.320 |
| mirin on XRCC4 KO vs KO | ifn_signaling | 13 | 5 | 2 | 6 | 0 | 0 | 0 | +0.211 |
| mirin on XRCC4 KO vs KO | apm | 19 | 0 | 4 | 15 | 0 | 1 | 1 | -0.141 |
| mirin on XRCC4 KO vs KO | nhej | 8 | 1 | 0 | 7 | 0 | 1 | 1 | -0.009 |
| mirin on XRCC4 KO vs KO | alt_nhej | 7 | 0 | 0 | 7 | 0 | 0 | 0 | +0.054 |

Gene-level table: `tables/panel_de.tsv`.

### Genes called UP or DOWN in at least one contrast

| gene | category | KO vs WT | mirin vs WT | KO+mirin vs WT | mirin on KO |
|---|---|---|---|---|---|
| MRE11A | qc | FLAT (+0.090) | UP (+0.368) | FLAT (+0.143) | FLAT (+0.053) |
| cGAS/MB21D1 | cgas_sting | FLAT (-0.229) | FLAT (-0.290) | DOWN (-0.500) | FLAT (-0.271) |
| STING/TMEM173 | cgas_sting | DOWN (-0.654) | FLAT (+0.048) | DOWN (-0.719) | FLAT (-0.066) |
| IRF3 | cgas_sting | FLAT (-0.239) | FLAT (+0.173) | FLAT (+0.061) | UP (+0.300) |
| SAMHD1 | cgas_sting | FLAT (+0.139) | FLAT (+0.434) | UP (+0.430) | FLAT (+0.291) |
| ENPP1 | cgas_sting | DOWN (-0.755) | FLAT (-0.204) | DOWN (-1.144) | FLAT (-0.389) |
| RNASEH2C | cgas_sting | FLAT (+0.047) | DOWN (-0.681) | DOWN (-0.877) | DOWN (-0.924) |
| IL6 | ifn_ligand | DOWN (-2.121) | DOWN (-0.606) | DOWN (-2.212) | FLAT (-0.090) |
| ISG15 | ifn_isg | DOWN (-1.866) | UP (+0.861) | DOWN (-0.970) | UP (+0.896) |
| ISG20 | ifn_isg | DOWN (-1.139) | FLAT (+0.296) | DOWN (-0.824) | FLAT (+0.314) |
| MX1 | ifn_isg | UP (+0.405) | UP (+0.575) | UP (+0.506) | FLAT (+0.101) |
| MX2 | ifn_isg | DOWN (-0.598) | UP (+0.609) | DOWN (-0.493) | FLAT (+0.104) |
| OAS1 | ifn_isg | DOWN (-1.142) | UP (+0.871) | DOWN (-0.691) | FLAT (+0.452) |
| OAS2 | ifn_isg | DOWN (-1.249) | UP (+2.405) | DISCORDANT (+0.074) | UP (+1.323) |
| OAS3 | ifn_isg | DOWN (-0.672) | UP (+1.378) | FLAT (+0.141) | UP (+0.813) |
| OASL | ifn_isg | DOWN (-1.595) | UP (+0.522) | DOWN (-1.446) | FLAT (+0.149) |
| IFIT1 | ifn_isg | DOWN (-1.145) | UP (+1.250) | DOWN (-0.486) | UP (+0.659) |
| IFIT2 | ifn_isg | DOWN (-1.156) | UP (+0.907) | FLAT (-0.310) | UP (+0.845) |
| IFIT3 | ifn_isg | DOWN (-1.601) | UP (+1.034) | DOWN (-0.724) | UP (+0.877) |
| IFIT5 | ifn_isg | DOWN (-0.793) | UP (+1.265) | DOWN (-0.483) | FLAT (+0.310) |
| RSAD2 | ifn_isg | DOWN (-0.932) | UP (+2.815) | UP (+0.690) | UP (+1.623) |
| USP18 | ifn_isg | DOWN (-0.468) | UP (+0.613) | FLAT (-0.148) | FLAT (+0.320) |
| IFI6 | ifn_isg | DOWN (-1.379) | UP (+1.516) | FLAT (-0.015) | UP (+1.364) |
| IFI16 | ifn_isg | FLAT (+0.146) | UP (+0.545) | UP (+0.456) | FLAT (+0.310) |
| IFI27 | ifn_isg | DOWN (-1.525) | UP (+1.435) | FLAT (-0.364) | UP (+1.161) |
| IFI35 | ifn_isg | DOWN (-0.600) | FLAT (+0.275) | DOWN (-1.050) | DOWN (-0.451) |
| IFI44 | ifn_isg | DOWN (-1.524) | UP (+0.582) | DOWN (-1.048) | UP (+0.476) |
| IFI44L | ifn_isg | DOWN (-0.970) | UP (+2.639) | FLAT (-0.110) | UP (+0.860) |
| IFITM1 | ifn_isg | DOWN (-1.917) | UP (+0.403) | DOWN (-1.852) | FLAT (+0.066) |
| IFITM2 | ifn_isg | FLAT (+0.322) | FLAT (-0.437) | FLAT (-0.195) | DOWN (-0.516) |
| IFITM3 | ifn_isg | FLAT (-0.022) | FLAT (-0.208) | DOWN (-0.597) | FLAT (-0.575) |
| DDX58 | ifn_isg | DOWN (-1.440) | FLAT (+0.382) | DOWN (-1.131) | FLAT (+0.309) |
| IFIH1 | ifn_isg | DOWN (-0.583) | UP (+0.849) | FLAT (-0.022) | UP (+0.561) |
| DDX60 | ifn_isg | DOWN (-1.719) | FLAT (-0.070) | DOWN (-1.532) | FLAT (+0.186) |
| HERC5 | ifn_isg | DOWN (-1.324) | UP (+0.582) | DOWN (-0.978) | FLAT (+0.346) |
| XAF1 | ifn_isg | DOWN (-0.860) | UP (+1.573) | FLAT (-0.178) | UP (+0.682) |
| EIF2AK2 | ifn_isg | FLAT (-0.115) | UP (+0.699) | FLAT (+0.226) | UP (+0.342) |
| SAMD9L | ifn_isg | DOWN (-1.024) | UP (+1.057) | DOWN (-0.810) | FLAT (+0.214) |
| TRIM22 | ifn_isg | FLAT (-0.199) | UP (+1.475) | UP (+0.655) | UP (+0.855) |
| PLSCR1 | ifn_isg | FLAT (-0.334) | UP (+1.215) | FLAT (+0.215) | UP (+0.548) |
| STAT1 | ifn_signaling | FLAT (-0.114) | UP (+0.763) | FLAT (+0.167) | FLAT (+0.281) |
| STAT2 | ifn_signaling | DOWN (-0.369) | FLAT (+0.384) | FLAT (+0.005) | UP (+0.374) |
| IRF1 | ifn_signaling | DOWN (-0.841) | FLAT (+0.239) | DOWN (-0.375) | UP (+0.466) |
| IRF7 | ifn_signaling | DOWN (-0.447) | UP (+0.664) | FLAT (+0.326) | UP (+0.774) |
| IRF9 | ifn_signaling | DOWN (-1.040) | UP (+0.875) | FLAT (-0.267) | UP (+0.773) |
| JAK1 | ifn_signaling | FLAT (+0.096) | UP (+0.470) | UP (+0.505) | UP (+0.409) |
| IFNAR2 | ifn_signaling | FLAT (-0.141) | DOWN (-0.363) | DOWN (-0.577) | DOWN (-0.436) |
| SOCS3 | ifn_signaling | FLAT (-0.219) | FLAT (-0.192) | DOWN (-0.726) | DOWN (-0.507) |
| GBP1 | ifn_isg | DOWN (-1.378) | UP (+0.637) | DOWN (-1.580) | DISCORDANT (-0.201) |
| GBP2 | ifn_isg | DOWN (-0.548) | DOWN (-0.303) | DOWN (-0.994) | DOWN (-0.446) |
| GBP4 | ifn_isg | FLAT (-0.598) | UP (+0.723) | DISCORDANT (-0.122) | UP (+0.476) |
| GBP5 | ifn_isg | FLAT (+0.265) | UP (+0.494) | FLAT (+0.338) | FLAT (+0.074) |
| CCL5 | ifn_isg | FLAT (-0.211) | UP (+0.711) | FLAT (+0.263) | FLAT (+0.475) |
| IDO1 | ifn_isg | DOWN (-0.599) | FLAT (+0.469) | DOWN (-0.536) | DISCORDANT (+0.063) |
| HLA-A | apm | DOWN (-0.640) | FLAT (+0.103) | DOWN (-0.581) | FLAT (+0.059) |
| HLA-B | apm | DOWN (-0.854) | FLAT (-0.165) | DOWN (-1.094) | FLAT (-0.240) |
| HLA-C | apm | DOWN (-0.568) | FLAT (+0.261) | FLAT (-0.377) | FLAT (+0.192) |
| HLA-E | apm | DOWN (-0.528) | FLAT (+0.191) | FLAT (-0.307) | FLAT (+0.221) |
| B2M | apm | DOWN (-0.445) | FLAT (+0.001) | DOWN (-0.699) | FLAT (-0.254) |
| NLRC5 | apm | DISCORDANT (+0.008) | FLAT (+0.207) | FLAT (-0.300) | DOWN (-0.308) |
| TAP1 | apm | FLAT (-0.402) | FLAT (+0.217) | DOWN (-0.539) | FLAT (-0.137) |
| PSMB8 | apm | FLAT (-0.101) | FLAT (-0.045) | DOWN (-0.457) | DOWN (-0.355) |
| PSMB9 | apm | FLAT (-0.427) | FLAT (+0.227) | DOWN (-0.489) | FLAT (-0.062) |
| PSMB10 | apm | DOWN (-0.450) | FLAT (-0.113) | DOWN (-1.029) | DOWN (-0.579) |
| ERAP1 | apm | FLAT (+0.056) | FLAT (+0.157) | DOWN (-0.377) | DOWN (-0.433) |
| ERAP2 | apm | FLAT (-0.260) | UP (+0.539) | FLAT (-0.025) | FLAT (+0.234) |
| CANX | apm | FLAT (-0.164) | FLAT (-0.257) | DOWN (-0.469) | FLAT (-0.305) |
| PSME1 | apm | DOWN (-0.465) | FLAT (-0.037) | DOWN (-0.358) | FLAT (+0.107) |
| Ku80/XRCC5 | nhej | UP (+0.403) | FLAT (+0.101) | FLAT (+0.209) | FLAT (-0.193) |
| Ku70/XRCC6 | nhej | FLAT (+0.241) | FLAT (+0.109) | UP (+0.352) | FLAT (+0.111) |
| DNA-PKcs/PRKDC | nhej | UP (+0.454) | FLAT (+0.131) | FLAT (+0.268) | FLAT (-0.186) |
| XLF/NHEJ1 | nhej | UP (+0.345) | FLAT (+0.189) | FLAT (+0.283) | FLAT (-0.062) |
| Artemis/DCLRE1C | nhej | FLAT (-0.247) | FLAT (+0.247) | FLAT (+0.104) | UP (+0.351) |
| RAD50 | alt_nhej | FLAT (-0.095) | UP (+0.310) | FLAT (-0.041) | FLAT (+0.054) |
| POLB | alt_nhej | FLAT (+0.293) | UP (+0.528) | UP (+0.386) | FLAT (+0.093) |

## Limits

- n=2. A Hallmark nominal p is a permutation of the gene rank, not a patient-level or even a well-powered sample test.
- HeLa + transfected TALEN breaks. Not a lung tumor, not ICI, not cGAS-STING stimulation as the experiment's intent.
- XRCC4 protein is gone by the authors' blot; XRCC4 mRNA is a weak QC.
- Mirin blocks MRE11 exonuclease activity. MRE11A mRNA is not expected to fall.
- This file does not re-litigate CLDN4 exclusion, TISMO, or the concordant-4 cohort.
