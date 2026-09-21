# CLDN4 and interferon programs in public TROP2-ADC and TACSTD2-loss omics

Public matrices were searched for TROP2 antibody-drug conjugate treatment, TACSTD2 knockdown or knockout, and SKB264 / MK-2870 / sacituzumab tirumotecan. Accessible count matrices were recomputed. No value below is taken from a paper.

Full precision is in `results/contrasts.tsv`, `results/scores.tsv`, `results/correlations.tsv`, `results/shift_alignment.tsv`, and `results/context_genes.tsv`. The search log is `results/catalog.tsv`.

## Answer

One cell-line knockdown lowers CLDN4 and interferon together. That co-decrease does not repeat on cisplatin, in Trop2-knockout xenografts, or under sacituzumab govitecan, and the immunocompetent 4T1 tumors do not show Cldn4 tracking effector or interferon scores.

The only immunocompetent bulk tumor transcriptome is 4T1 Trop2 knockout in BALB/c (`GSE334497`). Tacstd2 falls (mean normalized count 159 to 11; difference of mean log2(count+1) = −3.82; exact Mann-Whitney p = 0.00794). Cldn4’s mean falls from 2605 to 1674 (delta −0.82) with p = 0.42, and the ten tumors overlap (knockout 464–3408, wild type 1148–5056). Hallmark interferon-gamma moves by +0.13 (187 genes, p = 0.31). A 10-gene effector score and a 9-gene leukocyte score are higher in knockout tumors (deltas +0.73 and +0.35, both p = 0.056). Across the same 10 tumors, Cldn4 does not correlate with those scores (Spearman rho = −0.18 effector, −0.15 leukocyte, −0.21 interferon-gamma; p = 0.63, 0.68, 0.56).

Elsewhere the signs disagree, or one side does not move:

| Dataset | Contrast | CLDN4 delta | IFN-gamma delta | Sample Spearman, CLDN4 vs IFN-gamma |
|---|---|---:|---:|---|
| GSE245459 | SKOV3 shTACSTD2 vs shNC | −1.92 (p = 0.10) | −0.52 (p = 0.10) | rho = 0.77, exact p = 0.10, n = 6 |
| GSE245459 | Same shRNA on cisplatin | −0.14 (p = 0.10) | +0.20 (p = 0.10) | rho = −0.60, exact p = 0.24, n = 6 |
| GSE334497 | 4T1 Trop2 KO vs WT, BALB/c | −0.82 (p = 0.42) | +0.13 (p = 0.31) | rho = −0.21, p = 0.56, n = 10 |
| GSE289287 | T-47D Trop2 KO vs WT, NRG | +0.29 (p = 0.11) | +0.10 (p = 0.40) | rho = 0.32, exact p = 0.50, n = 7 |
| GSE312098 | CX-1 IMMU132 vs control | −0.86 (p = 0.10) | +0.33 (p = 0.10) | rho = −0.89, exact p = 0.033, n = 6 |
| GSE304294 | KYSE30 IMMU132 vs control | +0.91 (p = 0.20) | +0.22 (p = 0.20) | rho = 0.90, exact p = 0.083, n = 5; treated n = 2 |
| GSE311016 | Five CRC PDX, paired IMMU132 | −0.44 (Wilcoxon p = 0.0625) | +0.029 (Wilcoxon p = 0.81) | paired-delta rho = −0.90, exact p = 0.083, n = 5 |

Delta is the difference of mean log2(expression+1), treatment minus control. For n = 3 versus 3 the smallest two-sided exact Mann-Whitney p is 0.10, so that p means the three treated values sit entirely on one side of the three controls. Scores with fewer than 5 genes above an expression floor of 1 are left out of this table. Cell-line interferon scores are tumor-cell transcript programs. In the cell-line matrices and in the CRC PDX matrix, 0 or 1 leukocyte gene cleared that floor.

## What was reanalyzed

**GSE245459, SKOV3 TACSTD2 shRNA.** Series text says “knockout”; sample names are shNC and sh. Without cisplatin, TACSTD2 FPKM falls from 5.40 to 0.034 (delta −2.63, p = 0.10) and CLDN4 from 3.15 to 0.10 (delta −1.92, p = 0.10). Replicates are completely separated: CLDN4 shNC 3.44, 2.98, 3.02 versus sh 0.060, 0.041, 0.20. Hallmark interferon-alpha (75 genes) delta −0.61 and interferon-gamma (125 genes) delta −0.52, both p = 0.10. CLDN4 versus interferon-alpha Spearman rho = 0.89, exact p = 0.033, n = 6. That correlation is the separation of the two groups of three. Effector genes stay below the floor. On cisplatin, shNC CLDN4 is already 0.46, 0.33, 0.30 and shDDP is 0.27, 0.27, 0.17, so the knockdown delta shrinks to −0.14 while interferon-gamma moves the other way (+0.20). TACSTD2 itself still falls (delta −4.03).

**GSE334497, 4T1 Trop2 knockout tumors in BALB/c.** Column IDs are the GEO library names: knockout KO162, KO164, KO165, KO172, RESUB-KO163R versus wild type RESUB-171R, RESUB-170R, RESUB-169R, RESUB-168R, control170. Mouse Hallmark lists were matched through Ensembl. Interferon-gamma uses 187/188 genes, interferon-alpha 93/94, effector 10/11 (`Gzmk` absent), leukocyte 9/10 (`Ms4a1` absent). Tacstd2 is separated (knockout 0.9–17, wild type 66–361). Cldn7 is lower in knockout tumors (mean 286 to 139, delta −0.82, p = 0.016) and is reported as context, not as a substitute for Cldn4. The effector genes are Cd8a, Cd8b1, Cd3d, Cd3e, Gzma, Gzmb, Prf1, Nkg7, Ifng, Ctsw. This is bulk RNA of whole tumors, so a higher effector score is a tissue-level transcript change.

**GSE289287, T-47D Trop2 knockout xenografts in NRG mice.** Recomputed from protein-coding DESeq2 normalized counts in the deposited table (wild type animals 2808, 2810, 2812; knockout 2807, 2815, 2817, 2818). TACSTD2 falls from 1723 to 159 (delta −3.47, p = 0.057). The deposited DESeq2 log2 fold change for TACSTD2 is −3.26 and for CLDN4 is +0.28, matching the recomputed signs and magnitudes (CLDN4 delta +0.29, mean 909 to 1111, p = 0.11). Interferon-gamma delta is +0.10 (171 genes, p = 0.40). NRG mice lack adaptive lymphocytes, and the matrix is human genes. In vitro Trop2-knockout RNA-seq is described in the series design and is not among the deposited samples. The same study’s DSG2-knockout xenografts, analyzed from their own normalized-count table, lower CLDN4 (909 to 469, delta −1.00, p = 0.017) and raise interferon-gamma (delta +0.67, 180 genes, p = 0.017) while TACSTD2 stays in place (delta −0.22, p = 0.12). A CLDN4 decrease with an interferon increase in this xenograft system occurs in the DSG2 comparison.

**GSE312098, CX-1 cells, 48 h, 3 µg/ml IMMU132.** Libraries X_1–X_3 are control and X_4–X_6 are IMMU132. CLDN4 FPKM 92 to 50 (delta −0.86, p = 0.10). Interferon-gamma delta +0.33 (127 genes, p = 0.10). Spearman rho = −0.89, exact p = 0.033, n = 6. TACSTD2 rises (delta +0.68). The single effector gene above the floor is NKG7. GSK2606414 alone moves CLDN4 by −0.055 (p = 0.70). The IMMU132 plus GSK arm is a combination and is in the tables as such.

**GSE304294, KYSE30 cells, 1 day IMMU132.** Libraries OX2_1 and OX2_2 versus OX1_1–OX1_3. CLDN4 FPKM 74 to 139 (delta +0.91) and interferon-gamma delta +0.22. With two treated samples the exact Mann-Whitney p cannot fall below 0.20. IACS-010759 alone moves interferon-gamma by −0.008 (p = 1).

**GSE311016, five CRC PDX models, day 29, paired control versus IMMU132.** Human FPKM columns C_* and T_* match the GEO treatment labels. CLDN4 paired log2 deltas are all negative: −0.087, −0.94, −0.56, −0.042, −0.55 (models 114, 36, 82, 83, 196). Mean delta −0.44. Wilcoxon exact p = 0.0625, which is the smallest two-sided value for five pairs. Unpaired Mann-Whitney p = 0.42 because the models differ in baseline. Interferon-gamma paired deltas are −0.11, +0.48, +0.15, −0.55, +0.17 (mean +0.029, Wilcoxon p = 0.81, 154 genes). The paired-delta Spearman of −0.90 (exact p = 0.083) ranks those small interferon changes against the CLDN4 changes. Leukocyte genes are below the floor. The host strain is not in the GEO fields read for this series. The matrix contains human genes.

## Catalog of what else exists

**SKB264 / MK-2870 / sacituzumab tirumotecan.** GEO queries for those names returned no expression series. The preclinical SKB264 paper (PMID 36620535) reports pharmacology and does not name an RNA-seq accession. No public SKB264 transcriptome was found.

**Sacituzumab govitecan single-cell studies from the colorectal TROP2 paper (Nature, 2026).** These are real ADC-treatment matrices and were not recomputed:

- `E-MTAB-16433`, subcutaneous PDOX, vehicle versus sacituzumab govitecan. Processed counts about 0.75 GB.
- `E-MTAB-16843`, organoid time course versus a non-targeting IgG1-SN-38 control.
- `E-MTAB-16849`, liver-metastasis time course. One processed file is about 23 GB.

BioStudies lists the files as public. `ftp.ebi.ac.uk` accepted a TCP connection and then closed the TLS handshake, and FTP port 21 timed out, so no counts were read. `E-MTAB-16836` is FOLFIRI, and `E-MTAB-16583` / `E-MTAB-16585` are baseline cohorts.

**Named in papers, and not treatment matrices.** `GSE235812` is 36 breast PDX samples labeled BCX P0/P1. It is the accession cited for baseline RNA-seq by PMID 37567892 and by the datopotamab deruxtecan PDX paper PMID 39585341. `GSE278664` is pretreatment high-grade serous ovarian RNA-seq from NCT02203513. `GSE309617` is carboplatin resistance in TNBC PDX; sacituzumab govitecan appears later as a drug-screen result. `GSE302284` is residual EGFR-inhibitor lung tumors. `GSE128724` is sorted mouse prostate Trop2-positive cells. `GSE78736` and `GSE69715` are Huh7.5 versus HepG2 and HCC tumor versus liver from a paper on TACSTD2, CLDN1, and occludin; the deposited arrays do not include a TACSTD2 knockdown.

**LINCS L1000 CRISPR consensus.** The Harmonizome TACSTD2 knockout consensus lists CLDN4 among genes with a positive standardized value. That list was not recomputed, and the standardized values are not fold changes.

## Methods

Expression delta = mean log2(expression+1) in the treated or knockout samples minus the same mean in controls. Group tests are two-sided exact Mann-Whitney tests on those log2 values, with the asymptotic test only if the exact test is unavailable. Paired PDX tests are two-sided exact Wilcoxon tests on the five paired log2 differences. A program score is the mean log2(expression+1) of genes in the set whose maximum in the contrast samples exceeds 1. Hallmark interferon-alpha and interferon-gamma gene lists were downloaded from MSigDB on 21 Sep 2026 (`resources/hallmark_interferon_*_{human,mouse}.txt`; Liberzon et al., Cell Systems 2015). The human effector list is CD8A, CD8B, CD3D, CD3E, GZMA, GZMB, GZMK, PRF1, NKG7, GNLY, IFNG, CTSW. The mouse list is the same symbols in mouse case, with Cd8b1 and without Gnly. The leukocyte list is PTPRC, CD3D, CD3E, CD4, CD8A, CD68, ITGAM, MS4A1, FCGR3A, NCR1, plus mouse Adgre1 in place of FCGR3A. Scores built from fewer than 5 genes are stored and flagged `interpretable_program = False`. Spearman p-values for n ≤ 8 are exact permutation tests. For n = 3 versus 3, p = 0.10 is the smallest exact Mann-Whitney p; for n = 5 versus 5, p = 0.00794 is the smallest. Duplicate gene symbols keep the row with the highest mean. GSE289287 uses protein-coding rows only. Mouse symbols were mapped with the Ensembl REST lookup (`results/mouse_symbol_to_ensembl.tsv`).
