# GEO 2022-2023 human lung ICI series - verified catalog
- Candidate GSE series returned by the exhaustive search: **94**
- Passing text relevance (human + lung + ICI + expression): **42**
- With per-sample outcome annotation in GEO metadata: **30**
- Open, processed supplementary download <2 GB: **41**
- Leftover / recovered series given a deep verdict this pass: **3**

## Deeply verified series
| GSE | n | category | lung | usable for TACSTD2/CLDN4 | note |
|---|---|---|---|---|---|
| GSE207422 | 39 | bulk_tumor_rnaseq | True | True | PRIMARY: NSCLC neoadjuvant anti-PD-1+chemo, 24 baseline bulk tumours, whole-transcriptome log2TPM, MPR/NMPR + RECIST + residual tumour %. Both genes present. |
| GSE162520 | 92 | targeted_panel | True | False | NSCLC PD-1/PD-L1 cohort with OS/PFS (n=92) but expression is a ~2.5k-gene targeted panel WITHOUT TACSTD2/CLDN4. |
| GSE161537 | 82 | targeted_panel | True | False | NSCLC immunotherapy cohort with RECIST/OS/PFS (n=82) but ~2.5k-gene targeted panel WITHOUT TACSTD2/CLDN4. |
| GSE216297 | 286 | platelet_rnaseq | True | False | NSCLC baseline platelet (TEP) RNA for nivolumab response (n=286); platelet not tumour, matrix in .RData, no per-sample response label in GEO metadata. |
| GSE243238 | 20 | bulk_tumor_rnaseq | False | True | CROSS-CHECK (NOT lung): acral melanoma, ICI-treated, bulk raw counts, clinical-benefit label. Both genes present. |
| GSE235500 | 34 | single_cell | True | False | Tumour-infiltrating Treg scRNA-seq (checkpoint blockade); sorted immune cells, not bulk tumour epithelium. |
| GSE235603 | 86 | single_cell | True | False | SuperSeries of Treg scRNA-seq; sorted immune cells. |
| GSE185206 | 108 | single_cell | True | False | scRNA/TCR-seq of T cells during checkpoint blockade. |
| GSE160903 | 268 | single_cell | True | False | scRNA of T cells (mouse+human); mechanistic. |
| GSE214992 | 50 | cell_line | True | False | NSCLC cell lines +/- CD8 T-cell killing; not patient outcome. |
| GSE190731 | 23 | xenograft | True | False | EGFR-mutant NSCLC xenografts, durvalumab/oleclumab; not patient outcome. |
| GSE248830 | 44 | targeted_panel | True | False | Brain-metastasis (breast+lung) NanoString-scale panel; not an ICI response cohort. |
| GSE221733 | 93 | spatial_dsp | True | True | RECOVERED leftover: NSCLC GeoMx DSP CTA, immunotherapy-treated, 41 pts with Responder/Non-responder. Missed by first-pass ICI regex (immunotherapy != immunotherap\b). TACSTD2 on panel; CLDN4 absent. Analysed patient-level PanCK+ TACSTD2. |
| GSE221322 | 96 | spatial_protein | True | False | Sister DSP protein panel of GSE221733; no TACSTD2/CLDN4 proteins. |
| GSE248378 | 29 | bulk_tumor_rnaseq | True | True | LEFTOVER PRIMARY: neoadjuvant durvalumab +/- SBRT, 29 post-Rx non-MPR tumours, FPKM with TACSTD2+CLDN4. Recurrence joined from Nat Commun source data (DurvaNNN). 45-M-PO arm-discordant. |
| GSE193049 | 7 | balf_rnaseq | True | True | LEFTOVER: BALF (not tumour) RNA, PD-1 responders vs non-responders n=7; both genes present. Analysed with that caveat. |
| GSE189045 | 23 | exosomal_mirna | True | False | Serum exosomal miRNA with anti-PD-1/PD-L1 response labels; no mRNA so TACSTD2/CLDN4 cannot be measured. |
| GSE250262 | 51 | nanostring_io360 | True | False | NSCLC tumour NanoString IO360 RCC; panel has EPCAM not TACSTD2/CLDN4; GEO metadata has no ICI response/survival label. |
| GSE185204 | 40 | single_cell | True | False | scRNA-seq of CD3+ / multimer+ T cells during ICB; not bulk tumour. |
| GSE186446 | 32 | bulk_or_other | True | False | Regional featureCounts from 3 ICB patients (T-cell lineage paper); both genes present but no per-sample ICI outcome. |
| GSE235048 | 15 | pbmc_rnaseq | True | False | PBMC TPM from NSCLC on various ICI regimens; TACSTD2 present, no response label (treatment only). |
| GSE228419 | 14 | sorted_tcell | True | False | Sorted CD8 T cells (PBMC/TIL); both genes present; no ICI outcome. |
| GSE164146 | 12 | single_cell | True | False | Treg/CD4 scRNA-seq in lung cancer; not bulk tumour epithelium. |
| GSE212622 | 12 | infection | False | False | False positive: murine trypanosome lung infection, not ICI. |
| GSE224099 | 12 | sorted_tcell | False | False | Melanoma (not lung) pre-ICI T cells; not a lung tumour cohort. |
| GSE193719 | 10 | cell_line | True | False | NSCLC cell lines +/- miR-455-5p; no patient ICI outcome. |
| GSE189804 | 9 | cell_line | True | False | A549 NRF2 KO RNA-seq; no patient ICI outcome. |
| GSE217451 | 9 | cell_line | True | False | H1650 hMENA siRNA; no patient ICI outcome. |
| GSE224216 | 9 | cell_line | True | False | H2030 hMENA siRNA RNA-seq; no patient ICI outcome. |
| GSE150255 | 8 | cell_line | True | False | A549/HCC827 IFN-γ time course; no patient ICI outcome. |
| GSE194350 | 8 | cell_line | True | False | A549 CRISPR/in-vivo MEN1 screen RNA-seq; no patient ICI outcome. |
| GSE195770 | 8 | ipf | False | False | False positive: idiopathic pulmonary fibrosis, not lung cancer ICI. |
| GSE218402 | 8 | cell_line | True | False | A549 + PBMC norepinephrine/adenosine; no patient ICI outcome. |
| GSE224246 | 8 | cell_line | True | False | H1650 ATAC-seq (hMENA); chromatin, not mRNA outcome cohort. |
| GSE238006 | 8 | cell_line | True | False | SCLC cell-line decitabine RNA-seq; no patient ICI outcome. |
| GSE250254 | 8 | sorted_tcell | True | False | Sorted CD4 memory T cells (COPD/NSCLC); no ICI outcome. |
| GSE229353 | 7 | single_cell | True | False | CD45+ scRNA-seq, chemo vs anti-PD-1+chemo; sorted immune cells. |
| GSE178521 | 6 | cell_line | True | False | NCI-H460 PRMT5 knockdown; no patient ICI outcome. |
| GSE192591 | 6 | sorted_tcell | True | False | Blood CD8 T cells +/- IL-27; no tumour ICI outcome. |
| GSE192790 | 6 | cell_line | True | False | A549 tumoursphere RNA-seq; no patient ICI outcome. |
| GSE197236 | 6 | cell_line | True | False | A549 radioresistance lines; not ICI. |
| GSE198099 | 6 | single_cell | True | False | NSCLC TME scRNA-seq, n=2 patients; no ICI outcome label. |
| GSE213590 | 6 | cell_line | True | False | PC-9 c-Jun overexpression; no patient ICI outcome. |
| GSE213902 | 6 | single_cell | True | False | PBMC CD3+ scRNA/TCR on chemo-IO; circulating T cells, not tumour TACSTD2/CLDN4. |
| GSE223779 | 6 | organoid | True | False | ALK+ tumour organoids; treatment field is sample origin, not ICI outcome. |
| GSE193707 | 4 | bulk_tumor_rnaseq | True | False | Mediastinal LN after neoadjuvant chemo (not ICI); n=4. |

## All relevant series (sorted by n)
| GSE | n | date | category | outcome keys |
|---|---|---|---|---|
| GSE216297 | 286 | 2023/07/10 | platelet_rnaseq | treatment |
| GSE160903 | 268 | 2022/01/11 | single_cell | treatment(1) |
| GSE185206 | 108 | 2023/02/13 | single_cell | nan |
| GSE162520 | 92 | 2022/07/07 | targeted_panel | os (event);os (month);pfs (event);pfs (month) |
| GSE235603 | 86 | 2023/09/11 | single_cell | nan |
| GSE161537 | 82 | 2022/07/07 | targeted_panel | best response on immunotherapy (recist);immunotherapy line;os (event);os (month) |
| GSE214992 | 50 | 2022/10/12 | cell_line | treatment |
| GSE248830 | 44 | 2023/12/03 | targeted_panel | smoking status;treatment after surgery of bm |
| GSE185204 | 40 | 2023/02/14 | single_cell | nan |
| GSE207422 | 39 | 2023/03/21 | bulk_tumor_rnaseq | pathologic_response;treatment group |
| GSE235500 | 34 | 2023/09/11 | single_cell | nan |
| GSE186446 | 32 | 2023/02/14 | bulk_or_other | nan |
| GSE248378 | 29 | 2023/11/26 | bulk_tumor_rnaseq | treatment |
| GSE190731 | 23 | 2022/02/23 | xenograft | treatment |
| GSE243238 | 20 | 2023/10/20 | bulk_tumor_rnaseq | clinical benefit_(cb);ici treatment |
| GSE235048 | 15 | 2023/08/03 | pbmc_rnaseq | treatment |
| GSE228419 | 14 | 2023/04/02 | sorted_tcell | nan |
| GSE164146 | 12 | 2023/01/01 | single_cell | nan |
| GSE212622 | 12 | 2022/10/25 | infection | infection status |
| GSE224099 | 12 | 2023/07/17 | sorted_tcell | treatment |
| GSE193719 | 10 | 2022/12/31 | cell_line | nan |
| GSE189804 | 9 | 2023/03/22 | cell_line | treatment |
| GSE217451 | 9 | 2023/08/30 | cell_line | treatment |
| GSE224216 | 9 | 2023/08/30 | cell_line | treatment |
| GSE150255 | 8 | 2023/05/01 | cell_line | ifn-g treatment time |
| GSE194350 | 8 | 2022/11/03 | cell_line | nan |
| GSE195770 | 8 | 2023/01/01 | ipf | htll-280 status;pd-l1 status |
| GSE218402 | 8 | 2023/04/04 | cell_line | treatment |
| GSE224246 | 8 | 2023/08/30 | cell_line | treatment |
| GSE238006 | 8 | 2023/07/25 | cell_line | treatment |
| GSE250254 | 8 | 2023/12/17 | sorted_tcell | treatment |
| GSE193049 | 7 | 2022/04/07 | balf_rnaseq | nan |
| GSE229353 | 7 | 2023/05/16 | single_cell | treatment |
| GSE178521 | 6 | 2022/02/08 | cell_line | treatment |
| GSE192591 | 6 | 2022/02/18 | sorted_tcell | treatment |
| GSE192790 | 6 | 2022/01/02 | cell_line | treatment |
| GSE197236 | 6 | 2022/02/27 | cell_line | nan |
| GSE198099 | 6 | 2023/12/31 | single_cell | nan |
| GSE213590 | 6 | 2023/09/20 | cell_line | treatment |
| GSE213902 | 6 | 2023/04/05 | single_cell | treatment |
| GSE223779 | 6 | 2023/08/09 | organoid | treatment |
| GSE193707 | 4 | 2022/12/31 | bulk_tumor_rnaseq | response to neoadjuvant chemotherapy |
