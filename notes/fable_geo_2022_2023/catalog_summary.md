# GEO 2022-2023 human lung ICI series - verified catalog
- Candidate GSE series returned by the exhaustive search: **94**
- Passing text relevance (human + lung + ICI + expression): **42**
- With per-sample outcome annotation in GEO metadata: **30**
- Open, processed supplementary download <2 GB: **41**

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
| GSE186446 | 32 | 2023/02/14 | single_cell | nan |
| GSE248378 | 29 | 2023/11/26 | bulk_or_other | treatment |
| GSE190731 | 23 | 2022/02/23 | xenograft | treatment |
| GSE243238 | 20 | 2023/10/20 | bulk_tumor_rnaseq | clinical benefit_(cb);ici treatment |
| GSE235048 | 15 | 2023/08/03 | bulk_or_other | treatment |
| GSE228419 | 14 | 2023/04/02 | bulk_or_other | nan |
| GSE164146 | 12 | 2023/01/01 | single_cell | nan |
| GSE212622 | 12 | 2022/10/25 | bulk_or_other | infection status |
| GSE224099 | 12 | 2023/07/17 | single_cell | treatment |
| GSE193719 | 10 | 2022/12/31 | cell_line | nan |
| GSE189804 | 9 | 2023/03/22 | single_cell | treatment |
| GSE217451 | 9 | 2023/08/30 | cell_line | treatment |
| GSE224216 | 9 | 2023/08/30 | bulk_or_other | treatment |
| GSE150255 | 8 | 2023/05/01 | bulk_or_other | ifn-g treatment time |
| GSE194350 | 8 | 2022/11/03 | cell_line | nan |
| GSE195770 | 8 | 2023/01/01 | bulk_or_other | htll-280 status;pd-l1 status |
| GSE218402 | 8 | 2023/04/04 | bulk_or_other | treatment |
| GSE224246 | 8 | 2023/08/30 | bulk_or_other | treatment |
| GSE238006 | 8 | 2023/07/25 | bulk_or_other | treatment |
| GSE250254 | 8 | 2023/12/17 | single_cell | treatment |
| GSE193049 | 7 | 2022/04/07 | bulk_or_other | nan |
| GSE229353 | 7 | 2023/05/16 | single_cell | treatment |
| GSE178521 | 6 | 2022/02/08 | cell_line | treatment |
| GSE192591 | 6 | 2022/02/18 | single_cell | treatment |
| GSE192790 | 6 | 2022/01/02 | bulk_or_other | treatment |
| GSE197236 | 6 | 2022/02/27 | cell_line | nan |
| GSE198099 | 6 | 2023/12/31 | single_cell | nan |
| GSE213590 | 6 | 2023/09/20 | single_cell | treatment |
| GSE213902 | 6 | 2023/04/05 | single_cell | treatment |
| GSE223779 | 6 | 2023/08/09 | single_cell | treatment |
| GSE193707 | 4 | 2022/12/31 | bulk_or_other | response to neoadjuvant chemotherapy |
