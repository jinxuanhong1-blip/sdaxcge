# Public mouse lung GEMM catalog for Cldn4

Public inventory only. No private KL matrices. GEO, SRA, and ArrayExpress were searched on 21 Sep 2026 for mouse lung tumors whose genotype is Kras (K), Kras/Trp53 (KP), Kras/Lkb1 (KL), Kras/Trp53/Lkb1 (KPL), or Stk11/Lkb1, with whole-transcriptome scRNA or bulk so Cldn4 is in the feature space. 10x is listed first.

GEO esearch (`Mus musculus` + lung + Kras or Lkb1/Stk11 + expression profiling) returned 325 series. Sample titles, characteristics, and supplementary filenames were read from GEO SOFT. Where a series looked closed (DEG-only name, or a RAW.tar with no count file beside it), the archive listing, the file header, the paper data-availability statement, figshare supplements, and the linked SRA BioProject were checked.

**How to read the columns**

- **n mice** is biological mice when the GEO sample title names them. A library that pools mice, or a pair of GEX+HTO rows, is not counted as one mouse per GSM.
- **T/NK** is `yes, mixed` when the digest is whole tumor or total viable cells. Bulk whole tumor is `yes, genes` (T/NK transcripts are in the sample; they are not a cell fraction). Epithelial, AT2, or CD45− sorts are `no`. CD45+ sorts are `T/NK only` (Cldn4 is the wrong compartment).
- **Open matrix** means a processed count, TPM, FPKM, or normalized matrix is on GEO or in the RAW.tar. SRA FASTQ alone is not an open matrix.
- Whole-transcriptome 10x and bulk RNA-seq include Cldn4 by assay design. Rows marked **Cldn4 checked** are files opened in this pass.

ChIP, ATAC, barcode-seq, methylation, human series, EML4-ALK (GSE176185), TC1/MCT2 (GSE338088), and Lewis lung without a Kras/Lkb1 genotype are not in the tables. E-MTAB-13704 is a lung GEMM aPD-L1 RNA-seq with open counts, but the study record does not label Kras, Lkb1, or Stk11, so it is not a K/KP/KL set.

---

## 1. 10x / droplet scRNA — autochthonous lung GEMM, open matrix

| Accession | Genotype | n mice | T/NK | Open matrix | Download |
|---|---|---|---|---|---|
| GSE180963 | K (`KrasG12D/+`) vs KL (`KrasG12D/+;Lkb1fl/fl`). Two samples were pooled at library build and split by label. | 1 K + 1 KL | yes, mixed | yes. Cell Ranger matrix inside `GSE180963_RAW.tar` | [RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180963/suppl/GSE180963_RAW.tar) |
| GSE165641 | KL (`KrasLSL-G12D/+;Lkb1fl/fl`), 10 weeks after Ad-Cre. Rep2 title says Lkb2; characteristics say Lkb1. | 2 KL | yes, mixed | yes. Normalized Rdata plus 10x RAW.tar | [Rdata](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165641/suppl/GSE165641_processed_normalized_matrix_data.Rdata.gz) |
| GSE179501 | KT; Lkb1-XTR. Total viable cells. Restored n=2 (CM2260, CM2319) vs non-restored n=2 (CM2324, CM2328). | 4 | yes, mixed | yes. 10x MTX + barcodes + features | [MTX](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179501/suppl/GSE179501_XTR_scRNAseq_matrix.mtx.gz) |
| GSE179502 | Same Lkb1-XTR model. FACS neoplastic cells. Restored vs non-restored. | 6 | no | yes. 10x MTX | [MTX](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179502/suppl/GSE179502_XTR_sorted_scRNAseq_matrix.mtx.gz) |
| GSE154977 | KP (`Kras` mutant, `Trp53` null), 30-week tumors. 2 no-drug + 2 cisplatin 72 h. AT2-lineage FACS. | 4 KP | no | yes. rawCount and normTPM h5 | [raw h5](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154977/suppl/GSE154977_mmLung10x_cis_dSp_rawCount.h5) |
| GSE154978 | Same KP study. Tigit-sorted pre/post transplant. | 6 libraries (mouse n not 1:1 with GSM) | no | yes. rawCount and normTPM h5 | [raw h5](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154978/suppl/GSE154978_mmLung10x_Tigit_dSp_rawCount.h5) |
| GSE264739 | KP (`LSL-KrasG12D;Trp53fl/fl`) n=3 vs KPP (same + Rosa26 LSL-Plk1) n=3. Lung tumor digest. | 6 | yes, mixed (not marked as a sort) | yes. 10x matrices in RAW.tar | [RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264739/suppl/GSE264739_RAW.tar) |
| GSE266323 | KP/KPC LUAD. Malat1 CRISPRa (d10) n=2 vs Tomato n=2. GEX libraries only (ADT rows are not RNA). | 4 | yes, mixed | yes. Per-library MTX triples | [dTom_1 MTX](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE266nnn/GSE266323/suppl/GSE266323_dTom_1_MMT_matrix.mtx.gz) |
| GSE322632 | KP lung ± sgLkb1, plus non-transduced lung. Hashed CD45+ and extravascular CD45−. Paper n=2 per arm; GEO has 3 GEX libraries. | 2 per arm (hashed) | yes, mixed, but CD45-heavy | yes. Seurat RDS (~1.9 Gb) | [RDS](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE322nnn/GSE322632/suppl/GSE322632_LKB1_SeuratObj.RDS) |
| GSE322633 | sgLkb1 KP ± sgLif or sgLifr. CITE-seq RDS. | 3 libraries | yes, mixed | yes. Seurat RDS | [RDS](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE322nnn/GSE322633/suppl/GSE322633_sgLif_SeuratObj.RDS) |
| GSE322634 | sgLkb1 KP, isotype vs anti-LIF. CITE-seq. | 2 libraries | yes, mixed | yes. Seurat RDS | [RDS](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE322nnn/GSE322634/suppl/GSE322634_antiLif_SeuratObj.RDS) |
| GSE201247 | Whole lung. ATTAC WT n=2, Kras adenoma n=2, ATTAC;Kras adenoma n=2. Not KL. | 6 | yes, mixed | yes. Filtered 10x h5 inside RAW.tar | [RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE201nnn/GSE201247/suppl/GSE201247_RAW.tar) |
| GSE317576 | KP whole lung after influenza or PBS, plus tumor-free controls. Summary states KP; GSM titles do not repeat the genotype. | 6 libraries | yes, mixed | yes. 10x RAW.tar | [RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE317nnn/GSE317576/suppl/GSE317576_RAW.tar) |
| GSE133604 | KP lung ± Asf1a KO ± anti-PD-1. One library per arm. | 4 libraries (mice per arm not stated) | yes, mixed (not marked as a sort) | yes. 10x RAW.tar plus shared genes.tsv | [RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE133nnn/GSE133604/suppl/GSE133604_RAW.tar) |
| GSE149813 | K (`KrasLSL-G12D/WT`), 7 weeks. YFP− and YFP+ distal epithelium from 2 mice. | 2 | no | yes. 10x RAW.tar | [RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149813/suppl/GSE149813_RAW.tar) |
| GSE158109 | KP (`Kras-LSLG12D/+;P53f/f`) lung tumors ± oxaliplatin/cyclophosphamide ± CAR-T. One CAR-product library is not a tumor. | 4 tumor libraries | yes, mixed | yes. counts.txt.gz + meta | [counts](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE158nnn/GSE158109/suppl/GSE158109_counts.txt.gz) |
| GSE228965 | KrasLA1/+ lung tumors ± Dicer1S2D, plus one Dicer-only lung. | 5 libraries | yes, mixed (primary tumor, not marked as a sort) | yes. 10x RAW.tar | [RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE228nnn/GSE228965/suppl/GSE228965_RAW.tar) |
| GSE281744 | CCSP-rtTa; TetO-KrasG12D; Trp53LSL-R172H. Whole lobe, KRAS ON vs ON-then-OFF. Each library pools 2 mice. | 4 mice in 2 libraries | yes, mixed | yes. 10x RAW.tar (Cell Ranger 6) | [RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE281nnn/GSE281744/suppl/GSE281744_RAW.tar) |
| GSE183875 | One `Kras mutant` lung (CK, 6 weeks) vs 2 control lungs. | 1 K | yes, mixed | yes. 10x RAW.tar | [RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183875/suppl/GSE183875_RAW.tar) |
| GSE188436 | KP (`KrasFSF-G12D;p53frt`) vs KPF1F2 (same + Foxa1/2 flox). Pre/post depletion. | 4 libraries | not marked; tumor-cell design | yes. 10x RAW.tar | [RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE188nnn/GSE188436/suppl/GSE188436_RAW.tar) |
| GSE231674 | 10x multiome RNA. KT, KPT, KFT. Two libraries each. | 6 libraries | not marked | yes. RAW.tar | [RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE231nnn/GSE231674/suppl/GSE231674_RAW.tar) |
| GSE231675 | SpcPT vs SpcT lung. Two libraries each. | 4 libraries | not marked | yes. RAW.tar | [RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE231nnn/GSE231675/suppl/GSE231675_RAW.tar) |
| GSE246481 | shKras GEMM lung time course and KP lung ± MRTX, plus subcutaneous KP. HTO rows duplicate GEX rows. | 9 GEX libraries (HTO; mouse n > libraries) | tumor-cell libraries; T/NK not established | yes. `shKrasGEMM_norm_allgene.h5` plus RAW.tar | [h5](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246481/suppl/GSE246481_shKrasGEMM_norm_allgene.h5) |
| GSE277777 | KP (`Kras-G12D;Trp53`) with Slc4a11 or Hopx tracers. Cell type is tumor cells. GEX+HTO pairs. | many hashed mice; not 1 GSM = 1 mouse | no | yes. Per-batch h5ad plus RAW.tar | [series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277777) |
| GSE319598 | KcP Cas9 lung, vehicle vs sotorasib. eGFP+ tumor cells (7 libraries) plus 2 pooled eGFP− stromal libraries. | 7 tumor libraries + 2 stromal pools | T/NK only in the stromal pools | yes. `adata_tumor_final.h5ad` (tumor cells) plus RAW.tar | [h5ad](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE319nnn/GSE319598/suppl/GSE319598_adata_tumor_final.h5ad) |
| GSE253461 | KPY (`KrasLSL-G12D;p53fl/fl`) AT2, organoids, and in vivo sorted AT2 / mesenchyme / endothelium / immune at 4–16 weeks. | in vivo mouse n not in GSM titles | T/NK only in the immune libraries, separate from AT2 | yes. 10x RAW.tar | [RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253461/suppl/GSE253461_RAW.tar) |
| GSE282851 | KT (`KrasLSL-G12D;R26LSL-Tomato`) FACS neoplastic cells, in vivo CRISPR. | 10 libraries, not 10 mice | no | yes. DGE MTX | [wt MTX](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE282nnn/GSE282851/suppl/GSE282851_wt_DGE.mtx.gz) |
| GSE281964 | KP vs KPH (Hnf4a flox), pre/post depletion. Lung tumor. | 4 libraries | not marked | yes. RAW.tar | [RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE281nnn/GSE281964/suppl/GSE281964_RAW.tar) |
| GSE294571 | KN vs KNH (KrasFSF-G12D, Nkx2-1 flox ± Hnf4a). | 2 libraries | not marked | yes. RAW.tar | [RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE294nnn/GSE294571/suppl/GSE294571_RAW.tar) |
| GSE234471 | Lung tumor, vehicle vs dTAGV-1 (KRAS G12V degradation). | 2 libraries | not marked | yes. RAW.tar | [RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE234nnn/GSE234471/suppl/GSE234471_RAW.tar) |
| GSE295452 | snRNA. KrasLSL-G12D MADM p53 lung: WT lung vs GFP tumor. | 2 libraries | yes, nuclei of the whole sample | yes. RAW.tar plus cell-type table | [RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE295nnn/GSE295452/suppl/GSE295452_RAW.tar) |
| GSE300862 | SPC KP vs SPC KPI lung tumors. | 2 libraries | not marked | yes. RAW.tar | [RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE300nnn/GSE300862/suppl/GSE300862_RAW.tar) |
| GSE203447 | KRAS-driven lung, 8 weeks. mCherry+ (1) and mCherry− (2). | 3 libraries | yes, mixed lung cells | yes. Annotated Seurat rds | [rds](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE203nnn/GSE203447/suppl/GSE203447_seu_preprocessed_s1-3_CCycle_removed_annotated.rds.gz) |
| GSE297023 | Kras(L/+);Tom;Cas9 tumor-bearing lung, young vs old, ± Pten. Parse split-pool. 16 sublibraries of mixed samples, not 16 mice. | hashed; not 16 mice | yes, mixed | yes. h5ad (~9.9 Gb) | [h5ad](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE297nnn/GSE297023/suppl/GSE297023_Aging_Parse_All_Samples.h5ad) |
| GSE295824 | Orthotopic Sox2-high lines including Nkx2-1−/−;Lkb1−/− (SNL, p53, Pten). Not a classic KL GEMM. 4 lines × 4 reps. | 16 libraries | yes, mixed tumor | yes. RAW.tar | [RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE295nnn/GSE295824/suppl/GSE295824_RAW.tar) |

Sister plate series, not 10x: **GSE154989** is Smart-seq2 of the same KP/K timecourse (FACS epithelium, T/NK no). GEO has 7,282 cell-level GSMs. The public sample table is 30 mice (K and KP). Open rawCount and normTPM h5: [raw h5](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154989/suppl/GSE154989_mmLungPlate_fQC_dSp_rawCount.h5). Super-series GSE152607 points at these children.

---

## 2. 10x / droplet — Kras or Lkb1 syngeneic lung lines (not autochthonous GEMM)

| Accession | Genotype | n mice | T/NK | Open matrix | Download |
|---|---|---|---|---|---|
| GSE267321 | Subcutaneous LKR13. K (KrasG12D) n=2, KK (Keap1 KO) n=2, KLK (Keap1/Lkb1 KO) n=2. Non-malignant cells. | 6 tumors | yes, mixed | yes. Normalized CSV. **Cldn4 checked in a prior public score: nearly empty (6/7956 cells).** | [CSV](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE267nnn/GSE267321/suppl/GSE267321_Normalized_expression_matrix_02122026.csv.gz) |
| GSE285606 | Subcutaneous 344SQ (`KrasG12D;Trp53R172H`) parental vs PD1R1. Both arms are IgG, not on-treatment PD-1. | 2 + 2 | yes, mixed | yes. RAW.tar | [RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE285nnn/GSE285606/suppl/GSE285606_RAW.tar) |
| GSE157881 | HKP1 (Kras/p53) orthotopic lung. CD45− and CD45+ at 0 Gy and 4 Gy. Mice pooled. No PD-1 on these libraries. | 4 pooled libraries | T/NK in the CD45+ libraries; Cldn4 in the CD45− libraries | yes. RAW.tar | [RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157881/suppl/GSE157881_RAW.tar) |
| GSE194166 | Subcutaneous KP vs KPL (STK11/LKB1). GEX + TCR. Sorted immune. | 5 GEX libraries | T/NK only | yes. RAW.tar. Cldn4 not expected. | [RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE194nnn/GSE194166/suppl/GSE194166_RAW.tar) |
| GSE232730 | KPL-3M (`KrasG12D;Tp53+/-;Lkb1-/-`) tumors. All four libraries are CD45+. | 4 libraries | T/NK only | yes. RAW.tar. Cldn4 sparse. | [RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE232nnn/GSE232730/suppl/GSE232730_RAW.tar) |
| GSE275877 | Subcutaneous LKR13-H (Kras, Msh2-KO). Isotype, PD-1 response, combo response, acquired resistance. | 5 libraries, 1 per arm | yes, mixed | unfiltered droplet dump (barcodes + features; ~34M barcodes in the prior public note). Not a QC matrix. SRA is the raw reads. | [features](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE275nnn/GSE275877/suppl/GSE275877_features.tsv.gz) |
| GSE127465 | Mouse arm is inDrop, not 10x. Orthotopic KP1.9 vs healthy lung. CD45+ only. 2 tumor-bearing + 2 healthy. | 4 | T/NK only | yes. Mouse normalized MTX. Cldn4 is the wrong compartment. | [series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127465) |

---

## 3. In vivo bulk — open matrix, K / KP / KL / KPL / Lkb1 lung

| Accession | Genotype | n mice | T/NK | Open matrix | Download |
|---|---|---|---|---|---|
| GSE6135 | Ji et al. microarray. Mouse tumors: K, KP, K;Ink4a, KL (Lkb1 L/L, L/−, L/+), plus one KL metastasis. Human A549/H2126 arrays are a second platform. | 25 mouse tumors; some mice have two nodules | yes, genes | yes. Mouse series matrix (GPL8321). Cldn4 probe 1418283_at. | [series matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE6nnn/GSE6135/matrix/GSE6135-GPL8321_series_matrix.txt.gz) |
| GSE137396 | Autochthonous nodules. KL (`KrasG12D;Lkb1fl/fl`) n=5 vs KP (`KrasG12D;Trp53fl/fl`) n=5. | 5 + 5 nodules | yes, genes | yes. Raw and normalized gene tables | [raw](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137396/suppl/GSE137396_Raw_genetable_GEMMnodule.txt.gz) |
| GSE193895 | Whole lung tumors. K/Keap1: 5 mice × 2 tumors. KL: 5 mice × 2 tumors. K/Keap1/Lkb1: 5 mice, 9 tumors. Plus 6 Cre-naive lungs. | 5 KL mice (10 tumors) | yes, genes | yes. Genewise counts | [counts](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE193nnn/GSE193895/suppl/GSE193895_Tumor_KrasModel_GEMMs_RNAseq_GenewiseCounts.txt.gz) |
| GSE175479 | Lung tumors. K n=4, KL n=4, K;AMPKα1/α2 n=4. | 4 K + 4 KL | yes, genes | yes. Raw gene counts | [counts](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE175nnn/GSE175479/suppl/GSE175479_Raw_gene_count_matrix.txt.gz) |
| GSE164758 | Primary tumors, untreated. K n=9, KL n=9, KP n=8, KPL n=15. | 41 tumors | yes, genes | yes. FPKM | [FPKM](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/GSE164758/suppl/GSE164758_primary_tumors_fpkm.txt.gz) |
| GSE277929 | KL primary lung tumors. Vehicle and trametinib/entinostat arms. | 35 tumors | yes, genes | yes. FPKM | [FPKM](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277929/suppl/GSE277929_fpkm.txt.gz) |
| GSE133714 | Affymetrix MoGene. Tumors: K n=4, KL (Stk11) n=4, K/Keap1 n=4, K/Keap1/Lkb1 n=4, K/Nrf2 n=3, plus WT lung n=4. | 4 KL tumors | yes, genes | yes. Series matrix (29 Mb). Cldn4 is on Mouse Gene 2.0 ST. | [series matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE133nnn/GSE133714/matrix/GSE133714_series_matrix.txt.gz) |
| GSE114601 | KP GEMM lung nodules. Vehicle, JQ1, anti-PD-1, combo. 2 mice per arm. | 8 | yes, genes | yes. Raw and normalized counts | [raw](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE114nnn/GSE114601/suppl/GSE114601_counts.raw.csv.gz) |
| GSE179500 | FACS neoplastic cells. KT, Lkb1-XTR restored vs not, some also Trp53-null. Several mice have more than one library. | 32 libraries, fewer mice | no | yes. DESeq2 normalized counts | [norm counts](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179500/suppl/GSE179500_20210620_XTR_Bulk_RNAseq_DESeq2norm_Counts.txt.gz) |
| GSE322570 | KP lung, sgNeo (Lkb1 WT) n=3 vs sgLkb1 n=3. Same study as GSE322632. | 3 + 3 | yes, genes (not marked as a sort) | yes. Raw counts | [counts](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE322nnn/GSE322570/suppl/GSE322570_raw_count.txt.gz) |
| GSE253613 | KL (`KrasG12D;Lkb1-/-`) lung tumors, G6pd WT vs KO, 12 weeks. | 15 tumors | yes, genes | yes. Raw counts | [counts](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253613/suppl/GSE253613_raw_counts.txt.gz) |
| GSE118246 | Lung tumors: KP adenocarcinoma n=8, Lkb1/Pten SCC, Sox2;Lkb1 SCC, Lenti-Sox2-Cre;Lkb1, plus 3 normal lungs. | 8 KP tumors among 34 libraries | yes, genes | yes. Mouse counts | [counts](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE118nnn/GSE118246/suppl/GSE118246_Mouse_counts.txt.gz) |
| GSE204752 | Kras-G12C;Trp53-KO GEMM lung. Control n=8 vs sotorasib-resistant n=8. | 16 | yes, genes | yes. counts.txt.gz | [counts](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE204nnn/GSE204752/suppl/GSE204752_counts.txt.gz) |
| GSE139347 | KP lung adenocarcinoma n=4, Eml4-Alk n=4, normal lung n=3. | 4 KP | yes, genes | yes. counts CSV | [counts](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE139nnn/GSE139347/suppl/GSE139347_lung_tumor_RNAseq_counts.csv.gz) |
| GSE161609 | K (`KrasG12D`) whole tumors ± sgMga. Sample text says tumor cells plus immune and stroma. Some mice contribute two tumors. | 13 tumors | yes, genes | yes. counts | [counts](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161609/suppl/GSE161609_Kras_tumor.counts.txt.gz) |
| GSE94927 | GEMM lung tumors: K, KP, Kras/Smarca4, Kras/p53/Smarca4. | 12 tumors | yes, genes | yes. Raw-count xls | [xls](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE94nnn/GSE94927/suppl/GSE94927_RawCount_Matrix.xls.gz) |
| GSE306017 | KP n=4 vs KP;Usp22-KO n=4 lung tumors. | 8 | yes, genes | yes. Counts matrix | [counts](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE306nnn/GSE306017/suppl/GSE306017_Counts_Matrix.txt.gz) |
| GSE250082 | CD45− cells from KP vs KP;Usp25−/− lungs, 10 weeks. | 2 + 2 | no | yes. All-gene table (21,932 genes) with per-sample values. **Cldn4 checked** (padj 1, not a Cldn4 hit). SRA also has 4 runs (PRJNA1052004). | [xls](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE250nnn/GSE250082/suppl/GSE250082_KP25vsKP_compare.xls.gz) |
| GSE271713 | Lung, Kras vs Kras;Prtn3. Replicate labels only. | 2 + 2 | yes, genes | yes. counts.txt.gz | [counts](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE271nnn/GSE271713/suppl/GSE271713_counts.txt.gz) |
| GSE110108 | KrasG12D;Lkb1 lung SCC, vehicle vs mTOR inhibitor. Microarray. | 6 | yes, genes | yes. Series matrix | [series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110108) |
| GSE111313 | Lkb1 in lung SCC development. Microarray. | 10 | yes, genes | yes. Series matrix | [series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111313) |
| GSE111335 | Lkb1/Pten lung SCC and adenosquamous. Microarray. | 9 | yes, genes | yes. Series matrix | [series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111335) |
| GSE54351 | Lkb1/Pten lung SCC, epithelial fraction. Microarray. Sister GSE54352 is stroma (T/NK genes, not tumor Cldn4). | 9 epithelial | no (epithelial fraction) | yes. Series matrix | [series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54351) |

Cell-line bulk that is genotype-labeled and open, not an autochthonous tumor:

| Accession | Genotype | n | T/NK | Open matrix | Download |
|---|---|---|---|---|---|
| GSE137244 | Nodule-derived lines. Five B6AL10 libraries, KL-named lines (KL155, KL47, KLC, KLD, KLE), and one normal lung. This is the cell-line series, not the nodule bulk (that is GSE137396). | 11 libraries | no (cultured lines) | yes. Raw, normalized, and FPKM CSV | [raw](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137244/suppl/GSE137244_counts.raw.csv.gz) |
| GSE338923 | Lacun3 murine LUAD cells, parental vs STK11-deficient. In vitro. | 8 libraries | no | yes. Raw counts | [counts](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE338nnn/GSE338923/suppl/GSE338923_Lacun3_STK11_RNAseq_raw_counts.tsv.gz) |

---

## 4. Looked closed — mirrors checked until the matrix or the end of public files

| Accession | What GEO showed | What the extra search found | Usable for Cldn4? |
|---|---|---|---|
| GSE244452 | Only `KPvsKL_deg_all.txt.gz`. Subcutaneous KL vs KP cell-line tumors, 3 vs 3. | The file is a 7,450-gene table with columns KP1–3 and KL1–3, not a significant-only list. **Cldn4 checked** (present). Figshare 26674682 is a cytokine figure, not the matrix. SRA has 6 runs (PRJNA1023194) if the full transcriptome is required. Paper text says FPKM was uploaded to this GEO series. | yes, for genes in that table, including Cldn4. Not the full transcriptome. [file](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE244nnn/GSE244452/suppl/GSE244452_KPvsKL_deg_all.txt.gz) |
| GSE182228 | Only `RAW.tar` (3.0 Mb) plus SRA. Subcutaneous LKB1-deficient LUAD, 3 mice × 4 arms (vehicle, palbociclib, anti-PD-1, combo). | The tar is 12 per-mouse FPKM files, not FASTQ. **Cldn4 checked** in vehicle S388 (FPKM 0 in that mouse; the row exists). Nature Communications data availability names this accession. No second matrix on figshare. | yes. Whole-transcriptome FPKM inside the tar. [RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE182nnn/GSE182228/suppl/GSE182228_RAW.tar). T/NK: yes, genes. |
| GSE250082 | Only a “compare” xls. | The sheet is 21,932 genes with KO1, KO2, WT1, WT2. Cldn4 is present. SRA PRJNA1052004 has 4 runs. | yes. See table 3. CD45−, so T/NK is no. |
| GSE303269 | `RAW.tar` (73 Mb), KP vs KP;Cic lung adenocarcinomas, 10 libraries. | filelist is stringtie GTF only, not a gene-count matrix. SRA PRJNA1294735 has 10 runs. No count table on GEO. | no ready matrix. Cldn4 would require parsing GTF FPKM attributes or counting the FASTQ. |
| GSE275877 | Barcode and feature files only. | Unfiltered 10x dump. No QC matrix in the paper supplement that GEO indexed. | not a usable matrix without heavy filtering. |
| E-MTAB-13704 | ArrayExpress / BioStudies lung GEMM RNA-seq, aPD-L1 combinations. | Study text does not name Kras, Lkb1, or Stk11. Prior public note: SDRF genotype field is wild type, and a raw-count CSV exists in the study files. | not added as a K/KP/KL set. |
| E-MTAB-5311 | LLC ± exercise, open Expression Atlas counts. | Lewis lung, not a Kras/Lkb1 GEMM. | out of genotype scope. |

---

## Sets that are public but not a Cldn4 tumor matrix

- GSE192806: Nkx2.1-Cre;Lkb1 lung at 5 weeks vs Lkb1 flox. Not a tumor GEMM. RAW.tar is open.
- GSE149909: KrasLSL-G12D AT2 organoids infected in vitro. Two pooled samples. Not tumor-bearing mice.
- GSE274477: sorted AT2, KRAS WT vs G12V, 2 runs each. Open feature-count TSV. T/NK no. Early lesion, not an advanced KL/KP tumor.
- LaFave et al. KP/KPT lung (GSE145192, GSE134812, GSE145194) is scATAC, not RNA. Cldn4 is not an expression measurement there.
