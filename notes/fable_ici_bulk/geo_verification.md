# GEO page verification (verified 2026-08-16)

Task: TACSTD2 & CLDN4 vs ICI response in **open human bulk lung** datasets.
Download **processed matrices only** (<2GB, no FASTQ).

| GSE | Platform | Assay | Processed file (size) | Notes on eligibility |
|-----|----------|-------|-----------------------|----------------------|
| GSE126044 | GPL16791 Illumina HiSeq 2500 (Homo sapiens) | bulk RNA-seq | `GSE126044_counts.txt.gz` (547 KB) | NSCLC anti-PD-1 (Cho et al). Human bulk lung. **Eligible** if response labels present. |
| GSE135222 | GPL16791 Illumina HiSeq 2500 (Homo sapiens) | bulk RNA-seq | `GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz` (1.6 MB) | NSCLC anti-PD-(L)1 (Jung et al 2019). Human bulk lung. **Eligible** (PFS/DCB). |
| GSE136961 | GPL24014 Ion Torrent S5 XL (Homo sapiens) | targeted RNA panel | `GSE136961_TPM.tsv.gz` (68 KB), `raw_count.tsv.gz` (29 KB) | Need to verify tissue = lung and ICI + whether TACSTD2/CLDN4 on panel. |
| GSE166449 | GPL11154 Illumina HiSeq 2000 (Homo sapiens) | bulk RNA-seq | `GSE166449_Raw_gene_TPM_matrix.txt.gz` (1.6 MB) | Need to verify tissue = lung and ICI response labels. |
| GSE93157 | GPL19965 nCounter PanCancer Immune Profiling Panel | NanoString (~730 genes) | `GSE93157_raw_data_values.txt.gz` (83 KB) | Prat et al 2017 (NSCLC/melanoma/HNSCC, anti-PD-1). Targeted immune panel — TACSTD2/CLDN4 likely NOT on panel. Verify. |
| GSE207422 | GPL24676 Illumina NovaSeq 6000 (Homo sapiens) | bulk RNA-seq (+scRNA separate) | `GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz` (5.4 MB) + `_metadata.xlsx` (12 KB) | NSCLC bulk. **Eligible**. Use bulk only (ignore 175 MB scRNA). |

All processed files are well under 2 GB; no FASTQ used.

## Eligibility outcome (after inspecting matrices)

| GSE | Design | Response variable | TACSTD2/CLDN4 present? | Verdict |
|-----|--------|-------------------|------------------------|---------|
| GSE126044 | 16 NSCLC, pre-tx anti-PD-1, bulk RNA-seq raw counts | responder / non-responder (binary) | **Yes** | **USABLE** |
| GSE135222 | 27 NSCLC, anti-PD-(L)1, bulk RNA-seq TPM | PFS event + pfs.time (survival) | **Yes** (ENSG00000184292 / ENSG00000189143) | **USABLE (survival)** |
| GSE136961 | 21 NSCLC anti-PD-1, Ion Torrent **Oncomine Immune Response** targeted panel (395 immune genes) | DCB/NDB + survival | **NO** (epithelial genes not on immune panel) | **UNUSABLE** — target genes not measured |
| GSE166449 | 22 lung cancer, pre-tx immunotherapy, bulk RNA-seq TPM | responder / non-responder (binary) | **Yes** | **USABLE** |
| GSE93157 | 65 pts (NSCLC/melanoma/HNSCC) anti-PD-1, **NanoString PanCancer Immune 730** panel | best RECIST resp + PFS | **NO** (epithelial genes not on immune panel) | **UNUSABLE** — target genes not measured |
| GSE207422 | NSCLC neoadjuvant anti-PD-1+chemo; 24 **bulk** RNA-seq (log2TPM) + 15 scRNA (ignored) | pathologic response MPR / NMPR (binary) | **Yes** | **USABLE (bulk only)** |

Confirmed by grep: GSE136961 and GSE93157 contain no TACSTD/CLDN/TROP-family symbols at all. Both are targeted immune-oncology panels that simply do not include the epithelial markers TACSTD2 (TROP2) or CLDN4, so they cannot be analyzed for this question regardless of sample size.

## Additional matching cohorts (hunted via NCBI eutils, round 2)

Goal per user direction: prefer cohorts like theirs (neoadjuvant IO, **durvalumab**, TCGA/OncoSG),
and use tumor-only / purity-corrected expression. `esearch` on GEO DataSets for
`durvalumab lung` and `neoadjuvant NSCLC immunotherapy RNA-seq`.

| GSE | Cohort | Platform / file | TACSTD2/CLDN4? | Response label in GEO? | Use |
|-----|--------|-----------------|----------------|------------------------|-----|
| **GSE253564** | Neoadjuvant **durvalumab ± SBRT**, early NSCLC, **pre-treatment** resected tumor | NovaSeq bulk RNA-seq, `Pre-treatment_Samples_Pubs_FPKMs` (n=32) | **Yes** | No (outcome only in paper) | **Correlation-only** (TACSTD2/CLDN4 vs CD8/NK, purity-corrected) |
| **GSE248378** | Neoadjuvant **durvalumab ± SBRT**, NSCLC, **post-treatment** resected tumor | NovaSeq bulk RNA-seq, `Durva_Post_FPKMs` (n=29) | **Yes** | No (outcome only in paper) | **Correlation-only** |
| GSE110390 | Durvalumab NSCLC + urothelial (IFN-γ signature) | NextSeq, **21-gene** `*_21gene_expression` panel | **No** (21-gene IFN-γ panel) | best RECIST in SRA | **Skip** (panel without TACSTD2) |
| GSE243013 | Anti-PD1 NSCLC single-cell atlas | scRNA-seq | — | — | Not bulk; out of scope |

cBioPortal (public) large-n tumor cohorts used for the mechanistic anti-correlation test:
- **TCGA-LUAD** `luad_tcga_pan_can_atlas_2018` (510 RNA-seq samples)
- **TCGA-LUSC** `lusc_tcga_pan_can_atlas_2018` (484)
- **OncoSG LUAD** `luad_oncosg_2020` (169 with mRNA)

The two durvalumab series carry only tissue / histology / treatment-arm in their GEO
`series_matrix`; no responder label is deposited (outcomes are in the trial paper). They are
therefore used to test the **anti-correlation with CD8/NK** (which needs no response label),
not a responder comparison.

