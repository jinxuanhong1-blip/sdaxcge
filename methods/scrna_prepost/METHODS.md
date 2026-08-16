# Methods: public scRNA pre vs post ICI (additive extra)

Zhejiang paired TROP2 IHC (user A7) is taken as given and is not re-scored here. This extra is a public single-cell catalog plus one sample-unit test.

## Search

NCBI GEO (GDS) was queried on 2026-08-16 for human series matching lung/NSCLC × scRNA × ICI/PD-1/neoadjuvant, and lung × scRNA × pre/post/paired. Must-try accessions named in the request (GSE207422, GSE337519, GSE205335, GSE241934, GSE146100) plus 2023–2026 leftovers with Pre/Post in the title or design were opened (SOFT + supplementary file lists). Inclusion required a **public processed tumor-cell matrix** with **both** a pre-ICI/chemo-IO timepoint and a post timepoint. Unmatched patients were allowed. T-sorted, blood, post-only, pre-only, and DAC-controlled matrices were recorded and not tested.

## GSE207422 (only testable cohort)

Hu et al., *Genome Medicine* 2023 (PMID 36869384). Author UMI matrix `GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz` (24,292 genes × 92,330 cells) and sample metadata `GSE207422_NSCLC_scRNAseq_metadata.xlsx`. Three pre-treatment biopsies and twelve post-treatment resections are **different patients**.

Lineages were scored from canonical marker modules (T, NK, B, plasma, myeloid, neutrophil, mast, epithelial, fibroblast, endothelial) on log1p(CP10K). Malignant-like cells = epithelial lineage and not a clear normal-lung program (alveolar SFTPA/SFTPC/AGER, club SCGB1A1/SCGB3A2, or ciliated TPPP3/FOXJ1/CAPS log1p ≥ 1). GEO does not deposit barcode-level author annotations or CopyKAT objects.

**Unit = sample / patient** (one library per patient). Samples with fewer than 10 malignant-like cells are dropped, matching the source paper. For each remaining sample, TACSTD2 and CLDN4 were summarized as (i) mean log1p(CP10K), (ii) percent positive (UMI > 0), (iii) pseudobulk CPM. Pre vs post used two-sided Mann–Whitney U (unpaired). A residual-tumor sensitivity compared the three pre samples to NMPR post only. An epithelial-wide sensitivity dropped the normal-lung filter.

## GSE337519

Series design states one paired pre/post patient after three cycles of chemo-IO. GEO deposits a single 10x library (GSM9856929). No timepoint or hashing labels. The library is QC’d (cell counts, malignant-like TACSTD2/CLDN4 means) and is **not** entered into the pre vs post test.

## Combined estimate

Direction count and Stouffer’s Z across **independent tumor-cell cohorts** that actually have both timepoints. With one such cohort, Stouffer Z is the signed normal deviate of that Mann–Whitney p (up after treatment = positive). No random-effects meta-analysis is claimed.

## What this extra is not

It is not a paired IHC analysis and does not revisit the Zhejiang cohort.
