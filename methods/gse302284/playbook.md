# Playbook — public 10x H5 scoring of TACSTD2 / CLDN4 / TJ / IFN

Methods only. No results, no numbers.

## Estimand

For a deposited treatment contrast, estimate the direction of TACSTD2, CLDN4, a pre-specified tight-junction (TJ) signature, and IFN / MHC-I / APM between treated and control libraries.

If the accession has no TROP2-ADC / SKB264 / sacituzumab libraries, say so and score the contrast that exists (for GSE302284: osimertinib vs vehicle). Do not substitute another accession.

## Experimental unit

The library (GEO sample) is the biological replicate. Cells inside a library are not independent biological replicates. When n_library = 1 vs 1, there is no sample-level *p*-value. Report pseudobulk log2FC as the primary effect size. Cell-level tests are exploratory and must be labeled as such.

## Matrices

Prefer author-processed Cell Ranger filtered H5 from GEO supplementary files. Do not pull controlled FASTQ. Record URL, bytes, and MD5.

## QC (pre-specified)

Start from the filtered H5. Then, per library:

1. Drop cells in the bottom 10% of UMI counts.
2. Drop cells in the bottom 10% of detected genes.
3. Drop cells with mitochondrial UMI fraction > 25%.

Optional doublet calling is not required for reuse. If skipped, say so.

## Compartments

- Cell-line / PDX libraries aligned to GRCh38: all QC cells.
- Mixed patient tissue: restrict tumor-cell tests to epithelial cells (EPCAM / keratin high relative to PTPRC). Report how many cells remain.

## Gene sets (fixed)

- Targets: TACSTD2, CLDN4.
- TJ signature: CLDN1, CLDN4, CLDN7, F11R, PARD3, OCLN, TJP1. Score = mean of per-cell log1p(CP10k) over genes present.
- PRIORITY6: IFI27, OAS2, IFIT1, MX1, ISG15, HLA-A.
- MHC1_APM and IFN_ISG: lists in `gene_sets.py`.

Do not add or drop genes after seeing the values.

## Effects to report

For each gene and each set:

- n_library treat / control
- n_cells treat / control after QC (and compartment)
- percent positive (UMI > 0)
- mean log1p(CP10k)
- pseudobulk CPM and log2FC
- exploratory cell-level two-sided Mann–Whitney U, rank-biserial, Welch *t*

For co-expression, Spearman ρ of TACSTD2 vs CLDN4 / TJ / IFN / APM **within** each library. Always pair ρ with n and *p*.

## Honesty

- Do not treat osimertinib vs vehicle as a TROP2-ADC analog.
- Do not treat mRNA as protein.
- Do not fabricate accessions or *p*-values.
- Missing genes stay missing.
