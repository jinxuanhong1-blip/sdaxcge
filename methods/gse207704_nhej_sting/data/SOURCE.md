# GSE207704 processed matrix

`GSE207704_CLDN4_RNAseq.txt.gz` is the GEO supplementary cufflinks FPKM table (not FASTQ).

- Series: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207704
- File: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207704/suppl/GSE207704_CLDN4_RNAseq.txt.gz
- Paper: Murakami et al., Breast Cancer Research 2023 (CLDN4 CRISPR knockout in T47D and MCF7). DOI: 10.1186/s13058-023-01646-z

The deposit has one FPKM per genotype per line. GEO itself lists 2 biological replicates per genotype (GSM6310640–GSM6310647; SRA SRR20029118–SRR20029125), already collapsed in this file. The sweep re-quantified those public SRA objects with kallisto 0.51.1 against Ensembl 90 cDNA. FASTQ is not stored in the repo. Counts and PyDESeq2 tables are under `tables/`.

`ifn_apm_sets.json` holds Hallmark IFN-γ, Hallmark IFN-α, and the 21-gene MHC-I/APM list, with the same membership as the earlier CLDN4-loss prerank.
