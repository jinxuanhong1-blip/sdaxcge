# GSE334497 source

- Accession: [GSE334497](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE334497)
- Supplementary matrix: `GSE334497_normalized_counts.csv.gz`
  - GEO FTP: `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE334nnn/GSE334497/suppl/GSE334497_normalized_counts.csv.gz`
  - Author-normalized counts, GRCm38 Ensembl gene IDs
- Groups, from GEO SOFT `Sample_description` library names:
  - Trop2 KO (n=5): KO162, KO164, KO165, KO172, RESUB-KO163R
  - WT (n=5): RESUB-171R, RESUB-170R, RESUB-169R, RESUB-168R, control170
- Design: 4T1 tumors, 3 weeks in BALB/c, frozen sections. Wu *et al.*, *J Immunother Cancer* 2026;14:e012265. https://jitc.bmj.com/content/14/4/e012265
- `ensembl_to_symbol.tsv`: NCBI `gene2ensembl` (taxid 10090) joined to `Mus_musculus.gene_info` official symbols, restricted to IDs in this matrix. Runtime anchors: Tacstd2 `ENSMUSG00000051397`, Cldn4 `ENSMUSG00000047501`, Cldn7 `ENSMUSG00000018569`, Cldn1 `ENSMUSG00000022512`, Cxcl9 `ENSMUSG00000029417`.
