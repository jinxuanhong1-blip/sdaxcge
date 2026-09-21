# GSE180581 inputs

- Series: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180581
- Counts: `GSE180581_all_samples_counts.xlsx` (STAR raw counts, hg38). Stored here as `GSE180581_counts_ensembl.csv.gz`.
- Symbol map: NCBI `Homo_sapiens.gene_info` Ensembl cross-references (`ensembl_to_symbol.tsv`). Counts for one symbol are summed.
- Gene sets: Enrichr `MSigDB_Hallmark_2020`, `Reactome_2022`, `KEGG_2021_Human`, `GO_Biological_Process_2023`, plus the repository MHC-I/APM list and a curated cytosolic DNA-sensor list. Frozen in `genesets.json`.

GEO sample titles for GSM5465261–GSM5465266 say `293T_KuDNA-PKcs`. Sample characteristics say cell line `293T_DNA-PKcs` and treatment `siControl` or `siDNA-PKcs`. The series growth protocol describes monoallelic DNA-PKcs knockout. Those six libraries are the DNA-PKcs arm.

Design, from the series matrix: each subunit is a monoallelic HEK293T knockout plus 50 nM matching siRNA for 72 h, versus the same monoallelic line plus siControl. Three libraries per arm. No exogenous DNA is in the protocol.

Data note: Anisenko et al., Data in Brief 2021 (PMID 34849385). The authors' own differential-expression summary is Anisenko et al., Biochimie 2022 (PMID 35430316).
