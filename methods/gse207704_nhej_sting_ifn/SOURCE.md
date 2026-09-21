# Source

- Accession: GSE207704 (PMID 37059993)
- File: `GSE207704_CLDN4_RNAseq.txt.gz`
- URL: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207704/suppl/GSE207704_CLDN4_RNAseq.txt.gz
- Design on the GEO record: T47D WT GSM6310640–41, T47D CLDN4−/− GSM6310642–43, MCF7 WT GSM6310644–45, MCF7 CLDN4−/− GSM6310646–47
- The processed file has one FPKM column per group. Replicates are already pooled.
- WT is parental, “No treatment”. There is no siRNA negative-control arm.
- Reactome 2022 sets in `reactome_sets.gmt` were taken from the Enrichr library `Reactome_2022` (gene-set names and R-HSA ids are in the GMT description field).
- The matrix is downloaded by `scripts/gse207704_nhej_sting_ifn/analyze.py` and is not committed.
