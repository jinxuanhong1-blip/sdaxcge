# methods/scrna_twoway

Additive extra: public lung scRNA series that have **both** a timepoint and an MPR/RECIST label, then malignant TACSTD2/CLDN4 in that 2×2.

```bash
python3 methods/scrna_twoway/hunt_geo.py
python3 methods/scrna_twoway/analyze_gse207422.py \
  --matrix data/scrna_twoway/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz \
  --metadata data/scrna_twoway/GSE207422_NSCLC_scRNAseq_metadata.xlsx \
  --gene-chr methods/scrna_twoway/gene_chr.tsv \
  --outdir results/scrna_twoway
```

GEO downloads (not committed):

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207422
