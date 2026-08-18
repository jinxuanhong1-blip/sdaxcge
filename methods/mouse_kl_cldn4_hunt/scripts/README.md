# Re-run scored analyses

Put processed files in `/tmp/geo_dl` (no FASTQ):

```
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154989/suppl/GSE154989_mmLungPlate_fQC_dSp_normTPM.h5
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154989/suppl/GSE154989_mmLungPlate_fQC_geneTable.csv.gz
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154989/suppl/GSE154989_mmLungPlate_fQC_smpTable.csv.gz
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE6nnn/GSE6135/matrix/GSE6135-GPL8321_series_matrix.txt.gz
https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL8nnn/GPL8321/annot/GPL8321.annot.gz
https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/E-MTAB-5311/E-MTAB-5311-raw-counts.tsv
https://www.ebi.ac.uk/biostudies/files/E-MTAB-5311/E-MTAB-5311.sdrf.txt
```

Then:

```
python3 score_gse154989.py
python3 score_gse6135.py
python3 score_emtab5311.py
```
