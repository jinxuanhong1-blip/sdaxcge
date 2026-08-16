# Reproduce A3 GSE205335

Install the Python requirements, download only the three open processed/metadata
files listed below, and run `analyze.py`. Do not fetch controlled EGA raw data.

```bash
python3 -m pip install -r scripts/w200/A3_GSE205335/requirements.txt
mkdir -p /tmp/gse205335
curl -L --fail -o /tmp/gse205335/GSE205335_Lung_IO_CellIdentity.txt.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_CellIdentity.txt.gz
curl -L --fail -o /tmp/gse205335/GSE205335_Lung_IO_UMI_matrix.rds.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz
curl -L --fail -o /tmp/gse205335/GSE205335_family.soft.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/soft/GSE205335_family.soft.gz

python3 scripts/w200/A3_GSE205335/analyze.py \
  --matrix /tmp/gse205335/GSE205335_Lung_IO_UMI_matrix.rds.gz \
  --identities /tmp/gse205335/GSE205335_Lung_IO_CellIdentity.txt.gz \
  --soft /tmp/gse205335/GSE205335_family.soft.gz \
  --outdir results/w200/A3_GSE205335
```

GEO's matrix supplement is double-gzipped. The script removes the outer layer
in a temporary directory, parses the sparse RDS, extracts TACSTD2/EPCAM/PTPRC,
and never writes the full matrix into the repository.
