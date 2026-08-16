# scRNA TLS / B vs TACSTD2 / CLDN4 meta

Additive public-only slice. User A6/B6 immune-cold is taken as given.

```bash
python3 methods/scrna_tls_meta/download.py
# GSE253013 RDS walk (9.3 GB; skip if the file is absent)
python3 methods/scrna_tls_meta/extract_gse253013_rds.py \
  --rds data/scrna_tls_meta/GSE253013/GSE253013_all_luad_garnett_temp.rds.gz \
  --outdir data/scrna_tls_meta/GSE253013/extracted
python3 methods/scrna_tls_meta/assemble_gse253013.py \
  --extracted data/scrna_tls_meta/GSE253013/extracted
python3 methods/scrna_tls_meta/extract.py
python3 methods/scrna_tls_meta/analyze.py
```

Outputs: `results/scrna_tls_meta/`. Methods-only playbook: `playbook.md`.
