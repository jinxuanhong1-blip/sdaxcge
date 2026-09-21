# Public mouse K / KP / KL lung scRNA

Processed GEO matrices only. Private 8-KL counts are not used.

```bash
bash methods/public_mouse_klkp_epi_t/scripts/download.sh /tmp/public_mouse_scrna
python3 methods/public_mouse_klkp_epi_t/scripts/preprocess.py \
  --data /tmp/public_mouse_scrna \
  --out methods/public_mouse_klkp_epi_t
Rscript methods/public_mouse_klkp_epi_t/scripts/to_seurat.R \
  methods/public_mouse_klkp_epi_t /tmp/rlibs
```

`tables/` holds the public per-mouse scores, the Cldn4 threshold sweep, label checks, and per-cell metadata. `objects/*.h5ad` and `objects/*.rds` are regenerated locally and are not committed. Those objects are the QC'd cells restricted to the marker and program genes used for labels and scores.

Accessions scored: GSE165641, GSE180963, GSE154977, GSE179502, GSE179501 (unsorted sister of GSE179502), GSE154989, GSE267321, GSE127465, GSE136246, GSE149813.

Found and not scored: GSE277777 (tumor-cell objects; combined h5ad is 8 GB) and GSE319598 (stroma pooled across mice).
