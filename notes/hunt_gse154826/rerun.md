# Re-run GSE154826 TACSTD2 vs LCAM

Author cell barcodes and LCAM definitions are vendored under `data/leader_metadata/` (from https://github.com/effiken/Leader_et_al). Raw 10x matrices are not in git (~3.4 GB).

```bash
# 1. GEO tarballs
mkdir -p data/geo/tars
cd data/geo/tars
curl -sS https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154826/suppl/ \
  | grep -oE 'GSE154826_amp_batch_ID_[0-9]+\.tar\.gz' | sort -u \
  | xargs -P 4 -I{} curl -sS -O --retry 5 \
      https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154826/suppl/{}
cd ../../..

# 2. extract annotated cells + gate transcriptomes + immune pbulk
python3 scripts/01_extract_batches.py

# 3. LCAM / TACSTD2 / CD8 / CITE-seq
python3 scripts/02_analyze.py
```

Outputs land in `results/hunt_gse154826/{tables,figures}/`.

Python: pandas, numpy, scipy, matplotlib, pyarrow. No Seurat/scanpy required for the reported numbers.
