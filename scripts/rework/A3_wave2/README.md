# A3 wave-2 — GSE207422 malignant TACSTD2 recompute

Public GEO UMI only. Author CopyKAT IDs are not deposited.

```bash
python3 scripts/rework/A3_wave2/download.py
python3 scripts/rework/A3_wave2/analyze.py \
  --matrix data/GSE207422/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz \
  --metadata data/GSE207422/geo_scRNAseq_sample_metadata.tsv \
  --gene-chr scripts/rework/A3_wave2/gene_chr.tsv \
  --drmref scripts/rework/A3_wave2/drmref_cell_annotation.tsv.gz \
  --outdir results/rework/A3_wave2
```

Outputs land in `results/rework/A3_wave2/`.
