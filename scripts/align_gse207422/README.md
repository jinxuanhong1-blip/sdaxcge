# GSE207422 USER-ALIGN (TACSTD2 / TROP2)

Replicates the user's unpublished slide-5 scRNA design on the public GEO
processed matrix. Outputs go only to `results/align_gse207422/`.

## Data

- GEO: [GSE207422](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207422)
- Hu et al., *Genome Medicine* 2023. ~92k cells, neoadjuvant PD-1 + chemo NSCLC.
- Processed UMI matrix 175 MB gzipped (<2 GB). No raw FASTQ.

## Run

```bash
python3 scripts/align_gse207422/00_download.py
python3 scripts/align_gse207422/01_build_anndata.py
python3 scripts/align_gse207422/02_process_annotate.py
python3 scripts/align_gse207422/03_analyze_tacstd2.py
```

Large `.h5ad` / `.txt.gz` files are gitignored under `data/gse207422/`.
