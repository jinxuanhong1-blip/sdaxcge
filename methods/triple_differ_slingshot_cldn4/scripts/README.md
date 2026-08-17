# Scripts

1. `download.py` — public GEO processed files only. GSE123902+GSE131907+GSE205335.
   No GSE148071. No 7-pool. No FASTQ / EGA / 2.86 GB TPM / 36.5 GB Laughney H5.
2. `extract_epithelium.py` — author epithelium (GSE131907, GSE205335) plus
   marker epithelium (GSE123902) → joint h5ad, cap ≤350 cells/unit.
3. `analyze.py` — try Harmony on the triple; if Harmony dies, per-dataset
   graphs. AT2-rooted Slingshot if R is present, else scanpy DPT. Root is
   never CLDN4-high. Writes the stacked sample-level CLDN4 vs pseudotime table.
4. `run_slingshot.R` — real Slingshot (Street 2018) on a reduced-dimension
   matrix + cluster + start cluster.
5. `gene_sets.py` — locked AT2 / barrier (CLDN4 excluded) / comparator sets.
