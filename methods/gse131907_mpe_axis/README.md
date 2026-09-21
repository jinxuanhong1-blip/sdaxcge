# GSE131907 MPE vs primary: CLDN4, TACSTD2, ELF3, EPCAM

Public fallback after GSA-Human **HRA006761** (controlled access, DAC HDAC002197). Numbers are in [RESULTS.md](RESULTS.md).

The run reads the Kim et al. 2020 processed UMI matrix (`GSE131907`). It does not use the 2.9 GB log2TPM text or EGA FASTQ. HRA006761 is not downloaded.

```bash
bash methods/gse131907_mpe_axis/scripts/download.sh /tmp/gse131907
python3 methods/gse131907_mpe_axis/scripts/extract_genes.py \
  --datadir /tmp/gse131907 --outdir /tmp/gse131907/mpe_axis
python3 methods/gse131907_mpe_axis/scripts/analyze_mpe.py \
  --datadir /tmp/gse131907 --subset /tmp/gse131907/mpe_axis \
  --outdir methods/gse131907_mpe_axis/results
```

`extract_genes.py` streams the full matrix once, keeps 15 genes, and stores each cell’s full-transcriptome UMI sum. `analyze_mpe.py` log-normalizes with that full library size: `log1p(UMI / nCount_full * 10000)`.
