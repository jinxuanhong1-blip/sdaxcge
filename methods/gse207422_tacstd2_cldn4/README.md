# GSE207422 TACSTD2 and CLDN4, scored separately

Public Hu et al. Genome Medicine 2023 neoadjuvant NSCLC scRNA-seq (GEO GSE207422).

Slide 5–6 use this cohort for TROP2-high → higher tumor fraction, lower CD8 T/NK, and an MPR association. This folder scores **TACSTD2**, **CLDN4**, and **dual-high** on those endpoints with the same rules. CLDN4 is not edited to match TACSTD2.

```bash
python3 methods/gse207422_tacstd2_cldn4/scripts/download.py --outdir data/GSE207422
python3 methods/gse207422_tacstd2_cldn4/scripts/analyze.py \
  --matrix data/GSE207422/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz \
  --metadata data/GSE207422/GSE207422_NSCLC_scRNAseq_metadata.xlsx \
  --outdir methods/gse207422_tacstd2_cldn4
```

The UMI matrix is not committed (see `data/GSE207422/.gitignore`). Numbers are written to `FINDING.md` by the analysis script.
