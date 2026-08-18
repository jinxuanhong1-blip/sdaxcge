# Methods — GSE127465 mouse Cldn4

## Data

GEO series [GSE127465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127465) mouse supplementary files only:

- `GSE127465_mouse_counts_normalized_15939x28205.mtx.gz` (cells × genes, author-normalized)
- `GSE127465_mouse_cell_metadata_15939x12.tsv.gz`
- `GSE127465_gene_names_mouse_28205.tsv.gz`

Human MTX, human metadata, human gene names, RAW.tar, and SRA are not used.

## Design (taken from the deposit / paper)

CD45+ inDrops from lungs of 2 tumor-free and 2 KP1.9 tumor-bearing mice (Zilionis et al., *Immunity* 2019). After author QC: 15,939 cells. Major types are author labels (Neutrophils, B cells, MoMacDC, T cells, NK cells, pDC, Basophils).

## Scores

- Cldn4 = mean log1p of the author-normalized value in **all** cells of that mouse. No epithelial subset is defined by the authors.
- T fraction = author `T cells` / cells.
- T/NK fraction = (`T cells` + `NK cells`) / cells.
- Myeloid fraction = (`Neutrophils` + `MoMacDC` + `pDC` + `Basophils`) / cells.
- Leftover epithelium (audit only) = Epcam>0 or (Cdh1>0 and Krt8>0).

## Statistics

Spearman on **mice** (n=4). Fisher-z 95% CI. Tumor-only n=2: points only. p-values are descriptive. Tacstd2 is audited, not a gate.

## Reproduce

```bash
python3 methods/gse127465_mouse_cldn4/scripts/download.py --out /tmp/gse127465_mouse
python3 methods/gse127465_mouse_cldn4/scripts/analyze.py \
  --data /tmp/gse127465_mouse \
  --outdir methods/gse127465_mouse_cldn4/results \
  --finding methods/gse127465_mouse_cldn4/FINDING.md
```
