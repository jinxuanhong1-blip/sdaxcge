# Leftover public mouse lung-tumor scRNA hunt (Seurat, Cldn4-only)

Additive catalog of **public** mouse lung-tumor **scRNA** with **Cldn4** plus immune or GEMM labels. Assigned accessions were cataloged, not re-scored as the only analysis. Private 8 KL mice were not used. FASTQ was not downloaded. Scoring used **R 4.3.3 + Seurat 5.0.1**.

## Assigned (not re-scored here)

GSE179502, GSE154989, GSE267321, GSE154977, GSE165641, GSE180963, GSE127465, GSE274477.

## Catalog

`methods/seurat_mouse_kl_hunt/catalog.tsv` — leftover hunt plus assigned/no-go rows.

## Scored leftovers (Seurat)

Best unused processed leftovers (not the assigned eight):

1. **GSE179501** — Lkb1-XTR **total-viable** sister of assigned GSE179502. Combined 10x MTX. **n=4 mice** (2 Restored, 2 Non-Restored).
2. **GSE201247** — Kras / ATTAC whole-lung 10x h5. **n=6 mice** (2 WT, 2 Kras, 2 ATTAC;Kras).
3. **GSE266323** — KP LUAD TME, Malat1 CRISPRa vs Tomato. **n=4 mice**.

Per-accession folders: `methods/seurat_mouse_kl_hunt/GSE179501/`, `GSE201247/`, `GSE266323/` each with `mouse_level_scores.tsv`.

## Methods

**Cldn4 only.** Seurat `NormalizeData` (log1p, scale 1e4) and `AddModuleScore`.

**T/NK genes:** Cd3d/e/g, Cd2, Cd8a/b1, Cd4, Nkg7, Gzma/b, Prf1, Klrb1c, Ncr1, Klrd1, Klrc1, Ifng (present genes only).

**IFN/MHC genes:** Stat1/2, Irf1/7/9, Isg15, Ifit1/2/3, Mx1, Oasl2, Rsad2, Ifih1, Ddx58, Ifnb1, B2m, H2-K1/D1/Q4/Q6/Q7, H2-Aa/Ab1/Eb1, Tap1/2, Psmb8/9, Nlrc5, Ciita.

**Epithelial cells:** Epcam>0 or (Krt8>0 and Ptprc==0). **T/NK cells:** Cd3e/Cd3d/Nkg7/Ncr1>0.

**Unit:** one biological mouse. QC: 200–10000 features and ≥500 UMI. High vs low is a median split of epithelial Cldn4 among mice with ≥10 epithelial cells. Spearman + two-sided Mann–Whitney. Honest n; no multiple-testing theater.

## Results

Filled after Seurat scoring (see per-accession `mouse_level_scores.tsv` and `honest_n.tsv`).

## Re-run

```bash
# processed GEO files in /tmp/geo_dl (URLs in catalog.tsv)
Rscript methods/seurat_mouse_kl_hunt/scripts/score_leftovers.R
```
