# Public concordant lung atlas

Harmony UMAP of the four concordant public lung scRNA cohorts, colored by cluster cell class.

Cohorts: GSE123902, GSE131907, GSE205335, GSE189357.

## Coordinates that are in the repo

Cell-level UMAP coordinates and cluster labels are already committed:

`methods/seurat_concordant4_cldn4/results/objects/seurat_harmony_embeddings.rds`

- Branch: `cursor/seurat-concordant4-cldn4-802a`
- Blob: `01e4d80f28d8040e61ed90f66e3da4da26306bae` (4,130,411 bytes)
- Contents: `umap` (22,653 × 2), `harmony` (22,653 × 30), `meta` (22,653 × 11)

`data/umap_coordinates.tsv.gz` is that UMAP joined to `cell_class`, `dataset`, `unit_id`, and `seurat_cluster`. Values are copied from the object. None were filled in.

`meta` columns in the RDS are: `orig.ident`, `nCount_RNA`, `nFeature_RNA`, `dataset`, `unit_id`, `unit_type`, `cell_class_src`, `author_malignant`, `RNA_snn_res.0.6`, `seurat_clusters`, `cell_class`.

Figure: `fig_atlas_cellclass.pdf`, `fig_atlas_cellclass.svg`, `fig_atlas_cellclass.png`.

Counts taken from the coordinate table:

| cell class | cells |
|---|---:|
| malignant | 7,424 |
| T | 6,358 |
| myeloid | 4,220 |
| NK | 1,814 |
| B | 1,552 |
| other | 1,285 |

| cohort | units | cells in the UMAP |
|---|---:|---:|
| GSE123902 | 13 donors | 4,453 |
| GSE131907 | 21 samples | 7,350 |
| GSE205335 | 22 patients | 7,700 |
| GSE189357 | 9 patients | 3,150 |

65 units. The plotted cloud is 22,653 cells after the original QC and ≤350 cells/unit cap. The unit count is the inferential n used with this atlas; the cell count is how many points are drawn.

## CLDN4 / TACSTD2 feature panels

Per-cell CLDN4 and TACSTD2 are not in the committed embeddings, so they are not plotted.

The file that holds the RNA assay for these same cells, and that is not in the repository:

`methods/seurat_concordant4_cldn4/results/objects/seurat_harmony_integrated.rds`

On `cursor/seurat-concordant4-cldn4-802a`, `FINDING.md` describes that object as local (218 MB). The committed file is only the slim list of UMAP, Harmony, and metadata. `saveRDS` for the slim list does not store gene expression.

Also absent, and listed in `.gitignore` on `cursor/concordant4-atlas-umap-annotate-3655`:

`methods/concordant4_atlas_umap_annotate/results/tables/cell_obs.tsv`

That table is not in any commit. The Python atlas script would have written a `cldn4` column into it. It would still not have written TACSTD2, and the UMAP coordinates in that pipeline sit in `obsm["X_umap"]`, not in the table.

Full GEO count matrices were not downloaded for this figure.

## Reproduce

```bash
python3 nature_figures/fig_atlas_public/plot_fig_atlas.py
```

Reads `data/umap_coordinates.tsv.gz` only. Needs matplotlib and pandas.
