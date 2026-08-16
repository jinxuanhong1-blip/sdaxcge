# Honest skip: huge raw files not downloaded

This leftover B6 analysis uses **only small processed Space Ranger outputs**.
Huge raw / image / Loupe files were **not fetched** (bandwidth, disk, and they
are not required for CLDN4 vs immune-niche scoring).

Dataset page: https://www.10xgenomics.com/datasets/human-lung-cancer-ffpe-2-standard
CDN prefix: `https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma/`

This is **not** ArrayExpress **E-MTAB-13530**.

## Used (downloaded)

| File | Size | Why |
|---|---|---|
| `*_filtered_feature_bc_matrix.h5` | 24 MB | Spot x gene counts |
| `*_spatial.tar.gz` (lowres/hires PNG, positions, scalefactors) | 31 MB | Coordinates + tissue image for `sc.read_visium` |
| `*_metrics_summary.csv` | <1 KB | Published QC numbers |

The CytAssist TIFF inside `spatial/` (~21 MB) was copied locally but **not**
used by the script (scanpy uses the PNG pair).

## Skipped on purpose

| File | ~Size | Reason |
|---|---|---|
| `*_tissue_image.tif` | **3.0 GB** | Full-res H&E; lowres/hires PNGs are enough for spatial plots |
| `*_cloupe.cloupe` | **1.6 GB** | Loupe-only; not used by scanpy/squidpy |
| `*_molecule_info.h5` | **377 MB** | Molecule-level; not needed for spot-level niches |
| FASTQ / input files | multi-GB | Re-alignment not in scope |
| `*_raw_feature_bc_matrix.h5` | 31 MB | Unfiltered matrix; filtered H5 is the analysis input |
| `*_analysis.tar.gz` | 31 MB | Vendor clustering; we recompute Leiden |
| `*_image.tif` | 20 MB | Redundant with spatial PNG/TIFF |

Re-fetch only the used files with `scripts/fetch_b6_visium_extra.sh`.
Raw files are **not** in this repo (`data/` is gitignored).
