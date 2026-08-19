# Inventory — leftover OPEN lung ST (not CosMx demo)

Additive public spatial only. CLDN4-only. Skip controlled. Skip if CLDN4 absent. No private 8-KL.

| Accession / source | Platform | Tissue | Coords? | CLDN4 | Decision |
|---|---|---|---|---|---|
| 10x Visium HD Human Lung Cancer (IF FFPE) | Visium HD 16 µm bins | vendor lung cancer | yes (`feature_slice.h5`) | **4 UMIs total** | skip_CLDN4_too_sparse |
| [GSE301973](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301973) | Visium HD (StarDist cells) | EGFR-mutant NSCLC, 2 slides | yes | yes | **analyze** |
| [GSE328481](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE328481) | Stereo-XCR-seq | LUAD, 11 sections | yes (`obsm['spatial']`) | yes | **analyze** |
| [Zenodo 7306132](https://zenodo.org/records/7306132) (STopover / Na–Choi LUAD) | Visium | LUAD, 11 sections | yes (`tissue_positions_list.csv`) | yes | **analyze** |
| [GSE189357](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189357) | scRNA-seq | early LUAD atlas | no | — | skip_not_spatial (spatial sibling is already-used GSE189487) |
| [GSE200916](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200916) | Visium (10x barcodes) | MPLC | **no** `tissue_positions` | not run | skip_no_coords |
| [Figshare 32384337](https://doi.org/10.6084/m9.figshare.32384337) | Visium CytAssist (deposit) | KRAS G12C LUAD n=4 | plots / cluster tables only | not a spot matrix | skip_no_coords |
| [Zenodo 14624390](https://doi.org/10.5281/zenodo.14624390) Tavernari | Visium (claimed) | LUAD plasticity | files API empty | — | skip_controlled (`access_right=restricted`) |
| [Zenodo 14620362](https://zenodo.org/records/14620362) USZ TLS | Visium | kidney 3 + lung 5 | yes (Space Ranger) | not extracted | skip_incomplete (2.1 GB; STopover already covers Zenodo LUAD Visium) |
| [Zenodo 13337961](https://zenodo.org/records/13337961) | Visium CytAssist | lepidic/solid LUAD | rar downloads are JSON stubs | — | skip_no_files |
| [CNP0005129](https://db.cngb.org/data_resources/project/CNP0005129) | Stereo-seq | LUAD STAS | apply-for-data | — | skip_controlled |
| GSE276934 Visium HD | Visium HD | pulmonary fibrosis | 49 GB RAW | — | skip_size / not NSCLC leftover |

TSV: `tables/hunt_catalog.tsv`.
