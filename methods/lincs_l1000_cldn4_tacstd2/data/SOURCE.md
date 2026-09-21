# LINCS / CMap L1000 source

Public Broad LINCS 2020 beta release and the GEO phase I / phase II signature catalogs. No API key. The Level 5 gctx files are multi-gigabyte; this folder keeps a 20-signature slice, not the archive.

| File | URL | Role |
|---|---|---|
| `siginfo_beta.txt` | https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/siginfo_beta.txt | 1,201,944 signature rows. Gene symbol is `cmap_name`. |
| `geneinfo_beta.txt` | https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/geneinfo_beta.txt | 12,328 genes (978 landmark, 9,196 best inferred, 2,154 inferred). |
| `cellinfo_beta.txt` | https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/cellinfo_beta.txt | Lineage for the 10 TACSTD2 cell lines. |
| `level5_beta_trt_xpr_n142901x12328.gctx` | https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/level5/level5_beta_trt_xpr_n142901x12328.gctx | CRISPR Level 5 moderated z. 6,516,812,529 bytes. Columns are signatures, rows are genes. |
| `GSE92742_Broad_LINCS_sig_info.txt.gz` | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_sig_info.txt.gz | Phase I signature catalog (473,647). |
| `GSE92742_Broad_LINCS_pert_info.txt.gz` | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_pert_info.txt.gz | Phase I perturbagen catalog (51,383). |
| `GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz` | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz | Phase II signature catalog, 2017-03-06 (118,050). |

`download.py` reads the CRISPR gctx with HTTP `Range` requests (`h5py` fileobj driver) and writes only the TACSTD2 `trt_xpr` rows. The run that produced this slice fetched **175,024,881** bytes in **42** requests.

Committed slice:

| File | Contents |
|---|---|
| `tacstd2_xpr_level5_z.tsv.gz` | 12,328 genes × 20 signatures. sha256 `b1db7915370ee714af098eb2d8db295ec3688485ad5fc6130ad6ce2d1fd82d57` |
| `tacstd2_xpr_sig_meta.tsv` | siginfo + cellinfo columns for those 20 rows, plus gctx row index |
| `slice_manifest.tsv` | byte count and checksum |

ShRNA Level 5 (`level5_beta_trt_sh_n238351x12328.gctx`, ~11.7 GB) was not downloaded. CLDN4 and TACSTD2 are absent from `trt_sh`, `trt_sh.cgs`, and `trt_sh.css` in `siginfo_beta.txt`, so that file has no columns to slice.

Data use: Broad CMap / LINCS, https://clue.io/connectopedia/data_use_policy . Level 5 moderated z-scores are the 2020 beta release (Subramanian et al., *Cell* 2017, for the L1000 assay; LINCS 2020 build for these signature ids).

Phase I gene field is `pert_iname`. LINCS 2020 gene field is `cmap_name`. Counts are exact string matches.
