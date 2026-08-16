# Data sources (fable_tisch slice)

All objects are public TISCH2 NSCLC downloads
(`https://tisch.compbio.cn/static/data/<DS>/<DS>_expression.h5` plus
`<DS>_CellMetainfo_table.tsv`). Values in the h5 are the TISCH2/MAESTRO
`log2(TPM/10 + 1)` matrix. We did **not** download the ~12 GB extended
LuCA atlas.

| TISCH2 id | GEO / ArrayExpress | Paper | On-disk h5 | Cells used | Malignant definition |
|---|---|---|---|---|---|
| NSCLC_GSE131907 | GSE131907 | Kim et al. 2020, *Nat Commun* | 2.69 GB (subsetted to tumor-tissue cells before stats) | 65,702 tumor-tissue / 203,298 total | proxy: tumor-tissue `Epithelial` (TISCH2 has no Malignant call) |
| NSCLC_EMTAB6149 | E-MTAB-6149 | Lambrechts et al. 2018, *Nat Med* | 416 MB | 40,218 | TISCH2 `Malignant` |
| NSCLC_GSE127465 | GSE127465 | Zilionis et al. 2019, *Immunity* | 304 MB | 31,179 | TISCH2 `Malignant` |
| NSCLC_GSE148071 | GSE148071 | Wu et al. 2021, *Nat Commun* | 1.15 GB | 82,267 | TISCH2 `Malignant` |

CELLxGENE lung-cancer collections were listed as a fallback
(`ad10cef8-…`, `0bebef1a-…`, `edb893ee-…`, `62e8f058-…`) but were not
needed once the four TISCH2 objects above were in hand.

Raw h5/tsv files live in `data_fable_tisch/` (gitignored). Reproduce with
`bash scripts/fable_tisch/run_all.sh`.
