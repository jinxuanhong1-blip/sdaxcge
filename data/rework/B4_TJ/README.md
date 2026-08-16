Raw GEO matrices are downloaded by `scripts/rework/B4_TJ/download.py` and are gitignored.

Pinned after download:

- `download_manifest.json` — URLs, SHA256, byte sizes
- `kegg_hsa04530_genes.txt` — parsed KEGG hsa04530 symbols
- `symbol_to_ensembl.json` — mygene.info symbol → Ensembl (current)

GSE135222 is hg19; `scripts/rework/B4_TJ/gene_sets.py` carries GRCh37 Ensembl fallbacks (CLDN7 = ENSG00000181885).
