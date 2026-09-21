# Geneformer / scGPT: CLDN4-high vs low malignant cells (GSE131907)

Additive CLDN4-only embedding of author-malignant cells from GSE131907. The sample is the unit. This is not the concordant-4 n=65 T/NK correlation and not a spatial exclusion test.

```bash
python3 -m venv /tmp/venv
/tmp/venv/bin/pip install -r methods/geneformer_cldn4_gse131907/requirements.txt
GF_CACHE=/tmp/gf_cache /tmp/venv/bin/python methods/geneformer_cldn4_gse131907/analyze.py
```

`analyze.py` expects these files under `GF_CACHE` (default `/tmp/gf_cache`):

| file | source |
|---|---|
| `GSE131907_Lung_Cancer_cell_annotation.txt.gz` | GEO GSE131907 |
| `GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz` | GEO GSE131907 |
| `gf_v1/` (`config.json`, `model.safetensors`) | `ctheodoris/Geneformer` Geneformer-V1-10M |
| `token_dictionary_gc30M.pkl`, `gene_median_dictionary_gc30M.pkl`, `gene_name_id_dict_gc30M.pkl` | Geneformer gc30M dictionaries |
| `scgpt/best_model.pt`, `scgpt/vocab.json` | `wanglab/scGPT-human` |

Optional, benchmarked and skipped if slower than 1.5 s/cell: Geneformer-V2-104M_CLcancer plus the gc104M dictionaries.

Design: `METHODS.md`. Numbers: `FINDING.md` (written by the script).
