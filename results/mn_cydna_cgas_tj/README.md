# Micronuclei, cytosolic DNA, and cGAS: claudin / tight-junction panel

Public GEO only. The script downloads count matrices to `/tmp/geo_mn` and writes the tables in this directory.

```bash
python3 scripts/mn_cydna_cgas_tj/analyze.py
```

Requires `pandas`, `numpy`, `scipy`, `statsmodels`, and `pydeseq2` (see `scripts/mn_cydna_cgas_tj/requirements.txt`).

## Tables

| File | Contents |
| --- | --- |
| `inventory.tsv` | Series and papers, including those with no matrix |
| `search_log.md` | E-utilities queries |
| `VERDICT.txt` | Calls and the numbers they rest on |
| `de_panel.tsv` | Full pre-specified panel, every contrast |
| `de_claudin_tj.tsv` | Claudin, tight-junction, and TACSTD2 rows |
| `de_key_genes.tsv` | Short gene list used while reading the run |

## Test

PyDESeq2 Wald test, contrast case versus control, Cook's filter and independent filtering off. `padj_panel` is Benjamini-Hochberg across the panel in that contrast. `padj_genome` is the PyDESeq2 adjusted p-value. A call of UP or DOWN needs `|log2FC| >= 0.5` and `padj_panel < 0.05`, and the two group means must sum to at least 1.

GSE100771 is a Mann-Whitney test on log2(CPM+1) across cells from one dish. GSE98183 75 bp PE is a one-library ratio with no p-value.

Contrasts are not merged across accessions. GSE296527 reuses the GSE279172 high-density wild-type libraries.
