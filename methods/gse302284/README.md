# methods/gse302284 — reusable GSE302284 / 10x H5 scoring

Public GEO [GSE302284](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE302284): EGFR-mutant NSCLC DTP scRNA (patient residual + DFCI282 / PC9 osimertinib vs vehicle).

SKB264 / CLDN4 thesis is taken as given. These scripts score TACSTD2, CLDN4, a tight-junction signature, and IFN/MHC-I/APM on whatever contrast is actually deposited.

## What GEO contains

Six Cell Ranger filtered H5 libraries. **No sacituzumab / SKB264 / TROP2-ADC treated RNA.** The treatment contrast that exists is osimertinib vs vehicle (n_library = 1 vs 1 in each model).

## Run

```bash
python3 -m pip install -r methods/gse302284/requirements.txt
python3 methods/gse302284/download.py --cache-dir /tmp/gse302284
python3 methods/gse302284/analyze.py --cache-dir /tmp/gse302284
```

Outputs: `results/C_GSE302284/` (WRITEUP, figures, tables). H5 files stay in the cache and are gitignored.

## Reuse

| File | Role |
|---|---|
| `gene_sets.py` | Fixed human TJ / PRIORITY6 / IFN_ISG / MHC1_APM lists |
| `download.py` | GEO RAW tar + provenance JSON |
| `analyze.py` | 10x H5 load, author-like QC, pseudobulk log2FC, cell-level MWU (exploratory), Spearman, figures |

`gene_sets.py` and the H5 helpers in `analyze.py` can be pointed at another 10x GEO series. Do not invent a TROP2-ADC contrast if the accession does not have those libraries.

## Stats contract

- Biological unit = **library**. GSE302284 has n = 1 vs 1.
- Primary effect size = pseudobulk log2FC = log2((CPM_treat+1)/(CPM_ctrl+1)).
- Cell-level MWU *p* and rank-biserial are exploratory (pseudoreplication).
- Report n_library, n_cells, direction, log2FC, and *p* together. No fabricated stats.
