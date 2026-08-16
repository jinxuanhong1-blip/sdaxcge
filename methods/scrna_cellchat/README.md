# methods/scrna_cellchat

CellChat-style ligand–receptor comparison of **TACSTD2-high vs TACSTD2-low** epithelial/malignant cells against **T/NK** on public **GSE207422** (one UMI matrix; MPR labels present).

A3/B6 are taken as given. This folder is the communication add-on, not a claim audit.

```bash
python3 methods/scrna_cellchat/scripts/download.py --out methods/scrna_cellchat/data
python3 methods/scrna_cellchat/scripts/analyze.py \
  --matrix methods/scrna_cellchat/data/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz \
  --meta methods/scrna_cellchat/data/GSE207422_NSCLC_scRNAseq_metadata.xlsx \
  --out methods/scrna_cellchat/results
```

| Path | Role |
| --- | --- |
| [METHODS.md](METHODS.md) | Dataset, lineage, Hill probability, permutation, split rule |
| [RESULTS.md](RESULTS.md) | n, significant LR counts, kept split (written after the run) |
| [db/](db/) | CellChatDB v2 protein pairs + complex table |
| [scripts/analyze.py](scripts/analyze.py) | Stream UMI → lineage → CellChat-like LR |
| [results/](results/) | Tables, JSON, figures from the public matrix |

The 176 MB GEO matrix is not stored in git.
