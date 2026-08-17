# methods/gse148071_cellchat_cldn4

Additive CellChat-style ligand–receptor comparison of **CLDN4-high vs CLDN4-low** epithelial/malignant cells against **T/NK** on public **GSE148071** (Wu et al. 2021; 42 advanced NSCLC biopsies).

**CLDN4 only.** TACSTD2 is not used to define groups. Epithelial cells are a **putative** malignant compartment (no GEO labels).

```bash
python3 methods/gse148071_cellchat_cldn4/scripts/download.py --out methods/gse148071_cellchat_cldn4/data
python3 methods/gse148071_cellchat_cldn4/scripts/analyze.py \
  --data methods/gse148071_cellchat_cldn4/data \
  --out methods/gse148071_cellchat_cldn4/results
```

| Path | Role |
| --- | --- |
| [FINDING.md](FINDING.md) | Honest n, kept split, ligand table |
| [METHODS.md](METHODS.md) | Dataset, lineage, Hill probability, permutation, split rule |
| [db/](db/) | CellChatDB v2 protein pairs + complex table |
| [scripts/analyze.py](scripts/analyze.py) | Stream 42 UMI matrices → lineage → CellChat-like LR on CLDN4 |
| [results/ligand_table.tsv](results/ligand_table.tsv) | Significant differential pairs (kept split) |
| [results/fig_extra_ligand_table.png](results/fig_extra_ligand_table.png) | Extra figure: drawn ligand table |

The 172 MB GEO tar is not stored in git.
