# methods/scrna_cellchat_cldn4

Additive CellChat-style ligand–receptor comparison of **CLDN4-high vs CLDN4-low** epithelial/malignant cells against **T/NK** on public **GSE207422**.

Prior TACSTD2 CellChat (`methods/scrna_cellchat`) is given and is not re-run. **CLDN4 only.**

```bash
python3 methods/scrna_cellchat_cldn4/scripts/download.py --out methods/scrna_cellchat_cldn4/data
python3 methods/scrna_cellchat_cldn4/scripts/analyze.py \
  --matrix methods/scrna_cellchat_cldn4/data/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz \
  --meta methods/scrna_cellchat_cldn4/data/GSE207422_NSCLC_scRNAseq_metadata.xlsx \
  --out methods/scrna_cellchat_cldn4/results
```

| Path | Role |
| --- | --- |
| [FINDING.md](FINDING.md) | Honest n, kept split, ligand table |
| [METHODS.md](METHODS.md) | Dataset, lineage, Hill probability, permutation, split rule |
| [db/](db/) | CellChatDB v2 protein pairs + complex table |
| [scripts/analyze.py](scripts/analyze.py) | Stream UMI → lineage → CellChat-like LR on CLDN4 |
| [results/ligand_table.tsv](results/ligand_table.tsv) | Significant differential pairs (kept split) |
| [results/fig_extra_ligand_table.png](results/fig_extra_ligand_table.png) | Extra figure: drawn ligand table |

The 176 MB GEO matrix is not stored in git.
