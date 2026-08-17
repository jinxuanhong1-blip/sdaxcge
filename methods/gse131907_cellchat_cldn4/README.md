# methods/gse131907_cellchat_cldn4

Additive **CellChat-style** ligand–receptor comparison of **malignant CLDN4-high vs CLDN4-low** against **T/NK** on public **GSE131907**.

GSE207422 CellChat is a different folder / agent. This slice is **GSE131907 + CLDN4 only**.

```bash
python3 methods/gse131907_cellchat_cldn4/scripts/download.py --out /tmp/gse131907
python3 methods/gse131907_cellchat_cldn4/scripts/analyze.py \
  --matrix /tmp/gse131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz \
  --ann /tmp/gse131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz \
  --series /tmp/gse131907/GSE131907_series_matrix.txt.gz \
  --out methods/gse131907_cellchat_cldn4/results
```

| Path | Role |
| --- | --- |
| [FINDING.md](FINDING.md) | Honest n + ligand table |
| [METHODS.md](METHODS.md) | Author labels, Hill probability, permutation |
| [db/](db/) | CellChatDB v2 protein pairs + complex table |
| [scripts/analyze.py](scripts/analyze.py) | Stream UMI → CLDN4 tertile → CellChat-like LR |
| [results/ligand_table.tsv](results/ligand_table.tsv) | Primary outgoing Mal → T/NK table |

The 0.38 GB GEO UMI matrix is not stored in git.
