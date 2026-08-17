# methods/pair_207422_131907_cellchat_cldn4

ADDITIVE **CLDN4-only** pairwise merge of public **GSE207422** (Hu et al.,
*Genome Medicine* 2023) and **GSE131907** (Kim et al., *Nat Commun* 2020).

GSE207422 is included. No dual-high. Patient-level malignant CLDN4 vs T/NK
first (honest n, Q4 vs Q1), then CellChat-style outgoing CLDN4-high → T/NK.

```bash
python3 methods/pair_207422_131907_cellchat_cldn4/analyze.py
```

| Path | Role |
| --- | --- |
| [FINDING.md](FINDING.md) | Honest n, combo ρ, Q4 vs Q1, ligand table |
| [METHODS.md](METHODS.md) | Pair rules; no cell pooling across studies |
| [results/combo_rho.tsv](results/combo_rho.tsv) | Primary combo Spearman + Q4 vs Q1 |
| [results/combo_lr_table.tsv](results/combo_lr_table.tsv) | Full outgoing join |
| [results/ligand_table.tsv](results/ligand_table.tsv) | Significant-in-either + key pairs |
| [figures/](figures/) | Forests, scatters, Q4 boxes, extra LR figures |

Single-cohort CellChat and the Q4 T/NK extract are given and are not re-run.
