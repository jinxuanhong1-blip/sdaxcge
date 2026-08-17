# Results — GSE131907 CellChat-style CLDN4 → T/NK

| File | Role |
| --- | --- |
| `ligand_table.tsv` | **Primary.** Outgoing malignant CLDN4-high vs low → T/NK (kept tumor tertile). 2,116 CellChatDB v2 protein pairs. |
| `ligand_table_tumor_tertile.tsv` | Same table, named by split |
| `ligand_table_tlung_tertile.tsv` / `ligand_table_mets_tertile.tsv` | Sensitivities |
| `n_cells_*.tsv` / `dominance_*.tsv` | Honest cell / sample / patient n and top-sample fraction |
| `lr_pairs_*.tsv` | All directed tests + permutation p |
| `contrast_*_incoming.tsv` | T/NK → Mal (not the ligand table) |
| `fig_top_outgoing_tumor_tertile.png` | Significant outgoing ΔP |
| `summary.json` / `inventory.json` | Machine-readable n and skips |

See [../FINDING.md](../FINDING.md) for the written table and n.
