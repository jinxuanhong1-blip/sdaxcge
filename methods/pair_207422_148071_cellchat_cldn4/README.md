# methods/pair_207422_148071_cellchat_cldn4

ADDITIVE **CLDN4-only** pairwise merge of **GSE207422 + GSE148071**.

1. Patient-level malignant/epithelial CLDN4 vs T/NK (per cohort + combo rho).
2. CellChat-style outgoing CLDN4-high → T/NK (given ligand tables; consensus LR).

**Include 207422.** No dual-high. TACSTD2 is not a gate. Honest n = 25 + 12 = **37**, not 42+15.

Prior CellChat folders are taken as given (PR #324, PR #348). Matrices are not re-downloaded.

```bash
python3 methods/pair_207422_148071_cellchat_cldn4/scripts/analyze.py
```

| Path | Role |
| --- | --- |
| [FINDING.md](FINDING.md) | Combo rho, honest n, consensus ligand table |
| [METHODS.md](METHODS.md) | Eligibility, Spearman / DL RE, consensus rule |
| [results/combo_rho.tsv](results/combo_rho.tsv) | Patient-level ρ (primary + sensitivity) |
| [results/lr_table.tsv](results/lr_table.tsv) | Consensus LR (both sig, same sign) |
| [results/lr_table_outgoing.tsv](results/lr_table_outgoing.tsv) | Outgoing Mal → T/NK only |
| [figures/](figures/) | Extra figures (scatter, forest, n, ligand table) |
| [data/SOURCE.txt](data/SOURCE.txt) | Where the given tables came from |
