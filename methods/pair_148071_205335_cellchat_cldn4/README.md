# methods/pair_148071_205335_cellchat_cldn4

ADDITIVE **CLDN4-only** pairwise merge of **GSE148071 + GSE205335**.
GSE131907 is not included. No dual-high TACSTD2×CLDN4 score.

1. Patient-level malignant CLDN4 vs same-patient T/NK (honest n, Q4 vs Q1).
2. Combo Spearman (DerSimonian–Laird on the two cohort rhos).
3. CellChat-style outgoing CLDN4-high → T/NK on the same within-cohort tails.

```bash
python3 -m pip install -r methods/pair_148071_205335_cellchat_cldn4/requirements.txt
python3 methods/pair_148071_205335_cellchat_cldn4/scripts/download.py --out /tmp/gse148071
python3 methods/pair_148071_205335_cellchat_cldn4/scripts/analyze.py \
  --data-148071 /tmp/gse148071 \
  --out methods/pair_148071_205335_cellchat_cldn4
```

`--skip-cellchat` writes combo ρ / Q4 vs Q1 from the locked extracts
without streaming the GSE148071 matrices.

| Path | Role |
| --- | --- |
| [FINDING.md](FINDING.md) | Honest n, combo ρ, ligand table |
| [METHODS.md](METHODS.md) | Pair rule, quartiles, Hill probability |
| [data/](data/) | Locked patient extracts + GSE205335 per-patient LR |
| [results/combo_rho.tsv](results/combo_rho.tsv) | DL combo ρ |
| [results/ligand_table.tsv](results/ligand_table.tsv) | Outgoing Mal→T/NK Q4 vs Q1 |
