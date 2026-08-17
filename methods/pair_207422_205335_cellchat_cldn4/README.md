# methods/pair_207422_205335_cellchat_cldn4

ADDITIVE **CLDN4-only** pairwise merge of **GSE207422 + GSE205335** (both ICI-adjacent).
Include 207422. No dual-high.

Patient-level malignant CLDN4 vs same-patient T/NK (Spearman + Q4 vs Q1), then
CellChat-style and LIANA-style outgoing CLDN4-high → T/NK. Honest n.

```bash
python3 -m pip install -r methods/pair_207422_205335_cellchat_cldn4/requirements.txt
python3 methods/pair_207422_205335_cellchat_cldn4/download.py --out /tmp/geo
python3 methods/pair_207422_205335_cellchat_cldn4/analyze.py \
  --matrix /tmp/geo/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz \
  --out methods/pair_207422_205335_cellchat_cldn4
```

`--skip-matrix` writes combo ρ and the GSE205335-only LR reuse without streaming the UMI.

| Path | Role |
| --- | --- |
| [FINDING.md](FINDING.md) | Combo ρ + LR tables |
| [METHODS.md](METHODS.md) | Quartiles, Hill P, mean-of-means |
| [data/](data/) | Locked patient extracts + given 205335 per-patient LR |
| [results/combo_rho.tsv](results/combo_rho.tsv) | Fisher-z combo |
| [results/ligand_table_cellchat_outgoing.tsv](results/ligand_table_cellchat_outgoing.tsv) | Combo CellChat-style outgoing |
| [results/ligand_table_liana_outgoing.tsv](results/ligand_table_liana_outgoing.tsv) | Combo LIANA-style outgoing |
