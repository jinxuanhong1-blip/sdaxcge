# methods/gse205335_cldn4_extra

ADDITIVE **GSE205335-only**. Malignant **CLDN4** vs same-patient **T/NK** at
**Q4 vs Q1**, plus CellChat-style ligand–receptor scores on the same split.

Prior TACSTD2 A3 (`results/w200/A3_GSE205335`) and the multi-cohort CLDN4 Q4
meta (`methods/cldn4_malig_q4_tnk`) are taken as given. This folder does not
re-rank those pools. **CLDN4 only** — TACSTD2 is not used to define groups.

```bash
python3 -m pip install -r methods/gse205335_cldn4_extra/requirements.txt
python3 methods/gse205335_cldn4_extra/download.py --out /tmp/gse205335
python3 methods/gse205335_cldn4_extra/analyze.py \
  --patients methods/gse205335_cldn4_extra/data/GSE205335_patients.tsv \
  --matrix /tmp/gse205335/GSE205335_Lung_IO_UMI_matrix.rds.gz \
  --identities /tmp/gse205335/GSE205335_Lung_IO_CellIdentity.txt.gz \
  --soft /tmp/gse205335/GSE205335_family.soft.gz \
  --out methods/gse205335_cldn4_extra
```

`--skip-cellchat` writes the Q4 vs T/NK table from the given patient extract
without the 500 MB RDS.

| Path | Role |
| --- | --- |
| [FINDING.md](FINDING.md) | Honest n, Q4 vs T/NK table, ligand table |
| [METHODS.md](METHODS.md) | Quartiles, author labels, Hill probability |
| [data/GSE205335_patients.tsv](data/GSE205335_patients.tsv) | Given patient extract (PR #320 / #279) |
| [db/](db/) | CellChatDB v2 protein pairs (same parse as GSE207422 CellChat) |
| [results/q4q1_tnk.tsv](results/q4q1_tnk.tsv) | Q4 vs Q1 T/NK tests |
| [results/ligand_table.tsv](results/ligand_table.tsv) | Differential CellChat-style pairs |
