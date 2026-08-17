# methods/markermalig15_cldn4_hiend

ADDITIVE **CLDN4-only** high-end on the marker-malignant combo that already differs (GSE253013 + GSE291670, **n=15**, ρ=−0.714 given — not re-audited).

No dual-high. CellChat-style + NicheNet-style on CLDN4-high vs low marker-malignant → T/NK. Extra figures.

```bash
# GSE291670 RAW tar and GSE253013 RDS are not stored in git.
python3 methods/markermalig15_cldn4_hiend/scripts/extract_gse253013.py \
  --rds /tmp/markermalig15_data/GSE253013/GSE253013_all_luad_garnett_temp.rds.gz \
  --outdir /tmp/markermalig15_data/GSE253013/extracted \
  --panel methods/markermalig15_cldn4_hiend/data/lr_gene_panel.txt
python3 methods/markermalig15_cldn4_hiend/scripts/assemble_gse253013.py \
  --extracted /tmp/markermalig15_data/GSE253013/extracted
python3 methods/markermalig15_cldn4_hiend/scripts/convert_prior.py \
  --lt-rds /tmp/nichenet_priors/ligand_target_matrix_nsga2r_final.rds \
  --out /tmp/nichenet_priors/ligand_target.parquet
python3 methods/markermalig15_cldn4_hiend/scripts/analyze.py
```

| Path | Role |
| --- | --- |
| [FINDING.md](FINDING.md) | Honest n=15, LR table, NicheNet activity |
| [METHODS.md](METHODS.md) | Combo, CLDN4-only split, Hill probability, prior |
| [data/combo/samples_n15.tsv](data/combo/samples_n15.tsv) | The 15 samples from existing combo notes |
| [results/lr_table.tsv](results/lr_table.tsv) | CellChat-style differential LR table |
