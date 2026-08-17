# CXCL13+ trio — CLDN4-only high-end CellChat (patient-level)

Additive analysis on the CXCL13+ trio that already differs
(GSE148071 + GSE207422 + GSE253013, **n=60 ρ=−0.425 given**, not re-audited).

Outgoing CellChat-style probability from **CLDN4-high malignant** cells
toward **CXCL13+ T** and toward **T/NK**. Patient is the unit. Honest n.
No dual-high. Not GSE207422-only.

Write-up: `FINDING.md`. Ligand table: `results/ligand_table.tsv`.

```bash
python3 methods/cxcl13_trio_cldn4_hiend/scripts/download.py
# GSE253013 only (after the 9.3 GB RDS is present):
python3 methods/cxcl13_trio_cldn4_hiend/scripts/extract_gse253013_rds.py \
  --rds /tmp/cxcl13_trio_data/GSE253013/GSE253013_all_luad_garnett_temp.rds.gz \
  --outdir /tmp/cxcl13_trio_data/GSE253013/extracted
python3 methods/cxcl13_trio_cldn4_hiend/scripts/assemble_gse253013.py \
  --extracted /tmp/cxcl13_trio_data/GSE253013/extracted
python3 methods/cxcl13_trio_cldn4_hiend/scripts/analyze.py
```
