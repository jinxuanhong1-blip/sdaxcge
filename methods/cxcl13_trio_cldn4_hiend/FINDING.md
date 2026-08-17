# FINDING — CellChat-style outgoing CLDN4-high malignant on the CXCL13+ trio

**CLDN4 only. No dual-high. Not GSE207422-only.** Patient is the unit.

The CXCL13+ trio that already differs is taken as given and is **not re-audited**:
GSE148071 + GSE207422 + GSE253013, **n=60, ρ=−0.425, p=0.00121**
(`methods/scrna_cldn4_combo`, PR #290, family `tls/cxcl13pos/mean`).

This folder adds outgoing CellChat-style probability from CLDN4-high vs
CLDN4-low putative malignant epithelium toward CXCL13+ T and toward T/NK.
Results land in `results/ligand_table.tsv` after `scripts/analyze.py`.

Analysis running / see `results/` for the LR table.
