# FINDING — Seurat + CellChat concordant-four CLDN4-high vs low senders

ADDITIVE. **Thesis already correct. Ligands stay.**
CLDN4 only. No dual-high. Concordant four only
(GSE123902 + GSE131907 + GSE205335 + GSE189357).
Do **not** add GSE148071. This is **not** a Python reimplementation.

Analysis not yet written — `scripts/run_cellchat.R` overwrites this file
after R CellChat `computeCommunProb` on per-patient CLDN4-high vs low
senders → T/NK. Honest unit = patient. If Seurat or CellChat cannot
install, the run stops and this file is the no-go.
