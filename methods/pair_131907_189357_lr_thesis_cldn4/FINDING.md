# FINDING — CLDN4-only thesis-aligned ligand test on GSE131907 + GSE189357

ADDITIVE. **CLDN4 only. No dual-high.** Patient is the unit.
Do **not** add GSE148071. Do **not** re-audit the T/NK Spearman.

Thesis (already correct):

- CLDN4-high → more barrier/inhibitory outgoing to T/NK (F11R, NECTIN2–TIGIT, CDH1, LGALS9)
- CLDN4-low / KD-like → more IFN / T-recruit / MHC-I outgoing (CXCL9/10–CXCR3, CCL5, HLA–CD8)

Given PR #459 pair: GSE131907+GSE189357 %pos n=30, ρ=−0.542, Q4 vs Q1 r=−0.619.

Family tables are written by `scripts/analyze.py` to:

- `results/family_barrier_inhibitory.tsv`
- `results/family_ifn_recruit_mhci.tsv`

This page is replaced after the run.
