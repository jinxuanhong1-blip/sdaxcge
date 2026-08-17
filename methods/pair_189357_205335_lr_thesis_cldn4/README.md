# pair_189357_205335_lr_thesis_cldn4

ADDITIVE CLDN4-only ligand test aligned to the thesis, on the PR #459 pair
that differs (GSE189357+GSE205335). No dual-high. Patient is the unit.

Given (not re-audited): %pos n=31, ρ=−0.478, Q4 vs Q1 r=−0.750.

Done when both family tables exist:

- `results/family_barrier_inhibitory.tsv`
- `results/family_ifn_recruit_mhci.tsv`

Write-up: `FINDING.md`.

```bash
python3 methods/pair_189357_205335_lr_thesis_cldn4/scripts/download.py
python3 methods/pair_189357_205335_lr_thesis_cldn4/scripts/analyze.py
```
