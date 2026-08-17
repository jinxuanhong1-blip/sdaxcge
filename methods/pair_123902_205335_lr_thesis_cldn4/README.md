# pair_123902_205335_lr_thesis_cldn4

ADDITIVE CLDN4-only ligand test aligned to the thesis, on the PR #459 pair
that differs (GSE123902+GSE205335). No dual-high. No GSE148071.
Do not re-audit the T/NK Spearman. Patient/donor is the unit.

Given (PR #459, not re-audited): %pos n=35, ρ=−0.522, Q4 vs Q1 r=−0.802.

Done when both family tables exist:

- `results/family_barrier_inhibitory.tsv`
- `results/family_ifn_recruit_mhci.tsv`

Write-up: `FINDING.md`.

```bash
python3 methods/pair_123902_205335_lr_thesis_cldn4/scripts/download.py
python3 methods/pair_123902_205335_lr_thesis_cldn4/scripts/analyze.py
```
