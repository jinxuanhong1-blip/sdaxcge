# Triple that differs — thesis-aligned CLDN4-only ligand test

ADDITIVE. **CLDN4 only.** No dual-high. No GSE148071. No 7-pool.

Given cut (PR #459, not re-audited): GSE123902+GSE131907+GSE205335 %pos n=56, ρ=−0.522, Q4 vs Q1 r=−0.735.

Patient/donor/locked-sample is the unit. Done when both family tables exist:

- `results/family_barrier_inhibitory.tsv`
- `results/family_ifn_recruit_mhci.tsv`

```bash
python3 methods/triple_differ_lr_thesis_cldn4/scripts/download.py
python3 methods/triple_differ_lr_thesis_cldn4/scripts/analyze.py
```
