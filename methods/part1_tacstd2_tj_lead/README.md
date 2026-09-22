# PART 1 polish — Tacstd2-high → TJ/junction lead ranks

Re-ranks frozen public tables from:

- PR #741 `methods/concordant4_tacstd2_malignant_deg/` (concordant-4)
- PR #736 `methods/paper_funnel_luad_tacstd2_tj/` (OncoSG + GSE31210)

Does **not** re-download matrices or invent NES/FDR. **No CLDN4 pin.**

```bash
pip install -r methods/part1_tacstd2_tj_lead/requirements.txt
python3 methods/part1_tacstd2_tj_lead/rank_from_provenance.py
```

See `FINDING.md` for the Part 1 ending and best TJ/junction ranks.
