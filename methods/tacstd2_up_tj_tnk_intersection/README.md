# TACSTD2-high UP ∩ TJ ∩ low T/NK

Story gap after PR #741. Among genes up in TACSTD2-high malignant cells
on the #741 list, which are tight-junction members and which of those
individually predict low T/NK?

```bash
python3 methods/tacstd2_up_tj_tnk_intersection/analyze.py
```

Writes `FINDING.md`, `tables/`, and `figures/`. Does not refit the locked
CLDN4 %pos vs T/NK result except as a reference row on the same units.
