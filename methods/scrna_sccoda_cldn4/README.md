# scCODA-style composition: malignant CLDN4-high vs low

Additive public-only analysis. **CLDN4 only.** Patient is the unit.

Question: do **T/NK** or **B** fractions drop in patients with high malignant CLDN4
on GSE207422 and already-extracted public lung ICI scRNA?

Write-up: [`FINDING.md`](FINDING.md). Composition table:
[`results/composition_table.tsv`](results/composition_table.tsv).

```bash
python3 scripts/run_all.py
```

scCODA HMC is not faked if `sccoda` / tensorflow are missing. The engine is
ALR + patient permutation p, with Dirichlet-multinomial LRT as a diagnostic.
