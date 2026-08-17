# bfrac32_cldn4_hiend

ADDITIVE CLDN4-only high-end ligand–receptor on the **given** B-fraction combo
(GSE131907 + GSE241934 IIT, n=32, ρ=−0.513; PR #290). Do not re-audit that Spearman.

Patient is the unit. No dual-high. No GSE207422-only B-frac redo.

Write-up: [`FINDING.md`](FINDING.md). Primary table: [`results/lr_table.tsv`](results/lr_table.tsv).

```bash
python3 methods/bfrac32_cldn4_hiend/scripts/download.py --outdir /tmp/bfrac32_cldn4_hiend
python3 methods/bfrac32_cldn4_hiend/scripts/analyze.py --data /tmp/bfrac32_cldn4_hiend --out methods/bfrac32_cldn4_hiend
```
