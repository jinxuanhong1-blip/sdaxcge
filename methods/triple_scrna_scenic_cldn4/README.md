# Triple-merge CLDN4-only AUCell (GSE131907 + GSE148071 + GSE205335)

ADDITIVE. Malignant cells only. **CLDN4 only** — no TACSTD2 gate, no dual-high.
A10 ELF3–CLDN4 is given. This is **not** a redo of the 131907+205335-only SCENIC slice.

Patient-level question: do CLDN4-high tumors have different IFN / MHC / TJ
regulon activity than CLDN4-low tumors?

```bash
python3 methods/triple_scrna_scenic_cldn4/download.py
python3 methods/triple_scrna_scenic_cldn4/analyze.py
```

Raw matrices stay under `/tmp` and are not committed.

Primary table: `tables/regulon_table.tsv`. Write-up: `FINDING.md`.
