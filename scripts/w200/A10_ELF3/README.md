# A10 ELF3 — public ChIP + lung co-expression vs TACSTD2 / CLDN4

Honest public-data test only. No motif scan. No private ChIP.

```bash
python3 scripts/w200/A10_ELF3/chip.py
python3 scripts/w200/A10_ELF3/download.py
python3 scripts/w200/A10_ELF3/analyze.py
```

Raw matrices stay in `--cache-dir` (default `/tmp/a10_elf3_data`). Outputs: `results/w200/A10_ELF3/`.
