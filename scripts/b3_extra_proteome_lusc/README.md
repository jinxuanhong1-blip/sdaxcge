# B3 extra — CPTAC LUAD / remaining proteome + TCGA-LUSC

Public-only pipeline. See `methods/B3_extra_proteome_lusc.md`.

```bash
python3 scripts/b3_extra_proteome_lusc/download.py
python3 scripts/b3_extra_proteome_lusc/analyze.py
```

Caches under `data/b3_extra_proteome_lusc/` (gitignored).
LSCC protein is not downloaded.
