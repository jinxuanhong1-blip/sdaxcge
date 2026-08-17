# Scripts

- `download.py` — public GEO processed files + CellPhoneDB v5 tables → `/tmp`.
- `analyze.py` — propeller table, CellPhoneDB-style focus table, figures, `FINDING.md`.

```bash
python3 methods/winpair_131907_205335_propeller_cpdb_cldn4/scripts/download.py
python3 methods/winpair_131907_205335_propeller_cpdb_cldn4/scripts/analyze.py
```

`--skip-cpdb` writes the propeller table only.
