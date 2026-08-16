# A4 — TISMO all 64 ICB slices, Tacstd2 paired sign test

Recompute the user claim (Tacstd2 up in 49/64 TISMO models, p=5.8e-5)
from the public TISMO in-vivo ICB expression API.

```bash
python3 scripts/w200/A4_all64/download.py
python3 scripts/w200/A4_all64/analyze.py
```

Outputs land in `results/w200/A4_all64/`.
The downloaded per-sample table is written to
`results/w200/A4_all64/tacstd2_vivo_icb.csv` so the analysis can be
re-run offline.
