# B5 wave 2 — lung-first then all-open ICI, CLDN4 and TACSTD2

Second-wave rework of user claim B5 (11-cohort ICI meta, CLDN4-high OR=0.42).

Public processed matrices only. No EGA / dbGaP FASTQ. No fabricated cohorts.

## Locked result (honest)

CLDN4 ESTIMATEScore-residual, median split, DerSimonian–Laird RE:

| Pool | Endpoint | k | n | RE OR [95% CI] | p |
|---|---|---:|---:|---|---:|
| lung first | DCB (PFS≥6 mo) | 3 | 139 | 0.94 [0.44–1.99] | 0.87 |
| lung first | RECIST ORR | 2 | 51 | 0.64 [0.14–2.96] | 0.57 |
| all open ICI | DCB | 8 | 510 | 0.95 [0.66–1.36] | 0.77 |
| all open ICI | RECIST ORR | 13 | 921 | 0.85 [0.63–1.14] | 0.27 |

User 0.42 / k=11 is **not reproduced**. IMvigor210 CLDN4 raw median ORR = 1.47 p=0.214 (opposite).

Full write-up: `WRITEUP.md`. Tables: `tables/`. Forests: `figures/`.

```bash
python3 -m pip install -r results/rework/B5_wave2/requirements.txt
python3 results/rework/B5_wave2/download.py
python3 results/rework/B5_wave2/analyze.py
```
