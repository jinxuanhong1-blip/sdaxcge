# results/hunt_sclc_ici

Public-data hunt: **TACSTD2 (TROP2) / CLDN4 vs SCLC-I (inflamed)**.

Start at **[REPORT.md](REPORT.md)**. It is written to be honest about what is and is not public.

| Path | What |
|---|---|
| `REPORT.md` | Verdict + methods + numbers |
| `IMPOWER133_NOTE.md` | Why IMpower133 RNA-seq is not used (EGA) |
| `tables/` | All numeric outputs (George, Chan, IMpower133 leftovers) |
| `figures/` | George + Chan plots |
| `data/download_data.sh` | Re-fetch public raw files (not committed; large) |
| `../../src/hunt_sclc_ici/` | Code |

```bash
bash results/hunt_sclc_ici/data/download_data.sh
python3 -m src.hunt_sclc_ici.run_hunt
```

`impower133_template.py` will **refuse** to run without real EGA files. That is intentional.
