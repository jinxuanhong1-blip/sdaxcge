# Scripts

| Step | What |
| --- | --- |
| `00_download.py` | GEO processed UMI / identity / SOFT + NicheNet-v2 ligand–target RDS |
| `01_convert_prior.py` | RDS → `/tmp/multinichenet_cldn4/prior_ligand_target.parquet` |
| `02_extract_panels.py` | Stream / extract the ligand–receptor–program panel; keep malignant + T/NK |
| `03_analyze.py` | Patient-level MultiNicheNet-style ranks, figures, `FINDING.md` |

Helpers: `gene_sets.py`, `lib_io.py`, `lib_nichenet.py`, `lib_stats.py`.
