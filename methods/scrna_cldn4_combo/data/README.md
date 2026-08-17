# Patient-level extracts (no new 9 GB download)

All tables are copies of processed patient scores already computed in this
repository. GSE253013 is the existing 9-patient tumor extract — the 9 GB
GEO RDS is not re-downloaded.

| Path | Source PR | Used for |
|---|---|---|
| `tnk/` | #279 `methods/scrna_meta_tnk` | CLDN4 vs T/NK (and CD8 / B where the table has it) |
| `tls/patient_scores.tsv` | #274 | CLDN4 vs B fraction and CXCL13+ T fraction |
| `exh/patient_scores.tsv` | #260 | CLDN4 vs T/NK cytotoxicity (cyto-high T score) |
| `leftover/` | #280 | Extra n (GSE267108, GSE274595, E-MTAB-13526) |
| `defgrid/grid_spearman.tsv` | #271 | TACSTD2-primary GSE207422 grid, **taken as given** (not re-audited) |
| `given_a3.json` | #279 | Locked A3 TACSTD2 numbers |

TACSTD2-primary combinations in PR #279 / #271 / #274 are given and are not
re-cut here. This slice is CLDN4-first.
