# Existing scRNA patient-level malignant scores (no new download)

Copies of processed scores already computed in this repository. GSE253013
is the existing 9-patient tumor extract — the 9 GB GEO RDS is not
re-downloaded. GSE207422 A3 TACSTD2 is taken as given; CLDN4 is scored
from the same locked 12-patient DRMref table.

| Path | Source PR | Used for |
|---|---|---|
| `tnk/` | #279 `methods/scrna_meta_tnk` | Malignant CLDN4 vs T/NK (and B/CD8 where the table has it) |
| `tls/patient_scores.tsv` | #274 | Malignant-compartment CLDN4 vs B and CXCL13+ T |
| `given_a3.json` | #279 | Locked A3 TACSTD2 numbers (not re-audited) |

**Dropped here:** all-epithelial rows (`author_epi`, `marker_epi`, residual
Epi in GSE241934, leftover epi extracts). TLS rows with
`compartment != malignant` (epithelial / gate / insufficient) are dropped.
No dual-high TACSTD2×CLDN4 score is built.
