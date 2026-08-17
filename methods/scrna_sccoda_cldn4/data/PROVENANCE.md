# Provenance (already-extracted public tables only)

No new GEO matrices were downloaded. CLDN4 is the only exposure gene.
TACSTD2 columns that happen to sit on the source tables are ignored.

| File | Source branch / PR | Used for |
|---|---|---|
| `GSE207422_lineage_counts.tsv` | PR #318 `cursor/gse207422-cldn4-dualhigh-951b` | Hu marker T/NK/B counts |
| `GSE207422_per_patient.tsv` | PR #318 | A3-malignant CLDN4, timing, MPR |
| `GSE207422_drmref_patients.tsv` | PR #290 `cursor/scrna-cldn4-combo-a36d` | DRMref malignant CLDN4 (sensitivity) |
| `composition_long_pr283.tsv` | PR #283 `cursor/scrna-sccoda-tacstd2-53c7` | GSE241934 author_major + GSE207422 drmref_collapsed counts |
| `composition_all_pr284.tsv` | PR #284 `cursor/scrna-sccoda-composition-30fb` | GSE291670 6-part counts |
| `GSE241934_IIT_patients.tsv` | PR #290 | IIT epithelial CLDN4 |
| `GSE241934_Real_patients.tsv` | PR #290 | RWC/REAL epithelial CLDN4 |
| `GSE291670_patients.tsv` | PR #290 | marker-malignant CLDN4 |
| `GSE205335_patients.tsv` | PR #290 | author-malignant CLDN4; T/NK; B+plasma |
| `GSE253013_patients.tsv` | PR #290 | tumor malignant-like CLDN4 + T/NK |
| `tls_patient_scores.tsv` | PR #290 | GSE253013 B counts (marker_coarse TLS was 0) |

## Included (public lung ICI scRNA already extracted)

- GSE207422 Hu 2023 neoadjuvant PD-1
- GSE241934 NEOTIDE IIT + EGFR-WT real-world (RWC/REAL)
- GSE291670 neoadjuvant camrelizumab
- GSE205335 ICI with RECIST
- GSE253013 9-patient tumor extract (9.3 GB GEO RDS not downloaded)

## Explicitly not included

- GSE131907 — treatment-naive LUAD atlas, not ICI
- GSE325414 — already scored for Spearman elsewhere; not used here (B counts not in the extract)
- GSE267108 / GSE274595 leftover — not ICI; malignant CLDN4 empty for almost every sample
- Author CopyKAT barcodes for GSE207422 — not public on GEO
