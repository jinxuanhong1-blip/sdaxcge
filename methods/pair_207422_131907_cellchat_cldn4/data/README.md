# Source tables (taken as given; not re-downloaded)

Patient-level malignant CLDN4 vs T/NK:

- `tnk/GSE207422_drmref_patients.tsv` — locked A3 12-patient DRMref table (PR #279 / `cldn4_malig_q4_tnk`).
- `tnk/GSE131907_samples.tsv` — Kim atlas sample extract (same PR). Author `Malignant cells` only.
- `given_a3.json` — A3 TACSTD2 slide is given; this pair uses **CLDN4 only**.

CellChat-style outgoing Mal → T/NK (not re-run):

- `cellchat/GSE207422_outgoing.tsv` — kept split `median_post` (`methods/scrna_cellchat_cldn4`).
- `cellchat/GSE131907_outgoing.tsv` — kept split tumor-origin tertile (`methods/gse131907_cellchat_cldn4`).
- Incoming contrasts are stored for the extra IFNG row only.

GEO UMI matrices are not stored here.
