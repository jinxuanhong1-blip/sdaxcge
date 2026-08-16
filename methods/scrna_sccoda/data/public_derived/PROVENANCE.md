# Public-derived tables (not fabricated)

These TSVs are **re-derived from public GEO** by sibling open analyses in this
repository. They are bundled so GSE253013 (9.3 GB RDS) and GSE207422 DRMref
labels can enter the combinatorial grid without re-downloading over-budget or
third-party Seurat objects.

| File | Source accession | What it is | What it is not |
|---|---|---|---|
| `GSE207422_drmref_celltype_counts.tsv` | GSE207422 + DRMref public RDS | DRMref 16-type counts for 12 post-treatment samples | Hu/Zhang CopyKAT malignant IDs |
| `GSE207422_drmref_malig_tacstd2.tsv` | same | sample-level TACSTD2 in DRMref "Malignant cells" | author Seurat object |
| `GSE253013_per_patient_metrics.tsv` | GSE253013 Garnett RDS (streamed) | patient-level marker composition + TACSTD2 | a full cell×gene matrix |
| `GSE253013_author_label_per_patient.tsv` | GSE253013 `cell_type` slot | author epithelial / T counts only | a complete annotation |
| `GSE253013_sample_metadata.tsv` | GEO series matrix | GSM lanes, tumor vs ANT | MPR / ICI labels (none public) |

GSE207422 Hu-marker composition and GSE241934 author composition are
**recomputed in this folder** from the GEO matrices downloaded by
`scripts/00_download.py`.
