# GSE179994 CLDN4-only — no usable processed matrix

**Cohort:** Liu / Zhang (PKU), *Nat Cancer* 2022, PMID 35121991. Public GEO **GSE179994**. LUAD biopsies on pembrolizumab + carboplatin + pemetrexed. Interactive viewer mentioned in the paper (`http://nsclcpd1.cancer-pku.cn/`) is not a downloadable UMI/h5/TISCH extract.

**Task:** Additive **CLDN4-only** (no dual-high, no TACSTD2 gate). Malignant (author label or EPCAM+ epithelial) CLDN4 mean and %pos vs same-patient T/NK fraction, CD8, and CXCL13+ if present. Patient is the unit. Q4 vs Q1 only if `n_patients ≥ 16`. Public processed UMI / h5 / TISCH extract only. Skip files >2 GB.

## Verdict

**STOP. No usable processed all-cell UMI / h5 / TISCH extract <2 GB.**

The intended contrast cannot be estimated. **Analysis n (patients) = 0.** Do not cite T-cell barcode counts as n.

| Quantity | Value | Meaning |
|---|---:|---|
| Analysis n (patients with malignant CLDN4 and a same-patient all-cell T/NK denominator) | **0** | Endpoint not estimable |
| Patients in author T-cell metadata | 36 | T-cell barcodes only; not analysis n |
| Samples in author T-cell metadata | 47 | T-cell barcodes only |
| T-cell barcodes in metadata | 150,849 | **Not n** |
| Q4 vs Q1 | not done | Requires analysis `n_patients ≥ 16` |

## What was checked

### GEO series supplementary (all files <2 GB)

| File | Size | Role |
|---|---:|---|
| `GSE179994_all.Tcell.rawCounts.rds.gz` | 421.0 MB | T-cell UMI RDS. Under the 2 GB cap, but **T cells only**. |
| `GSE179994_Tcell.metadata.tsv.gz` | 0.89 MB | Author `celltype` / `cluster` / patient / sample. |
| `GSE179994_RAW.tar` | 44.1 MB | Two files: one blood T-cell Seurat RDS (`GSM5444629_P019.blood.post.Tcell.Seurat.object.rds.gz`) and its TCR TSV. Not an all-cell matrix. |
| `GSE179994_all.scTCR.tsv.gz` | 6.5 MB | scTCR. No expression. |
| `GSE179994_PBMC.bulkTCR.tsv.gz` | 2.6 MB | Bulk TCR. No expression. |

No all-cell UMI, MTX, h5, or h5ad is on the series record. Series text: “Raw data not provided for this record. Processed data are available on Series record.”

### Why the 421 MB T-cell RDS is not usable

Author metadata (`GSE179994_Tcell.metadata.tsv.gz`) covers 150,849 barcodes from 36 patients / 47 samples. `celltype` has exactly three values: **CD8** (59,656), **CD4** (58,706), **NA** (32,487). Clusters are T-cell states (Non-exhausted, Tex, Treg, naive/memory CD4, Prolif., XCL1). There is **no** epithelial, malignant, tumor, or EPCAM label.

Without a malignant / EPCAM+ epithelial compartment:

- CLDN4 mean and %pos in malignant cells cannot be computed.
- Same-patient T/NK **fraction** is undefined: the public matrix has no non-T denominator (myeloid, B, stroma, epithelium). A T-cell-only matrix would give T/NK ≈ 1 by construction.
- CD8 and CXCL13+ among T cells could be counted, but they are not paired to malignant CLDN4.

The RDS was **not** downloaded. The filename, GEO description, and full author metadata already establish T-cell-only content. Downloading 421 MB of T cells would not create a malignant compartment.

### Per-sample GEO records

Example `GSM5444604` (and the series sample block): “Supplementary data files not provided.” “Raw data not provided for this record.” Processed data sit only on the series record (T-cell files above).

### TISCH extract

HTTP HEAD of the usual TISCH2 paths:

- `NSCLC_GSE179994` expression.h5 and CellMetainfo → **404**
- `LUAD_GSE179994` expression.h5 and CellMetainfo → **404**
- `NSCLC_GSE179994_aPD1` / `LUAD_GSE179994_aPD1` → **404**

GSE179994 is not listed in the TISCH NSCLC or LUAD galleries. A prior TISCH NSCLC ICI pool in this repo (17 catalog datasets) also does not include GSE179994.

### Other public processed extracts (not used)

- CELLxGENE has no GSE179994 collection (404 on the accession path).
- The PKU web server is a viewer, not a public UMI/h5/TISCH download under the allowed input rule.
- Third-party reprocessings (e.g. CellResponse RDS) are outside “public processed UMI/h5/TISCH extract only.”

## Figures

- `figures/fig_geo_file_inventory.png` — GEO supplementary sizes vs the 2 GB skip line. Largest expression file is the T-cell RDS.
- `figures/fig_author_celltypes_t_only.png` — author `celltype` is CD4 / CD8 / NA only.
- `figures/fig_tcell_metadata_patient_census.png` — barcode census by patient in T-cell metadata. Caption states analysis n for malignant CLDN4 = 0.

## Reproduce

```bash
python3 methods/gse179994_cldn4/audit.py
```

Writes `tables/feasibility.json`, inventory TSVs, and the figures above. Downloads only the 0.89 MB metadata and 44 MB RAW.tar (listing). Does not fetch the T-cell RDS.

## Note

GSE179994 remains usable for **T-cell / TCR** questions (clonal Tex, pre/on-treatment). It is not usable for tumor-cell-intrinsic CLDN4 vs same-patient T/NK composition from a public processed all-cell matrix under 2 GB.
