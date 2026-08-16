# E-MTAB-13526 extra n — malignant TACSTD2/CLDN4 vs T/NK

Additive only. Public Cvejic / De Zuani NSCLC scRNA atlas
([E-MTAB-13526](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13526);
De Zuani et al., *Nat Commun* 2024, PMID 38821935). **Not an ICI cohort.**
Treatment-naive tumor resections. GSE207422 and other user analyses are not
re-run.

## Why a subset

Author-annotated objects on ArrayExpress:

| File | Size |
|---|---:|
| `10X_Lung_Tumour_Annotated_v2.h5ad` | 58.7 GB |
| `10X_Lung_Healthy_Background_Annotated_v2.h5ad` | 45.5 GB |

Those files are not downloaded. Instead we use the deposited processed 10x
Cell Ranger matrices for **CD235a− tumor** lanes only (RBC-depleted, not
CD45-sorted): 15 lanes, 12 patients. CD45+/MDSC-sorted tumor lanes cannot
support a malignant-vs-T/NK contrast (no epithelium; T/NK fractions inflated
by FACS). Background and donor lanes are omitted.

Matrices are unfiltered Cell Ranger 3.1.0 raw feature-barcode outputs
(33,538 genes × 6,794,880 whitelist barcodes). Author QC is applied while
streaming: UMI 400–100,000, genes 180–6,000, mitochondrial fraction ≤ 20%.

## Reproduce

```bash
python3 methods/emtab13526_tacstd2/download.py
python3 methods/emtab13526_tacstd2/extract.py
python3 methods/emtab13526_tacstd2/analyze.py
```

Raw matrices stay in `data/emtab13526/raw/` (gitignored). Patient-level
tables, Spearman n/ρ/p, epithelial-restriction stats, and the figure are
written into this folder.

## Analysis

- Unit = **patient** (lanes from the same patient are pooled).
- Lineage = argmax of mean log1p marker scores (same panel as the GSE253013 extra).
- Malignant-like = epithelial with near-zero normal-lung markers
  (*SFTPA2*, *AGER*, *SCGB1A1*, *SCGB3A1*, *TPPP3*).
- T/NK fraction among all QC cells in CD235a− tumor lanes.
- Eligible: ≥10 malignant-like and ≥20 T/NK cells.
- Association: Spearman rank correlation.
- Epithelial restriction: paired Wilcoxon of %positive malignant-like vs T/NK.

No public MPR / R / ICI labels exist on this accession.
