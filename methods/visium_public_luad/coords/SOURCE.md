# Visium barcode → array coordinate maps

These CSVs are the public 10x Genomics barcode inclusion lists (barcode, `array_row`, `array_col` only). They are **not** section-specific `tissue_positions` from GSE277206 or Zenodo 13337961 (those depositions did not include spatial folders).

| File | Chemistry | Public source used to extract the list |
|---|---|---|
| `visium_v5_cytassist_11mm.csv` | CytAssist 11 mm (14,336 spots) | 10x `CytAssist_11mm_FFPE_Human_Lung_Cancer_spatial.tar.gz` (spatial-exp/2.0.1) |
| `visium_v4_cytassist_6p5mm.csv` | CytAssist 6.5 mm (4,992 spots) | 10x `CytAssist_FFPE_Human_Breast_Cancer_spatial.tar.gz` (spatial-exp/2.0.0) |

Barcode sequences on a given slide chemistry are fixed. Pixel columns from those demo slides were discarded because they belong to the demo H&E, not the LUAD sections.
