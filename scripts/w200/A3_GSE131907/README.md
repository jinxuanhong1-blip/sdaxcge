# W200-A3 analog — GSE131907 epithelial TACSTD2 vs immune

A3 (user claim) is **GSE207422**: malignant-cell `TACSTD2` higher in NMPR than MPR, and
negatively correlated with T/NK presence (ρ ≈ −0.40 to −0.50).

This folder applies the **same contrast that GSE131907 can support**: author-labeled
epithelial/tumor-epithelial `TACSTD2` vs the T/NK compartment. GSE131907 is a treatment-naive
LUAD atlas (Kim et al., Nat Commun 2020; PMID 32385277). It has **no ICI / MPR labels**, so
the NMPR>MPR half of A3 is not testable here.

## Data policy

- Use the processed GEO UMI matrix (`GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz`, 0.38 GB).
- Download the cell annotation (1.8 MB) and series matrix (6 KB).
- **Skip** the 2.86 GB `normalized_log2TPM_matrix.txt.gz`.
- **Skip** EGA raw FASTQ (`EGAD00001005054`); GEO states raw data are not open.

```bash
python3 scripts/w200/A3_GSE131907/download.py
python3 scripts/w200/A3_GSE131907/analyze.py
```

Outputs land in `results/w200/A3_GSE131907/`. The UMI matrix is written under `/tmp/gse131907/`
and is not committed.
