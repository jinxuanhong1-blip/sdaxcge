# Concordant-4 max-effect barrier communication

Patient-level CellChat / CellPhoneDB / LIANA-family scores for F11R, NECTIN2, CDH1, and LGALS9 from CLDN4-high vs CLDN4-low malignant cells to T/NK.

The PR 616 CellChat family sum is +0.0037 because population-size scaling compresses the probability. This folder recomputes that Hill score and the other standard magnitudes, and keeps the largest absolute family delta that is positive in all four cohorts. IFN/recruit chemokines and HLA–CD8 are separate families.

```bash
bash methods/concordant4_ccc_maxeffect_cldn4/scripts/download.sh /tmp/concordant4_raw
python3 methods/concordant4_ccc_maxeffect_cldn4/scripts/prepare_131907.py
gcc -O3 -o /tmp/extract_131907 methods/concordant4_ccc_maxeffect_cldn4/scripts/extract_131907.c -lz
/tmp/extract_131907 \
  /tmp/concordant4_raw/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz \
  /tmp/concordant4_cache/gse131907/genes_request.txt \
  /tmp/concordant4_cache/gse131907/keep_idx.i32 \
  208506 \
  /tmp/concordant4_cache/gse131907
python3 methods/concordant4_ccc_maxeffect_cldn4/scripts/run_maxeffect.py
```

GSE205335 needs R with the Matrix package (`export_gse205335.R` peels the double-gzipped RDS).
Narrative: `FINDING.md`. Formulas: `METHODS.md`.
