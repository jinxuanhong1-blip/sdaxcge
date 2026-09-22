# Concordant-4 max-effect TACSTD2 (TROP2+) barrier communication

Patient-level CellChat / CellPhoneDB / LIANA-family scores for F11R, NECTIN2, CDH1, and LGALS9 from **TACSTD2-high vs TACSTD2-low** malignant cells to T/NK.

Same score menu as the CLDN4 max-effect folder (PR 713). Gate gene is TACSTD2 for the Part1 PPT **TROP2+ barrier face**. Does not replace the locked CLDN4 table.

```bash
bash methods/concordant4_ccc_maxeffect_tacstd2/scripts/download.sh /tmp/concordant4_raw
python3 methods/concordant4_ccc_maxeffect_tacstd2/scripts/prepare_131907.py
gcc -O3 -o /tmp/extract_131907 methods/concordant4_ccc_maxeffect_tacstd2/scripts/extract_131907.c -lz
/tmp/extract_131907 \
  /tmp/concordant4_raw/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz \
  /tmp/concordant4_cache/gse131907/genes_request.txt \
  /tmp/concordant4_cache/gse131907/keep_idx.i32 \
  208506 \
  /tmp/concordant4_cache/gse131907
python3 methods/concordant4_ccc_maxeffect_tacstd2/scripts/run_maxeffect.py
```

GSE205335 needs R with the Matrix package (`export_gse205335.R` peels the double-gzipped RDS).
Narrative: `FINDING.md`. Formulas: `METHODS.md`. PPT one-liner: `results/PPT_SLIDE.md` (written after the run).
