# Concordant-4 scVI / scANVI (patient batch)

CLDN4-only. Four public cohorts already known to point the same way:

GSE123902 + GSE131907 + GSE205335 + GSE189357.

Not GSE148071, GSE127465, GSE154826, GSE207422, or GSE200563.

The locked unit stays the donor (GSE123902), tumor-bearing sample (GSE131907), or patient (GSE205335, GSE189357). n = 65. The integration cell count is not n.

```bash
python3 methods/scanvi_concordant4_cldn4/download.py --out /tmp/geo_c4
python3 methods/scanvi_concordant4_cldn4/prepare.py --raw /tmp/geo_c4
python3 methods/scanvi_concordant4_cldn4/integrate_model.py
```

`prepare.py` calls `extract_gse205335.R` (R + Matrix; the UMI file is double-gzipped). The GSE131907 streamer is `src/extract_cols.c` (`gcc -O3 -lz`). Mixed models are `mixed_model.R` (lme4). scVI batch key is the patient unit. scANVI annotates malignant cells on the two cohorts that have no author labels.
