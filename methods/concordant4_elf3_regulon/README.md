# Concordant-4 ELF3–TACSTD2–CLDN4–CLDN7

Public malignant-cell coexpression and DoRothEA regulon activity for ELF3,
GRHL2, TACSTD2, CLDN4, and CLDN7, plus the patient-level association with
T/NK fraction. Cohorts are only GSE123902, GSE131907, GSE205335, and
GSE189357 (n = 65 units).

```bash
bash methods/concordant4_elf3_regulon/scripts/download.sh /tmp/geo_c4
python3 methods/concordant4_elf3_regulon/scripts/analyze.py --geo /tmp/geo_c4
```

GSE205335 needs R with the Matrix package (`extract_gse205335.R`).
See `METHODS.md` and `FINDING.md`.
