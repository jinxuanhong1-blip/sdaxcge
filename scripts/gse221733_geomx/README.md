# GSE221733 GeoMx CTA RNA (CLDN4 panel check)

Public GeoMx DSP RNA from Monkman et al., *Immunology* 2023 (PMID 37022147).

```bash
bash scripts/gse221733_geomx/00_fetch.sh
python3 scripts/gse221733_geomx/analyze.py
```

CLDN4 is not on the deposited Cancer Transcriptome Atlas (CTA) panel. The analyzer inventories the panel, writes ROI metadata, and skips the CLDN4 vs CD8 / ICB-response tests. It does not substitute TACSTD2 or EPCAM for CLDN4.
