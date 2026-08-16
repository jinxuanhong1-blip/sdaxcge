# Zenodo/figshare extra: leftover lung ICI × TACSTD2/CLDN4

## Bottom line

No leftover open, compact, processed lung-ICI matrix supports a TACSTD2 or
CLDN4 comparison against ICI outcome.

1. [Zenodo 2635194](https://doi.org/10.5281/zenodo.2635194) is the leftover
   compact NSCLC anti-PD-1 expression cohort (55 discovery, 36 validation).
   Its NanoString panel has **neither TACSTD2 nor CLDN4**.
2. [Zenodo 11198494](https://doi.org/10.5281/zenodo.11198494) is a leftover
   processed lung PD-1 spatial record missed by PR42. The compact GeoMx table
   measures **TACSTD2** in 285 ROIs and **does not measure CLDN4**. TACSTD2 is
   above LOQ in only 28/285 ROIs, and 18 of those 28 detections come from one
   slide. The GeoMx table has **no ICI response labels**; the paper used
   immunotherapy-naive NSCLC to describe spatial niches.

No surrogate gene, below-LOQ quantitative test, or ICI-response p-value was
substituted.

## What is reported

- `result.json`: combined no-valid-ICI-target result.
- `target_coverage.csv` and `cohort_summary.csv`: NanoString gene and response
  counts only.
- `geomx_tacstd2_rois.csv`: extracted TACSTD2 rows from the leftover GeoMx
  table (original `geomx.csv` MD5 `8421f72c459e8aa6eb06fa48ba974501`).
- `geomx_detection_by_roi_type.csv` and `geomx_detection_by_slide.csv`:
  above-LOQ counts. `q_norm` medians are descriptive only.
- `geomx_tacstd2_above_loq.svg`: detection counts by ROI type and slide.
- `repository_audit.csv`: leftover dispositions, including restricted TROP2
  H-score data and records already in PR42.
- `source_data/`: four small CC BY 4.0 NanoString files from Zenodo 2635194.
- `analyze_zenodo_extra.py`: standard-library reproduction.

The 111 MB GeoMx source table is not committed.

## Other leftovers checked in this update

- [Zenodo 18543127](https://doi.org/10.5281/zenodo.18543127) / 18494664:
  resected-NSCLC TROP-2 H-scores exist, but the files are restricted and the
  cohort is not an ICI-outcome series.
- [Zenodo 16916391](https://doi.org/10.5281/zenodo.16916391): open ACTS-30
  durvalumab plus chemoradiotherapy objects; the total Seurat file is 23.34 GB
  and was not downloaded.
- RNAscope files in 11198494 (`ACD_stemimmunity.csv`) contain EPCAM and immune
  genes, not TACSTD2 or CLDN4.

## Reproduce

```bash
python3 results/w200/Zenodo_extra/analyze_zenodo_extra.py
```

The script fails if a NanoString checksum changes, if expression and metadata
rows do not align, or if the extracted GeoMx TACSTD2 table is not 285 ROIs.

## Interpretation limits

Missing assay features are not zero expression. Below-LOQ GeoMx values are not
quantitative TACSTD2 expression. Slide-concentrated detections are not a
spatial-niche ICI effect. This audit does not test whether TACSTD2 or CLDN4
predict immunotherapy outcome.

## Sources

- Wiesweg et al., Zenodo 2635194, CC BY 4.0; paper
  [10.1093/annonc/mdz049](https://doi.org/10.1093/annonc/mdz049).
- Chen, Nieman, Spurrell et al., Zenodo 11198494, CC BY 4.0; paper
  [10.1038/s41590-024-01792-2](https://doi.org/10.1038/s41590-024-01792-2).
