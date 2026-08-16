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
3. Both leftover RNAscope tables in 11198494 (`ACD_stemimmunity.csv` and
   `ACD_tumor.csv`) use the same 12-gene panel and measure **neither target**.
4. The closest new lung-IO spatial leftover,
   [Zenodo 21807253](https://doi.org/10.5281/zenodo.21807253) (neoadjuvant
   immunotherapy NSCLC spatial transcriptomics), is **restricted**.

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
- `rnascope_panel.csv`: verified 12-gene leftover RNAscope panel for both ACD
  files.
- `ici_module195_genes.csv`: 195 ICI-module genes from Zenodo 21967935; neither
  requested target is present.
- `repository_audit.csv`: leftover dispositions, including restricted TROP2
  H-score data, records already in PR42, and the Grok 4.6 hunt additions.
- `source_data/`: four small CC BY 4.0 NanoString files from Zenodo 2635194.
- `analyze_zenodo_extra.py`: standard-library reproduction.

The 111 MB GeoMx source table is not committed.

## Other leftovers checked in this update

- [Zenodo 21807253](https://doi.org/10.5281/zenodo.21807253): neoadjuvant-IO
  NSCLC spatial transcriptomics; **restricted**.
- [Zenodo 18825477](https://doi.org/10.5281/zenodo.18825477): processed ICI
  scRNA plus Visium; **restricted**.
- [Zenodo 15053161](https://doi.org/10.5281/zenodo.15053161): n=1 Visium after
  chemo-immunotherapy; **restricted**.
- [Zenodo 18870185](https://doi.org/10.5281/zenodo.18870185): NSCLC CosMx WTA;
  embargoed until 2029 and not an ICI-outcome series.
- [Zenodo 21967935](https://doi.org/10.5281/zenodo.21967935): open ICI-module
  supplementary zip; 195 genes, no expression matrix, no TACSTD2/CLDN4.
- [Zenodo 14511579](https://doi.org/10.5281/zenodo.14511579): `v4 ICB
  datasets.xlsx` lists skin, liver, breast, and kidney only; no lung row.
- [Zenodo 18543127](https://doi.org/10.5281/zenodo.18543127) / 18494664:
  resected-NSCLC TROP-2 H-scores exist, but the files are restricted and the
  cohort is not an ICI-outcome series.
- [Zenodo 16916391](https://doi.org/10.5281/zenodo.16916391): open ACTS-30
  durvalumab plus chemoradiotherapy objects; the total Seurat file is 23.34 GB
  and was not downloaded.
- Both RNAscope files in 11198494 share the 12-gene panel
  IDO1/CCR7/FOXP3/CD4/PDCD1/CCL22/EPCAM/CCL19/CD8A/CXCL10/CD3E/TCF7.

## Reproduce

```bash
python3 results/w200/Zenodo_extra/analyze_zenodo_extra.py
```

The script fails if a NanoString checksum changes, if expression and metadata
rows do not align, if the extracted GeoMx TACSTD2 table is not 285 ROIs, if
either leftover RNAscope panel is not the verified 12-gene list, or if the
195-gene ICI module table is not 195 unique symbols without TACSTD2/CLDN4.

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
- Goldman, Zenodo 21967935, CC BY 4.0; 195-gene ICI module dictionary only.
