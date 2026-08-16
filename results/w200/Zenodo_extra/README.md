# Zenodo/figshare extra: lung ICI × TACSTD2/CLDN4

## Bottom line

The leftover open, compact, processed lung-ICI expression cohort is
[Zenodo 2635194](https://doi.org/10.5281/zenodo.2635194): advanced NSCLC
treated with anti-PD-1, split into 55 discovery and 36 validation samples.
Its deposited NanoString panel measures 784 features including controls, but
it measures **neither TACSTD2 nor CLDN4**. Therefore this cohort cannot support
a TACSTD2/CLDN4 response comparison, correlation, effect size, or p-value.

No surrogate was substituted. EPCAM, CDH1, MUC1, and CEACAM6 are present, but
they are not interchangeable with TACSTD2 or CLDN4.

## What is reported

- `result.json`: machine-readable no-target result.
- `target_coverage.csv`: exact-symbol coverage in both deposited matrices.
- `cohort_summary.csv`: sample, response, and histology counts only.
- `repository_audit.csv`: dispositions for relevant Zenodo/figshare leftovers
  and explicit exclusions of the two PR42 records.
- `source_manifest.csv`: source-file sizes and verified repository MD5 hashes.
- `source_data/`: the four small CC BY 4.0 files from Zenodo 2635194.
- `analyze_zenodo_extra.py`: standard-library-only reproduction and integrity
  checks.

## Other exact lung-ICI leftovers

- [Zenodo 16916391](https://doi.org/10.5281/zenodo.16916391) is the open
  ACTS-30 neoadjuvant durvalumab plus chemoradiotherapy scRNA/spatial cohort.
  The target-capable total Seurat object is 23.34 GB inside a 29.15 GB archive;
  the repository offers no gene-level table or smaller epithelial object.
  It was not downloaded, and no target result is claimed.
- [Zenodo 15495774](https://doi.org/10.5281/zenodo.15495774) is a processed
  peripheral-blood anti-PD1 CyTOF/IMC cohort. Its antibody panel does not
  measure TACSTD2/TROP2 or CLDN4.
- [figshare collection 6590041](https://doi.org/10.6084/m9.figshare.c.6590041)
  accompanies a neoadjuvant PD-1 plus chemotherapy scRNA study, but its four
  deposited supplements are clinical/serum tables and a PDF—not a processed
  expression matrix.

The complete record-level audit also distinguishes restricted data, wrong
compartments/species, non-ICI atlases, and records already handled in PR42.

## Reproduce

From the repository root:

```bash
python3 results/w200/Zenodo_extra/analyze_zenodo_extra.py
```

The script fails if a source checksum changes, if expression and metadata row
counts differ, or if their row identifiers do not align.

## Interpretation limits

This is a coverage audit and an honest negative result, not evidence that
TACSTD2 or CLDN4 lack association with ICI outcome. Absence from a targeted
assay is missing measurement, not zero expression. The discovery/validation
response counts are descriptive and were not reinterpreted as a target
analysis.

## Source and license

Wiesweg et al., “RNA datasets to derive predictors for immune checkpoint
inhibitor therapy of non-small cell lung cancer,”
[Zenodo 2635194](https://doi.org/10.5281/zenodo.2635194), CC BY 4.0. The
corresponding publication is DOI
[10.1093/annonc/mdz049](https://doi.org/10.1093/annonc/mdz049).
