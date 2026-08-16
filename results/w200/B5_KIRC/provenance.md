# Provenance — B5_KIRC

All inputs are public. No controlled-access (EGA/dbGaP) data were used.

## Downloads

| File | Source | Retrieved |
|---|---|---|
| `data/braun2020/braun2020_supp.xlsx` (125,569,016 bytes; sha256 `a1f8683676d416b255291dcaa007aa0fa48c736c23e6a7a566664ba528bf205b`) | https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41591-020-0839-y/MediaObjects/41591_2020_839_MOESM2_ESM.xlsx | 2026-08-16 |
| `data/javelin101/javelin_supp.xlsx` (125,767,390 bytes; sha256 `673c5fd66813e3267b3a514074ac22d090a258b1094b45119d592405d5f663ed`) | https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41591-020-1044-8/MediaObjects/41591_2020_1044_MOESM3_ESM.xlsx | 2026-08-16 |
| `data/gse67501/GSE67501_series_matrix.txt.gz` (1,029,286 bytes; sha256 `e48b9c3b01a160ff59233d185a17df1c1210ba244b637a0a85ea08aaf4b5bd2a`) | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE67nnn/GSE67501/matrix/GSE67501_series_matrix.txt.gz | 2026-08-16 |

Papers:

- Braun DA et al. *Nat Med* 2020;26:909–918. doi:10.1038/s41591-020-0839-y
- Motzer RJ et al. *Nat Med* 2020;26:1733–1741. doi:10.1038/s41591-020-1044-8
- Ascierto ML et al. *Cancer Immunol Res* 2016;4:726–733. GEO GSE67501. PMID 27491898

## Processing notes

- Braun expression used as distributed (ComBat-corrected, UQ-normalized log2 TPM
  + constant). No re-normalization.
- JAVELIN TPM used as distributed (log2 TPM, floor 0.01 per the paper).
- GSE67501: Illumina HumanHT-12 v4 probe **ILMN_2132458** (CLDN4), author-scaled
  log2 intensities from the series matrix.
- Braun `PFS_CNSR` / `OS_CNSR`: coded as 1 = event (BMS-style; 80%+ of PD cases
  have CNSR=1). JAVELIN `PFS_P_CNSR`: 1 = censored (CDISC; longer PFS when
  CNSR=1). Documented in `summary.json`.

## Not used (not open expression)

IMmotion150/151 RNA (EGA), CheckMate-214 RNA, HCRN GU16-260, any institutional
RCC ICI matrix.
