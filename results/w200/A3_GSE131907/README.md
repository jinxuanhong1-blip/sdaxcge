# W200-A3 · GSE131907 — epithelial TACSTD2 vs T/NK

**Task:** A3 analog (malignant/epithelial TACSTD2 vs immune) on GSE131907, using processed data only.

## Verdict

**A3-style epithelial TACSTD2 vs T/NK anti-correlation is not supported** in this atlas
(tumor-site Spearman ρ = 0.17, p = 0.33, n = 36; tLung ρ = 0.09, p = 0.79, n = 11).
TACSTD2 **is** epithelial-restricted (median 76.1% vs 2.9% T/NK positive; Wilcoxon p = 1.4e-7).

The A3 **NMPR > MPR** half **cannot** be tested: GSE131907 is treatment-naive and has no ICI /
pathologic-response labels. Processed UMI (0.38 GB) was in budget. Raw EGA FASTQ and the
2.86 GB log2TPM text matrix were **not** downloaded.

## What was used

- Kim et al., *Nat Commun* 2020 (PMID 32385277); GEO [GSE131907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907).
- 208,506 cells / 58 samples / 44 LUAD patients; author `Cell_type` and `Cell_subtype`.
- Expression: processed raw UMI matrix (0.38 GB gzip). Selected genes streamed; other rows not parsed.
- Epithelial = `Cell_type == Epithelial cells`. In primary tLung these are tS1/tS2/tS3; author `Malignant cells` appear in metastases / PE / tL-B only.
- T/NK = `T lymphocytes` + `NK cells`.
- Eligible sample: ≥20 epithelial and ≥20 T/NK cells.
- Genes missing from the UMI matrix: none.

## Honest limits vs claim A3

- A3 source claim is GSE207422 (neoadjuvant ICI, MPR vs NMPR). GSE131907 is a treatment-naive LUAD atlas.
- No MPR / RECIST / PD-1 labels in the GEO series matrix.
- Primary tLung n=11 samples (eligible n=11); underpowered for a ρ≈−0.45 claim.
- Metric is mean log1p(raw UMI), not author log2(TPM+1), because the 2.86 GB normalized text file was skipped.

## Results

- tumor_sites: epi TACSTD2 mean_log1p vs T/NK fraction: ρ=0.168, p=0.326, n=36
- tumor_sites: epi TACSTD2 %pos vs T/NK fraction: ρ=-0.009, p=0.96, n=36
- tLung: epi TACSTD2 mean_log1p vs T/NK fraction: ρ=0.091, p=0.79, n=11
- metastasis/PE/tL-B: malignant TACSTD2 mean_log1p vs T/NK fraction: ρ=0.082, p=0.724, n=21
- paired tumor samples: epi TACSTD2 > T/NK TACSTD2 (Wilcoxon signed-rank): n=36, p=1.35e-07; median %pos epithelial=76.1 vs T/NK=2.9.

## Files

| File | Role |
|---|---|
| `feasibility.json` / `file_manifest.tsv` | Size budget and skip decisions |
| `sample_metadata.tsv` | GEO patient / stage / origin |
| `cell_type_composition.tsv` / `sample_composition.tsv` | Author labels |
| `per_sample_tacstd2.tsv` | Sample-level epithelial and T/NK TACSTD2 |
| `association_statistics.tsv` | Spearman tests |
| `summary.json` / `paired_compartment_test.json` / `sanity_checks.json` / `audit.json` | Verdict + extra tests |
| `fig_epi_tacstd2_vs_tnk_fraction.png` | A3-style correlation panels |
| `fig_tacstd2_epithelial_vs_tnk.png` | Tumor-restricted TACSTD2 |

## Reproduce

```bash
python3 scripts/w200/A3_GSE131907/download.py
python3 scripts/w200/A3_GSE131907/analyze.py
```
