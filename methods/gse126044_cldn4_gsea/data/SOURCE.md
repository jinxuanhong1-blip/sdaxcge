# GSE126044 source

Public processed counts only. No FASTQ / SRA.

| Item | Value |
|---|---|
| GEO | [GSE126044](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126044) |
| Title | Genome-wide identification of differentially methylated promoters and enhancers associated with response to anti-PD-1 therapy in non-small cell lung cancer |
| Paper | Cho et al., *Exp Mol Med* 2020, PMID 32879421 |
| File | `GSE126044_counts.txt.gz` from NCBI GEO FTP `suppl/` |
| URL | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/suppl/GSE126044_counts.txt.gz |
| Matrix | integer gene × sample counts, HGNC symbols, 18,747 rows × 16 columns |
| Samples | 16 pre-treatment NSCLC biopsies (11 fresh, 5 FFPE), anti-PD-1 |
| Labels | GEO `patient_response`: 5 responder / 11 non-responder |

Clinical map (`GSE126044_clinical.csv`) is transcribed from the GEO series matrix (`source_name_ch1`, `characteristics_ch1`). All 5 FFPE samples are non-responders.

CLDN4 vs ESTIMATE ImmuneScore (Spearman r = −0.524, p = 0.037, n = 16) is already in `methods/bulk_immune/results/demo/GSE126044_correlation_spearman.tsv` and is **not** re-tested here.
