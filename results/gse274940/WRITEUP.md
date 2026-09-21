# GSE274940: EpH4 Cldn-null is not a CLDN4-only knockout, and Cldn4 mRNA is not depleted

Public bulk RNA-seq, mouse mammary EpH4, day 7. Three wild-type libraries versus three Cldn-null libraries (GSM8462258–GSM8462263). Counts are the GEO file `GSE274940_raw_counts.csv.gz` (RSEM expected counts, mm10). This series is the RNA-seq source data for Kashihara et al., Science Advances 2025 (PMID 41171911, doi:10.1126/sciadv.adx7431).

**Cldn4 mRNA is not depleted.** Mean raw counts are 15,074 (WT) and 16,733 (Cldn-null), 111% of the WT mean. After library-size normalization the Cldn-null mean is 128% of WT CPM (log2 fold-change +0.32, 95% Welch CI −0.65 to +1.30, p = 0.37, BH q = 0.64, n = 3 vs 3). The interval excludes a halving. Every null library still has Cldn4 in the same range as wild type.

**This is not a CLDN4-only knockout.** The cell line is a multiplex CRISPR null for the eight claudins EpH4 expresses: Cldn3, Cldn4, Cldn7, Cldn8, Cldn9, Cldn12, Cldn23, and Cldn25. The IFN and tight-junction numbers below are that eight-gene contrast. They are not a Cldn4 knockdown result, and they do not go in a CLDN4-only column.

The authors report frameshift indels and loss of claudin protein on immunoblot (their Fig. 1). This reanalysis does not rescore the blot. A frameshift can remove protein while leaving mRNA, which is what the deposited counts show for Cldn4.

## Targeted claudins

Calls use the pre-specified eight-gene list. `not_depleted`: wild-type mean raw count ≥ 10 and CPM remaining ≥ 80%. `depleted`: remaining < 50% and Welch p < 0.05. `reduced`: remaining < 80% and Welch p < 0.05. Genome-wide BH q is also shown. No gene in the transcriptome has q < 0.05 (13,661 genes tested; smallest q = 0.073), so the q column is not a discovery list. The Cldn4 question is answered by the estimate and the confidence interval.

| Gene | Mean raw WT | Mean raw null | CPM remaining | log2FC (null − WT) | 95% CI | Welch p | BH q | Call |
|---|---:|---:|---:|---:|---|---:|---:|---|
| Cldn3 | 27,669 | 3,257 | 14% | −2.93 | −4.18 to −1.69 | 0.0067 | 0.19 | depleted |
| Cldn4 | 15,074 | 16,733 | 128% | +0.32 | −0.65 to +1.30 | 0.37 | 0.64 | not depleted |
| Cldn7 | 6,982 | 2,305 | 37% | −1.45 | −2.01 to −0.88 | 0.0021 | 0.14 | depleted |
| Cldn8 | 1,060 | 594 | 64% | −0.66 | −1.24 to −0.08 | 0.035 | 0.28 | reduced |
| Cldn9 | 2,347 | 1,319 | 63% | −0.68 | −1.40 to +0.03 | 0.056 | 0.32 | lower, not significant |
| Cldn12 | 1,540 | 729 | 53% | −0.92 | −1.52 to −0.33 | 0.013 | 0.22 | reduced |
| Cldn23 | 843 | 1,044 | 140% | +0.49 | +0.19 to +0.78 | 0.011 | 0.21 | not depleted |
| Cldn25 | 0 | 0 | — | — | — | — | — | undetectable in both |

Cldn4 expected counts, library by library: WT 12,380 / 14,085 / 18,757; Cldn-null 16,155 / 12,982 / 21,063. The middle null library (12,982) sits inside the wild-type range. CPM for those six libraries is 500 / 553 / 709 versus 765 / 518 / 966.

Cldn3 and Cldn7 mRNA do fall, in every library. Any downstream shift in this contrast is therefore a multi-claudin effect dominated, at the RNA level, by Cldn3 (about 14% of WT CPM) and Cldn7 (about 37%), not by Cldn4.

Cldn23 mRNA is higher, not lower. Cldn25 has zero expected counts in all six libraries, so this matrix cannot confirm that Cldn25 was expressed or that its mRNA was removed.

Cldn1, which was not a CRISPR target, rises from 0/0/0 to 15/85/32 counts (mean CPM 0 versus 1.9). That is detectable, and it is about 400-fold below Cldn4 in the null cells. No other non-targeted claudin is expressed. There is no compensatory claudin program on the scale of the endogenous Cldn3/Cldn4 RNAs.

## IFN and tight junctions

Same contrast, same n. log2FC is Cldn-null minus WT on log2(median-of-ratios normalized count + 1). The sample-level test is a mean z-score per library, then a Welch t-test across the three versus three libraries. That is the test that uses the biological n. The competitive test compares gene-level t statistics inside the set with the rest of the transcriptome; its n is a number of genes. Set q values are Benjamini-Hochberg across the sets that were tested.

| Set | Genes tested | Sample Δ (null − WT) | Sample p | Sample q | Median log2FC | Competitive p | Competitive q |
|---|---:|---:|---:|---:|---:|---:|---:|
| TJ scaffolds, claudins removed | 15 | −0.36 | 0.28 | 0.80 | −0.17 | 0.20 | 0.61 |
| Adherens junction | 5 | +0.15 | 0.80 | 0.80 | +0.12 | 0.75 | 0.75 |
| IFN-I / ISG | 26 | −0.23 | 0.68 | 0.80 | −0.21 | 0.36 | 0.70 |
| IFN-II | 17 | −0.63 | 0.28 | 0.80 | −0.24 | 0.020 | 0.12 |
| MHC-I / APM | 15 | −0.16 | 0.80 | 0.80 | −0.02 | 0.58 | 0.70 |
| Epithelial identity, KO claudins removed | 9 | −0.18 | 0.67 | 0.80 | +0.05 | 0.58 | 0.70 |

No set moves at the sample level. The IFN-II competitive p of 0.020 is a downward rank shift (median log2FC −0.24), it does not survive set-level BH (q = 0.12), and the three-versus-three test is p = 0.28. That is not an interferon opening.

Genes that are actually expressed stay near zero:

| Gene | log2FC | 95% CI | Welch p | BH q | Mean CPM, WT → null |
|---|---:|---|---:|---:|---|
| Ocln | +0.21 | −0.52 to +0.93 | 0.39 | 0.65 | 66 → 77 |
| Tjp1 | +0.07 | −0.55 to +0.69 | 0.77 | 0.89 | 69 → 72 |
| F11r | −0.58 | −1.08 to −0.08 | 0.033 | 0.28 | 366 → 246 |
| Cdh1 | −0.64 | −1.29 to +0.00 | 0.051 | 0.31 | 751 → 485 |
| Tacstd2 | −0.61 | −1.15 to −0.07 | 0.038 | 0.28 | 2,092 → 1,377 |
| Stat1 | −0.07 | −1.08 to +0.94 | 0.86 | 0.93 | 21 → 20 |
| Isg15 | −0.12 | −0.97 to +0.74 | 0.68 | 0.84 | 3.5 → 3.3 |
| Ifit1 | −0.15 | −2.74 to +2.44 | 0.87 | 0.94 | 6.0 → 5.0 |
| B2m | +0.19 | −0.88 to +1.26 | 0.62 | 0.80 | 108 → 127 |
| H2-K1 | +0.12 | −0.45 to +0.69 | 0.49 | 0.72 | 56 → 61 |
| Cd274 | −0.23 | −2.62 to +2.17 | 0.80 | 0.91 | 6.3 → 5.8 |

Nominal single-gene shifts that do not survive transcriptome-wide BH, reported because they sit in these sets:

- Ifitm3 is up (log2FC +0.94, CPM 20 → 38, p = 0.0019, q = 0.14). Ifitm2 is slightly up (log2FC +0.33, p = 0.020, q = 0.24).
- Cldnd1 mRNA falls (log2FC −4.57, CPM 97 → 4.1, p = 0.0012, q = 0.13). Cldnd1 was not one of the eight CRISPR targets. Marveld3 is modestly down (log2FC −0.53, p = 0.016, q = 0.23).
- Psmb8 is down (log2FC −4.99, CPM 28 → 1.0, p = 0.013, q = 0.22). The MHC-I set as a whole is flat, and B2m and H2-K1 do not move.

Ifng, Cxcl9, Cxcl10, and Cxcl11 are at or near zero counts in every library. Mx1 is essentially absent (wild-type counts 0/0/0). Checkpoint receptors other than a low Tigit signal are absent; Tigit itself is log2FC +1.06 with a CI that includes zero (p = 0.12). Ddx58 is not in the count matrix. This is an epithelial monolayer, not an immune co-culture, so those chemokines and receptors were not available to score.

Tacstd2 mRNA is nominally lower. Cldn4 mRNA is not lower in the same samples, and Cldn3/Cldn7 mRNA are, so that Tacstd2 shift cannot be read as an effect of losing Cldn4.

## What this series is not used for

- It is not a CLDN4-only knockdown or knockout transcriptome.
- It is not lung, and it is not an immune co-culture.
- It does not show that removing Cldn4 opens or closes IFN, MHC-I, or the rest of the tight junction. Cldn4 mRNA was not removed.
- Author immunoblots are the protein-level evidence. These tables are the RNA-level evidence.

## Methods

Expected counts were summed across Ensembl ids that share a gene symbol. Size factors are the DESeq2 median-of-ratios factors (WT 1.03, 1.07, 1.11; null 0.86, 1.08, 0.92). Differential expression is a Welch t-test on log2(normalized count + 1), genes with a group-mean raw count of at least 10 in the higher group (13,661 genes). The Mann-Whitney p-value is stored in the gene table and is not used: with n = 3 vs 3 its two-sided floor is 0.1. Gene-set membership is `scripts/gse274940/gene_sets_mouse.tsv`. Mouse IFN symbols are orthologs of the sets used in the earlier CLDN4/TACSTD2 knockdown wave (Oas1a, Oasl1/Oasl2, Herc6, Trim30a, H2-K1/H2-D1), not a title-case of the human symbols. Tight-junction scaffolds are scored with every Cldn gene removed, so the set is not the CRISPR deletion.

Counts: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274940/suppl/GSE274940_raw_counts.csv.gz
sha256 `ff1ebb53a6732bf1b053389be6bdd467d6f740601cb932e651b1f9af739b949f`

Reproduce with `python3 scripts/gse274940/analyze_gse274940.py`.
