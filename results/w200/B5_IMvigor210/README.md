# B5 — IMvigor210 CLDN4-high vs ORR / OS

Public processed RNA from **IMvigor210** (atezolizumab, metastatic urothelial carcinoma). Question: is **CLDN4-high** associated with objective response or overall survival?

**Honest headline: no. Primary tests are null.**

## OR / n / p (primary, prespecified)

| Endpoint | Contrast | Effect | 95% CI | n | p |
|---|---|---|---|---|---|
| **ORR** (CR/PR vs SD/PD) | CLDN4-high vs low (median) | **Fisher OR 1.465** | 0.821–2.637 | **298** (39/149 vs 29/149) | **0.214** |
| **OS** | CLDN4-high vs low (median) | **Cox HR 0.890** | 0.688–1.152 | **347** (231 deaths) | **0.377** |

2×2 ORR (NE excluded):

| CLDN4 | SD/PD | CR/PR | ORR |
|---|---|---|---|
| low | 120 | 29 | 19.5% |
| high | 110 | 39 | 26.2% |

Direction is a modest **higher** ORR and slightly **better** OS in CLDN4-high. Intervals include the null. Continuous CLDN4 (per SD) is also null: ORR OR 1.13 (n=298, p=0.40); OS HR 0.98 (n=347, p=0.75).

Do not read this as evidence that CLDN4-high predicts ICI benefit or resistance in urothelial cancer.

## What was opened

Processed counts + clinical annotation from the official Genentech package, not a guessed GEO series.

| Item | Value | How verified |
|---|---|---|
| Package | IMvigor210CoreBiologies **1.0.0** | Downloaded from `research-pub.gene.com/IMvigor210CoreBiologies/packageVersions/` |
| Object | `cds` CountDataSet | 31,286 genes × **348** RNA samples |
| SHA256 | `cfdd3176d7b34de5b04fb9416bfd2b20fa4b6e238aaad5f20b048a34329ea178` | `sha256sum` on the 122,127,298-byte tarball |
| Trials | NCT02108652, NCT02951767 | Package `cds` help |
| Paper | Mariathasan et al., *Nature* 2018 | DOI `10.1038/nature25501` |
| Raw RNA study | **EGAS00001002556** | Stated on the package landing page and in `?cds` `@source` |
| Raw RNA dataset | **EGAD00001003977** | Listed on the EGA study page (FFPE capture RNA, pre-anti-PD-L1) |
| Raw DNA dataset | **EGAD00001004218** | Listed on the EGA study page; **not used** |

**No GEO / GSE / SRA / PRJ accession is claimed.** None is printed on the official package page, the `cds` help, or the EGA study page. Raw EGA is controlled-access; this run did not download FASTQ/BAM.

CLDN4 is present once: Entrez **1364**, IGIS feature, length 1831 nt.

## Prespecified analysis

1. **Exposure.** `log2(DESeq size-factor-normalized CLDN4 count + 1)` using the package `sizeFactor`. High = ≥ median (8.869). Median is prespecified, not optimized.
2. **ORR.** Package `binaryResponse`: CR/PR vs SD/PD. `NE` dropped (49 patients after dedup). Fisher exact OR is primary; logistic OR is reported beside it.
3. **OS.** Package `os` (months) + `censOS` (1 = death). Unadjusted Cox HR + log-rank.
4. **Dedup.** 348 RNA samples, 347 `ANONPT_ID`. Patient 10285 has two bladder RNA aliquots, both `NE`, identical OS (0.49 mo), CLDN4 9.86 vs 9.61. First sample kept. ORR is unaffected.

Sensitivities (all still compatible with null): tertile T3 vs T1, Q4 vs Q1, median of `log2(CPM+1)` (Pearson r = 0.982 vs DESeq log), continuous per SD, multivariable complete-case (IC level, ECOG, sex, liver mets, prior platinum, log1p TMB; n=210).

## Exploratory subgroups (not primary)

IC1 is the only unadjusted p < 0.05: Fisher OR 5.59 (1.47–31.8), **17/63 vs 3/49**, n=112, p=0.0054. OS HR 0.68 (0.45–1.02), p=0.061.

That cell is small (3 responders in CLDN4-low IC1), the contrast was not prespecified, and ~15 subgroups were scanned. It does **not** overturn the primary null. Treat as a hypothesis, not a finding.

## Limitations (read these)

- Single-arm PD-L1 blockade. Even a real association would be prognostic, not predictive of atezolizumab vs chemo.
- FFPE capture RNA-seq, not poly(A). Counts are the published processed matrix, not a re-alignment from EGA.
- Median split is crude; the continuous model is the better-powered test and is closer to null.
- TMB missing in 76/347 patients, so multivariable n drops to 210.
- Immune phenotype missing in 64/347.
- No multiple-testing correction on subgroups.
- Package LICENSE is a PDF; this repo keeps derived sample-level scores + clinical fields, not the full 31k-gene matrix.

## Reproduce

```bash
# official tarball (HTTP; HTTPS from some networks resets)
curl -O http://research-pub.gene.com/IMvigor210CoreBiologies/packageVersions/IMvigor210CoreBiologies_1.0.0.tar.gz
sha256sum IMvigor210CoreBiologies_1.0.0.tar.gz
# expect cfdd3176d7b34de5b04fb9416bfd2b20fa4b6e238aaad5f20b048a34329ea178

# needs R with Biobase; scripts/01_extract_cds.R uses a tiny DESeq S4 stub
# if the real DESeq package is absent (Bioconductor DESeq is archived)
bash results/w200/B5_IMvigor210/scripts/run.sh IMvigor210CoreBiologies_1.0.0.tar.gz
```

## Files

- `tables/headline_OR_n_p.tsv` — OR / n / p
- `tables/orr_binary.tsv`, `orr_continuous.tsv`, `orr_2x2_primary.tsv`
- `tables/os_cox.tsv`, `multivariable.tsv`, `subgroups.tsv`, `recist4_by_cldn4.tsv`
- `figures/orr_bar_median.png`, `km_os_median.png`, `cldn4_by_binary_response.png`, `cldn4_by_recist4.png`, forests
- `data/sample_level.tsv` — one row per RNA sample (CLDN4 + clinical only)
- `data/analysis_cohort.tsv` — one row per patient after dedup + analysis flags
- `data/provenance.json` — accessions and checksum
- `scripts/01_extract_cds.R`, `02_analyze_cldn4.R`
