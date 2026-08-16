# B5 analog: Miao 2018 RCC ICI — CLDN4 and TACSTD2 versus response

## Bottom line

Pretreatment tumor **CLDN4 and TACSTD2 do not separate ICI responders from
nonresponders in this public RNA subset**. Among 33 RNA-profiled tumors
(8 RECIST CR/PR, 25 SD/PD), the responder-higher AUCs were 0.510 for CLDN4
(bootstrap 95% CI 0.305–0.715; two-sided Mann–Whitney p=0.950) and 0.412 for
TACSTD2 (0.170–0.670; p=0.475). These are null results. They are not evidence
that either gene predicts response, and they are not evidence of a protective
or adverse effect.

The cohort is too small, too mixed, and too batch-split to support a biomarker
claim. A non-significant p-value here is lack of support, not proof of no
effect.

## Cohort and endpoint

- Paper: Miao et al., *Science* 2018;359:801–806
  ([PubMed 29301960](https://pubmed.ncbi.nlm.nih.gov/29301960/)).
- Public RNA+response table used here: BHK Lab / ORCESTRA `ICB_Miao1` TSV
  extract, Zenodo [10.5281/zenodo.7058399](https://doi.org/10.5281/zenodo.7058399),
  CC-BY-4.0. This is a re-annotated public extract, not a new alignment of the
  original FASTQs.
- Disease / treatment: kidney cancer treated with PD-1/PD-L1 blockade
  (nivolumab n=27, nivolumab+ipilimumab n=4, atezolizumab n=2).
- Analysis set: all 33 samples with RNA (`rna=tpm`). One sample per ID.
- Genes: **CLDN4** `ENSG00000189143.8` and **TACSTD2** `ENSG00000184292.5`,
  each a single protein-coding row.
- Expression scale: the source matrix is already log2(TPM) with a floor at
  ≈log2(0.001) (−9.96578428466209). No further transform was applied.
- Primary endpoint (B5 analog): RECIST objective response, CR/PR = responder
  versus SD/PD = nonresponder. This matches the Braun B5 rule and uses the
  public `recist` field, not a fitted cutpoint.

Secondary endpoints, all prespecified: BHK curated `response` R versus NR
(n=10/18; 5 SD samples are NA); Miao-style `response_category` clinical
benefit versus no clinical benefit, excluding intermediate (n=12/13); overall
survival Cox models. **PFS was not tested** (see data integrity).

## Primary result

| Gene | R / NR | Median R | Median NR | Δ median (R−NR) | 95% bootstrap CI | AUC (R higher) | 95% bootstrap CI | MW p |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| CLDN4 | 8 / 25 | 1.930 | 1.190 | +0.741 | −1.802 to +2.476 | 0.510 | 0.305–0.715 | 0.950 |
| TACSTD2 | 8 / 25 | −1.830 | −0.808 | −1.022 | −3.547 to +2.333 | 0.412 | 0.170–0.670 | 0.475 |

Medians are on the source log2(TPM) scale. CLDN4 was above the floor in 8/8
responders and 24/25 nonresponders; TACSTD2 in 8/8 and 24/25 (Fisher p=1.0
for both). Detection does not associate with response. The two genes are
moderately co-expressed (Spearman ρ=0.488, p=0.004), as expected for
epithelial / tight-junction transcripts, but that correlation is not a
response signal.

## Sensitivities (still null)

| Contrast | Gene | n | AUC (high group higher) | MW p |
|---|---|---:|---:|---:|
| RECIST ORR, ID prefix PD_ | CLDN4 | 5 / 12 | 0.567 | 0.712 |
| RECIST ORR, ID prefix PD_ | TACSTD2 | 5 / 12 | 0.233 | 0.102 |
| RECIST ORR, ID prefix RCC_ | CLDN4 | 3 / 13 | 0.385 | 0.590 |
| RECIST ORR, ID prefix RCC_ | TACSTD2 | 3 / 13 | 0.667 | 0.420 |
| RECIST ORR, nivolumab only | CLDN4 | 4 / 23 | 0.609 | 0.517 |
| RECIST ORR, nivolumab only | TACSTD2 | 4 / 23 | 0.576 | 0.657 |
| RECIST ORR, prefix-residualized | CLDN4 | 8 / 25 | 0.480 | 0.883 |
| RECIST ORR, prefix-residualized | TACSTD2 | 8 / 25 | 0.400 | 0.413 |
| Curated R vs NR, pooled | CLDN4 | 10 / 18 | 0.489 | 0.943 |
| Curated R vs NR, pooled | TACSTD2 | 10 / 18 | 0.439 | 0.615 |
| CB vs NCB, pooled | CLDN4 | 12 / 13 | 0.532 | 0.807 |
| CB vs NCB, pooled | TACSTD2 | 12 / 13 | 0.474 | 0.849 |

The lowest p-value in the table is TACSTD2 inside the PD_ prefix (p=0.102),
with the opposite point estimate in RCC_. That is what a batch-split n=3–5
responder analysis looks like. It is not a subgroup finding.

Exploratory OS (23 events / 33): CLDN4 HR 1.01 per log2 unit (95% CI
0.89–1.14, p=0.886); TACSTD2 HR 0.94 (0.76–1.15, p=0.528). Prefix-adjusted
Cox models are the same conclusion.

## Data integrity (why this cohort cannot carry a claim)

1. **Expression batch dominates biology.** On the 2,000 most variable genes
   detected in >50% of samples, PC1 explains 68.4% of variance and separates
   sample-ID prefixes PD_ (n=17) from RCC_ (n=16) at Mann–Whitney p=1.1×10⁻6.
   Prefix is also confounded with regimen: all 4 nivo+ipi and both atezolizumab
   cases sit in PD_. Residualizing each gene on prefix mean does not create a
   response association.
2. **Responders are few and regimen-mixed.** Only 8 tumors are CR/PR. Four of
   those eight are combination or PD-L1 monotherapy (3 nivo+ipi, 1 atezo).
   Nivolumab-only ORR is 4 versus 23.
3. **The curated PFS event field is not usable.** All 16 RECIST PD tumors are
   coded `event_occurred_pfs=0`, while some CR/PR tumors are coded as events.
   That is incompatible with a progression endpoint. PFS was therefore not
   analyzed. OS event coding is at least directionally plausible (23/33 events)
   and is reported only as exploratory.
4. **This is not the Braun CheckMate RNA set.** Braun 2020 (B5_Braun) is a
   larger nivolumab RCC transcriptome. Miao 2018 is a separate, much smaller
   public RNA slice of a WES-first paper.

## Reproduction

From this directory:

```bash
python3 -m pip install -r requirements.txt
python3 analyze.py
```

The script downloads `ICB_Miao1.zip` from Zenodo into `.cache/` when absent,
checks SHA-256 `121c7a843aa3ceb95c8ec58354b92d131d998d56d029762ea8a3ab46bca07ffe`
and the Zenodo MD5, then regenerates the tables and figure.

Outputs:

- `cldn4_tacstd2_vs_response.png` / `.svg`: primary ORR plus prefix strata
- `summary.csv`: all prespecified contrasts, effect sizes, and CIs
- `sample_level.csv`: one row per RNA sample
- `os_cox.csv`: exploratory overall-survival models
- `results.json`: source, integrity checks, and primary numbers
- `provenance.tsv`: URL, DOI, license, and checksums

## Limitations

This is a post hoc, two-gene association test in a sequenced subset of a
published ICI series. No multiplicity correction was applied beyond naming
both genes as prespecified; even the smallest p-value would not survive a
six-test view of two genes × three response labels. The bootstrap intervals
are descriptive. Response association in treated patients cannot distinguish
a predictive ICI interaction from a general prognostic association. Tumor
purity, biopsy site, IMDC risk, and prior VEGF therapy are not in the public
TSV used here and were not imputed.

## References

- Miao D, Margolis CA, Gao W, et al. Genomic correlates of response to immune
  checkpoint therapies in clear cell renal cell carcinoma. *Science*.
  2018;359:801–806.
  [doi:10.1126/science.aan5951](https://doi.org/10.1126/science.aan5951)
- Haibe-Kains B, et al. ICB Data TSV Files. Zenodo.
  [doi:10.5281/zenodo.7058399](https://doi.org/10.5281/zenodo.7058399)
  (BHK Lab / ORCESTRA `ICB_Miao1`; CC-BY-4.0)
