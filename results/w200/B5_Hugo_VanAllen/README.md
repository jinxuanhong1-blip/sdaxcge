# Hugo 2016 and Van Allen 2015 melanoma ICI: TACSTD2 and CLDN4 versus response

## Bottom line

Neither public melanoma ICI RNA cohort supports TACSTD2 or CLDN4 as a
response marker. Pretreatment **CLDN4 is null in both studies**. Pretreatment
**TACSTD2 point estimates go in opposite directions** (lower in Hugo anti–PD-1
responders; higher in Van Allen anti–CTLA-4 responders) and neither contrast
is significant. These are melanoma datasets, not lung, and they are not a
replication of the same drug.

## Cohorts and endpoints

### Hugo 2016 (GSE78220)

- Source: Hugo et al., *Cell* 2016, GEO
  [GSE78220](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78220).
- Treatment: pembrolizumab (anti–PD-1) in advanced melanoma.
- Public RNA: 28 biopsies, GEO FPKM matrix `GSE78220_PatientFPKM.xlsx`.
- Genes: **TACSTD2** and **CLDN4** by official symbol (one row each).
- Response: GEO `anti-pd-1 response` recoded to RECIST CR / PR / PD. There is
  no SD in this series.
- Primary analysis: one pretreatment sample per patient. Pt16 is on-treatment
  and was excluded. Pt27A and Pt27B are two pretreatment biopsies from the same
  CR patient; their `log2(FPKM+1)` values were averaged.
- Binary endpoint: CR/PR = responder (n=14); PD = nonresponder (n=12).
- Expression scale: `log2(FPKM + 1)`. Mann–Whitney p-values are unchanged by
  this monotonic transform.

### Van Allen 2015

- Source: Van Allen et al., *Science* 2015. Raw RNA is dbGaP
  [phs000452](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000452)
  (controlled). The public analysis uses the lab-released STAR/RSEM TPM matrix
  (Gencode v19) for 42 pretreatment tumors
  ([vanallenlab/VanAllen_CTLA4_Science_RNASeq_TPM](https://github.com/vanallenlab/VanAllen_CTLA4_Science_RNASeq_TPM)).
- Treatment: ipilimumab (anti–CTLA-4).
- Genes: **TACSTD2** `ENSG00000184292`, **CLDN4** `ENSG00000189143`.
- Response: cBioPortal study `skcm_dfci_2015`. The field named
  `DURABLE_CLINICAL_BENEFIT` is RECIST (CR/PR/SD/PD/X), not the paper's
  clinical-benefit label.
- 40/42 RNA patients have a cBioPortal clinical row. **Pat20 and Pat91 have
  public TPM but no public RECIST** in that table. Science Table S2 was not
  available without a browser challenge; those two patients were excluded
  rather than labeled from a third-party recode.
- Primary endpoint: RECIST CR/PR vs PD (7 vs 26). SD (n=6) and X (n=1) were
  excluded from the ORR contrast.
- Secondary endpoint: paper clinical benefit (CR/PR, or SD with OS > 12
  months) vs no benefit (PD, or SD with OS ≤ 12 months): 12 vs 27. The X
  patient and the two unlabeled RNA samples were excluded.
- Expression scale: `log2(TPM + 1)`.

The two matrices are not on a common scale. They are not meta-analyzed.

## Results

Primary analyses are pretreatment, patient-level, RECIST CR/PR vs PD.
Medians are on the log2(x+1) scale used for the test. AUC > 0.5 means
responders have higher expression.

| Analysis | Gene | R / NR | Median R | Median NR | AUC (R higher) | 95% bootstrap CI | MW p |
|---|---|---:|---:|---:|---:|---:|---:|
| Hugo pre, patient | TACSTD2 | 14 / 12 | 0.503 | 0.729 | 0.351 | 0.143–0.577 | 0.208 |
| Hugo pre, patient | CLDN4 | 14 / 12 | 0.477 | 0.369 | 0.542 | 0.310–0.762 | 0.738 |
| Van Allen pre, RECIST | TACSTD2 | 7 / 26 | 0.390 | 0.104 | 0.698 | 0.500–0.868 | 0.118 |
| Van Allen pre, RECIST | CLDN4 | 7 / 26 | 0.029 | 0.050 | 0.360 | 0.179–0.566 | 0.261 |

Hugo Pt28 (PR) is an extreme high-expresser of both genes
(`log2(FPKM+1)` TACSTD2 9.33, CLDN4 7.22). Dropping that one patient does
not create a supported association; TACSTD2 moves further toward
responder-lower (AUC 0.301, CI 0.103–0.526, p=0.097) and CLDN4 stays null
(AUC 0.506, p=0.978).

Keeping both Pt27 biopsies as independent samples (15 vs 12) does not change
the Hugo conclusion (TACSTD2 AUC 0.361, p=0.232; CLDN4 AUC 0.522, p=0.864).

Van Allen paper clinical-benefit (12 vs 27) is the same story: TACSTD2 AUC
0.645 (0.443–0.824), p=0.157; CLDN4 AUC 0.448 (0.267–0.633), p=0.611.

Both genes are detected in every Hugo pretreatment patient. Van Allen CLDN4
is low (detected in 4/7 RECIST responders and 18/26 nonresponders; Fisher
p=0.661). Detection rate is not associated with response.

TACSTD2 and CLDN4 are correlated with each other (Hugo patients Spearman
ρ=0.532, p=0.005; Van Allen ρ=0.471, p=0.002). That is expected for two
epithelial/junction genes. It is not evidence that either one predicts ICI
response.

## What this does not show

1. It does not show that high TACSTD2 or high CLDN4 marks ICI response in
   melanoma. The only numerically higher TACSTD2 signal (Van Allen RECIST)
   has n_R=7, p=0.118, and a bootstrap CI that starts at 0.50. Hugo goes the
   other way.
2. It does not show the opposite either. Hugo TACSTD2 responder-lower is also
   non-significant, and one outlier (Pt28) is a high-TACSTD2 responder.
3. It is not a lung ICI result and should not be cited as one.
4. Anti–PD-1 and anti–CTLA-4 are different mechanisms. A null in one does not
   “replicate” a null in the other; it also does not rescue a weak trend in
   the other.

## Important limitations

1. Both primary analyses are small. Van Allen has only seven RECIST
   responders with public RNA and labels.
2. Hugo GEO characteristic rows after the mutation fields are not
   rectangular. Response, survival, and prior MAPKi were taken from the
   aligned early characteristic rows; on- vs pretreatment was taken from the
   FPKM column suffix (`.baseline` vs `.OnTx`), not from the shifted
   “biopsy time” row.
3. Public Hugo metadata do not support a credible multivariable model
   (purity, biopsy site, subtype). Prior MAPKi is available but splits the
   already small table further.
4. Van Allen raw counts remain in dbGaP. This analysis uses the public TPM
   release plus cBioPortal RECIST. Two RNA samples are unlabeled.
5. FPKM/TPM are suitable for a within-gene rank test. This is not a
   raw-count differential-expression analysis.
6. Association with RECIST is not a test of treatment-specific prediction
   versus prognosis in untreated melanoma.

## Reproduction

From this directory:

```bash
python3 -m pip install -r requirements.txt
python3 analyze.py
```

The script downloads the four public inputs when absent, checks gene identity
and sample counts, and regenerates tables and figures. `provenance.tsv`
records input URLs and SHA-256 hashes. Downloaded matrices are gitignored.

Outputs:

- `tacstd2_cldn4_vs_response.png` / `.svg`: primary RECIST boxplots
- `summary.csv`: tests, effect sizes, confidence intervals, detection
- `sample_level.csv`: parsed public annotations and expression
- `hugo_patient_level.csv`: Hugo patient collapse (Pt27 averaged)
- `gene_correlation.csv`: TACSTD2–CLDN4 Spearman
- `provenance.tsv`: source URLs and checksums

## References

- Hugo W, Zaretsky JM, Sun L, et al. Genomic and Transcriptomic Features of
  Response to Anti-PD-1 Therapy in Metastatic Melanoma. *Cell*.
  2016;165:35–44.
  [doi:10.1016/j.cell.2016.02.065](https://doi.org/10.1016/j.cell.2016.02.065)
- NCBI GEO accession
  [GSE78220](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78220)
- Van Allen EM, Miao D, Schilling B, et al. Genomic correlates of response to
  CTLA-4 blockade in metastatic melanoma. *Science*. 2015;350:207–211.
  [doi:10.1126/science.aad0095](https://doi.org/10.1126/science.aad0095)
- Van Allen lab public TPM matrix:
  [vanallenlab/VanAllen_CTLA4_Science_RNASeq_TPM](https://github.com/vanallenlab/VanAllen_CTLA4_Science_RNASeq_TPM)
- cBioPortal study `skcm_dfci_2015`
  ([datahub](https://github.com/cBioPortal/datahub/tree/master/public/skcm_dfci_2015))
