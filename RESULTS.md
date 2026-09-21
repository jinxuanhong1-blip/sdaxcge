# CLDN4 / TACSTD2 in public ICI NSCLC bulk

OAK and POPLAR linear RNA were not scored. Both TPM matrices and the clinical tables are controlled-access EGA (checked 2026-09-21). The public bulk series that do contain CLDN4 or TACSTD2 and a response label were swept. The thesis slide’s 11-cohort pooled OR 0.42 (95% CI 0.18–0.95) for high CLDN4 and worse ICI response was not re-estimated. These public series are a different collection of endpoints and sample sizes, and they were not pooled into one odds ratio.

Odds ratio is the odds of the labeled benefit in the high group versus the low group. OR below 1 means high expression has worse benefit. The pre-specified split is the median, with ties at the cut called high. Quantiles 0.25, 0.33, 0.67, and 0.75, and a Q3-versus-Q1 tertile contrast, were searched afterwards and were not multiplicity-adjusted. A row is called stable when both groups have at least 5 samples and no cell is zero (no Haldane correction). Fisher p-values are two-sided. Intervals are Woolf.

## Pre-specified CLDN4 median, unadjusted

| Cohort | Endpoint | n | Benefit high / non-benefit high | Benefit low / non-benefit low | OR | 95% CI | Fisher p |
|---|---|---:|---|---|---:|---|---:|
| GSE126044 | Deposited responder label | 16 | 1 / 7 | 4 / 4 | 0.143 | 0.012–1.76 | 0.282 |
| GSE166449 | Deposited responder label | 22 | 3 / 8 | 4 / 7 | 0.656 | 0.108–4.00 | 1.000 |
| GSE135222 | PFS ≥ 180 days | 27 | 2 / 12 | 5 / 8 | 0.267 | 0.041–1.73 | 0.209 |
| GSE207422 | Major pathologic response | 24 | 2 / 10 | 7 / 5 | 0.143 | 0.021–0.958 | 0.089 |
| GSE207422 | RECIST CR/PR versus SD | 24 | 7 / 5 | 10 / 2 | 0.280 | 0.042–1.88 | 0.371 |
| GSE190265 | PFS ≥ 6 months | 43 | 8 / 14 | 6 / 15 | 1.43 | 0.395–5.16 | 0.747 |
| GSE190266 | Still progression-free at a 6-month cap | 69 | 11 / 24 | 6 / 28 | 2.14 | 0.688–6.65 | 0.265 |

Five of these seven medians have OR below 1. The two larger anti-PD-1 monotherapy series, Dijon France3 (GSE190265) and France4 (GSE190266), have OR above 1. No unadjusted CLDN4 median has a Fisher p-value below 0.05. The GSE207422 pathologic-response row is the only one whose Woolf interval excludes 1; its Fisher p-value is 0.089. That cohort is neoadjuvant toripalimab plus platinum, and the endpoint is major pathologic response, including MPR (pCR), versus NMPR.

## Searched CLDN4 cuts

Among stable unadjusted CLDN4 rows with OR below 1, the odds ratio closest to 0.42 is a searched cut, not the median:

- GSE126044, upper third versus the rest (quantile 0.67), n = 16 (5 high, 11 low), 1 and 4 versus 4 and 7, OR 0.438 (95% CI 0.035–5.40), Fisher p = 1.0.

That interval runs from well below 0.42 to well above 1. It is the nearest stable public OR in this sweep. It is not a reproduction of the slide’s pooled estimate.

The strongest stable unadjusted CLDN4 odds ratios below 1 are tied at 0.143:

| Cohort | Cut | 2×2 (benefit/non-benefit, high then low) | 95% CI | Fisher p |
|---|---|---|---|---:|
| GSE207422 MPR | Median (pre-specified) | 2/10 vs 7/5 | 0.021–0.958 | 0.089 |
| GSE207422 MPR | Quantile 0.67 (searched) | 1/7 vs 8/8 | 0.014–1.44 | 0.178 |
| GSE126044 deposited responder | Median (pre-specified) | 1/7 vs 4/4 | 0.012–1.76 | 0.282 |

Full rows are in `results/public_ici_nsclc_cldn4/sweep.tsv`. The machine-readable picks are in `summary.json`.

## Keratin residual

Each scored gene was also residualized by ordinary least squares on an intercept plus KRT8, KRT18, and KRT19, then split the same way. Exact symbols only. At the median, unadjusted and keratin-residual CLDN4 agree in GSE126044 (OR 0.143) and GSE190266 (OR 2.14). They move in the other cohorts:

| Cohort | Endpoint | Keratin-residual CLDN4 median OR | 95% CI | Fisher p |
|---|---|---:|---|---:|
| GSE166449 | Deposited responder | 1.52 | 0.250–9.30 | 1.000 |
| GSE135222 | PFS ≥ 180 days | 0.614 | 0.108–3.49 | 0.678 |
| GSE207422 | MPR | 0.700 | 0.133–3.68 | 1.000 |
| GSE207422 | RECIST CR/PR versus SD | 0.0909 | 0.0088–0.943 | 0.069 |
| GSE190265 | PFS ≥ 6 months | 0.933 | 0.261–3.34 | 1.000 |

The GSE207422 ORR residual is 6/6 versus 11/1. It is the smallest stable CLDN4 odds ratio in the whole sweep, including keratin rows. The interval excludes 1 and the Fisher p-value is 0.069. It is a residual of a chemotherapy-combination RECIST contrast, and the deposited RECIST column contains CR, PR, and SD only, so disease control cannot be contrasted.

## TACSTD2 in the same public files

Unadjusted median:

| Cohort | Endpoint | 2×2 | OR | 95% CI | Fisher p |
|---|---|---|---:|---|---:|
| GSE126044 | Deposited responder | 2/6 vs 3/5 | 0.556 | 0.065–4.76 | 1.000 |
| GSE166449 | Deposited responder | 5/6 vs 2/9 | 3.75 | 0.540–26.0 | 0.361 |
| GSE135222 | PFS ≥ 180 days | 3/11 vs 4/9 | 0.614 | 0.108–3.49 | 0.678 |
| GSE207422 | MPR | 3/9 vs 6/6 | 0.333 | 0.059–1.88 | 0.400 |
| GSE207422 | RECIST CR/PR versus SD | 8/4 vs 9/3 | 0.667 | 0.113–3.93 | 1.000 |
| GSE190265 | PFS ≥ 6 months | 8/14 vs 6/15 | 1.43 | 0.395–5.16 | 0.747 |

GSE190265 CLDN4 and TACSTD2 share this 2×2. The high sets are not the same people: Spearman correlation is 0.35 and 14 of 22 high calls overlap. GSE190266 has no TACSTD2 column. The deposited TPM header stops at MTMR14 (16,383 genes).

The only median row in the sweep with Fisher p below 0.05 is keratin-residual TACSTD2 in GSE135222: PFS ≥ 180 days, 1/13 versus 6/7, OR 0.0897 (95% CI 0.0089–0.902), Fisher p = 0.033. The unadjusted TACSTD2 median in that series is OR 0.614. This is a keratin residual of a PFS binary.

## Published OAK / POPLAR TACSTD2 summary (not recomputed)

Bessede et al., *Clinical Cancer Research* 2024 (10.1158/1078-0432.CCR-23-2566), report TACSTD2 on the atezolizumab arm of OAK plus POPLAR, n = 405. The high/low cut is a maxstat threshold fit on PFS, not a median. They report median PFS 2.5 versus 4.1 months (P < 0.001), median OS 12.6 versus 16.3 months (P = 0.007), and durable clinical benefit 15.7% versus 26.2% (P = 0.009). Durable clinical benefit is objective response or stable disease lasting at least 12 months. The same associations were not seen on docetaxel. The main text does not give the two group sizes, so those percentages were not turned into an odds ratio here. The main text that was read does not report a CLDN4 result. The FASTQ they used came from EGA DAC EGAC00001002120. Those files were not downloaded in this session.

The same paper’s BIP tumor-transcript set is described as 72 ICI-treated NSCLC tumors, without a numeric CLDN4 or TACSTD2 odds ratio in the main text read here. Plasma TROP2 protein in 74 BIP patients is reported as durable clinical benefit 17.5% versus 52.9% (P = 0.003). That is a protein measurement, and it was not recomputed.

## Cohorts that were scored

| Cohort | What was joined | Expression | Benefit counts |
|---|---|---|---|
| GSE126044 | 16 pre-treatment NSCLC biopsies | Counts converted to log2(CPM+1) with the full-library size | 5 responders, 11 non-responders, GEO labels |
| GSE166449 | 22 pre-treatment advanced lung cancers | Values already on a log scale (CLDN4 max 3.13); not logged again | 7 responders, 15 non-responders, from the sample title |
| GSE135222 | 27 tumor RNA samples, anti-PD-1/PD-L1 | Ensembl counts to log2(count+1) | 7 with PFS ≥ 180 days, 20 events before 180 days |
| GSE207422 | 24 pre-treatment biopsies with a patient id and a response | Deposited log2 TPM | MPR 9 vs NMPR 15; CR/PR 17 vs SD 7 |
| GSE190265 | 43 Dijon anti-PD-1 tumors at diagnosis | log2(TPM+1). Sample id is the extra leading field | 14 with PFS ≥ 6, 29 events before 6. Times 0.1 to 46.4 |
| GSE190266 | 69 of 70 Dijon anti-PD-1 lung tumors | log2(TPM+1). Decimals use a comma | 17 at the cap, 52 events before it |

GSE135222’s public flag was checked before use: samples with flag 1 have mean PFS 65 days (n = 21) and samples with flag 0 have mean 348 days (n = 6), so flag 1 was treated as the event. Benefit is time ≥ 180 days. Early censoring would have been left missing; all 27 samples received a label.

GSE190266’s time field is the GEO string `pfs_time (6 months)`. Observed times run from 0.03 to 6.0. All 17 samples at 6.0 have event 0. Sample A16-1489-1 (time 0.03, event 0) was left out as an early censor. The binary is whether the deposited time reached that 6-month cap.

GSE207422 has 29 metadata rows. Five have no patient id and no response and were dropped. The remaining 24 patient ids are unique. There is no PD in the RECIST column.

## Inspected and not scored

| Series | Why it is not in the odds-ratio table |
|---|---|
| OAK EGAD00001008391 and POPLAR EGAD00001008390 | `access_type = controlled` on the EGA metadata API. No credential was used |
| TCCIA OAK / POPLAR rows | The app’s analyses are circRNA ID and host-gene back-splice sums. Those sums were not used as CLDN4 or TACSTD2 mRNA. Rows 21–25 of the Repository were not captured |
| GSE93157 | Immune-panel matrix. EPCAM is present. CLDN4 and TACSTD2 are absent |
| GSE161537 | Targeted RNA, n = 82, RECIST present. CLDN3, KRT8, KRT18, KRT19, EPCAM, and CD274 are present. CLDN4 and TACSTD2 are absent |
| GSE111414 | PBMC CD8+ T cells from NSCLC patients on PD-1 blockade. CLDN4 and TACSTD2 are absent from the expression file |
| GSE248378 | Resected tumors after neoadjuvant durvalumab ± radiation. CLDN4 and TACSTD2 are present. Series-matrix characteristics are histology and Arm1/Arm2. No response field |
| GSE248249 | Affymetrix acquired-resistance FFPE series. The series matrix uses probe ids. A scan of the first 200,000 lines found no CLDN4 or TACSTD2 symbol. Probes were not remapped. There is no responder-versus-non-responder contrast |
| GSE249000 | Mouse CT26 acquired-resistance RNA-seq, not patient NSCLC |
| GSE246922 | Mouse LLC1 and KP cell-line RNA-seq, not patient NSCLC |
| SU2C-MARK | dbGaP phs002822.v1.p1, controlled. The Nature Genetics data note describes 152 RNA-seq profiles. Not downloaded |
| PLOS One 10.1371/journal.pone.0260500 | 40 nivolumab-treated NSCLC tumors. The data-availability statement points at the paper and supporting information. The supporting items listed on the article page are figures plus DEG and GSEA tables, not a sample-by-gene matrix |
| Bessede BIP RNA and plasma | Published summaries only. Individual expression was not public in the files used here |

A GEO DataSets query for “NSCLC immunotherapy RNA-seq response” on 2026-09-21 returned 19 hits. The human tumor bulk series with a usable response label and CLDN4 in the same public file are the six cohorts above. Blood, single-cell, cell-line, and mouse hits were not scored.

## DAC checklist for OAK / POPLAR

The contrast that is still blocked is median CLDN4 and TACSTD2 versus best confirmed response and OS, inside atezolizumab and inside chemotherapy, with a keratin adjustment. Ensembl ids to match in the TPM header: CLDN4 ENSG00000189143, TACSTD2 ENSG00000184292, KRT8 ENSG00000170421, KRT18 ENSG00000111057, KRT19 ENSG00000171345.

1. Create an individual EGA account.
2. Open the Genentech DAC page [EGAC00001002120](https://ega-archive.org/dacs/EGAC00001002120). Study page: [EGAS00001005013](https://ega-archive.org/studies/EGAS00001005013).
3. Request these datasets: EGAD00001008390 (POPLAR log2(TPM+1), 192 tumors), EGAD00001008391 (OAK log2(TPM+1), 699 tumors), EGAD00001008548 (POPLAR clinical), EGAD00001008549 (OAK clinical). Add EGAD00001008550 if the PD-L1 / TMB / STK11 / KEAP1 / EGFR table is needed.
4. Download the academic Data Access Agreement and the Data Access Request Form from [http://research-pub.gene.com/seqdata/dataAccess.html](http://research-pub.gene.com/seqdata/dataAccess.html).
5. Sign them and email them to DevSci DAC, `devsci-dac-d@gene.com`. Policy [EGAP00001002123](https://metadata.ega-archive.org/policies/EGAP00001002123) is “Genentech Data Access Agreement for Academic Institutions”: three-year term, then destroy the data; somatic mutation and/or gene-expression use; a collaborator outside the institution needs a separate application.
6. Treat EGA “Request access” and the signed Genentech agreement as one application. Approval has to come from Genentech.
7. After approval, download with pyEGA3:

```bash
pyega3 -cf CREDENTIALS_FILE fetch EGAD00001008390 --output-dir data/oak_poplar
pyega3 -cf CREDENTIALS_FILE fetch EGAD00001008391 --output-dir data/oak_poplar
pyega3 -cf CREDENTIALS_FILE fetch EGAD00001008548 --output-dir data/oak_poplar
pyega3 -cf CREDENTIALS_FILE fetch EGAD00001008549 --output-dir data/oak_poplar
```

8. Do not commit those files. The agreement forbids transferring the data, or identifiable derivatives, outside the approved institution.
9. Patil et al. also point to vivli.org for additional clinical data. That request is separate and was not made.
10. This session did not submit a DAC request and did not use a credential file.

| Dataset | File | Size (bytes) |
|---|---|---:|
| EGAD00001008390 POPLAR log2(TPM+1) | EGAF00005797825 | 50,248,881 |
| EGAD00001008391 OAK log2(TPM+1) | EGAF00005797824 | 283,985,625 |
| EGAD00001008548 POPLAR clinical | EGAF00005797821; EGAF00006143195; EGAF50000088676 | 28,009; 44,439; 29,943 |
| EGAD00001008549 OAK clinical | EGAF00005797822; EGAF00006143195 | 94,237; 44,439 |
| EGAD00001008550 OAK biomarkers | EGAF00007237036; EGAF50000086634 | 146,254; 26,338 |

Counts, CPM, and FASTQ remain controlled and are not required once the TPM matrices are in hand: OAK EGAF00006072294 and EGAF00006072295, POPLAR EGAF00006072292 and EGAF00006072293, FASTQ EGAD00001007703 (891 samples, 1,782 files). Checksums are in `results/oak_poplar_access/access_status.json`.

TCCIA (https://shiny.hiplot.cn/TCCIA/ and https://shiny.zhoulab.ac.cn/TCCIA/) lists Patil OAK as EGAD00001008549, n = 699, and Patil POPLAR as EGAD00001008548, n = 192. Those ids are the clinical datasets. Zenodo 10.5281/zenodo.7969298 is a circRNA ID catalog. GEO and cBioPortal returned no OAK/POPLAR TPM matrix.

## Reproduce

```bash
python3 scripts/check_oak_poplar_access.py
python3 scripts/public_ici_nsclc_cldn4_or.py
```

The access script exits if any listed EGA dataset is not `controlled`. The sweep script reads public GEO files only, writes `results/public_ici_nsclc_cldn4/sweep.tsv` and `summary.json`, and refuses to binarize a PFS flag if the event class does not have the shorter mean time. Downloaded matrices stay in `data/` and are gitignored.
