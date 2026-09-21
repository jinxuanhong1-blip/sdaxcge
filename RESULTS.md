# OAK / POPLAR — CLDN4 and TACSTD2 vs atezolizumab

No ORR, DCR, or OS result. The linear RNA and the clinical tables are controlled-access EGA. This session did not download them and did not compute a high-vs-low contrast.

The slide claim this was meant to sit beside is page 12 of the thesis deck: 11 ICI cohorts, pooled OR 0.42 (95% CI 0.18–0.95) for high CLDN4 and lower response. That pooled figure is not re-estimated here, and it is not treated as failed. OAK (GO28915) and POPLAR (GO28753) would have been a separate, randomized atezolizumab-vs-chemotherapy test. That test did not run.

## What was asked

Pretreatment tumor RNA from the two randomized 2L NSCLC trials in Patil et al., *Cancer Cell* 2022 (10.1016/j.ccell.2022.02.002):

- CLDN4 high vs low, and TACSTD2 high vs low
- ORR, DCR, and OS
- inside atezolizumab and inside chemotherapy, so the arm contrast is visible
- keratin adjustment (KRT8 / KRT18 / KRT19) if the same matrix contains those genes

## What is public, and what was opened

| Place checked | What came back |
|---|---|
| EGA metadata API, 2026-09-21 | `access_type = controlled` for every TPM, count, CPM, FASTQ, and clinical dataset below |
| EGA file catalog | File accessions and sizes only. Bytes were not fetched |
| TCCIA https://shiny.hiplot.cn/TCCIA/ | OAK and POPLAR are in the Repository. The app analyzes circRNA, not linear TPM |
| TCCIA mirror https://shiny.zhoulab.ac.cn/TCCIA/ | HTTP 200, same Shiny app |
| http://tccia.zmu.zhoulab.com/ and http://biotrainee.vip:18888/TCCIA/ | No HTTP response in this session |
| Zenodo 10.5281/zenodo.7969298 (record 10091408) | CircRNA ID catalogs (`TCCIA_Ensemble_circRNAs.tsv.gz`, `TCCIA_4_methods_circRNAs.tsv.gz`). Not a sample-by-gene TPM |
| GEO `esearch` for OAK/POPLAR atezolizumab RNA-seq | 0 datasets |
| cBioPortal studies keyword `OAK` | empty |

No EGA account and no DAC approval were used. `pyega3` was not pointed at a credential file.

## TCCIA is not the TPM matrix

The live Repository table lists both trials. Rows 1–20 of 25 were read; rows 21–25 were not. The two NSCLC rows are:

| Cohort label in TCCIA | ID shown | n | Treatment text | Drug text | Sampling | Survival fields |
|---|---|---:|---|---|---|---|
| Patil et al 2022 (NSCLC) - OAK | EGAD00001008549 | 699 | anti-PD-L1; Chemotherapy | Atezolizumab | Pre | OS, PFS |
| Patil et al 2022 (NSCLC) - POPLAR | EGAD00001008548 | 192 | anti-PD-L1; Chemotherapy | Atezolizumab | Pre | OS, PFS |

Those IDs are the **clinical** datasets, not EGAD00001008391 / EGAD00001008390. The cohort-centred controls expose two analysis types only: `CircRNA ID` and `Map to Host Gene`. Host-gene mode sums back-splice junction counts of circRNAs that overlap a gene. That is not CLDN4 or TACSTD2 mRNA, and it was not used as a stand-in.

TCCIA’s own data-availability note (JITC 2024, 10.1136/jitc-2023-008040) says processed circRNA matrices are requested from the project leader. The site also says analysis modules keep immunotherapy-related samples by default, with a Setting-page switch for other arms. The chemotherapy contrast was not queried through that switch.

Wang et al. state that Patil et al. contributed two cohorts and that the raw RNA-seq was obtained after DAC approval. IOBR signatures inside TCCIA were computed by the authors from expression they are not redistributing as a public TPM.

## DAC path

Study: [EGAS00001005013](https://ega-archive.org/studies/EGAS00001005013)

DAC: [EGAC00001002120](https://ega-archive.org/dacs/EGAC00001002120) (Genentech DAC)

Policy: [EGAP00001002123](https://metadata.ega-archive.org/policies/EGAP00001002123), title “Genentech Data Access Agreement for Academic Institutions”, `dac_accession_id` EGAC00001002120

DAC contact on the EGA record: DevSci DAC, `devsci-dac-d@gene.com`, Genentech

Form page linked from public descriptions of this DAC: [http://research-pub.gene.com/seqdata/dataAccess.html](http://research-pub.gene.com/seqdata/dataAccess.html). Academic applicants download the academic Data Access Agreement and the Data Access Request Form, sign them, and email them to the address on the EGA DAC page. The agreement term on the public policy text is three years, then the data are destroyed. Use is limited to the approved project, somatic mutation and/or gene-expression analysis, and a collaborator outside the institution needs a separate application.

EGA-side request, same DAC: create an individual EGA account, open the DAC page, use Request access, select the datasets, and send the request. Approval still has to come from Genentech. Patil et al. also point to vivli.org for additional clinical data; that request is separate and was not made.

Files to request for the contrast that was not run:

| Dataset | What it is | File | Size (bytes) |
|---|---|---|---:|
| EGAD00001008390 | POPLAR log2(TPM+1), 192 tumors, GO28753 | EGAF00005797825 csv | 50,248,881 |
| EGAD00001008391 | OAK log2(TPM+1), 699 tumors, GO28915 | EGAF00005797824 csv | 283,985,625 |
| EGAD00001008548 | POPLAR arm, histology, OS, PFS, best confirmed overall response | EGAF00005797821 csv; EGAF00006143195 csv; EGAF50000088676 | 28,009; 44,439; 29,943 |
| EGAD00001008549 | OAK arm, histology, OS, PFS, best confirmed overall response | EGAF00005797822 csv; EGAF00006143195 csv | 94,237; 44,439 |
| EGAD00001008550 | OAK PD-L1 22C3, TMB, STK11 / KEAP1 / EGFR | EGAF00007237036 txt; EGAF50000086634 | 146,254; 26,338 |

Same study, still controlled, not required once the TPM matrices are in hand: OAK counts EGAF00006072294, OAK CPM EGAF00006072295, POPLAR counts EGAF00006072292, POPLAR CPM EGAF00006072293, and the FASTQ set EGAD00001007703 (891 samples; the public catalog lists 1,782 files). Checksums for the TPM and clinical files are in `results/oak_poplar_access/access_status.json`. The FASTQ dataset is recorded there as a file count only.

After approval, download with pyEGA3 and the EGA credentials file:

```bash
pyega3 -cf CREDENTIALS_FILE fetch EGAD00001008390 --output-dir data/oak_poplar
pyega3 -cf CREDENTIALS_FILE fetch EGAD00001008391 --output-dir data/oak_poplar
pyega3 -cf CREDENTIALS_FILE fetch EGAD00001008548 --output-dir data/oak_poplar
pyega3 -cf CREDENTIALS_FILE fetch EGAD00001008549 --output-dir data/oak_poplar
```

Do not commit those files. The agreement forbids transferring the data, or identifiable derivatives, outside the approved institution.

## Analysis that was not run

Once those four datasets are local, the contrast is:

- Expression: `CLDN4` and `TACSTD2` from the log2(TPM+1) matrices. Ensembl IDs to confirm against the matrix header: ENSG00000189143 and ENSG00000184292.
- Split: within each trial and arm, median high vs low. A tertile sensitivity (Q3 vs Q1) if the arm is large enough.
- Response: best confirmed overall response. ORR = CR or PR versus not. DCR = CR, PR, or SD versus PD. Report the exact strings in the clinical file before coding them.
- OS: time and event as released, Cox within arm, then a gene-by-arm interaction on the pooled patients of that trial.
- Keratin: the TPM matrix is transcriptome-wide, so KRT8, KRT18, and KRT19 should be in the same file. Adjust CLDN4 (and, separately, TACSTD2) for those three before the high/low call, and fit a second model with them as covariates. Skip the adjustment if any of the three symbols is absent, and say which one.

No odds ratio, hazard ratio, response count, or keratin coefficient is reported, because the inputs were not available.

## Reproduce the access check

```bash
python3 scripts/check_oak_poplar_access.py
```

The script reads public metadata, writes `results/oak_poplar_access/access_status.json`, and exits if any listed dataset is not `controlled`. It does not compute associations. The Repository excerpt is `results/oak_poplar_access/tccia_repository_rows_1_to_20.tsv`.
