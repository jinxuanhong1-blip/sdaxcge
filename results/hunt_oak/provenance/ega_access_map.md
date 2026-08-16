# EGA access map — EGAS00001005013 (the OAK/POPLAR RNA-seq Bessede used)

Study: **EGAS00001005013** — "Intratumoral plasma cells predict outcomes to PD-L1
blockade in non-small cell lung cancer" (Patil et al., *Cancer Cell* 2022). This is
the source of the 891-tumor POPLAR+OAK transcriptome collection that Bessede et al.
2024 (*Clin Cancer Res* 30:779-785) reused for the TACSTD2/TROP2 analysis.

Every sub-dataset is **controlled access**, gated by the Genentech Data Access
Committee **EGAC00001002120**. A Data Access Committee (DAC) request/approval is
required; nothing here is downloadable openly. Checked 2026-08-16.

| EGA dataset ID   | Contents                                             | Trial   | Access     |
|------------------|------------------------------------------------------|---------|------------|
| EGAD00001007703  | Full tumor transcriptomes, 891 samples (FASTQ)       | both    | Controlled |
| EGAD00001008390  | log2(TPM+1) matrix, 192 samples                      | POPLAR  | Controlled |
| EGAD00001008391  | log2(TPM+1) matrix, 699 samples                      | OAK     | Controlled |
| EGAD00001008628  | raw count table, 699 samples                         | OAK     | Controlled |
| EGAD00001008629  | counts-per-million table, 699 samples                | OAK     | Controlled |
| EGAD00001008630  | raw count table, 192 samples                         | POPLAR  | Controlled |
| EGAD00001008631  | counts-per-million table, 192 samples                | POPLAR  | Controlled |
| EGAD00001008548  | clinical (arm, histology, OS, PFS, best response)    | POPLAR  | Controlled |
| EGAD00001008549  | clinical (arm, histology, OS, PFS, best response)    | OAK     | Controlled |

Additional OAK/POPLAR clinical data: available only via request through
**vivli.org** (Roche/Genentech data-sharing) — also gated, not open.

Bessede et al. data-availability statement (verbatim gist): raw and processed
transcriptomic data + limited clinical data are at EGA under EGAS00001005013;
additional clinical data via vivli.org; the immunofluorescence data are not public
and require ethics-committee approval + request to the corresponding author.

## What this means for an open re-test of "TACSTD2-high -> atezo resistance" in OAK
The per-sample OAK (and POPLAR) `TACSTD2` expression joined to atezolizumab outcome
— i.e. the exact table needed to reproduce Bessede on OAK — is **not** available in
any open form. The log2(TPM+1), counts, CPM, and clinical tables that would contain
it are all behind the Genentech DAC. Bessede's own supplements
(Supplementary Tables S1/S3/S4) are patient-*characteristic* summary tables, not
per-sample expression + outcome matrices, so they cannot be used to recompute the
TACSTD2 hazard ratios either.

**Verdict for OAK itself: no open leftover exists that can test TACSTD2 vs atezo
benefit.** The honest fallback is to test the same hypothesis in independent,
genuinely-open atezolizumab / PD-(L)1 cohorts (see ../README.md).
