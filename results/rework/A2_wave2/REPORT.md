# REWORK A2 wave 2 — do TACSTD2/CLDN4 vs MPR/PFS/CYT/CD8 survive MKI67, ESTIMATE, and epithelial residuals?

**Self-contained. Public data only. Written to be read without the rest of the repo.**

**Verdict: PARTLY — leftover MPR gap is real but not residual-robust; post-treatment TACSTD2–CYT/CD8 negative correlation is residual-robust to ESTIMATE and remains after MKI67/epithelial; no new PACIFIC/ADRIA matrix.**

The user slide (ρ=−0.65 / purity −0.46, p=2e-4) still has **no open PACIFIC per-patient RNA**. ADRIATIC (NCT03703297) published group-level RNA-seq *medians* in 2025, not a matrix. The only open human durvalumab NSCLC whole-transcriptome matrices remain **GSE253564** (pre, n=32) and **GSE248378** (post, n=29) from NCT02904954.

## Why this rework exists

Wave 1 (`results/rework/A2_defs/`) recomputed TACSTD2 vs five immune definitions with partial Spearman on ESTIMATEScore ranks. It did **not** residualise on MKI67 or an epithelial score, and it did not test MPR/PFS. Leftover PR 175 found pre-treatment TACSTD2 *higher* in MPR (p=0.0155) that died after MKI67 OLS residual (p=0.25). This wave retests that leftover and the A2 immune claim after each covariate separately.

## Analysis set

| Series | Timepoint | n | Usable endpoints |
|---|---|---|---|
| GSE253564 | pre-treatment FPKM | 32 | MPR 11 vs 21; pathologic depth; PFS (7 events); CYT; CD8A |
| GSE248378 | post-treatment FPKM | 29 | **MPR not estimable** (0 MPR in matrix); recurrence 9 vs 20; PFS (9 events); CYT; CD8A |

Clinical join: Table S1 (PMID 38401548, `mmc2.xlsx`) + Nat Commun source data (PMID 38114518). Titles join 32/32 and 29/29. Arm conflicts: 1 (documented `45-M-PO` / durva045).

## Pre-specified features and covariates

| Item | Definition | Honest limitation |
|---|---|---|
| TACSTD2, CLDN4 | log2(FPKM+1) | Bulk, not protein |
| CYT | mean log2(FPKM+1) of GZMA+PRF1 (Rooney 2015) | Two-gene cytolytic score |
| CD8 | log2(FPKM+1) CD8A | Single gene, not a CD8 panel |
| MKI67 | log2(FPKM+1) | Single-gene proliferation; leftover covariate |
| ESTIMATE | ESTIMATEScore = Stromal+Immune (v1.0.13 port) | Expression-derived; cosine TumorPurity out of bounds on FPKM |
| EpiScore | mean log2(FPKM+1) of EPCAM, KRT7, KRT8, KRT18, KRT19, CDH1, MUC1 | Locked lung-epithelial markers. **TACSTD2 and CLDN4 excluded** (circularity). Present pre: EPCAM, KRT7, KRT8, KRT18, KRT19, CDH1, MUC1; post: EPCAM, KRT7, KRT8, KRT18, KRT19, CDH1, MUC1. Single-gene EPCAM is sensitivity only.

Residualisation is **separate, never joint**. Primary MPR/PFS tests use leftover-style OLS residual of the gene on the covariate, then Mann–Whitney / Cox. Immune tests also report A2-matched partial Spearman (Pearson on ranks, n−3 df). No multiple-testing correction — every p is shown raw.

## Direct answer

1. **PACIFIC / ADRIATIC per-patient matrix: still none open.** This re-hunt verified that live. Entrez returned 8 unique GSE accessions after dropping GSM/GPL children. The bare keyword PACIFIC without durvalumab is polluted by PacBio / Pacific-Islander series and was not used as a hit list. None of the GSE hits is a new PACIFIC or ADRIATIC per-patient whole-transcriptome matrix. Open usable lung durvalumab RNA remains the two NCT02904954 series (2 rows marked `yes` in `hunt_triage.tsv`).
2. **Pre-treatment TACSTD2 is higher in MPR, not lower** (median 7.302 vs 5.888; p=0.0155; Cliff δ=0.532). After MKI67 residual p=0.25 (δ=0.255); after ESTIMATE p=0.0473; after EpiScore p=0.0622. The leftover MKI67 kill is reproduced. This is not evidence that TROP2-high tumours are ICI-resistant.
3. **CLDN4 vs MPR is null** (p=0.341; MKI67 residual p=0.267).
4. **Post-treatment TACSTD2 Cox PFS** unadjusted HR=2.236 (1.162–4.303), p=0.016, 9 events. After MKI67 HR=2.073 p=0.00955; after ESTIMATE HR=1.675 p=0.174; after EpiScore HR=2.761 p=0.00459. Small-event finding.
5. **TACSTD2 vs CYT/CD8 (the A2 claim).** Post CYT partial | ESTIMATE ρ=-0.710 p=2.29e-05 (older hunt −0.71 reproduced). Post CD8A | ESTIMATE ρ=-0.577 p=0.0013. Post CYT | MKI67 ρ=-0.781; | EpiScore ρ=-0.826. Pre CYT | ESTIMATE ρ=-0.165 p=0.375 (collapses, as in wave 1). Pre CYT | MKI67 ρ=-0.590; | EpiScore ρ=-0.445.
6. **User coefficient match is still post ESTIMATE ImmuneScore**, not a new PACIFIC number: raw ρ=-0.665, partial | ESTIMATEScore ρ=-0.474 p=0.0109. User p=2e-4 is the raw-p order, not the adjusted p. p=2e-4 attached to adj ρ=−0.46 does not reproduce at n=29.

## Primary result — MPR / PFS after each residual

### GSE253564 pre-treatment MPR (11 vs 21)

| Gene | Covariate | n_MPR / n_no | median_MPR | median_no | Cliff δ | p |
|---|---|---|---|---|---|---|
| TACSTD2 | none | 11/21 | 7.302 | 5.888 | 0.532 | 0.0155 |
| TACSTD2 | MKI67 | 11/21 | 0.826 | 0.564 | 0.255 | 0.25 |
| TACSTD2 | ESTIMATEScore | 11/21 | 1.081 | -0.078 | 0.437 | 0.0473 |
| TACSTD2 | EpiScore | 11/21 | 0.818 | -0.600 | 0.411 | 0.0622 |
| TACSTD2 | EPCAM_sensitivity | 11/21 | 0.888 | -0.456 | 0.351 | 0.113 |
| CLDN4 | none | 11/21 | 5.975 | 5.659 | 0.212 | 0.341 |
| CLDN4 | MKI67 | 11/21 | 0.504 | 0.104 | 0.247 | 0.267 |
| CLDN4 | ESTIMATEScore | 11/21 | 0.350 | 0.083 | 0.065 | 0.781 |
| CLDN4 | EpiScore | 11/21 | 0.230 | -0.026 | 0.056 | 0.812 |
| CLDN4 | EPCAM_sensitivity | 11/21 | -0.043 | 0.097 | -0.091 | 0.692 |

### GSE248378 post-treatment PFS (9 events) and recurrence

| Gene | Covariate | Cox HR (95% CI) | Cox p | log-rank p | recurrence MWU p |
|---|---|---|---|---|---|
| TACSTD2 | none | 2.236 (1.162–4.303) | 0.016 | 0.166 | 0.0562 |
| TACSTD2 | MKI67 | 2.073 (1.194–3.599) | 0.00955 | 0.137 | 0.212 |
| TACSTD2 | ESTIMATEScore | 1.675 (0.797–3.520) | 0.174 | 0.52 | 0.465 |
| TACSTD2 | EpiScore | 2.761 (1.368–5.572) | 0.00459 | 0.418 | 0.195 |
| CLDN4 | none | 1.410 (0.863–2.301) | 0.17 | 0.251 | 0.311 |
| CLDN4 | MKI67 | 1.552 (0.829–2.903) | 0.169 | 0.268 | 0.334 |
| CLDN4 | ESTIMATEScore | 1.216 (0.749–1.975) | 0.43 | 0.798 | 0.588 |
| CLDN4 | EpiScore | 2.790 (0.794–9.805) | 0.11 | 0.436 | 0.556 |

Post MPR is **not estimable** (0 vs 29). Do not pool this accession as an MPR replication.

## TACSTD2 / CLDN4 vs CYT and CD8A

Two estimators: leftover-style Spearman of the OLS residual vs the *raw* immune score, and A2-matched partial Spearman. Both are in `immune_correlations.tsv`. Partial Spearman is the A2 comparator.

| Dataset | Gene | Immune | Covariate | partial ρ | partial p | OLS-resid ρ | OLS-resid p |
|---|---|---|---|---|---|---|---|
| GSE253564 | TACSTD2 | CYT | MKI67 | -0.590 | 0.000473 | -0.578 | 0.000535 |
| GSE253564 | TACSTD2 | CYT | ESTIMATEScore | -0.165 | 0.375 | -0.243 | 0.181 |
| GSE253564 | TACSTD2 | CYT | EpiScore | -0.445 | 0.0121 | -0.408 | 0.0206 |
| GSE253564 | TACSTD2 | CD8A | MKI67 | -0.682 | 2.35e-05 | -0.651 | 5.39e-05 |
| GSE253564 | TACSTD2 | CD8A | ESTIMATEScore | -0.243 | 0.189 | -0.301 | 0.0942 |
| GSE253564 | TACSTD2 | CD8A | EpiScore | -0.497 | 0.00448 | -0.455 | 0.00884 |
| GSE253564 | CLDN4 | CYT | MKI67 | -0.370 | 0.0403 | -0.413 | 0.0188 |
| GSE253564 | CLDN4 | CYT | ESTIMATEScore | -0.287 | 0.118 | -0.268 | 0.139 |
| GSE253564 | CLDN4 | CYT | EpiScore | -0.381 | 0.0344 | -0.516 | 0.0025 |
| GSE253564 | CLDN4 | CD8A | MKI67 | -0.422 | 0.018 | -0.456 | 0.00865 |
| GSE253564 | CLDN4 | CD8A | ESTIMATEScore | -0.353 | 0.0516 | -0.291 | 0.106 |
| GSE253564 | CLDN4 | CD8A | EpiScore | -0.387 | 0.0316 | -0.474 | 0.00614 |
| GSE248378 | TACSTD2 | CYT | MKI67 | -0.781 | 9.52e-07 | -0.711 | 1.52e-05 |
| GSE248378 | TACSTD2 | CYT | ESTIMATEScore | -0.710 | 2.29e-05 | -0.547 | 0.00212 |
| GSE248378 | TACSTD2 | CYT | EpiScore | -0.826 | 6.35e-08 | -0.738 | 4.81e-06 |
| GSE248378 | TACSTD2 | CD8A | MKI67 | -0.694 | 4.15e-05 | -0.633 | 0.000225 |
| GSE248378 | TACSTD2 | CD8A | ESTIMATEScore | -0.577 | 0.0013 | -0.488 | 0.00722 |
| GSE248378 | TACSTD2 | CD8A | EpiScore | -0.725 | 1.29e-05 | -0.648 | 0.000143 |
| GSE248378 | CLDN4 | CYT | MKI67 | -0.433 | 0.0215 | -0.364 | 0.0526 |
| GSE248378 | CLDN4 | CYT | ESTIMATEScore | -0.407 | 0.0317 | -0.379 | 0.0427 |
| GSE248378 | CLDN4 | CYT | EpiScore | -0.580 | 0.00121 | -0.568 | 0.00131 |
| GSE248378 | CLDN4 | CD8A | MKI67 | -0.414 | 0.0284 | -0.343 | 0.0686 |
| GSE248378 | CLDN4 | CD8A | ESTIMATEScore | -0.332 | 0.0843 | -0.328 | 0.0828 |
| GSE248378 | CLDN4 | CD8A | EpiScore | -0.531 | 0.00364 | -0.557 | 0.00171 |

Unadjusted (covariate = none) Spearman rows are in the same TSV.

## Hunt — PACIFIC / ADRIATIC / open durvalumab

Live Entrez `gds` + PubMed queries were run at analysis time (`hunt_manifest.tsv`). A GEO accession is listed only if Entrez returned it. Near-misses are not hits.

| Resource | ID | Open per-patient matrix? | Status |
|---|---|---|---|
| PACIFIC (durvalumab after CRT, stage III NSCLC) | `NCT02125461` | no | NOT_OPEN |
| ADRIATIC (consolidation durvalumab, LS-SCLC) | `NCT03703297` | no | NOT_OPEN |
| SUBMARINE / WJOG11518L PACIFIC-regimen NanoString | `PMID 37364849` | no | NOT_DEPOSITED |
| Neoadjuvant durvalumab ± SBRT pre-treatment FPKM | `GSE253564` | yes | ANALYZED |
| Neoadjuvant durvalumab ± SBRT post-treatment FPKM | `GSE248378` | yes | ANALYZED |
| Durvalumab Study 1108 IFNγ 21-gene panel | `GSE110390` | no | EXCLUDED_UNUSABLE |
| Durvalumab ± tremelimumab scRNA (2 NSCLC tumors) | `GSE131933` | no | EXCLUDED_WRONG_ASSAY |
| Adjuvant durvalumab esophageal/GEJ RNA-seq | `GSE183924` | yes_but_not_lung | EXCLUDED_WRONG_DISEASE |
| VCN-01 + durvalumab HNSCC RNA-seq | `GSE333537` | yes_but_not_lung | EXCLUDED_WRONG_DISEASE |
| CANTABRICO ES-SCLC GeoMx CTA (chemo-IO) | `GSE261345` | no | EXCLUDED_WRONG_ASSAY |
| Durvalumab/oleclumab NSCLC xenograft array | `GSE190731` | no | EXCLUDED_NOT_PATIENT |
| OAK/POPLAR atezolizumab (Bessede 2024 TROP2 paper) | `EGAS00001005013` | no | CONTROLLED_ACCESS |
| POSEIDON / MYSTIC / AEGEAN / PACIFIC-2/4/5/6/8/9 | `NCT03164616 / NCT02453282 / NCT03800134 / PACIFIC-x` | no | NOT_OPEN |

## Honest interpretation

1. Do not cite a single ρ. Pre-treatment CYT/CD8A vs TACSTD2 is mostly ESTIMATE/purity; post-treatment CYT/CD8A is not.
2. The leftover TACSTD2-higher-in-MPR result is the opposite of “TROP2-high = ICI-cold/resistant” and is not independent of proliferation.
3. n=32/29 and 7/9 PFS events. Only large effects are detectable. A null CLDN4 result does not rule out a modest effect.
4. ESTIMATE TumorPurity cosine is out of bounds for 23/32 pre and 6/29 post FPKM samples; ESTIMATEScore (not the cosine) is the covariate.
5. EpiScore is a locked marker mean, not a published deconvolution epithelial fraction. It is reported so the residual request can be tested without inventing a commercial assay.
6. ADRIATIC RNA-seq exists as a trial biomarker analysis (group medians). That is not a downloadable per-patient matrix.
7. p=2e-4 attached to adj ρ=−0.46 does not reproduce at these sample sizes (needs n≈55+).

## Reproduce

```bash
pip install -r requirements.txt
python3 scripts/rework_A2_wave2.py
```

## Files

- `REPORT.md` — this write-up
- `residual_tests.tsv` — MPR / recurrence / depth / PFS after each residual
- `immune_correlations.tsv` — vs CYT / CD8A / ImmuneScore, both estimators
- `hunt_manifest.tsv` — live Entrez log (search counts + GSE hits)
- `hunt_triage.tsv` — GSE hits plus locked PACIFIC/ADRIA triage
- `sample_table_GSE253564.tsv`, `sample_table_GSE248378.tsv`
- `summary.json`, `provenance.json`
- `figures/GSE253564_TACSTD2_MPR_residuals.png`
- `figures/TACSTD2_partial_rho_CYT_CD8.png`
- `figures/GSE248378_TACSTD2_recurrence.png`

## Public extras (added on top of A2)

The user durvalumab result is taken as given. Additional leftover public lung IO RNA *besides* GSE253564 — TACSTD2/CLDN4 vs MPR/DCB/ORR and vs CD8/GEP after an ESTIMATEScore residual — is in `extra/` (`REPORT.md`, `paper_tables.md`, four figures). GSE253564 is not re-scored there.

## Data

- GSE253564 / GSE248378 FPKM and series matrices (NCBI GEO FTP)
- ESTIMATE `SI_geneset.gmt` v1.0.13
- PMID 38401548 Table S1 (`mmc2.xlsx` via Europe PMC OA)
- PMID 38114518 source data (`41467_2023_44195_MOESM6_ESM.xlsx`)
- Extra leftover matrices under `data/A2_wave2/extra/` (GSE207422, GSE126044, GSE166449, GSE135222, GSE329813)

