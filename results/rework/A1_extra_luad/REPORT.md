# EXTRA LUAD RNA tables — TACSTD2 / CLDN4 vs CD8 / GEP after purity residual

The user's **TCGA + OncoSG** purity-corrected TACSTD2–immune result is **taken as given**. OncoSG numbers below are from PR #139 and were not re-run. This folder **adds** public non-TCGA / East-Asian LUAD RNA cohorts.

## Extra table (this run)

| Cohort | Gene | Feature | n | Unadj ρ | Unadj p | Partial ρ | Partial p | Purity |
|---|---|---|---:|---:|---:|---:|---:|---|
| GSE31210 | TACSTD2 | CD8 | 226 | -0.289 | 9.76e-06 | -0.107 | 0.108 | ESTIMATE TumorPurity (ssGSEA Yoshihara 141+141) |
| GSE31210 | TACSTD2 | GEP18 | 226 | -0.248 | 0.00017 | -0.013 | 0.846 | ESTIMATE TumorPurity (ssGSEA Yoshihara 141+141) |
| GSE31210 | CLDN4 | CD8 | 226 | -0.341 | 1.45e-07 | -0.127 | 0.0574 | ESTIMATE TumorPurity (ssGSEA Yoshihara 141+141) |
| GSE31210 | CLDN4 | GEP18 | 226 | -0.287 | 1.13e-05 | -0.005 | 0.94 | ESTIMATE TumorPurity (ssGSEA Yoshihara 141+141) |
| GSE72094 | TACSTD2 | CD8 | 442 | -0.227 | 1.39e-06 | -0.085 | 0.0741 | ESTIMATE TumorPurity (ssGSEA Yoshihara 141+141) |
| GSE72094 | TACSTD2 | GEP18 | 442 | -0.216 | 4.47e-06 | -0.046 | 0.331 | ESTIMATE TumorPurity (ssGSEA Yoshihara 141+141) |
| GSE72094 | CLDN4 | CD8 | 442 | -0.219 | 3.35e-06 | -0.051 | 0.282 | ESTIMATE TumorPurity (ssGSEA Yoshihara 141+141) |
| GSE72094 | CLDN4 | GEP18 | 442 | -0.194 | 4.17e-05 | 0.023 | 0.637 | ESTIMATE TumorPurity (ssGSEA Yoshihara 141+141) |
| GSE68465 | TACSTD2 | CD8 | 443 | -0.152 | 0.00138 | -0.170 | 0.000324 | ESTIMATE TumorPurity (ssGSEA Yoshihara 141+141) |
| GSE68465 | TACSTD2 | GEP18 | 443 | -0.120 | 0.0114 | -0.164 | 0.000556 | ESTIMATE TumorPurity (ssGSEA Yoshihara 141+141) |
| GSE68465 | CLDN4 | CD8 | 443 | -0.177 | 0.000174 | -0.092 | 0.054 | ESTIMATE TumorPurity (ssGSEA Yoshihara 141+141) |
| GSE68465 | CLDN4 | GEP18 | 443 | -0.215 | 4.72e-06 | -0.113 | 0.0173 | ESTIMATE TumorPurity (ssGSEA Yoshihara 141+141) |
| GSE19804 | TACSTD2 | CD8 | 60 | -0.024 | 0.854 | 0.196 | 0.138 | ESTIMATE TumorPurity (ssGSEA Yoshihara 141+141) |
| GSE19804 | TACSTD2 | GEP18 | 60 | -0.106 | 0.421 | 0.108 | 0.413 | ESTIMATE TumorPurity (ssGSEA Yoshihara 141+141) |
| GSE19804 | CLDN4 | CD8 | 60 | -0.465 | 0.000185 | -0.274 | 0.0356 | ESTIMATE TumorPurity (ssGSEA Yoshihara 141+141) |
| GSE19804 | CLDN4 | GEP18 | 60 | -0.287 | 0.0262 | 0.083 | 0.532 | ESTIMATE TumorPurity (ssGSEA Yoshihara 141+141) |
| CPTAC_LUAD_RNA | TACSTD2 | CD8 | 108 | -0.206 | 0.0306 | -0.188 | 0.053 | published WES_purity (Gillette 2020 freeze) |
| CPTAC_LUAD_RNA | TACSTD2 | CD8 | 110 | -0.206 | 0.0306 | -0.284 | 0.00275 | published ESTIMATE_ESTIMATEScore (RNA impurity axis) |
| CPTAC_LUAD_RNA | TACSTD2 | GEP18 | 108 | 0.015 | 0.88 | 0.042 | 0.671 | published WES_purity (Gillette 2020 freeze) |
| CPTAC_LUAD_RNA | TACSTD2 | GEP18 | 110 | 0.015 | 0.88 | -0.046 | 0.638 | published ESTIMATE_ESTIMATEScore (RNA impurity axis) |
| CPTAC_LUAD_RNA | CLDN4 | CD8 | 108 | -0.260 | 0.00606 | -0.226 | 0.019 | published WES_purity (Gillette 2020 freeze) |
| CPTAC_LUAD_RNA | CLDN4 | CD8 | 110 | -0.260 | 0.00606 | -0.282 | 0.00293 | published ESTIMATE_ESTIMATEScore (RNA impurity axis) |
| CPTAC_LUAD_RNA | CLDN4 | GEP18 | 108 | -0.008 | 0.933 | 0.039 | 0.693 | published WES_purity (Gillette 2020 freeze) |
| CPTAC_LUAD_RNA | CLDN4 | GEP18 | 110 | -0.008 | 0.933 | 0.044 | 0.648 | published ESTIMATE_ESTIMATEScore (RNA impurity axis) |

## OncoSG reference (already done)

| Cohort | Gene | Feature | n | Unadj ρ | Partial ρ | Partial p | Purity |
|---|---|---|---:|---:|---:|---:|---|
| OncoSG PR #139 | TACSTD2 | CD8A | 169 | −0.380 | −0.309 | 4.69e-05 | published clinical PURITY |
| OncoSG PR #139 | TACSTD2 | GEP18 (17/18) | 169 | −0.414 | −0.349 | 3.62e-06 | published clinical PURITY |

## Cohorts

| Cohort | Region | Platform | Tumor n | Purity used |
|---|---|---|---:|---|
| GSE31210 Okayama / Kohno | East Asia (Japan) | GPL570 U133 Plus 2.0 | 226 | ESTIMATE TumorPurity (ssGSEA, Yoshihara 141+141) |
| GSE19804 Taiwan never-smoker female | East Asia (Taiwan) | GPL570 U133 Plus 2.0 | 60 | ESTIMATE TumorPurity (ssGSEA) |
| GSE72094 | USA | GPL15048 HuRSTA | 442 | ESTIMATE TumorPurity (ssGSEA) |
| GSE68465 Director's Challenge | USA / multi-site | GPL96 U133A | 443 | ESTIMATE TumorPurity (ssGSEA) |
| CPTAC LUAD RNA (Gillette *Cell* 2020) | USA | RNA-seq RSEM log2 UQ | 110 | published WES_purity **and** published ESTIMATEScore |

Primary tumors only. GSE19804 paired normals dropped. GSE31210 restricted to `tissue: primary lung tumor`.

**Coverage notes.** GSE68465 GEP18 is **15/18** (`CD274`, `CD276`, `TIGIT` absent on U133A). Other GEO cohorts and CPTAC are 18/18. CPTAC WES complete-case n=108 (2 tumors lack `WES_purity`). ESTIMATE gene coverage on GPL570 is 136/141 stromal and 138/141 immune.

## Definitions

- **TACSTD2 / CLDN4 / CD8:** microarray = max-mean probe collapse to HUGO; CPTAC = Ensembl `ENSG00000184292` / `ENSG00000189143` / `ENSG00000153563` (version suffix stripped).
- **CD8** = `CD8A`.
- **GEP18** = unweighted within-cohort z-mean of the Ayers 2017 18-gene list. Missing genes are not imputed (coverage in `cohort_coverage.tsv`).
- **GEO purity:** ESTIMATE ssGSEA (Barbie/GSVA tau=0.25, ranks 1..10000) on Yoshihara 2013 Stromal141 + Immune141, then `TumorPurity = cos(0.6049872018 + 0.0001467884 × ESTIMATEScore)`. Gene list: `data/signatures/estimate_yoshihara_2013.tsv` (Supp Data 1).
- **CPTAC purity:** published `WES_purity` (DNA). Published `ESTIMATE_ESTIMATEScore` is a second RNA axis. The Yoshihara cosine is **not** applied to the freeze scores (they wrap; see PR #99).
- **Partial Spearman:** Pearson of rank residuals; df = n − 3.

## 2023–2026 GEO hunt

NCBI GEO (GDS) search `lung adenocarcinoma[Title] AND 2023:2026[Publication Date] AND Expression profiling by high throughput sequencing` returned 122 records. The first page is cell-line, n=1–36, scRNA, or plasma EV studies (GSE293914 n=1; GSE330179 n=18; GSE332750 n=4; GSE299604 n=11; GSE317139 n=36). GSE226481 (China, early- vs late-onset LUAD) has **n=14** tumors — too small for a purity residual. GSE40419 (Seo Korea RNA-seq) and GSE140343 (paired LUAD RNA-seq) have **no processed series matrix** (SRA/GTF only). No 2023–2026 GEO LUAD bulk tumor matrix with n≥50 was used.

## Files

- `EXTRA_TABLE_tacstd2_cldn4_cd8_gep.tsv` — all extra n/ρ/p
- `extra_correlations.tsv` — same
- `cohort_coverage.tsv`
- `oncosg_reference_pr139.tsv`
- `samples_*.tsv`
- `figures/bar_extra_partial_rho.png`
- `summary.json` / `provenance.json`

## Reproduce

```
pip install -r requirements.txt
python scripts/rework_A1_extra_luad.py
```

Public GEO series matrices + CPTAC S3 freeze v1.2 + Yoshihara Supp Data 1. Raw downloads under `data/extra_luad/` (gitignored).
