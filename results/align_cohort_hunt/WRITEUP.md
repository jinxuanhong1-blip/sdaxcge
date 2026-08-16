# ALIGN cohort hunt: exact datasets and purity-aware TACSTD2 correlations

## Bottom line

The public cohorts can be identified without inventing accessions:

1. **OncoSG LUAD** is Chen et al., *Nature Genetics* 2020, OncoSG dataset **GIS031**, mirrored in cBioPortal as **`luad_oncosg_2020`**. Raw RNA FASTQs are controlled at **EGAS00001002941 / EGAD00001004421**; open processed RNA and sample purity are available from OncoSG/cBioPortal.
2. The durvalumab NSCLC slide clue **“TACSTD2 vs immune, ρ≈−0.65”** resolves to **GSE248378**, the post-neoadjuvant resection RNA-seq matrix from **NCT02904954** (durvalumab alone or durvalumab + SBRT). Using the canonical 141-gene ESTIMATE immune signature reproduces **Spearman ρ = −0.6576, p = 1.06×10⁻⁴, n = 29**. The companion pretreatment accession is **GSE253564**, not PACIFIC.
3. **GSE76628** is real, but it is **not KL lung cancer**. It is a 78-sample Ad-VEGF-A164 flank-tissue time course in **athymic nude mice**, with DC101/G6 interventions, from a gastric-cancer stromal-signature paper. It can be analyzed alone, but it must not be relabeled as a KRAS/LKB1 tumor cohort.
4. No open PACIFIC tumor-RNA accession was found. PACIFIC is **NCT02125461**; its publications direct researchers to AstraZeneca's request-based data-sharing process.

No slide file was present in the repository, so “exact” means an accession-level match to the cohort descriptions and a numerical reproduction of the reported correlation.

## Downloaded processed data

All files are below 2 GB and stored in `data/`:

| Cohort | Local file | Public content |
|---|---|---|
| OncoSG GIS031 | `OncoSG_GIS031_selected_expression_zscores.tsv` | cBioPortal RNA-seq V2 RSEM all-sample z-scores for TACSTD2 and prespecified signature genes, 169 tumors |
| OncoSG GIS031 | `OncoSG_GIS031_sample_purity.tsv` | Public `PURITY` clinical field, 302 samples (169 overlap RNA) |
| GSE76628 | `GSE76628_series_matrix.txt.gz` | Full processed MAS5 array matrix and sample metadata, 78 samples |
| GPL1261 | `GPL1261-55999.txt.gz` | Mouse 430 2.0 probe-to-gene annotation |
| GSE110390 | `GSE110390_CP1108_NSCLC_21gene_expression.txt.gz` | Public 21-gene NSCLC matrix, 97 tumors |
| GSE110390 | `GSE110390_series_matrix.txt.gz` | Series metadata |
| GSE248378 | `GSE248378_Durva_Post_FPKMs.txt.gz` | Post-treatment bulk RNA FPKM, 29 tumors |
| GSE253564 | `GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz` | Pretreatment bulk RNA FPKM, 32 tumors |

The OncoSG downloader records every API URL and the exact molecular profile in `provenance/oncosg_api_manifest.json`. `provenance/sha256sums.txt` provides checksums.

## Analysis

### Methods

- Spearman correlations use pairwise-complete samples.
- CD8 score: mean per-gene z-score of available `CD8A`, `CD8B`, `CCL5`, `GZMK`, `LCK`, `TRAC` orthologs.
- NK score: mean per-gene z-score of available `NKG7`, `KLRD1`, `GNLY`, `PRF1`, `CTSW`, `NCR1` orthologs.
- Pan-immune score: mean per-gene z-score of available `PTPRC`, `CD3D`, `CD3E`, `LST1`, `CD74`, `HLA-DRA` orthologs.
- ESTIMATE ImmuneScore is canonical ssGSEA using the 141-gene `ImmuneSignature` from ESTIMATE. The exact GMT is archived in `provenance/ESTIMATE_SI_geneset.gmt`.
- “Purity-aware” is a partial Spearman correlation: rank-transform target, score, and covariates; regress target and score ranks on covariates; correlate residuals.
- OncoSG uses its **reported sample purity**. GSE248378/GSE253564 have no public measured-purity field, so the adjustment uses a prespecified epithelial-expression proxy (`EPCAM`, `KRT7/8/18/19`, `MUC1`). This is labeled a proxy and is not misrepresented as ABSOLUTE or histologic purity.
- GSE76628 has no tumor compartment, so tumor purity is undefined. Its adjusted model uses an epithelial-content proxy plus experimental day and treatment indicators.

Full machine-readable results are in `analysis/correlations.tsv`.

### OncoSG LUAD (GIS031)

| Outcome | Raw ρ | p | Purity-adjusted partial ρ | p | n |
|---|---:|---:|---:|---:|---:|
| CD8 score | −0.378 | 4.18×10⁻⁷ | −0.304 | 5.98×10⁻⁵ | 169 |
| NK score | −0.477 | 5.47×10⁻¹¹ | −0.428 | 6.34×10⁻⁹ | 169 |
| Pan-immune score | −0.312 | 3.72×10⁻⁵ | −0.229 | 0.00271 | 169 |

The inverse TACSTD2–immune association weakens after measured-purity adjustment but remains statistically clear.

### Durvalumab NSCLC

**GSE248378 post-treatment is the numerical slide match.**

| Outcome | Raw ρ | p | Epithelial-proxy-adjusted partial ρ | p | n |
|---|---:|---:|---:|---:|---:|
| ESTIMATE ImmuneScore | **−0.658** | **1.06×10⁻⁴** | **−0.594** | **6.75×10⁻⁴** | 29 |
| CD8 score | −0.784 | 4.81×10⁻⁷ | −0.783 | 5.08×10⁻⁷ | 29 |
| NK score | −0.823 | 4.35×10⁻⁸ | −0.807 | 1.23×10⁻⁷ | 29 |
| Pan-immune score | −0.560 | 0.00158 | −0.517 | 0.00406 | 29 |

The −0.65 value is therefore reproducible from the open post-treatment matrix with the established ESTIMATE immune score. It is not supportable from GSE110390's open processed matrix because that matrix contains only 21 immune genes and omits TACSTD2.

The same NCT02904954 trial's **pretreatment** matrix, GSE253564, gives ESTIMATE immune ρ = −0.712 (p = 4.91×10⁻⁶; proxy-adjusted partial ρ = −0.710, p = 5.28×10⁻⁶; n = 32). It is useful corroboration but is not the −0.65 slide value.

### GSE76628 alone

| Outcome | Raw ρ | p | Adjusted partial ρ | p | n |
|---|---:|---:|---:|---:|---:|
| CD8 score | −0.226 | 0.0471 | −0.254 | 0.0248 | 78 |
| NK score | −0.334 | 0.00279 | −0.179 | 0.116 | 78 |
| Pan-immune score | −0.141 | 0.218 | −0.158 | 0.167 | 78 |

Adjustment includes epithelial content, collection day, and anti-VEGFR treatment. The raw Tacstd2–NK association does not survive this design-aware adjustment.

The CD8 result requires an especially strong caveat: these are **Foxn1-null athymic nude mice**, so conventional T-cell biology is abnormal, and the samples are flank tissue rather than tumors. A “CD8 signature” here is an expression summary, not evidence of ordinary CD8 T-cell infiltration. These data cannot validate a KL-lung mechanism.

## PACIFIC and other open durvalumab lung RNA

See `durvalumab_rna_catalog.tsv`.

- **PACIFIC / NCT02125461:** no public sample-level RNA accession identified; request-based sponsor data only.
- **GSE110390:** 97 advanced NSCLC biopsies from CP1108/NCT01693562, but only a 21-gene processed matrix is open and TACSTD2 is absent.
- **GSE131933:** open longitudinal NSCLC scRNA-seq after durvalumab or durvalumab+tremelimumab; mechanistic and very small, not PACIFIC.
- **GSE248378/GSE253564:** open post/pre bulk RNA from NCT02904954; the relevant sample-level matrices.
- **NeoCOAST / NCT03794544:** open aggregate differential-expression supplement, but no comprehensive open sample-level matrix or RNA accession found.

## Reproduce

```bash
python3 results/align_cohort_hunt/fetch_oncosg.py
MPLBACKEND=Agg python3 results/align_cohort_hunt/analyze.py
sha256sum results/align_cohort_hunt/data/* > results/align_cohort_hunt/provenance/sha256sums.txt
```

Versions used in this run: Python 3.12, pandas 3.0.5, SciPy 1.18.0, matplotlib 3.11.1.
