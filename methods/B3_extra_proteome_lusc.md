# Methods — B3 extra: CPTAC LUAD / remaining open proteome, and TCGA-LUSC

This extra is **additive** to user B3. The TCGA-LUAD RNA result (TJ-high tumors have lower CD8 / GEP) is taken as given and is not re-estimated here. CPTAC LSCC CLDN4 protein versus ImmuneScore is already on the public record and is not re-run as an audit.

Question asked here:

1. In **CPTAC LUAD** protein, what are the Spearman associations of **CLDN4** and a **structural TJ protein score** with ESTIMATE ImmuneScore, CD8, and related ESTIMATE axes?
2. Do any **other open CPTAC freeze proteomes** (same public S3 layout, LSCC excluded) show the same protein-versus-ImmuneScore numbers?
3. In **TCGA-LUSC** RNA, using the same TJ / CD8 / GEP definitions as B3, what are n / ρ / p for TJ-high versus CD8 and GEP?

Public files only. No ICI labels are present in CPTAC surgical cohorts; OS/PFS are not used.

## Data

### CPTAC LUAD (primary protein extra)

Open CPTAC pan-cancer freeze v1.2 (LinkedOmics / AWS S3; HTTP 200):

`https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LUAD/`

| File | Role |
|---|---|
| `LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt` | TMT gene abundance, already log2 vs reference |
| `LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt` | matched tumor RNA (for CD8 / GEP / RNA TJ) |
| `LUAD_phenotype.txt` | freeze ESTIMATE, CIBERSORT, xCell, WES purity |

Cohort: treatment-naive surgical LUAD, Gillette et al., *Cell* 2020. Sample IDs are harmonized (`C3L.00001` = `C3L-00001`). Protein values are **not** logged again.

### Remaining open proteomes

The same freeze prefix exposes protein + phenotype for BRCA, COAD, CCRCC, GBM, HNSCC, OV, PDAC, and UCEC (HEAD HTTP 200). Each is scored with the same CLDN4 / TJ-15 protein rules versus whatever ImmuneScore / CD8 columns exist in that phenotype table. **LSCC is not downloaded.**

Other lung proteomes (APOLLO / Soltis, Xu early-stage LUAD, raw MassIVE MSV000086793) were not used: they are not compact open gene-abundance + ImmuneScore tables on this public prefix.

### TCGA-LUSC (histology extra)

| File | Source |
|---|---|
| `TCGA.LUSC.sampleMap/HiSeqV2` | UCSC Xena S3, `log2(norm_count+1)` |
| ESTIMATE RNAseqV2 | MD Anderson / Yoshihara 2013 (`lung_squamous_cell_carcinoma_RNAseqV2.txt`) |

Primary solid tumors only (barcode `*-01`) present in both matrices.

## Gene sets (locked)

Same structural TJ used for B3 (PR #69 / PR #90). Not a search over claudins.

- **TJ-15:** CLDN1, CLDN3, CLDN4, CLDN7, OCLN, TJP1, TJP2, TJP3, F11R, JAM2, JAM3, MARVELD2, MARVELD3, CGN, CGNL1.
- **TJ-7** (claim-page module, sensitivity): CLDN1, CLDN4, CLDN7, F11R, TJP1, TJP2, OCLN.
- **CD8:** CD8A, CD8B.
- **CYT:** GZMA, PRF1 (Rooney).
- **GEP18:** Ayers 2017 18-gene T-cell-inflamed list.

CPTAC matrices are Ensembl-indexed. Symbols are mapped to GENCODE Ensembl IDs (`scripts/b3_extra_proteome_lusc/genes.py`) and matched by ID prefix (version suffix allowed). Xena LUSC is symbol-indexed.

## Scores

- **Mean z-score:** each gene is z-scored across samples in that matrix; the sample score is the mean of available z-values.
- TJ-15 requires ≥5 member genes with ≥8 non-missing samples; TJ-7 requires ≥4; GEP18 requires ≥10; CD8/CYT require both genes.
- **CLDN4 protein** is the single TMT row. Missing values are left missing (MNAR-typical for a small tetraspan). No RNA imputation. In CPTAC LUAD, CLDN1 protein is observed in 9 tumors and CLDN7 in 17; they enter the TJ-15 mean only where measured.
- CPTAC immune endpoints are the freeze phenotype columns: ESTIMATE ImmuneScore / StromalScore, CIBERSORT CD8, xCell CD8, xCell immune score, WES purity when present.
- CPTAC freeze ESTIMATE scores are **not** converted with the Yoshihara cosine purity formula (those scores sit outside the Affymetrix calibration range). WES purity is the DNA residual covariate when present.
- TCGA-LUSC purity = `cos(0.6049872018 + 0.0001467884 · ESTIMATEScore)` (Yoshihara 2013), same formula as B3.

## Statistics

- **Spearman ρ** on pairwise-complete samples. Two-sided p from SciPy.
- **95% CI:** percentile bootstrap, 2,000 resamples, seed `20260816`.
- **Partial Spearman:** rank-transform X, Y, and the covariate; residualize X and Y ranks on the covariate rank; Spearman of residuals (n−3 implicit in the residual correlation).
- **TJ-high vs TJ-low:** cohort-internal median split of the TJ (or CLDN4) score; two-sided Mann–Whitney U on the immune endpoint.
- Remaining-proteome rows are descriptive extras. They are not pooled with LUAD and are not an LSCC audit.

## What is not done

- No re-fit of TCGA-LUAD B3.
- No CPTAC LSCC protein re-analysis.
- No ICI response model (these surgical proteomes have no ICI labels).
- No gene-set fishing to enlarge |ρ|.

## Reproduce

```bash
python3 -m pip install -r requirements.txt
python3 scripts/b3_extra_proteome_lusc/download.py
python3 scripts/b3_extra_proteome_lusc/analyze.py
```

Outputs: `results/B3_extra_proteome_lusc/` (tables, figures, `RESULTS.md`).
