# GSE81089 Swedish NSCLC RNA-seq — CLDN4 vs CD8A / ImmuneScore / CD274 / IFN-γ / MHC-I

**Additive CLDN4-only.** User thesis (given, not re-audited): CLDN4-high is a tight-junction / epithelial programme. No slide was re-scored. No TACSTD2∩CLDN4 dual-high.

Public Uppsala / SciLifeLab resected NSCLC RNA-seq (Djureinovic et al., *JCI Insight* 2016, PMID 29282718; GEO [GSE81089](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81089)). Unit is the **patient / tumor array** (one tumor GSM per patient). Processed FPKM (Cufflinks 2.1.1, Ensembl v73 GRCh37) and featureCounts 1.4.0-p1 are on the series record. No FASTQ.

Primary histology codes are the GEO overall_design integers: `1` = squamous cell cancer (LUSC), `2` = AC unspecified (LUAD), `3` = large cell / NOS. Mixed NSCLC is not a LUAD finding.

## Honest n

| Item | Public? | n | Note |
|---|---|---:|---|
| GSM / FPKM columns in the series | yes | **218** | 199 tumor + 19 paired normal |
| **Tumor patients (primary NSCLC n)** | yes | **199** | one tumor array per patient; **this is the n used below** |
| Paired normal arrays | yes | 19 | inventoried; not mixed into tumor Spearman |
| **histology = 2 AC unspecified (LUAD)** | yes | **108** | GEO code 2; this is the LUAD n |
| **histology = 1 squamous (LUSC)** | yes | **67** | GEO code 1; this is the LUSC n |
| histology = 3 large cell / NOS | yes | 24 | **not** folded into LUAD or LUSC; n&lt;40 |
| Stage (pTNM code) | yes | 199 | 1=1a … 7=IV; not a covariate |
| OS / vital date / dead | yes | 199 | deposited; not this claim |
| ICI / treatment | no | 0 | resected 2006–2010 atlas |
| Tumor % / ABSOLUTE purity | no | 0 | only public proxy is an RNA epithelial score |
| CLDN4 `ENSG00000189143` finite (tumors) | yes | **199** | Cufflinks FPKM |
| CD8A / CD274 / TACSTD2 finite (tumors) | yes | **199** | TACSTD2 is companion only |
| ESTIMATE ImmuneSignature genes present | yes | **141 / 141** | Yoshihara 2013; mean-z, not official ssGSEA |
| Ayers IFN-γ 6-gene / MHC-I 6-gene | yes | **6 / 6** | IFNG STAT1 CXCL9 CXCL10 IDO1 HLA-DRA · HLA-A/B/C B2M TAP1 TAP2 |

Do not write n=218 for the tumor tests. Do not write n=199 for LUAD (108) or LUSC (67).

## One-row table

| dataset | histology | n | CLDN4–CD8A ρ (p) | residual | CLDN4–CD274 ρ (p) | residual | CLDN4–ImmuneScore ρ (p) | residual | CLDN4–IFN-γ ρ (p) | residual | CLDN4–MHC-I ρ (p) | residual |
|---|---|---:|---|---|---|---|---|---|---|---|---|---|
| GSE81089 Uppsala RNA-seq | NSCLC tumors | **199** | **−0.102 (0.15)** | −0.038 (0.60) | **−0.103 (0.15)** | −0.054 (0.45) | **+0.118 (0.097)** | +0.179 (0.012) | **−0.017 (0.81)** | +0.018 (0.80) | **+0.152 (0.032)** | +0.125 (0.079) |
| GSE81089 | LUAD (`AC unspecified`) | **108** | −0.011 (0.91) | +0.068 (0.49) | +0.005 (0.96) | +0.026 (0.79) | +0.133 (0.17) | +0.225 (0.020) | +0.098 (0.31) | +0.119 (0.22) | +0.168 (0.083) | +0.142 (0.14) |
| GSE81089 | LUSC (`squamous`) | **67** | −0.244 (0.047) | −0.203 (0.10) | +0.016 (0.90) | +0.047 (0.71) | −0.075 (0.55) | +0.011 (0.93) | −0.147 (0.24) | −0.139 (0.27) | +0.001 (0.99) | +0.016 (0.90) |

Full numeric row: `tables/one_row.tsv`. Residual = partial Spearman on the pan-epithelial mean-z (EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7; 6/6 present).

## Locked design (before ρ)

| Item | Choice |
|---|---|
| Matrix | Public Cufflinks FPKM. No FASTQ / SRA. |
| Transform | `log2(FPKM+1)` |
| Sensitivity | featureCounts `log2(CPM+1)` on the same GSM |
| Unit | tumor GSM = patient |
| ImmuneScore | Yoshihara 2013 ESTIMATE ImmuneSignature **mean-z** (141/141). Not official R ssGSEA |
| IFN-γ | Ayers 2017 6-gene mean-z |
| MHC-I | HLA-A / HLA-B / HLA-C / B2M / TAP1 / TAP2 mean-z |
| Epithelial residual | KRT / EPCAM / CDH1 mean-z |
| Extra | CLDN4 Q4 vs Q1 MWU; LUAD vs LUSC split |
| Companion | TACSTD2 vs the same axes (same-run only). **No dual-high** |
| HOLDS (same as PR 229 / PR 313) | n≥40, ρ_adj&lt;0, p_adj&lt;0.05 |

## CLDN4 vs partners (primary, NSCLC tumors n=199)

Spearman is rank-based. Bootstrap 95% CI, 2,000 resamples, seed `20260817`.

| pair | n | ρ | 95% CI | p | ρ_adj | p_adj |
|---|---:|---:|---|---:|---:|---:|
| CLDN4 vs **CD8A** | **199** | **−0.102** | −0.232 to +0.040 | 0.15 | −0.038 | 0.60 |
| CLDN4 vs **CD274** | 199 | **−0.103** | −0.243 to +0.040 | 0.15 | −0.054 | 0.45 |
| CLDN4 vs **ImmuneScore** | 199 | **+0.118** | −0.026 to +0.256 | 0.097 | +0.179 | 0.012 |
| CLDN4 vs **IFN-γ Ayers6** | 199 | **−0.017** | −0.151 to +0.126 | 0.81 | +0.018 | 0.80 |
| CLDN4 vs **MHC-I** | 199 | **+0.152** | +0.009 to +0.285 | 0.032 | +0.125 | 0.079 |
| CLDN4 Q4 vs Q1 on CD8A | 50 vs 50 | — | — | MWU 0.034 | — | — |
| CLDN4 Q4 vs Q1 on CD274 | 50 vs 50 | — | — | MWU 0.19 | — | — |
| CLDN4 Q4 vs Q1 on ImmuneScore | 50 vs 50 | — | — | MWU 0.24 | — | — |
| CLDN4 Q4 vs Q1 on IFN-γ | 50 vs 50 | — | — | MWU 0.52 | — | — |
| CLDN4 Q4 vs Q1 on MHC-I | 50 vs 50 | — | — | MWU 0.071 | — | — |

CD8A / CD274 / IFN-γ stay near zero after the epithelial residual. Crude MHC-I is a small positive and loses p&lt;0.05 after residual. ImmuneScore residual is a small positive (CLDN4-high is not ImmuneScore-low once epithelium is removed). Under the locked HOLDS rule, no residual immune-low pair holds.

Q4 vs Q1 CD8A (50 vs 50, rank-biserial −0.25) is the only extra cut with MWU p&lt;0.05 on mixed tumors; it is not a residual Spearman HOLDS and is not quoted as LUAD.

## Histology split (honest n)

Do not quote the mixed row as a LUAD result. Large-cell n=24 is underpowered and is in `tables/spearman_cldn4_vs_partners.tsv` only.

| cohort | n | CLDN4–CD8A ρ (p) / residual | CLDN4–CD274 | CLDN4–ImmuneScore | CLDN4–IFN-γ | CLDN4–MHC-I |
|---|---:|---|---|---|---|---|
| NSCLC tumors | **199** | −0.102 (0.15) / −0.038 (0.60) | −0.103 (0.15) / −0.054 (0.45) | +0.118 (0.097) / +0.179 (0.012) | −0.017 (0.81) / +0.018 (0.80) | +0.152 (0.032) / +0.125 (0.079) |
| LUAD AC unspecified | **108** | −0.011 (0.91) / +0.068 (0.49) | +0.005 (0.96) / +0.026 (0.79) | +0.133 (0.17) / +0.225 (0.020) | +0.098 (0.31) / +0.119 (0.22) | +0.168 (0.083) / +0.142 (0.14) |
| LUSC squamous | **67** | −0.244 (0.047) / −0.203 (0.10) | +0.016 (0.90) / +0.047 (0.71) | −0.075 (0.55) / +0.011 (0.93) | −0.147 (0.24) / −0.139 (0.27) | +0.001 (0.99) / +0.016 (0.90) |

LUSC CD8A is the only crude p&lt;0.05 histology cell (n=67, CI crosses 0). Residual p=0.10. LUAD n=108 is null on CD8A / CD274 / IFN-γ.

## What CLDN4 *does* track here

| pair | n | ρ | 95% CI | p | ρ_adj |
|---|---:|---:|---|---:|---:|
| CLDN4 vs epithelial mean-z | 199 | **+0.629** | 0.531 to 0.712 | 2.6×10⁻²³ | — |
| CLDN4 vs TACSTD2 (companion) | 199 | **+0.271** | 0.138 to 0.396 | 1.1×10⁻⁴ | +0.157 (0.027) |
| CLDN4 vs epithelial (LUAD) | 108 | +0.497 | 0.336 to 0.632 | 4.5×10⁻⁸ | — |
| CLDN4 vs TACSTD2 (LUAD) | 108 | +0.472 | 0.310 to 0.612 | 2.5×10⁻⁷ | +0.306 (0.0013) |
| CLDN4 vs epithelial (LUSC) | 67 | +0.567 | 0.374 to 0.720 | 5.8×10⁻⁷ | — |

CLDN4 sits on the epithelial / TACSTD2 side of this RNA-seq. That is the same direction as a tight-junction / tumour-cell programme. TACSTD2 is a companion column only.

## featureCounts sensitivity (same 199 tumors)

`log2(CPM+1)` from the deposited featureCounts matrix. Same genes, same GSM.

| pair | n | ρ (p) | ρ_adj (p) |
|---|---:|---|---|
| CLDN4 vs CD8A | 199 | −0.091 (0.20) | −0.046 (0.52) |
| CLDN4 vs CD274 | 199 | −0.091 (0.20) | −0.050 (0.48) |
| CLDN4 vs ImmuneScore | 199 | +0.115 (0.11) | +0.168 (0.018) |
| CLDN4 vs IFN-γ | 199 | −0.009 (0.90) | +0.015 (0.83) |
| CLDN4 vs MHC-I | 199 | +0.090 (0.21) | +0.081 (0.26) |
| CLDN4 vs epithelial | 199 | +0.649 (3.2×10⁻²⁵) | — |

Quantification choice does not move CD8A / CD274 / IFN-γ. Crude MHC-I p=0.032 on FPKM is not reproduced on CPM (p=0.21).

## What this does not test

- ICI response (not an ICI series).
- Official ESTIMATE ssGSEA ImmuneScore from the R package (this is ImmuneSignature mean-z).
- Pathologist CD8 / PD-L1 IHC.
- A LUAD-wide CLDN4–immune law. This is one public Swedish HiSeq 2500 series.
- OS (labels are on GEO; not scored here).
- TACSTD2∩CLDN4 dual-high (out of scope).

## Reproduce

```bash
python3 -m pip install -r methods/gse81089_cldn4/requirements.txt
python3 methods/gse81089_cldn4/analyze.py
```

Downloads (not committed) go to `$GSE81089_CLDN4_DATA` (default `/tmp/gse81089_cldn4`):

- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81089/matrix/GSE81089_series_matrix.txt.gz`
- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81089/suppl/GSE81089_FPKM_cufflinks.tsv.gz`
- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81089/suppl/GSE81089_readcounts_featurecounts.tsv.gz`
- `https://ftp.ensembl.org/pub/release-73/gtf/homo_sapiens/Homo_sapiens.GRCh37.73.gtf.gz` (symbol map)
- `https://download.r-forge.r-project.org/src/contrib/estimate_1.0.13.tar.gz` (ImmuneSignature GMT only)

## Files

- `analyze.py` — download, Ensembl map, Spearman / partial Spearman, figures
- `tables/one_row.tsv` — n / ρ / p / residual
- `tables/label_inventory.tsv` — public vs missing fields; honest n
- `tables/histology_counts.tsv` — GEO histology codes
- `tables/gene_id_confirm.tsv` — named Ensembl v73 IDs
- `tables/immunescore_gene_coverage.tsv` — 141/141 ImmuneSignature genes
- `tables/spearman_cldn4_vs_partners.tsv` — NSCLC + LUAD + LUSC + large-cell; FPKM + CPM
- `tables/highlow_cldn4.tsv` — CLDN4 Q4 vs Q1
- `tables/sample_annotation.tsv` — 218 GSM
- `tables/summary.json`
- `figures/fig1_cldn4_vs_partners.png`
- `figures/fig2_cldn4_forest.png`
- `figures/fig3_forest_luad_lusc.png`
- `figures/fig4_q4q1_boxplots.png`
