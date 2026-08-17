# GSE37745 NSCLC array — CLDN4 vs CD8A / CD274

**Additive only.** Public Botling / Uppsala resected NSCLC Affymetrix U133 Plus 2.0 series (Botling et al., *Clin Cancer Res* 2013, PMID 23032747; GEO [GSE37745](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37745)). Unit is the **array**. One titled patient per array. No slide was re-scored. No ICI arm.

Primary pairs are **CLDN4 vs CD8A** and **CLDN4 vs CD274**. Histology is mixed NSCLC and is not collapsed into one law.

## Honest n

| Item | Public? | n | Note |
|---|---|---:|---|
| Arrays in the series matrix | yes | **196** | 54,675 probes × 196 GSM; **0** missing RMA values |
| Unique GSM / unique titles | yes | **196** | all unique |
| Unique patient IDs from title | yes | **196** | `Patient N`; **1 array / patient**; no duplicate pair |
| Patients (GEO design text) | yes | **196** | “196 consecutive NSCLC patients, operated between 1995 and 2005” |
| Histology: adeno (LUAD) | yes | **106** | GEO `histology: adeno` |
| Histology: squamous (LUSC) | yes | **66** | GEO `histology: squamous` |
| Histology: large (LCC) | yes | **24** | GEO `histology: large`; kept in the NSCLC pool; **UNDERPOWERED** alone (n&lt;40) |
| Stage | yes | 196 | 1a:40, 1b:90, 2a:6, 2b:29, 3a:21, 3b:6, 4:4 |
| OS (dead + days) | yes | **196** | 145 dead / 51 alive; days finite for every array. **Not this claim.** |
| Recurrence known | partial | 96 | 100 of 196 are `not known`; not used |
| Adjuvant known | partial | 100 | 96 of 196 are `not known` |
| ICI / PD-1 treatment | no | 0 | resected 1995–2005 atlas |
| Tumor % / ABSOLUTE purity | no | 0 | only public proxy is an RNA epithelial mean-z |
| CLDN4 finite (`201428_at`) | yes | **196** | named Plus-2 probe; official GPL570 annot = CLDN4 / Entrez 1364 |
| CD8A finite (`205758_at`) | yes | **196** | named Plus-2 probe; official GPL570 annot = CD8A / Entrez 925 |
| CD274 finite (`223834_at`, `227458_at`) | yes | **196** | two Plus-2 probes; see probe note |
| **Primary pairwise n (CLDN4 + CD8A)** | yes | **196** | all NSCLC |
| **Primary pairwise n (CLDN4 + CD274)** | yes | **196** | all NSCLC |

Do not write n=196 as if it were one histology. The computable public n is **196 arrays = 196 patients**. The histology split is **106 / 66 / 24**.

## One-row table

| dataset | histology | platform | n arrays | n patients | CLDN4 | CD8A | CD274 | CLDN4–CD8A ρ (p) | adj ρ (p) | CLDN4–CD274 ρ (p) | adj ρ (p) | verdict (gene-level) | ICI |
|---|---|---|---:|---:|---|---|---|---|---|---|---|---|---|
| GSE37745 Botling | NSCLC mixed (106 adeno / 66 squam / 24 large) | GPL570 U133 Plus 2.0 RMA | **196** | **196** | `201428_at` | `205758_at` | two probes | **−0.185 (0.0094)** | −0.069 (0.34) | **−0.164 (0.021)** collapse | −0.046 (0.53) | **NO_EVIDENCE** after epithelium | none |

Full numbers: `tables/spearman_cldn4_vs_cd8a_cd274.tsv`, `tables/label_inventory.tsv`, `tables/spearman_named_probe_sensitivity.tsv`.

## CLDN4 vs CD8A / CD274 (primary, all NSCLC)

RMA as deposited (Bioconductor `affy` standard RMA in the series matrix). Spearman is rank-based. Bootstrap 95% CI, 2,000 resamples, seed `20260817`. Partial Spearman residualises on the pan-epithelial mean-z (EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7; **6/6 present**). Gene-level partners use max-mean collapse of official GEO GPL570 symbols. HOLDS rule (same as PR 229 / GSE4573): n≥40, ρ_adj&lt;0, p_adj&lt;0.05.

| pair | n | ρ | 95% CI | p | ρ_adj | p_adj | verdict |
|---|---:|---:|---|---:|---:|---:|---|
| CLDN4 vs **CD8A** | **196** | **−0.185** | −0.325 to −0.037 | 0.0094 | −0.069 | 0.34 | **NO_EVIDENCE** |
| CLDN4 vs **CD274** (collapse `227458_at`) | **196** | **−0.164** | −0.299 to −0.003 | 0.021 | −0.046 | 0.53 | **NO_EVIDENCE** |
| CLDN4 Q4 vs Q1 on CD8A | 49 vs 49 | — | — | MWU 0.028 | — | — | crude only (rank-biserial −0.26) |
| CLDN4 Q4 vs Q1 on CD274 | 49 vs 49 | — | — | MWU 0.036 | — | — | crude only (rank-biserial −0.25) |

Crude NSCLC-wide negatives exist. After eating the epithelial score they do **not** hold. Do not call a mixed-NSCLC CLDN4–CD8A or CLDN4–CD274 law from the gene-level residual.

## Histology split (do not pool)

| subset | n | pair | ρ | p | ρ_adj | p_adj | verdict |
|---|---:|---|---:|---:|---:|---:|---|
| LUAD | **106** | CLDN4 vs CD8A | −0.119 | 0.23 | +0.064 | 0.51 | **NO_EVIDENCE** |
| LUAD | **106** | CLDN4 vs CD274 (collapse) | **+0.227** | 0.019 | **+0.296** | 0.0022 | **OPPOSITE** |
| LUSC | **66** | CLDN4 vs CD8A | **−0.333** | 0.0063 | **−0.292** | 0.018 | **HOLDS** |
| LUSC | **66** | CLDN4 vs CD274 (collapse) | **−0.384** | 0.0014 | **−0.350** | 0.0043 | **HOLDS** |
| LCC | **24** | CLDN4 vs CD8A | −0.203 | 0.34 | −0.212 | 0.33 | **UNDERPOWERED** |
| LCC | **24** | CLDN4 vs CD274 (collapse) | −0.020 | 0.93 | −0.007 | 0.98 | **UNDERPOWERED** |

The mixed-NSCLC crude negative is histology-confounded. **LUSC n=66** is the only stratum that meets HOLDS for both partners. **LUAD n=106** is null for CD8A and **positive** for collapse-CD274. **LCC n=24** is below the n≥40 rule; do not write a large-cell claim.

## CD274 probe note (honest; do not pick one quietly)

GPL570 has two CD274 probes. Max-mean collapse keeps `227458_at` (mean RMA 7.96). The named / commonly cited probe is `223834_at` (mean RMA 6.48). They disagree on the mixed pool.

| subset | n | CD274 probe | ρ | p | ρ_adj | p_adj | verdict |
|---|---:|---|---:|---:|---:|---:|---|
| NSCLC_all | 196 | named `223834_at` | −0.285 | 5.2×10⁻⁵ | **−0.181** | **0.012** | **HOLDS** |
| NSCLC_all | 196 | collapse `227458_at` | −0.164 | 0.021 | −0.046 | 0.53 | **NO_EVIDENCE** |
| LUAD | 106 | named `223834_at` | +0.020 | 0.84 | +0.126 | 0.20 | NO_EVIDENCE |
| LUAD | 106 | collapse `227458_at` | +0.227 | 0.019 | +0.296 | 0.0022 | OPPOSITE |
| LUSC | 66 | named `223834_at` | −0.433 | 2.8×10⁻⁴ | −0.411 | 6.6×10⁻⁴ | HOLDS |
| LUSC | 66 | collapse `227458_at` | −0.384 | 0.0014 | −0.350 | 0.0043 | HOLDS |

CD8A has one clean probe (`205758_at`); named = collapse. CLDN4 named `201428_at` is the collapse probe (the other CLDN4 probe `1569421_at` is low-mean). Headline gene-level CD274 uses collapse. The named-probe residual HOLDS on n=196; that is a **probe disagreement**, not a second independent cohort.

## What CLDN4 *does* track here

| pair | n | ρ | 95% CI | p | ρ_adj | p_adj |
|---|---:|---:|---|---:|---:|---:|
| CLDN4 vs epithelial mean-z | 196 | **+0.481** | 0.362 to 0.591 | 9.8×10⁻¹³ | — | — |
| CLDN4 vs TACSTD2 | 196 | **+0.269** | 0.130 to 0.393 | 1.4×10⁻⁴ | +0.229 | 0.0013 |

CLDN4 sits on the epithelial / TACSTD2 side of this array. That is compatible with a tight-junction / tumour-cell program. It is the covariate that eats the mixed-NSCLC CD8A / collapse-CD274 residuals.

## TACSTD2 extras (companion only; not this claim)

| pair | n | ρ | p | ρ_adj | p_adj | verdict |
|---|---:|---:|---:|---:|---:|---|
| TACSTD2 vs CD8A | 196 | −0.085 | 0.23 | −0.049 | 0.50 | NO_EVIDENCE |
| TACSTD2 vs CD274 (collapse) | 196 | +0.070 | 0.33 | +0.114 | 0.11 | NO_EVIDENCE |

Do not treat CLDN4 as interchangeable with TACSTD2 on this matrix.

## What this does not test

- ICI response, PFS, or PD-1 treatment (1995–2005 surgical cohort; no ICI labels).
- Pathologist CD8 or PD-L1 IHC.
- OS as a CLDN4 claim (OS **is** deposited for 196/196; it is not scored here).
- A single NSCLC-wide CLDN4–immune law. Gene-level residual is null in the mixed pool; LUSC n=66 is the only HOLDS stratum; LUAD CD274 collapse is opposite.
- LCC n=24.

## Reproduce

```bash
python3 -m pip install -r methods/gse37745_cldn4/requirements.txt
python3 methods/gse37745_cldn4/analyze.py
```

Downloads (not committed) go to `$GSE37745_CLDN4_DATA` (default `/tmp/gse37745_cldn4`):

- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE37nnn/GSE37745/matrix/GSE37745_series_matrix.txt.gz`
- `https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz`

## Files

- `analyze.py` — download, probe map, Spearman / partial Spearman, figures
- `tables/label_inventory.tsv` — public vs missing fields; honest n
- `tables/probe_confirm.tsv` — named CLDN4 / CD8A / CD274 / TACSTD2 probes
- `tables/probe_means.tsv` — Plus-2 probe means; CD274 collapse ≠ named
- `tables/spearman_cldn4_vs_cd8a_cd274.tsv` — primary pairs by histology
- `tables/spearman_named_probe_sensitivity.tsv` — named vs collapse, by histology
- `tables/spearman_tacstd2_companion.tsv` — TACSTD2 vs CD8A / CD274
- `tables/highlow_cldn4_cd8a_cd274.tsv` — CLDN4 Q4 vs Q1
- `tables/sample_annotation.tsv` — 196 arrays
- `tables/summary.json`
- `figures/fig1_cldn4_vs_cd8a_cd274.png`
- `figures/fig2_cldn4_forest_by_histology.png`
