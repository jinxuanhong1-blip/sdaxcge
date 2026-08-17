# GSE8894 NSCLC array — CLDN4 vs CD8A / CD274

**Additive only.** Public Lee / Kim resected NSCLC Affymetrix U133 Plus 2.0 series (Lee et al., *Clin Cancer Res* 2008, PMID 19010856; GEO [GSE8894](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE8894)). Unit is the **array**. No slide was re-scored. No ICI arm. This slice is **CLDN4 vs CD8A and CD274**. TACSTD2 is a same-run companion only.

## Honest n

| Item | Public? | n | Note |
|---|---|---:|---|
| Arrays in the series matrix | yes | **138** | 54,675 probes × 138 GSM; 0 missing GCRMA values |
| Unique GSM / unique `SML*` titles | yes | **138** | all unique |
| Tissue | yes | **138** | every array is `frozen tissue of primary lung tumor`; 0 adjacent/normal |
| Histology ADC (GEO `cell type`) | yes | **63** | `Adenocarcinoma` |
| Histology SCC (GEO `cell type`) | yes | **75** | `squamous cell carcinoma` |
| Paper / secondary 69 ADC + 69 SCC | text only | — | **not the GEO counts**. Do not write 69/69 for the tests |
| Recurrence event (status 1/0) | yes | 69 / 69 | deposited; **not** the primary claim |
| RFS months | yes | 138 | deposited; **not** the primary claim |
| Age | incomplete | 136 | missing on GSM225774 and GSM225864 |
| Gender | yes | 138 | Male 104, Female 34 |
| Stage | no | 0 | not a GEO characteristic |
| ICI / treatment | no | 0 | resected 2008 atlas, not an ICI series |
| Tumor % / ABSOLUTE purity | no | 0 | only public proxy used here is an RNA epithelial mean-z |
| CLDN4 finite (`201428_at`) | yes | **138** | max-mean collapse chose the named probe |
| CD8A finite (`205758_at`) | yes | **138** | named U133 Plus 2.0 probe |
| CD274 finite | yes | **138** | collapse chose `227458_at` (higher mean); named `223834_at` also present |
| **Primary pairwise n (CLDN4 + CD8A / CD274)** | yes | **138** | this is the n used below |

Do not write n=69 ADC / 69 SCC. The computable public n is **138 arrays (63 ADC + 75 SCC)**.

## One-row table

| dataset | histology | platform | n | CLDN4–CD8A ρ (p) | adj ρ (p) | verdict | CLDN4–CD274 ρ (p) | adj ρ (p) | verdict | ICI |
|---|---|---|---:|---|---|---|---|---|---|---|
| GSE8894 Lee | NSCLC mixed | GPL570 U133 Plus 2.0 GCRMA | **138** | **−0.318 (1.4×10⁻⁴)** | **−0.183 (0.032)** | **HOLDS** | −0.194 (0.022) | −0.100 (0.24) | **NO_EVIDENCE** | no |
| GSE8894 | ADC only | same | **63** | −0.241 (0.058) | −0.181 (0.16) | NO_EVIDENCE | −0.043 (0.74) | +0.026 (0.84) | NO_EVIDENCE | no |
| GSE8894 | SCC only | same | **75** | −0.161 (0.17) | −0.116 (0.33) | NO_EVIDENCE | +0.005 (0.97) | −0.026 (0.83) | NO_EVIDENCE | no |

Full numbers: `tables/one_row.tsv`, `tables/spearman_cldn4_vs_cd8a_cd274.tsv`, `tables/label_inventory.tsv`.

## CLDN4 vs CD8A / CD274 (primary)

GCRMA as deposited. Spearman is rank-based. Bootstrap 95% CI, 2,000 resamples, seed `20260817`. Partial Spearman residualises on the pan-epithelial mean-z (EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7; 6/6 present). HOLDS rule (same as PR 229 / GSE4573): n≥40, ρ_adj<0, p_adj<0.05.

| pair | stratum | n | ρ | 95% CI | p | ρ_adj | p_adj | verdict |
|---|---|---:|---:|---|---:|---:|---:|---|
| CLDN4 vs **CD8A** | all NSCLC | **138** | **−0.318** | −0.471 to −0.144 | 1.4×10⁻⁴ | **−0.183** | **0.032** | **HOLDS** |
| CLDN4 vs **CD274** | all NSCLC | **138** | −0.194 | −0.357 to −0.028 | 0.022 | −0.100 | 0.24 | **NO_EVIDENCE** |
| CLDN4 vs CD8A | ADC | 63 | −0.241 | −0.497 to +0.051 | 0.058 | −0.181 | 0.16 | NO_EVIDENCE |
| CLDN4 vs CD274 | ADC | 63 | −0.043 | −0.302 to +0.207 | 0.74 | +0.026 | 0.84 | NO_EVIDENCE |
| CLDN4 vs CD8A | SCC | 75 | −0.161 | −0.394 to +0.077 | 0.17 | −0.116 | 0.33 | NO_EVIDENCE |
| CLDN4 vs CD274 | SCC | 75 | +0.005 | −0.218 to +0.230 | 0.97 | −0.026 | 0.83 | NO_EVIDENCE |
| CLDN4 Q4 vs Q1 on CD8A | all | 35 vs 35 | — | — | MWU 3.7×10⁻⁴ | — | — | rank-biserial −0.50 |
| CLDN4 Q4 vs Q1 on CD274 | all | 35 vs 35 | — | — | MWU 0.053 | — | — | rank-biserial −0.27 |

Pooled NSCLC: higher CLDN4 tracks lower CD8A, and the inverse **survives** the epithelial residual (HOLDS). The same residual **removes** the CD274 association.

## Histology is a real mix, not 69/69

GEO `cell type` is **63 ADC + 75 SCC**. CLDN4 is higher in ADC; CD8A and CD274 are higher in SCC (MWU p = 1.7×10⁻¹⁰, 0.013, 3.9×10⁻⁴). That mix can manufacture a pooled inverse. Residualising ranks on a binary ADC indicator:

| pair | n | crude ρ (p) | ρ \| histology (p) | verdict |
|---|---:|---|---|---|
| CLDN4 vs CD8A | 138 | −0.318 (1.4×10⁻⁴) | **−0.247 (0.0037)** | still HOLDS |
| CLDN4 vs CD274 | 138 | −0.194 (0.022) | −0.036 (0.67) | NO_EVIDENCE |

The CD8A inverse is not only ADC-versus-SCC mixing. It is also **not** a within-histology law: ADC n=63 and SCC n=75 are both ≥40 and both **NO_EVIDENCE** after the epithelial residual. Do not promote the pooled HOLDS as an ADC or SCC result.

## What CLDN4 *does* track here

| pair | n | ρ | 95% CI | p |
|---|---:|---:|---|---:|
| CLDN4 vs epithelial mean-z | 138 | **+0.601** | 0.478 to 0.706 | 6.8×10⁻¹⁵ |
| CLDN4 vs TACSTD2 | 138 | +0.228 | 0.043 to 0.393 | 0.0072 |

CLDN4 sits on the epithelial / TACSTD2 side of this array. That is compatible with a tight-junction / tumour-cell program.

## TACSTD2 companion (not this claim)

| pair | n | ρ (p) | ρ_adj (p) |
|---|---:|---|---|
| TACSTD2 vs CD8A | 138 | −0.193 (0.023) | −0.177 (0.039) |
| TACSTD2 vs CD274 | 138 | −0.096 (0.26) | −0.082 (0.34) |

Same direction as CLDN4 vs CD8A on the pooled set. Not re-claimed as a TACSTD2 result.

## CD274 probe note

Max-mean collapse used `227458_at`. The named probe `223834_at` on the same 138 arrays is weaker: ρ = −0.144 (p = 0.092); ρ_adj = −0.135 (p = 0.12). Neither probe supports a CLDN4–CD274 HOLDS call.

## What this does not test

- ICI response (not an ICI series).
- Pathologist CD8 / PD-L1 IHC.
- Stage (not on GEO).
- A 69/69 ADC/SCC split (GEO is 63/75).
- Recurrence or RFS as a CLDN4 endpoint (labels are public; they are not this claim).
- A histology-specific CLDN4–CD8A law. The HOLDS call is **pooled NSCLC n=138** after an epithelial residual.

## Reproduce

```bash
python3 -m pip install -r methods/gse8894_cldn4/requirements.txt
python3 methods/gse8894_cldn4/analyze.py
```

Downloads (not committed) go to `$GSE8894_CLDN4_DATA` (default `/tmp/gse8894_cldn4`):

- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE8nnn/GSE8894/matrix/GSE8894_series_matrix.txt.gz`
- `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570&targ=self&form=text&view=data`

## Files

- `analyze.py` — download, probe map, Spearman / partial Spearman, figures
- `tables/label_inventory.tsv` — public vs missing fields; honest n
- `tables/probe_confirm.tsv` — named vs collapse probes for CLDN4 / CD8A / CD274 / TACSTD2
- `tables/gene_coverage.tsv`
- `tables/spearman_cldn4_vs_cd8a_cd274.tsv` — primary pairs, strata, histology residual, named-probe sensitivity
- `tables/highlow_cldn4.tsv` — CLDN4 Q4 vs Q1 on CD8A / CD274
- `tables/one_row.tsv`
- `tables/sample_annotation.tsv` — 138 arrays
- `tables/summary.json`
- `figures/fig1_cldn4_vs_cd8a.png`
- `figures/fig2_cldn4_vs_cd274.png`
- `figures/fig3_cldn4_forest.png`
