# GSE101929 NSCLC array — CLDN4 vs CD8A / CD274 / ImmuneScore

**Additive CLDN4-only.** Public Mitchell / Ryan NCI-Maryland African-American vs European-American resected NSCLC Affymetrix U133 Plus 2.0 series (Mitchell et al., *Clin Cancer Res* 2017, PMID 29196495; GEO [GSE101929](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101929), SuperSeries [GSE102287](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102287)). Unit is the **array**. No slide was re-scored. No ICI arm. No dual-high. TACSTD2 is not scored.

Primary question: **CLDN4 vs CD8A**, **CLDN4 vs CD274**, and **CLDN4 vs ImmuneScore** on **tumor** arrays. Ancestry (AA / EA) and paper histology are extras, not the headline n.

## Honest n

Do not write n=66 for the pairwise tests. The series matrix has 66 arrays because it mixes 32 tumors with 34 adjacent normals. Do not write n=41 (patients) or n=22 / 19 (paper mRNA-cohort headcount) as the correlation n.

| Item | Public? | n | Note |
|---|---|---:|---|
| Arrays in the series matrix | yes | **66** | 54,675 probes × 66 GSM; **0** missing GCOS values |
| Unique GSM / unique titles | yes | **66** | `patient <ID>T` / `…N`; all unique |
| Unique patients (`individual`) | yes | **41** | paper mRNA cohort 22 AA + 19 EA; GEO matches |
| **Tumor arrays (`tumor_normal status: T`)** | yes | **32** | **this is the n used below** (16 AA + 16 EA) |
| Adjacent-normal arrays (`N`) | yes | **34** | inventoried; not mixed into the NSCLC pairwise tests |
| Matched T+N patients | yes | **25** | paper text is 11 AA + 14 EA matched pairs; unmatched T or N remain |
| AA tumor / EA tumor | yes | **16 / 16** | ancestry extra only; both **UNDERPOWERED** for HOLDS (n≥40) |
| Histology per array | **no** | 0 | not a GEO characteristic |
| Paper Table 1 histology (aggregate) | paper only | 32 | ADC 27 / SCC 1 / other 4; **not** mapped to GSM IDs |
| Stage | yes | 32 | I / II / III on every array; not a covariate here |
| OS days / months / years + lung-cancer death | yes | 66 | deposited; not this claim |
| ICI / treatment | no | 0 | resected atlas, not an ICI series |
| Numeric tumor % / ABSOLUTE purity | no | 0 | macro-dissected; only public proxy is RNA epithelial mean-z |
| CLDN4 finite on tumors (`201428_at`) | yes | **32** | max-mean keeps the named probe |
| CD8A finite (`205758_at`) | yes | **32** | only CD8A probe on GPL570 |
| CD274 finite (`227458_at`) | yes | **32** | max-mean of two clean probes |
| ESTIMATE ImmuneSignature genes present | yes | **137 / 141** | missing ARHGAP15, IFI30, KLRK1, CD302 (same four as GSE50081) |
| **Primary pairwise n (tumors)** | yes | **32** | complete-case CLDN4 + CD8A + CD274 + ImmuneScore |

The computable public n for this claim is **32 tumor arrays**. Use 16 / 16 only when the row is ancestry-stratified. Do not pool the 34 normals.

## One-row table

| dataset | histology | platform | n tumors | CLDN4 | CD8A | CD274 | ImmuneScore | CLDN4–CD8A ρ (p) | CLDN4–CD274 ρ (p) | CLDN4–ImmuneScore ρ (p) | adj ρ epi (p) | verdict | OS / ICI |
|---|---|---|---:|---|---|---|---|---|---|---|---|---|---|
| GSE101929 Mitchell NCI-Maryland | NSCLC mixed; per-array histology **not deposited** | GPL570 U133 Plus 2.0 GCOS default | **32** | `201428_at` | `205758_at` | `227458_at` | ESTIMATE ImmuneSignature mean-z, 137/141 | **−0.436 (0.013)** | **−0.207 (0.26)** | **−0.547 (0.0012)** | −0.277 (0.13) / −0.190 (0.31) / −0.296 (0.11) | **UNDERPOWERED** (n=32 < 40) | OS deposited; ICI not |

Full numbers: `tables/one_row.tsv`, `tables/spearman.tsv`, `tables/label_inventory.tsv`.

## CLDN4 vs CD8A / CD274 / ImmuneScore (primary, tumors n=32)

GCOS default as deposited. The paper’s Partek RMA matrix is not the GEO series table. Spearman is rank-based. Bootstrap 95% CI, 2,000 resamples, seed `20260817`. Partial Spearman residualises on the pan-epithelial mean-z (EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7; **6/6** present) or, as an extra, on a binary AA indicator. HOLDS rule (same as PR 229 / GSE4573): n≥40, ρ_adj<0, p_adj<0.05.

ImmuneScore here is **not** the official R `estimate::estimateScore` ssGSEA. It is the mean of gene-wise z-scores (z within the analysis set) of the Yoshihara 2013 ImmuneSignature (141 genes from `estimate` 1.0.13 `SI_geneset.gmt`).

| pair | subset | n | ρ | 95% CI | p | ρ_adj epi (p) | ρ_adj ancestry (p) | verdict |
|---|---|---:|---:|---|---:|---|---|---|
| CLDN4 vs **CD8A** | tumors | **32** | **−0.436** | −0.670 to −0.119 | 0.013 | −0.277 (0.13) | −0.445 (0.012) | **UNDERPOWERED** |
| CLDN4 vs **CD274** | tumors | **32** | **−0.207** | −0.523 to +0.140 | 0.26 | −0.190 (0.31) | −0.203 (0.27) | **UNDERPOWERED** |
| CLDN4 vs **ImmuneScore** | tumors | **32** | **−0.547** | −0.746 to −0.252 | 0.0012 | −0.296 (0.11) | −0.550 (0.0013) | **UNDERPOWERED** |

Crude CLDN4–CD8A and CLDN4–ImmuneScore are negative on these 32 tumors. After the epithelial residual both intervals include 0 (p=0.13 and 0.11). CLDN4–CD274 is null before and after. n=32 is below the locked HOLDS threshold. Do not promote the crude p-values.

Named-probe CD8A (`201428_at` vs `205758_at`) is identical to the collapsed row (single CD8A probe). Named-probe CD274 `223834_at` is still null (ρ=−0.117, p=0.53; epi adj −0.191, p=0.30).

### Q4 vs Q1 (descriptive; 8 vs 8)

| subset | endpoint | n Q4 vs Q1 | MWU p | rank-biserial |
|---|---|---|---:|---:|
| tumors | CD8A | 8 vs 8 | 0.021 | −0.69 |
| tumors | CD274 | 8 vs 8 | 0.28 | −0.34 |
| tumors | ImmuneScore | 8 vs 8 | 0.0047 | −0.81 |

Quartile n=8 vs 8 is not a HOLDS test.

## Extra: ancestry (not the headline n)

GEO `race` is AA or EA on every array. Tumor split is 16 / 16. Both rows are **UNDERPOWERED**. Residualising the pooled 32 on an AA indicator does not replace the epithelial residual.

| pair | subset | n | ρ (p) | ρ_adj epi (p) | verdict |
|---|---|---:|---|---|---|
| CLDN4 vs CD8A | AA tumors | 16 | −0.588 (0.017) | −0.482 (0.069) | UNDERPOWERED |
| CLDN4 vs CD274 | AA tumors | 16 | −0.247 (0.36) | −0.190 (0.50) | UNDERPOWERED |
| CLDN4 vs ImmuneScore | AA tumors | 16 | −0.459 (0.074) | −0.274 (0.32) | UNDERPOWERED |
| CLDN4 vs CD8A | EA tumors | 16 | −0.229 (0.39) | −0.008 (0.98) | UNDERPOWERED |
| CLDN4 vs CD274 | EA tumors | 16 | −0.144 (0.59) | −0.120 (0.67) | UNDERPOWERED |
| CLDN4 vs ImmuneScore | EA tumors | 16 | −0.597 (0.015) | −0.281 (0.31) | UNDERPOWERED |

CLDN4, CD8A, CD274, and ImmuneScore do not differ by ancestry on the 32 tumors (MWU p=0.38 / 0.92 / 0.84 / 0.90). That is not an immune claim.

## Extra: histology (not deposited per array)

GEO has no histology characteristic. Mitchell et al. Table 1 reports aggregate counts only for the 16 EA + 16 AA tumor patients used in the paper: adenocarcinoma 13+14=27, squamous 0+1=1, other 3+1=4. Those counts cannot be joined to GSM IDs from public files. No ADC / SCC correlation row is computed. Do not invent per-array histology.

## What CLDN4 *does* track here

| pair | n | ρ | 95% CI | p |
|---|---:|---:|---|---:|
| CLDN4 vs epithelial mean-z | 32 | **+0.610** | 0.235 to 0.832 | 2.1×10⁻⁴ |
| CLDN4 vs ancestry (AA) | 32 | −0.162 | — | 0.37 |

CLDN4 sits on the epithelial side of this array. The crude CD8A / ImmuneScore negatives shrink after that residual. That is compatible with a tight-junction / tumour-cell program. It is not a HOLDS immune finding on n=32.

## What this does not test

- Dual-high TACSTD2+CLDN4 (not scored).
- Official ESTIMATE ssGSEA ImmuneScore from the R package (this is ImmuneSignature mean-z).
- Per-array ADC vs SCC (labels are not on GEO).
- ICI response (not an ICI series).
- Pathologist CD8 / PD-L1 IHC.
- OS (labels are on GEO; not scored here).
- A NSCLC-wide CLDN4–immune law. This is one public U133 Plus 2.0 series with 32 tumors.

## Methods

- **Matrix:** GEO `GSE101929_series_matrix.txt.gz` (Affymetrix GCOS default as deposited). Probe × sample values used as published.
- **Annotation:** official GPL570 `Gene symbol` (NCBI platform table, first `///` token). **Max-mean** probe collapse to HUGO.
- **CD8** = `CD8A`. **CD274** = collapsed max-mean (`227458_at`); `223834_at` is a named-probe sensitivity row.
- **ImmuneScore:** unweighted mean of gene-wise z for the Yoshihara ImmuneSignature, z-scored within the analysis set (137/141 genes).
- **Epithelial residual:** unweighted mean of gene-wise z for EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7.
- **Ancestry residual:** binary AA indicator on the pooled 32 tumors only. Not applied inside AA or EA.
- **Partial Spearman:** Pearson of rank residuals; df = n − 2 − k.
- **Not done:** no official ESTIMATE ssGSEA. No dual-high. No ICI / OS model. No imputed histology.

## Files

- `tables/one_row.tsv` — headline row
- `tables/spearman.tsv` — all n / ρ / p / CI / residuals
- `tables/label_inventory.tsv` — honest n
- `tables/sample_annotation.tsv` — per-array genes + ancestry + tissue
- `tables/probe_confirm.tsv` / `gene_coverage.tsv` / `immunescore_gene_coverage.tsv` / `highlow_cldn4.tsv` / `ancestry_contrast.tsv` / `summary.json`
- `figures/fig1_cldn4_vs_cd8a_cd274_immunescore.png` / `fig2_forest.png`

## Reproduce

```bash
python3 -m pip install -r methods/gse101929_cldn4/requirements.txt
export GSE101929_CLDN4_DATA=/tmp/gse101929_cldn4
python3 methods/gse101929_cldn4/analyze.py
```

Downloads (not committed) go to `$GSE101929_CLDN4_DATA`:

- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE101nnn/GSE101929/matrix/GSE101929_series_matrix.txt.gz`
- `https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz`
- `https://download.r-forge.r-project.org/src/contrib/estimate_1.0.13.tar.gz` (ImmuneSignature GMT only)
