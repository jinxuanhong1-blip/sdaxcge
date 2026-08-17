# GSE3141 NSCLC array — CLDN4 vs CD8A / CD274

**Additive only.** Public Bild / Nevins resected NSCLC Affymetrix U133 Plus 2.0 series (Bild et al., *Nature* 2006, PMID 16273092; GEO [GSE3141](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3141)). Unit is the **array**. No slide was re-scored. No ICI arm.

Primary question: does **CLDN4** track **CD8A** or **CD274** (PD-L1) on this mixed NSCLC matrix? TACSTD2 is a same-run companion only.

## Honest n

| Item | Public? | n | Note |
|---|---|---:|---|
| Arrays in the series matrix | yes | **111** | 54,675 probes × 111 GSM; **0** missing MAS5 values |
| Unique GSM / unique titles | yes | **111** | Duke case IDs; all unique |
| Source | yes | **111** | every row is `frozen tissue of primary lung tumor`; **no** adjacent-normal rows |
| Histology LUAD (`Cell type: A`) | yes | **58** | GEO characteristic, not inferred |
| Histology LUSC (`Cell type: S`) | yes | **53** | GEO characteristic, not inferred |
| Histology other | yes | **0** | no large-cell / NOS / normal codes |
| OS months numeric | yes | **110** | GSM70230 is `Surv(months): >72` (alive LUSC); not a number |
| OS event (`STATUS` 0/1) | yes | **111** | deposited on every array |
| Stage | no | 0 | not a GEO characteristic |
| ICI / treatment | no | 0 | 2006 oncogenic-pathway atlas, not an ICI series |
| Tumor % / ABSOLUTE purity | no | 0 | only public proxy used here is an RNA epithelial mean-z |
| CLDN4 finite (max-mean) | yes | **111** | collapse picks `201428_at` (mean 2138) over `1569421_at` (mean 46) |
| CD8A finite | yes | **111** | named U133 probe `205758_at` |
| CD274 finite (max-mean) | yes | **111** | collapse picks `227458_at` (mean 1010) over `223834_at` (mean 369) |
| **Primary pairwise n (CLDN4 + CD8A)** | yes | **111** | this is the n used below |
| **Primary pairwise n (CLDN4 + CD274)** | yes | **111** | this is the n used below |

Do not write n=110 for the correlations. The computable public n is **111 arrays**. The `>72` token is an OS-time encoding only.

## One-row table

| dataset | histology | platform | n arrays | n LUAD / LUSC | CLDN4 | CD8A | CD274 | CLDN4–CD8A ρ (p) | adj ρ (p) | CLDN4–CD274 ρ (p) | adj ρ (p) | verdict | ICI |
|---|---|---|---:|---|---|---|---|---|---|---|---|---|---|
| GSE3141 Bild | NSCLC mixed | GPL570 U133 Plus 2.0 MAS5 | **111** | 58 / 53 | `201428_at` | `205758_at` | `227458_at` | **−0.181 (0.058)** | −0.010 (0.92) | **−0.078 (0.41)** | +0.048 (0.62) | **NO_EVIDENCE** | not deposited |

Full numbers: `tables/one_row.tsv`, `tables/spearman_cldn4_vs_cd8a_cd274.tsv`, `tables/label_inventory.tsv`.

## CLDN4 vs CD8A / CD274 (primary)

MAS5 as deposited (linear signal). Spearman is rank-based, so a log2 transform would not change ρ. Bootstrap 95% CI, 2,000 resamples, seed `20260817`. Partial Spearman residualises on the pan-epithelial mean-z (EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7; **6/6** present). HOLDS rule (same as PR 229 / GSE4573 CLDN4): n≥40, ρ_adj<0, p_adj<0.05.

| pair | subset | n | ρ | 95% CI | p | ρ_adj | p_adj | verdict |
|---|---|---:|---:|---|---:|---:|---:|---|
| CLDN4 vs **CD8A** | all tumors | **111** | **−0.181** | −0.370 to +0.017 | 0.058 | −0.010 | 0.92 | **NO_EVIDENCE** |
| CLDN4 vs **CD274** | all tumors | **111** | **−0.078** | −0.262 to +0.122 | 0.41 | +0.048 | 0.62 | **NO_EVIDENCE** |
| CLDN4 vs CD8A | LUAD | **58** | −0.151 | −0.390 to +0.096 | 0.26 | −0.090 | 0.51 | NO_EVIDENCE |
| CLDN4 vs CD274 | LUAD | **58** | −0.199 | −0.450 to +0.084 | 0.13 | −0.057 | 0.67 | NO_EVIDENCE |
| CLDN4 vs CD8A | LUSC | **53** | −0.247 | −0.528 to +0.044 | 0.075 | +0.007 | 0.96 | NO_EVIDENCE |
| CLDN4 vs CD274 | LUSC | **53** | +0.002 | −0.279 to +0.281 | 0.99 | +0.103 | 0.47 | NO_EVIDENCE |
| CLDN4 Q4 vs Q1 on CD8A | all tumors | 28 vs 28 | — | — | MWU 0.081 | — | — | null (rank-biserial −0.27) |
| CLDN4 Q4 vs Q1 on CD274 | all tumors | 28 vs 28 | — | — | MWU 0.70 | — | — | null (rank-biserial −0.06) |

The all-tumor CLDN4–CD8A crude interval includes both a modest negative and a near-zero effect (p=0.058). After the epithelial residual the association is gone. CLDN4–CD274 is null on the collapsed probe in every pre-specified subset.

## What CLDN4 *does* track here

| pair | n | ρ | 95% CI | p |
|---|---:|---:|---|---:|
| CLDN4 vs epithelial mean-z | 111 | **+0.621** | 0.479 to 0.729 | 3.5×10⁻¹³ |
| CLDN4 vs TACSTD2 | 111 | +0.199 | 0.004 to 0.379 | 0.036 |

CLDN4 sits on the epithelial side of this array (ρ = +0.62). The crude TACSTD2 correlation does not survive the same epithelial residual (ρ_adj = −0.022, p=0.82). That is compatible with a tight-junction / tumour-cell program. It is not evidence that CLDN4-high NSCLC is CD8-low or PD-L1-low on GSE3141.

## Probe sensitivity (not the claim)

Max-mean collapse is the primary map. CD274 has two clean GPL570 probes; the lower-mean probe is **not** interchangeable with the collapse.

| pair | n | ρ | p | ρ_adj | p_adj |
|---|---:|---:|---:|---:|---:|
| `201428_at` vs CD8A `205758_at` | 111 | −0.181 | 0.058 | −0.010 | 0.92 |
| `201428_at` vs CD274 `227458_at` (collapse) | 111 | −0.078 | 0.41 | +0.048 | 0.62 |
| `201428_at` vs CD274 `223834_at` (alt) | 111 | −0.241 | 0.011 | −0.138 | 0.15 |

The alt CD274 probe is crude-negative and still **NO_EVIDENCE** after the epithelial residual. Do not promote it over the max-mean row.

## TACSTD2 extras (companion; not this claim)

| pair | n | ρ | p | ρ_adj | p_adj | verdict |
|---|---:|---:|---:|---:|---:|---|
| TACSTD2 vs CD8A | 111 | −0.257 | 0.0065 | −0.178 | 0.063 | NO_EVIDENCE |
| TACSTD2 vs CD274 | 111 | −0.025 | 0.80 | +0.043 | 0.66 | NO_EVIDENCE |

TACSTD2–CD8A is the stronger crude inverse on this matrix; it still misses the locked HOLDS rule after epithelial correction. Do not treat CLDN4 as interchangeable with TACSTD2 here.

## What this does not test

- ICI response, PFS, or a treatment arm (not on GEO).
- Pathologist CD8 / PD-L1 IHC.
- Stage-adjusted models (stage is not deposited).
- OS as a CLDN4 claim (time is deposited, including one `>72` token; that is inventory only).
- A NSCLC-wide CLDN4–immune law. This is one public Plus-2 series (n=111 mixed).

## Reproduce

```bash
python3 -m pip install -r methods/gse3141_cldn4/requirements.txt
python3 methods/gse3141_cldn4/analyze.py
```

Downloads (not committed) go to `$GSE3141_CLDN4_DATA` (default `/tmp/gse3141_cldn4`):

- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE3nnn/GSE3141/matrix/GSE3141_series_matrix.txt.gz`
- `https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz`

## Files

- `analyze.py` — download, probe map, Spearman / partial Spearman, figures
- `tables/label_inventory.tsv` — public vs missing fields; honest n
- `tables/probe_confirm.tsv` — named CLDN4 / CD8A / CD274 / TACSTD2 probes
- `tables/gene_coverage.tsv` — collapse probe used
- `tables/spearman_cldn4_vs_cd8a_cd274.tsv` — primary, histology split, probe sensitivity
- `tables/highlow_cldn4.tsv` — CLDN4 Q4 vs Q1 on CD8A and CD274 (28 vs 28)
- `tables/one_row.tsv` — headline row
- `tables/sample_annotation.tsv` — 111 arrays
- `tables/summary.json`
- `figures/fig1_cldn4_vs_cd8a_cd274.png`
- `figures/fig2_cldn4_forest.png`
