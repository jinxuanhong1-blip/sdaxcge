# GSE41271 LUAD array — CLDN4 vs CD8A / CD274

**Additive only.** Public MD Anderson / Girard NSCLC Illumina WG-6 v3 series (Sato / Girard / Wistuba / Minna; GEO [GSE41271](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41271); GPL6884). Unit is the **array**. No slide was re-scored. No ICI arm.

Primary slice is **LUAD**. The series is mixed NSCLC. Do not write n=275 for the tests.

## Honest n

| Item | Public? | n | Note |
|---|---|---:|---|
| Arrays in the series matrix | yes | **275** | 48,803 beads × 275 GSM; 0 missing values on the three genes |
| Unique GSM / unique titles | yes | **275** | titles like `15-T`; all unique |
| Histology Adenocarcinoma (LUAD) | yes | **183** | `histology: Adenocarcinoma` |
| Histology Squamous | yes | 80 | dropped from primary |
| Other histology | yes | 12 | LCC-NE 3, adenosquamous 2, sarcomatoid 2, sarcomatoid-adeno 1, sarcomatoid-squamous 1, pleomorphic 1, NSCLC 1, SCLC & NSCLC 1 |
| GEO design text ADC | text only | 183 | series summary: “mainly adenocarcinomas (n = 183)” — matches the label count |
| CLDN4 finite (`ILMN_2132458`) | yes | **183** | single CLDN4 bead; official GPL6884 annot = CLDN4 / Entrez 1364 |
| CD8A finite (max-mean) | yes | **183** | three beads; collapse keeps `ILMN_2353732` |
| CD274 finite (`ILMN_1701914`) | yes | **183** | single CD274 bead; Entrez 29126 |
| **Primary pairwise n (LUAD, CLDN4 + CD8A + CD274)** | yes | **183** | this is the n used below |
| Tobacco history | yes | 272 | 3 arrays lack the token (2 squamous, 1 LUAD); later SOFT rows shift if parsed by row index |
| Stage / vital A–D / recurrence Y–N | yes | 275 | keyed parse; **not** used in the correlations |
| ICI / PD-1 treatment | **no** | 0 | surgical 1997–2005 tumors |

Do not write n=275 for the correlations. Do not add adenosquamous or sarcomatoid-adenocarcinoma to LUAD. The computable public LUAD n is **183 arrays**.

## One-row table

| dataset | histology | platform | n arrays | n LUAD | CLDN4 | CD8A | CD274 | CLDN4–CD8A ρ (p) | adj ρ (p) | CLDN4–CD274 ρ (p) | adj ρ (p) | verdict | ICI |
|---|---|---|---:|---:|---|---|---|---|---|---|---|---|---|
| GSE41271 Girard / MDACC | LUAD | GPL6884 WG-6 v3 | 275 | **183** | `ILMN_2132458` | max-mean `ILMN_2353732` | `ILMN_1701914` | **−0.142 (0.054)** | −0.083 (0.27) | **−0.097 (0.19)** | −0.075 (0.31) | **NO_EVIDENCE** | not an ICI series |

Full numbers: `tables/spearman_cldn4_vs_cd8a_cd274.tsv`, `tables/label_inventory.tsv`.

## CLDN4 vs CD8A / CD274 (primary)

Deposited Illumina processed intensity (already on a log-like scale). Spearman is rank-based. Bootstrap 95% CI, 2,000 resamples, seed `20260817`. Partial Spearman residualises on the pan-epithelial mean-z (EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7; 6/6 present). HOLDS rule (same as PR 229 / GSE4573): n≥40, ρ_adj<0, p_adj<0.05.

| pair | n | ρ | 95% CI | p | ρ_adj | p_adj | verdict |
|---|---:|---:|---|---:|---:|---:|---|
| CLDN4 vs **CD8A** | **183** | **−0.142** | −0.276 to −0.001 | 0.054 | −0.083 | 0.27 | **NO_EVIDENCE** |
| CLDN4 vs **CD274** | **183** | **−0.097** | −0.239 to +0.056 | 0.19 | −0.075 | 0.31 | **NO_EVIDENCE** |
| CLDN4 Q4 vs Q1 on CD8A | 46 vs 47 | — | — | MWU 0.058 | — | — | null (rank-biserial −0.23) |
| CLDN4 Q4 vs Q1 on CD274 | 46 vs 47 | — | — | MWU 0.29 | — | — | null (rank-biserial −0.13) |

Unadjusted CLDN4–CD8A is a weak negative (p=0.054). The bootstrap CI just excludes zero. After the epithelial residual it is not significant. That is **not** HOLDS. CLDN4–CD274 is null on both the crude and residual tests.

## CD8A bead sensitivity

Three official CD8A beads. Max-mean keeps the highest-mean bead (`ILMN_2353732`). That bead is also the most anti-correlated with CLDN4. The other two are weaker. None HOLDS after epithelial residual.

| CD8A probe | mean (LUAD) | max-mean? | ρ vs CLDN4 (p) | ρ_adj (p) | verdict |
|---|---:|---|---|---|---|
| `ILMN_2353732` | 8.11 | **yes** | −0.142 (0.054) | −0.083 (0.27) | NO_EVIDENCE |
| `ILMN_1768482` | 7.85 | no | −0.120 (0.11) | −0.055 (0.46) | NO_EVIDENCE |
| `ILMN_1760374` | 4.60 | no | −0.074 (0.32) | −0.011 (0.89) | NO_EVIDENCE |

## What CLDN4 *does* track here

| pair | n | ρ | 95% CI | p | ρ_adj (p) |
|---|---:|---:|---|---:|---|
| CLDN4 vs epithelial mean-z | 183 | **+0.373** | 0.247 to 0.495 | 2.0×10⁻⁷ | — |
| CLDN4 vs TACSTD2 | 183 | **+0.283** | 0.142 to 0.419 | 1.0×10⁻⁴ | +0.135 (0.069) |

CLDN4 sits on the epithelial / TACSTD2 side of this LUAD array. That is compatible with a tight-junction / tumour-cell program. It is not evidence that CLDN4-high LUAD is CD8A-low or CD274-low after eating epithelium.

TACSTD2 vs CD8A companion (not the claim): ρ = −0.150, p = 0.043; ρ_adj = −0.077, p = 0.30 (NO_EVIDENCE after epithelium).

## What this does not test

- ICI response, PFS, or a locked OS model (dates are on GEO; this slice does not score survival).
- Pathologist CD8 / PD-L1 IHC.
- ESTIMATE TumorPurity (not computed here; the residual is an RNA epithelial mean-z).
- A LUAD-wide CLDN4–immune law. This is one public Illumina series.
- n=275 mixed NSCLC. Squamous and the 12 other histologies are out of the primary n.

## Reproduce

```bash
python3 -m pip install -r methods/gse41271_cldn4/requirements.txt
python3 methods/gse41271_cldn4/analyze.py
```

Downloads (not committed) go to `$GSE41271_CLDN4_DATA` (default `/tmp/gse41271_cldn4`):

- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE41nnn/GSE41271/matrix/GSE41271_series_matrix.txt.gz`
- `https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6884/annot/GPL6884.annot.gz`

## Files

- `analyze.py` — download, keyed histology parse, max-mean collapse, Spearman / partial Spearman, figures
- `tables/label_inventory.tsv` — public vs missing fields; honest n
- `tables/histology_counts.tsv` — 10 GEO histology labels
- `tables/probe_confirm.tsv` — CLDN4 / CD8A / CD274 / TACSTD2 beads
- `tables/cd8a_probe_sensitivity.tsv` — three CD8A beads vs CLDN4
- `tables/gene_coverage.tsv`
- `tables/spearman_cldn4_vs_cd8a_cd274.tsv`
- `tables/one_row.tsv`
- `tables/highlow_cldn4.tsv` — CLDN4 Q4 vs Q1 (46 vs 47)
- `tables/sample_annotation.tsv` — 183 LUAD arrays
- `tables/summary.json`
- `figures/fig1_cldn4_vs_cd8a.png`
- `figures/fig2_cldn4_vs_cd274.png`
