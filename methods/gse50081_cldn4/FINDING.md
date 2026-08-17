# GSE50081 LUAD array — CLDN4 vs CD8A / CD274 / ImmuneScore

**Additive only.** This slice is **CLDN4** on the public Der UHN181 early-stage NSCLC Affymetrix U133 Plus 2.0 series (Der et al., *J Thorac Oncol* 2014, PMID 24305008; GEO [GSE50081](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50081)). Unit is the **array**. No slide was re-scored. No ICI arm.

Primary histology is the GEO string `histology: adenocarcinoma`. Mixed NSCLC is not a LUAD finding.

## Honest n

| Item | Public? | n | Note |
|---|---|---:|---|
| Arrays in the series matrix | yes | **181** | 54,675 probes × 181 GSM; RMA log2 as deposited |
| Unique GSM / unique titles | yes | **181** | all unique |
| **histology = adenocarcinoma (LUAD primary)** | yes | **127** | strict GEO string; **this is the n used below** |
| histology = NSCLC-favor adenocarcinoma | yes | 1 | **not** folded into LUAD |
| histology = adenosquamous carcinoma | yes | 2 | **not** folded into LUAD |
| histology = squamous cell carcinoma (LUSC extra) | yes | 42 | strict GEO string; extra row only |
| histology = squamous cell carcinoma X2 | yes | 1 | **not** folded into LUSC extra |
| histology = large cell / mixed large cell | yes | 7 / 1 | not LUAD |
| Stage | yes | 127 LUAD | 1A 36, 1B 56, 2A 7, 2B 28; not a covariate |
| OS time / status | yes | 127 | deposited; not this claim |
| ICI / treatment | no | 0 | resected stage I–II atlas |
| Tumor % / ABSOLUTE purity | no | 0 | only public proxy is an RNA epithelial score |
| CLDN4 finite (`201428_at`) | yes | **127** | named Plus-2 probe; max-mean keeps this probe (alt `1569421_at` is low-mean) |
| CD8A finite (`205758_at`) | yes | **127** | single mapped probe |
| CD274 finite | yes | **127** | two probes; primary = max-mean `227458_at` |
| **Primary pairwise n (LUAD CLDN4 + partners)** | yes | **127** | this is the n used below |
| ESTIMATE ImmuneSignature genes present | yes | **137 / 141** | missing ARHGAP15, IFI30, KLRK1, CD302 |

Do not write n=181 for the LUAD tests. Do not write n=128 by adding “NSCLC-favor adenocarcinoma”.

## One-row table

| dataset | histology | platform | n arrays | CLDN4 | CD8A | CD274 | ImmuneScore | CLDN4–CD8A ρ (p) | CLDN4–CD274 ρ (p) | CLDN4–ImmuneScore ρ (p) | adj ImmuneScore | verdict | ICI |
|---|---|---|---:|---|---|---|---|---|---|---|---|---|---|
| GSE50081 Der UHN181 | LUAD (`adenocarcinoma`) | GPL570 U133 Plus 2.0 RMA log2 | **127** | `201428_at` | `205758_at` | `227458_at` (max-mean) | ESTIMATE ImmuneSignature mean-z, 137/141 | **−0.150 (0.093)** | **−0.076 (0.39)** | **−0.193 (0.030)** | −0.046 (0.61) | **NO_EVIDENCE** after epithelium | none |

Full numbers: `tables/spearman_cldn4_vs_partners.tsv`, `tables/label_inventory.tsv`.

## CLDN4 vs CD8A / CD274 / ImmuneScore (primary, LUAD n=127)

RMA log2 as deposited. Spearman is rank-based. Bootstrap 95% CI, 2,000 resamples, seed `20260817`. Partial Spearman residualises on the pan-epithelial mean-z (EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7; 6/6 present). HOLDS rule (same as PR 229 / PR 313): n≥40, ρ_adj<0, p_adj<0.05.

ImmuneScore here is **not** the official R `estimate::estimateScore` ssGSEA. It is the mean of gene-wise z-scores (z within the analysis set) of the Yoshihara 2013 ImmuneSignature (141 genes from `estimate` 1.0.13 `SI_geneset.gmt`).

| pair | n | ρ | 95% CI | p | ρ_adj | p_adj | verdict |
|---|---:|---:|---|---:|---:|---:|---|
| CLDN4 vs **CD8A** | **127** | **−0.150** | −0.327 to +0.032 | 0.093 | −0.039 | 0.66 | **NO_EVIDENCE** |
| CLDN4 vs **CD274** | **127** | **−0.076** | −0.247 to +0.105 | 0.39 | −0.030 | 0.74 | **NO_EVIDENCE** |
| CLDN4 vs **ImmuneScore** | **127** | **−0.193** | −0.365 to −0.007 | 0.030 | −0.046 | 0.61 | **NO_EVIDENCE** (adj) |
| CLDN4 Q4 vs Q1 on CD8A | 32 vs 32 | — | — | MWU 0.17 | — | — | null (rank-biserial −0.20) |
| CLDN4 Q4 vs Q1 on CD274 | 32 vs 32 | — | — | MWU 0.46 | — | — | null (rank-biserial −0.11) |
| CLDN4 Q4 vs Q1 on ImmuneScore | 32 vs 32 | — | — | MWU 0.028 | — | — | crude only (rank-biserial −0.32) |

CLDN4-high is not CD8A-low or CD274-low on this LUAD array. The crude ImmuneScore association is a small negative (CI just excludes 0) and **vanishes after the epithelial residual**. Under the locked HOLDS rule, nothing holds.

## CD274 probe note

GPL570 has two CD274 probes. Max-mean collapse keeps `227458_at` (higher mean). Named `223834_at` is the sensitivity. Both are null vs CLDN4.

| CD274 probe | role | n | mean RMA (LUAD) | ρ vs CLDN4 | p |
|---|---|---:|---:|---:|---:|
| `227458_at` | max-mean (primary) | 127 | 5.85 | −0.076 | 0.39 |
| `223834_at` | named | 127 | 4.66 | +0.018 | 0.84 |

Probes vs each other: ρ = 0.57. Do not treat a CD274-probe choice as a finding.

## Extra histology (not LUAD)

Do not quote the mixed row as a LUAD result. LUSC n=42 is the smallest extra that still meets n≥40; CIs are wide.

| cohort | n | CLDN4–CD8A ρ (p) | CLDN4–CD274 ρ (p) | CLDN4–ImmuneScore ρ (p) | ImmuneScore ρ_adj (p) |
|---|---:|---|---|---|---|
| LUAD adenocarcinoma | **127** | −0.150 (0.093) | −0.076 (0.39) | −0.193 (0.030) | −0.046 (0.61) |
| LUSC squamous extra | 42 | −0.240 (0.13) | +0.050 (0.75) | −0.232 (0.14) | −0.063 (0.70) |
| all NSCLC mixed | 181 | −0.118 (0.11) | −0.102 (0.17) | −0.104 (0.16) | +0.033 (0.66) |

Mixed ImmuneScore is diluted to null. That is composition, not a confirmation.

## What CLDN4 *does* track here

| pair | n | ρ | 95% CI | p |
|---|---:|---:|---|---:|
| CLDN4 vs epithelial mean-z | 127 | **+0.543** | 0.393 to 0.670 | 4.2×10⁻¹¹ |
| CLDN4 vs TACSTD2 | 127 | **+0.462** | 0.299 to 0.606 | 4.6×10⁻⁸ |

CLDN4 sits on the epithelial / TACSTD2 side of this array. That is compatible with a tight-junction / tumour-cell program. It is not evidence that CLDN4-high LUAD is immune-low after an epithelial residual.

## What this does not test

- ICI response (not an ICI series).
- Official ESTIMATE ssGSEA ImmuneScore from the R package (this is ImmuneSignature mean-z).
- Pathologist CD8 / PD-L1 IHC.
- A LUAD-wide CLDN4–immune law. This is one public U133 Plus 2.0 series.
- OS (labels are on GEO; not scored here).

## Reproduce

```bash
python3 -m pip install -r methods/gse50081_cldn4/requirements.txt
python3 methods/gse50081_cldn4/analyze.py
```

Downloads (not committed) go to `$GSE50081_CLDN4_DATA` (default `/tmp/gse50081_cldn4`):

- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50081/matrix/GSE50081_series_matrix.txt.gz`
- `https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz`
- `https://download.r-forge.r-project.org/src/contrib/estimate_1.0.13.tar.gz` (ImmuneSignature GMT only)

## Files

- `analyze.py` — download, probe map, Spearman / partial Spearman, figures
- `tables/label_inventory.tsv` — public vs missing fields; honest n
- `tables/histology_counts.tsv` — GEO histology strings
- `tables/probe_confirm.tsv` — named vs max-mean probes
- `tables/probe_all_mapped.tsv` — all mapped CLDN4 / CD8A / CD274 / TACSTD2 probes
- `tables/cd274_probe_sensitivity.tsv` — both CD274 probes vs CLDN4
- `tables/immunescore_gene_coverage.tsv` — 137/141 ImmuneSignature genes
- `tables/spearman_cldn4_vs_partners.tsv` — LUAD + LUSC extra + mixed
- `tables/highlow_cldn4.tsv` — CLDN4 Q4 vs Q1
- `tables/sample_annotation.tsv` — 181 arrays
- `tables/summary.json`
- `figures/fig1_cldn4_vs_partners.png`
- `figures/fig2_cldn4_forest.png`
