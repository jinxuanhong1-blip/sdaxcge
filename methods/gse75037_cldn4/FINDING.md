# GSE75037 LUAD array — CLDN4 vs CD8A / CD274 / ImmuneScore

**Additive CLDN4-only.** Public Girard / Gazdar / Lam matched-pair LUAD Illumina WG-6 v3 series (GEO [GSE75037](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75037); GPL6884; PMID 27354471). Unit is the **array**. No slide was re-scored. No ICI arm. **No dual-high** TACSTD2+CLDN4 gate.

Primary set is **tumor only**. The series is 83 LUAD + 83 matched adjacent non-malignant lung. Mixed tumor+normal is not a LUAD finding.

## Honest n

| Item | Public? | n | Note |
|---|---|---:|---|
| Arrays in the series matrix | yes | **166** | 48,803 beads × 166 GSM; MBCB + quantile + log2 as deposited |
| Unique GSM / unique titles | yes | **166** | titles like `05L4_N` / `05L4_T`; all unique |
| GEO design text (matched pairs) | text only | 83 | series title/design: 83 LUAD + 83 matched adjacent non-malignant |
| **histology Adenocarcinoma AND source Lung cancer (tumor primary)** | yes | **83** | strict GEO strings; **this is the n used below** |
| histology Non-malignant AND source Non-malignant lung | yes | 83 | matched adjacent; **not** mixed into the LUAD tests |
| Complete patient pairs (title prefix) | yes | 83 | 83 unique prefixes; title T/N matches histology/source (mismatch = 0) |
| Title suffix used as the tumor gate | no | 0 | suffix is consistent here but is not the gate |
| Stage (tumors) | yes | 83 | IA 24, IB 26, IIA 3, IIB 17, IIIA 9, III 1, IIIB 1, IV 2; not a covariate |
| EGFR / KRAS / LKB1 (tumors) | yes | 83 / 83 / 78 | EGFR MUT 20, KRAS MUT 35, LKB1 MUT 15 (5 tumors lack LKB1); not this claim |
| Race / smoker / pack-years | yes | 164 / 166 / 102 | 2 arrays lack a race token and shift later SOFT rows |
| ICI / treatment | **no** | 0 | resected matched-pair LUAD atlas |
| Tumor % / ABSOLUTE purity | **no** | 0 | only public proxy is an RNA epithelial score |
| CLDN4 finite (`ILMN_2132458`) | yes | **83** | single CLDN4 bead; official GPL6884 annot = CLDN4 / Entrez 1364 |
| CD8A finite (max-mean) | yes | **83** | three beads; collapse keeps `ILMN_2353732` |
| CD274 finite (`ILMN_1701914`) | yes | **83** | single CD274 bead; Entrez 29126 |
| **Primary pairwise n (tumor, CLDN4 + partners)** | yes | **83** | this is the n used below |
| ESTIMATE ImmuneSignature genes present | yes | **141 / 141** | all Yoshihara 2013 genes map after first-symbol collapse |
| Dual-high TACSTD2+CLDN4 gate | no | 0 | CLDN4-only slice |

Do not write n=166 for the LUAD tests. Do not write n=84 (some secondary papers inflate the pair count). The computable public tumor n is **83 arrays**.

## One-row table

| dataset | histology | platform | n arrays | n tumor | CLDN4 | CD8A | CD274 | ImmuneScore | CLDN4–CD8A ρ (p) | adj ρ (p) | CLDN4–CD274 ρ (p) | adj ρ (p) | CLDN4–ImmuneScore ρ (p) | adj ρ (p) | verdict | dual-high | ICI |
|---|---|---|---:|---:|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GSE75037 Girard / Gazdar / Lam | LUAD tumor only | GPL6884 WG-6 v3 MBCB log2 | 166 | **83** | `ILMN_2132458` | max-mean `ILMN_2353732` | `ILMN_1701914` | ESTIMATE ImmuneSignature mean-z, 141/141 | **−0.041 (0.71)** | −0.041 (0.71) | **+0.009 (0.93)** | +0.010 (0.93) | **−0.038 (0.73)** | −0.039 (0.73) | **NO_EVIDENCE** | no | none |

Full numbers: `tables/one_row.tsv`, `tables/spearman_cldn4_vs_partners.tsv`, `tables/label_inventory.tsv`.

## CLDN4 vs CD8A / CD274 / ImmuneScore (primary, tumor n=83)

Deposited Illumina intensity (MBCB background-correction, quantile-normalized, log2; Ding et al., *Nucleic Acids Res* 2008). Spearman is rank-based. Bootstrap 95% CI, 2,000 resamples, seed `20260817`. Partial Spearman residualises on the pan-epithelial mean-z (EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7; 6/6 present). HOLDS rule (same as PR 229 / GSE4573): n≥40, ρ_adj<0, p_adj<0.05.

ImmuneScore here is **not** the official R `estimate::estimateScore` ssGSEA. It is the mean of gene-wise z-scores (z within the analysis set) of the Yoshihara 2013 ImmuneSignature (141 genes from `estimate` 1.0.13 `SI_geneset.gmt`). All 141 genes are present after collapse.

| pair | n | ρ | 95% CI | p | ρ_adj | p_adj | verdict |
|---|---:|---:|---|---:|---:|---:|---|
| CLDN4 vs **CD8A** | **83** | **−0.041** | −0.247 to +0.178 | 0.71 | −0.041 | 0.71 | **NO_EVIDENCE** |
| CLDN4 vs **CD274** | **83** | **+0.009** | −0.195 to +0.215 | 0.93 | +0.010 | 0.93 | **NO_EVIDENCE** |
| CLDN4 vs **ImmuneScore** | **83** | **−0.038** | −0.254 to +0.177 | 0.73 | −0.039 | 0.73 | **NO_EVIDENCE** |
| CLDN4 Q4 vs Q1 on CD8A | 21 vs 21 | — | — | MWU 0.92 | — | — | null (rank-biserial +0.02) |
| CLDN4 Q4 vs Q1 on CD274 | 21 vs 21 | — | — | MWU 0.48 | — | — | null (rank-biserial +0.13) |
| CLDN4 Q4 vs Q1 on ImmuneScore | 21 vs 21 | — | — | MWU 1.00 | — | — | null (rank-biserial −0.00) |

CLDN4-high is not CD8A-low, CD274-low, or ImmuneScore-low on this LUAD tumor array. The epithelial residual does not change the signs or the verdict. Under the locked HOLDS rule, nothing holds.

## CD8A bead sensitivity

Three official CD8A beads. Max-mean keeps the highest-mean bead (`ILMN_2353732`). The lowest-mean bead is the most anti-correlated and is still null. None HOLDS after epithelial residual.

| CD8A probe | mean (tumor) | max-mean? | ρ vs CLDN4 (p) | ρ_adj (p) | verdict |
|---|---:|---|---|---|---|
| `ILMN_2353732` | 8.11 | **yes** | −0.041 (0.71) | −0.041 (0.71) | NO_EVIDENCE |
| `ILMN_1768482` | 7.84 | no | −0.052 (0.64) | −0.052 (0.64) | NO_EVIDENCE |
| `ILMN_1760374` | 4.93 | no | −0.182 (0.10) | −0.182 (0.10) | NO_EVIDENCE |

Do not treat a CD8A-bead choice as a finding.

## What CLDN4 *does* track here

| pair | n | ρ | 95% CI | p | ρ_adj (p) |
|---|---:|---:|---|---:|---|
| CLDN4 vs TACSTD2 | 83 | **+0.398** | 0.216 to 0.550 | 2.0×10⁻⁴ | +0.398 (2.1×10⁻⁴) |
| CLDN4 vs epithelial mean-z | 83 | −0.007 | −0.222 to +0.204 | 0.95 | — |

CLDN4 sits with TACSTD2 on this tumor array (companion only; not a dual-high gate). Within tumors, CLDN4 does **not** track the pan-epithelial mean-z. That residual is still reported because it is the locked covariate, not because epithelium explains the immune pairs (those pairs are already null).

TACSTD2 vs CD8A companion (not the claim): ρ = −0.037, p = 0.74; ρ_adj = −0.036, p = 0.75 (NO_EVIDENCE).

## Extra tissue (not LUAD tumor)

Do not quote the mixed row as a LUAD result. Adjacent lung is inventoried only. The mixed n=166 ImmuneScore association is tumor-versus-normal composition: CLDN4 is higher in tumor, ImmuneScore is higher in adjacent lung.

| cohort | n | CLDN4–CD8A ρ (p) | CLDN4–CD274 ρ (p) | CLDN4–ImmuneScore ρ (p) | ImmuneScore ρ_adj (p) |
|---|---:|---|---|---|---|
| LUAD tumor | **83** | −0.041 (0.71) | +0.009 (0.93) | −0.038 (0.73) | −0.039 (0.73) |
| adjacent non-malignant extra | 83 | −0.054 (0.63) | −0.104 (0.35) | +0.072 (0.52) | +0.105 (0.35) |
| all arrays mixed | 166 | −0.332 (1.3×10⁻⁵) | −0.211 (0.006) | −0.451 (1.0×10⁻⁹) | −0.187 (0.016) |

Mixed ImmuneScore looks like HOLDS. That is composition, not a confirmation. Do not write n=166 for the claim.

## What this does not test

- ICI response (not an ICI series).
- Official ESTIMATE ssGSEA ImmuneScore from the R package (this is ImmuneSignature mean-z).
- Pathologist CD8 / PD-L1 IHC.
- A TACSTD2+CLDN4 dual-high quadrant (this slice is CLDN4-only).
- A LUAD-wide CLDN4–immune law. This is one public Illumina series.
- n=166 mixed tumor+normal. Adjacent arrays are out of the primary n.

## Reproduce

```bash
python3 -m pip install -r methods/gse75037_cldn4/requirements.txt
python3 methods/gse75037_cldn4/analyze.py
```

Downloads (not committed) go to `$GSE75037_CLDN4_DATA` (default `/tmp/gse75037_cldn4`):

- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75037/matrix/GSE75037_series_matrix.txt.gz`
- `https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6884/annot/GPL6884.annot.gz`
- `https://download.r-forge.r-project.org/src/contrib/estimate_1.0.13.tar.gz` (ImmuneSignature GMT only)

## Files

- `analyze.py` — download, keyed histology/source parse, max-mean collapse, Spearman / partial Spearman, figures
- `tables/label_inventory.tsv` — public vs missing fields; honest n
- `tables/histology_counts.tsv` — Adenocarcinoma 83 / Non-malignant 83
- `tables/probe_confirm.tsv` — CLDN4 / CD8A / CD274 / TACSTD2 beads
- `tables/probe_all_mapped.tsv` — all mapped CLDN4 / CD8A / CD274 / TACSTD2 / epithelial beads
- `tables/cd8a_probe_sensitivity.tsv` — three CD8A beads vs CLDN4
- `tables/gene_coverage.tsv`
- `tables/immunescore_gene_coverage.tsv` — 141/141 ImmuneSignature genes
- `tables/spearman_cldn4_vs_partners.tsv` — tumor + adjacent extra + mixed
- `tables/primary_pairs.tsv` — the three requested tumor pairs
- `tables/one_row.tsv` — one-row summary
- `tables/highlow_cldn4.tsv` — CLDN4 Q4 vs Q1 (21 vs 21)
- `tables/sample_annotation.tsv` — 166 arrays
- `tables/summary.json`
- `figures/fig1_cldn4_vs_partners.png`
- `figures/fig2_cldn4_forest.png`
- `figures/fig3_epithelial_residual.png`
