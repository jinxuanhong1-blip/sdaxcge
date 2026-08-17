# GSE7670 paired LUAD array — CLDN4 vs CD8A / CD274 / ImmuneScore

**Additive only. CLDN4-only. No dual-high.** Public Su / Taipei VGH adjacent-normal–tumor matched lung array (Su et al., *BMC Genomics* 2007, PMID 17540040; GEO [GSE7670](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7670)). Platform is Affymetrix HG-U133A (GPL96). Unit is the **array**. No slide was re-scored. No ICI arm.

Primary pairs are **tumor CLDN4 vs CD8A**, **tumor CLDN4 vs ImmuneScore**, and **tumor CLDN4 vs CD274**. Matched adjacent-normal arrays are a **companion only** and are not mixed into those tests. TACSTD2 is a same-run companion only. Dual-high (CLDN4-high and TACSTD2-high) was **not** done.

## Honest n

| Item | Public? | n | Note |
|---|---|---:|---|
| Arrays in the series matrix | yes | **66** | 22,283 probes × 66 GSM; **0** missing MAS5. Do **not** use 66 for the LUAD tests |
| Unique GSM / unique titles | yes | **66** | all unique |
| Patients (GEO `overall_design` text) | text only | **27** | “pairwise samples from 27 patients” at Taipei VGH |
| Patient tumor arrays | yes | **27** | 26 LUAD + **1 LCC** (pair 19; title typo `lage cell`) |
| **LUAD patient tumors (primary n)** | yes | **26** | source `*T` and adenocarcinoma. Excludes pair 19, mixtures, cell lines |
| LUAD matched adjacent-normal | yes | **26** | companion only; not in the tumor pairwise tests |
| LCC pair 19 (tumor + normal) | yes | 2 | inventoried; **not LUAD**; not in primary n |
| Tissue mixtures (Taichung VGH) | yes | 2 | `Adj. N. mixture` + `T. mixture`; not an individual tumor |
| Commercial normal lung | yes | 2 | Stratagene 735020 + Clontech 636524 |
| Cell lines | yes | 8 | A549, H661, H1299, NL-20, CL1-0, CL1-1, CL1-5, CL1-5-F4 |
| Stage | no | 0 | design text says early and late stages; not a per-array characteristic |
| ICI / treatment | no | 0 | 2007 surgical atlas, not an ICI series |
| Tumor % / ABSOLUTE / ESTIMATE | no | 0 | ESTIMATE skipped; RNA epithelial mean-z is the purity proxy |
| CLDN4 finite (`201428_at`) | yes | **26** | named GPL96 / Entrez 1364 |
| CD8A finite (`205758_at`) | yes | **26** | named GPL96 / Entrez 925 |
| **CD274 finite (Entrez 29126)** | **no** | **0** | no GPL96 probe; Plus-2 `223834_at` / `227458_at` are not on U133A |
| ImmuneScore A1 genes present | yes | **7/8** | EOMES **absent** on U133A |
| Epithelial genes present | yes | **6/6** | EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7 |
| **Primary pairwise n (CLDN4 + CD8A / ImmuneScore)** | yes | **26** | this is the n used below |
| **Primary pairwise n (CLDN4 + CD274)** | **no** | **0** | ABSENT — do not invent a PD-L1 ρ |

Do not write n=66, n=27, or n=54. The computable public LUAD tumor n is **26 arrays**. The GEO 27-patient sentence includes one large-cell pair.

PDCD1LG2 (PD-L2, `220049_s_at`) **is** on U133A. It is not CD274 and is **not** substituted.

## One-row table

| dataset | histology | platform | n LUAD tumors | n pairs | CLDN4 | CD8A | CD274 | ImmuneScore | CLDN4–CD8A ρ (p) | adj ρ (p) | CLDN4–ImmuneScore ρ (p) | adj ρ (p) | CLDN4–CD274 | verdict | ICI |
|---|---|---|---:|---:|---|---|---|---|---|---|---|---|---|---|---|
| GSE7670 Su | LUAD paired (pair 19 LCC out) | GPL96 U133A MAS5 | **26** | 26 | `201428_at` | `205758_at` | **ABSENT** | A1 7/8 mean-z (no EOMES; not ESTIMATE) | **−0.303 (0.13)** | +0.170 (0.42) | **−0.227 (0.26)** | −0.243 (0.24) | **ABSENT (n=0)** | **UNDERPOWERED** | no |

Full numbers: `tables/one_row.tsv`, `tables/spearman_cldn4_vs_genes.tsv`, `tables/label_inventory.tsv`.

## CLDN4 vs CD8A / ImmuneScore / CD274 (primary)

MAS5 as deposited (target intensity 500). Spearman is rank-based, so a log2 transform would not change ρ. Bootstrap 95% CI, 2,000 resamples, seed `20260817`. Partial Spearman residualises on the pan-epithelial mean-z (EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7; **6/6** present). ImmuneScore is the A1 8-gene T-cell effector mean-z (`CD8A GZMA GZMB IFNG EOMES CXCL9 CXCL10 TBX21`) on genes present (**7/8**; EOMES absent). This is **not** Yoshihara ESTIMATE ImmuneScore.

HOLDS rule (same as PR 229 / GSE4573): n≥40, ρ_adj<0, p_adj<0.05. **n=26 cannot HOLDS.**

| pair | n | ρ | 95% CI | p | ρ_adj | p_adj | verdict |
|---|---:|---:|---|---:|---:|---:|---|
| CLDN4 vs **CD8A** | **26** | **−0.303** | −0.686 to +0.178 | 0.13 | +0.170 | 0.42 | **UNDERPOWERED** |
| CLDN4 vs **ImmuneScore** | **26** | **−0.227** | −0.623 to +0.239 | 0.26 | −0.243 | 0.24 | **UNDERPOWERED** |
| CLDN4 vs **CD274** | **0** | — | — | — | — | — | **ABSENT** |
| CLDN4 vs epithelial mean-z | 26 | **+0.826** | 0.597 to 0.938 | 2.0×10⁻⁷ | — | — | tracks epithelium |
| CLDN4 Q4 vs Q1 on CD8A | 7 vs 7 | — | — | MWU 0.16 | — | — | rank-biserial −0.47 |
| CLDN4 Q4 vs Q1 on ImmuneScore | 7 vs 7 | — | — | MWU 0.46 | — | — | rank-biserial −0.27 |
| CLDN4 Q4 vs Q1 on CD274 | — | — | — | — | — | — | **ABSENT** |

The crude CLDN4–CD8A interval includes both a moderate negative and a positive effect. After the epithelial residual the sign **flips** to a non-significant positive. CLDN4–ImmuneScore stays weakly negative and non-significant. Neither pair can HOLDS at n=26.

Including the LCC pair (n=27 patient tumors) does not change the call: CD8A ρ=−0.224 (p=0.26), ImmuneScore ρ=−0.168 (p=0.40); both still UNDERPOWERED.

## What CLDN4 *does* track here

| pair | n | ρ | 95% CI | p |
|---|---:|---:|---|---:|
| CLDN4 vs epithelial mean-z | 26 | **+0.826** | 0.597 to 0.938 | 2.0×10⁻⁷ |
| CLDN4 vs TACSTD2 | 26 | +0.447 | 0.038 to 0.754 | 0.022 |

CLDN4 sits on the epithelial / TACSTD2 side of this array. The TACSTD2 correlation does not survive the epithelial residual (ρ_adj=+0.083, p=0.69). That is compatible with a tight-junction / tumour-cell program. It is not evidence that CLDN4-high LUAD is CD8-low or ImmuneScore-low on GSE7670.

CD8A vs ImmuneScore ρ=+0.589 (p=0.0015) is a positive control and is partly by construction (CD8A is one of the seven genes used).

## Companion: tumor vs matched normal (not the claim)

26 LUAD pairs. Wilcoxon signed-rank, two-sided. These arrays are **not** entered into the primary Spearman.

| endpoint | n pairs | T > N | median T | median N | p | role |
|---|---:|---:|---:|---:|---:|---|
| CLDN4 | 26 | 23 | 1197 | 542 | 1.6×10⁻⁶ | companion |
| CD8A | 26 | 9 | 408 | 455 | 0.28 | companion |
| ImmuneScore (z on 52 T+N) | 26 | 7 | −0.30 | +0.02 | 0.022 | companion |
| CD274 | 0 | — | — | — | — | **ABSENT** |

CLDN4 is higher in the tumor half of the pair. That is a tumour-versus-normal fact. It is **not** a CLDN4-versus-immune fact and is not the primary claim.

## TACSTD2 companion (not this claim; not dual-high)

| pair | n | ρ (p) | ρ_adj (p) |
|---|---:|---|---|
| TACSTD2 vs CD8A | 26 | −0.461 (0.018) | −0.300 (0.15) |
| CLDN4 vs TACSTD2 | 26 | +0.447 (0.022) | +0.083 (0.69) |

No CLDN4-high / TACSTD2-high quadrant was cut. Dual-high is out of scope.

## What this does not test

- ICI response, PFS, or a treatment arm (not on GEO).
- Pathologist CD8 / PD-L1 IHC.
- A PD-L1 RNA claim (CD274 is not on the chip).
- A PD-L2 (PDCD1LG2) substitute for CD274 (ρ=−0.200, p=0.33; not used).
- ESTIMATE ImmuneScore / xCell / MCP-counter.
- Stage-adjusted models (stage is not deposited per array).
- n=66, n=27, or a 27-pair LUAD law (pair 19 is LCC).
- A dual-high CLDN4+TACSTD2 class.
- A NSCLC-wide CLDN4–immune law. This is one public U133A paired series with honest LUAD tumor n=26.

## Reproduce

```bash
python3 -m pip install -r methods/gse7670_cldn4/requirements.txt
python3 methods/gse7670_cldn4/analyze.py
```

Downloads (not committed) go to `$GSE7670_CLDN4_DATA` (default `/tmp/gse7670_cldn4`):

- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE7nnn/GSE7670/matrix/GSE7670_series_matrix.txt.gz`
- `https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz`

## Files

- `analyze.py` — download, sample class, probe map, Spearman / partial Spearman, paired Wilcoxon, figures
- `tables/label_inventory.tsv` — public vs missing fields; honest n
- `tables/probe_confirm.tsv` — named CLDN4 / CD8A / TACSTD2 / PDCD1LG2; CD274 Plus-2 probes marked absent
- `tables/gene_coverage.tsv`
- `tables/spearman_cldn4_vs_genes.tsv` — primary, LCC-in sensitivity, companions
- `tables/highlow_cldn4.tsv` — CLDN4 Q4 vs Q1 (7 vs 7)
- `tables/paired_tumor_vs_normal.tsv` — companion 26 LUAD pairs
- `tables/one_row.tsv` — headline row
- `tables/sample_annotation.tsv` — 66 arrays
- `tables/luad_pairs.tsv` — 26 LUAD pair IDs
- `tables/summary.json`
- `figures/fig1_cldn4_vs_cd8a_immunescore.png`
- `figures/fig2_cldn4_forest.png`
- `figures/fig3_companion_paired_cldn4.png`
- `figures/fig4_epithelial_residual.png`
