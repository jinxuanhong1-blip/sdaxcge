# GSE309652 NSCLC metabolic-subtype array — CLDN4 vs CD8A / ImmuneScore / CD274 / IFN

**Additive only. CLDN4-only. No dual-high.** Public Lee / National Cancer Center Hospital Korea stage-IV NSCLC NanoString series (GEO [GSE309652](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE309652); GPL31904 nCounter Human Metabolic Pathways Panel, 768 genes). Unit is the **array**. No slide was re-scored.

Primary question: does **CLDN4** track **CD8A**, an **ImmuneScore**, **CD274**, or an **IFN** score on this mixed NSCLC matrix? TACSTD2 is a same-run companion only and is also absent.

**CLDN4 is not on the panel.** The primary pairwise n is **0**. This is an honest ABSENT, not a null Spearman.

## Honest n

| Item | Public? | n | Note |
|---|---|---:|---|
| Arrays in the series matrix | yes | **72** | 768 genes × 72 GSM (GSM9271195–GSM9271266); **0** missing nSolver values |
| Unique GSM / unique titles | yes | **72** | pathology IDs; all unique |
| Source | yes | **72** | archival **tumor** (design text); `source_name` = histology; **no** adjacent-normal / blood arrays |
| Histology LUAD (`tissue: Adenocarcinoma`) | yes | **34** | GEO characteristic, not inferred |
| Histology LUSC (`tissue: Squamous cell carcinoma`) | yes | **27** | GEO characteristic, not inferred |
| Histology Other | yes | **11** | GEO code `Other`; not recoded |
| Anti-PD-(L)1 ICI cohort | yes | **72** | stage IV; pretreatment archival tumor |
| ICI response R / NR | yes | **72** | **24 R / 48 NR** |
| PD-L1 IHC numeric | yes | **71** | GEO `pdl1` percent; one non-numeric |
| SUVmax | yes | **69** | PET; not this claim |
| OS / PFS time | no | 0 | design text mentions survival; times are **not** GEO characteristics |
| Tumor % / ESTIMATE purity | no | 0 | not deposited |
| CLDN4 finite | **no** | **0** | 0/768 symbols; 0 `CLDN*` genes; GPL31904 SOFT platform table has no CLDN4 |
| CD8A finite | yes | **72** | symbol `CD8A` |
| CD274 finite | yes | **72** | symbol `CD274` |
| IFNG finite | yes | **72** | symbol `IFNG` |
| ImmuneScore A1 8-gene | partial | 72 | **7/8** (`CXCL10` missing). ESTIMATE is **not** computed |
| IFN Ayers-6 | partial | 72 | **4/6** (`CXCL10`, `HLA-DRA` missing). `HLA-DRB1` is on the panel and is **not** substituted |
| Epithelial mean-z (6-gene) | **no** | **0** | EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7 all absent. `KRT1` is on the panel and is **not** used as a substitute |
| TACSTD2 companion | **no** | **0** | not on GPL31904 |
| **Primary pairwise n (CLDN4 + any partner)** | **no** | **0** | this is the n used below |

Do not write n≈72 as if CLDN4 were measured. The computable public n for a CLDN4 test is **0 arrays**. The series itself is **72** tumor arrays.

No PMID is attached on the GEO record (public 30 Sep 2025).

## One-row table

| dataset | histology | platform | n arrays | n LUAD / LUSC / Other | CLDN4 | CD8A | ImmuneScore | CD274 | IFN | epi residual | CLDN4–CD8A | CLDN4–ImmuneScore | CLDN4–CD274 | CLDN4–IFN | verdict | ICI |
|---|---|---|---:|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GSE309652 Lee NCC Korea | NSCLC mixed | GPL31904 Metabolic Pathways 768-gene nSolver | **72** | **34 / 27 / 11** | **ABSENT** | present n=72 | A1 **7/8** (no CXCL10); ESTIMATE no | present n=72 | Ayers-6 **4/6**; IFNG n=72 | **0/6 NOT_COMPUTABLE** | n=0 | n=0 | n=0 | n=0 | **ABSENT** | stage IV anti-PD-(L)1; 24 R / 48 NR; CLDN4 vs response not testable |

Full numeric row: `tables/one_row.tsv`. Coverage: `tables/gene_coverage.tsv`, `tables/label_inventory.tsv`.

## CLDN4 vs CD8A / ImmuneScore / CD274 / IFN (primary)

nSolver counts as deposited (linear). Spearman would be rank-based. HOLDS rule (same as PR 229 / GSE4573 CLDN4): n≥40, residual Spearman ρ<0, p<0.05. Epithelial residual = mean-z of EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7.

| pair | subset | n | ρ | p | ρ_adj epithelial | verdict |
|---|---|---:|---|---|---|---|
| CLDN4 vs **CD8A** | all tumors | **0** | — | — | — | **ABSENT** |
| CLDN4 vs **ImmuneScore** (A1 7/8) | all tumors | **0** | — | — | — | **ABSENT** |
| CLDN4 vs **CD274** | all tumors | **0** | — | — | — | **ABSENT** |
| CLDN4 vs **IFN** (Ayers 4/6) | all tumors | **0** | — | — | — | **ABSENT** |
| CLDN4 vs CD8A | LUAD | **0** | — | — | — | **ABSENT** |
| CLDN4 vs CD8A | LUSC | **0** | — | — | — | **ABSENT** |
| CLDN4 Q4 vs Q1 | all tumors | 0 vs 0 | — | — | — | **ABSENT** |

The LUAD/LUSC split is labeled and is reported. It does not create a CLDN4 test. Other n=11 is not a third histology claim.

Epithelial residual is **not computable** (0/6 locked genes). A residualised ρ is not invented.

## What *is* on the panel (inventory; not this claim)

CD8A, CD274, and IFNG are finite on all 72 arrays. The A1 T-cell effector list is 7/8 (`CD8A GZMA GZMB IFNG EOMES CXCL9 TBX21`; `CXCL10` absent). Ayers-6 is 4/6 (`IFNG STAT1 CXCL9 IDO1`; `CXCL10` and `HLA-DRA` absent). These scores are written so the missing members are visible. They are **not** a CLDN4 result.

| pair | n | ρ | 95% CI | p | role |
|---|---:|---:|---|---:|---|
| CD8A vs CD274 | 72 | +0.420 | 0.217 to 0.591 | 2.4×10⁻⁴ | panel sanity |
| CD8A vs IFNG | 72 | +0.362 | 0.142 to 0.562 | 0.0018 | panel sanity |
| CD8A vs ImmuneScore (7/8) | 72 | +0.765 | 0.620 to 0.861 | 5.1×10⁻¹⁵ | CD8A is inside the score |
| CD274 vs PD-L1 IHC | 71 | +0.072 | −0.164 to +0.287 | 0.55 | RNA vs IHC; not CLDN4 |

LUAD vs LUSC on the deposited immune genes is null (CD8A MWU p=0.97; CD274 p=0.86; IFNG p=0.13). That is inventory, not a CLDN4–histology law.

No `CLDN*` symbol is on GPL31904. TACSTD2 is absent and is not substituted.

## What this does not test

- A CLDN4–CD8A / ImmuneScore / CD274 / IFN association (CLDN4 was never measured).
- Dual-high TACSTD2/CLDN4 (both genes absent; dual-high was not used).
- Pathologist CLDN4 IHC.
- ESTIMATE ImmuneScore or a 6-gene epithelial residual (genes not on the metabolic panel).
- CLDN4 vs ICI response / PFS / OS (response is deposited; CLDN4 is not).
- A NSCLC-wide CLDN4–immune law.

## Reproduce

```bash
python3 -m pip install -r methods/gse309652_cldn4/requirements.txt
python3 methods/gse309652_cldn4/analyze.py
```

Downloads (not committed) go to `$GSE309652_CLDN4_DATA` (default `/tmp/gse309652_cldn4`):

- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE309nnn/GSE309652/matrix/GSE309652_series_matrix.txt.gz`
- `https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL31nnn/GPL31904/soft/GPL31904_family.soft.gz`

## Files

- `analyze.py` — download, symbol map, honest-n inventory, figures
- `tables/label_inventory.tsv` — public vs missing fields; honest n
- `tables/probe_confirm.tsv` — requested / locked genes; CLDN4 n_finite=0
- `tables/gene_coverage.tsv` — primary / A1 / Ayers-6 / epithelial / claudin family
- `tables/spearman_cldn4_vs_cd8a_immunescore_cd274_ifn.tsv` — ABSENT rows, including LUAD/LUSC
- `tables/spearman_panel_inventory.tsv` — CD8A/CD274/IFNG sanity (not the claim)
- `tables/histology_contrast.tsv` — LUAD vs LUSC on deposited immune genes
- `tables/highlow_cldn4.tsv` — empty Q4/Q1 (no CLDN4)
- `tables/one_row.tsv` — headline row
- `tables/sample_annotation.tsv` — 72 arrays
- `tables/summary.json`
- `figures/fig1_gene_coverage.png`
- `figures/fig2_immune_by_histology.png`
- `figures/fig3_one_row.png`
