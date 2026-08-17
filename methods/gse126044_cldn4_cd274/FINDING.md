# GSE126044 — CLDN4 vs CD274 / HLA (CD274 extra)

**Additive. CLDN4 only.** No TACSTD2 gate. Patient is the unit.

Public Cho et al. anti-PD-1 NSCLC biopsy RNA-seq ([GSE126044](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126044); Cho et al., *Nat Commun* 2020). Author count matrix → `log2(CPM+1)`.

**Already known (not re-audited as a claim).** CLDN4 vs ESTIMATE ImmuneScore on this same matrix is in leftover OPEN ICI bulk ([PR #296](https://github.com/jinxuanhong1-blip/sdaxcge/pull/296)): n=16, ρ=−0.57, p=0.022. This page recomputes ImmuneScore only so CD274 / HLA can be residualized on the same covariate. The ImmuneScore row is a **given**, not a new finding.

**Extra.** `CD274` (PD-L1) was not on the leftover immune-axis table (that table had CD8A / IFN / 6-gene MHC / ImmuneScore). HLA-A/B/C and MHC-I mean-z here are **HLA-A/B/C only (3/3)**, not B2M/TAP.

This page does **not** re-audit R vs NR (CLDN4 Cliff δ=−0.53, p=0.115 is given in PR #296 / PR #149). Correlations use all **16** patients.

Reproduce: `python3 methods/gse126044_cldn4_cd274/analyze.py`

---

## Honest n

| Item | Public? | n | Note |
|---|---|---:|---|
| GEO series GSM | yes | **16** | GSE126044 series matrix |
| Count-matrix columns | yes | **16** | `GSE126044_counts.txt.gz` |
| Unique GSM / unique titles | yes | **16 / 16** | 1:1 |
| Merged patients (unit) | yes | **16** | title `RNA-seq_*` stripped to matrix ID |
| GEO responder / non-responder | yes | **5 / 11** | used only as a plot color; not the test |
| Fresh / FFPE (`sample:` characteristic) | yes | **11 / 5** | all 5 FFPE are NR (given) |
| LUAD vs LUSC | **no** | 0 | GEO says NSCLC; no histology table |
| Tumor % / ABSOLUTE purity | **no** | 0 | RNA ImmuneScore is the public proxy |
| CLDN4 / CD274 / HLA-A / HLA-B / HLA-C finite | yes | **16 / 16 / 16 / 16 / 16** | symbols on the count matrix |
| ImmuneScore finite (Yoshihara immune list mean-z) | yes | **16** | coverage 140/141 |
| **Primary pairwise n (CLDN4 + CD274)** | yes | **16** | this is the n used below |

Do not write n>16. Do not drop FFPE to chase a p-value. Do not treat 5 vs 11 as the CD274 n.

---

## One-row table

| Cohort | n | scale | CD274 ρ (p) **extra** | HLA-A ρ (p) | HLA-B ρ (p) | HLA-C ρ (p) | MHC-I mean-z ρ (p) | CD274 \| ImmuneScore ρ (p) | MHC-I \| ImmuneScore ρ (p) | ImmuneScore ρ (p) **given** |
|---|---:|---|---|---|---|---|---|---|---|---|
| GSE126044 | **16** | log2(CPM+1) | **-0.635 (0.00818)** | -0.244 (0.362) | -0.488 (0.055) | -0.474 (0.0639) | -0.459 (0.0738) | -0.532 (0.0411) | +0.159 (0.572) | -0.568 (0.0218) |

Numeric row: `tables/one_row.tsv`. Axis-level tests: `tables/spearman.tsv`.

---

## How to read the row

**CD274 extra (n=16).** CLDN4 is **CD274-low** (ρ=-0.635, p=0.00818). After residualizing ranks on the already-known ImmuneScore, the partial is -0.532 (0.0411). CD8A residual (sensitivity, not primary): -0.580 (0.0234).

**HLA / MHC-I (3/3).** HLA-B/C trend negative; HLA-A is weaker. MHC-I mean-z ρ=-0.459, p=0.0738. After ImmuneScore the MHC-I residual is +0.159 (0.572). Leftover PR #296 used a **6-gene** MHC set (HLA-A/B/C+B2M+TAP1+TAP2, raw ρ=−0.44, p=0.087; partial | ImmuneScore +0.19, p=0.49). That is a different score and is not re-used as the MHC-I number here.

**ImmuneScore given.** Recomputed Yoshihara immune-list mean-z matches the leftover row (ρ=-0.568, p=0.0218, n=16, coverage 140/141). Cited so the CD274 residual has a named covariate. Not a new ImmuneScore claim.

Honest n=16 is small. A large |ρ| is required to reach p<0.05. Do not pool this ρ with GSE218989 or other ICI bulk CD274 rows (signs conflict; see PR #326).

---

## Methods

- **Matrix:** GEO `GSE126044_counts.txt.gz`. Columns matched to series titles after stripping the `RNA-seq_` prefix (`Dis_01` …).
- **Scale:** `log2(CPM+1)` from library-size CPM. Symbols collapsed by mean after uppercasing.
- **CD274 extra** = HUGO `CD274` (aliases PDCD1LG1 / PDL1 / B7-H1 checked; matrix uses `CD274`).
- **MHC-I mean-z** = gene-wise z of `HLA-A`, `HLA-B`, `HLA-C` only (HLA-A,HLA-B,HLA-C; 3/3). Not B2M/TAP.
- **ImmuneScore** = Yoshihara 2013 immune-signature **mean z** on the log2-CPM matrix (`estimate_gene_sets.csv`; same list as PR #296). Coverage 140/141. This is not the Affymetrix-calibrated ESTIMATE R purity transform.
- **Partial Spearman:** Pearson of rank residuals; df = n − 3. Primary residual covariate = **ImmuneScore** (already-known infiltrate axis). CD8A residual is sensitivity.
- **Unit:** patient = GSM. All 16 patients enter every pairwise test. No FFPE drop. No LUAD/LUSC split (none deposited).

---

## What was not done

- No TACSTD2 column on the claim table.
- No R vs NR re-test and no search for a cutoff that hits p=0.019.
- No Fisher-z meta-analysis with other ICI bulk CD274 rows.
- No invented histology or purity.
- ImmuneScore is **given**, not sold as a new finding.

---

## Files

- `FINDING.md` — this page
- `analyze.py` — download-once, score, write tables/figures
- `estimate_gene_sets.csv` — Yoshihara stromal/immune lists (tidyestimate / PR #296)
- `tables/one_row.tsv` — the one-row table
- `tables/spearman.tsv` / `inventory.tsv` / `summary.json`
- `processed/GSE126044_patient.tsv`
- `figures/fig1_cldn4_vs_cd274_hla.png`, `fig2_forest.png`
