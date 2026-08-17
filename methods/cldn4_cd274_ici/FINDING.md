# CLDN4 vs CD274 and HLA-A/B/C — four lung ICI bulk series

**Additive. CLDN4 only.** No TACSTD2 gate and no dual-high score. Patient is the unit.

Public GEO processed matrices: [GSE218989](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218989), [GSE126044](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126044), [GSE166449](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166449), [GSE182328](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182328). This page does **not** re-audit ICI response / DCB / Akkermansia contrasts already reported in PR #292 and PR #296.

**Question.** Does CLDN4 track PD-L1 (`CD274`) or classical MHC-I (`HLA-A`, `HLA-B`, `HLA-C`) on these ICI-treated tumor matrices? Partial Spearman residualizes ranks on `CD8A` when it is present (it is, on all four).

**Scale.** TPM matrices: `log2(TPM+1)`. Count matrices: `log2(CPM+1)`. MHC-I mean-z = gene-wise z of HLA-A/B/C only (3/3; not B2M/TAP). Symbols are on every matrix. No histology split: GSE218989 has none on GEO; the other three are deposited as NSCLC / LUAD without a usable LUAD vs LUSC table for these tests.

Reproduce: `python3 methods/cldn4_cd274_ici/analyze.py` (GEO files cached under `/tmp/cldn4_cd274_ici/`).

---

## One row per cohort (CLDN4)

| Cohort | n | scale | CD274 ρ (p) | HLA-A ρ (p) | HLA-B ρ (p) | HLA-C ρ (p) | MHC-I mean-z ρ (p) | CD274 \| CD8A ρ (p) | MHC-I \| CD8A ρ (p) |
|---|---:|---|---|---|---|---|---|---|---|
| **GSE218989** | **355** | log2(TPM+1) | **+0.148 (0.0053)** | +0.094 (0.076) | −0.026 (0.62) | +0.007 (0.90) | +0.021 (0.70) | **+0.219 (3.3×10⁻⁵)** | +0.046 (0.39) |
| **GSE126044** | **16** | log2(CPM+1) | **−0.635 (0.0082)** | −0.244 (0.36) | −0.488 (0.055) | −0.474 (0.064) | −0.459 (0.074) | **−0.580 (0.023)** | −0.180 (0.52) |
| **GSE166449** | **22** | log2(TPM+1) | −0.020 (0.93) | −0.001 (1.00) | −0.057 (0.80) | +0.049 (0.83) | +0.001 (1.00) | +0.175 (0.45) | +0.243 (0.29) |
| **GSE182328** | **44** | log2(CPM+1) | −0.169 (0.27) | −0.040 (0.80) | −0.121 (0.43) | +0.124 (0.42) | +0.005 (0.97) | −0.032 (0.84) | +0.148 (0.34) |

Full numeric row: `tables/one_row.tsv`. Axis-level tests: `tables/spearman.tsv`.

Do **not** pool these four ρ values. Signs conflict and n is 16 / 22 / 44 / 355.

---

## How to read the row

**GSE218989 (n=355).** CLDN4 is weakly **CD274-high** (ρ=+0.148, p=0.0053). After CD8A the residual is larger, not smaller (partial ρ=+0.219, p=3.3×10⁻⁵). Classical MHC-I is **null** (mean-z ρ=+0.021, p=0.70; HLA-B/C ~0; HLA-A +0.094, p=0.076). This matches the already-published GSE218989 MHC-I row in PR #292 (ρ=+0.021, p=0.70) and does not make CLDN4 an MHC-I-high marker. CD274-high after CD8A residual is compatible with a combination-candidate axis; it is not an ICI-response test (GEO R vs NR for CLDN4 remains the given null, p=0.40).

**GSE126044 (n=16).** CLDN4 is **CD274-low** (ρ=−0.635, p=0.0082). The CD8A partial stays negative (ρ=−0.580, p=0.023). HLA-B/C trend negative (p=0.055 / 0.064); MHC-I mean-z ρ=−0.459, p=0.074; after CD8A the MHC-I residual is null (p=0.52). Honest n=16. This is the opposite CD274 sign from GSE218989. Do not average them.

**GSE166449 (n=22).** All six tests are null (|ρ|≤0.06 raw; partials p≥0.29). Pembrolizumab LUAD TPM. GEO titles are 7 R / 15 NR; those labels are not used here.

**GSE182328 (n=44).** CD274 ρ=−0.169, p=0.27; MHC-I mean-z ρ=+0.005, p=0.97. CD8A partials remain null. GEO has Akkermansia detectable vs not only — not RECIST / DCB / PFS / MPR. The leftover 6-gene MHC set in PR #296 (HLA-A/B/C+B2M+TAP1+TAP2, ρ=−0.21, p=0.16) is a different score and is not re-used.

---

## What was not done

- No TACSTD2 column on the claim table.
- No Fisher-z meta-analysis and no “CLDN4 vs PD-L1 holds in ICI bulk” sentence.
- No invented histology, purity, or response labels.
- Partial is **CD8A**, not ESTIMATE ImmuneScore (that residual is already in PR #296).

---

## Files

- `FINDING.md` — this page
- `analyze.py` — download-once, score, write tables/figures
- `tables/one_row.tsv` — the four-row table
- `tables/spearman.tsv` / `inventory.tsv` / `summary.json`
- `processed/*_patient.tsv`
- `figures/fig1_forest_raw.png`, `fig2_partial_cd8a.png`, `fig3_scatters.png`
