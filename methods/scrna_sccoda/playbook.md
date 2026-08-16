# scCODA / Dirichlet-multinomial composition playbook
# scCODA / 狄利克雷-多项 组成分析手册

**Focus / 主题.** Patient-level cell-type composition ~ **malignant TACSTD2**
(and, separately, ~ **MPR**) in three public lung-tumor scRNA series.
**Combinatorial question:** which **cohort × annotation** recovers **T/NK or TLS
down** in TACSTD2-high.

> **Additive methods + public-only numbers.** Everything lives under
> `methods/scrna_sccoda/`. No private / EGA / GSA raw. No fabricated p-values.
> The replicate unit is the **patient/sample**, never the cell.

---

## 0. Why not Spearman-on-fractions as the primary test

Cell-type fractions sum to 1. A rise in epithelium *mechanically* lowers T/NK
even if T/NK counts are unchanged (Aitchison; Gloor et al.). scCODA
(Büttner, Ostner et al., *Nat Commun* 2021) models the count vector with a
**Dirichlet-multinomial** and an ALR link relative to a reference cell type.

This folder implements that estimand as a **frequentist DM-GLM** (default).
`sccoda` + tensorflow is an optional Bayesian companion; if it does not
import, the run writes `not_available` and does not invent spike-and-slab
inclusion probabilities.

Naive MWU / Spearman on fractions are written to `companion_fraction_naive`
so a reader can see how often the compositional and non-compositional tests
disagree.

---

## 1. Estimands (pre-registered)

Two **separate** models. Do not put TACSTD2 and MPR in one formula on these
n.

| Estimand | Exposure | Primary compartments | Family |
|---|---|---|---|
| TACSTD2-high (primary) | median split of **malignant / author-Epi** TACSTD2 mean log1p, within cohort, among samples with ≥10 such cells | TNK, TLS | `primary_tacstd2` |
| TACSTD2 continuous | within-cohort z-score of the same score | all types | `secondary_tacstd2_continuous` |
| MPR | MPR (incl. pCR) vs NMPR; pre-treatment / unlabeled samples dropped | all types | `mpr` |

**Orientation.** A **negative** ALR slope on TNK or TLS means that compartment
is **down** in TACSTD2-high (or in MPR, for the MPR family).

**Recovery rule (primary only).** Patient-level **permutation p** of the ALR
slope (exact enumeration when `C(n, n_high) ≤ 2000`, else 499 shuffles).
This p **does not grow with n_cells**. DM-GLM LRT is written as
`diagnostic_dm_lrt` only — large libraries make it anti-conservative.
Recovery = TACSTD2-high, compartment ∈ {TNK, TLS} (or T/NK/B/plasma),
**effect < 0** and **BH q < 0.10** inside `primary_tacstd2`.
Bonferroni on that grid is a sensitivity column.

---

## 2. Cohorts (public only)

| Cohort | Paper | n patients | MPR public? | Author labels? | Matrix used |
|---|---|---|---|---|---|
| **GSE207422** | Hu et al. *Genome Med* 2023 | 15 (12 post) | yes | **no** (GEO sample sheet only) | 175 MB UMI txt.gz |
| **GSE241934 IIT** | Zhang/Sun/Zhong *Cell Rep Med* 2024 (NEOTIDE) | 11 EGFR-mut | yes | yes (`major.cell.type`, `cell.type`) | 444 MB MTX |
| **GSE241934 RWC** | same, real-world EGFR-WT | 34 | yes | yes | 1.2 GB MTX |
| **GSE253013** | Sze/Xiang *Cancer Res* 2024 | 9 tumor | **no** | Garnett `cell_type` (limited table) | 9.3 GB RDS **not** re-downloaded |

GSE241934 design is FACS 7AAD− CD235a− (non-erythrocyte), not CD45 sort, so
epithelium is present but **immune-rich**. Composition is still compositional
inside that capture.

GSE253013 has **no public MPR / ICI label**. The MPR family is skipped there.
The 9.3 GB RDS exceeds the 2 GB processed-file budget; patient-level tables
are public-derived (see `data/public_derived/PROVENANCE.md`).

Skipped (catalog only): GSE207422 GSA `HRA001033`, GSE241934 GSA bulk
`HRA007419`.

---

## 3. Annotation axis (the other half of the grid)

| Scheme | Cohorts | What it is |
|---|---|---|
| `hu_markers` | GSE207422, GSE241934 | argmax of Hu canonical log1p marker scores |
| `hu_collapsed` | same | T+NK→TNK, B+plasma→TLS |
| `author_major` | GSE241934 | paper `major.cell.type` |
| `author_collapsed` | GSE241934 | T+NK→TNK, B→TLS |
| `author_fine_immune` | GSE241934 | fine T/B/Plasma/Tfh vs major fallback for NA |
| `drmref` / `drmref_collapsed` | GSE207422 post only | DRMref 16-type (not CopyKAT) |
| `marker_coarse` | GSE207422, GSE253013 | TNK / TLS / epithelial / other |
| `author_garnett_limited` | GSE253013 | author epi + T only; residual = other |

TLS in dissociated scRNA is a **proxy** (B ± plasma ± Tfh), not a spatial
structure. Say so in any figure legend.

Malignant TACSTD2:

- GSE207422: epithelial **and** near-zero normal-lung markers
  (`SFTPA2, AGER, SCGB1A1, SCGB3A1, TPPP3`). Not CopyKAT.
- GSE241934: author **Epi** cells (paper did not deposit CopyKAT).
- GSE253013: sibling marker malignant-like score.

---

## 4. Model

For sample *i* and cell type *k*,

```
y_i ~ DirMult(n_i, π_i)
log(π_{ik} / π_{i,ref}) = β_{0k} + β_{1k} x_i
```

Reference = least-variable type with mean fraction ≥ 2%, preferring
`stromal` / `myeloid` / `other` / `mast` so TNK/TLS remain testable.
Shared precision φ; weak ridge on slopes (default 1.0) because n is small.
Failures → `not_estimable`, not p = 0.

Companions: ALR linear model (HC3) and fraction MWU/Spearman.

Optional: `sccoda.CompositionalAnalysis(..., reference_cell_type=...)`.

---

## 5. Honest n and FDR

- **n = patients** who pass the eligibility filter for that contrast.
- Cells are the DM library size, not replicates.
- BH **inside** each `fdr_family`. Report `n_tests_in_family`.
- Do not BH the companions together with the primary family.
- With n=9–15, a q<0.10 recovery is still a **small-n** result. Replicate
  across the grid; do not promote a single cohort×annotation hit.

---

## 6. How to run

```bash
pip install -r methods/scrna_sccoda/requirements.txt
python3 methods/scrna_sccoda/scripts/selftest_dm.py
python3 methods/scrna_sccoda/scripts/00_download.py
python3 methods/scrna_sccoda/scripts/01_extract.py --which all
python3 methods/scrna_sccoda/scripts/02_build_composition.py
python3 methods/scrna_sccoda/scripts/03_fit_grid.py
python3 methods/scrna_sccoda/scripts/04_figures.py
```

Large GEO files are gitignored. Results TSVs/JSON/PNG are the committed
artifacts.

---

## 7. What not to claim

- Author CopyKAT malignant IDs (not public for GSE207422).
- A TLS *structure* from dissociated scRNA.
- Cell-level Wilcoxon / Spearman as the composition test.
- MPR on GSE253013.
- That scCODA HMC ran if `summary.json` says `not_available`.
- A pooled “TACSTD2 excludes T/NK” claim unless the recovery table shows
  the same direction after FDR in more than one independent cohort.

---

## 8. References (verify DOIs at use)

- Büttner, Ostner et al. scCODA. *Nat Commun* 12, 6876 (2021).
- Aitchison J. *The Statistical Analysis of Compositional Data* (1986).
- Gloor et al. Microbiome datasets are compositional. *Front Microbiol* (2017).
- Squair et al. Confronting false discoveries in single-cell DE. *Nat Commun* (2021).
- Hu et al. GSE207422. *Genome Med* 15, 14 (2023).
- Zhang et al. GSE241934 / NEOTIDE. *Cell Rep Med* (2024), PMID 38897205.
- Sze et al. GSE253013. *Cancer Res* (2024), PMID 38335304.
