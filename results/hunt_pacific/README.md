# Hunt: PACIFIC/durvalumab NSCLC data for purity-adjusted TACSTD2 (TROP2) vs immune

**Question.** Find PACIFIC/durvalumab NSCLC supplements, NanoString, or open RNA where a
purity-adjusted correlation between `TACSTD2` and immune signatures can be computed.
User hypothesis: negative Spearman ρ.

**Honesty policy.** Every accession/ID below was verified live against NCBI eutils, GEO FTP,
PubMed, GDC, or the cBioPortal API before inclusion. Datasets that looked relevant but were
unusable are reported as such, not silently dropped.

## 1. What exists (and what does not)

| Resource | Verdict |
|---|---|
| **PACIFIC trial (NCT02125461)** | **No open per-patient RNA.** Individual patient-level data (including any omics) are controlled-access via AstraZeneca's data-sharing program on [Vivli](https://vivli.org/ourmember/astrazeneca/). Published PACIFIC biomarker supplements are PD-L1 IHC subgroup analyses only. |
| **SUBMARINE (WJOG11518L)**, [PMID 37364849](https://pubmed.ncbi.nlm.nih.gov/37364849/) | **Closest PACIFIC-regimen NanoString dataset that exists**: nCounter IO360 on ~85 pre-CRT tumors from stage III NSCLC treated with dCCRT + durvalumab. **Not deposited in GEO** (searches for `WJOG11518L`/`SUBMARINE` return no relevant series; the single hit, GSE278471, is an archaeal study — false positive). Author request would be required. |
| **GSE253564** | Open. Pre-treatment tumor RNA-seq FPKM (n=32) from a randomized phase II neoadjuvant durvalumab ± SBRT trial (PMID 38401548). TACSTD2 present. **Analyzed.** |
| **GSE248378** | Open. Post-durvalumab resected tumors, RNA-seq FPKM (n=29), same trial (PMIDs 38114518, 38401548). TACSTD2 present. **Analyzed.** |
| **GSE110390** | Durvalumab Study 1108 NSCLC — only the 21 IFNγ-signature genes were deposited; **TACSTD2 absent (verified) → unusable.** |
| **TCGA-LUAD / TCGA-LUSC (PanCanAtlas)** | Open RNA where purity adjustment is fully rigorous: DNA-based ABSOLUTE purity (GDC PanCanAtlas file `TCGA_mastercalls.abs_tables_JSedit.fixed.txt`) is independent of expression. Treatment-naive — a biological reference, not durvalumab tissue. **Analyzed.** |

Caveat on the two open durvalumab sets: they are neoadjuvant, stages I–III (resectable), not the
PACIFIC stage III unresectable CRT population. They are the only open durvalumab-NSCLC tumor
transcriptomes found.

## 2. Results — TACSTD2 vs immune, Spearman ρ

Full tables: `tacstd2_immune_correlations.csv` (GEO sets), `tacstd2_immune_correlations_tcga.csv`
(TCGA); per-sample scores in `scores_*.csv`; run logs in `analysis_log*.txt`.

### GSE253564 — pre-treatment, durvalumab trial (n=32)

| Immune score | ρ (unadj.) | p | ρ adj. ESTIMATE | p | ρ adj. stromal-only | p |
|---|---|---|---|---|---|---|
| ESTIMATE ImmuneScore | −0.72 | <0.001 | −0.53* | 0.002 | −0.71 | <0.001 |
| IFNγ-6 (Ayers) | −0.29 | 0.11 | +0.12 | 0.54 | −0.24 | 0.19 |
| GEP-18 T-cell-inflamed | −0.44 | 0.012 | −0.02 | 0.89 | −0.40 | 0.024 |
| CYT (GZMA/PRF1) | −0.46 | 0.008 | −0.17 | 0.37 | −0.45 | 0.011 |
| CD8A | −0.54 | 0.002 | −0.24 | 0.19 | −0.51 | 0.003 |

### GSE248378 — post-durvalumab resected tumors (n=29)

| Immune score | ρ (unadj.) | p | ρ adj. ESTIMATE | p | ρ adj. stromal-only | p |
|---|---|---|---|---|---|---|
| ESTIMATE ImmuneScore | −0.67 | <0.001 | −0.47* | 0.011 | −0.59 | <0.001 |
| IFNγ-6 (5/6 genes; IFNG missing) | −0.49 | 0.007 | −0.24 | 0.22 | −0.39 | 0.039 |
| GEP-18 T-cell-inflamed | −0.65 | <0.001 | −0.43 | 0.023 | −0.57 | 0.002 |
| CYT (GZMA/PRF1) | −0.81 | <0.001 | −0.71 | <0.001 | −0.77 | <0.001 |
| CD8A | −0.71 | <0.001 | −0.58 | 0.001 | −0.66 | <0.001 |

\* Partially circular: the covariate (ESTIMATE score) contains the ImmuneScore.

### TCGA (ABSOLUTE DNA purity — no circularity)

| Dataset | Immune score | ρ (unadj.) | ρ purity-adj. | adj. p |
|---|---|---|---|---|
| LUSC (n=479) | IFNγ-6 | −0.19 | **−0.22** | 8.5e−7 |
| LUSC | GEP-18 | −0.18 | **−0.23** | 3.0e−7 |
| LUSC | CYT | −0.13 | **−0.17** | 1.7e−4 |
| LUSC | CD8A | −0.21 | **−0.25** | 4.0e−8 |
| LUAD (n=497) | IFNγ-6 | −0.02 | −0.02 | 0.68 |
| LUAD | GEP-18 | −0.03 | −0.03 | 0.53 |
| LUAD | CYT | −0.11 | −0.12 | 0.008 |
| LUAD | CD8A | −0.10 | −0.10 | 0.020 |

## 3. Honest interpretation vs the target (ρ negative)

- **Direction supported, magnitude context-dependent.** All 26 unadjusted TACSTD2–immune
  correlations across four datasets are negative (one adjusted estimate flips to a
  non-significant +0.12).
- **Strongest and adjustment-robust**: post-durvalumab resected tumors (GSE248378; CYT adj.
  ρ = −0.71) and TCGA-LUSC, where DNA-purity-adjusted ρ is negative for all four scores and
  slightly *strengthens* after adjustment (p ≤ 1.7e−4).
- **Weakest**: TCGA-LUAD (near zero for IFNγ/GEP-18; small negative for CYT/CD8A) and the
  pre-treatment durvalumab cohort under the most conservative (full-ESTIMATE) adjustment,
  where T-cell scores drop to non-significance at n=32. Under stromal-only adjustment they
  remain negative and mostly significant. With n≈30 the two GEO sets are hypothesis-level
  evidence, not confirmation.
- In both GEO sets TACSTD2 correlates *positively* with the purity proxy (ρ ≈ 0.54–0.59), i.e.
  part of the raw negative TACSTD2–immune correlation is a purity artifact; the TCGA-LUSC result
  shows a genuine residual negative association beyond purity.

## 4. Methods (brief)

log2(FPKM+1) (GEO) or log2(RSEM+1) (TCGA, cBioPortal API). Immune scores: ESTIMATE ImmuneScore
(ssGSEA, official 141-gene sets from `estimate` R pkg v1.0.13 on R-Forge), Ayers IFNγ-6 and
18-gene T-cell-inflamed GEP (unweighted mean approximation of the weighted NanoString original),
Rooney CYT, CD8A. Purity: GEO sets — ESTIMATE; because the Affymetrix-calibrated cos-formula is
miscalibrated on RNA-seq, partial Spearman uses the rank of the raw ESTIMATE score (invariant to
monotone transforms), with a stromal-only sensitivity covariate; TCGA — ABSOLUTE purity.
Partial Spearman = Pearson on rank residuals, p from t with n−3 df.
Reproduce: `scripts/fetch_data.sh` then `scripts/hunt_pacific_tacstd2.py` and
`scripts/tcga_tacstd2_validation.py`.

## 5. Routes to actual PACIFIC(-regimen) data

1. **Vivli request** to AstraZeneca for PACIFIC (NCT02125461) individual-level data.
2. **Author request** for SUBMARINE (WJOG11518L) nCounter IO360 data (corresponding author
   K. Haratani, Kindai University; PMID 37364849) — the closest PACIFIC-regimen NanoString set.
