# Data provenance & methods — TACSTD2/CLDN4 in EGFR/ALK lung cancer

This slice studies **TACSTD2 (TROP2)** and **CLDN4 (Claudin-4)** — two antibody–drug
conjugate (ADC) / tight-junction targets — against the anti-tumour immune
microenvironment and immune-checkpoint-inhibitor (ICI) outcome in **EGFR/ALK-mutant
lung cancer**.

All raw inputs are **public** and downloaded on the fly by `scripts/fable_egfr_alk/run_all.sh`
into `data/fable_egfr_alk/` (git-ignored). Only compact, processed outputs are committed
under `results/fable_egfr_alk/` (total < 2 MB, far below the 2 GB cap).

## Datasets

| ID | Source | n | Assay | Key labels | Role |
|----|--------|---|-------|-----------|------|
| TCGA-LUAD | UCSC Xena GDC hub (`gdc-hub.s3`) | 528 primary tumours | STAR TPM (log2 TPM+1) + WXS mutations + survival | EGFR-mut / ALK-driven / WT (derived) | **Immune axis (primary, EGFR/ALK-stratified)** |
| TCGA-LUSC | UCSC Xena GDC hub | 501 primary tumours | STAR TPM + WXS + survival | same | Immune-axis replication |
| GSE126044 | GEO (Cho et al. 2020) | 16 | RNA-seq counts | responder / non-responder to anti-PD-1 | ICI response axis |
| GSE135222 | GEO (Jung et al. 2019, Nat Commun) | 27 | RNA-seq (TPM/FPKM) | PFS event + time under anti-PD-1/PD-L1 | ICI PFS axis |

Exact download URLs are encoded in `run_all.sh`.

## Driver-group definition (TCGA)

- **EGFR-mut** — ≥1 non-silent somatic mutation in `EGFR` (Xena WXS/MC3-style calls).
  71/528 LUAD (13.4%), consistent with the literature (~11–15%).
- **ALK-driven** — EML4-ALK fusion **proxy**: samples with aberrantly high `ALK`
  expression (> median + 3·MAD **and** > 95th percentile). Fusion is not directly in
  the SNV file; ALK is normally silent in lung tissue and fusion drives strong
  3′-kinase-domain expression, so the outlier flag is a recognised surrogate.
  25/528 LUAD (4.7%), consistent with ~3–7% ALK-rearranged LUAD.
  ALK **point** mutations are deliberately *not* used to define drivers (ALK is a
  large gene; most non-silent SNVs in this smoking-enriched cohort are passengers);
  they are retained only as a transparency column.
- **EGFR_ALK_wt** — WXS-profiled, neither of the above.

## Expression handling

- TCGA: Xena STAR `star_tpm` values are already `log2(TPM+1)`; used directly.
- GSE126044: raw counts → log2(CPM+1).
- GSE135222: supplied TPM/FPKM → log2(x+1).
- Ensembl IDs matched to symbols after stripping the `.NN` version suffix.

## Immune signatures (transparent, marker-based)

Each signature score = **mean of per-sample z-scored genes** (log scale), so it is
robust and needs no external reference. Definitions in `common.py`:

- `CYT` — cytolytic activity (GZMA, PRF1; Rooney et al. 2015).
- `IFNG_6gene` — IFN-γ program (IFNG, STAT1, IDO1, CXCL9, CXCL10, HLA-DRA; after Ayers et al. 2017).
- `CD8_effector` — CD8A/B, GZMA/B/K, PRF1, NKG7, IFNG.
- `Checkpoint` — CD274, PDCD1, CTLA4, LAG3, HAVCR2, TIGIT.
- `Broad_immune` — T-cell + APC + myeloid markers (mixes in CD68/ITGAX myeloid signal).

## Statistics

- Associations: **Spearman** correlation (rank-based, outlier-robust).
- Group differences: **Mann–Whitney U** (+ rank-biserial effect size), **Kruskal–Wallis** across groups.
- Multiplicity: **Benjamini–Hochberg** FDR within each target × group.
- Survival: Kaplan–Meier + **log-rank** on median split (`lifelines`).
- Verification (`03_verify.py`): 5,000× bootstrap 95% CIs; 10,000× permutation p-values;
  cross-cohort sign concordance (binomial test); alternative immune scorings; ICI positive control.
- Seed: `20260816` (see `common.RANDOM_SEED`).

## Key limitation (stated honestly)

The two ICI RNA-seq cohorts carry **no EGFR/ALK genotype**, so EGFR/ALK-mutant
*specificity* is established in TCGA (immune microenvironment), while ICI *response*
is assessed in overall NSCLC-ICI cohorts. The direction of every effect is
consistent across the two axes, but the ICI cohorts (n=16, n=27) are underpowered
for target-level significance. This is the honest scope of a "parallel slice".
