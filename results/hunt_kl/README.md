# Hunt: public KL / immune-resistant lung scRNA with a TROP2-high, immune-low tumor subset

Goal: find **public** KRAS/LKB1 (KL) or immune-resistant mouse/human lung single-cell RNA
datasets in which one could identify a **TROP2 (gene `TACSTD2` / mouse `Tacstd2`)-high,
immune-low tumor subset**. The user's private 8-KL-mouse cohort is explicitly off limits and was
not used. Everything here is from GEO.

## TL;DR (read this first — it is the honest bottom line)

1. **`GSE76628` (one of the two example accessions in the request) is NOT a lung KL scRNA dataset.**
   GEO record `GSE76628` is *"Stromal-Based Signatures for the Classification of Gastric Cancer"* —
   mouse **gastric** cancer, **microarray** (Affymetrix, expression profiling by array), not
   single-cell and not lung/KL. It does not match the brief. The user likely transposed the
   accession; I did not silently substitute a "corrected" one. See `Reference-dataset check` below.
2. **`GSE180963` is exactly on target**: mouse `KrasG12D` vs `KrasG12D;Lkb1fl/fl` (KL) lung GEMM,
   scRNA-seq of dissociated tumor nodules, explicitly framed around the LKB1-loss **immune-desert**
   TME. Small (n=2 series-level).
3. **No public dataset states a "TROP2-high immune-low subset" in its abstract.** TROP2/`TACSTD2`
   is an expression-level feature, essentially never mentioned in KL/immune-cold lung abstracts.
   The claim can only be established by **reanalysis**, so I ran a targeted check on the most
   tractable KL scRNA matrices instead of asserting it.
4. **What I verified with real data:** in **`GSE179502`** (sorted neoplastic KL/`Lkb1` lung tumor
   cells) TROP2/`Tacstd2` is heterogeneously expressed in **54.7%** of tumor cells with a clear
   **HIGH subset (~14% of tumor cells)**. So a TROP2-high tumor subset **demonstrably exists** in a
   real KL model. The paired "immune-low" property could not be tested there because the cells were
   pre-sorted to tumor only; that half needs a tumor+immune co-capture dataset.

## Recommended datasets (see `curated_shortlist.tsv` for the full table + URLs)

Tier 1 — directly usable for the KL + TROP2 question:

| Accession | Organism | Model | Assay | Why |
|---|---|---|---|---|
| `GSE180963` | mouse | KrasG12D vs KL lung GEMM | scRNA whole nodules | on-target reference; tumor+immune co-captured; small n |
| `GSE179502` | mouse | KT;Lkb1 restorable, sorted tumor cells | scRNA (10x) | **verified**: distinct Tacstd2-high tumor subset (~14%) |
| `GSE179501` | mouse | KT;Lkb1 restorable, total viable | scRNA (10x) | tumor+immune; verified but noisy (see caveat) |
| `GSE280232` | human | resectable KRAS-mut, STK11-comut vs wt | scRNA + scTCR | strongest **human** KL analog with immune outcomes |

Tier 2 — immune-resistant lung scRNA or KL immune compartment (adjacent):
`GSE194166` (KPL immune-cell scRNA, immune-cold TME), `GSE274934` (ALK+ NSCLC, poor-immunotherapy,
CITE-seq), `GSE198099` (immune-suppressive NSCLC), `GSE146100` (immunotherapy-response LUAD).

Tier 3 — KL / immune-cold **bulk** context for signature validation (not single-cell):
`GSE274351` (K/KP/KL LCM), `GSE277929` (T-cell exclusion in KL), `GSE244452`, `GSE175479`,
`GSE295685`, `GSE302247`.

## How this list was built (reproducible)

1. **Search** — 10 GEO (`db=gds`) queries via NCBI E-utilities covering the KL genotype
   (`Kras` AND `Lkb1`/`Stk11`), lung, single-cell, immune-cold/immunotherapy-resistant, and
   TROP2/`TACSTD2` axes → 164 unique GSE-level records (`verification/uids.json`).
2. **Rank** — each record scored on genotype (KL/LKB1/KRAS), lung, single-cell, immune-cold vs
   generic-immune, and TROP2 mentions → `candidates_all.tsv` (all 164, sorted).
3. **Assay truthing** — for the top hits I pulled the GEO `Series_type` and `overall_design` and
   reclassified them, because several series titles say "scRNA" while the relevant subseries are
   actually **bulk** or **immune-sorted only** (e.g. `GSE244452`, `GSE175479`, `GSE302247`,
   `GSE277929` are bulk; `GSE194166` is immune-only). This is captured in `curated_shortlist.tsv`.
4. **Data-level verification** — downloaded the KL scRNA count matrices for `GSE179501` and
   `GSE179502` and measured `Tacstd2`/immune-marker expression directly
   (`verification/verify_GSE179501.txt`, `verify_GSE179502.txt`).

## Reference-dataset check (verbatim GEO metadata)

- `GSE76628` — *Stromal-Based Signatures for the Classification of Gastric Cancer [part II]*;
  organism **Mus musculus**; type **Expression profiling by array**; 78 samples; 2016. → gastric,
  microarray, not lung/KL/scRNA. **Mismatch — flagged, not used.**
- `GSE180963` — *Single cell RNA sequencing of tumor sections from GEMM harboring KrasG12D/+ or
  KrasG12D/+Lkb1fl/fl (KL) mutation*; **Mus musculus**; high-throughput sequencing; 2022; summary
  explicitly discusses LKB1-loss immune-desert/"cold" TME. **On target.**

## Verification results (real numbers, honest caveats)

**`GSE179502`** (sorted neoplastic KL tumor cells; 16,017 cells): 73.1% `Epcam+`, 0.7% `Ptprc+`
(clean tumor epithelium). **`Tacstd2`/TROP2 positive in 54.7%** of tumor cells; top-quartile
"TROP2-high" subset = **13.7%** of tumor cells. → a TROP2-high tumor subset **exists** in KL tumor
cells. Immune-low pairing not testable here (no immune cells retained).

**`GSE179501`** (total-viable KL cells; 20,351 cells): immune-dominated (63.2% `Ptprc+`, 13% `Epcam+`,
3.4% `Tacstd2+`). `Tacstd2`-high vs `Tacstd2`-negative immune-score difference is **weak**
(2.11 vs 2.39; Pearson r = **−0.087**), and 63% of `Epcam+` cells are also `Ptprc+`, indicating
ambient-RNA/doublet contamination. → consistent with, but does not on its own prove, a TROP2-high
immune-low subset; **proper clustering + doublet QC required**.

### Method caveats (do not over-read the verification)
- Normalization is CP10k + log1p with **crude single-gene thresholds**; no clustering, no doublet
  removal, no batch integration, no immune deconvolution. These are **feasibility/plausibility
  checks**, not a validated subset call.
- "Immune-low" was proxied by lineage markers (`Ptprc`, `Cd3e`, `Cd8a`, `Cd68`, `Nkg7`). A rigorous
  analysis should use spatial or neighborhood/co-capture data and a full immune signature.
- I did **not** re-run every dataset; Tier-2/3 assignments rely on GEO metadata, not reanalysis.

## What is NOT claimed
- I did **not** confirm a "TROP2-high immune-low tumor subset" as a published finding anywhere; no
  abstract states it. The strongest honest statement is: a TROP2-high KL **tumor** subset is real
  and reproducible (`GSE179502`), and the right public datasets exist to test the immune-low pairing
  (`GSE180963`, `GSE280232`).
- No private data were accessed.

## Files
- `curated_shortlist.tsv` — hand-verified tiered shortlist with assay truthing and URLs.
- `candidates_all.tsv` — all 164 scored GEO hits (machine-readable).
- `verification/verify_GSE179502.txt`, `verify_GSE179501.txt` — real expression measurements.
- `verification/recs.json`, `uids.json` — raw GEO metadata + search UIDs for reproducibility.

## Suggested next step
Take `GSE180963` (mouse KL, tumor+immune co-capture) and/or `GSE280232` (human KRAS/STK11), run
standard Seurat/Scanpy QC + clustering, annotate the tumor epithelial compartment, score
`Tacstd2`/`TACSTD2`, and test whether the TROP2-high tumor cluster is depleted of neighboring
immune cells (or has low cytotoxic infiltration by deconvolution). Use `GSE179502` as the positive
control that a TROP2-high tumor state exists in KL.
