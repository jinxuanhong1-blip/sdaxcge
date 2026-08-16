# Methods — A9 extra-cohort TJ recurrence (supplement)

Additive public-data supplement. User A9 (intersection **CLDN1 / CLDN4 / CLDN7 / F11R / PARD3**) is taken as given and is **not** re-tested here. This analysis does not use that five-gene list as an inclusion filter, a ranking target, or a success criterion.

## Question

In three extra public LUAD RNA cohorts, which genes from a published tight-junction (TJ) catalog are higher in TACSTD2-high than in TACSTD2-low tumors, and which of those calls recur?

This is a supplement table of recurrent TJ-catalog genes. It is not a check of the slide gene list.

## Cohorts (public only)

| Cohort | Accession / freeze | Assay | Samples used | Why this cohort |
|--------|--------------------|-------|--------------|-----------------|
| OncoSG LUAD | cBioPortal `luad_oncosg_2020` (Chen et al. 2020) | RNA-seq V2 RSEM, all-sample z-scores | RNA sample list `luad_oncosg_2020_rna_seq_v2_mrna` | Independent East-Asian LUAD RNA-seq |
| CPTAC LUAD RNA | LinkedOmics / S3 `data_freeze_v1.2_reorganized/LUAD`, tumor RSEM | log2(RSEM coding UQ 1500) | All tumor columns in the tumor matrix | Independent proteogenomic LUAD RNA |
| GSE31210 LUAD | GEO GSE31210, GPL570 (Okayama et al. *Cancer Res* 2012) | Affymetrix U133 Plus 2.0 | `tissue: primary lung tumor` only | East-Asian stage I–II LUAD microarray |

GSE72094 was considered as the “one East-Asian GEO” example and **not used**: it is Moffitt (US) LUAD (Schabath et al.), not an East-Asian series. GSE31210 is the East-Asian GEO.

No dbGaP, EGA, or other controlled-access data.

OncoSG public cBioPortal RNA is z-scored. Spearman is rank-invariant. The tertile high-minus-low difference is in z units, not log2 RNA. Direction is still used for the call.

## TJ catalog (pre-specified)

Universe = union of:

1. **GO:0005923** bicellular tight junction (Enrichr GO_Cellular_Component_2021; same GMT as the sibling A8/A9 GSEA)
2. **R-HSA-420029** Tight Junction Interactions (Enrichr Reactome_2022; same GMT)
3. **JAM2** only, as a documented canonical JAM present in KEGG Tight junction but missing from (1) and (2)

TACSTD2 is the splitter and is not a tested gene.

Retired catalog symbols are queried under the current HGNC name and renamed
back to the catalog symbol: **BVES → POPDC1**, **DDX58 → RIGI**. CPTAC RNA
rows are Ensembl gene IDs; they are mapped with the public HGNC custom table
(approved symbol + previous + alias → Ensembl gene ID; version suffix stripped).

KEGG Tight junction actin / myosin / MAPK / tubulin members are **not** added. They are pathway genes, not a TJ-gene catalog.

The A9 slide genes are annotated in the output table when they occur in this catalog. They do not define who is tested.

## Split and tests

- Within each cohort, TACSTD2-high = upper tertile; TACSTD2-low = lower tertile. Mid tertile is unused.
- Association: Spearman ρ of each TJ gene vs continuous TACSTD2.
- High vs low: Welch two-sample t-test on the cohort native scale. Effect = mean(high) − mean(low).
- Multiple testing: Benjamini–Hochberg FDR across TJ-catalog genes **present in that cohort**.
- **Call** `up_in_tacstd2_high` only if ρ > 0 **and** high−low > 0 **and** Welch BH-FDR < 0.05.

Median split is a labeled sensitivity and does not define recurrence.
The tertile 3/3 set and the median 3/3 set are both reported; the larger
median set is not substituted for the pre-specified tertile table.

GSE31210 probes: GPL570 `.annot` gene symbols; probes mapping to more than one symbol (`///`) are dropped; remaining probes for a gene are collapsed by **max-mean** probe. The GEO series matrix is MAS5 linear intensity; the compact extract is **log2(MAS5 + 1)** so the high-vs-low difference is a log2 fold-change. Spearman is rank-invariant to that transform.

## Recurrence (honest counts)

Let P_c be the set of catalog genes that pass in cohort c.

| Count | Definition |
|-------|------------|
| n_pass(c) | \|P_c\| |
| pairwise | \|P_a ∩ P_b\| for each pair |
| recur_2of3 | genes in at least two of the three P_c |
| recur_3of3 | genes in P_OncoSG ∩ P_CPTAC ∩ P_GSE31210 |

Empty intersections are reported as empty. They are not replaced with a weaker cutoff, a different split, or a smaller gene list.

A gene that is absent from a platform (no unique probe / not in the RNA matrix) is `absent` in that cohort and cannot pass there. Absence is not a no-call on measured data.

## What is not done

- No tuning of tertiles, FDR, or lineage to enlarge the overlap.
- No TCGA in this supplement (TCGA was the A9 discovery setting).
- No protein, spatial, or single-cell analysis.
- No claim that TROP2 binds any recurrent gene.

## Reproduce

```bash
python3 scripts/w200/A9_extra_tj_intersection/download.py
python3 scripts/w200/A9_extra_tj_intersection/analyze.py
```

Outputs: `results/w200/A9_extra_tj_intersection/`
(`README.md`, `summary.json`, `tables/supplement_tj_recurrence.tsv`).
