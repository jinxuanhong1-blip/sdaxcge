# methods/pair_131907_148071_milo_cldn4

Additive **CLDN4-only** Milo-style kNN neighbourhood DA on the public
**GSE131907 + GSE148071 pair** (Kim et al. 2020; Wu et al. 2021).

This is **not** the triple (131907+148071+205335) and **not** the
131907+205335 pair. TACSTD2 is not a gate. Dual-high is not used.

**miloR is not used.** Same documented fallback as the single-cohort CLDN4
Milo folders: sample-level Spearman / Welch on neighbourhood proportions,
SpatialFDR with k-distance weights.

The independent unit is the **patient / sample**, not the cell and not the
overlapping neighbourhood. Graphs are built **per dataset / site**. Do not
pool tLung with mBrain, and do not pool neighbourhoods across GEO series
into one SpatialFDR.

Pairwise sample-level tests use GSE148071 biopsies + GSE131907 **tLung**
only. mBrain is a same-atlas sensitivity and is not added into the pair *n*.

```bash
pip install -r methods/pair_131907_148071_milo_cldn4/requirements.txt
python3 methods/pair_131907_148071_milo_cldn4/scripts/download.py
python3 methods/pair_131907_148071_milo_cldn4/scripts/run_pair.py
```

`FINDING.md` is written from the public-count run. The nhood table is
`tables/nhoods.tsv`.
