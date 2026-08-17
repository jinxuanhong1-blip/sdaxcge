# Harvested GSE19804 tumor scores

**Additive recut only.** This table is the 60-tumor extract already used for the
continuous CLDN4–CD8 purity-partial in the extra LUAD RNA table
(PR #235, `results/rework/A1_extra_luad/samples_GSE19804.tsv`,
branch `cursor/a1-luad-wave2-6484`). It is not a new GEO download.

| item | value |
|---|---|
| GEO | [GSE19804](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19804) |
| GDS | GDS3837 (120 arrays) |
| platform | GPL570 Affymetrix U133 Plus 2.0 |
| paper | Lu et al., *Cancer Epidemiol Biomarkers Prev* 2010, PMID 20802022 |
| arrays on series | 120 (60 tumor + 60 paired adjacent normal) |
| rows here | **60 tumors** (paired normals already dropped) |
| columns | `sample`, `TACSTD2`, `CLDN4`, `CD8A`, `GEP18`, ESTIMATE stromal / immune / score / TumorPurity |

Purity is ESTIMATE TumorPurity (ssGSEA Yoshihara 2013 141+141, then the
published cosine). Continuous Spearman / partial on this same table is
**already known** and is not the new claim. This folder only adds CLDN4
**Q4 vs Q1 vs CD8**.
