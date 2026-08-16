# Hunt GSE131907 — epithelial TACSTD2/CLDN4 vs CD8/NK/TLS; TJ/keratin GSEA

Pipeline is in `scripts/hunt_gse131907/`. This file is overwritten with the
honest verdict when `analyze.py` finishes.

**Not spatial.** GSE131907 is dissociated 10x. Neighborhood = sample-level
co-occurrence. TLS = B / GC-B / 12-chemokine proxies.

**Epithelial-restricted.** TACSTD2 and CLDN4 used for association and GSEA
grouping are scored only in author epithelial cells.
