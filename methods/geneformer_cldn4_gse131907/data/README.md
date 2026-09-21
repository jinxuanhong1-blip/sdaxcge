`a8_sets.json` is the gene-set file already used by the concordant-4 CLDN4 analysis (Hallmark, KEGG tight junction, GOBP tight-junction organization, the custom MHC-I/APM list, and the focal genes). It is copied here so this script does not depend on another branch being checked out.

GEO matrices and model weights are not stored in the repo. `analyze.py` reads them from `GF_CACHE` (default `/tmp/gf_cache`).
