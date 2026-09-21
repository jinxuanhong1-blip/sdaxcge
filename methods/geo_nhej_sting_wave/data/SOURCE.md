# Where the matrices come from

NCBI GEO, database `gds`, queried the day the script is run.

Series matrix and supplementary tables are read from `https://ftp.ncbi.nlm.nih.gov/geo/series/`. Platform gene maps are `*.annot.gz` on the same FTP host.

Matrices are not stored in this repository. `analyze.py` downloads each open table, scores it, and deletes the file.

The strict series query in the request is:

`(CLDN4 OR Cldn4) AND ("DNA repair" OR NHEJ OR PRKDC OR STING OR cGAS OR interferon) AND gse[Entry Type]`

NCBI maps bare `STING` toward the MeSH term "bites and stings" and bare `cGAS` toward chromogranin A. Those two tokens are counted, and the download queries use `STING1`, `TMEM173`, `"cGAS"`, and `"cyclic GMP-AMP"` so the MeSH expansion is not treated as a hit.
