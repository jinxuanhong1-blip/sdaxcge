#!/usr/bin/env bash
# =============================================================================
# geo_query_edirect.sh  -  find GEO Series with NCBI EDirect (esearch/efetch).
# -----------------------------------------------------------------------------
# EDirect is NCBI's official command-line toolkit. Install once:
#   sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
#   export PATH=${HOME}/edirect:${PATH}
# Set an API key to raise the rate limit:
#   export NCBI_API_KEY=xxxxxxxx
#
# All accessions returned are REAL (queried live from GEO). Never fabricate them.
# =============================================================================
set -euo pipefail

# ---- 1. gene-perturbation searches ----------------------------------------
# Use fielded queries so you only get sequencing Series in the right organism.
search () {
  local term="$1"
  echo "########## ${term} ##########"
  esearch -db gds -query "${term} AND \"expression profiling by high throughput sequencing\"[DataSet Type] AND gse[Entry Type]" \
    | esummary \
    | xtract -pattern DocumentSummary -element Accession taxon n_samples title
  echo
}

search 'CLDN4[Description] AND (knockdown OR knockout OR CRISPR OR shRNA OR siRNA)'
search 'TACSTD2[Description] AND (knockdown OR knockout OR shRNA OR siRNA)'
search 'TROP2 AND (knockdown OR knockout)'

# ---- 2. inspect one Series' processed supplementary files -----------------
# Pick an accession from above, then list its suppl/ dir + file sizes so you can
# choose a processed count matrix under your size budget (e.g. < 2 GB).
inspect_suppl () {
  local acc="$1"                       # e.g. GSE207704
  local stub="${acc%???}"              # GSE207704 -> GSE207
  local url="https://ftp.ncbi.nlm.nih.gov/geo/series/${stub}nnn/${acc}/suppl/"
  echo "### suppl files for ${acc} ###"
  curl -s "$url" \
    | grep -oP '(?<=href=")[^"]+' \
    | grep -viP '^/geo|^http|Parent' \
    | while read -r f; do
        sz=$(curl -sI "${url}${f}" | awk 'tolower($1)=="content-length:"{print $2}' | tr -d '\r')
        printf '  %s  (%s bytes)\n' "$f" "${sz:-?}"
      done
}

# Example (uncomment to run against a real series you selected):
# inspect_suppl GSE207704
