#!/usr/bin/env bash
# NicheNet v2 ligand-target matrix and ligand-receptor network (Zenodo 7074291).
set -euo pipefail
PRI="${1:-/tmp/nichenet_prior}"
mkdir -p "$PRI"
dl() {
  local dest="$1" url="$2"
  if [[ -s "$dest" ]] && [[ "$(stat -c%s "$dest")" -gt 1000 ]]; then
    echo "HAVE $dest"
    return 0
  fi
  curl -fL --retry 6 --retry-delay 8 --retry-all-errors -o "${dest}.partial" "$url"
  mv "${dest}.partial" "$dest"
  echo "OK $dest ($(stat -c%s "$dest") bytes)"
}
dl "$PRI/ligand_target_matrix_nsga2r_final.rds" \
  "https://zenodo.org/records/7074291/files/ligand_target_matrix_nsga2r_final.rds?download=1"
dl "$PRI/lr_network_human_21122021.rds" \
  "https://zenodo.org/records/7074291/files/lr_network_human_21122021.rds?download=1"
