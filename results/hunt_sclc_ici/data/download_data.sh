#!/usr/bin/env bash
# Download the PUBLIC datasets used in this analysis.
# Raw data is intentionally NOT committed to git (size + redistribution terms).
# Run this script from results/hunt_sclc_ici/data/ to reproduce.
set -euo pipefail
cd "$(dirname "$0")"

echo "[1/2] George et al. 2015 SCLC bulk RNA-seq (cBioPortal datahub, PUBLIC)"
BASE="https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/sclc_ucologne_2015"
for f in data_mrna_seq_fpkm.txt data_clinical_sample.txt data_clinical_patient.txt; do
  [ -f "$f" ] || curl -sSL -o "$f" "$BASE/$f"
  echo "  - $f"
done

echo "[2/2] Chan et al. 2021 SCLC atlas (cellxgene, PUBLIC) — combined object (~1.5 GB)"
# Full collection: https://cellxgene.cziscience.com/collections/62e8f058-9c37-48bc-9200-e767f318a8ec
CHAN_URL="https://datasets.cellxgene.cziscience.com/a9d92e38-9a6e-401b-8484-74bb15122341.h5ad"
[ -f chan_combined.h5ad ] || curl -sSL -o chan_combined.h5ad "$CHAN_URL"
echo "  - chan_combined.h5ad"

echo "Done. NOTE: IMpower133 per-patient RNA-seq is NOT downloadable here —"
echo "it is CONTROLLED ACCESS via EGA study EGAS00001004888 (see IMPOWER133_NOTE.md)."
