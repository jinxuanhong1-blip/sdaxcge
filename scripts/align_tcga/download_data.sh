#!/usr/bin/env bash
# Download all PUBLIC OPEN-ACCESS matrices used by the USER-ALIGN replication.
# Nothing here is controlled-access. Re-running is idempotent.
set -euo pipefail
DEST="$(cd "$(dirname "$0")/../../results/align_tcga/data" && pwd)"
cd "$DEST"

echo "[1/3] TCGA-LUAD expression (UCSC Xena, RSEM log2(norm+1)) ..."
curl -fsSL -o LUAD.HiSeqV2.gz \
  "https://tcga.xenahubs.net/download/TCGA.LUAD.sampleMap/HiSeqV2.gz"

echo "[2/3] TCGA-LUSC expression (UCSC Xena, RSEM log2(norm+1)) ..."
curl -fsSL -o LUSC.HiSeqV2.gz \
  "https://tcga.xenahubs.net/download/TCGA.LUSC.sampleMap/HiSeqV2.gz"

echo "[3/3] TCGA tumour purity (Aran et al. Nat Commun 2015, Supp Data 1) ..."
curl -fsSL -o Aran_CPE_purity.xlsx \
  "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fncomms9971/MediaObjects/41467_2015_BFncomms9971_MOESM1236_ESM.xlsx"

echo "Done. Files in: $DEST"
ls -la "$DEST"
