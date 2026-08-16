#!/usr/bin/env bash
# Download public inputs for results/w200/Xena_LUAD/.
# All files land in data/w200/raw/ and are checksummed into MANIFEST.tsv.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RAW="${1:-$ROOT/data/w200/raw}"
MANIFEST="$RAW/MANIFEST.tsv"
mkdir -p "$RAW"

fetch () {
  local out="$RAW/$1" url="$2"
  if [[ -s "$out" ]]; then echo "[skip] $1"; return 0; fi
  echo "[get ] $1"
  for attempt in 1 2 3 4 5; do
    if curl -fsSL --retry 3 --retry-delay 2 --connect-timeout 30 --max-time 3600 -o "$out.part" "$url"; then
      mv "$out.part" "$out"; return 0
    fi
    echo "       retry $attempt after failure"; sleep $((2 ** attempt))
  done
  echo "[FAIL] $1 <- $url" >&2; rm -f "$out.part"; return 1
}

# UCSC Xena GDC hub (STAR, log2(TPM+1), Ensembl IDs)
fetch TCGA-LUAD.star_tpm.tsv.gz \
  https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LUAD.star_tpm.tsv.gz
fetch gencode.v36.annotation.gtf.gene.probemap \
  https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap
fetch TCGA-LUAD.clinical.tsv.gz \
  https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LUAD.clinical.tsv.gz

# UCSC Xena legacy TCGA hub (RSEM HiSeqV2) — matched to MDACC ESTIMATE RNAseqV2
fetch TCGA.LUAD.HiSeqV2.gz \
  https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap/HiSeqV2.gz
fetch LUAD_clinicalMatrix \
  https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap/LUAD_clinicalMatrix

# Purity
fetch TCGA_mastercalls.abs_tables_JSedit.fixed.txt \
  https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5
fetch Aran2015_purity_S1.xlsx \
  'https://static-content.springer.com/esm/art%3A10.1038%2Fncomms9971/MediaObjects/41467_2015_BFncomms9971_MOESM1236_ESM.xlsx'
fetch MDACC_estimate_LUAD_RNAseqV2.txt \
  https://bioinformatics.mdanderson.org/estimate/tables/lung_adenocarcinoma_RNAseqV2.txt

# PanImmune (Thorsson 2018) — methylation leukocyte fraction is RNA-independent
fetch TCGA_all_leuk_estimate.masked.20170107.tsv \
  https://api.gdc.cancer.gov/data/6f75c9d7-5134-4ed1-b8f3-72856c98a4e8
fetch TCGA.Kallisto.fullIDs.cibersort.relative.tsv \
  https://api.gdc.cancer.gov/data/b3df502e-3594-46ef-9f94-d041a20a0b9a
fetch Scores_160_Signatures.tsv.gz \
  https://api.gdc.cancer.gov/data/80a82092-161d-4615-9d96-e858f113618d

{
  printf 'file\tbytes\tsha256\n'
  for f in "$RAW"/*; do
    [[ "$(basename "$f")" == "MANIFEST.tsv" ]] && continue
    [[ -f "$f" ]] || continue
    printf '%s\t%s\t%s\n' "$(basename "$f")" "$(stat -c%s "$f")" "$(sha256sum "$f" | cut -d' ' -f1)"
  done
} > "$MANIFEST"
echo "--- manifest ---"; cat "$MANIFEST"
