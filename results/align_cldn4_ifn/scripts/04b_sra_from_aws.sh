#!/usr/bin/env bash
# Download GSE207704 SRA from AWS Open Data and quantify with salmon.
# ENA FASTQ 404'd and NCBI fasterq-dump over the network segfaulted;
# local SRA files from sra-pub-run-odp are the remaining path.
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
WD=/tmp/gse207704
REF=$WD/ref
mkdir -p "$WD/sra" "$WD/fastq" "$WD/quant" "$ROOT/data/GSE207704_salmon" "$ROOT/logs"
export MAMBA_ROOT_PREFIX=/tmp/mamba
MM="/tmp/bin/micromamba run -n rnaseq"

TX=$REF/gencode.v44.transcripts.fa.gz
if [ ! -s "$TX" ]; then
  curl -fsSL --retry 4 -o "$TX" \
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.transcripts.fa.gz || exit 1
fi
if [ ! -d "$REF/salmon_idx" ]; then
  $MM salmon index -t "$TX" -i "$REF/salmon_idx" -k 23 --gencode -p 4 \
      > "$ROOT/logs/salmon_index.log" 2>&1 || exit 1
fi

# GSM mapping (from SRA runinfo)
# SRR20029125 T47D WT rep1
# SRR20029124 T47D WT rep2
# SRR20029123 T47D CLDN4-/- rep1
# SRR20029122 T47D CLDN4-/- rep2
# SRR20029121 MCF7 WT rep1
# SRR20029120 MCF7 WT rep2
# SRR20029119 MCF7 CLDN4-/- rep1
# SRR20029118 MCF7 CLDN4-/- rep2
runs="SRR20029125 SRR20029124 SRR20029123 SRR20029122 SRR20029121 SRR20029120 SRR20029119 SRR20029118"

for r in $runs; do
  sra=$WD/sra/$r
  if [ ! -s "$sra" ]; then
    echo "[dl] $r from AWS ODP"
    curl -fL --retry 4 --max-time 3600 \
      -o "$sra" "https://sra-pub-run-odp.s3.amazonaws.com/sra/$r/$r" \
      || { echo "[FAIL dl] $r"; continue; }
  fi
  fq=$WD/fastq/$r.fastq.gz
  if [ ! -s "$fq" ]; then
    echo "[dump] $r"
    ( cd "$WD/fastq" && $MM fasterq-dump --split-spot -e 4 -O . "$sra" \
        && $MM pigz -p 4 "$r.fastq" ) \
      || { echo "[FAIL dump] $r"; continue; }
  fi
  if [ ! -s "$WD/quant/$r/quant.sf" ]; then
    echo "[quant] $r"
    $MM salmon quant -i "$REF/salmon_idx" -l A -r "$fq" -p 4 \
        --validateMappings --seqBias --gcBias -o "$WD/quant/$r" \
        > "$WD/quant/$r.log" 2>&1 || echo "[FAIL quant] $r"
  fi
  if [ -s "$WD/quant/$r/quant.sf" ]; then
    cp "$WD/quant/$r/quant.sf" "$ROOT/data/GSE207704_salmon/$r.quant.sf"
    echo "[ok] $r"
    rm -f "$fq" "$sra"   # reclaim disk after success
  fi
done
echo "[done] $(ls "$ROOT/data/GSE207704_salmon"/*.quant.sf 2>/dev/null | wc -l) quants"
