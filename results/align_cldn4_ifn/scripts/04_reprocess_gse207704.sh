#!/usr/bin/env bash
# Re-quantify GSE207704 (T47D / MCF7 CLDN4-/- vs WT) from raw reads.
#
# Rationale: the GEO supplementary file only contains cuffdiff FPKM collapsed
# to ONE column per genotype, so the 2 biological replicates per arm cannot be
# separated and no replicate-level statistics are possible. Re-quantifying the
# 8 SRA runs with salmon restores per-sample values.
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
WD=/tmp/gse207704
REF=$WD/ref
mkdir -p "$WD/fastq" "$REF" "$WD/quant"
export MAMBA_ROOT_PREFIX=/tmp/mamba
MM="/tmp/bin/micromamba run -n rnaseq"

TX=$REF/gencode.v44.transcripts.fa.gz
if [ ! -s "$TX" ]; then
  echo "[ref] downloading GENCODE v44 transcripts"
  curl -fsSL --retry 4 -o "$TX" \
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.transcripts.fa.gz || exit 1
fi

if [ ! -d "$REF/salmon_idx" ]; then
  echo "[ref] building salmon index (k=23; reads are 51 bp)"
  $MM salmon index -t "$TX" -i "$REF/salmon_idx" -k 23 --gencode -p 4 \
      > "$ROOT/results_logs_index.txt" 2>&1 || exit 1
fi

runs="SRR20029125 SRR20029124 SRR20029123 SRR20029122 SRR20029121 SRR20029120 SRR20029119 SRR20029118"
for r in $runs; do
  fq=$WD/fastq/$r.fastq.gz
  if [ ! -s "$fq" ]; then
    d1=${r: -2}; d2=${r: -1}
    for url in \
      "https://ftp.sra.ebi.ac.uk/vol1/fastq/${r:0:6}/0${d2}/$r/$r.fastq.gz" \
      "https://ftp.sra.ebi.ac.uk/vol1/fastq/${r:0:6}/${d1}/$r/$r.fastq.gz" \
      "https://ftp.sra.ebi.ac.uk/vol1/fastq/${r:0:6}/$r/$r.fastq.gz" ; do
      echo "[dl] $r <- $url"
      curl -fsSL --retry 3 --max-time 3600 -o "$fq" "$url" && break
      rm -f "$fq"
    done
  fi
  if [ ! -s "$fq" ]; then
    echo "[dl] ENA failed for $r, falling back to fasterq-dump"
    ( cd "$WD/fastq" && $MM fasterq-dump --split-spot -e 4 -O . "$r" \
        && $MM pigz -p 4 "$r.fastq" ) || { echo "[FAIL] $r"; continue; }
  fi
  if [ ! -s "$WD/quant/$r/quant.sf" ]; then
    echo "[quant] $r"
    $MM salmon quant -i "$REF/salmon_idx" -l A -r "$fq" -p 4 \
        --validateMappings --seqBias --gcBias -o "$WD/quant/$r" \
        > "$WD/quant/$r.log" 2>&1 || echo "[FAIL quant] $r"
  fi
  rm -f "$fq"          # keep disk small; each fastq is ~1.5 GB
done

mkdir -p "$ROOT/data/GSE207704_salmon"
for r in $runs; do
  [ -s "$WD/quant/$r/quant.sf" ] && cp "$WD/quant/$r/quant.sf" \
      "$ROOT/data/GSE207704_salmon/$r.quant.sf"
  [ -s "$WD/quant/$r/aux_info/meta_info.json" ] && cp "$WD/quant/$r/aux_info/meta_info.json" \
      "$ROOT/data/GSE207704_salmon/$r.meta_info.json"
done
echo "[done] $(ls "$ROOT/data/GSE207704_salmon" | wc -l) files"
