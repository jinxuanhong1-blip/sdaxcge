#!/usr/bin/env bash
# Quantify GSE207704 (CLDN4 KO T47D/MCF7) with kallisto.
# SRA objects: sra-pub-run-odp S3 copies of SRR20029118–SRR20029125.
# sratoolkit 3.4.1 fasterq-dump segfaults on these files here; fastq-dump is
# run in 2 million spot chunks and concatenated.
# Transcriptome: Ensembl 116 GRCh38 cDNA.
set -euo pipefail

ROOT="${1:-/tmp/gse207704}"
KALLISTO="${KALLISTO:-kallisto}"
THREADS="${THREADS:-4}"
mkdir -p "$ROOT"/{ref,fastq,quant,logs}
cd "$ROOT"

if [[ ! -s ref/cdna.fa.gz ]]; then
  curl -fL --retry 4 -o ref/cdna.fa.gz \
    "https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
fi

if [[ ! -s ref/hs_cdna.idx ]]; then
  "$KALLISTO" index -i ref/hs_cdna.idx ref/cdna.fa.gz | tee logs/index.log
fi

# ENA FTP timed out from this environment. Full SRA objects are on NCBI S3.
mkdir -p sra chunks
FASTQ_DUMP="${FASTQ_DUMP:-fastq-dump}"
declare -A SPOTS=(
  [SRR20029118]=44087640
  [SRR20029119]=47026612
  [SRR20029120]=45458325
  [SRR20029121]=46699152
  [SRR20029122]=44387540
  [SRR20029123]=44424212
  [SRR20029124]=45214824
  [SRR20029125]=44062012
)
CHUNK=2000000
for srr in "${!SPOTS[@]}"; do
  if [[ -s fastq/${srr}.fastq.gz ]]; then
    continue
  fi
  if [[ ! -s sra/${srr} ]]; then
    curl -fL --retry 4 -o "sra/${srr}.partial" \
      "https://sra-pub-run-odp.s3.amazonaws.com/sra/${srr}/${srr}"
    mv "sra/${srr}.partial" "sra/${srr}"
  fi
  total=${SPOTS[$srr]}
  start=1
  rm -f "fastq/${srr}.fastq.gz"
  while (( start <= total )); do
    end=$(( start + CHUNK - 1 ))
    if (( end > total )); then end=$total; fi
    odir="chunks/${srr}/${start}"
    mkdir -p "$odir"
    "$FASTQ_DUMP" --gzip -N "$start" -X "$end" -O "$odir" "sra/${srr}"
    gzip -t "${odir}/${srr}.fastq.gz"
    cat "${odir}/${srr}.fastq.gz" >> "fastq/${srr}.fastq.gz"
    start=$(( end + 1 ))
  done
  gzip -t "fastq/${srr}.fastq.gz"
done

# Fragment length is not in the GEO record. -l 200 -s 30 is the kallisto
# single-end default prior. Gene-level estimated counts, not TPM, are the
# differential-expression input, so this prior does not set the log2FC sign.
for srr in SRR20029125 SRR20029124 SRR20029123 SRR20029122 SRR20029121 SRR20029120 SRR20029119 SRR20029118; do
  if [[ ! -s quant/${srr}/abundance.tsv ]]; then
    "$KALLISTO" quant -i ref/hs_cdna.idx -o "quant/${srr}" \
      --single -l 200 -s 30 -t "$THREADS" "fastq/${srr}.fastq.gz" \
      | tee "logs/${srr}.quant.log"
  fi
done

echo "quantification complete: $ROOT/quant"
