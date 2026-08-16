#!/usr/bin/env bash
# Download processed expression matrices / DE tables for each series.
set -u
cd "$(dirname "$0")/../data" || exit 1
B=https://ftp.ncbi.nlm.nih.gov/geo

get () { # url outfile
  [ -s "$2" ] && { echo "have $2"; return; }
  for i in 1 2 3 4; do
    curl -fsSL --retry 3 --max-time 900 -o "$2" "$1" && { echo "ok  $2"; return; }
    echo "retry $i $1"; sleep $((4 ** i))
  done
  echo "FAIL $1"
}

get $B/series/GSE207nnn/GSE207704/suppl/GSE207704_CLDN4_RNAseq.txt.gz     GSE207704_CLDN4_RNAseq.txt.gz
get $B/series/GSE22nnn/GSE22493/suppl/GSE22493_RAW.tar                    GSE22493_RAW.tar
get $B/series/GSE50nnn/GSE50927/suppl/GSE50927_Cldn4lungWTvsKOgenes.csv.gz GSE50927_Cldn4lungWTvsKOgenes.csv.gz
get $B/series/GSE50nnn/GSE50927/suppl/GSE50927_VILIwtkohiGenes.csv.gz     GSE50927_VILIwtkohiGenes.csv.gz
get $B/series/GSE50nnn/GSE50927/suppl/GSE50927_VILIwtkoloGenes.csv.gz     GSE50927_VILIwtkoloGenes.csv.gz
get $B/series/GSE274nnn/GSE274940/suppl/GSE274940_raw_counts.csv.gz       GSE274940_raw_counts.csv.gz
get $B/series/GSE334nnn/GSE334497/suppl/GSE334497_normalized_counts.csv.gz GSE334497_normalized_counts.csv.gz
get $B/series/GSE289nnn/GSE289287/suppl/GSE289287_DESeq2-Trop2KO_tumors_vs_WT.tsv.gz GSE289287_Trop2KO_tumors_vs_WT.tsv.gz
get $B/series/GSE245nnn/GSE245459/suppl/GSE245459_fpkm.anno.txt.gz        GSE245459_fpkm.anno.txt.gz

mkdir -p GSE22493_RAW && tar -xf GSE22493_RAW.tar -C GSE22493_RAW 2>/dev/null
ls -la
