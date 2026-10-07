#!/usr/bin/env bash
# Stream HD4246_SG_mets_counts.txt (genes x cells, tab-delimited, ~4.96 GB) via
# parallel HTTP range requests. Nothing but the panel gene rows and per-cell
# column sums (library sizes) is ever written to disk.
#
# Each worker owns the lines whose first byte lies in [start, end). It reads
# OVERLAP extra bytes so the last owned line is complete, drops the leading
# partial line, and stops at the first line starting at or after `end`.
set -euo pipefail

URL=${URL:-"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/849/E-MTAB-16849/Files/HD4246_SG_mets_counts.txt"}
SIZE=${SIZE:-4960521535}
NW=${NW:-16}
OVERLAP=${OVERLAP:-4194304}
OUT=${OUT:-/tmp/ae/stream}
GENES="B2M,HLA-A,HLA-B,HLA-C,HLA-E,HLA-F,TAP1,TAP2,PSMB8,PSMB9,NLRC5,TAPBP,STAT1,IRF1,CXCL9,CXCL10,CXCL11,GBP1,GBP2,GBP4,IRF7,ISG15,MX1,OAS1,IFIT1,IFIT3,IDO1,STING1,TMEM173,CGAS,MB21D1,IRF3,CLDN4,TACSTD2"

mkdir -p "$OUT"
CHUNK=$(( (SIZE + NW - 1) / NW ))

{ curl -s --retry 3 -r 0-2097151 "$URL" || true; } | head -n 1 > "$OUT/header.tsv"

worker() {
  local i=$1
  local s=$(( i * CHUNK ))
  # Start one byte early so a line beginning exactly at s is not discarded as the leading fragment.
  (( i > 0 )) && s=$(( s - 1 ))
  local e=$(( (i + 1) * CHUNK )); (( e > SIZE )) && e=$SIZE
  local r=$(( e + OVERLAP - 1 )); (( r > SIZE - 1 )) && r=$(( SIZE - 1 ))
  { curl -s --retry 3 -r "$s-$r" "$URL" || true; } | mawk -F'\t' \
    -v start="$s" -v endb="$e" -v wantstr="$GENES" \
    -v genef="$OUT/genes_$i.tsv" -v sumf="$OUT/colsum_$i.tsv" -v logf="$OUT/log_$i.txt" '
    BEGIN { n = split(wantstr, a, ","); for (k = 1; k <= n; k++) want[a[k]] = 1; off = start }
    {
      len = length($0) + 1
      if (NR == 1) { off += len; next }
      if (off >= endb) exit
      nl++
      if (nc == 0) nc = NF
      if (NF != nc) bad++
      if ($1 in want) print $0 > genef
      for (j = 2; j <= NF; j++) cs[j] += $j
      off += len
    }
    END {
      for (j = 2; j <= nc; j++) printf "%d\t%.0f\n", j - 1, cs[j] > sumf
      printf "worker_start=%d end=%d lines=%d ncol=%d bad_nf=%d\n", start, endb, nl, nc, bad + 0 > logf
    }'
  echo "worker $i done"
}
export -f worker
export URL SIZE CHUNK OVERLAP OUT GENES

seq 0 $(( NW - 1 )) | xargs -P "$NW" -I{} bash -c 'worker {}'
echo "ALL_DONE"
