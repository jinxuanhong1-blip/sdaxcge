#!/usr/bin/env bash
# Processed matrices for the noskip blood/PBMC ICI catalog (each file < 2 GB).
set -euo pipefail
GEO_DIR="${GEO_DIR:-/tmp/geo}"
NS="$GEO_DIR/noskip"
mkdir -p "$NS"
cd "$NS"
base="https://ftp.ncbi.nlm.nih.gov/geo"

fetch() { local dest="$1" url="$2"; [[ -s "$dest" ]] && echo "have $dest" && return; echo "GET $dest"; curl -fsS -o "$dest" "$url"; }

fetch GSE111414_gene_counts.csv.gz "$base/series/GSE111nnn/GSE111414/suppl/GSE111414_gene_counts.csv.gz"
fetch GSE111414_sm.txt.gz          "$base/series/GSE111nnn/GSE111414/matrix/GSE111414_series_matrix.txt.gz"
fetch GSE202417_sm.txt.gz          "$base/series/GSE202nnn/GSE202417/matrix/GSE202417_series_matrix.txt.gz"
fetch GSE141479_sm.txt.gz          "$base/series/GSE141nnn/GSE141479/matrix/GSE141479_series_matrix.txt.gz"
fetch GSE235048_24CFlow_TPM.txt.gz "$base/series/GSE235nnn/GSE235048/suppl/GSE235048_24CFlow_TPM.txt.gz"
fetch GSE235048_sm.txt.gz          "$base/series/GSE235nnn/GSE235048/matrix/GSE235048_series_matrix.txt.gz"
fetch GSE216297_TEP_Count_Matrix.RData.gz "$base/series/GSE216nnn/GSE216297/suppl/GSE216297_TEP_Count_Matrix.RData.gz"
fetch GSE216297_sm.txt.gz          "$base/series/GSE216nnn/GSE216297/matrix/GSE216297_series_matrix.txt.gz"
fetch GSE310370_sm.txt.gz          "$base/series/GSE310nnn/GSE310370/matrix/GSE310370_series_matrix.txt.gz"
fetch GSE100860_RAW.tar            "$base/series/GSE100nnn/GSE100860/suppl/GSE100860_RAW.tar"
fetch GSE100860_sm.txt.gz          "$base/series/GSE100nnn/GSE100860/matrix/GSE100860_series_matrix.txt.gz"
fetch GSE213902_MQpool1_features.tsv.gz "$base/series/GSE213nnn/GSE213902/suppl/GSE213902_MQpool1_features.tsv.gz"
fetch GSE213902_MQpool2_features.tsv.gz "$base/series/GSE213nnn/GSE213902/suppl/GSE213902_MQpool2_features.tsv.gz"
fetch GSE213902_MQpool1_matrix.mtx.gz   "$base/series/GSE213nnn/GSE213902/suppl/GSE213902_MQpool1_matrix.mtx.gz"
fetch GSE213902_MQpool2_matrix.mtx.gz   "$base/series/GSE213nnn/GSE213902/suppl/GSE213902_MQpool2_matrix.mtx.gz"
fetch GSE213902_sm.txt.gz          "$base/series/GSE213nnn/GSE213902/matrix/GSE213902_series_matrix.txt.gz"
echo "done"
