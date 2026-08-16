#!/usr/bin/env bash
# =============================================================================
# 04 - pySCENIC: infer transcription-factor regulons in lung epithelium
#
# Question this feeds: which TFs define the TACSTD2/CLDN4-high epithelial state,
# and are the T-cell-recruitment chemokine genes (CXCL9/10/11, CCL5, CXCL16)
# INSIDE any active regulon or anticorrelated with the state's regulons?
#
# What SCENIC does / does not do. Pipeline = GRNBoost2 (TF-target co-expression)
# -> cisTarget (prune to targets with the TF motif) -> AUCell (per-cell regulon
# activity). It finds regulons whose activity CORRELATES with a cell state; it
# does NOT prove a TF regulates (or represses) a gene. Only TFs with motifs in
# the (genome-specific) cisTarget databases get regulons; repression is not
# modeled. GRNBoost2 is stochastic -> run multiple seeds. See playbook §SCENIC.
#
# PRACTICALITY on 10x lung ICI: cisTarget databases are ~GBs and GRNBoost2 is
# heavy on 10^5 cells. Subsample epithelium for GRN inference, then AUCell-score
# ALL cells. Budget CPU/RAM accordingly.
# =============================================================================
set -euo pipefail

# ---- inputs / resources (EDIT) ----------------------------------------------
LOOM_IN="data/epithelium_for_scenic.loom"     # from 04b_scenic_downstream.py step 0
OUTDIR="results/scenic"; mkdir -p "$OUTDIR"
N_WORKERS="${N_WORKERS:-8}"
SEED="${SEED:-777}"

# Resources for GRCh38 (download once; pin versions):
#   - TF list:      allTFs_hg38.txt                (aertslab)
#   - motif ranks:  hg38_*.genes_vs_motifs.rankings.feather  (cisTarget DBs)
#   - motif2tf:     motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
TF_LIST="${TF_LIST:-resources/allTFs_hg38.txt}"
RANKINGS="${RANKINGS:-resources/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather}"
MOTIF2TF="${MOTIF2TF:-resources/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl}"

# ---- 1. GRN inference (stochastic; run a few seeds & take consensus) --------
pyscenic grn \
  --num_workers "$N_WORKERS" \
  --seed "$SEED" \
  --method grnboost2 \
  -o "$OUTDIR/adj.seed${SEED}.tsv" \
  "$LOOM_IN" "$TF_LIST"

# ---- 2. cisTarget: prune to motif-supported regulons ------------------------
pyscenic ctx \
  "$OUTDIR/adj.seed${SEED}.tsv" \
  "$RANKINGS" \
  --annotations_fname "$MOTIF2TF" \
  --expression_mtx_fname "$LOOM_IN" \
  --mode "dask_multiprocessing" \
  --min_genes 10 \
  --num_workers "$N_WORKERS" \
  --output "$OUTDIR/regulons.seed${SEED}.csv"

# ---- 3. AUCell: score regulon activity per cell -----------------------------
pyscenic aucell \
  "$LOOM_IN" \
  "$OUTDIR/regulons.seed${SEED}.csv" \
  --num_workers "$N_WORKERS" \
  --seed "$SEED" \
  --output "$OUTDIR/aucell.seed${SEED}.loom"

echo "[done] SCENIC seed=${SEED} -> $OUTDIR . Run >=3 seeds, then 04b to analyze."
