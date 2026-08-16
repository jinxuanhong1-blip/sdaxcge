#!/usr/bin/env bash
# End-to-end reproduction of the fable_tisch slice.
#   bash scripts/fable_tisch/run_all.sh [DATA_DIR]
set -euo pipefail
cd "$(dirname "$0")/../.."
DATA_DIR="${1:-data_fable_tisch}"

bash scripts/fable_tisch/01_download_data.sh "$DATA_DIR"

for ds in NSCLC_GSE131907 NSCLC_EMTAB6149; do
  python scripts/fable_tisch/02_trop2_cldn4_tumor_vs_immune.py \
    --dataset "$ds" --data-dir "$DATA_DIR" --out-root results/fable_tisch
done

python scripts/fable_tisch/03_cross_dataset_summary.py --out-root results/fable_tisch
