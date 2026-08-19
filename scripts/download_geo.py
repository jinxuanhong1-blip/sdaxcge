#!/usr/bin/env python3
"""Download processed Visium matrices + coordinates from GEO (no images)."""

from __future__ import annotations

import argparse
import os
import sys
import time
import urllib.error
import urllib.request

FTP = "https://ftp.ncbi.nlm.nih.gov/geo/samples"


def gsm_bucket(gsm: str) -> str:
    # GSM5702473 -> GSM5702nnn
    return gsm[:7] + "nnn"


def download(url: str, dest: str, retries: int = 5) -> None:
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        print(f"  exists {dest} ({os.path.getsize(dest)} bytes)")
        return
    tmp = dest + ".part"
    delay = 4
    last_err = None
    for attempt in range(1, retries + 1):
        try:
            print(f"  GET {url} -> {dest} (try {attempt})")
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=300) as r, open(tmp, "wb") as f:
                while True:
                    chunk = r.read(1024 * 1024)
                    if not chunk:
                        break
                    f.write(chunk)
            os.replace(tmp, dest)
            print(f"  wrote {dest} ({os.path.getsize(dest)} bytes)")
            return
        except (urllib.error.URLError, TimeoutError, OSError) as e:
            last_err = e
            print(f"  fail: {e}")
            if os.path.exists(tmp):
                os.remove(tmp)
            time.sleep(delay)
            delay *= 2
    raise RuntimeError(f"failed {url}: {last_err}")


def gsm_file(gsm: str, filename: str, dest_dir: str) -> str:
    url = f"{FTP}/{gsm_bucket(gsm)}/{gsm}/suppl/{filename}"
    dest = os.path.join(dest_dir, filename)
    download(url, dest)
    return dest


# Leftover series selected for full run (processed matrix + coordinates, open).
SERIES = {
    "GSE189487": {
        "dir": "data/GSE189487",
        "files": [
            ("GSM5702473", "GSM5702473_TD1_barcodes.tsv.gz"),
            ("GSM5702473", "GSM5702473_TD1_features.tsv.gz"),
            ("GSM5702473", "GSM5702473_TD1_matrix.mtx.gz"),
            ("GSM5702473", "GSM5702473_TD1_tissue_positions_list.csv.gz"),
            ("GSM5702474", "GSM5702474_TD2_barcodes.tsv.gz"),
            ("GSM5702474", "GSM5702474_TD2_features.tsv.gz"),
            ("GSM5702474", "GSM5702474_TD2_matrix.mtx.gz"),
            ("GSM5702474", "GSM5702474_TD2_tissue_positions_list.csv.gz"),
            ("GSM5702475", "GSM5702475_TD3_barcodes.tsv.gz"),
            ("GSM5702475", "GSM5702475_TD3_features.tsv.gz"),
            ("GSM5702475", "GSM5702475_TD3_matrix.mtx.gz"),
            ("GSM5702475", "GSM5702475_TD3_tissue_positions_list.csv.gz"),
            ("GSM5702476", "GSM5702476_TD5_barcodes.tsv.gz"),
            ("GSM5702476", "GSM5702476_TD5_features.tsv.gz"),
            ("GSM5702476", "GSM5702476_TD5_matrix.mtx.gz"),
            ("GSM5702476", "GSM5702476_TD5_tissue_positions_list.csv.gz"),
            ("GSM5702477", "GSM5702477_TD6_barcodes.tsv.gz"),
            ("GSM5702477", "GSM5702477_TD6_features.tsv.gz"),
            ("GSM5702477", "GSM5702477_TD6_matrix.mtx.gz"),
            ("GSM5702477", "GSM5702477_TD6_tissue_positions_list.csv.gz"),
            ("GSM5702478", "GSM5702478_TD8_barcodes.tsv.gz"),
            ("GSM5702478", "GSM5702478_TD8_features.tsv.gz"),
            ("GSM5702478", "GSM5702478_TD8_matrix.mtx.gz"),
            ("GSM5702478", "GSM5702478_TD8_tissue_positions_list.csv.gz"),
        ],
    },
    "GSE322553": {
        "dir": "data/GSE322553",
        "files": [
            # NSCLC-only leftover slides (CRC/HCC skipped at analysis time)
            ("GSM9554209", "GSM9554209_NSCLC_936_1_filtered_feature_bc_matrix.h5.gz"),
            ("GSM9554209", "GSM9554209_NSCLC_936_1_tissue_positions_list.csv.gz"),
            ("GSM9554210", "GSM9554210_NSCLC_41485_2_filtered_feature_bc_matrix.h5.gz"),
            ("GSM9554210", "GSM9554210_NSCLC_41485_2_tissue_positions_list.csv.gz"),
        ],
    },
    "GSE273378": {
        "dir": "data/GSE273378",
        "files": [
            ("GSM8427428", "GSM8427428_LM_SD_1216_1_barcodes.tsv.gz"),
            ("GSM8427428", "GSM8427428_LM_SD_1216_1_features.tsv.gz"),
            ("GSM8427428", "GSM8427428_LM_SD_1216_1_matrix.mtx.gz"),
            ("GSM8427428", "GSM8427428_LM_SD_1216_1_tissue_positions_list.csv.gz"),
            ("GSM8427429", "GSM8427429_LM_SD_16_barcodes.tsv.gz"),
            ("GSM8427429", "GSM8427429_LM_SD_16_features.tsv.gz"),
            ("GSM8427429", "GSM8427429_LM_SD_16_matrix.mtx.gz"),
            ("GSM8427429", "GSM8427429_LM_SD_16_tissue_positions_list.csv.gz"),
            ("GSM8427430", "GSM8427430_LM_SD_11_barcodes.tsv.gz"),
            ("GSM8427430", "GSM8427430_LM_SD_11_features.tsv.gz"),
            ("GSM8427430", "GSM8427430_LM_SD_11_matrix.mtx.gz"),
            ("GSM8427430", "GSM8427430_LM_SD_11_tissue_positions_list.csv.gz"),
            ("GSM8427431", "GSM8427431_LM_SD_2_barcodes.tsv.gz"),
            ("GSM8427431", "GSM8427431_LM_SD_2_features.tsv.gz"),
            ("GSM8427431", "GSM8427431_LM_SD_2_matrix.mtx.gz"),
            ("GSM8427431", "GSM8427431_LM_SD_2_tissue_positions_list.csv.gz"),
            ("GSM8427432", "GSM8427432_LM_SD_3_barcodes.tsv.gz"),
            ("GSM8427432", "GSM8427432_LM_SD_3_features.tsv.gz"),
            ("GSM8427432", "GSM8427432_LM_SD_3_matrix.mtx.gz"),
            ("GSM8427432", "GSM8427432_LM_SD_3_tissue_positions_list.csv.gz"),
            ("GSM8427433", "GSM8427433_LM_SD_4_barcodes.tsv.gz"),
            ("GSM8427433", "GSM8427433_LM_SD_4_features.tsv.gz"),
            ("GSM8427433", "GSM8427433_LM_SD_4_matrix.mtx.gz"),
            ("GSM8427433", "GSM8427433_LM_SD_4_tissue_positions_list.csv.gz"),
            ("GSM8427434", "GSM8427434_LM_SD_5_barcodes.tsv.gz"),
            ("GSM8427434", "GSM8427434_LM_SD_5_features.tsv.gz"),
            ("GSM8427434", "GSM8427434_LM_SD_5_matrix.mtx.gz"),
            ("GSM8427434", "GSM8427434_LM_SD_5_tissue_positions_list.csv.gz"),
            ("GSM8427435", "GSM8427435_LM_SD_6_barcodes.tsv.gz"),
            ("GSM8427435", "GSM8427435_LM_SD_6_features.tsv.gz"),
            ("GSM8427435", "GSM8427435_LM_SD_6_matrix.mtx.gz"),
            ("GSM8427435", "GSM8427435_LM_SD_6_tissue_positions_list.csv.gz"),
            ("GSM8427436", "GSM8427436_LM_SD_7_barcodes.tsv.gz"),
            ("GSM8427436", "GSM8427436_LM_SD_7_features.tsv.gz"),
            ("GSM8427436", "GSM8427436_LM_SD_7_matrix.mtx.gz"),
            ("GSM8427436", "GSM8427436_LM_SD_7_tissue_positions_list.csv.gz"),
            ("GSM8427437", "GSM8427437_LM_SD_1216_8_barcodes.tsv.gz"),
            ("GSM8427437", "GSM8427437_LM_SD_1216_8_features.tsv.gz"),
            ("GSM8427437", "GSM8427437_LM_SD_1216_8_matrix.mtx.gz"),
            ("GSM8427437", "GSM8427437_LM_SD_1216_8_tissue_positions_list.csv.gz"),
            ("GSM8427438", "GSM8427438_LM_SD_9_barcodes.tsv.gz"),
            ("GSM8427438", "GSM8427438_LM_SD_9_features.tsv.gz"),
            ("GSM8427438", "GSM8427438_LM_SD_9_matrix.mtx.gz"),
            ("GSM8427438", "GSM8427438_LM_SD_9_tissue_positions_list.csv.gz"),
            ("GSM8427439", "GSM8427439_LM_SD_10_barcodes.tsv.gz"),
            ("GSM8427439", "GSM8427439_LM_SD_10_features.tsv.gz"),
            ("GSM8427439", "GSM8427439_LM_SD_10_matrix.mtx.gz"),
            ("GSM8427439", "GSM8427439_LM_SD_10_tissue_positions_list.csv.gz"),
            ("GSM8427440", "GSM8427440_LM_SD_1216_12_barcodes.tsv.gz"),
            ("GSM8427440", "GSM8427440_LM_SD_1216_12_features.tsv.gz"),
            ("GSM8427440", "GSM8427440_LM_SD_1216_12_matrix.mtx.gz"),
            ("GSM8427440", "GSM8427440_LM_SD_1216_12_tissue_positions_list.csv.gz"),
            ("GSM8427441", "GSM8427441_LM_SD_13_barcodes.tsv.gz"),
            ("GSM8427441", "GSM8427441_LM_SD_13_features.tsv.gz"),
            ("GSM8427441", "GSM8427441_LM_SD_13_matrix.mtx.gz"),
            ("GSM8427441", "GSM8427441_LM_SD_13_tissue_positions_list.csv.gz"),
            ("GSM8427442", "GSM8427442_LM_SD_1216_14_barcodes.tsv.gz"),
            ("GSM8427442", "GSM8427442_LM_SD_1216_14_features.tsv.gz"),
            ("GSM8427442", "GSM8427442_LM_SD_1216_14_matrix.mtx.gz"),
            ("GSM8427442", "GSM8427442_LM_SD_1216_14_tissue_positions_list.csv.gz"),
            ("GSM8427443", "GSM8427443_LM_SD_15_barcodes.tsv.gz"),
            ("GSM8427443", "GSM8427443_LM_SD_15_features.tsv.gz"),
            ("GSM8427443", "GSM8427443_LM_SD_15_matrix.mtx.gz"),
            ("GSM8427443", "GSM8427443_LM_SD_15_tissue_positions_list.csv.gz"),
        ],
    },
}


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("series", nargs="*", default=list(SERIES))
    args = p.parse_args()
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    for acc in args.series:
        if acc not in SERIES:
            print(f"unknown series {acc}", file=sys.stderr)
            return 2
        spec = SERIES[acc]
        dest_dir = os.path.join(root, spec["dir"])
        print(f"=== {acc} -> {dest_dir} ===")
        for gsm, filename in spec["files"]:
            gsm_file(gsm, filename, dest_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
