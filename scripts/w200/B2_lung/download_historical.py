#!/usr/bin/env python3
"""Optional historical CCLE extracts (22Q2 CCLE_expression + 2018 RPKM GCT).

These are named sensitivities only. The primary B2 number uses DepMap 24Q4.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import sys
import urllib.request
from pathlib import Path

UA = "sdaxcge-w200-B2-lung/1.0"
CCLE22_EXPR = "https://ndownloader.figshare.com/files/34989919"
CCLE22_SI = "https://ndownloader.figshare.com/files/35020903"
CCLE2018_GCT = "https://data.broadinstitute.org/ccle/CCLE_RNAseq_genes_rpkm_20180929.gct.gz"
CCLE2018_ANN = "https://data.broadinstitute.org/ccle/Cell_lines_annotations_20181226.txt"


def get(url: str):
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    return urllib.request.urlopen(req, timeout=600)


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"GET {url} -> {dest}", flush=True)
    with get(url) as resp, dest.open("wb") as out:
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)


def extract_22q2(out_csv: Path) -> None:
    print("streaming 22Q2 CCLE_expression.csv", flush=True)
    with get(CCLE22_EXPR) as resp:
        reader = csv.reader(io.TextIOWrapper(resp, encoding="utf-8", newline=""))
        header = next(reader)
        tac, cld = "TACSTD2 (4070)", "CLDN4 (1364)"
        if tac not in header or cld not in header:
            raise SystemExit(f"22Q2 missing gene columns: {header[:5]}")
        i_t, i_c = header.index(tac), header.index(cld)
        with out_csv.open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["ModelID", "TACSTD2", "CLDN4", "TACSTD2_column", "CLDN4_column"])
            n = 0
            for row in reader:
                if not row:
                    continue
                w.writerow([row[0], row[i_t], row[i_c], tac, cld])
                n += 1
    print(f"22Q2 extract n={n}", flush=True)


def extract_2018(out_csv: Path, gct_cache: Path) -> None:
    if not gct_cache.exists() or gct_cache.stat().st_size < 1000:
        download(CCLE2018_GCT, gct_cache)
    need = {"TACSTD2": None, "CLDN4": None}
    with gzip.open(gct_cache, "rt") as tf:
        tf.readline()
        tf.readline()
        header = tf.readline()
        samples = header.rstrip("\n").split("\t")[2:]
        for line in tf:
            parts = line.rstrip("\n").split("\t")
            desc = parts[1]
            if desc in need and need[desc] is None:
                need[desc] = parts[2:]
                print("found", desc, parts[0], flush=True)
            if all(need.values()):
                break
    if any(v is None for v in need.values()):
        raise SystemExit(f"2018 GCT missing genes: {need}")
    with out_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["CCLE_ID", "TACSTD2_RPKM", "CLDN4_RPKM"])
        for i, s in enumerate(samples):
            w.writerow([s, need["TACSTD2"][i], need["CLDN4"][i]])
    print(f"2018 extract n={len(samples)}", flush=True)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--cache-dir", default="/tmp/ccle2018")
    p.add_argument("--out-dir", default="results/w200/B2_lung")
    args = p.parse_args()
    cache = Path(args.cache_dir)
    out = Path(args.out_dir)
    cache.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)

    out22 = out / "ccle22q2_tacstd2_cldn4_all_models.csv"
    if not out22.exists():
        extract_22q2(out22)
    si = cache / "sample_info_22q2.csv"
    if not si.exists():
        download(CCLE22_SI, si)
    out18 = out / "ccle2018_rpkm_tacstd2_cldn4.csv"
    if not out18.exists():
        extract_2018(out18, cache / "CCLE_RNAseq_genes_rpkm_20180929.gct.gz")
    ann = cache / "Cell_lines_annotations_20181226.txt"
    if not ann.exists():
        download(CCLE2018_ANN, ann)
    print("historical extracts ready")
    return 0


if __name__ == "__main__":
    sys.exit(main())
