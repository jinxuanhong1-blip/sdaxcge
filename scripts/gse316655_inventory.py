#!/usr/bin/env python3
"""Inventory GSE316655 and score CLDN4 in the deposited 10x matrices.

The GEO series is anti-LILRB2 versus isotype scRNA-seq of FACS-sorted human
CD45+ cells (Liu et al., Sci Immunol 2026, eadt7832). This script does not
test immune exclusion. It records file sizes and how many cells carry CLDN4
or CLDN18 UMIs.
"""

from __future__ import annotations

import argparse
import gzip
import urllib.request
from collections import defaultdict
from pathlib import Path

FILELIST_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE316nnn/GSE316655/suppl/filelist.txt"
)

SAMPLES = [
    {
        "gsm": "GSM9457798",
        "library": "ND_anti_B2",
        "title": "ND tumor LILRB2 treatment",
        "tissue": "tumor",
        "treatment": "anti-LILRB2",
        "line_in_protocol": "SK-MEL-5",
        "srx": "SRX31826394",
        "biosample": "SAMN54677775",
    },
    {
        "gsm": "GSM9457799",
        "library": "ND_CTR",
        "title": "ND tumor Isotype",
        "tissue": "tumor",
        "treatment": "isotype",
        "line_in_protocol": "SK-MEL-5",
        "srx": "SRX31826395",
        "biosample": "SAMN54677774",
    },
    {
        "gsm": "GSM9457800",
        "library": "PB_anti_B2",
        "title": "PB LILRB2 treatment",
        "tissue": "peripheral blood",
        "treatment": "anti-LILRB2",
        "line_in_protocol": "SK-MEL-5",
        "srx": "SRX31826396",
        "biosample": "SAMN54677773",
    },
    {
        "gsm": "GSM9457801",
        "library": "PB_CTR",
        "title": "PB Isotype",
        "tissue": "peripheral blood",
        "treatment": "isotype",
        "line_in_protocol": "SK-MEL-5",
        "srx": "SRX31826397",
        "biosample": "SAMN54677772",
    },
]

HUMAN = [
    "CLDN4",
    "CLDN18",
    "CLDN3",
    "CLDN7",
    "LILRB2",
    "LILRB5",
    "CD8A",
    "NKG7",
    "GZMB",
    "PRF1",
    "IFNG",
    "CD274",
    "IDO1",
    "NEAT1",
    "PTPRC",
    "KRT8",
    "EPCAM",
]
MOUSE = ["Cldn4", "Cldn18", "Ptprc"]


def strip_symbol(symbol: str) -> tuple[str, str]:
    if symbol.startswith("GRCh38_"):
        return symbol[len("GRCh38_") :], "GRCh38"
    if symbol.startswith("mm10___"):
        return symbol[len("mm10___") :], "mm10"
    if symbol.startswith("mm10_"):
        return symbol[len("mm10_") :], "mm10"
    return symbol, "other"


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise SystemExit(f"no rows for {path}")
    cols = list(rows[0].keys())
    lines = ["\t".join(cols)]
    for row in rows:
        lines.append("\t".join(str(row[c]) for c in cols))
    path.write_text("\n".join(lines) + "\n")


def fetch_filelist(dest: Path) -> str:
    if dest.exists() and dest.stat().st_size > 0:
        return dest.read_text()
    req = urllib.request.Request(FILELIST_URL, headers={"User-Agent": "Mozilla/5.0"})
    text = urllib.request.urlopen(req, timeout=60).read().decode()
    dest.write_text(text)
    return text


def parse_filelist(text: str) -> list[dict]:
    rows = []
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        kind, name, time, size, typ = line.split("\t")
        rows.append(
            {
                "kind": kind,
                "name": name,
                "geo_time": time,
                "size_bytes": size,
                "type": typ,
            }
        )
    return rows


def score_sample(mtx_dir: Path, sample: dict) -> list[dict]:
    gsm = sample["gsm"]
    feat_path = mtx_dir / f"{gsm}_features.tsv.gz"
    bar_path = mtx_dir / f"{gsm}_barcodes.tsv.gz"
    mtx_path = mtx_dir / f"{gsm}_matrix.mtx.gz"
    genes = []
    with gzip.open(feat_path, "rt") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            symbol, genome = strip_symbol(parts[1] if len(parts) > 1 else parts[0])
            genes.append((symbol, genome))
    n_barcodes = sum(1 for _ in gzip.open(bar_path, "rt"))
    wanted = {}
    for index, (symbol, genome) in enumerate(genes, start=1):
        if genome == "GRCh38" and symbol in HUMAN:
            wanted[index] = symbol
        elif genome == "mm10" and symbol in MOUSE:
            wanted[index] = f"mm10:{symbol}"
    sums = defaultdict(float)
    nz = defaultdict(int)
    with gzip.open(mtx_path, "rt") as handle:
        line = handle.readline()
        while line.startswith("%"):
            line = handle.readline()
        n_genes, n_cells, nnz = map(int, line.split())
        for line in handle:
            row_s, _col_s, val_s = line.split()
            name = wanted.get(int(row_s))
            if name is None:
                continue
            val = float(val_s)
            sums[name] += val
            if val > 0:
                nz[name] += 1
    if n_barcodes != n_cells:
        raise SystemExit(f"{gsm}: barcodes {n_barcodes} != matrix columns {n_cells}")
    present_names = set(wanted.values())
    rows = []
    for symbol in HUMAN + [f"mm10:{s}" for s in MOUSE]:
        rows.append(
            {
                "gsm": gsm,
                "library": sample["library"],
                "tissue": sample["tissue"],
                "treatment": sample["treatment"],
                "n_cells": n_cells,
                "n_genes": n_genes,
                "nnz": nnz,
                "gene": symbol,
                "on_reference": "yes" if symbol in present_names else "no",
                "nonzero_cells": nz[symbol],
                "umi_sum": int(sums[symbol]),
                "pct_cells": f"{(100.0 * nz[symbol] / n_cells):.4f}",
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/gse316655_cldn_lilrb"),
    )
    parser.add_argument(
        "--mtx-dir",
        type=Path,
        default=None,
        help="Directory with {GSM}_{features,barcodes,matrix} files already downloaded",
    )
    parser.add_argument(
        "--filelist",
        type=Path,
        default=Path("/tmp/gse316655/filelist.txt"),
    )
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    text = fetch_filelist(args.filelist)
    write_tsv(args.out / "file_inventory.tsv", parse_filelist(text))
    write_tsv(args.out / "sample_design.tsv", SAMPLES)
    if args.mtx_dir is None:
        print("file inventory written; pass --mtx-dir to score UMIs")
        return
    scored = []
    for sample in SAMPLES:
        scored.extend(score_sample(args.mtx_dir, sample))
    write_tsv(args.out / "gene_detection.tsv", scored)
    print(f"wrote {len(scored)} gene rows to {args.out / 'gene_detection.tsv'}")


if __name__ == "__main__":
    main()
