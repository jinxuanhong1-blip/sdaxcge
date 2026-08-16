#!/usr/bin/env python3
"""Honest TACSTD2 / CLDN4 probe of leftover E-MTAB-15784 (TIL–organoid, not ICI).

Streams each MTX and extracts only a small gene panel. Does not claim ICI
response associations: this experiment is autologous TIL cytotoxicity against
patient-derived NSCLC organoids, with no PD-1/PD-L1 treatment arm.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path

TARGETS = [
    "TACSTD2",
    "CLDN4",
    "CLDN3",
    "CLDN7",
    "TJP1",
    "OCLN",
    "EPCAM",
    "KRT19",
    "KRT7",
    "CD3D",
    "CD8A",
    "CD4",
    "PDCD1",
    "CD274",
    "CD47",
]


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open()


def load_features(path: Path) -> list[str]:
    names = []
    with open_text(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            # 10x-style: id, name[, type]; some files are name only
            if len(parts) >= 2 and parts[1] and parts[1] != parts[0]:
                names.append(parts[1])
            else:
                names.append(parts[0])
    return names


def count_barcodes(path: Path) -> int:
    with open_text(path) as handle:
        return sum(1 for _ in handle)


def stream_gene_sums(mtx_path: Path, gene_rows: dict[str, int], n_cells: int) -> dict[str, dict]:
    """gene_rows maps symbol -> 1-based MTX row index."""
    wanted = {idx: sym for sym, idx in gene_rows.items()}
    umi = {sym: 0.0 for sym in gene_rows}
    pos = {sym: 0 for sym in gene_rows}
    nnz = 0
    with gzip.open(mtx_path, "rt") as handle:
        # skip comments / banner
        line = handle.readline()
        while line.startswith("%"):
            line = handle.readline()
        # header: rows cols nnz
        _nrows, ncols, _nnz = [int(x) for x in line.split()]
        if n_cells and ncols != n_cells:
            # still proceed; some writers swap orientation
            pass
        for line in handle:
            row_s, _col_s, val_s = line.split()
            row = int(row_s)
            if row in wanted:
                val = float(val_s)
                if val > 0:
                    sym = wanted[row]
                    umi[sym] += val
                    pos[sym] += 1
                    nnz += 1
    out = {}
    n = n_cells or ncols
    for sym in gene_rows:
        out[sym] = {
            "n_positive_cells": pos[sym],
            "detection_rate": (pos[sym] / n) if n else None,
            "mean_umi_all_cells": (umi[sym] / n) if n else None,
            "mean_umi_positive": (umi[sym] / pos[sym]) if pos[sym] else 0.0,
            "total_umi": umi[sym],
        }
    out["_n_cells"] = n
    out["_target_nnz"] = nnz
    return out


def parse_sdrf(path: Path) -> dict[str, dict]:
    # SDRF repeats "Derived Array Data File"; csv.DictReader would keep only the last.
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")

    def col(name: str) -> int:
        return header.index(name)

    by_prefix = {}
    for line in lines[1:]:
        if not line.strip():
            continue
        cols = line.split("\t")
        matrix = next(c for c in cols if c.endswith("_matrix.mtx.gz"))
        prefix = matrix.replace("_matrix.mtx.gz", "")
        source = cols[col("Source Name")]
        rec = {
            "source": source,
            "patient": prefix.split("_")[0],
            "disease": cols[col("Characteristics[disease]")],
            "sample_type": cols[col("Characteristics[sample type]")],
            "stage": cols[col("Characteristics[developmental stage]")],
            "sex": cols[col("Characteristics[sex]")],
            "sdrf_source": source,
        }
        # LCP90_Tumor_Organoid SDRF points at LCP89_T_ORG_* files (depositor
        # copy-paste). Keep the first matching matrix owner; do not overwrite.
        if prefix in by_prefix:
            rec["note"] = (
                f"SDRF source {source} reuses processed files of {by_prefix[prefix]['sdrf_source']}"
            )
            by_prefix[prefix].setdefault("sdrf_collisions", []).append(source)
            continue
        by_prefix[prefix] = rec
    return by_prefix


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("results/w200/AE_leftover/downloads/E-MTAB-15784"),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("results/w200/AE_leftover"),
    )
    args = parser.parse_args()

    sdrf = parse_sdrf(args.data_dir / "E-MTAB-15784.sdrf.txt")
    matrices = sorted(args.data_dir.glob("*_matrix.mtx.gz"))
    rows_out = []
    missing_genes = defaultdict(list)

    for mtx in matrices:
        prefix = mtx.name.replace("_matrix.mtx.gz", "")
        feat = args.data_dir / f"{prefix}_features.tsv.gz"
        bc = args.data_dir / f"{prefix}_barcodes.tsv.gz"
        if not feat.exists() or not bc.exists():
            raise FileNotFoundError(f"missing features/barcodes for {prefix}")
        names = load_features(feat)
        n_cells = count_barcodes(bc)
        name_to_row = {}
        for i, name in enumerate(names, start=1):
            # keep first occurrence
            name_to_row.setdefault(name, i)
        gene_rows = {}
        for sym in TARGETS:
            if sym in name_to_row:
                gene_rows[sym] = name_to_row[sym]
            else:
                missing_genes[prefix].append(sym)
        print(f"{prefix}: {n_cells} cells, {len(names)} genes, targets={list(gene_rows)}", flush=True)
        stats = stream_gene_sums(mtx, gene_rows, n_cells)
        meta = sdrf.get(prefix, {})
        for sym in TARGETS:
            rec = {
                "accession": "E-MTAB-15784",
                "sample_prefix": prefix,
                "source": meta.get("source"),
                "patient": meta.get("patient"),
                "disease": meta.get("disease"),
                "sample_type": meta.get("sample_type"),
                "stage": meta.get("stage"),
                "sex": meta.get("sex"),
                "n_cells": n_cells,
                "n_genes": len(names),
                "gene": sym,
                "present_in_features": sym in gene_rows,
            }
            if sym in gene_rows:
                rec.update(stats[sym])
            else:
                rec.update(
                    {
                        "n_positive_cells": None,
                        "detection_rate": None,
                        "mean_umi_all_cells": None,
                        "mean_umi_positive": None,
                        "total_umi": None,
                    }
                )
            rows_out.append(rec)

    tsv = args.out_dir / "E-MTAB-15784_tacstd2_cldn4.tsv"
    args.out_dir.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows_out[0].keys())
    with tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows_out)

    # compact per-sample TACSTD2/CLDN4 view
    key_genes = ["TACSTD2", "CLDN4", "EPCAM", "CD3D", "CD8A"]
    compact = [r for r in rows_out if r["gene"] in key_genes]
    compact_path = args.out_dir / "E-MTAB-15784_key_genes.tsv"
    with compact_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(compact)

    # honest summary numbers
    def mean(xs):
        xs = [x for x in xs if x is not None]
        return sum(xs) / len(xs) if xs else None

    by_type = defaultdict(lambda: defaultdict(list))
    for r in rows_out:
        if r["present_in_features"] and r["detection_rate"] is not None:
            by_type[r["sample_type"]][r["gene"]].append(r["detection_rate"])

    summary = {
        "accession": "E-MTAB-15784",
        "n_matrices": len(matrices),
        "n_rows": len(rows_out),
        "missing_genes_by_sample": dict(missing_genes),
        "mean_detection_rate_by_sample_type": {
            st: {g: mean(vals) for g, vals in genes.items()} for st, genes in by_type.items()
        },
        "caveats": [
            "Not an ICI experiment. No PD-1/PD-L1/CTLA-4 treatment arm or response labels.",
            "Sample types are Tumor_Organoid, Normal Organoid, TILs, and Co-Culture.",
            "Counts are raw UMI from deposited MTX; no cell-type annotation file was deposited.",
            "Detection in TIL libraries is expected to be near-zero for epithelial genes (TACSTD2/CLDN4).",
            "Co-culture mixes organoid epithelium and TILs, so bulk-per-cell means are diluted.",
        ],
    }
    (args.out_dir / "E-MTAB-15784_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n"
    )
    print(json.dumps(summary["mean_detection_rate_by_sample_type"], indent=2))
    print(f"Wrote {tsv}")


if __name__ == "__main__":
    main()
