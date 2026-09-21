#!/usr/bin/env python3
"""Stream GSE131907 UMI text to a slim genes x cells matrix.

Python only parses the wide text file. CellChat itself runs in R.
Keeps malignant and T/NK cells from locked samples (n_malignant > 0).
Library size is the full-transcriptome UMI sum of those cells.
"""
import argparse
import gzip
import os
import sys

import numpy as np
from scipy import sparse
from scipy.io import mmwrite


MALIG = {"Malignant cells", "tS1", "tS2", "tS3"}
TNK = {"T lymphocytes", "NK cells"}


def read_tsv(path):
    rows = []
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            rows.append(dict(zip(header, parts)))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--umi", required=True)
    ap.add_argument("--annot", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    wanted = {}
    with open(args.genes) as handle:
        for line in handle:
            g = line.strip()
            if g:
                wanted[g.upper()] = g
    print(f"wanted genes {len(wanted)}", flush=True)

    locked = set()
    for row in read_tsv(args.samples):
        try:
            nmal = float(row.get("n_malignant", "0"))
        except ValueError:
            nmal = 0.0
        if nmal > 0:
            locked.add(row["sample"])
    print(f"locked samples {len(locked)}", flush=True)

    annot = {}
    with gzip.open(args.annot, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            rec = dict(zip(header, parts))
            annot[rec["Index"]] = rec

    with gzip.open(args.umi, "rt") as handle:
        umi_header = handle.readline().rstrip("\n").split("\t")
    cell_ids = umi_header[1:]
    col_of = {cid: i for i, cid in enumerate(cell_ids)}

    keep = []
    for cid, i in col_of.items():
        rec = annot.get(cid)
        if rec is None or rec.get("Sample") not in locked:
            continue
        subtype = rec.get("Cell_subtype", "")
        ctype = rec.get("Cell_type", "")
        mal = subtype in MALIG
        tnk = (ctype in TNK) and (not mal)
        if not mal and not tnk:
            continue
        keep.append((cid, i, rec["Sample"], ctype, subtype, mal, tnk))
    keep.sort(key=lambda x: (x[2], x[0]))
    print(f"kept cells {len(keep)}", flush=True)
    if not keep:
        sys.exit("no malignant/TNK cells matched annotation Index to UMI columns")

    src_cols = np.array([k[1] for k in keep], dtype=np.int64)
    lib = np.zeros(len(keep), dtype=np.float64)
    gene_rows = {}

    with gzip.open(args.umi, "rt") as handle:
        handle.readline()
        n_seen = 0
        for line in handle:
            n_seen += 1
            tab = line.find("\t")
            if tab < 0:
                continue
            gene = line[:tab].strip()
            key = gene.upper()
            vals = np.fromstring(line[tab + 1:], sep="\t", dtype=np.float32)
            if vals.size != len(cell_ids):
                if vals.size < len(cell_ids):
                    vals = np.pad(vals, (0, len(cell_ids) - vals.size))
                else:
                    vals = vals[: len(cell_ids)]
            took = vals[src_cols].astype(np.float64, copy=False)
            lib += took
            if key in wanted:
                gene_rows[wanted[key]] = took.astype(np.float32, copy=False)
            if n_seen % 5000 == 0:
                print(f"streamed {n_seen} genes, kept {len(gene_rows)}", flush=True)
    print(f"stream done genes_seen {n_seen} kept {len(gene_rows)}", flush=True)
    if not gene_rows:
        sys.exit("none of the wanted genes were found")

    genes = sorted(gene_rows)
    mat = np.vstack([gene_rows[g] for g in genes])
    sm = sparse.csr_matrix(mat)
    mmwrite(os.path.join(args.out, "matrix.mtx"), sm)
    with open(os.path.join(args.out, "genes.tsv"), "w") as handle:
        handle.write("\n".join(genes) + "\n")
    with open(os.path.join(args.out, "cells.tsv"), "w") as handle:
        handle.write("cell\tsample\tcell_type\tcell_subtype\tmalignant\ttnk\tlibsize\n")
        for (cid, _i, sample, ctype, subtype, mal, tnk), lib_i in zip(keep, lib):
            handle.write(
                f"{cid}\t{sample}\t{ctype}\t{subtype}\t{int(mal)}\t{int(tnk)}\t{lib_i:.6f}\n"
            )
    print(f"wrote {args.out} {len(genes)} x {len(keep)} nnz {sm.nnz}", flush=True)


if __name__ == "__main__":
    main()
