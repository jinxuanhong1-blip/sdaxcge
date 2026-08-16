#!/usr/bin/env python3
"""Join protein/RNA target genes to CIBERSORT, xCell, or LinkedOmics immune clusters.

Sample IDs are normalized (C3L.00001 == C3L-00001). The join never invents
ICI response labels — CPTAC molecular phenotypes are treatment-naive.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from lib_io import (  # noqa: E402
    lookup_gene_rows,
    normalize_sample_id,
    read_cibersort,
    read_expression_matrix,
    read_linkedomics_phenotype,
    read_xcell,
)


RNA_IMMUNE_MARKERS = {
    "CD8A": "CD8_T",
    "CD4": "CD4",
    "FOXP3": "Treg",
    "CD274": "PDL1",
    "CXCL13": "TLS_proxy",
    "CD68": "Macrophage",
    "NCAM1": "NK_proxy",
}


def _gene_vector(matrix: pd.DataFrame, gene: str, prefix: str) -> pd.Series:
    rows = lookup_gene_rows(matrix, gene)
    if rows.empty:
        return pd.Series(dtype=float, name=f"{prefix}_{gene}")
    s = rows.iloc[0].copy()
    s.index = [normalize_sample_id(i) for i in s.index]
    s.name = f"{prefix}_{gene}"
    return s


def marker_scores(rna: pd.DataFrame) -> pd.DataFrame:
    cols = {}
    for gene, label in RNA_IMMUNE_MARKERS.items():
        rows = lookup_gene_rows(rna, gene)
        if rows.empty:
            continue
        s = rows.iloc[0].copy()
        s.index = [normalize_sample_id(i) for i in s.index]
        cols[f"rna_{label}"] = s
    if not cols:
        return pd.DataFrame()
    return pd.DataFrame(cols)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--protein", type=Path, action="append", default=[], help="genes x samples")
    ap.add_argument("--rna", type=Path, action="append", default=[], help="genes x samples")
    ap.add_argument("--genes", nargs="+", default=["TACSTD2", "CLDN4"])
    ap.add_argument("--cibersort", type=Path, default=None)
    ap.add_argument("--xcell", type=Path, default=None)
    ap.add_argument("--phenotype", type=Path, default=None, help="LinkedOmics .tsi (e.g. Immune.Cluster.rna)")
    ap.add_argument("--cohort", default="")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)

    pieces: list[pd.DataFrame] = []
    for i, path in enumerate(args.protein):
        mat = read_expression_matrix(path)
        tag = args.cohort or f"p{i}"
        for gene in args.genes:
            pieces.append(_gene_vector(mat, gene, f"protein_{tag}").to_frame())
    rna_frames = []
    for i, path in enumerate(args.rna):
        mat = read_expression_matrix(path)
        tag = args.cohort or f"r{i}"
        for gene in args.genes:
            pieces.append(_gene_vector(mat, gene, f"rna_{tag}").to_frame())
        rna_frames.append(mat)
    if rna_frames:
        pieces.append(marker_scores(rna_frames[0]))
    if args.cibersort:
        pieces.append(read_cibersort(args.cibersort).add_prefix("cibersort_"))
    if args.xcell:
        pieces.append(read_xcell(args.xcell).add_prefix("xcell_"))
    if args.phenotype:
        ph = read_linkedomics_phenotype(args.phenotype)
        pieces.append(ph.add_prefix("pheno_"))

    if not pieces:
        raise SystemExit("nothing to join")
    out = pieces[0]
    for extra in pieces[1:]:
        out = out.join(extra, how="outer")
    out.index.name = "sample_id"
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.reset_index().to_csv(args.out, sep="\t", index=False)
    print(f"wrote {args.out}  n={out.shape[0]} p={out.shape[1]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
