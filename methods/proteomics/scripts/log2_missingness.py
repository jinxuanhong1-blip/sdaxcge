#!/usr/bin/env python3
"""log2-scale check + missingness, with CLDN4 absence as the worked example."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from lib_io import (  # noqa: E402
    TARGET_GENES,
    apply_log2_if_needed,
    detect_log_scale,
    gene_presence,
    read_expression_matrix,
    read_gdc_star_gene_counts,
    read_pride_protein_table,
)


CLAUDIN_PANEL = [
    "TACSTD2",
    "CLDN1",
    "CLDN2",
    "CLDN3",
    "CLDN4",
    "CLDN5",
    "CLDN7",
    "CLDN12",
    "CLDN18",
    "EPCAM",
    "KRT7",
]


def gene_wise_missingness(matrix: pd.DataFrame) -> pd.Series:
    return matrix.isna().mean(axis=1)


def summarize_matrix(path: Path, kind: str) -> dict:
    if kind == "pride":
        mat = read_pride_protein_table(path)
    elif kind == "star":
        star = read_gdc_star_gene_counts(path)
        # presence of the two Ensembl IDs
        info = {
            "path": str(path),
            "kind": "gdc_star",
            "n_genes": int(star.shape[0]),
            "targets": {},
        }
        for symbol, meta in TARGET_GENES.items():
            ens = meta["ensembl"]
            info["targets"][symbol] = {
                "ensembl": ens,
                "present": ens in star.index,
                "counts": None if ens not in star.index else float(star.loc[ens, "counts"]),
                "tpm": None if ens not in star.index else float(star.loc[ens, "tpm"]),
                "fpkm": None if ens not in star.index else float(star.loc[ens, "fpkm"]),
            }
        vals = star["tpm"].dropna()
        info["log2_check"] = detect_log_scale(vals, name=f"{path.name}:tpm")
        return info
    else:
        mat = read_expression_matrix(path)

    logged, log_info = apply_log2_if_needed(mat, name=path.name)
    na = gene_wise_missingness(mat)
    targets = {g: gene_presence(mat, g) for g in TARGET_GENES}
    panel = {g: gene_presence(mat, g) for g in CLAUDIN_PANEL}
    return {
        "path": str(path),
        "kind": kind,
        "n_genes": int(mat.shape[0]),
        "n_samples": int(mat.shape[1]),
        "log2_check": log_info,
        "missingness": {
            "median_gene_na_frac": float(na.median()),
            "p90_gene_na_frac": float(na.quantile(0.9)),
            "max_gene_na_frac": float(na.max()),
            "n_genes_all_na": int((na >= 1.0).sum()),
            "n_genes_absent_targets_note": (
                "CLDN4 is often *missing as a row* after NArm, not a row of NAs. "
                "That is worse than MNAR: the protein never enters the matrix."
            ),
        },
        "targets": targets,
        "claudin_panel": panel,
        "double_log_warning": (
            None
            if not log_info.get("already_log")
            else "Values already look log-scaled. Do not apply log2 again."
        ),
        "logged_preview_used": log_info.get("already_log", True),
        "n_logged_equal_raw": int(logged.shape[0]),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("tables", nargs="+", type=Path, help="Protein/RNA/STAR/PRIDE tables")
    ap.add_argument("--kind", default="auto", choices=["auto", "matrix", "pride", "star"])
    ap.add_argument("--out-json", type=Path, required=True)
    ap.add_argument("--out-tsv", type=Path, default=None)
    args = ap.parse_args(argv)

    reports = []
    rows = []
    for path in args.tables:
        kind = args.kind
        if kind == "auto":
            name = path.name.lower()
            if "star" in name or name.endswith("gene_counts.tsv"):
                kind = "star"
            elif "proteingroups" in name or name.startswith("report.pg"):
                kind = "pride"
            else:
                kind = "matrix"
        rep = summarize_matrix(path, kind)
        reports.append(rep)
        if "targets" in rep:
            for gene, info in rep["targets"].items():
                row = {"table": path.name, "gene": gene, **info}
                rows.append(row)

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(reports, indent=2), encoding="utf-8")
    if args.out_tsv:
        pd.json_normalize(rows).to_csv(args.out_tsv, sep="\t", index=False)
    print(f"wrote {args.out_json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
