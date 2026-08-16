#!/usr/bin/env python3
"""Protein–RNA concordance for TACSTD2 / CLDN4 (and any extra genes)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from lib_io import (  # noqa: E402
    apply_log2_if_needed,
    lookup_gene_rows,
    normalize_sample_id,
    read_expression_matrix,
)


def _vector(matrix: pd.DataFrame, gene: str) -> pd.Series | None:
    rows = lookup_gene_rows(matrix, gene)
    if rows.empty:
        return None
    s = rows.iloc[0].copy()
    s.index = [normalize_sample_id(i) for i in s.index]
    return s


def concordance(protein: pd.DataFrame, rna: pd.DataFrame, gene: str, cohort: str) -> dict:
    p = _vector(protein, gene)
    r = _vector(rna, gene)
    if p is None and r is None:
        return {
            "cohort": cohort,
            "gene": gene,
            "status": "absent_both",
            "n_paired": 0,
            "spearman_rho": np.nan,
            "spearman_p": np.nan,
            "note": "gene missing from protein and RNA tables",
        }
    if p is None:
        return {
            "cohort": cohort,
            "gene": gene,
            "status": "protein_absent",
            "n_rna": int(r.notna().sum()) if r is not None else 0,
            "n_paired": 0,
            "spearman_rho": np.nan,
            "spearman_p": np.nan,
            "note": (
                "RNA present, protein row absent (typical CLDN4 TMT/NArm). "
                "Concordance is undefined — do not impute protein from RNA."
            ),
        }
    if r is None:
        return {
            "cohort": cohort,
            "gene": gene,
            "status": "rna_absent",
            "n_protein": int(p.notna().sum()),
            "n_paired": 0,
            "spearman_rho": np.nan,
            "spearman_p": np.nan,
            "note": "protein present, RNA row absent",
        }
    paired = pd.concat([p.rename("protein"), r.rename("rna")], axis=1, join="inner").dropna()
    if paired.shape[0] < 5:
        return {
            "cohort": cohort,
            "gene": gene,
            "status": "too_few_pairs",
            "n_paired": int(paired.shape[0]),
            "spearman_rho": np.nan,
            "spearman_p": np.nan,
            "note": "fewer than 5 paired non-missing samples",
        }
    rho, pval = stats.spearmanr(paired["protein"], paired["rna"])
    return {
        "cohort": cohort,
        "gene": gene,
        "status": "ok",
        "n_paired": int(paired.shape[0]),
        "spearman_rho": float(rho),
        "spearman_p": float(pval),
        "protein_median": float(paired["protein"].median()),
        "rna_median": float(paired["rna"].median()),
        "note": "Spearman on paired samples; IDs harmonized to C3L/C3N",
    }


def paired_table(protein: pd.DataFrame, rna: pd.DataFrame, gene: str, cohort: str) -> pd.DataFrame:
    p = _vector(protein, gene)
    r = _vector(rna, gene)
    if p is None or r is None:
        return pd.DataFrame()
    out = pd.concat([p.rename("protein"), r.rename("rna")], axis=1, join="inner")
    out["gene"] = gene
    out["cohort"] = cohort
    out.index.name = "sample_id"
    return out.reset_index()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--protein", action="append", nargs=2, metavar=("COHORT", "PATH"), required=True)
    ap.add_argument("--rna", action="append", nargs=2, metavar=("COHORT", "PATH"), required=True)
    ap.add_argument("--genes", nargs="+", default=["TACSTD2", "CLDN4"])
    ap.add_argument("--out-stats", type=Path, required=True)
    ap.add_argument("--out-pairs", type=Path, required=True)
    args = ap.parse_args(argv)

    prot = {c: apply_log2_if_needed(read_expression_matrix(p), name=f"{c}_protein")[0] for c, p in args.protein}
    rna = {c: apply_log2_if_needed(read_expression_matrix(p), name=f"{c}_rna")[0] for c, p in args.rna}

    stats_rows = []
    pair_frames = []
    for cohort, pmat in prot.items():
        if cohort not in rna:
            raise SystemExit(f"no RNA matrix for cohort {cohort}")
        for gene in args.genes:
            stats_rows.append(concordance(pmat, rna[cohort], gene, cohort))
            pair_frames.append(paired_table(pmat, rna[cohort], gene, cohort))

    args.out_stats.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(stats_rows).to_csv(args.out_stats, sep="\t", index=False)
    pd.concat([f for f in pair_frames if not f.empty], ignore_index=True).to_csv(args.out_pairs, sep="\t", index=False)
    print(f"wrote {args.out_stats} and {args.out_pairs}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
