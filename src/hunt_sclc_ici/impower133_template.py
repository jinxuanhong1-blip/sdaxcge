"""IMpower133 analysis TEMPLATE — runs ONLY on real EGA-controlled data.

IMpower133 per-patient RNA-seq is controlled access (EGA EGAS00001004888).
This template intentionally REFUSES to run on placeholder/simulated inputs so
that no fabricated IMpower133 result can be produced. Obtain the data via EGA,
place the files as described below, and re-run.

Expected inputs (place under results/hunt_sclc_ici/data/impower133/):
  - expression.tsv : genes (rows) x samples (cols), log2(TPM+1) or FPKM.
                     First column = gene symbol. (EGAD00001006927/6928)
  - subtypes.tsv   : columns [sample_id, subtype in {SCLC-A,N,P,I}]
                     (EGAD00001006926)
  - clinical.tsv   : optional, columns [sample_id, arm, os_months, os_event]
                     (arm in {EP+atezo, EP+placebo}) if you also want the
                     treatment-interaction analysis.

Given those, it reuses the exact George scoring/statistics so IMpower133 results
are directly comparable.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

from . import bulk_george as bg
from . import signatures as sig

IMP_DIR = bg.DATA / "impower133"
REQUIRED = ["expression.tsv", "subtypes.tsv"]


class Impower133NotPublic(FileNotFoundError):
    """Raised when EGA-controlled IMpower133 files are absent. Not a bug."""


def _guard():
    missing = [f for f in REQUIRED if not (IMP_DIR / f).exists()]
    if missing:
        msg = (
            "IMpower133 controlled-access data not found. NOT running. "
            f"Missing: {', '.join(missing)}. "
            "This template refuses to fabricate results. "
            "See results/hunt_sclc_ici/IMPOWER133_NOTE.md."
        )
        print("=" * 72)
        print(msg)
        print("=" * 72)
        raise Impower133NotPublic(msg)


def run():
    _guard()
    expr = pd.read_csv(IMP_DIR / "expression.tsv", sep="\t", index_col=0)
    # collapse duplicate symbols and ensure log-scale (heuristic: large values -> log)
    if expr.max().max() > 50:
        expr = np.log2(expr + 1.0)
    z = bg.zscore_genes(expr)

    scores = pd.DataFrame(index=expr.columns)
    for name, genes in sig.MEANZ_SIGNATURES.items():
        scores[name] = bg.mean_z_signature(z, genes, name)
    scores["CYT_log2"] = bg.cyt_score(expr)
    if "CD8A" in expr.index:
        scores["CD8A_log2"] = expr.loc["CD8A"]

    subt = pd.read_csv(IMP_DIR / "subtypes.tsv", sep="\t")
    subt = subt.set_index(subt.columns[0])["subtype"]
    subt = subt.reindex(expr.columns)
    scores["subtype"] = subt

    corr = bg.correlation_table(expr, scores)
    grp = bg.group_comparison(expr, subt)

    outdir = bg.TABLES
    corr.to_csv(outdir / "impower133_correlations.tsv", sep="\t", index=False)
    grp.to_csv(outdir / "impower133_group_comparison.tsv", sep="\t", index=False)
    scores.to_csv(outdir / "impower133_sample_scores.tsv", sep="\t")
    print("IMpower133 real-data analysis complete; tables written to", outdir)


if __name__ == "__main__":
    try:
        run()
    except Impower133NotPublic:
        sys.exit(2)
