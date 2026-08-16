#!/usr/bin/env python3
"""Download TCGA LUAD + LUSC from UCSC Xena and write prepared matrices.

Xena HiSeqV2 is already log2(RSEM+1) on current gene symbols. A linear
reconstruction ``2**x - 1`` is written for CIBERSORT / CYT / TIDE (TIDE will
re-log and row-centre it).

Outputs under --dest:

    TCGA_NSCLC.log2rsem.tsv.gz
    TCGA_NSCLC.linear.tsv.gz
    TCGA_NSCLC.pheno.tsv
    prepare_tcga.log.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from bulkimmune.genes import collapse_duplicates  # noqa: E402
from bulkimmune.preprocess import detect_scale  # noqa: E402
from bulkimmune.tcga import load_tcga_nsclc  # noqa: E402


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dest", type=Path, default=ROOT / "data")
    p.add_argument("--cohorts", nargs="+", default=["LUAD", "LUSC"])
    p.add_argument("--no-fetch", action="store_true")
    args = p.parse_args()
    args.dest.mkdir(parents=True, exist_ok=True)

    expr, pheno = load_tcga_nsclc(
        args.dest / "raw",
        cohorts=tuple(args.cohorts),
        fetch=not args.no_fetch,
    )
    expr = collapse_duplicates(expr, method="max_mean")
    # Xena HiSeqV2 is log2(RSEM + 1). Recover a linear scale for methods that
    # need it; zeros stay zeros.
    linear = np.clip(np.power(2.0, expr) - 1.0, a_min=0.0, a_max=None)

    expr.to_csv(args.dest / "TCGA_NSCLC.log2rsem.tsv.gz", sep="\t", compression="gzip")
    linear.to_csv(args.dest / "TCGA_NSCLC.linear.tsv.gz", sep="\t", compression="gzip")
    pheno.to_csv(args.dest / "TCGA_NSCLC.pheno.tsv", sep="\t")

    log = {
        "n_genes": int(expr.shape[0]),
        "n_samples": int(expr.shape[1]),
        "cohorts": pheno["cohort"].value_counts().to_dict(),
        "sex": pheno["sex"].value_counts(dropna=False).to_dict(),
        "scale_log": detect_scale(expr),
        "scale_linear": detect_scale(linear),
        "targets_present": {g: bool(g in expr.index) for g in ("TACSTD2", "CLDN4")},
        "note": (
            "TCGA is treatment-naive surgical NSCLC, not an ICI cohort. "
            "Use it for TACSTD2/CLDN4 ~ immune-score correlations; do not "
            "quote it as ICI-response evidence."
        ),
    }
    (args.dest / "prepare_tcga.log.json").write_text(json.dumps(log, indent=2, default=str))
    print(json.dumps(log, indent=2, default=str))


if __name__ == "__main__":
    main()
