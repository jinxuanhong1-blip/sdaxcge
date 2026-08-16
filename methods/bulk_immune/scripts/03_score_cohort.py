#!/usr/bin/env python3
"""Score one prepared cohort (GEO or TCGA) with every available method.

Example:

    python scripts/03_score_cohort.py \\
        --name GSE126044 \\
        --log data/GSE126044.log2cpm.tsv.gz \\
        --linear data/GSE126044.counts.tsv.gz \\
        --linear-is-counts \\
        --pheno data/GSE126044.pheno.tsv
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from bulkimmune.pipeline import score_all  # noqa: E402
from bulkimmune.preprocess import cpm  # noqa: E402


def _read_matrix(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", index_col=0, compression="infer")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--name", required=True)
    p.add_argument("--log", type=Path, required=True)
    p.add_argument("--linear", type=Path, required=True)
    p.add_argument("--linear-is-counts", action="store_true")
    p.add_argument("--pheno", type=Path, default=None)
    p.add_argument("--resources", type=Path, default=ROOT / "resources")
    p.add_argument("--outdir", type=Path, default=None)
    p.add_argument("--skip-xcell", action="store_true")
    p.add_argument("--skip-tide", action="store_true")
    p.add_argument("--lm22", type=Path, default=None)
    p.add_argument("--permutations", type=int, default=0)
    args = p.parse_args()

    outdir = args.outdir or (ROOT / "results" / args.name)
    outdir.mkdir(parents=True, exist_ok=True)

    expr_log = _read_matrix(args.log)
    expr_linear = _read_matrix(args.linear)
    if args.linear_is_counts:
        expr_linear = cpm(expr_linear, log=False)

    result = score_all(
        expr_log,
        expr_linear,
        resource_dir=args.resources,
        run_xcell=not args.skip_xcell,
        run_tide=not args.skip_tide,
        lm22_path=args.lm22,
        permutations=args.permutations,
    )

    result["scores"].to_csv(outdir / "scores.tsv", sep="\t")
    result["targets"].to_csv(outdir / "targets.tsv", sep="\t")
    for name, frame in result["blocks"].items():
        frame.to_csv(outdir / f"block_{name}.tsv", sep="\t")
    for name, frame in result["reports"].items():
        frame.to_csv(outdir / f"report_{name}.tsv", sep="\t")

    if args.pheno is not None:
        pheno = pd.read_csv(args.pheno, sep="\t", index_col=0)
        pheno.to_csv(outdir / "pheno.tsv", sep="\t")

    meta = {
        "name": args.name,
        "n_samples": int(expr_log.shape[1]),
        "n_genes_log": int(expr_log.shape[0]),
        "n_scores": int(result["scores"].shape[1]),
        "targets": list(result["targets"].columns),
        "blocks": {k: list(v.columns) for k, v in result["blocks"].items()},
        "skipped": result["skipped"],
        "signature_table": result["library"].table().to_dict(orient="records"),
    }
    (outdir / "score_meta.json").write_text(json.dumps(meta, indent=2, default=str))
    print(json.dumps({"name": args.name, "n_scores": meta["n_scores"], "skipped": meta["skipped"]}, indent=2))


if __name__ == "__main__":
    main()
