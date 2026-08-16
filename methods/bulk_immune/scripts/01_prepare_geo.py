#!/usr/bin/env python3
"""Download GSE126044 + GSE135222, map identifiers, write prepared matrices.

Outputs under --dest (default methods/bulk_immune/data/):

    GSE126044.counts.tsv.gz
    GSE126044.log2cpm.tsv.gz
    GSE126044.pheno.tsv
    GSE135222.fpkm.tsv.gz
    GSE135222.log2fpkm.tsv.gz
    GSE135222.pheno.tsv
    prepare_geo.log.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from bulkimmune.genes import HGNC, collapse_duplicates, strip_ensembl_version  # noqa: E402
from bulkimmune.geo import (  # noqa: E402
    GSE126044_COUNTS,
    GSE126044_MATRIX,
    GSE135222_EXPR,
    GSE135222_MATRIX,
    download,
    load_gse126044,
    load_gse135222,
)
from bulkimmune.preprocess import cpm, detect_scale, log2p1  # noqa: E402


def _maybe_hgnc(resource_dir: Path) -> HGNC | None:
    path = resource_dir / "hgnc_complete_set.txt"
    if not path.exists():
        print("[prepare] no HGNC table; Ensembl IDs will be version-stripped only")
        return None
    return HGNC.from_file(path)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dest", type=Path, default=ROOT / "data")
    p.add_argument("--resources", type=Path, default=ROOT / "resources")
    args = p.parse_args()
    args.dest.mkdir(parents=True, exist_ok=True)
    hgnc = _maybe_hgnc(args.resources)
    log: dict = {}

    print("[prepare] GSE126044")
    cpath = download(GSE126044_COUNTS, args.dest / "raw" / "GSE126044_counts.txt.gz")
    mpath = download(GSE126044_MATRIX, args.dest / "raw" / "GSE126044_series_matrix.txt.gz")
    counts, pheno = load_gse126044(cpath, mpath)
    if hgnc is not None:
        mapped, renamed = hgnc.harmonise(list(counts.index))
        new_index = [m if m else o for m, o in zip(mapped, counts.index)]
        counts.index = new_index
        log["GSE126044_renamed_symbols"] = renamed
    counts = collapse_duplicates(counts, method="max_mean")
    log2cpm = cpm(counts, log=True, prior_count=1.0)
    counts.to_csv(args.dest / "GSE126044.counts.tsv.gz", sep="\t", compression="gzip")
    log2cpm.to_csv(args.dest / "GSE126044.log2cpm.tsv.gz", sep="\t", compression="gzip")
    pheno.to_csv(args.dest / "GSE126044.pheno.tsv", sep="\t")
    log["GSE126044"] = {
        "n_genes": int(counts.shape[0]),
        "n_samples": int(counts.shape[1]),
        "response": pheno["response"].value_counts().to_dict(),
        "tissue_preservation": pheno["tissue_preservation"].value_counts().to_dict(),
        "scale": detect_scale(counts),
        "targets_present": {
            g: bool(g in log2cpm.index) for g in ("TACSTD2", "CLDN4", "GZMA", "PRF1")
        },
    }
    print("  ", log["GSE126044"])

    print("[prepare] GSE135222")
    epath = download(GSE135222_EXPR, args.dest / "raw" / "GSE135222_expr.tsv.gz")
    mpath = download(GSE135222_MATRIX, args.dest / "raw" / "GSE135222_series_matrix.txt.gz")
    expr, pheno2 = load_gse135222(epath, mpath)
    raw_ids = list(expr.index)
    stripped = strip_ensembl_version(raw_ids)
    if hgnc is not None:
        symbols = hgnc.ensembl_to_symbol(stripped)
        keep_idx, keep_sym = [], []
        for i, (ens, sym) in enumerate(zip(stripped, symbols)):
            if sym:
                keep_idx.append(i)
                keep_sym.append(sym)
        expr = expr.iloc[keep_idx]
        expr.index = keep_sym
        log["GSE135222_unmapped_ensembl"] = len(raw_ids) - len(keep_idx)
    else:
        expr.index = stripped
    expr = collapse_duplicates(expr, method="max_mean")
    log2 = log2p1(expr)
    expr.to_csv(args.dest / "GSE135222.fpkm.tsv.gz", sep="\t", compression="gzip")
    log2.to_csv(args.dest / "GSE135222.log2fpkm.tsv.gz", sep="\t", compression="gzip")
    pheno2.to_csv(args.dest / "GSE135222.pheno.tsv", sep="\t")
    log["GSE135222"] = {
        "n_genes": int(expr.shape[0]),
        "n_samples": int(expr.shape[1]),
        "pfs_events": int(pheno2["pfs_event"].sum()),
        "pfs_event_rate": float(pheno2["pfs_event"].mean()),
        "sex": pheno2["sex"].value_counts().to_dict(),
        "scale": detect_scale(expr),
        "targets_present": {
            g: bool(g in log2.index) for g in ("TACSTD2", "CLDN4", "GZMA", "PRF1")
        },
    }
    print("  ", log["GSE135222"])

    (args.dest / "prepare_geo.log.json").write_text(json.dumps(log, indent=2, default=str))
    print("wrote", args.dest)


if __name__ == "__main__":
    main()
