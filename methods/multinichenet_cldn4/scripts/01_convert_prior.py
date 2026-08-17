#!/usr/bin/env python3
"""Convert the NicheNet-v2 ligand–target RDS to parquet (no R required)."""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import rdata

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_io import CACHE, log  # noqa: E402

LT_RDS = CACHE / "ligand_target_matrix_nsga2r_final.rds"
LT_OUT = CACHE / "prior_ligand_target.parquet"


def convert_rds(path: Path):
    log(f"reading {path} ({path.stat().st_size} bytes)")
    parsed = rdata.parser.parse_file(path)
    return rdata.conversion.convert(parsed)


def as_matrix(obj) -> pd.DataFrame:
    if isinstance(obj, pd.DataFrame):
        return obj
    if hasattr(obj, "to_pandas"):
        return obj.to_pandas()
    if hasattr(obj, "coords") and hasattr(obj, "values"):
        rows = [str(x) for x in obj.coords[obj.dims[0]].values]
        cols = [str(x) for x in obj.coords[obj.dims[1]].values]
        return pd.DataFrame(obj.values, index=rows, columns=cols)
    return pd.DataFrame(obj)


def main() -> int:
    if LT_OUT.exists() and LT_OUT.stat().st_size > 10_000_000:
        log(f"HAVE {LT_OUT} ({LT_OUT.stat().st_size} bytes)")
        return 0
    mat = as_matrix(convert_rds(LT_RDS))
    mat.index = mat.index.astype(str)
    mat.columns = mat.columns.astype(str)
    log(f"ligand_target {mat.shape} index_ex={list(mat.index[:3])} cols_ex={list(mat.columns[:3])}")
    # NicheNet-v2: targets in rows, ligands in columns
    if mat.shape[0] < mat.shape[1]:
        log("warning: more columns than rows; leaving orientation as converted")
    mat.to_parquet(LT_OUT)
    log(f"wrote {LT_OUT} ({LT_OUT.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
