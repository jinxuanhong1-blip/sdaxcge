#!/usr/bin/env python3
"""Convert NicheNet-v2 ligand-target RDS to parquet (no R required)."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import rdata


def convert_rds(path: Path):
    print(f"reading {path} ({path.stat().st_size} bytes)", flush=True)
    parsed = rdata.parser.parse_file(path)
    return rdata.conversion.convert(parsed)


def as_dataframe(obj) -> pd.DataFrame:
    if isinstance(obj, pd.DataFrame):
        return obj
    if hasattr(obj, "to_pandas"):
        return obj.to_pandas()
    if isinstance(obj, dict) and all(hasattr(v, "__len__") for v in obj.values()):
        return pd.DataFrame(obj)
    raise TypeError(f"unsupported RDS type: {type(obj)}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--lt-rds", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    lt = convert_rds(args.lt_rds)
    if isinstance(lt, pd.DataFrame):
        mat = lt
    elif hasattr(lt, "to_pandas"):
        mat = lt.to_pandas()
    elif hasattr(lt, "coords") and hasattr(lt, "values"):
        rows = [str(x) for x in lt.coords[lt.dims[0]].values]
        cols = [str(x) for x in lt.coords[lt.dims[1]].values]
        mat = pd.DataFrame(lt.values, index=rows, columns=cols)
    else:
        mat = pd.DataFrame(lt)
    mat.index = mat.index.astype(str)
    mat.columns = mat.columns.astype(str)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    mat.to_parquet(args.out)
    print(f"ligand_target {mat.shape} -> {args.out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
