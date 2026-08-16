#!/usr/bin/env python3
"""Convert NicheNet-v2 RDS priors to TSV/parquet (no R required)."""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import rdata

CACHE = Path("/tmp/scrna_nichenet")
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

LR_RDS = CACHE / "lr_network_human_21122021.rds"
LT_RDS = CACHE / "ligand_target_matrix_nsga2r_final.rds"
LR_OUT = DATA / "lr_network.tsv"
LT_OUT = CACHE / "prior_ligand_target.parquet"


def convert_rds(path: Path):
    print(f"reading {path} ({path.stat().st_size} bytes)", flush=True)
    parsed = rdata.parser.parse_file(path)
    obj = rdata.conversion.convert(parsed)
    return obj


def as_dataframe(obj) -> pd.DataFrame:
    if isinstance(obj, pd.DataFrame):
        return obj
    if hasattr(obj, "to_pandas"):
        return obj.to_pandas()
    # rdata may return a dict of columns or a numpy matrix with dimnames
    if isinstance(obj, dict):
        if all(hasattr(v, "__len__") for v in obj.values()):
            return pd.DataFrame(obj)
    raise TypeError(f"unsupported RDS type: {type(obj)}")


def main() -> int:
    lr = as_dataframe(convert_rds(LR_RDS))
    # NicheNet lr_network: from, to, (optional source/database)
    cols = [c for c in lr.columns if str(c).lower() in {"from", "to", "source", "database"}]
    if "from" not in lr.columns or "to" not in lr.columns:
        # sometimes first two columns
        lr = lr.rename(columns={lr.columns[0]: "from", lr.columns[1]: "to"})
    lr = lr.loc[:, [c for c in ["from", "to"] if c in lr.columns]].drop_duplicates()
    lr.to_csv(LR_OUT, sep="\t", index=False)
    print(f"lr_network {lr.shape} -> {LR_OUT}", flush=True)
    print(lr.head().to_string(), flush=True)

    lt = convert_rds(LT_RDS)
    print(f"ligand_target raw type={type(lt)}", flush=True)
    if isinstance(lt, pd.DataFrame):
        mat = lt
    elif hasattr(lt, "to_pandas"):
        mat = lt.to_pandas()
    elif hasattr(lt, "coords") and hasattr(lt, "values"):
        rows = [str(x) for x in lt.coords[lt.dims[0]].values]
        cols = [str(x) for x in lt.coords[lt.dims[1]].values]
        mat = pd.DataFrame(lt.values, index=rows, columns=cols)
    else:
        try:
            mat = pd.DataFrame(lt)
        except Exception as exc:  # noqa: BLE001
            print(f"DataFrame() failed: {exc}", flush=True)
            print(dir(lt)[:40], flush=True)
            raise
    mat.index = mat.index.astype(str)
    mat.columns = mat.columns.astype(str)
    print(f"ligand_target {mat.shape} index_ex={list(mat.index[:3])} cols_ex={list(mat.columns[:3])}", flush=True)
    mat.to_parquet(LT_OUT)
    print(f"wrote {LT_OUT} ({LT_OUT.stat().st_size} bytes)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
