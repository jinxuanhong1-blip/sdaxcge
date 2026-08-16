#!/usr/bin/env python3
"""Extract panel genes + library size from GSE205335 processed dgCMatrix RDS.

GEO file is gzipped (sometimes double-gzipped). Uses the `rdata` package;
does not invent accessions or labels.
"""
from __future__ import annotations

import gzip
import shutil
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import rdata
from scipy import sparse

DATA = Path("/tmp/scrna_data")
RDS_GZ = DATA / "GSE205335_Lung_IO_UMI_matrix.rds.gz"
PANEL = Path("/workspace/scripts/scrna/gene_panel.tsv")
OUT = DATA / "gse205335_panel_cells.tsv.gz"


def _ungzip_until_rds(src: Path, dest: Path) -> None:
    """Decompress gzip wrapper(s) until the file looks like an RDS (magic 0x1f8b or 'X\\n')."""
    current = src
    tmp_paths = []
    for _ in range(3):
        with open(current, "rb") as fh:
            magic = fh.read(2)
        if magic == b"\x1f\x8b":
            nxt = dest.with_suffix(dest.suffix + f".pass{len(tmp_paths)}")
            with gzip.open(current, "rb") as zin, open(nxt, "wb") as zout:
                shutil.copyfileobj(zin, zout, 16 * 1024 * 1024)
            tmp_paths.append(nxt)
            current = nxt
            continue
        break
    if current != dest:
        shutil.copyfile(current, dest)
    for p in tmp_paths:
        if p.exists() and p != dest:
            p.unlink(missing_ok=True)


def main() -> None:
    panel = list(pd.read_csv(PANEL, sep="\t")["gene"])
    with tempfile.TemporaryDirectory(prefix="gse205335-") as tmp:
        rds_path = Path(tmp) / "matrix.rds"
        print("decompressing RDS...", flush=True)
        _ungzip_until_rds(RDS_GZ, rds_path)
        print(f"rds bytes={rds_path.stat().st_size}", flush=True)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message='Missing constructor for R class "dgCMatrix"')
            obj = rdata.read_rds(rds_path)

    expected = {"i", "p", "Dim", "Dimnames", "x"}
    if not expected.issubset(vars(obj)):
        raise TypeError(f"RDS is not dgCMatrix; keys={list(vars(obj))}")

    genes = np.asarray(obj.Dimnames[0], dtype=str)
    barcodes = np.asarray(obj.Dimnames[1], dtype=str)
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False)
    print(f"matrix {matrix.shape[0]} genes x {matrix.shape[1]} cells", flush=True)

    gene_index = {g: i for i, g in enumerate(genes)}
    present = [g for g in panel if g in gene_index]
    missing = [g for g in panel if g not in gene_index]
    print(f"panel present={len(present)} missing={missing}", flush=True)

    total_umi = np.asarray(matrix.sum(axis=0)).ravel()
    df = pd.DataFrame({"barcode": barcodes, "total_umi": total_umi.astype(np.int64)})
    for g in present:
        df[g] = np.asarray(matrix.getrow(int(gene_index[g])).toarray()).ravel().astype(np.int32)
    df.to_csv(OUT, sep="\t", index=False, compression="gzip")
    print(f"wrote {OUT} rows={len(df)}", flush=True)


if __name__ == "__main__":
    main()
