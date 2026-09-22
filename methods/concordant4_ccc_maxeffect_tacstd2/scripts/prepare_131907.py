#!/usr/bin/env python3
"""Write the GSE131907 column map for extract_131907.

Keeps malignant and T/NK cells from locked samples (n_malignant > 0).
Column order matches the UMI matrix header. Compact slot 0..n_kept-1.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
RAW = Path("/tmp/concordant4_raw")
OUT = Path("/tmp/concordant4_cache/gse131907")
MALIG = {"Malignant cells", "tS1", "tS2", "tS3"}
TNK = {"T lymphocytes", "NK cells"}


def main() -> None:
    ann = pd.read_csv(RAW / "GSE131907" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz", sep="\t")
    ann["Index"] = ann["Index"].astype(str)
    ann = ann.drop_duplicates("Index").set_index("Index")
    samples = pd.read_csv(HERE / "data" / "GSE131907_samples.tsv", sep="\t")
    locked = set(samples.loc[samples["n_malignant"] > 0, "sample"].astype(str))
    mtx = RAW / "GSE131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    with gzip.open(mtx, "rb") as fh:
        header = fh.readline().rstrip(b"\n").split(b"\t")
    cells = [h.decode() for h in header[1:]]
    if len(cells) != len(ann):
        raise RuntimeError(f"header cells {len(cells)} vs annotation {len(ann)}")
    subtype = ann["Cell_subtype"].astype(str)
    ctype = ann["Cell_type"].astype(str)
    sample = ann["Sample"].astype(str)
    keep_map = np.full(len(cells), -1, dtype=np.int32)
    rows = []
    slot = 0
    for i, cell in enumerate(cells):
        if sample.at[cell] not in locked:
            continue
        if subtype.at[cell] in MALIG:
            comp = "MAL"
        elif ctype.at[cell] in TNK:
            comp = "TNK"
        else:
            continue
        keep_map[i] = slot
        rows.append((cell, sample.at[cell], comp))
        slot += 1
    OUT.mkdir(parents=True, exist_ok=True)
    keep_map.tofile(OUT / "keep_idx.i32")
    pd.DataFrame(rows, columns=["cell", "sample", "comp"]).to_csv(OUT / "cells.tsv", sep="\t", index=False)
    genes = sorted(
        {
            "EPCAM", "KRT8", "KRT18", "KRT19", "PTPRC", "CD3D", "CD3E", "CD8A", "NKG7", "GNLY",
            "KLRD1", "CLDN4", "TACSTD2", "F11R", "JAM1", "NECTIN2", "PVRL2", "CDH1", "LGALS9", "ITGAL",
            "ITGB2", "TIGIT", "CD96", "ITGAE", "ITGB7", "KLRG1", "HAVCR2", "CD44", "CXCL9",
            "CXCL10", "CXCL11", "CCL5", "CXCL16", "CXCR3", "CCR5", "CCR1", "CXCR6", "HLA-A",
            "HLA-B", "HLA-C",
        }
    )
    (OUT / "genes_request.txt").write_text("\n".join(genes) + "\n")
    print(f"kept {slot} / {len(cells)}")


if __name__ == "__main__":
    main()
