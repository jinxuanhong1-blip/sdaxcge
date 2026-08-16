#!/usr/bin/env python3
"""Stream GSE207422 genes×cells UMI TSV; keep panel genes + n_genes + total UMI.

Does not materialize the full dense matrix.
"""
from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from pathlib import Path as _P
import sys as _sys
_sys.path.insert(0, str(_P(__file__).resolve().parent))
from gene_sets import all_panel_genes

DATA = Path("/tmp/scrna_cytotrace")
MATRIX = DATA / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
OUT = DATA / "gse207422_panel_cells.tsv.gz"


def main() -> None:
    panel = set(all_panel_genes())
    print(f"panel={len(panel)}", flush=True)
    with gzip.open(MATRIX, "rb") as fh:
        header = fh.readline().decode("ascii").rstrip("\n")
        cells = header.split("\t")[1:]
        n_cells = len(cells)
        print(f"cells={n_cells}", flush=True)
        total = np.zeros(n_cells, dtype=np.float64)
        n_genes = np.zeros(n_cells, dtype=np.int32)
        kept: dict[str, np.ndarray] = {}
        n_seen = 0
        for raw in fh:
            n_seen += 1
            tab = raw.find(b"\t")
            gene = raw[:tab].decode("ascii")
            arr = np.fromstring(raw[tab + 1 :].rstrip(b"\r\n"), sep="\t", dtype=np.float32)
            if arr.size != n_cells:
                raise ValueError(f"{gene}: got {arr.size}, expected {n_cells}")
            total += arr
            n_genes += arr > 0
            if gene in panel:
                kept[gene] = arr.copy()
            if n_seen % 2000 == 0:
                print(f"  genes={n_seen} kept={len(kept)}", flush=True)
    missing = sorted(panel - set(kept))
    print(f"done genes={n_seen} kept={len(kept)} missing={missing}", flush=True)
    df = pd.DataFrame({
        "barcode": cells,
        "total_umi": total.astype(np.int64),
        "n_genes": n_genes,
    })
    for g in sorted(kept):
        df[g] = kept[g].astype(np.int32)
    df.to_csv(OUT, sep="\t", index=False, compression="gzip")
    print(f"wrote {OUT} rows={len(df)}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
