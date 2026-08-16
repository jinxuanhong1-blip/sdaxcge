#!/usr/bin/env python3
"""Stream GSE207422 genes x cells UMI txt.gz; keep panel genes + per-cell total UMI.

Does not load the full dense matrix. Output is cells x panel genes.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path("/tmp/scrna_data")
MATRIX = DATA / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
PANEL = Path("/workspace/scripts/scrna/gene_panel.tsv")
OUT = DATA / "gse207422_panel_cells.tsv.gz"


def main() -> None:
    panel = set(pd.read_csv(PANEL, sep="\t")["gene"])
    print(f"panel size {len(panel)}", flush=True)

    with gzip.open(MATRIX, "rb") as fh:
        header = fh.readline().decode("ascii").rstrip("\n")
        cells = header.split("\t")[1:]
        n_cells = len(cells)
        print(f"cells {n_cells}", flush=True)

        total = np.zeros(n_cells, dtype=np.float64)
        kept: dict[str, np.ndarray] = {}
        n_genes = 0
        for raw in fh:
            n_genes += 1
            tab = raw.find(b"\t")
            gene = raw[:tab].decode("ascii")
            arr = np.fromstring(raw[tab + 1 :].rstrip(b"\r\n"), sep="\t", dtype=np.float32)
            if arr.size != n_cells:
                raise ValueError(f"{gene}: got {arr.size} values, expected {n_cells}")
            total += arr
            if gene in panel:
                kept[gene] = arr.copy()
            if n_genes % 2000 == 0:
                print(f"  genes={n_genes} kept={len(kept)}", flush=True)

    missing = sorted(panel - set(kept))
    print(f"done genes={n_genes} kept={len(kept)} missing={missing}", flush=True)

    df = pd.DataFrame({"barcode": cells, "total_umi": total.astype(np.int64)})
    for g in sorted(kept):
        df[g] = kept[g].astype(np.int32)
    df.to_csv(OUT, sep="\t", index=False, compression="gzip")
    print(f"wrote {OUT} rows={len(df)}", flush=True)


if __name__ == "__main__":
    main()
