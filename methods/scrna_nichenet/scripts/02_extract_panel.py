#!/usr/bin/env python3
"""Stream GSE207422 dense UMI matrix; keep ligands, receptors, T/NK programs.

Line-by-line gzip parse: only requested gene rows are tokenized.
"""
from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
CACHE = Path("/tmp/scrna_nichenet")
sys.path.insert(0, str(ROOT / "scripts"))
from gene_sets import CYTOTOXICITY, EXHAUSTION, EXTRA_TNK, LINEAGE  # noqa: E402

MATRIX = CACHE / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
LR = DATA / "lr_network.tsv"
OUT = CACHE / "extracted_panel.parquet"


def panel_genes() -> set[str]:
    lr = pd.read_csv(LR, sep="\t")
    genes = set(lr["from"].astype(str)) | set(lr["to"].astype(str))
    for block in (CYTOTOXICITY, EXHAUSTION, EXTRA_TNK, LINEAGE):
        genes.update(block)
    return genes


def main() -> int:
    want = panel_genes()
    print(f"panel requested: {len(want)}", flush=True)
    kept: dict[str, np.ndarray] = {}
    n_genes = 0
    barcodes: list[str] | None = None
    with gzip.open(MATRIX, "rt") as fh:
        header = fh.readline().rstrip("\n")
        barcodes = header.split("\t")[1:]
        print(f"cells: {len(barcodes)}", flush=True)
        for line in fh:
            n_genes += 1
            tab = line.find("\t")
            if tab < 0:
                continue
            gene = line[:tab]
            if gene in want:
                vals = np.fromstring(line[tab + 1 :], sep="\t", dtype=np.int32)
                kept[gene] = vals
            if n_genes % 2000 == 0:
                print(f"genes {n_genes}; kept {len(kept)}", flush=True)
    missing = sorted(want - set(kept))
    print(f"total genes {n_genes}; present {len(kept)}; missing {len(missing)}", flush=True)
    if missing[:20]:
        print("missing example:", missing[:20], flush=True)
    df = pd.DataFrame(kept, index=barcodes)
    df.index.name = "barcode"
    df.to_parquet(OUT)
    print(f"wrote {OUT} shape={df.shape} bytes={OUT.stat().st_size}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
