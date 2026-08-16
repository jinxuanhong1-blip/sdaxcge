#!/usr/bin/env python3
"""Stream the GSE207422 dense UMI matrix (genes x cells, txt.gz) and extract
panel-gene rows plus per-cell total UMI counts.
Output: /tmp/data_scrna_ici/gse207422_panel_cells.tsv.gz
"""
import gzip
from collections import defaultdict
import numpy as np
import pandas as pd

DATA = "/tmp/data_scrna_ici"
MATRIX = f"{DATA}/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
PANEL = "/workspace/scripts/fable_scrna_ici/gene_panel.tsv"
OUT = f"{DATA}/gse207422_panel_cells.tsv.gz"

panel = set(pd.read_csv(PANEL, sep="\t")["gene"])

with gzip.open(MATRIX, "rt") as fh:
    cells = fh.readline().rstrip("\n").split("\t")[1:]
n_cells = len(cells)
print(f"cells: {n_cells}", flush=True)

total = np.zeros(n_cells, dtype=np.float64)
kept = {}
n_genes = 0
dtypes = defaultdict(lambda: np.int32)
dtypes["Gene"] = str
reader = pd.read_csv(MATRIX, sep="\t", index_col=0, chunksize=1500,
                     dtype=dtypes)
for chunk in reader:
    n_genes += chunk.shape[0]
    total += chunk.to_numpy().sum(axis=0)
    hits = chunk.index.intersection(panel)
    for g in hits:
        kept[g] = chunk.loc[g].to_numpy()
    print(f"processed {n_genes} genes, kept {len(kept)}", flush=True)

missing = sorted(panel - kept.keys())
print(f"total genes: {n_genes}; panel present: {len(kept)}; missing: {missing}")

df = pd.DataFrame({"barcode": cells, "total_umi": total.astype(np.int64)})
for g in sorted(kept):
    df[g] = kept[g]
df.to_csv(OUT, sep="\t", index=False, compression="gzip")
print(f"wrote {OUT}")
