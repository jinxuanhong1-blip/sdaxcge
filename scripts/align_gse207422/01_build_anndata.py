"""Stream the GSE207422 dense UMI text matrix into a sparse AnnData (.h5ad).

Input : data/gse207422/scRNAseq_UMI_matrix.txt.gz  (genes x cells, tab-delimited)
        data/gse207422/scRNAseq_metadata.xlsx      (per-sample clinical metadata)
Output: data/gse207422/gse207422_raw.h5ad          (cells x genes, raw counts)
"""
import gzip
import time
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

MAT = "data/gse207422/scRNAseq_UMI_matrix.txt.gz"
META = "data/gse207422/scRNAseq_metadata.xlsx"
OUT = "data/gse207422/gse207422_raw.h5ad"

t0 = time.time()
with gzip.open(MAT, "rt") as fh:
    header = fh.readline().rstrip("\n").split("\t")
    cells = header[1:]
    n_cells = len(cells)
    print(f"cells: {n_cells}")

    genes = []
    rows = []  # per-gene sparse row vectors (1 x n_cells)
    for i, line in enumerate(fh):
        tab = line.index("\t")
        genes.append(line[:tab])
        vals = np.fromstring(line[tab + 1:], sep="\t", dtype=np.float32)
        if vals.shape[0] != n_cells:
            raise ValueError(f"row {i} has {vals.shape[0]} values, expected {n_cells}")
        rows.append(sp.csr_matrix(vals))
        if (i + 1) % 2000 == 0:
            print(f"  parsed {i+1} genes  ({time.time()-t0:.0f}s)")

print(f"parsed all {len(genes)} genes in {time.time()-t0:.0f}s; stacking...")
X = sp.vstack(rows, format="csr")  # genes x cells
del rows
X = X.T.tocsr()  # cells x genes
print("matrix shape (cells x genes):", X.shape)

# per-cell sample id is the prefix before the last underscore, e.g. BD_immune01_612637
sample = ["_".join(c.split("_")[:-1]) for c in cells]
obs = pd.DataFrame({"cell": cells, "Sample": sample}).set_index("cell")

meta = pd.read_excel(META)
meta = meta.set_index("Sample")
for col in meta.columns:
    vals = meta.loc[obs["Sample"].values, col].values
    obs[col] = pd.Series(vals, index=obs.index).astype(str)

adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))
adata.var_names_make_unique()
print(adata)
print(obs["Sample"].value_counts())
adata.write(OUT)
print(f"wrote {OUT} in {time.time()-t0:.0f}s total")
