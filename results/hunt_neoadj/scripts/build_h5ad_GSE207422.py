#!/usr/bin/env python3
"""Load GSE207422 dense TSV UMI matrix (genes x cells) into a sparse AnnData (cells x genes)."""
import time
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

DATA = "results/hunt_neoadj/data"
SRC = f"{DATA}/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
OUT = f"{DATA}/GSE207422_scRNA.h5ad"

t0 = time.time()
blocks = []
genes = []
cells = None
reader = pd.read_csv(SRC, sep="\t", index_col=0, chunksize=1000,
                     dtype=str)  # read as str; convert per-chunk to int
for i, chunk in enumerate(reader):
    if cells is None:
        cells = chunk.columns.to_numpy()
    genes.extend(chunk.index.tolist())
    arr = chunk.to_numpy(dtype=np.float32)  # genes_chunk x cells
    blocks.append(sp.csr_matrix(arr))
    print(f"chunk {i}: genes so far {len(genes)} elapsed {time.time()-t0:.0f}s", flush=True)

X = sp.vstack(blocks).tocsr()  # genes x cells
del blocks
X = X.T.tocsr()  # cells x genes
adata = ad.AnnData(X=X,
                   obs=pd.DataFrame(index=pd.Index(cells, name="cell")),
                   var=pd.DataFrame(index=pd.Index(genes, name="gene")))
adata.obs["sample"] = [c.rsplit("_", 1)[0] for c in adata.obs_names]
print("AnnData:", adata.shape, "samples:", adata.obs['sample'].nunique(), flush=True)
print(adata.obs['sample'].value_counts().to_string(), flush=True)
adata.write(OUT)
print("wrote", OUT, "elapsed", round(time.time()-t0), "s")
