#!/usr/bin/env python3
"""Pull a few extra GSE207422 gene rows used to separate malignant-like
epithelium from normal lung epithelium (Hu et al. Genome Med 2023 markers).
Does not recompute total UMI.
"""
import gzip
from collections import defaultdict
import numpy as np
import pandas as pd

MATRIX = "/tmp/data_scrna_ici/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
OUT = "/tmp/data_scrna_ici/gse207422_extra_genes.tsv.gz"
WANTED = {
    "SFTPA2", "SFTPD", "AGER", "SCGB1A1", "SCGB3A2", "TPPP3", "FOXJ1",
    "CAPS", "DST", "SERPINB9", "PCNA", "KRT15", "KRT6A", "TP63", "SOX2",
    "NKX2-1", "NAPSA", "WFDC2", "CAPS",
}

dtypes = defaultdict(lambda: np.int32)
dtypes["Gene"] = str
kept = {}
reader = pd.read_csv(MATRIX, sep="\t", index_col=0, chunksize=2000, dtype=dtypes)
n = 0
for chunk in reader:
    n += chunk.shape[0]
    hits = chunk.index.intersection(WANTED)
    for g in hits:
        kept[g] = chunk.loc[g].to_numpy()
    print(f"scanned {n} genes, kept {len(kept)}", flush=True)

print("present:", sorted(kept))
print("missing:", sorted(WANTED - set(kept)))
# barcode column from the already-extracted panel file
barcodes = pd.read_csv("/tmp/data_scrna_ici/gse207422_panel_cells.tsv.gz",
                       sep="\t", usecols=["barcode"])["barcode"]
df = pd.DataFrame({"barcode": barcodes})
for g in sorted(kept):
    df[g] = kept[g]
df.to_csv(OUT, sep="\t", index=False, compression="gzip")
print("wrote", OUT)
