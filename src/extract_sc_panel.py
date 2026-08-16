"""Stream-extract the A11 gene panel from GSE131907 (genes x ~208k cells).

The full matrix is 2.9 GB gzipped; we never load it. One pass writes a
cells x genes parquet of only the panel genes.
"""

import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C

SRC = C.DATA / "sc" / "GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz"
ANN = C.DATA / "sc" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
OUT = C.DATA / "derived"
OUT.mkdir(parents=True, exist_ok=True)

# Also pull a few extra lineage genes used for sanity checks.
EXTRA = ["EPCAM", "KRT8", "PTPRC", "COL1A1", "PECAM1"]


def main():
    wanted = set(C.all_genes()) | set(EXTRA)
    # Apply aliases in reverse: older symbols that might appear as keys.
    alias_rev = {v: k for k, v in C.ALIASES.items()}
    also = {alias_rev[g] for g in wanted if g in alias_rev}
    wanted |= also

    print("streaming", SRC, flush=True)
    genes = {}
    barcodes = None
    n_seen = 0
    with gzip.open(SRC, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        barcodes = header[1:]
        print(f"n_cells={len(barcodes)}", flush=True)
        for line in fh:
            n_seen += 1
            tab = line.find("\t")
            gene = line[:tab]
            if gene in wanted:
                vals = np.fromstring(line[tab + 1:], sep="\t", dtype=np.float32)
                genes[gene] = vals
                print(f"  captured {gene} ({len(genes)}/{len(wanted)}) at row {n_seen}",
                      flush=True)
            if n_seen % 5000 == 0:
                print(f"  scanned {n_seen} genes, captured {len(genes)}", flush=True)

    print(f"done scan: {n_seen} genes, captured {list(genes)}", flush=True)
    mat = pd.DataFrame(genes, index=barcodes)
    mat.index.name = "Index"

    # Harmonise aliases to current symbols.
    rename = {old: new for old, new in C.ALIASES.items() if old in mat.columns}
    if rename:
        mat = mat.rename(columns=rename)
        mat = mat.loc[:, ~mat.columns.duplicated()]

    missing = sorted(set(C.all_genes()) - set(mat.columns))
    print("missing from GSE131907:", missing, flush=True)

    ann = pd.read_csv(ANN, sep="\t")
    # Annotation Index should match matrix barcodes.
    print("ann", ann.shape, "overlap", ann["Index"].isin(mat.index).sum(), flush=True)
    ann = ann.set_index("Index")
    common = mat.index.intersection(ann.index)
    mat = mat.loc[common]
    ann = ann.loc[common]
    print("aligned", mat.shape, ann.shape, flush=True)
    print(ann["Cell_type"].value_counts().to_string())
    print(ann["Sample_Origin"].value_counts().to_string())

    mat.to_parquet(OUT / "sc_GSE131907_panel.parquet")
    ann.to_parquet(OUT / "sc_GSE131907_annot.parquet")
    print("wrote", OUT / "sc_GSE131907_panel.parquet", mat.shape, flush=True)


if __name__ == "__main__":
    main()
