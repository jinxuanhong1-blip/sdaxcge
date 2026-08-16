"""Shared helpers for the USER-ALIGN TROP2 tight-junction analysis.

All outputs live under results/align_tj/.
"""
import os
import gzip
import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data")
TABLES = os.path.join(BASE, "tables")
FIGURES = os.path.join(BASE, "figures")
for d in (DATA, TABLES, FIGURES):
    os.makedirs(d, exist_ok=True)

# The user's private hypothesis gene sets (what we are trying to align public data to)
USER_TJ_GENES = ["CLDN1", "CLDN4", "CLDN7", "F11R", "PARD3"]
USER_TFS = ["ELF3", "GRHL1", "KLF4", "TFAP2A"]
TROP2 = "TACSTD2"

# Themes the user expects enriched in TROP2-high tumors
THEME_KEYWORDS = {
    "keratinization": ["keratin"],
    "skin_barrier": ["skin", "epiderm", "cornif"],
    "tight_junction": ["tight junction", "apical junction", "cell-cell junction",
                        "bicellular tight junction", "cell junction"],
    "emt": ["epithelial mesenchymal transition", "epithelial-mesenchymal"],
}


def read_xena_matrix(path):
    """Read a UCSC Xena gene x sample matrix (gene symbols in first column)."""
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index.name = "gene"
    # collapse duplicate gene rows by mean if any
    if df.index.duplicated().any():
        df = df.groupby(level=0).mean()
    return df


def stream_umi_columns(path, keep_barcodes=None, goi=None, genome_wide=False,
                       gene_col=0, progress_every=4000):
    """Stream a gene x cell TSV (gzip OK).

    keep_barcodes : if set, only those cell columns are extracted
    goi           : if set and genome_wide is False, only those gene rows kept
    Always returns libsize over ALL genes for the kept cells.

    Returns (raw_df of kept genes x kept cells, libsize Series).
    """
    import gzip
    import numpy as np
    import pandas as pd
    opener = gzip.open if str(path).endswith(".gz") else open
    keep_set = set(keep_barcodes) if keep_barcodes is not None else None
    goi = set(goi) if goi is not None else None
    with opener(path, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        if keep_set is None:
            idxs = list(range(1, len(hdr)))
        else:
            idxs = [i for i, b in enumerate(hdr) if i > 0 and b in keep_set]
        names = [hdr[i] for i in idxs]
        print(f"  stream: header_cells={len(hdr)-1} kept_cells={len(names)}")
        libsize = np.zeros(len(idxs), dtype=np.float64)
        genes, rows = [], []
        n = 0
        for line in fh:
            n += 1
            parts = line.rstrip("\n").split("\t")
            gene = parts[gene_col]
            vals = np.fromiter((parts[i] for i in idxs), dtype=np.float32, count=len(idxs))
            libsize += vals
            take = genome_wide or (goi is not None and gene in goi)
            if take:
                genes.append(gene)
                rows.append(vals)
            if progress_every and n % progress_every == 0:
                print(f"  ...rows {n}, kept genes {len(genes)}")
    if not rows:
        raw = pd.DataFrame(columns=names)
    else:
        raw = pd.DataFrame(np.vstack(rows), index=genes, columns=names)
        raw = raw[~raw.index.duplicated()]
    lib = pd.Series(libsize, index=names, name="libsize")
    return raw, lib


def load_sc_matrix(path, genes_of_interest, chunksize=1000, sep="\t"):
    """Stream a big gene x cell UMI text matrix (genes in rows, cells in cols).

    Returns (goi_raw, libsize):
      goi_raw : DataFrame (genes_of_interest present) x cells, raw counts
      libsize : Series per cell = total UMI across ALL genes
    Only the genes_of_interest rows are retained in memory; library size is
    accumulated over every gene so normalization is correct.
    """
    import pandas as pd
    goi = set(genes_of_interest)
    libsize = None
    kept = []
    reader = pd.read_csv(path, sep=sep, index_col=0, chunksize=chunksize)
    ncells = None
    for chunk in reader:
        chunk = chunk.apply(pd.to_numeric, errors="coerce").fillna(0)
        if ncells is None:
            ncells = chunk.shape[1]
        s = chunk.sum(axis=0)
        libsize = s if libsize is None else libsize.add(s, fill_value=0)
        hit = chunk.index.intersection(goi)
        if len(hit):
            kept.append(chunk.loc[hit])
    goi_raw = pd.concat(kept) if kept else pd.DataFrame()
    goi_raw = goi_raw[~goi_raw.index.duplicated()]
    return goi_raw, libsize


def lognorm(goi_raw, libsize, target=1e4):
    """CP10K + log1p normalization of the retained genes-of-interest rows."""
    import numpy as np
    sf = (libsize / target).replace(0, np.nan)
    norm = goi_raw.div(sf, axis=1)
    return np.log1p(norm)


def module_score(lognorm_df, genes):
    """Mean log-normalized expression across available genes of a set (per cell)."""
    import numpy as np
    g = [x for x in genes if x in lognorm_df.index]
    if not g:
        return None, 0
    return lognorm_df.loc[g].mean(axis=0), len(g)


def assign_lineage(lognorm_df, lineage_markers):
    """Assign each cell to the lineage with the highest mean-marker z-score."""
    import numpy as np
    import pandas as pd
    scores = {}
    for lin, markers in lineage_markers.items():
        g = [m for m in markers if m in lognorm_df.index]
        if not g:
            continue
        z = lognorm_df.loc[g]
        z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1).replace(0, np.nan), axis=0)
        scores[lin] = z.mean(axis=0)
    S = pd.DataFrame(scores)
    return S.idxmax(axis=1), S


def tumor_samples(cols):
    """Keep TCGA primary tumor samples (barcode sample-type code 01)."""
    keep = []
    for c in cols:
        parts = c.split("-")
        if len(parts) >= 4:
            code = parts[3][:2]
            if code == "01":  # primary solid tumor
                keep.append(c)
    return keep
