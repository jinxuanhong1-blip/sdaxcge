"""Pass 1 over the GSE154826 raw 10x matrices.

For every amplification batch this writes three things:

  cellstats_<batch>.parquet  per-barcode totals + marker/ADT counts for every
                             cell-like barcode (>=200 GEX UMI) and every barcode
                             that the authors annotated
  gate_<batch>.npz           full transcriptome of the annotated cells that fall in
                             the authors' "epi_endo_fibro_doublet" gate, i.e. the
                             pool that any epithelial cell has to come from
  pbulk_<batch>.npz          sample x cluster x gene count sums over annotated cells

The published cell->sample and cell->cluster assignments (Leader et al.
input_tables/cell_metadata.csv) are used verbatim; nothing is re-clustered here.
"""

import gzip
import io
import os
import subprocess
import sys
import tarfile
from multiprocessing import Pool

import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_gse154826 as L

SCRATCH = "/tmp/gse_scratch"
MIN_UMI_STATS = 200


def read_batch_files(batch):
    """Return (features_df, barcodes array, triplet arrays) for one amp batch."""
    with tarfile.open(L.batch_tar(batch), "r:gz") as tf:
        names = tf.getnames()
        fname = [n for n in names if n.endswith("features.tsv")][0]
        bname = [n for n in names if n.endswith("barcodes.tsv")][0]
        mname = [n for n in names if n.endswith("matrix.mtx")][0]
        feats = pd.read_csv(tf.extractfile(fname), sep="\t", header=None,
                            names=["ens", "symbol", "ftype"], dtype=str)
        barcodes = pd.read_csv(tf.extractfile(bname), header=None)[0].values
        fh = tf.extractfile(mname)
        head = []
        while True:
            pos = fh.tell()
            line = fh.readline().decode()
            if line.startswith("%"):
                head.append(line)
                continue
            dims = [int(x) for x in line.split()]
            break
        trip = pd.read_csv(fh, sep=" ", header=None, names=["row", "col", "val"],
                           dtype=np.int32, engine="c")
    assert dims[0] == len(feats) and dims[1] == len(barcodes), (batch, dims)
    return feats, barcodes, trip


def process_batch(batch):
    out_stats = os.path.join(L.BATCHDIR, f"cellstats_{batch}.parquet")
    if os.path.exists(out_stats):
        return batch, "cached"

    cm = pd.read_csv(os.path.join(L.PROC, "cell_metadata_with_batch.csv"))
    cm = cm[cm.amp_batch_ID == batch]

    feats, barcodes, trip = read_batch_files(batch)
    is_gex = (feats.ftype == "Gene Expression").values
    symbols = feats.symbol.values
    n_feat, n_cell = len(feats), len(barcodes)

    row = trip.row.values - 1
    col = trip.col.values - 1
    val = trip.val.values.astype(np.int64)
    del trip

    gex_mask = is_gex[row]
    gex_umi = np.bincount(col[gex_mask], weights=val[gex_mask], minlength=n_cell)
    gex_genes = np.bincount(col[gex_mask], minlength=n_cell)

    stats = pd.DataFrame({"barcode": barcodes, "gex_umi": gex_umi.astype(np.int64),
                          "n_genes": gex_genes.astype(np.int32)})

    # marker / target panel counts per barcode
    sym_to_rows = {}
    for i, s in enumerate(symbols):
        if is_gex[i]:
            sym_to_rows.setdefault(s, []).append(i)
    for g in L.PANEL_GENES:
        idx = sym_to_rows.get(g)
        if idx is None:
            stats[f"g_{g}"] = 0
            continue
        m = np.isin(row, idx)
        stats[f"g_{g}"] = np.bincount(col[m], weights=val[m],
                                      minlength=n_cell).astype(np.int32)
    mito_rows = [i for i, s in enumerate(symbols) if is_gex[i] and s.startswith("MT-")]
    m = np.isin(row, mito_rows)
    stats["mito_umi"] = np.bincount(col[m], weights=val[m], minlength=n_cell).astype(np.int32)

    # ADT (CITE-seq / hashing) features, when the batch has them
    adt_rows = np.where(~is_gex)[0]
    if len(adt_rows):
        m = ~gex_mask
        adt = sp.coo_matrix((val[m], (row[m], col[m])), shape=(n_feat, n_cell)).tocsr()
        for r in adt_rows:
            stats[f"adt_{symbols[r]}"] = np.asarray(adt[r].todense()).ravel().astype(np.int32)
        del adt

    stats["amp_batch_ID"] = batch
    ann = cm.set_index("barcode")
    stats["sample_ID"] = stats.barcode.map(ann.sample_ID)
    stats["cluster_ID"] = stats.barcode.map(ann.cluster_ID)
    stats["annotated"] = stats.sample_ID.notna()
    keep = (stats.gex_umi >= MIN_UMI_STATS) | stats.annotated
    stats = stats[keep].reset_index(drop=True)
    assert stats.annotated.sum() == len(cm), (batch, stats.annotated.sum(), len(cm))

    # ---- full transcriptome of annotated gate cells -------------------------
    bc_to_col = pd.Series(np.arange(n_cell), index=barcodes)
    gate_bc = cm.loc[cm.is_gate, "barcode"].values
    gate_cols = bc_to_col.loc[gate_bc].values if len(gate_bc) else np.array([], int)
    gene_rows = np.where(is_gex)[0]
    row_to_gene = np.full(n_feat, -1, np.int32)
    row_to_gene[gene_rows] = np.arange(len(gene_rows))

    if len(gate_cols):
        col_to_new = np.full(n_cell, -1, np.int64)
        col_to_new[gate_cols] = np.arange(len(gate_cols))
        m = gex_mask & (col_to_new[col] >= 0)
        X = sp.coo_matrix((val[m], (col_to_new[col[m]], row_to_gene[row[m]])),
                          shape=(len(gate_cols), len(gene_rows))).tocsr()
    else:
        X = sp.csr_matrix((0, len(gene_rows)))
    sub = cm[cm.is_gate]
    np.savez_compressed(os.path.join(L.BATCHDIR, f"gate_{batch}.npz"),
                        data=X.data.astype(np.int32), indices=X.indices,
                        indptr=X.indptr, shape=np.array(X.shape),
                        barcode=sub.barcode.values.astype(str),
                        sample_ID=sub.sample_ID.values,
                        cluster_ID=sub.cluster_ID.values)
    del X

    # ---- sample x cluster pseudobulk over annotated cells -------------------
    keys = list(zip(cm.sample_ID.values, cm.cluster_ID.values))
    uniq = sorted(set(keys))
    key_idx = {k: i for i, k in enumerate(uniq)}
    cell_group = np.full(n_cell, -1, np.int64)
    cell_group[bc_to_col.loc[cm.barcode.values].values] = [key_idx[k] for k in keys]
    m = gex_mask & (cell_group[col] >= 0)
    flat = cell_group[col[m]] * len(gene_rows) + row_to_gene[row[m]]
    pb = np.bincount(flat, weights=val[m],
                     minlength=len(uniq) * len(gene_rows)).reshape(len(uniq), len(gene_rows))
    ncells = np.bincount([key_idx[k] for k in keys], minlength=len(uniq))
    np.savez_compressed(os.path.join(L.BATCHDIR, f"pbulk_{batch}.npz"),
                        counts=pb.astype(np.int64),
                        sample_ID=np.array([k[0] for k in uniq]),
                        cluster_ID=np.array([k[1] for k in uniq]),
                        n_cells=ncells,
                        genes=symbols[gene_rows].astype(str))

    stats.to_parquet(out_stats, index=False)
    return batch, f"{len(stats)} rows, {int(stats.annotated.sum())} annotated"


def main():
    annots = L.load_annots_list()
    cm = L.load_cell_metadata()
    geo = L.load_sample_annots()
    cm = cm.merge(geo[["sample_ID", "amp_batch_ID"]], on="sample_ID", how="left")
    lin = annots.set_index("cluster")["lineage"]
    cm["lineage"] = cm.cluster_ID.map(lin)
    cm["is_gate"] = cm.lineage == L.GATE_LINEAGE
    assert cm.amp_batch_ID.notna().all()
    cm.to_csv(os.path.join(L.PROC, "cell_metadata_with_batch.csv"), index=False)

    batches = sorted(cm.amp_batch_ID.unique().astype(int))
    print(f"{len(batches)} batches with annotated cells", flush=True)
    with Pool(3) as p:
        for b, msg in p.imap_unordered(process_batch, batches):
            print(f"batch {b}: {msg}", flush=True)


if __name__ == "__main__":
    main()
