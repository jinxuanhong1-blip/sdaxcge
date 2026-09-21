#!/usr/bin/env python3
"""Extract malignant-cell counts for the concordant-4 cohorts.

GSE123902 and GSE189357: marker malignant
    (EPCAM|KRT8|KRT18|KRT19) > 0 and PTPRC == 0.
GSE131907: author Cell_subtype == "Malignant cells" on the locked 21 samples.
GSE205335 is exported by extract_gse205335.R (author lineage, non-normal).

One matrix per cohort: genes x cells, raw counts, plus meta and library size.
"""
from __future__ import annotations

import argparse
import gzip
import subprocess
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread


EPI = ("EPCAM", "KRT8", "KRT18", "KRT19")
ROOT = Path(__file__).resolve().parents[1]


def locked_units(dataset: str) -> list[str]:
    df = pd.read_csv(ROOT / "data" / "locked_units.tsv", sep="\t")
    sub = df.loc[df["dataset"] == dataset, "unit_id"].astype(str)
    return sub.tolist()


def collapse_genes(X: sparse.spmatrix, genes: np.ndarray) -> tuple[sparse.csr_matrix, np.ndarray]:
    genes = np.array([str(g).strip().upper() for g in genes], dtype=object)
    keep = np.array([g not in {"", "NAN", "NA", "NONE"} and not str(g).startswith("ENSG") for g in genes])
    # Keep Ensembl ids only when no symbol was available; caller should pass symbols.
    # Here ENSG-looking strings are dropped if a real symbol list was provided.
    # If EVERY name looks like ENSG, keep them (should not happen for these cohorts).
    if keep.sum() == 0:
        keep[:] = True
    if not keep.all():
        X = X[keep]
        genes = genes[keep]
    if len(set(genes.tolist())) == len(genes):
        return X.tocsr(), genes
    uniq, inv = np.unique(genes, return_inverse=True)
    coo = X.tocoo()
    Xc = sparse.csr_matrix((coo.data, (inv[coo.row], coo.col)), shape=(len(uniq), X.shape[1]))
    Xc.sum_duplicates()
    return Xc, uniq


def save_cohort(out: Path, X: sparse.csr_matrix, genes: np.ndarray, meta: pd.DataFrame, lib: np.ndarray) -> None:
    out.mkdir(parents=True, exist_ok=True)
    X = X.tocsr().astype(np.float32)
    X.sort_indices()
    sparse.save_npz(out / "counts.npz", X)
    pd.Series(genes).to_csv(out / "genes.tsv", index=False, header=False)
    meta.to_csv(out / "meta.tsv", sep="\t", index=False)
    np.save(out / "libsize.npy", np.asarray(lib, dtype=np.float64))
    print(f"saved {out} genes={X.shape[0]} cells={X.shape[1]} nnz={X.nnz}", flush=True)


def marker_mask(X_cells_by_genes: np.ndarray, genes: np.ndarray) -> np.ndarray:
    idx = {g: i for i, g in enumerate(genes)}
    missing = [g for g in list(EPI) + ["PTPRC"] if g not in idx]
    if missing:
        raise SystemExit(f"marker genes missing: {missing}")
    epi = np.column_stack([X_cells_by_genes[:, idx[g]] > 0 for g in EPI]).any(axis=1)
    return epi & (X_cells_by_genes[:, idx["PTPRC"]] == 0)


def extract_gse123902(raw: Path, out_root: Path) -> None:
    print("==== GSE123902 ====", flush=True)
    tar_path = raw / "GSE123902" / "GSE123902_RAW.tar"
    csv_dir = raw / "GSE123902" / "csv"
    csv_dir.mkdir(parents=True, exist_ok=True)
    if not list(csv_dir.glob("*_dense.csv.gz")):
        with tarfile.open(tar_path) as tf:
            tf.extractall(csv_dir)
    locked = locked_units("GSE123902")
    # Prefer PRIMARY over METASTASIS. Ignore NORMAL. Matches the locked donor table.
    files = sorted(csv_dir.glob("*_dense.csv.gz"))
    chosen: dict[str, Path] = {}
    for fp in files:
        name = fp.name
        if "_NORMAL_" in name:
            continue
        # GSM..._MSK_{patient}_{PRIMARY_TUMOUR|METASTASIS}_dense.csv.gz
        parts = name.split("_")
        # patient token after MSK
        try:
            i = parts.index("MSK")
            patient = parts[i + 1]
        except ValueError:
            continue
        if patient not in locked:
            continue
        rank = 0 if "PRIMARY" in name else 1
        prev = chosen.get(patient)
        if prev is None or (rank == 0 and "PRIMARY" not in prev.name):
            chosen[patient] = fp
    missing = [u for u in locked if u not in chosen]
    if missing:
        raise SystemExit(f"GSE123902 missing files for {missing}")

    # Dense CSVs do not share a gene universe. Test the intersection so a
    # missing symbol is not filled in as zero. Library size still uses every
    # gene measured in that file.
    stored = []
    for patient in locked:
        fp = chosen[patient]
        print(f"  read {fp.name}", flush=True)
        df = pd.read_csv(fp, index_col=0)
        genes = np.array([str(c).upper() for c in df.columns], dtype=object)
        barcodes = df.index.astype(str).to_numpy()
        mat = np.asarray(df.to_numpy(), dtype=np.float32)
        del df
        mal = marker_mask(mat, genes)
        sub = mat[mal]
        lib = sub.sum(axis=1).astype(np.float64)
        Xs, genes_c = collapse_genes(sparse.csr_matrix(sub.T), genes)
        stored.append((patient, barcodes[mal], Xs, genes_c, lib, int(mal.sum()), int(mat.shape[0])))
        print(f"    {patient} malignant {int(mal.sum())} / {mat.shape[0]} genes {len(genes_c)}", flush=True)
        del mat, sub
    common = set(stored[0][3].tolist())
    for rec in stored[1:]:
        common &= set(rec[3].tolist())
    common_genes = np.array(sorted(common), dtype=object)
    print(f"  intersection genes {len(common_genes)}", flush=True)
    blocks = []
    metas = []
    libs = []
    for patient, barcodes, Xs, genes_c, lib, _nmal, _n in stored:
        pos = {g: i for i, g in enumerate(genes_c.tolist())}
        rows = [pos[g] for g in common_genes]
        blocks.append(Xs[rows])
        metas.append(pd.DataFrame({
            "barcode": [f"{patient}_{b}" for b in barcodes],
            "unit_id": patient,
            "dataset": "GSE123902",
        }))
        libs.append(lib)
    X = sparse.hstack(blocks, format="csr")
    meta = pd.concat(metas, ignore_index=True)
    save_cohort(out_root / "GSE123902", X, common_genes, meta, np.concatenate(libs))


def extract_gse189357(raw: Path, out_root: Path) -> None:
    print("==== GSE189357 ====", flush=True)
    tar_path = raw / "GSE189357" / "GSE189357_RAW.tar"
    ex = raw / "GSE189357" / "raw"
    ex.mkdir(parents=True, exist_ok=True)
    if not list(ex.glob("*_matrix.mtx.gz")):
        with tarfile.open(tar_path) as tf:
            tf.extractall(ex)
    locked = locked_units("GSE189357")
    blocks = []
    metas = []
    libs = []
    gene_ref = None
    for patient in locked:
        mtx = ex / f"GSM5699777_{patient}_matrix.mtx.gz"
        # GSM ids differ; find by patient token
        hits = sorted(ex.glob(f"*_{patient}_matrix.mtx.gz"))
        if len(hits) != 1:
            raise SystemExit(f"expected one mtx for {patient}, got {hits}")
        mtx = hits[0]
        feat = Path(str(mtx).replace("_matrix.mtx.gz", "_features.tsv.gz"))
        bc = Path(str(mtx).replace("_matrix.mtx.gz", "_barcodes.tsv.gz"))
        print(f"  read {mtx.name}", flush=True)
        X = mmread(gzip.open(mtx, "rb")).tocsr().astype(np.float32)
        feats = pd.read_csv(feat, sep="\t", header=None)
        if feats.shape[1] >= 2:
            symbols = feats.iloc[:, 1].astype(str).to_numpy()
        else:
            symbols = feats.iloc[:, 0].astype(str).to_numpy()
        barcodes = pd.read_csv(bc, sep="\t", header=None).iloc[:, 0].astype(str).to_numpy()
        X, symbols = collapse_genes(X, symbols)
        if gene_ref is None:
            gene_ref = symbols
            X_use = X
        elif np.array_equal(symbols, gene_ref):
            X_use = X
        else:
            idx = {g: i for i, g in enumerate(symbols.tolist())}
            rows = []
            data_rows = []
            # reindex to gene_ref; genes only in this sample are ignored until union pass
            # Build union incrementally: if new genes appear, expand gene_ref and past blocks.
            new = [g for g in symbols.tolist() if g not in set(gene_ref.tolist())]
            if new:
                gene_ref = np.concatenate([gene_ref, np.array(new, dtype=object)])
                blocks = [sparse.vstack([b, sparse.csr_matrix((len(new), b.shape[1]), dtype=np.float32)], format="csr") for b in blocks]
            pos = {g: i for i, g in enumerate(gene_ref.tolist())}
            coo = X.tocoo()
            row2 = np.fromiter((pos[symbols[r]] for r in coo.row), dtype=np.int32, count=coo.nnz)
            X_use = sparse.csr_matrix((coo.data, (row2, coo.col)), shape=(len(gene_ref), X.shape[1]))
        # marker gate on this sample's own matrix (symbols may have been collapsed)
        # Rebuild a cells x genes dense slice for the four markers only.
        sym_to_i = {g: i for i, g in enumerate(symbols.tolist())}
        # Use X (collapsed, this sample order) for the gate, then columns of X_use match cells.
        missing = [g for g in list(EPI) + ["PTPRC"] if g not in sym_to_i]
        if missing:
            raise SystemExit(f"{patient} missing {missing}")
        epi = np.zeros(X.shape[1], dtype=bool)
        for g in EPI:
            epi |= np.asarray(X[sym_to_i[g]].todense()).ravel() > 0
        pt = np.asarray(X[sym_to_i["PTPRC"]].todense()).ravel() == 0
        mal = epi & pt
        X_use = X_use[:, mal]
        lib = np.asarray(X[:, mal].sum(axis=0)).ravel().astype(np.float64)
        blocks.append(X_use.tocsr())
        metas.append(pd.DataFrame({
            "barcode": [f"{patient}_{b}" for b in barcodes[mal]],
            "unit_id": patient,
            "dataset": "GSE189357",
        }))
        libs.append(lib)
        print(f"    {patient} malignant {int(mal.sum())} / {len(barcodes)}", flush=True)
    # pad early blocks if gene_ref grew after them — already padded in-loop
    X = sparse.hstack(blocks, format="csr")
    meta = pd.concat(metas, ignore_index=True)
    lib = np.concatenate(libs)
    save_cohort(out_root / "GSE189357", X, gene_ref, meta, lib)


def extract_gse131907(raw: Path, out_root: Path) -> None:
    print("==== GSE131907 ====", flush=True)
    ann = raw / "GSE131907" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    mat = raw / "GSE131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    locked = set(locked_units("GSE131907"))
    print("  index header", flush=True)
    with gzip.open(mat, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
    cells = header[1:]
    col_of = {b: i for i, b in enumerate(cells)}
    keep_cols = []
    keep_units = []
    keep_bc = []
    n_ann = 0
    n_miss = 0
    with gzip.open(ann, "rt") as f:
        ah = f.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(ah)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            sample = p[ix["Sample"]]
            if sample not in locked:
                continue
            if p[ix["Cell_subtype"]] != "Malignant cells":
                continue
            n_ann += 1
            index = p[ix["Index"]]
            c = col_of.get(index)
            if c is None:
                n_miss += 1
                continue
            keep_cols.append(c)
            keep_units.append(sample)
            keep_bc.append(index)
    if n_miss:
        raise SystemExit(f"GSE131907 malignant barcodes missing from matrix: {n_miss}")
    order = np.argsort(np.asarray(keep_cols), kind="mergesort")
    keep_cols = [keep_cols[i] for i in order]
    keep_units = [keep_units[i] for i in order]
    keep_bc = [keep_bc[i] for i in order]
    print(f"  malignant cells {len(keep_cols)}", flush=True)
    out = out_root / "GSE131907"
    out.mkdir(parents=True, exist_ok=True)
    keep_path = out / "keep_cols.txt"
    keep_path.write_text("\n".join(str(c) for c in keep_cols) + "\n")
    cc = ROOT / "scripts" / "extract_subset.c"
    bin_path = out / "extract_subset"
    subprocess.check_call(["gcc", "-O3", "-o", str(bin_path), str(cc)])
    prefix = str(out / "stream")
    print("  stream matrix", flush=True)
    # gzip.GzipFile.fileno() is the compressed fd, so spawn gzip -dc instead.
    with open(prefix + ".log", "w") as log:
        src = subprocess.Popen(["gzip", "-dc", str(mat)], stdout=subprocess.PIPE)
        try:
            subprocess.check_call([str(bin_path), str(keep_path), prefix], stdin=src.stdout, stderr=log)
        finally:
            if src.stdout:
                src.stdout.close()
            rc = src.wait()
        if rc != 0:
            raise SystemExit(f"gzip -dc failed status {rc}")
    print(Path(prefix + ".log").read_text()[-500:], flush=True)
    shape = Path(prefix + ".shape.txt").read_text().strip().split("\t")
    n_genes, n_kept, nnz = int(shape[0]), int(shape[1]), int(shape[2])
    genes = pd.read_csv(prefix + ".genes.txt", header=None, dtype=str)[0].to_numpy()
    if len(genes) != n_genes:
        raise SystemExit(f"gene count {len(genes)} != {n_genes}")
    if n_kept != len(keep_cols):
        raise SystemExit(f"kept {n_kept} != {len(keep_cols)}")
    rec = np.fromfile(prefix + ".coo.bin", dtype=np.dtype([("gi", "<i4"), ("cj", "<i4"), ("v", "<f4")]))
    if len(rec) != nnz:
        raise SystemExit(f"nnz file {len(rec)} != header {nnz}")
    X = sparse.csr_matrix((rec["v"], (rec["gi"], rec["cj"])), shape=(n_genes, n_kept))
    del rec
    lib = np.fromfile(prefix + ".libsize.f64", dtype="<f8")
    meta = pd.DataFrame({"barcode": keep_bc, "unit_id": keep_units, "dataset": "GSE131907"})
    X, genes = collapse_genes(X, genes)
    # collapse may reorder; lib is per cell and stays aligned
    save_cohort(out, X, genes, meta, lib)


def compare_locked(out_root: Path) -> None:
    locked = pd.read_csv(ROOT / "data" / "locked_units.tsv", sep="\t")
    rows = []
    for ds in ["GSE123902", "GSE131907", "GSE189357", "GSE205335"]:
        meta_path = out_root / ds / "meta.tsv"
        if not meta_path.exists():
            print(f"missing {ds}", flush=True)
            continue
        meta = pd.read_csv(meta_path, sep="\t")
        got = meta.groupby("unit_id").size()
        exp = locked.loc[locked["dataset"] == ds].set_index("unit_id")["n_malignant"]
        for uid, e in exp.items():
            g = int(got.get(uid, 0))
            rows.append((ds, uid, int(e), g, g - int(e)))
    df = pd.DataFrame(rows, columns=["dataset", "unit_id", "locked_n_malignant", "export_n", "delta"])
    df.to_csv(out_root / "malignant_count_check.tsv", sep="\t", index=False)
    bad = df.loc[df["delta"] != 0]
    print(df.groupby("dataset")[["locked_n_malignant", "export_n"]].sum().to_string(), flush=True)
    if len(bad):
        print("COUNT MISMATCHES", flush=True)
        print(bad.head(30).to_string(index=False), flush=True)
    else:
        print("malignant counts match locked units", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default="/tmp/concordant4_raw")
    ap.add_argument("--out", default="/tmp/concordant4_malig")
    ap.add_argument("--cohort", default="all", choices=["all", "GSE123902", "GSE131907", "GSE189357", "check"])
    args = ap.parse_args()
    raw, out = Path(args.raw), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    if args.cohort in {"all", "GSE123902"}:
        extract_gse123902(raw, out)
    if args.cohort in {"all", "GSE189357"}:
        extract_gse189357(raw, out)
    if args.cohort in {"all", "GSE131907"}:
        extract_gse131907(raw, out)
    if args.cohort in {"all", "check"}:
        compare_locked(out)


if __name__ == "__main__":
    main()
