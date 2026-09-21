#!/usr/bin/env python3
"""Build a malignant-only AnnData for the locked concordant-4.

Gates match the locked patient table (PR #539 / #503), not a new malignant call:

- GSE123902: tumor/metastasis files only (NORMAL dropped). Marker malignant =
  (EPCAM|KRT8|KRT18|KRT19)>0 and PTPRC==0. Unit = donor.
- GSE131907: Cell_subtype == "Malignant cells" only. tS1/tS2/tS3 and nLung AT2
  are not in that locked gate and are not added. Unit = sample.
- GSE205335: lineage.sub == "Malignant cells" on non-normal tissue. Unit = patient.
- GSE189357: same marker gate. Unit = patient TD1–TD9.

This object is a cell-state matrix. It does not encode ICI time or RECIST.
"""
from __future__ import annotations

import argparse
import gc
import gzip
import io
import json
import re
import shutil
import tarfile
import tempfile
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import io as spio
from scipy import sparse

CAP = 180
SEED = 1
MIN_GENES = 200
MIN_COUNTS = 500
MAX_MITO = 20.0
EPI = ("EPCAM", "KRT8", "KRT18", "KRT19")


def _upper_unique(genes: list[str]) -> list[str]:
    seen: dict[str, int] = {}
    out = []
    for g in genes:
        u = str(g).upper().strip()
        if u in seen:
            seen[u] += 1
            u = f"{u}__{seen[u]}"
        else:
            seen[u] = 0
        out.append(u)
    return out


def _qc_cap(adata: ad.AnnData, unit_col: str, cap: int, seed: int) -> ad.AnnData:
    X = adata.X
    if sparse.issparse(X):
        counts = np.asarray(X.sum(axis=1)).ravel().astype(np.float64)
        n_genes = np.asarray((X > 0).sum(axis=1)).ravel()
    else:
        counts = np.asarray(X.sum(axis=1)).ravel().astype(np.float64)
        n_genes = np.asarray((X > 0).sum(axis=1)).ravel()
    mito = np.array([g.startswith("MT-") or g.startswith("MT.") for g in adata.var_names])
    if mito.any():
        mx = X[:, mito]
        msum = np.asarray(mx.sum(axis=1)).ravel().astype(np.float64)
        pct = np.divide(msum, counts, out=np.zeros(msum.shape, dtype=np.float64), where=counts > 0) * 100
    else:
        pct = np.zeros(adata.n_obs, dtype=np.float64)
    keep = (n_genes >= MIN_GENES) & (counts >= MIN_COUNTS) & (pct < MAX_MITO)
    print(
        f"  QC keep {int(keep.sum())}/{adata.n_obs} "
        f"(genes>={MIN_GENES}, counts>={MIN_COUNTS}, mito<{MAX_MITO})",
        flush=True,
    )
    adata = adata[keep].copy()
    adata.obs["n_counts"] = counts[keep]
    adata.obs["n_genes"] = n_genes[keep]
    rng = np.random.default_rng(seed)
    idx = []
    for unit, sub in adata.obs.groupby(unit_col, observed=True):
        ids = sub.index.to_numpy()
        if len(ids) > cap:
            ids = rng.choice(ids, size=cap, replace=False)
        idx.extend(ids.tolist())
    adata = adata[idx].copy()
    print(f"  after cap≤{cap}: {adata.n_obs} cells, {adata.obs[unit_col].nunique()} units", flush=True)
    return adata


def _finalize(adata: ad.AnnData, dataset: str, unit_col: str, malig_def: str) -> ad.AnnData:
    adata.obs["dataset"] = dataset
    adata.obs["unit_id"] = adata.obs[unit_col].astype(str)
    adata.obs["malig_def"] = malig_def
    adata.obs_names = [f"{dataset}:{c}" for c in adata.obs_names.astype(str)]
    adata.var_names = _upper_unique(list(adata.var_names))
    adata.var_names_make_unique()
    if not sparse.issparse(adata.X):
        adata.X = sparse.csr_matrix(adata.X.astype(np.float32))
    else:
        adata.X = adata.X.tocsr().astype(np.float32)
    adata.layers["counts"] = adata.X.copy()
    return adata


def _marker_mask(X: sparse.spmatrix, genes: list[str]) -> np.ndarray:
    idx = {g: i for i, g in enumerate(genes)}

    def col(g: str) -> np.ndarray:
        if g not in idx:
            return np.zeros(X.shape[1] if X.shape[0] == len(genes) else X.shape[0])
        # X is genes x cells
        v = X[idx[g]]
        return np.asarray(v.toarray()).ravel() if sparse.issparse(v) else np.asarray(v).ravel()

    epi = np.zeros(X.shape[1], dtype=bool)
    for g in EPI:
        epi |= col(g) > 0
    pt = col("PTPRC")
    return epi & (pt == 0)


def extract_gse123902(root: Path, cap: int) -> ad.AnnData:  # noqa: F811
    tar_path = root / "gse123902" / "GSE123902_RAW.tar"
    adatas = []
    with tarfile.open(tar_path) as tf:
        names = [m.name for m in tf.getmembers() if m.name.endswith(".csv.gz")]
        tumor = [n for n in names if "NORMAL" not in n.upper()]
        print(f"GSE123902 tumor files {len(tumor)} / {len(names)}", flush=True)
        for name in sorted(tumor):
            m = re.search(r"MSK_(LX[0-9A-Z]+)_(PRIMARY_TUMOUR|METASTASIS)", name)
            if not m:
                raise SystemExit(f"unparsed GSE123902 name {name}")
            donor, tissue = m.group(1), m.group(2)
            raw = tf.extractfile(name).read()
            df = pd.read_csv(io.BytesIO(gzip.decompress(raw)), index_col=0)
            df.columns = [str(c).upper() for c in df.columns]
            if df.columns.duplicated().any():
                df = df.T.groupby(level=0).max().T
            mat = sparse.csr_matrix(df.to_numpy(dtype=np.float32)).T
            mask = _marker_mask(mat, list(df.columns))
            n_mal = int(mask.sum())
            print(f"  {donor} {tissue}: cells={df.shape[0]} marker-malignant={n_mal}", flush=True)
            if n_mal == 0:
                continue
            sub = mat[:, mask].T.tocsr()
            obs = pd.DataFrame(
                {
                    "donor": donor,
                    "tissue": "PRIMARY" if "PRIMARY" in tissue else "METASTASIS",
                    "author_subtype": "marker_malignant",
                },
                index=pd.Index([f"{donor}_{j}" for j in range(n_mal)], name="cell"),
            )
            a = ad.AnnData(X=sub, obs=obs, var=pd.DataFrame(index=pd.Index(df.columns, name="gene")))
            adatas.append(a)
            del df, mat, sub
            gc.collect()
    if not adatas:
        raise SystemExit("GSE123902: no marker-malignant cells")
    shared = adatas[0].var_names
    for a in adatas[1:]:
        shared = shared.intersection(a.var_names)
    adatas = [a[:, shared].copy() for a in adatas]
    cat = ad.concat(adatas, join="inner", merge="same")
    cat = _qc_cap(cat, "donor", cap, SEED)
    return _finalize(cat, "GSE123902", "donor", "marker_malig")


def _stream_umi(umi_path: Path, keep_ids: list[str]) -> tuple[list[str], list[str], sparse.csr_matrix]:
    keep = set(keep_ids)
    with gzip.open(umi_path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        col_idx = [i for i, c in enumerate(cell_ids) if c in keep]
        if not col_idx:
            raise SystemExit("no GSE131907 malignant IDs in UMI header")
        ordered = [cell_ids[i] for i in col_idx]
        print(f"GSE131907 stream keep={len(col_idx)} / {len(cell_ids)}", flush=True)
        genes: list[str] = []
        data: list[np.ndarray] = []
        indices: list[np.ndarray] = []
        indptr = [0]
        nnz = 0
        for gi, line in enumerate(f, start=1):
            raw = line.rstrip("\n")
            if not raw:
                continue
            tab = raw.find("\t")
            gene = raw[:tab].upper()
            vals = np.fromstring(raw[tab + 1 :], sep="\t", dtype=np.float32)
            if vals.size != len(cell_ids):
                raise SystemExit(f"row {gi} {gene}: {vals.size} != {len(cell_ids)}")
            sub = vals[col_idx]
            nz = np.flatnonzero(sub)
            if nz.size:
                data.append(sub[nz])
                indices.append(nz.astype(np.int32))
                nnz += int(nz.size)
            indptr.append(nnz)
            genes.append(gene)
            if gi % 4000 == 0:
                print(f"  streamed {gi} genes nnz={nnz}", flush=True)
    data_a = np.concatenate(data) if data else np.array([], dtype=np.float32)
    indices_a = np.concatenate(indices) if indices else np.array([], dtype=np.int32)
    mat = sparse.csr_matrix(
        (data_a, indices_a, np.asarray(indptr, dtype=np.int64)),
        shape=(len(genes), len(ordered)),
        dtype=np.float32,
    )
    return ordered, genes, mat


def extract_gse131907(root: Path, cap: int) -> ad.AnnData:
    ann = pd.read_csv(
        root / "gse131907" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz",
        sep="\t",
        compression="gzip",
    )
    mal = ann.loc[ann["Cell_subtype"].astype(str).eq("Malignant cells")].copy()
    print(
        f"GSE131907 author malignant cells={len(mal)} samples={mal['Sample'].nunique()} "
        f"(tS1/tS2/tS3 not included)",
        flush=True,
    )
    cells, genes, mat = _stream_umi(
        root / "gse131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
        mal["Index"].astype(str).tolist(),
    )
    obs = mal.set_index("Index").reindex(cells)
    adata = ad.AnnData(
        X=mat.T.tocsr(),
        obs=pd.DataFrame(
            {
                "sample": obs["Sample"].astype(str).to_numpy(),
                "tissue": obs["Sample_Origin"].astype(str).to_numpy(),
                "author_subtype": obs["Cell_subtype"].astype(str).to_numpy(),
            },
            index=pd.Index(cells, name="cell"),
        ),
        var=pd.DataFrame(index=pd.Index(genes, name="gene")),
    )
    del mat
    gc.collect()
    adata = _qc_cap(adata, "sample", cap, SEED + 1)
    return _finalize(adata, "GSE131907", "sample", "author_malig")


def _parse_soft(path: Path) -> pd.DataFrame:
    opener = gzip.open if str(path).endswith(".gz") else open
    records = []
    current = None
    with opener(path, "rt", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    records.append(current)
                current = {"gsm": line.split(" = ", 1)[1]}
            elif current is not None and line.startswith("!Sample_description = "):
                current.setdefault("description", line.split(" = ", 1)[1])
            elif current is not None and line.startswith("!Sample_characteristics_ch1 = "):
                value = line.split(" = ", 1)[1]
                if ": " in value:
                    k, v = value.split(": ", 1)
                    current[k] = v
        if current is not None:
            records.append(current)
    meta = pd.DataFrame(records)
    read_end = meta["platform"].str.extract(r"Single Cell ([35])'")[0]
    meta["orig.ident"] = meta["description"].str.replace("_", "-", regex=False) + "-" + read_end + "P"
    return meta


def _load_rds_matrix(path: Path):
    import rdata

    with tempfile.TemporaryDirectory(prefix="gse205335-") as tmp:
        layer1 = Path(tmp) / "layer1"
        print("peel GSE205335 double-gzip", flush=True)
        with gzip.open(path, "rb") as src, layer1.open("wb") as dest:
            shutil.copyfileobj(src, dest, 16 * 1024 * 1024)
        magic = layer1.read_bytes()[:2]
        matrix_path = layer1
        if magic == b"\x1f\x8b":
            layer2 = Path(tmp) / "matrix.rds"
            with gzip.open(layer1, "rb") as src, layer2.open("wb") as dest:
                shutil.copyfileobj(src, dest, 16 * 1024 * 1024)
            matrix_path = layer2
            layer1.unlink()
        print(f"read RDS {matrix_path.stat().st_size} bytes", flush=True)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Missing constructor")
            obj = rdata.read_rds(matrix_path)
    # rdata may return a dict-like or an object with slots
    if isinstance(obj, dict):
        slots = obj
    else:
        slots = {k: getattr(obj, k) for k in ("i", "p", "Dim", "Dimnames", "x") if hasattr(obj, k)}
        if "x" not in slots and hasattr(obj, "value"):
            slots = obj.value
    expected = {"i", "p", "Dim", "Dimnames", "x"}
    if not expected.issubset(set(slots)):
        raise TypeError(f"RDS slots {sorted(slots)[:12]} not a dgCMatrix")
    genes = np.asarray(slots["Dimnames"][0], dtype=str)
    barcodes = np.asarray(slots["Dimnames"][1], dtype=str)
    dim = tuple(int(x) for x in slots["Dim"])
    print(f"build CSC {dim}", flush=True)
    matrix = sparse.csc_matrix((slots["x"], slots["i"], slots["p"]), shape=dim, copy=False)
    return matrix, genes, barcodes


def extract_gse205335(root: Path, cap: int) -> ad.AnnData:
    ident = pd.read_csv(root / "gse205335" / "GSE205335_Lung_IO_CellIdentity.txt.gz", sep="\t")
    meta = _parse_soft(root / "gse205335" / "GSE205335_family.soft.gz")
    cells = ident.merge(
        meta[["orig.ident", "patient", "tissue", "cancer subtype"]],
        on="orig.ident",
        how="left",
        validate="many_to_one",
    )
    if cells["patient"].isna().any():
        missing = cells.loc[cells["patient"].isna(), "orig.ident"].unique().tolist()
        raise SystemExit(f"GSE205335 samples missing patient map: {missing[:8]}")
    normal = cells["tissue"].astype(str).str.startswith("Normal")
    mal = cells.loc[cells["lineage.sub"].eq("Malignant cells") & ~normal].copy()
    print(
        f"GSE205335 malignant non-normal={len(mal)} patients={mal['patient'].nunique()} "
        f"normal-epi-like dropped={(cells['lineage.sub'].eq('Malignant cells') & normal).sum()}",
        flush=True,
    )
    matrix, genes, barcodes = _load_rds_matrix(root / "gse205335" / "GSE205335_Lung_IO_UMI_matrix.rds.gz")
    keep = set(mal["barcode"].astype(str))
    col_idx = np.array([i for i, b in enumerate(barcodes) if str(b) in keep], dtype=np.int64)
    print(f"GSE205335 subset {len(col_idx)} / {len(barcodes)}", flush=True)
    sub = matrix[:, col_idx].tocsr()
    del matrix
    gc.collect()
    ordered = [str(barcodes[i]) for i in col_idx]
    obs = mal.set_index("barcode").reindex(ordered)
    adata = ad.AnnData(
        X=sub.T.tocsr(),
        obs=pd.DataFrame(
            {
                "patient": obs["patient"].astype(str).to_numpy(),
                "tissue": obs["tissue"].astype(str).to_numpy(),
                "author_subtype": obs["lineage.sub"].astype(str).to_numpy(),
                "histology": obs["cancer subtype"].astype(str).to_numpy(),
            },
            index=pd.Index(ordered, name="cell"),
        ),
        var=pd.DataFrame(index=pd.Index([g.upper() for g in genes], name="gene")),
    )
    del sub
    gc.collect()
    adata = _qc_cap(adata, "patient", cap, SEED + 2)
    return _finalize(adata, "GSE205335", "patient", "author_malig")


def _read_10x_member(tf: tarfile.TarFile, sample: str) -> tuple[list[str], list[str], sparse.csr_matrix]:
    members = {m.name: m for m in tf.getmembers()}

    def find(token: str) -> tarfile.TarInfo:
        hits = [m for n, m in members.items() if sample in n and n.endswith(token)]
        if len(hits) != 1:
            raise SystemExit(f"{sample} {token}: {len(hits)} hits")
        return hits[0]

    feat_raw = gzip.decompress(tf.extractfile(find("features.tsv.gz")).read()).decode()
    genes = []
    for line in feat_raw.splitlines():
        parts = line.split("\t")
        genes.append((parts[1] if len(parts) > 1 else parts[0]).upper())
    bc_raw = gzip.decompress(tf.extractfile(find("barcodes.tsv.gz")).read()).decode().splitlines()
    mtx_bytes = gzip.decompress(tf.extractfile(find("matrix.mtx.gz")).read())
    mat = spio.mmread(io.BytesIO(mtx_bytes)).tocsr()
    if mat.shape[0] != len(genes):
        mat = mat.T.tocsr()
    if mat.shape[0] != len(genes) or mat.shape[1] != len(bc_raw):
        raise SystemExit(f"{sample} matrix {mat.shape} genes {len(genes)} bc {len(bc_raw)}")
    return genes, bc_raw, mat


def extract_gse189357(root: Path, cap: int) -> ad.AnnData:
    tar_path = root / "gse189357" / "GSE189357_RAW.tar"
    adatas = []
    with tarfile.open(tar_path) as tf:
        samples = sorted({re.search(r"(TD\d+)", m.name).group(1) for m in tf.getmembers() if "TD" in m.name})
        print(f"GSE189357 samples {samples}", flush=True)
        for sample in samples:
            genes, bcs, mat = _read_10x_member(tf, sample)
            # collapse duplicate gene symbols
            gser = pd.Series(np.arange(len(genes)), index=genes)
            if gser.index.duplicated().any():
                # sum duplicates via sparse grouping is heavier; max-keep first
                keep_g = ~pd.Index(genes).duplicated(keep="first")
                mat = mat[keep_g]
                genes = [g for g, k in zip(genes, keep_g) if k]
            mask = _marker_mask(mat, genes)
            n_mal = int(mask.sum())
            print(f"  {sample}: cells={mat.shape[1]} marker-malignant={n_mal}", flush=True)
            if n_mal == 0:
                continue
            sub = mat[:, mask].T.tocsr()
            obs = pd.DataFrame(
                {
                    "patient": sample,
                    "tissue": "tumor",
                    "author_subtype": "marker_malignant",
                },
                index=pd.Index([f"{sample}_{b}" for b, k in zip(bcs, mask) if k], name="cell"),
            )
            adatas.append(
                ad.AnnData(X=sub, obs=obs, var=pd.DataFrame(index=pd.Index(genes, name="gene")))
            )
            del mat, sub
            gc.collect()
    if not adatas:
        raise SystemExit("GSE189357: no marker-malignant cells")
    shared = adatas[0].var_names
    for a in adatas[1:]:
        shared = shared.intersection(a.var_names)
    adatas = [a[:, shared].copy() for a in adatas]
    cat = ad.concat(adatas, join="inner", merge="same")
    cat = _qc_cap(cat, "patient", cap, SEED + 3)
    return _finalize(cat, "GSE189357", "patient", "marker_malig")


def _sanitize(obs: pd.DataFrame) -> pd.DataFrame:
    out = obs.copy()
    for col in out.columns:
        if pd.api.types.is_bool_dtype(out[col]):
            out[col] = out[col].map({True: "True", False: "False"}).astype(str)
        elif out[col].dtype == object or str(out[col].dtype) == "category":
            out[col] = out[col].astype(str)
    return out


def concat_four(pieces: list[ad.AnnData]) -> ad.AnnData:
    shared = pieces[0].var_names
    for a in pieces[1:]:
        shared = shared.intersection(a.var_names)
    if "CLDN4" not in shared:
        raise SystemExit("CLDN4 missing from shared genes")
    if len(shared) < 8000:
        raise SystemExit(f"too few shared genes: {len(shared)}")
    cols = ["dataset", "unit_id", "malig_def", "tissue", "author_subtype"]
    kept = []
    for a in pieces:
        b = a[:, shared].copy()
        for c in cols:
            if c not in b.obs.columns:
                b.obs[c] = "NA"
        b.obs = _sanitize(b.obs[cols])
        kept.append(b)
    out = ad.concat(kept, join="inner", merge="same")
    out.layers["counts"] = out.X.copy()
    out.obs["unit_key"] = out.obs["dataset"].astype(str) + ":" + out.obs["unit_id"].astype(str)
    return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=Path, default=Path("/tmp/c4raw"))
    p.add_argument("--out", type=Path, default=Path("/tmp/c4_malignant.h5ad"))
    p.add_argument("--cap", type=int, default=CAP)
    args = p.parse_args()
    cache = args.data / "cache"
    cache.mkdir(parents=True, exist_ok=True)
    builders = [
        ("GSE123902", extract_gse123902),
        ("GSE131907", extract_gse131907),
        ("GSE205335", extract_gse205335),
        ("GSE189357", extract_gse189357),
    ]
    pieces = []
    for name, fn in builders:
        path = cache / f"{name}.h5ad"
        if path.is_file():
            print(f"reuse {path}", flush=True)
            pieces.append(ad.read_h5ad(path))
            continue
        print(f"=== extract {name} ===", flush=True)
        a = fn(args.data, args.cap)
        a.write_h5ad(path)
        print(f"cached {path} cells={a.n_obs} genes={a.n_vars}", flush=True)
        pieces.append(a)
        gc.collect()
    joint = concat_four(pieces)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    joint.write_h5ad(args.out)
    inv = {
        "n_cells": int(joint.n_obs),
        "n_genes": int(joint.n_vars),
        "n_units": int(joint.obs["unit_key"].nunique()),
        "by_dataset_cells": joint.obs["dataset"].value_counts().to_dict(),
        "by_dataset_units": joint.obs.groupby("dataset")["unit_id"].nunique().to_dict(),
        "cap": args.cap,
        "malig_def": joint.obs.groupby("dataset")["malig_def"].first().to_dict(),
        "not_an_ici_clock": True,
        "out": str(args.out),
    }
    args.out.with_suffix(".inventory.json").write_text(json.dumps(inv, indent=2))
    print(json.dumps(inv, indent=2), flush=True)


if __name__ == "__main__":
    main()
