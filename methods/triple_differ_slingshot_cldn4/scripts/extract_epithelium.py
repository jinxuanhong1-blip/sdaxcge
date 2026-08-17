#!/usr/bin/env python3
"""Build a joint epithelial AnnData for the triple that differs.

GSE123902 + GSE131907 + GSE205335. Do not add GSE148071. Not the 7-pool.

GSE131907: author Epithelial cells on tLung / nLung / tL/B / mLN / mBrain
(PE unlabeled epithelium is dropped). nLung AT2 is kept for the DPT root.
GSE205335: author lineage.total == Epithelial cells on non-normal tissues.
GSE123902: no public author cell-type table (36.5 GB H5 skipped). Marker
epithelium = (EPCAM|KRT8|KRT18|KRT19) > 0. NORMAL samples kept as units
for an AT2-like root. PRIMARY preferred over METASTASIS is NOT applied —
each sample is its own unit (like tLung vs nLung).
"""
from __future__ import annotations

import argparse
import gzip
import io
import json
import shutil
import sys
import tarfile
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from gene_sets import EPI_MARKERS  # noqa: E402

KEEP_ORIGINS_131907 = ("tLung", "nLung", "tL/B", "mLN", "mBrain")
EPI_TYPE_131907 = "Epithelial cells"
NORMAL_TISSUE_PREFIX = ("Normal ",)
CAP_PER_UNIT = 350
SEED = 1
COHORTS = ("GSE123902", "GSE131907", "GSE205335")


def parse_geo_soft(path: Path) -> pd.DataFrame:
    opener = gzip.open if str(path).endswith(".gz") else open
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    descriptions: list[str] = []
    titles: list[str] = []
    with opener(path, "rt", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    current["description"] = descriptions[0] if descriptions else ""
                    current["title"] = titles[0] if titles else ""
                    records.append(current)
                current = {"gsm": line.split(" = ", 1)[1]}
                descriptions, titles = [], []
            elif current is not None and line.startswith("!Sample_title = "):
                titles.append(line.split(" = ", 1)[1])
            elif current is not None and line.startswith("!Sample_description = "):
                descriptions.append(line.split(" = ", 1)[1])
            elif current is not None and line.startswith("!Sample_characteristics_ch1 = "):
                value = line.split(" = ", 1)[1]
                if ": " in value:
                    key, item = value.split(": ", 1)
                    current[key] = item
        if current is not None:
            current["description"] = descriptions[0] if descriptions else ""
            current["title"] = titles[0] if titles else ""
            records.append(current)
    metadata = pd.DataFrame(records)
    read_end = metadata["platform"].str.extract(r"Single Cell ([35])'")[0]
    metadata["orig.ident"] = (
        metadata["description"].str.replace("_", "-", regex=False) + "-" + read_end + "P"
    )
    return metadata.rename(
        columns={"tumor stage": "tumor_stage", "cancer subtype": "cancer_subtype"}
    )


def stream_umi_subset(
    umi_path: Path,
    keep_ids: list[str],
    gene_cap: int | None = None,
) -> tuple[list[str], list[str], sparse.csr_matrix]:
    keep_set = set(keep_ids)
    with gzip.open(umi_path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        col_idx = [i for i, c in enumerate(cell_ids) if c in keep_set]
        if not col_idx:
            raise SystemExit("no requested GSE131907 cell IDs in UMI header")
        ordered_cells = [cell_ids[i] for i in col_idx]
        print(
            f"GSE131907 UMI header cells={len(cell_ids)} keep={len(col_idx)}",
            flush=True,
        )
        genes: list[str] = []
        data: list[np.ndarray] = []
        indices: list[np.ndarray] = []
        indptr = [0]
        nnz_total = 0
        for gi, line in enumerate(f, start=1):
            raw = line.rstrip("\n")
            if not raw:
                continue
            tab0 = raw.find("\t")
            gene = raw[:tab0]
            rest = raw[tab0 + 1 :]
            vals = np.fromstring(rest, sep="\t", dtype=np.float32)
            if vals.size != len(cell_ids):
                raise SystemExit(
                    f"row {gi} {gene}: expected {len(cell_ids)} values, got {vals.size}"
                )
            sub = vals[col_idx]
            nz = np.flatnonzero(sub)
            if nz.size:
                data.append(sub[nz].astype(np.float32, copy=False))
                indices.append(nz.astype(np.int32, copy=False))
                nnz_total += int(nz.size)
            indptr.append(nnz_total)
            genes.append(gene)
            if gi % 2000 == 0:
                print(f"  streamed {gi} genes, nnz={nnz_total}", flush=True)
            if gene_cap is not None and gi >= gene_cap:
                break
    if data:
        data_a = np.concatenate(data)
        indices_a = np.concatenate(indices)
    else:
        data_a = np.array([], dtype=np.float32)
        indices_a = np.array([], dtype=np.int32)
    mat = sparse.csr_matrix(
        (data_a, indices_a, np.asarray(indptr, dtype=np.int64)),
        shape=(len(genes), len(ordered_cells)),
        dtype=np.float32,
    )
    return ordered_cells, genes, mat


def cap_barcodes(
    frame: pd.DataFrame,
    unit_col: str,
    cap: int,
    seed: int,
    protect: pd.Series | None = None,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    keep_idx: list[int] = []
    for _, sub in frame.groupby(unit_col, observed=True):
        if protect is not None:
            prot = sub.index[protect.loc[sub.index].to_numpy()]
        else:
            prot = sub.index[:0]
        prot = pd.Index(prot)
        rest = sub.index.difference(prot)
        n_rest = max(0, cap - len(prot))
        if len(rest) > n_rest:
            chosen = rng.choice(rest.to_numpy(), size=n_rest, replace=False)
            picked = prot.append(pd.Index(chosen))
        else:
            picked = prot.append(rest)
        keep_idx.extend(picked.tolist())
    return frame.loc[keep_idx].copy()


def _sanitize_obs(obs: pd.DataFrame) -> pd.DataFrame:
    out = obs.copy()
    for col in out.columns:
        if pd.api.types.is_bool_dtype(out[col]) or str(out[col].dtype) == "boolean":
            out[col] = out[col].map({True: "True", False: "False"}).astype(str)
        elif isinstance(out[col].dtype, pd.CategoricalDtype) or out[col].dtype == object:
            out[col] = out[col].astype(str).replace({"nan": "NA", "None": "NA", "<NA>": "NA"})
    return out


def _upper_var(adata):
    names = pd.Index([str(x).strip().upper() for x in adata.var_names])
    # keep first on collision
    seen = {}
    keep = []
    for i, n in enumerate(names):
        if n and n not in seen:
            seen[n] = i
            keep.append(i)
    adata = adata[:, keep].copy()
    adata.var_names = names[keep]
    adata.var_names_make_unique()
    return adata


def extract_gse131907(data: Path, gene_cap: int | None, cap: int) -> "ad.AnnData":
    import anndata as ad

    ann_path = data / "GSE131907" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    umi_path = data / "GSE131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    # also accept winpair-style nested folder
    if not ann_path.is_file():
        ann_path = data / "gse131907" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
        umi_path = data / "gse131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    ann = pd.read_csv(ann_path, sep="\t", compression="gzip")
    epi = ann.loc[
        (ann["Cell_type"] == EPI_TYPE_131907)
        & (ann["Sample_Origin"].isin(KEEP_ORIGINS_131907))
    ].copy()
    epi["Index"] = epi["Index"].astype(str)
    epi = epi.set_index("Index", drop=False)
    protect = epi["Cell_subtype"].astype(str).eq("AT2") & epi["Sample_Origin"].eq("nLung")
    catalog = {
        "n_ann_rows": int(len(ann)),
        "n_epithelial_kept_origins": int(len(epi)),
        "by_origin": epi["Sample_Origin"].value_counts().to_dict(),
        "by_subtype": epi["Cell_subtype"].fillna("NA").value_counts().to_dict(),
        "n_samples": int(epi["Sample"].nunique()),
        "n_nLung_AT2": int(protect.sum()),
    }
    print(json.dumps({"GSE131907_catalog": catalog}, indent=2), flush=True)
    epi = cap_barcodes(epi, "Sample", cap, SEED, protect)
    keep_ids = epi["Index"].astype(str).tolist()
    cells, genes, mat = stream_umi_subset(umi_path, keep_ids, gene_cap=gene_cap)
    X = mat.T.tocsr()
    obs = epi.set_index("Index").reindex(cells)
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(genes, name="gene")))
    adata = _upper_var(adata)
    adata.obs["dataset"] = "GSE131907"
    adata.obs["histology"] = "LUAD"
    adata.obs["unit_id"] = "GSE131907:" + adata.obs["Sample"].astype(str)
    adata.obs["patient_id"] = (
        adata.obs["Sample"]
        .astype(str)
        .str.extract(r"(?:LUNG_[NT]|EBUS_|BRONCHO_|EFFUSION_)?(\d+)", expand=False)
    )
    adata.obs["author_lineage"] = "Epithelial cells"
    adata.obs["author_subtype"] = adata.obs["Cell_subtype"].astype(str)
    adata.obs["tissue"] = adata.obs["Sample_Origin"].astype(str)
    adata.obs["is_normal_tissue"] = adata.obs["Sample_Origin"].eq("nLung").map(
        {True: "True", False: "False"}
    )
    adata.layers["counts"] = adata.X.copy()
    print(f"GSE131907 written cells={adata.n_obs} genes={adata.n_vars}", flush=True)
    return adata


def load_rds_matrix(path: Path):
    import rdata

    with tempfile.TemporaryDirectory(prefix="gse205335-") as tmp:
        matrix_path = path
        if str(path).endswith(".gz"):
            matrix_path = Path(tmp) / "matrix.rds"
            print(f"decompress {path.name}", flush=True)
            with gzip.open(path, "rb") as source, matrix_path.open("wb") as dest:
                shutil.copyfileobj(source, dest, 16 * 1024 * 1024)
        print("read GSE205335 RDS", flush=True)
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", message='Missing constructor for R class "dgCMatrix"'
            )
            obj = rdata.read_rds(matrix_path)
    expected = {"i", "p", "Dim", "Dimnames", "x"}
    if not expected.issubset(vars(obj)):
        raise TypeError("RDS is not the expected Matrix::dgCMatrix object")
    genes = np.asarray(obj.Dimnames[0], dtype=str)
    barcodes = np.asarray(obj.Dimnames[1], dtype=str)
    print(f"build CSC {tuple(obj.Dim)}", flush=True)
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False)
    return matrix, genes, barcodes


def extract_gse205335(data: Path, cap: int) -> "ad.AnnData":
    import anndata as ad

    ident_path = data / "GSE205335" / "GSE205335_Lung_IO_CellIdentity.txt.gz"
    soft_path = data / "GSE205335" / "GSE205335_family.soft.gz"
    rds_path = data / "GSE205335" / "GSE205335_Lung_IO_UMI_matrix.rds.gz"
    if not ident_path.is_file():
        ident_path = data / "gse205335" / "GSE205335_Lung_IO_CellIdentity.txt.gz"
        soft_path = data / "gse205335" / "GSE205335_family.soft.gz"
        rds_path = data / "gse205335" / "GSE205335_Lung_IO_UMI_matrix.rds.gz"
    ident = pd.read_csv(ident_path, sep="\t")
    meta = parse_geo_soft(soft_path)
    cells = ident.merge(
        meta[
            [
                "orig.ident",
                "gsm",
                "patient",
                "tissue",
                "recist",
                "cancer_subtype",
                "tumor_stage",
            ]
        ],
        on="orig.ident",
        how="left",
        validate="many_to_one",
    )
    if cells["patient"].isna().any():
        missing = cells.loc[cells["patient"].isna(), "orig.ident"].unique().tolist()
        raise SystemExit(f"GSE205335 identity samples missing GEO metadata: {missing}")
    cells["is_normal_tissue"] = cells["tissue"].astype(str).str.startswith(NORMAL_TISSUE_PREFIX)
    cells["is_epi"] = cells["lineage.total"].eq("Epithelial cells")
    epi = cells.loc[cells["is_epi"] & ~cells["is_normal_tissue"]].copy()
    catalog = {
        "n_ident": int(len(ident)),
        "n_epithelial_non_normal": int(len(epi)),
        "by_lineage_sub": epi["lineage.sub"].value_counts().to_dict(),
        "by_tissue": epi["tissue"].value_counts().to_dict(),
        "n_patients": int(epi["patient"].nunique()),
        "n_normal_epi_dropped": int((cells["is_epi"] & cells["is_normal_tissue"]).sum()),
    }
    print(json.dumps({"GSE205335_catalog": catalog}, indent=2), flush=True)
    epi = cap_barcodes(epi.set_index("barcode", drop=False), "patient", cap, SEED + 1)
    keep = set(epi["barcode"].astype(str))
    matrix, genes, barcodes = load_rds_matrix(rds_path)
    col_idx = np.array([i for i, b in enumerate(barcodes) if b in keep], dtype=np.int64)
    print(f"GSE205335 subset columns {len(col_idx)} / {len(barcodes)}", flush=True)
    sub = matrix[:, col_idx].tocsr()
    del matrix
    ordered = [str(barcodes[i]) for i in col_idx]
    obs = epi.set_index("barcode").reindex(ordered)
    adata = ad.AnnData(
        X=sub.T.tocsr(),
        obs=obs,
        var=pd.DataFrame(index=pd.Index(genes, name="gene")),
    )
    adata = _upper_var(adata)
    adata.obs["dataset"] = "GSE205335"
    adata.obs["histology"] = adata.obs["cancer_subtype"].astype(str)
    adata.obs["Sample"] = adata.obs["patient"].astype(str)
    adata.obs["Sample_Origin"] = adata.obs["tissue"].astype(str)
    adata.obs["Cell_type"] = "Epithelial cells"
    adata.obs["Cell_subtype"] = adata.obs["lineage.sub"].astype(str)
    adata.obs["unit_id"] = "GSE205335:" + adata.obs["patient"].astype(str)
    adata.obs["patient_id"] = adata.obs["patient"].astype(str)
    adata.obs["author_lineage"] = adata.obs["lineage.total"].astype(str)
    adata.obs["author_subtype"] = adata.obs["lineage.sub"].astype(str)
    adata.obs["is_normal_tissue"] = "False"
    adata.layers["counts"] = adata.X.copy()
    print(f"GSE205335 written cells={adata.n_obs} genes={adata.n_vars}", flush=True)
    return adata


def parse_123902_name(name: str) -> tuple[str, str, str]:
    stem = Path(name).name.replace(".csv.gz", "")
    parts = stem.split("_")
    gsm = parts[0] if parts else stem
    donor = parts[2] if len(parts) > 2 else stem
    if "NORMAL" in name:
        tissue = "NORMAL"
    elif "METASTASIS" in name:
        tissue = "METASTASIS"
    elif "PRIMARY" in name:
        tissue = "PRIMARY"
    else:
        tissue = "OTHER"
    return gsm, donor, tissue


def _read_123902_csv(handle) -> pd.DataFrame:
    df = pd.read_csv(handle, index_col=0)
    df.columns = [str(c).strip().upper() for c in df.columns]
    df.index = df.index.astype(str)
    return df


def extract_gse123902(data: Path, cap: int) -> "ad.AnnData":
    import anndata as ad

    tar_path = data / "GSE123902" / "GSE123902_RAW.tar"
    if not tar_path.is_file():
        tar_path = data / "gse123902" / "GSE123902_RAW.tar"
    if not tar_path.is_file():
        raise SystemExit(f"missing {tar_path}")

    blocks: list[tuple[pd.DataFrame, pd.DataFrame]] = []
    catalog_rows = []
    with tarfile.open(tar_path) as tf:
        members = [m for m in tf.getmembers() if m.name.endswith(".csv.gz")]
        members.sort(key=lambda m: m.name)
        print(f"GSE123902 tar members csv.gz={len(members)}", flush=True)
        for m in members:
            name = Path(m.name).name
            gsm, donor, tissue = parse_123902_name(name)
            if tissue == "OTHER":
                print(f"  skip {name} (unclassified)", flush=True)
                continue
            raw = gzip.GzipFile(fileobj=tf.extractfile(m))
            df = _read_123902_csv(io.TextIOWrapper(raw, encoding="utf-8", errors="replace"))
            n_cells, n_genes = df.shape
            epi_cols = [g for g in EPI_MARKERS if g in df.columns]
            if not epi_cols:
                raise SystemExit(f"{name}: none of {EPI_MARKERS} in columns")
            epi_any = (df[epi_cols] > 0).any(axis=1)
            epi = df.loc[epi_any].copy()
            sftpc = epi["SFTPC"] if "SFTPC" in epi.columns else pd.Series(0, index=epi.index)
            protect = (tissue == "NORMAL") & (sftpc.to_numpy() > 0)
            obs = pd.DataFrame(
                {
                    "barcode": epi.index.astype(str),
                    "gsm": gsm,
                    "patient": donor,
                    "tissue": tissue,
                    "file": name,
                    "Sample": f"{donor}_{tissue}",
                    "Sample_Origin": tissue,
                    "protect_at2like": protect,
                },
                index=epi.index.astype(str),
            )
            catalog_rows.append(
                {
                    "file": name,
                    "patient": donor,
                    "tissue": tissue,
                    "n_cells": int(n_cells),
                    "n_marker_epi": int(epi.shape[0]),
                    "n_genes": int(n_genes),
                    "n_sftpc_pos": int((sftpc > 0).sum()),
                }
            )
            print(
                f"  {name}: cells={n_cells} marker_epi={epi.shape[0]} "
                f"sftpc+={int((sftpc > 0).sum())} {donor} {tissue}",
                flush=True,
            )
            if epi.empty:
                continue
            blocks.append((epi, obs))

    if not blocks:
        raise SystemExit("GSE123902: no marker-epithelial cells")

    gene_sets = [set(e.columns) for e, _ in blocks]
    shared = set.intersection(*gene_sets)
    if len(shared) < 2000:
        raise SystemExit(f"GSE123902 too few shared genes across samples: {len(shared)}")
    genes = sorted(shared)
    mats = []
    obs_parts = []
    for epi, obs in blocks:
        sub = epi.loc[:, genes].to_numpy(dtype=np.float32)
        mats.append(sparse.csr_matrix(sub))
        obs_parts.append(obs)
    X = sparse.vstack(mats).tocsr()
    obs = pd.concat(obs_parts, axis=0)
    obs["unit_id"] = "GSE123902:" + obs["Sample"].astype(str)
    obs = cap_barcodes(obs, "unit_id", cap, SEED + 2, protect=obs["protect_at2like"])
    # recap X to remaining barcodes
    keep = set(obs.index.astype(str))
    # rebuild from blocks after cap — simpler: slice by reindex
    all_barcodes = []
    for _, o in obs_parts:
        all_barcodes.extend(o.index.astype(str).tolist())
    barcode_to_row = {b: i for i, b in enumerate(all_barcodes)}
    row_idx = [barcode_to_row[b] for b in obs.index.astype(str)]
    X = X[row_idx, :]
    adata = ad.AnnData(
        X=X,
        obs=obs,
        var=pd.DataFrame(index=pd.Index(genes, name="gene")),
    )
    adata.obs["dataset"] = "GSE123902"
    adata.obs["histology"] = "LUAD"
    adata.obs["patient_id"] = adata.obs["patient"].astype(str)
    adata.obs["author_lineage"] = "marker_epithelial"
    adata.obs["author_subtype"] = np.where(
        adata.obs["Sample_Origin"].eq("NORMAL") & adata.obs["protect_at2like"].astype(bool),
        "AT2_like",
        "marker_epithelial",
    )
    adata.obs["Cell_type"] = "marker_epithelial"
    adata.obs["Cell_subtype"] = adata.obs["author_subtype"].astype(str)
    adata.obs["is_normal_tissue"] = adata.obs["Sample_Origin"].eq("NORMAL").map(
        {True: "True", False: "False"}
    )
    adata.layers["counts"] = adata.X.copy()
    print(
        json.dumps(
            {
                "GSE123902_catalog": catalog_rows,
                "n_cells_capped": int(adata.n_obs),
                "n_units": int(adata.obs["unit_id"].nunique()),
                "n_shared_genes": int(len(genes)),
            },
            indent=2,
        ),
        flush=True,
    )
    print(f"GSE123902 written cells={adata.n_obs} genes={adata.n_vars}", flush=True)
    return adata


def concat_shared(objects: list) -> tuple:
    import anndata as ad

    shared = objects[0].var_names
    for a in objects[1:]:
        shared = shared.intersection(a.var_names)
    if len(shared) < 5000:
        raise SystemExit(f"too few shared genes: {len(shared)}")
    copies = []
    for a in objects:
        b = a[:, shared].copy()
        ds = str(b.obs["dataset"].iloc[0])
        b.obs_names = ds + ":" + b.obs_names.astype(str)
        copies.append(b)
    cols = sorted(set().union(*[set(c.obs.columns) for c in copies]))
    for frame in copies:
        for c in cols:
            if c not in frame.obs.columns:
                frame.obs[c] = "NA"
        frame.obs = _sanitize_obs(frame.obs[cols])
    out = ad.concat(copies, axis=0, join="inner", merge="same")
    out.obs = _sanitize_obs(out.obs)
    out.layers["counts"] = out.X.copy()
    return out, int(len(shared))


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=Path, default=Path("/tmp/triple_differ_raw"))
    p.add_argument(
        "--out",
        type=Path,
        default=Path("/tmp/triple_differ_raw/epithelium.h5ad"),
    )
    p.add_argument("--cap-per-unit", type=int, default=CAP_PER_UNIT)
    p.add_argument("--gene-cap", type=int, default=None, help="debug only")
    args = p.parse_args()

    import anndata as ad

    caches = {
        "GSE123902": args.data / "GSE123902_epithelium_cap.h5ad",
        "GSE131907": args.data / "GSE131907_epithelium_cap.h5ad",
        "GSE205335": args.data / "GSE205335_epithelium_cap.h5ad",
    }
    loaded = {}
    if caches["GSE123902"].is_file() and args.gene_cap is None:
        print(f"reuse {caches['GSE123902']}", flush=True)
        loaded["GSE123902"] = ad.read_h5ad(caches["GSE123902"])
    else:
        loaded["GSE123902"] = extract_gse123902(args.data, args.cap_per_unit)
        loaded["GSE123902"].obs = _sanitize_obs(loaded["GSE123902"].obs)
        loaded["GSE123902"].write_h5ad(caches["GSE123902"])
        print(f"cached {caches['GSE123902']}", flush=True)

    if caches["GSE131907"].is_file() and args.gene_cap is None:
        print(f"reuse {caches['GSE131907']}", flush=True)
        loaded["GSE131907"] = ad.read_h5ad(caches["GSE131907"])
    else:
        loaded["GSE131907"] = extract_gse131907(args.data, args.gene_cap, args.cap_per_unit)
        loaded["GSE131907"].obs = _sanitize_obs(loaded["GSE131907"].obs)
        loaded["GSE131907"].write_h5ad(caches["GSE131907"])
        print(f"cached {caches['GSE131907']}", flush=True)

    if caches["GSE205335"].is_file() and args.gene_cap is None:
        print(f"reuse {caches['GSE205335']}", flush=True)
        loaded["GSE205335"] = ad.read_h5ad(caches["GSE205335"])
    else:
        loaded["GSE205335"] = extract_gse205335(args.data, args.cap_per_unit)
        loaded["GSE205335"].obs = _sanitize_obs(loaded["GSE205335"].obs)
        loaded["GSE205335"].write_h5ad(caches["GSE205335"])
        print(f"cached {caches['GSE205335']}", flush=True)

    joint, n_shared = concat_shared([loaded[c] for c in COHORTS])
    del loaded
    args.out.parent.mkdir(parents=True, exist_ok=True)
    joint.write_h5ad(args.out)
    inv = {
        "n_cells": int(joint.n_obs),
        "n_genes_shared": n_shared,
        "by_dataset": joint.obs["dataset"].value_counts().to_dict(),
        "by_origin": joint.obs["Sample_Origin"].astype(str).value_counts().to_dict(),
        "by_subtype": joint.obs["author_subtype"].astype(str).value_counts().to_dict(),
        "n_units": int(joint.obs["unit_id"].nunique()),
        "cap_per_unit": args.cap_per_unit,
        "gse148071_added": False,
        "seven_pool": False,
        "out": str(args.out),
    }
    args.out.with_suffix(".inventory.json").write_text(json.dumps(inv, indent=2))
    print(json.dumps(inv, indent=2), flush=True)


if __name__ == "__main__":
    main()
