#!/usr/bin/env python3
"""Build a joint author-malignant AnnData for GSE131907 + GSE205335.

ADDITIVE. Author malignant cells only. No GSE148071. No GSE207422.
GSE131907: Cell_subtype == 'Malignant cells' (tLung has 0; those cells are tS*).
GSE205335: lineage.sub == 'Malignant cells' on non-normal tissues.
"""
from __future__ import annotations

import argparse
import gzip
import json
import shutil
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

NORMAL_TISSUE_PREFIX = ("Normal ",)
CAP_PER_UNIT = 200
SEED = 1
MIN_UMI = 200
MIN_GENES = 200


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


def cap_barcodes(frame: pd.DataFrame, unit_col: str, cap: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    keep_idx: list[int] = []
    for _, sub in frame.groupby(unit_col, observed=True):
        if len(sub) > cap:
            chosen = rng.choice(sub.index.to_numpy(), size=cap, replace=False)
            keep_idx.extend(chosen.tolist())
        else:
            keep_idx.extend(sub.index.tolist())
    return frame.loc[keep_idx].copy()


def extract_gse131907(data: Path, cap: int) -> "ad.AnnData":
    import anndata as ad

    ann_path = data / "gse131907" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    umi_path = data / "gse131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    ann = pd.read_csv(ann_path, sep="\t", compression="gzip")
    subtype_counts = ann["Cell_subtype"].fillna("NA").value_counts().to_dict()
    tlung = ann[ann["Sample_Origin"] == "tLung"]
    mal = ann.loc[ann["Cell_subtype"].astype(str).eq("Malignant cells")].copy()
    catalog = {
        "n_ann_rows": int(len(ann)),
        "n_author_malignant": int(len(mal)),
        "by_origin": mal["Sample_Origin"].value_counts().to_dict(),
        "n_samples": int(mal["Sample"].nunique()),
        "tLung_epithelial": int((tlung["Cell_type"] == "Epithelial cells").sum()),
        "tLung_author_malignant": int((tlung["Cell_subtype"] == "Malignant cells").sum()),
        "tLung_subtype_counts": tlung["Cell_subtype"].fillna("NA").value_counts().to_dict(),
        "author_subtype_counts_all": subtype_counts,
        "malignant_rule": "author_Cell_subtype==Malignant cells; tLung is 0 (tS1/tS2/tS3)",
    }
    print(json.dumps({"GSE131907_catalog": catalog}, indent=2), flush=True)
    mal["Index"] = mal["Index"].astype(str)
    mal = mal.set_index("Index", drop=False)
    mal = cap_barcodes(mal, "Sample", cap, SEED)
    keep_ids = mal["Index"].astype(str).tolist()
    cells, genes, mat = stream_umi_subset(umi_path, keep_ids)
    X = mat.T.tocsr()
    obs = mal.set_index("Index").reindex(cells)
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(genes, name="gene")))
    adata.obs["dataset"] = "GSE131907"
    adata.obs["histology"] = "LUAD"
    adata.obs["unit_id"] = "GSE131907:" + adata.obs["Sample"].astype(str)
    adata.obs["patient_id"] = adata.obs["Sample"].astype(str)
    adata.obs["author_lineage"] = adata.obs["Cell_type"].astype(str)
    adata.obs["author_subtype"] = adata.obs["Cell_subtype"].astype(str)
    adata.obs["tissue"] = adata.obs["Sample_Origin"].astype(str)
    adata.obs["is_normal_tissue"] = False
    adata.obs["malignant_rule"] = "author_Cell_subtype==Malignant cells"
    adata.layers["counts"] = adata.X.copy()
    print(f"GSE131907 written cells={adata.n_obs} genes={adata.n_vars}", flush=True)
    adata.uns["catalog"] = catalog
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
    cells["is_malignant"] = cells["lineage.sub"].astype(str).eq("Malignant cells")
    mal = cells.loc[cells["is_malignant"] & ~cells["is_normal_tissue"]].copy()
    catalog = {
        "n_ident": int(len(ident)),
        "n_author_malignant_non_normal": int(len(mal)),
        "by_lineage_sub": mal["lineage.sub"].value_counts().to_dict(),
        "by_tissue": mal["tissue"].value_counts().to_dict(),
        "by_subtype": mal["cancer_subtype"].value_counts().to_dict(),
        "n_patients": int(mal["patient"].nunique()),
        "n_normal_malignant_dropped": int((cells["is_malignant"] & cells["is_normal_tissue"]).sum()),
        "malignant_rule": "author_lineage.sub==Malignant cells; drop normal tissues",
    }
    print(json.dumps({"GSE205335_catalog": catalog}, indent=2), flush=True)
    mal = cap_barcodes(mal.set_index("barcode", drop=False), "patient", cap, SEED + 1)
    keep = set(mal["barcode"].astype(str))
    matrix, genes, barcodes = load_rds_matrix(rds_path)
    col_idx = np.array([i for i, b in enumerate(barcodes) if b in keep], dtype=np.int64)
    print(f"GSE205335 subset columns {len(col_idx)} / {len(barcodes)}", flush=True)
    sub = matrix[:, col_idx].tocsr()
    del matrix
    ordered = [str(barcodes[i]) for i in col_idx]
    obs = mal.set_index("barcode").reindex(ordered)
    adata = ad.AnnData(
        X=sub.T.tocsr(),
        obs=obs,
        var=pd.DataFrame(index=pd.Index(genes, name="gene")),
    )
    adata.obs["dataset"] = "GSE205335"
    adata.obs["histology"] = adata.obs["cancer_subtype"].astype(str)
    adata.obs["Sample"] = adata.obs["patient"].astype(str)
    adata.obs["Sample_Origin"] = adata.obs["tissue"].astype(str)
    adata.obs["Cell_type"] = adata.obs["lineage.total"].astype(str)
    adata.obs["Cell_subtype"] = adata.obs["lineage.sub"].astype(str)
    adata.obs["unit_id"] = "GSE205335:" + adata.obs["patient"].astype(str)
    adata.obs["patient_id"] = adata.obs["patient"].astype(str)
    adata.obs["author_lineage"] = adata.obs["lineage.total"].astype(str)
    adata.obs["author_subtype"] = adata.obs["lineage.sub"].astype(str)
    adata.obs["malignant_rule"] = "author_lineage.sub==Malignant cells"
    adata.layers["counts"] = adata.X.copy()
    print(f"GSE205335 written cells={adata.n_obs} genes={adata.n_vars}", flush=True)
    adata.uns["catalog"] = catalog
    return adata


def _sanitize_obs(obs: pd.DataFrame) -> pd.DataFrame:
    out = obs.copy()
    for col in out.columns:
        if pd.api.types.is_bool_dtype(out[col]) or str(out[col].dtype) == "boolean":
            out[col] = out[col].map({True: "True", False: "False"}).astype(str)
        elif pd.api.types.is_object_dtype(out[col]) or str(out[col].dtype) == "category":
            out[col] = out[col].astype(str).replace({"nan": "NA", "None": "NA", "<NA>": "NA"})
    return out


def concat_shared(a, b):
    import anndata as ad

    shared = a.var_names.intersection(b.var_names)
    if len(shared) < 5000:
        raise SystemExit(f"too few shared genes: {len(shared)}")
    a2 = a[:, shared].copy()
    b2 = b[:, shared].copy()
    a2.obs_names = "GSE131907:" + a2.obs_names.astype(str)
    b2.obs_names = "GSE205335:" + b2.obs_names.astype(str)
    cols = sorted(set(a2.obs.columns) | set(b2.obs.columns))
    for frame in (a2, b2):
        for c in cols:
            if c not in frame.obs.columns:
                frame.obs[c] = "NA"
        frame.obs = _sanitize_obs(frame.obs[cols])
    out = ad.concat([a2, b2], axis=0, join="inner", merge="same")
    out.obs = _sanitize_obs(out.obs)
    out.layers["counts"] = out.X.copy()
    return out, int(len(shared))


def qc_filter(adata, min_umi: int, min_genes: int):
    from scipy import sparse as sp

    X = adata.layers["counts"] if "counts" in adata.layers else adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    n_umi = np.asarray(X.sum(axis=1)).ravel()
    n_genes = np.asarray((X > 0).sum(axis=1)).ravel()
    adata.obs["n_umi"] = n_umi
    adata.obs["n_genes"] = n_genes
    keep = (n_umi >= min_umi) & (n_genes >= min_genes)
    print(
        f"QC keep {int(keep.sum())}/{adata.n_obs} (UMI>={min_umi}, genes>={min_genes})",
        flush=True,
    )
    return adata[keep].copy()


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=Path, default=Path("/tmp/winpair_cytotrace2_cldn4"))
    p.add_argument("--out", type=Path, default=Path("/tmp/winpair_cytotrace2_cldn4/malignant.h5ad"))
    p.add_argument("--cap-per-unit", type=int, default=CAP_PER_UNIT)
    args = p.parse_args()

    cache_a = args.data / "GSE131907_malignant_cap.h5ad"
    cache_b = args.data / "GSE205335_malignant_cap.h5ad"
    import anndata as ad

    if cache_a.is_file():
        print(f"reuse {cache_a}", flush=True)
        a = ad.read_h5ad(cache_a)
    else:
        a = extract_gse131907(args.data, args.cap_per_unit)
        a.obs = _sanitize_obs(a.obs)
        a.write_h5ad(cache_a)
        print(f"cached {cache_a}", flush=True)
    if cache_b.is_file():
        print(f"reuse {cache_b}", flush=True)
        b = ad.read_h5ad(cache_b)
    else:
        b = extract_gse205335(args.data, args.cap_per_unit)
        b.obs = _sanitize_obs(b.obs)
        b.write_h5ad(cache_b)
        print(f"cached {cache_b}", flush=True)
    joint, n_shared = concat_shared(a, b)
    del a, b
    joint = qc_filter(joint, MIN_UMI, MIN_GENES)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    joint.write_h5ad(args.out)
    inv = {
        "n_cells": int(joint.n_obs),
        "n_genes_shared": n_shared,
        "by_dataset": joint.obs["dataset"].value_counts().to_dict(),
        "by_origin": joint.obs["Sample_Origin"].astype(str).value_counts().to_dict(),
        "by_histology": joint.obs["histology"].astype(str).value_counts().to_dict(),
        "n_units": int(joint.obs["unit_id"].nunique()),
        "cap_per_unit": args.cap_per_unit,
        "gse148071_added": False,
        "gse207422_added": False,
        "dual_high": False,
        "malignant_only": True,
        "out": str(args.out),
    }
    args.out.with_suffix(".inventory.json").write_text(json.dumps(inv, indent=2))
    print(json.dumps(inv, indent=2), flush=True)


if __name__ == "__main__":
    main()
