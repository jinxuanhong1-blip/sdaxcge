#!/usr/bin/env python3
"""Extract subsampled malignant + T/NK AnnData for the winning pair.

A previous 53k-cell Harmony/UMAP run OOM'd on ~15 GB. This path caps
cells per unit *before* matrix subset:

- malignant: 150 / unit
- T/NK: 80 / unit

GSE131907 unit = sample with n_malignant ≥ 20 on tumor origins.
GSE205335 unit = patient (drop Normal* libraries).
GSE148071 and GSE207422 are not added.
"""
from __future__ import annotations

import argparse
import gzip
import json
import shutil
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
TUMOR_ORIGINS_131907 = ("tLung", "tL/B", "mLN", "PE", "mBrain")
MIN_MALIGNANT = 20
MAL_CAP = 150
TNK_CAP = 80
SEED = 1


def winning_gse131907_samples(given: Path) -> list[str]:
    df = pd.read_csv(given / "GSE131907_samples.tsv", sep="\t")
    keep = df.loc[
        df["origin"].isin(TUMOR_ORIGINS_131907)
        & (pd.to_numeric(df["n_malignant"], errors="coerce") >= MIN_MALIGNANT),
        "sample",
    ]
    return [str(x) for x in keep]


def winning_gse205335_patients(given: Path) -> list[str]:
    df = pd.read_csv(given / "GSE205335_patients.tsv", sep="\t")
    return [str(x) for x in df["patient"]]


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


def stream_umi_subset(umi_path: Path, keep_ids: list[str]) -> tuple[list[str], list[str], sparse.csr_matrix]:
    keep_set = set(keep_ids)
    with gzip.open(umi_path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        col_idx = [i for i, c in enumerate(cell_ids) if c in keep_set]
        if not col_idx:
            raise SystemExit("no requested GSE131907 cell IDs in UMI header")
        ordered_cells = [cell_ids[i] for i in col_idx]
        print(f"GSE131907 UMI header cells={len(cell_ids)} keep={len(col_idx)}", flush=True)
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
                raise SystemExit(f"row {gi} {gene}: expected {len(cell_ids)} values, got {vals.size}")
            sub = vals[col_idx]
            nz = np.flatnonzero(sub)
            if nz.size:
                data.append(sub[nz].astype(np.float32, copy=False))
                indices.append(nz.astype(np.int32, copy=False))
                nnz_total += int(nz.size)
            indptr.append(nnz_total)
            genes.append(gene)
            if gi % 4000 == 0:
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
    return metadata.rename(columns={"tumor stage": "tumor_stage", "cancer subtype": "cancer_subtype"})


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
            warnings.filterwarnings("ignore", message='Missing constructor for R class "dgCMatrix"')
            obj = rdata.read_rds(matrix_path)
    expected = {"i", "p", "Dim", "Dimnames", "x"}
    if not expected.issubset(vars(obj)):
        raise TypeError("RDS is not the expected Matrix::dgCMatrix object")
    genes = np.asarray(obj.Dimnames[0], dtype=str)
    barcodes = np.asarray(obj.Dimnames[1], dtype=str)
    print(f"build CSC {tuple(obj.Dim)}", flush=True)
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False)
    return matrix, genes, barcodes


def extract_gse131907(data: Path, given: Path, mal_cap: int, tnk_cap: int):
    import anndata as ad

    samples = set(winning_gse131907_samples(given))
    ann_path = data / "gse131907" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    umi_path = data / "gse131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    ann = pd.read_csv(ann_path, sep="\t", compression="gzip")
    ann["Index"] = ann["Index"].astype(str)
    tumor = ann["Sample"].astype(str).isin(samples)
    mal = ann.loc[tumor & ann["Cell_subtype"].eq("Malignant cells")].copy()
    tnk = ann.loc[tumor & ann["Cell_type"].isin(["T lymphocytes", "NK cells"])].copy()
    catalog = {
        "n_ann": int(len(ann)),
        "n_winning_samples": int(len(samples)),
        "n_malignant_catalog": int(len(mal)),
        "n_tnk_catalog": int(len(tnk)),
        "mal_by_sample": mal["Sample"].value_counts().to_dict(),
        "tnk_by_sample": tnk["Sample"].value_counts().to_dict(),
    }
    print(json.dumps({"GSE131907_catalog": catalog}, indent=2), flush=True)
    mal = cap_barcodes(mal.set_index("Index", drop=False), "Sample", mal_cap, SEED)
    tnk = cap_barcodes(tnk.set_index("Index", drop=False), "Sample", tnk_cap, SEED + 7)
    mal["compartment"] = "malignant"
    tnk["compartment"] = "TNK"
    keep = pd.concat([mal, tnk], axis=0)
    cells, genes, mat = stream_umi_subset(umi_path, keep["Index"].astype(str).tolist())
    X = mat.T.tocsr()
    obs = keep.set_index("Index").reindex(cells)
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(genes, name="gene")))
    adata.obs["dataset"] = "GSE131907"
    adata.obs["unit_id"] = "GSE131907:" + adata.obs["Sample"].astype(str)
    adata.obs["unit_type"] = "sample"
    adata.obs["author_lineage"] = adata.obs["Cell_type"].astype(str)
    adata.obs["author_subtype"] = adata.obs["Cell_subtype"].astype(str)
    adata.obs["tissue"] = adata.obs["Sample_Origin"].astype(str)
    adata.layers["counts"] = adata.X.copy()
    print(f"GSE131907 written cells={adata.n_obs} genes={adata.n_vars}", flush=True)
    return adata, catalog


def extract_gse205335(data: Path, given: Path, mal_cap: int, tnk_cap: int):
    import anndata as ad

    patients = set(winning_gse205335_patients(given))
    ident_path = data / "gse205335" / "GSE205335_Lung_IO_CellIdentity.txt.gz"
    soft_path = data / "gse205335" / "GSE205335_family.soft.gz"
    rds_path = data / "gse205335" / "GSE205335_Lung_IO_UMI_matrix.rds.gz"
    ident = pd.read_csv(ident_path, sep="\t")
    meta = parse_geo_soft(soft_path)
    cells = ident.merge(
        meta[["orig.ident", "gsm", "patient", "tissue", "recist", "cancer_subtype"]],
        on="orig.ident",
        how="left",
        validate="many_to_one",
    )
    if cells["patient"].isna().any():
        missing = cells.loc[cells["patient"].isna(), "orig.ident"].unique().tolist()
        raise SystemExit(f"GSE205335 identity samples missing GEO metadata: {missing}")
    cells["is_normal"] = cells["tissue"].astype(str).str.startswith("Normal")
    keep_pat = cells["patient"].astype(str).isin(patients) & ~cells["is_normal"]
    mal = cells.loc[keep_pat & cells["lineage.sub"].eq("Malignant cells")].copy()
    tnk = cells.loc[keep_pat & cells["lineage.total"].eq("T/NK cells")].copy()
    catalog = {
        "n_ident": int(len(ident)),
        "n_patients": int(len(patients)),
        "n_malignant_catalog": int(len(mal)),
        "n_tnk_catalog": int(len(tnk)),
        "n_normal_dropped": int(cells["is_normal"].sum()),
    }
    print(json.dumps({"GSE205335_catalog": catalog}, indent=2), flush=True)
    mal = cap_barcodes(mal.set_index("barcode", drop=False), "patient", mal_cap, SEED + 1)
    tnk = cap_barcodes(tnk.set_index("barcode", drop=False), "patient", tnk_cap, SEED + 8)
    mal["compartment"] = "malignant"
    tnk["compartment"] = "TNK"
    keep = pd.concat([mal, tnk], axis=0)
    keep_ids = set(keep["barcode"].astype(str))
    matrix, genes, barcodes = load_rds_matrix(rds_path)
    col_idx = np.array([i for i, b in enumerate(barcodes) if b in keep_ids], dtype=np.int64)
    print(f"GSE205335 subset columns {len(col_idx)} / {len(barcodes)}", flush=True)
    sub = matrix[:, col_idx].tocsr()
    del matrix
    ordered = [str(barcodes[i]) for i in col_idx]
    obs = keep.set_index("barcode").reindex(ordered)
    adata = ad.AnnData(
        X=sub.T.tocsr(),
        obs=obs,
        var=pd.DataFrame(index=pd.Index(genes, name="gene")),
    )
    adata.obs["dataset"] = "GSE205335"
    adata.obs["Sample"] = adata.obs["patient"].astype(str)
    adata.obs["Sample_Origin"] = adata.obs["tissue"].astype(str)
    adata.obs["Cell_type"] = adata.obs["lineage.total"].astype(str)
    adata.obs["Cell_subtype"] = adata.obs["lineage.sub"].astype(str)
    adata.obs["unit_id"] = "GSE205335:" + adata.obs["patient"].astype(str)
    adata.obs["unit_type"] = "patient"
    adata.obs["author_lineage"] = adata.obs["lineage.total"].astype(str)
    adata.obs["author_subtype"] = adata.obs["lineage.sub"].astype(str)
    adata.layers["counts"] = adata.X.copy()
    print(f"GSE205335 written cells={adata.n_obs} genes={adata.n_vars}", flush=True)
    return adata, catalog


def _sanitize_obs(obs: pd.DataFrame) -> pd.DataFrame:
    out = obs.copy()
    for col in out.columns:
        if pd.api.types.is_bool_dtype(out[col]) or str(out[col].dtype) == "boolean":
            out[col] = out[col].map({True: "True", False: "False"}).astype(str)
        elif out[col].dtype == object or str(out[col].dtype) == "category":
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


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=Path, default=Path("/tmp/pair_131907_205335"))
    p.add_argument("--given", type=Path, default=ROOT / "data" / "given")
    p.add_argument("--out", type=Path, default=Path("/tmp/pair_131907_205335/pair_subsample.h5ad"))
    p.add_argument("--mal-cap", type=int, default=MAL_CAP)
    p.add_argument("--tnk-cap", type=int, default=TNK_CAP)
    args = p.parse_args()

    cache_a = args.data / "GSE131907_cap.h5ad"
    cache_b = args.data / "GSE205335_cap.h5ad"
    import anndata as ad

    if cache_a.is_file():
        print(f"reuse {cache_a}", flush=True)
        a = ad.read_h5ad(cache_a)
        cat_a = {"reused": True}
    else:
        a, cat_a = extract_gse131907(args.data, args.given, args.mal_cap, args.tnk_cap)
        a.obs = _sanitize_obs(a.obs)
        a.write_h5ad(cache_a)
    if cache_b.is_file():
        print(f"reuse {cache_b}", flush=True)
        b = ad.read_h5ad(cache_b)
        cat_b = {"reused": True}
    else:
        b, cat_b = extract_gse205335(args.data, args.given, args.mal_cap, args.tnk_cap)
        b.obs = _sanitize_obs(b.obs)
        b.write_h5ad(cache_b)

    joint, n_shared = concat_shared(a, b)
    del a, b
    args.out.parent.mkdir(parents=True, exist_ok=True)
    joint.write_h5ad(args.out)
    inv = {
        "n_cells": int(joint.n_obs),
        "n_genes_shared": n_shared,
        "by_dataset": joint.obs["dataset"].value_counts().to_dict(),
        "by_compartment": joint.obs["compartment"].value_counts().to_dict(),
        "n_units": int(joint.obs["unit_id"].nunique()),
        "mal_cap_per_unit": args.mal_cap,
        "tnk_cap_per_unit": args.tnk_cap,
        "subsampled": True,
        "previous_oom_cells": 53296,
        "gse148071_added": False,
        "gse207422_added": False,
        "GSE131907": cat_a,
        "GSE205335": cat_b,
        "out": str(args.out),
    }
    args.out.with_suffix(".inventory.json").write_text(json.dumps(inv, indent=2))
    print(json.dumps(inv, indent=2), flush=True)


if __name__ == "__main__":
    main()
