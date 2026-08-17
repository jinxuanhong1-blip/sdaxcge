#!/usr/bin/env python3
"""Build a joint epithelial AnnData for GSE189357 + GSE205335.

ADDITIVE. Epithelium only. PR #459 is not re-scored.
GSE189357: marker epithelium (EPCAM|KRT8|KRT18|KRT19)>0 and PTPRC==0.
GSE205335: author lineage.total == Epithelial cells on non-normal tissues.
No dual-high gate. GSE131907 is not added.
"""
from __future__ import annotations

import argparse
import gzip
import json
import shutil
import tarfile
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

HERE = Path(__file__).resolve().parents[1]
DATA = HERE / "data"
EPI = ["EPCAM", "KRT8", "KRT18", "KRT19"]
NORMAL_TISSUE_PREFIX = ("Normal ",)
CAP_PER_UNIT = 350
SEED = 1


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


def _sanitize_obs(obs: pd.DataFrame) -> pd.DataFrame:
    out = obs.copy()
    for col in out.columns:
        if pd.api.types.is_bool_dtype(out[col]) or str(out[col].dtype) == "boolean":
            out[col] = out[col].map({True: "True", False: "False"}).astype(str)
        elif pd.api.types.is_categorical_dtype(out[col]) or out[col].dtype == object:
            out[col] = out[col].astype(str).replace({"nan": "NA", "None": "NA", "<NA>": "NA"})
    return out


def _unique_genes(names: list[str]) -> list[str]:
    seen: dict[str, int] = {}
    out = []
    for raw in names:
        gene = raw if raw else "NA"
        if gene in seen:
            seen[gene] += 1
            out.append(f"{gene}.{seen[gene]}")
        else:
            seen[gene] = 0
            out.append(gene)
    return out


def _read_10x_features(path: Path) -> list[str]:
    genes = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            parts = line.strip().split("\t")
            symbol = parts[1] if len(parts) > 1 else parts[0]
            genes.append(symbol.upper())
    return _unique_genes(genes)


def _read_10x_barcodes(path: Path) -> list[str]:
    with gzip.open(path, "rt") as handle:
        return [line.strip().split("\t")[0] for line in handle if line.strip()]


def _member(members: dict[str, tarfile.TarInfo], sample: str, kind: str) -> str:
    hits = [n for n in members if f"_{sample}_{kind}" in n]
    if not hits:
        raise SystemExit(f"GSE189357 {sample}: missing {kind}")
    return hits[0]


def extract_gse189357(tar_path: Path, scratch: Path, cap: int) -> "ad.AnnData":
    import anndata as ad

    units = pd.read_csv(DATA / "GSE189357_marker_units.tsv", sep="\t")
    keep = set(units.loc[units["eligible"].astype(str).str.lower() == "true", "patient"].astype(str))
    pieces = []
    catalog_rows = []
    with tarfile.open(tar_path) as tf:
        members = {m.name: m for m in tf.getmembers()}
        for i, sample in enumerate([f"TD{k}" for k in range(1, 10)]):
            if sample not in keep:
                print(f"  skip {sample} (not PR #459 eligible)", flush=True)
                continue
            print(f"GSE189357 {sample}", flush=True)
            dest = scratch / sample
            dest.mkdir(parents=True, exist_ok=True)
            local = {}
            for kind in ("features", "barcodes", "matrix"):
                name = _member(members, sample, kind)
                out = dest / Path(name).name
                if not (out.exists() and out.stat().st_size > 1000):
                    with tf.extractfile(members[name]) as src, out.open("wb") as fh:
                        shutil.copyfileobj(src, fh, 8 * 1024 * 1024)
                local[kind] = out
            genes = _read_10x_features(local["features"])
            barcodes = _read_10x_barcodes(local["barcodes"])
            mtx_plain = dest / f"{sample}_matrix.mtx"
            if not (mtx_plain.exists() and mtx_plain.stat().st_size > 1_000_000):
                with gzip.open(local["matrix"], "rb") as src, mtx_plain.open("wb") as fh:
                    shutil.copyfileobj(src, fh, 8 * 1024 * 1024)
            mat = mmread(mtx_plain).tocsc()
            if mat.shape[0] != len(genes) or mat.shape[1] != len(barcodes):
                raise SystemExit(
                    f"{sample}: mtx {mat.shape} vs genes {len(genes)} barcodes {len(barcodes)}"
                )
            idx = {g: i for i, g in enumerate(genes)}

            def col(g: str) -> np.ndarray:
                if g not in idx:
                    return np.zeros(mat.shape[1], dtype=np.float32)
                return np.asarray(mat[idx[g], :].todense()).ravel()

            epi = np.zeros(mat.shape[1], dtype=bool)
            for g in EPI:
                epi |= col(g) > 0
            epi &= col("PTPRC") == 0
            n_epi = int(epi.sum())
            catalog_rows.append(
                {
                    "patient": sample,
                    "n_cells": int(mat.shape[1]),
                    "n_marker_epithelium": n_epi,
                }
            )
            print(f"  cells={mat.shape[1]} marker_epi={n_epi}", flush=True)
            if n_epi == 0:
                del mat
                continue
            keep_idx = np.flatnonzero(epi)
            frame = pd.DataFrame(
                {
                    "barcode": [barcodes[j] for j in keep_idx],
                    "patient": sample,
                }
            ).set_index(pd.Index(keep_idx), drop=False)
            frame = cap_barcodes(frame, "patient", cap, SEED + i)
            col_idx = frame.index.to_numpy()
            sub = mat[:, col_idx].tocsr()
            del mat
            obs = pd.DataFrame(
                {
                    "dataset": "GSE189357",
                    "patient": sample,
                    "Sample": sample,
                    "unit_id": "GSE189357:" + sample,
                    "patient_id": sample,
                    "histology": "LUAD",
                    "tissue": "TUMOR",
                    "Sample_Origin": "TUMOR",
                    "author_lineage": "marker_epithelium",
                    "author_subtype": "marker_epithelium",
                    "is_normal_tissue": "False",
                    "malig_def": "marker_epi_PTPRC0",
                    "cancer_subtype": "LUAD",
                    "recist": "NA",
                },
                index=["GSE189357:" + sample + ":" + b for b in frame["barcode"].astype(str)],
            )
            adata = ad.AnnData(
                X=sub.T.tocsr().astype(np.float32),
                obs=obs,
                var=pd.DataFrame(index=pd.Index(genes, name="gene")),
            )
            adata.layers["counts"] = adata.X.copy()
            pieces.append(adata)
            print(f"  kept {adata.n_obs}", flush=True)
    if not pieces:
        raise SystemExit("GSE189357: no epithelial cells")
    shared = pieces[0].var_names
    for p in pieces[1:]:
        shared = shared.intersection(p.var_names)
    pieces = [p[:, shared].copy() for p in pieces]
    out = ad.concat(pieces, axis=0, join="inner", merge="same")
    out.obs = _sanitize_obs(out.obs)
    out.layers["counts"] = out.X.copy()
    out.uns["gse189357_catalog"] = catalog_rows
    print(f"GSE189357 written cells={out.n_obs} genes={out.n_vars} units={out.obs['unit_id'].nunique()}", flush=True)
    return out


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

    locked = pd.read_csv(DATA / "GSE205335_patients.tsv", sep="\t")
    keep = set(locked["patient"].astype(str))
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
    cells["patient"] = cells["patient"].astype(str)
    cells["is_normal_tissue"] = cells["tissue"].astype(str).str.startswith(NORMAL_TISSUE_PREFIX)
    cells["is_epi"] = cells["lineage.total"].eq("Epithelial cells")
    cells["in_locked"] = cells["patient"].isin(keep)
    epi = cells.loc[cells["is_epi"] & ~cells["is_normal_tissue"] & cells["in_locked"]].copy()
    catalog = {
        "n_ident": int(len(ident)),
        "n_locked_patients": int(len(keep)),
        "n_epithelial_non_normal_locked": int(len(epi)),
        "by_lineage_sub": epi["lineage.sub"].value_counts().to_dict(),
        "by_tissue": epi["tissue"].value_counts().to_dict(),
        "by_subtype": epi["cancer_subtype"].value_counts().to_dict(),
        "n_patients": int(epi["patient"].nunique()),
        "n_normal_epi_dropped": int((cells["is_epi"] & cells["is_normal_tissue"]).sum()),
        "n_unlocked_epi_dropped": int((cells["is_epi"] & ~cells["in_locked"] & ~cells["is_normal_tissue"]).sum()),
    }
    print(json.dumps({"GSE205335_catalog": catalog}, indent=2), flush=True)
    epi = cap_barcodes(epi.set_index("barcode", drop=False), "patient", cap, SEED + 11)
    keep_bc = set(epi["barcode"].astype(str))
    matrix, genes, barcodes = load_rds_matrix(rds_path)
    col_idx = np.array([i for i, b in enumerate(barcodes) if b in keep_bc], dtype=np.int64)
    print(f"GSE205335 subset columns {len(col_idx)} / {len(barcodes)}", flush=True)
    sub = matrix[:, col_idx].tocsr()
    del matrix
    ordered = [str(barcodes[i]) for i in col_idx]
    obs = epi.set_index("barcode").reindex(ordered)
    adata = ad.AnnData(
        X=sub.T.tocsr().astype(np.float32),
        obs=obs,
        var=pd.DataFrame(index=pd.Index(pd.Index(genes).astype(str).str.upper(), name="gene")),
    )
    adata.obs["dataset"] = "GSE205335"
    adata.obs["histology"] = adata.obs["cancer_subtype"].astype(str)
    adata.obs["Sample"] = adata.obs["patient"].astype(str)
    adata.obs["Sample_Origin"] = adata.obs["tissue"].astype(str)
    adata.obs["unit_id"] = "GSE205335:" + adata.obs["patient"].astype(str)
    adata.obs["patient_id"] = adata.obs["patient"].astype(str)
    adata.obs["author_lineage"] = adata.obs["lineage.total"].astype(str)
    adata.obs["author_subtype"] = adata.obs["lineage.sub"].astype(str)
    adata.obs["malig_def"] = "author_lineage_total_epithelial"
    adata.obs_names = "GSE205335:" + adata.obs_names.astype(str)
    adata.layers["counts"] = adata.X.copy()
    adata.uns["gse205335_catalog"] = catalog
    print(f"GSE205335 written cells={adata.n_obs} genes={adata.n_vars}", flush=True)
    return adata


def concat_shared(a, b):
    import anndata as ad

    # collapse duplicate gene symbols before intersect
    def collapse(ad_obj):
        names = pd.Index(ad_obj.var_names.astype(str).str.upper())
        if names.has_duplicates:
            X = ad_obj.X.tocsr()
            order = []
            seen = {}
            for i, g in enumerate(names):
                if g not in seen:
                    seen[g] = i
                    order.append(i)
            ad_obj = ad_obj[:, order].copy()
            ad_obj.var_names = names[order]
        else:
            ad_obj.var_names = names
        ad_obj.var_names_make_unique()
        return ad_obj

    a = collapse(a)
    b = collapse(b)
    shared = a.var_names.intersection(b.var_names)
    if len(shared) < 5000:
        raise SystemExit(f"too few shared genes: {len(shared)}")
    a2 = a[:, shared].copy()
    b2 = b[:, shared].copy()
    if not str(a2.obs_names[0]).startswith("GSE189357:"):
        a2.obs_names = "GSE189357:" + a2.obs_names.astype(str)
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
    p.add_argument("--data", type=Path, default=Path("/tmp/geo_pair_189357_205335"))
    p.add_argument("--out", type=Path, default=Path("/tmp/geo_pair_189357_205335/epithelium.h5ad"))
    p.add_argument("--scratch", type=Path, default=Path("/tmp/geo_pair_189357_205335/extract"))
    p.add_argument("--cap-per-unit", type=int, default=CAP_PER_UNIT)
    args = p.parse_args()

    cache_a = args.data / "GSE189357_epithelium_cap.h5ad"
    cache_b = args.data / "GSE205335_epithelium_cap.h5ad"
    import anndata as ad

    if cache_a.is_file():
        print(f"reuse {cache_a}", flush=True)
        a = ad.read_h5ad(cache_a)
    else:
        tar189 = args.data / "gse189357" / "GSE189357_RAW.tar"
        if not tar189.exists():
            raise SystemExit(f"missing {tar189}; run download.py")
        a = extract_gse189357(tar189, args.scratch, args.cap_per_unit)
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
    args.out.parent.mkdir(parents=True, exist_ok=True)
    joint.write_h5ad(args.out)
    inv = {
        "n_cells": int(joint.n_obs),
        "n_genes_shared": n_shared,
        "by_dataset": joint.obs["dataset"].value_counts().to_dict(),
        "by_origin": joint.obs["Sample_Origin"].astype(str).value_counts().to_dict(),
        "by_subtype": joint.obs["author_subtype"].astype(str).value_counts().to_dict(),
        "n_units": int(joint.obs["unit_id"].nunique()),
        "n_units_189357": int(joint.obs.loc[joint.obs["dataset"].eq("GSE189357"), "unit_id"].nunique()),
        "n_units_205335": int(joint.obs.loc[joint.obs["dataset"].eq("GSE205335"), "unit_id"].nunique()),
        "cap_per_unit": args.cap_per_unit,
        "gse131907_added": False,
        "dual_high": False,
        "out": str(args.out),
    }
    args.out.with_suffix(".inventory.json").write_text(json.dumps(inv, indent=2))
    print(json.dumps(inv, indent=2), flush=True)


if __name__ == "__main__":
    main()
