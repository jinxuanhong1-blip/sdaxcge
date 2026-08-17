#!/usr/bin/env python3
"""Extract marker-epithelial cells for GSE123902 + GSE189357.

Same gate as PR #459: (EPCAM|KRT8|KRT18|KRT19)>0 and PTPRC==0.

Given inferential units (not re-audited): 13 GSE123902 tumor/met donors +
9 GSE189357 patients = 22. GSE123902 NORMAL libraries are kept only as a
root pool (AT2-like, never CLDN4-high). No dual-high gate. Cap ≤350/unit.
"""
from __future__ import annotations

import argparse
import gzip
import json
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse
from scipy.io import mmread

HERE = Path(__file__).resolve().parents[1]
DATA = HERE / "data"
EPI = ["EPCAM", "KRT8", "KRT18", "KRT19"]
CAP = 350
RNG = np.random.default_rng(22)


def _upper_map(names: list[str]) -> dict[str, int]:
    return {n.upper(): i for i, n in enumerate(names)}


def _mal_mask_from_panel(panel: np.ndarray, names: list[str]) -> np.ndarray:
    idx = {g: i for i, g in enumerate(names)}

    def col(g: str) -> np.ndarray:
        if g not in idx:
            return np.zeros(panel.shape[0], dtype=float)
        return panel[:, idx[g]].astype(float)

    epi = np.zeros(panel.shape[0], dtype=bool)
    for g in EPI:
        epi |= col(g) > 0
    return epi & (col("PTPRC") == 0)


def _cap_indices(n: int, cap: int, protect: np.ndarray | None = None) -> np.ndarray:
    idx = np.arange(n)
    if n <= cap:
        return idx
    if protect is not None and protect.any():
        prot = idx[protect]
        rest = idx[~protect]
        n_keep_rest = max(0, cap - prot.size)
        if prot.size >= cap:
            return RNG.choice(prot, size=cap, replace=False)
        picked = RNG.choice(rest, size=n_keep_rest, replace=False) if n_keep_rest else np.array([], dtype=int)
        return np.sort(np.concatenate([prot, picked]))
    return np.sort(RNG.choice(idx, size=cap, replace=False))


def _dedup_genes(mat: sparse.spmatrix, genes: list[str]) -> tuple[sparse.spmatrix, list[str]]:
    """Sum duplicate gene symbols (cells x genes CSR)."""
    ser = pd.Series(np.arange(len(genes)), index=genes)
    if not ser.index.duplicated().any():
        return mat.tocsr(), genes
    groups = ser.groupby(level=0).apply(lambda s: list(s.values))
    uniq = list(groups.index)
    cols = []
    for g in uniq:
        ix = groups[g]
        if len(ix) == 1:
            cols.append(mat[:, ix[0]])
        else:
            cols.append(mat[:, ix].sum(axis=1))
    out = sparse.hstack(cols, format="csr")
    return out, uniq


def extract_gse123902(tar_path: Path, units: pd.DataFrame, cap: int) -> tuple[list[AnnData], list[dict]]:
    tumor = units[units["tissue"].isin(["PRIMARY", "METASTASIS"])].copy()
    tumor = tumor.sort_values(["patient", "tissue"]).drop_duplicates("patient", keep="first")
    keep_tumor = set(tumor.loc[tumor["eligible"] == True, "file"].astype(str))
    donor_of = dict(zip(tumor["file"].astype(str), tumor["patient"].astype(str)))
    tissue_of = dict(zip(tumor["file"].astype(str), tumor["tissue"].astype(str)))

    adatas: list[AnnData] = []
    catalog: list[dict] = []
    with tarfile.open(tar_path) as tf:
        for m in tf.getmembers():
            if not m.name.endswith(".csv.gz"):
                continue
            name = Path(m.name).name
            parts = name.split("_")
            donor = parts[2] if len(parts) > 2 else name
            tissue = (
                "NORMAL"
                if "NORMAL" in name
                else (
                    "METASTASIS"
                    if "METASTASIS" in name
                    else ("PRIMARY" if "PRIMARY" in name else "OTHER")
                )
            )
            is_given = name in keep_tumor
            is_root = tissue == "NORMAL"
            if not is_given and not is_root:
                print(f"  skip {name}", flush=True)
                continue
            print(f"  reading {name}", flush=True)
            raw = gzip.GzipFile(fileobj=tf.extractfile(m))
            df = pd.read_csv(raw, index_col=0)
            df.columns = [str(c).upper() for c in df.columns]
            if df.columns.duplicated().any():
                df = df.T.groupby(level=0).sum().T
            genes = [str(c) for c in df.columns]
            mat = sparse.csr_matrix(df.to_numpy(dtype=np.float32))
            panel_names = [g for g in ["PTPRC", *EPI, "SFTPC", "CLDN4"] if g in genes]
            gix = {g: i for i, g in enumerate(genes)}
            panel = np.zeros((mat.shape[0], len(panel_names)), dtype=np.float32)
            for j, g in enumerate(panel_names):
                panel[:, j] = np.asarray(mat[:, gix[g]].todense()).ravel()
            mal = _mal_mask_from_panel(panel, panel_names)
            n_mal = int(mal.sum())
            rec = {
                "cohort": "GSE123902",
                "file": name,
                "patient": donor,
                "tissue": tissue,
                "role": "given" if is_given else "root_pool",
                "n_cells_library": int(mat.shape[0]),
                "n_marker_epithelial": n_mal,
                "n_genes": int(mat.shape[1]),
                "n_kept": 0,
            }
            if n_mal == 0:
                catalog.append(rec)
                print(f"  WARN {name}: 0 marker-epithelial", flush=True)
                continue
            sub = mat[mal]
            barcodes = df.index.astype(str).to_numpy()[mal]
            protect = None
            if is_root and "SFTPC" in panel_names:
                sft = panel[mal, panel_names.index("SFTPC")]
                protect = sft > 0
            keep = _cap_indices(sub.shape[0], cap, protect)
            sub = sub[keep]
            barcodes = barcodes[keep]
            rec["n_kept"] = int(sub.shape[0])
            catalog.append(rec)
            obs = pd.DataFrame(
                {
                    "dataset": "GSE123902",
                    "patient": donor,
                    "tissue": tissue,
                    "role": "given" if is_given else "root_pool",
                    "unit_id": f"GSE123902:{donor}" if is_given else f"GSE123902:{donor}_NORMAL",
                    "file": name,
                    "barcode": barcodes,
                }
            )
            obs.index = [f"GSE123902:{donor}:{tissue}:{b}" for b in barcodes]
            var = pd.DataFrame(index=pd.Index(genes, name="gene"))
            adatas.append(AnnData(X=sub, obs=obs, var=var))
            print(
                f"  {name}: lib={mat.shape[0]} mal={n_mal} kept={sub.shape[0]} role={rec['role']}",
                flush=True,
            )
            del df, mat, sub
    return adatas, catalog


def _read_10x_features(tf: tarfile.TarFile, members: dict, sample: str) -> list[str]:
    name = next(n for n in members if f"_{sample}_features" in n or f"_{sample}_genes" in n)
    genes = []
    with gzip.GzipFile(fileobj=tf.extractfile(members[name])) as f:
        for line in f:
            p = line.decode().strip().split("\t")
            genes.append((p[1] if len(p) > 1 else p[0]).upper())
    return genes


def _read_10x_barcodes(tf: tarfile.TarFile, members: dict, sample: str) -> list[str]:
    name = next(n for n in members if f"_{sample}_barcodes" in n)
    out = []
    with gzip.GzipFile(fileobj=tf.extractfile(members[name])) as f:
        for line in f:
            out.append(line.decode().strip().split("\t")[0])
    return out


def _extract_mtx(tf: tarfile.TarFile, members: dict, sample: str, dest: Path) -> Path:
    name = next(n for n in members if f"_{sample}_matrix.mtx" in n)
    dest.mkdir(parents=True, exist_ok=True)
    out = dest / f"{sample}_matrix.mtx"
    if out.exists() and out.stat().st_size > 1_000_000:
        return out
    with gzip.GzipFile(fileobj=tf.extractfile(members[name])) as src, out.open("wb") as fh:
        while True:
            chunk = src.read(8 * 1024 * 1024)
            if not chunk:
                break
            fh.write(chunk)
    return out


def extract_gse189357(
    tar_path: Path, units: pd.DataFrame, scratch: Path, cap: int
) -> tuple[list[AnnData], list[dict]]:
    keep = set(units.loc[units["eligible"] == True, "patient"].astype(str))
    adatas: list[AnnData] = []
    catalog: list[dict] = []
    with tarfile.open(tar_path) as tf:
        members = {m.name: m for m in tf.getmembers()}
        for sample in [f"TD{i}" for i in range(1, 10)]:
            if sample not in keep:
                print(f"  skip {sample}", flush=True)
                continue
            print(f"  reading {sample}", flush=True)
            genes = _read_10x_features(tf, members, sample)
            barcodes = _read_10x_barcodes(tf, members, sample)
            mtx_path = _extract_mtx(tf, members, sample, scratch)
            mat = mmread(mtx_path).tocsc()  # genes x cells
            if mat.shape[0] != len(genes) or mat.shape[1] != len(barcodes):
                raise RuntimeError(
                    f"{sample}: features {len(genes)} barcodes {len(barcodes)} mtx {mat.shape}"
                )
            idx = _upper_map(genes)
            panel_names = [g for g in ["PTPRC", *EPI] if g in idx]
            cell_panel = np.zeros((mat.shape[1], len(panel_names)), dtype=np.float32)
            for j, g in enumerate(panel_names):
                cell_panel[:, j] = np.asarray(mat[idx[g], :].todense()).ravel()
            mal = _mal_mask_from_panel(cell_panel, panel_names)
            n_mal = int(mal.sum())
            rec = {
                "cohort": "GSE189357",
                "file": f"{sample}_matrix.mtx",
                "patient": sample,
                "tissue": "TUMOR",
                "role": "given",
                "n_cells_library": int(mat.shape[1]),
                "n_marker_epithelial": n_mal,
                "n_genes": int(mat.shape[0]),
                "n_kept": 0,
            }
            if n_mal == 0:
                catalog.append(rec)
                print(f"  WARN {sample}: 0 marker-epithelial", flush=True)
                continue
            keep_i = np.flatnonzero(mal)
            keep_i = keep_i[_cap_indices(keep_i.size, cap)]
            sub = mat[:, keep_i].T.tocsr().astype(np.float32)  # cells x genes
            sub, genes_u = _dedup_genes(sub, genes)
            rec["n_kept"] = int(sub.shape[0])
            catalog.append(rec)
            bc = [barcodes[i] for i in keep_i]
            obs = pd.DataFrame(
                {
                    "dataset": "GSE189357",
                    "patient": sample,
                    "tissue": "TUMOR",
                    "role": "given",
                    "unit_id": f"GSE189357:{sample}",
                    "file": rec["file"],
                    "barcode": bc,
                }
            )
            obs.index = [f"GSE189357:{sample}:{b}" for b in bc]
            var = pd.DataFrame(index=pd.Index(genes_u, name="gene"))
            adatas.append(AnnData(X=sub, obs=obs, var=var))
            print(f"  {sample}: lib={mat.shape[1]} mal={n_mal} kept={sub.shape[0]}", flush=True)
            del mat, sub
    return adatas, catalog


def concat_inner(adatas: list[AnnData]) -> AnnData:
    if not adatas:
        raise SystemExit("no cells extracted")
    common = set(adatas[0].var_names)
    for ad in adatas[1:]:
        common &= set(ad.var_names)
    common = sorted(common)
    if "CLDN4" not in common:
        raise SystemExit("CLDN4 missing from gene intersection")
    parts = [ad[:, common].copy() for ad in adatas]
    from anndata import concat

    out = concat(parts, join="inner", index_unique=None)
    out.var_names_make_unique()
    return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--tars", type=Path, default=Path("/tmp/geo_pair_123902_189357"))
    p.add_argument("--scratch", type=Path, default=Path("/tmp/geo_pair_123902_189357/extract"))
    p.add_argument("--out", type=Path, default=Path("/tmp/geo_pair_123902_189357/epithelium.h5ad"))
    p.add_argument("--cap", type=int, default=CAP)
    args = p.parse_args()

    u123 = pd.read_csv(DATA / "GSE123902_marker_units.tsv", sep="\t")
    u189 = pd.read_csv(DATA / "GSE189357_marker_units.tsv", sep="\t")
    p123 = args.tars / "GSE123902_RAW.tar"
    p189 = args.tars / "GSE189357_RAW.tar"
    if not p123.exists() or not p189.exists():
        raise SystemExit(f"missing tars in {args.tars}; run download.py")

    print("Extracting GSE123902 (given tumor/met + NORMAL root pool)", flush=True)
    a123, c123 = extract_gse123902(p123, u123, args.cap)
    print("Extracting GSE189357 (given patients)", flush=True)
    a189, c189 = extract_gse189357(p189, u189, args.scratch, args.cap)

    adata = concat_inner(a123 + a189)
    adata.obs["cell_id"] = adata.obs_names.astype(str)
    catalog = c123 + c189
    args.out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.out)
    inv = {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_given_units": int(adata.obs.loc[adata.obs["role"] == "given", "unit_id"].nunique()),
        "n_root_pool_units": int(adata.obs.loc[adata.obs["role"] == "root_pool", "unit_id"].nunique()),
        "cap_per_unit": args.cap,
        "gate": "(EPCAM|KRT8|KRT18|KRT19)>0 and PTPRC==0",
        "catalog": catalog,
    }
    inv_path = args.out.with_suffix(".inventory.json")
    inv_path.write_text(json.dumps(inv, indent=2))
    print(
        f"wrote {args.out} cells={adata.n_obs} genes={adata.n_vars} "
        f"given_units={inv['n_given_units']} root_units={inv['n_root_pool_units']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
