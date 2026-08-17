#!/usr/bin/env python3
"""Build a joint epithelial AnnData for GSE123902 + GSE131907.

ADDITIVE. Epithelium only. GSE148071 is not added. No dual-high gate.

GSE131907: author Epithelial cells on tLung / nLung / tL/B / mLN / mBrain
(PE unlabeled epithelium is dropped). nLung AT2 is kept for the Slingshot root.
Unit = Sample.

GSE123902: marker epithelium on PRIMARY_TUMOUR / METASTASIS dense CSVs
(Laughney 2020; author 36.5 GB H5 skipped). NORMAL samples dropped (PR #459).
Unit = donor. Marker gate = (EPCAM|KRT8|KRT18|KRT19)>0 and PTPRC==0.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gene_sets import EPI_MARKERS_123902  # noqa: E402

KEEP_ORIGINS_131907 = ("tLung", "nLung", "tL/B", "mLN", "mBrain")
EPI_TYPE_131907 = "Epithelial cells"
TUMOR_SITES_123902 = ("PRIMARY_TUMOUR", "METASTASIS")
CAP_PER_UNIT = 350
SEED = 1
NAME_RE = re.compile(
    r"(GSM\d+)_(MSK_LX[^_]+(?:B)?)_(PRIMARY_TUMOUR|METASTASIS|NORMAL)_dense\.csv\.gz"
)


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


def extract_gse131907(data: Path, gene_cap: int | None, cap: int) -> "ad.AnnData":
    import anndata as ad

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
    adata.obs["is_normal_tissue"] = adata.obs["Sample_Origin"].eq("nLung")
    adata.obs["donor"] = adata.obs["patient_id"].astype(str)
    adata.layers["counts"] = adata.X.copy()
    print(f"GSE131907 written cells={adata.n_obs} genes={adata.n_vars}", flush=True)
    return adata


def parse_123902_name(fname: str) -> dict:
    m = NAME_RE.match(fname)
    if not m:
        raise ValueError(fname)
    donor = m.group(2).replace("MSK_", "")
    return {"gsm": m.group(1), "patient": donor, "site": m.group(3), "file": fname}


def extract_tar_123902(raw_tar: Path, dest: Path) -> list[Path]:
    dest.mkdir(parents=True, exist_ok=True)
    existing = sorted(dest.glob("GSM*_dense.csv.gz"))
    if existing:
        return existing
    with tarfile.open(raw_tar) as tf:
        tf.extractall(dest)
    return sorted(dest.glob("**/GSM*_dense.csv.gz"))


def load_123902_epithelial(path: Path) -> tuple[list[str], list[str], sparse.csr_matrix, dict]:
    """Stream one dense CSV; keep marker-epithelial cells as sparse counts."""
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split(",")
        genes = header[1:]
        idx = {g: i for i, g in enumerate(genes)}
        marker_i = [idx[g] for g in EPI_MARKERS_123902 if g in idx]
        ptprc_i = idx.get("PTPRC")
        if not marker_i:
            raise SystemExit(f"{path.name}: no epithelial markers in header")
        barcodes: list[str] = []
        rows: list[np.ndarray] = []
        n_total = 0
        n_epi = 0
        for line in handle:
            n_total += 1
            parts = line.rstrip("\n").split(",")
            vals = np.fromstring(",".join(parts[1:]), sep=",", dtype=np.float32)
            if vals.size != len(genes):
                raise SystemExit(f"{path.name}: {vals.size} != {len(genes)}")
            epi = False
            for i in marker_i:
                if vals[i] > 0:
                    epi = True
                    break
            if ptprc_i is not None and vals[ptprc_i] > 0:
                epi = False
            if not epi:
                continue
            n_epi += 1
            barcodes.append(parts[0])
            rows.append(vals)
    if not rows:
        empty = sparse.csr_matrix((0, len(genes)), dtype=np.float32)
        return [], genes, empty, {"n_total": n_total, "n_epithelial": 0}
    mat = sparse.csr_matrix(np.vstack(rows), dtype=np.float32)
    return barcodes, genes, mat, {"n_total": n_total, "n_epithelial": n_epi}


def extract_gse123902(data: Path, cap: int) -> "ad.AnnData":
    import anndata as ad

    raw_tar = data / "gse123902" / "GSE123902_RAW.tar"
    unpacked = data / "gse123902" / "unpacked"
    files = extract_tar_123902(raw_tar, unpacked)
    catalog_rows = []
    adatas = []
    for path in files:
        meta = parse_123902_name(path.name)
        if meta["site"] not in TUMOR_SITES_123902:
            catalog_rows.append({**meta, "kept": False, "reason": "NORMAL dropped"})
            continue
        barcodes, genes, mat, counts = load_123902_epithelial(path)
        print(
            f"GSE123902 {path.name} total={counts['n_total']} epi={counts['n_epithelial']}",
            flush=True,
        )
        catalog_rows.append({**meta, "kept": True, **counts})
        if not barcodes:
            continue
        obs = pd.DataFrame(
            {
                "Index": barcodes,
                "gsm": meta["gsm"],
                "Sample": meta["patient"],
                "patient": meta["patient"],
                "donor": meta["patient"],
                "Sample_Origin": meta["site"],
                "tissue": meta["site"],
                "file": path.name,
            },
            index=barcodes,
        )
        adata = ad.AnnData(
            X=mat,
            obs=obs,
            var=pd.DataFrame(index=pd.Index(genes, name="gene")),
        )
        adatas.append(adata)
    if not adatas:
        raise SystemExit("no GSE123902 epithelial cells after marker gate")
    shared = adatas[0].var_names
    for a in adatas[1:]:
        shared = shared.intersection(a.var_names)
    adatas = [a[:, shared].copy() for a in adatas]
    out = ad.concat(adatas, axis=0, join="inner", merge="same")
    out.obs_names_make_unique()
    out.obs["dataset"] = "GSE123902"
    out.obs["histology"] = "LUAD"
    out.obs["unit_id"] = "GSE123902:" + out.obs["donor"].astype(str)
    out.obs["patient_id"] = out.obs["donor"].astype(str)
    out.obs["author_lineage"] = "marker_epithelial"
    out.obs["author_subtype"] = "marker_epithelial"
    out.obs["Cell_type"] = "marker_epithelial"
    out.obs["Cell_subtype"] = "marker_epithelial"
    out.obs["is_normal_tissue"] = "False"
    catalog = {
        "n_files": int(len(files)),
        "files": catalog_rows,
        "n_epithelial_tumor": int(out.n_obs),
        "n_donors": int(out.obs["donor"].nunique()),
        "by_site": out.obs["Sample_Origin"].value_counts().to_dict(),
        "by_donor": out.obs["donor"].value_counts().to_dict(),
    }
    print(json.dumps({"GSE123902_catalog": catalog}, indent=2, default=str), flush=True)
    out.obs = out.obs.copy()
    out.obs["_protect"] = False
    capped = cap_barcodes(out.obs, "unit_id", cap, SEED + 2)
    out = out[capped.index].copy()
    out.layers["counts"] = out.X.copy()
    print(f"GSE123902 written cells={out.n_obs} genes={out.n_vars}", flush=True)
    (data / "gse123902" / "extract_catalog.json").write_text(
        json.dumps(catalog, indent=2, default=str)
    )
    return out


def _sanitize_obs(obs: pd.DataFrame) -> pd.DataFrame:
    out = obs.copy()
    for col in out.columns:
        if pd.api.types.is_bool_dtype(out[col]) or str(out[col].dtype) == "boolean":
            out[col] = out[col].map({True: "True", False: "False"}).astype(str)
        elif pd.api.types.is_categorical_dtype(out[col]) or out[col].dtype == object:
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
    b2.obs_names = "GSE123902:" + b2.obs_names.astype(str)
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
    p.add_argument("--data", type=Path, default=Path("/tmp/pair_123902_131907"))
    p.add_argument(
        "--out",
        type=Path,
        default=Path("/tmp/pair_123902_131907/epithelium.h5ad"),
    )
    p.add_argument("--cap-per-unit", type=int, default=CAP_PER_UNIT)
    p.add_argument("--gene-cap", type=int, default=None, help="debug only")
    args = p.parse_args()

    cache_a = args.data / "GSE131907_epithelium_cap.h5ad"
    cache_b = args.data / "GSE123902_epithelium_cap.h5ad"
    import anndata as ad

    if cache_a.is_file() and args.gene_cap is None:
        print(f"reuse {cache_a}", flush=True)
        a = ad.read_h5ad(cache_a)
    else:
        a = extract_gse131907(args.data, args.gene_cap, args.cap_per_unit)
        a.obs = _sanitize_obs(a.obs)
        a.write_h5ad(cache_a)
        print(f"cached {cache_a}", flush=True)
    if cache_b.is_file() and args.gene_cap is None:
        print(f"reuse {cache_b}", flush=True)
        b = ad.read_h5ad(cache_b)
    else:
        b = extract_gse123902(args.data, args.cap_per_unit)
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
        "n_units_gse131907": int(
            joint.obs.loc[joint.obs["dataset"] == "GSE131907", "unit_id"].nunique()
        ),
        "n_units_gse123902": int(
            joint.obs.loc[joint.obs["dataset"] == "GSE123902", "unit_id"].nunique()
        ),
        "cap_per_unit": args.cap_per_unit,
        "gse148071_added": False,
        "dual_high": False,
        "out": str(args.out),
    }
    args.out.with_suffix(".inventory.json").write_text(json.dumps(inv, indent=2))
    print(json.dumps(inv, indent=2), flush=True)


if __name__ == "__main__":
    main()
