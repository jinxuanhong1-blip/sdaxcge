#!/usr/bin/env python3
"""Build a joint epithelial AnnData for GSE123902 + GSE205335.

ADDITIVE. Epithelium only. GSE148071 is not added. No dual-high gate.

GSE123902 (Laughney 2020): marker-epithelial cells
  (EPCAM|KRT8|KRT18|KRT19)>0 and PTPRC==0.
  Tumor/met: locked PR #459 donors, PRIMARY preferred.
  NORMAL libraries kept for the AT2-like root (never CLDN4-high).

GSE205335 (Ahn/Lee 2024): author lineage.total == Epithelial cells
  on non-normal tissues.
"""

from __future__ import annotations

import argparse
import gzip
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
ROOT = HERE.parent
DATA = ROOT / "data"

CAP_PER_UNIT = 350
SEED = 1
EPI_MARKERS = ("EPCAM", "KRT8", "KRT18", "KRT19")
NORMAL_TISSUE_PREFIX = ("Normal ",)


def locked_tumor_donors() -> set[str]:
    d = pd.read_csv(DATA / "GSE123902_marker_units.tsv", sep="\t")
    tumor = d[d["tissue"].isin(["PRIMARY", "METASTASIS"])].copy()
    tumor = tumor.sort_values(["patient", "tissue"]).drop_duplicates("patient", keep="first")
    el = tumor[tumor["eligible"].astype(str).str.lower() == "true"]
    return set(el["patient"].astype(str))


def tissue_of(name: str) -> str:
    if "NORMAL" in name:
        return "NORMAL"
    if "METASTASIS" in name:
        return "METASTASIS"
    if "PRIMARY" in name:
        return "PRIMARY"
    return "OTHER"


def marker_epi(df: pd.DataFrame) -> np.ndarray:
    cols = {c.upper(): c for c in df.columns}
    epi = np.zeros(len(df), dtype=bool)
    for g in EPI_MARKERS:
        if g in cols:
            epi |= df[cols[g]].to_numpy(dtype=float) > 0
    ptprc = df[cols["PTPRC"]].to_numpy(dtype=float) if "PTPRC" in cols else np.zeros(len(df))
    return epi & (ptprc == 0)


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
        elif pd.api.types.is_categorical_dtype(out[col]) or out[col].dtype == object:
            out[col] = out[col].astype(str).replace({"nan": "NA", "None": "NA", "<NA>": "NA"})
    return out


def extract_gse123902(data: Path, cap: int) -> "ad.AnnData":
    import anndata as ad

    tar_path = data / "GSE123902" / "GSE123902_RAW.tar"
    keep_tumor = locked_tumor_donors()
    print(f"GSE123902 locked tumor donors n={len(keep_tumor)}", flush=True)

    members = []
    with tarfile.open(tar_path) as tf:
        for m in tf.getmembers():
            if not m.name.endswith(".csv.gz"):
                continue
            name = Path(m.name).name
            parts = name.split("_")
            donor = parts[2] if len(parts) > 2 else name
            tissue = tissue_of(name)
            if tissue == "OTHER":
                continue
            if tissue in {"PRIMARY", "METASTASIS"} and donor not in keep_tumor:
                continue
            members.append((donor, tissue, m, name))
        members.sort(key=lambda x: (x[0], 0 if x[1] == "PRIMARY" else 1 if x[1] == "METASTASIS" else 2))

        chosen_tumor: dict[str, str] = {}
        adatas = []
        catalog_rows = []
        for donor, tissue, m, name in members:
            if tissue in {"PRIMARY", "METASTASIS"}:
                if donor in chosen_tumor:
                    continue
                chosen_tumor[donor] = tissue
            raw = gzip.GzipFile(fileobj=tf.extractfile(m))
            df = pd.read_csv(raw, index_col=0)
            df.columns = [str(c).upper() for c in df.columns]
            df.index = [f"{donor}:{tissue}:{i}" for i in df.index.astype(str)]
            keep = marker_epi(df)
            n_all = int(len(df))
            n_epi = int(keep.sum())
            catalog_rows.append(
                {
                    "donor": donor,
                    "tissue": tissue,
                    "file": name,
                    "n_cells": n_all,
                    "n_marker_epi": n_epi,
                }
            )
            print(f"  {name}: cells={n_all} marker_epi={n_epi}", flush=True)
            if n_epi < 1:
                continue
            sub = df.loc[keep]
            # protect NORMAL SFTPC+ for the root
            sft = sub["SFTPC"].to_numpy(dtype=float) if "SFTPC" in sub.columns else np.zeros(len(sub))
            obs = pd.DataFrame(
                {
                    "dataset": "GSE123902",
                    "patient_id": donor,
                    "Sample": donor,
                    "Sample_Origin": tissue,
                    "tissue": tissue,
                    "histology": "NA",
                    "author_lineage": "marker_epithelial",
                    "author_subtype": "NORMAL_AT2like" if tissue == "NORMAL" else "marker_malignant",
                    "Cell_type": "Epithelial cells",
                    "Cell_subtype": "NORMAL" if tissue == "NORMAL" else "tumor_marker_epi",
                    "is_normal_tissue": "True" if tissue == "NORMAL" else "False",
                    "is_tumor_unit": "False" if tissue == "NORMAL" else "True",
                    "unit_id": (
                        f"GSE123902:{donor}:NORMAL" if tissue == "NORMAL" else f"GSE123902:{donor}"
                    ),
                    "sftpc_count": sft,
                    "library": name,
                },
                index=sub.index.astype(str),
            )
            X = sparse.csr_matrix(sub.to_numpy(dtype=np.float32))
            adatas.append(
                ad.AnnData(
                    X=X,
                    obs=obs,
                    var=pd.DataFrame(index=pd.Index(sub.columns.astype(str), name="gene")),
                )
            )
            del df, sub

    if not adatas:
        raise SystemExit("GSE123902: no marker-epithelial cells")
    genes = adatas[0].var_names
    for a in adatas[1:]:
        genes = genes.intersection(a.var_names)
    adatas = [a[:, genes].copy() for a in adatas]
    adata = ad.concat(adatas, axis=0, join="inner", merge="same")
    del adatas
    protect = (adata.obs["Sample_Origin"] == "NORMAL") & (adata.obs["sftpc_count"].astype(float) > 0)
    adata.obs = adata.obs.copy()
    before = int(adata.n_obs)
    kept = cap_barcodes(adata.obs, "unit_id", cap, SEED, protect)
    adata = adata[kept.index].copy()
    adata.layers["counts"] = adata.X.copy()
    inv = {
        "n_before_cap": before,
        "n_after_cap": int(adata.n_obs),
        "by_tissue": adata.obs["Sample_Origin"].value_counts().to_dict(),
        "n_tumor_donors": int((adata.obs["is_tumor_unit"] == "True").sum() and adata.obs.loc[adata.obs["is_tumor_unit"] == "True", "patient_id"].nunique()),
        "n_normal_units": int(adata.obs.loc[adata.obs["Sample_Origin"] == "NORMAL", "unit_id"].nunique()),
        "catalog": catalog_rows,
    }
    print(json.dumps({"GSE123902_extract": {k: v for k, v in inv.items() if k != "catalog"}}, indent=2), flush=True)
    return adata


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
    genes_u = pd.Index([str(g).upper() for g in genes], name="gene")
    # collapse duplicate gene symbols after uppercasing
    if genes_u.has_duplicates:
        print(f"GSE205335 collapsing {int(genes_u.duplicated().sum())} duplicate symbols", flush=True)
    adata = ad.AnnData(
        X=sub.T.tocsr(),
        obs=obs,
        var=pd.DataFrame(index=genes_u),
    )
    if adata.var_names.has_duplicates:
        adata = adata[:, ~adata.var_names.duplicated()].copy()
    adata.obs["dataset"] = "GSE205335"
    adata.obs["histology"] = adata.obs["cancer_subtype"].astype(str)
    adata.obs["Sample"] = adata.obs["patient"].astype(str)
    adata.obs["Sample_Origin"] = adata.obs["tissue"].astype(str)
    adata.obs["tissue"] = adata.obs["tissue"].astype(str)
    adata.obs["Cell_type"] = "Epithelial cells"
    adata.obs["Cell_subtype"] = adata.obs["lineage.sub"].astype(str)
    adata.obs["unit_id"] = "GSE205335:" + adata.obs["patient"].astype(str)
    adata.obs["patient_id"] = adata.obs["patient"].astype(str)
    adata.obs["author_lineage"] = adata.obs["lineage.total"].astype(str)
    adata.obs["author_subtype"] = adata.obs["lineage.sub"].astype(str)
    adata.obs["is_normal_tissue"] = "False"
    adata.obs["is_tumor_unit"] = "True"
    adata.obs["sftpc_count"] = 0.0
    adata.obs["library"] = adata.obs["orig.ident"].astype(str)
    adata.layers["counts"] = adata.X.copy()
    print(f"GSE205335 written cells={adata.n_obs} genes={adata.n_vars}", flush=True)
    return adata


def concat_shared(a, b):
    import anndata as ad

    a.var_names = pd.Index([str(g).upper() for g in a.var_names], name="gene")
    b.var_names = pd.Index([str(g).upper() for g in b.var_names], name="gene")
    if a.var_names.has_duplicates:
        a = a[:, ~a.var_names.duplicated()].copy()
    if b.var_names.has_duplicates:
        b = b[:, ~b.var_names.duplicated()].copy()
    shared = a.var_names.intersection(b.var_names)
    if len(shared) < 5000:
        raise SystemExit(f"too few shared genes: {len(shared)}")
    a2 = a[:, shared].copy()
    b2 = b[:, shared].copy()
    a2.obs_names = "GSE123902:" + a2.obs_names.astype(str)
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
    p.add_argument("--data", type=Path, default=Path("/tmp/pair_123902_205335_traj"))
    p.add_argument(
        "--out",
        type=Path,
        default=Path("/tmp/pair_123902_205335_traj/epithelium.h5ad"),
    )
    p.add_argument("--cap-per-unit", type=int, default=CAP_PER_UNIT)
    args = p.parse_args()

    cache_a = args.data / "GSE123902_epithelium_cap.h5ad"
    cache_b = args.data / "GSE205335_epithelium_cap.h5ad"
    import anndata as ad

    if cache_a.is_file():
        print(f"reuse {cache_a}", flush=True)
        a = ad.read_h5ad(cache_a)
    else:
        a = extract_gse123902(args.data, args.cap_per_unit)
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
        "n_tumor_units": int(joint.obs.loc[joint.obs["is_tumor_unit"] == "True", "unit_id"].nunique()),
        "cap_per_unit": args.cap_per_unit,
        "gse148071_added": False,
        "out": str(args.out),
    }
    args.out.with_suffix(".inventory.json").write_text(json.dumps(inv, indent=2))
    print(json.dumps(inv, indent=2), flush=True)


if __name__ == "__main__":
    main()
