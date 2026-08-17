#!/usr/bin/env python3
"""Build a GSE205335-only epithelial AnnData (public processed UMI).

Author epithelium on non-normal tissues. Normal Lung / LN / Brain dropped.
GSE148071 and GSE131907 are not used. No dual-high gate.
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
CAP_PER_PATIENT = 500
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


def recist_arm(value: str) -> str:
    recist = str(value).strip().upper()
    if recist == "PR":
        return "PR"
    if recist in {"PD", "SD"}:
        return "PD_SD"
    return "NE_or_other"


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


def _sanitize_obs(obs: pd.DataFrame) -> pd.DataFrame:
    out = obs.copy()
    for col in out.columns:
        if pd.api.types.is_bool_dtype(out[col]) or str(out[col].dtype) == "boolean":
            out[col] = out[col].map({True: "True", False: "False"}).astype(str)
        elif str(out[col].dtype) == "object" or pd.api.types.is_string_dtype(out[col]):
            out[col] = out[col].astype(str).replace({"nan": "NA", "None": "NA", "<NA>": "NA"})
    return out


def extract_gse205335(data: Path, cap: int) -> "ad.AnnData":
    import anndata as ad

    ident_path = data / "GSE205335_Lung_IO_CellIdentity.txt.gz"
    soft_path = data / "GSE205335_family.soft.gz"
    rds_path = data / "GSE205335_Lung_IO_UMI_matrix.rds.gz"
    ident = pd.read_csv(ident_path, sep="\t")
    meta = parse_geo_soft(soft_path)
    keep_cols = [
        c
        for c in (
            "orig.ident",
            "gsm",
            "patient",
            "tissue",
            "recist",
            "cancer_subtype",
            "tumor_stage",
        )
        if c in meta.columns
    ]
    cells = ident.merge(meta[keep_cols], on="orig.ident", how="left", validate="many_to_one")
    if cells["patient"].isna().any():
        missing = cells.loc[cells["patient"].isna(), "orig.ident"].unique().tolist()
        raise SystemExit(f"GSE205335 identity samples missing GEO metadata: {missing}")
    cells["is_normal_tissue"] = cells["tissue"].astype(str).str.startswith(NORMAL_TISSUE_PREFIX)
    cells["is_epi"] = cells["lineage.total"].eq("Epithelial cells")
    cells["recist_arm"] = cells["recist"].map(recist_arm)
    epi = cells.loc[cells["is_epi"] & ~cells["is_normal_tissue"]].copy()
    catalog = {
        "n_ident": int(len(ident)),
        "n_epithelial_non_normal": int(len(epi)),
        "by_lineage_sub": epi["lineage.sub"].value_counts().to_dict(),
        "by_tissue": epi["tissue"].value_counts().to_dict(),
        "by_recist": epi.drop_duplicates("patient")["recist"].value_counts().to_dict(),
        "by_recist_arm": epi.drop_duplicates("patient")["recist_arm"].value_counts().to_dict(),
        "by_histology": epi.drop_duplicates("patient")["cancer_subtype"].value_counts().to_dict(),
        "n_patients": int(epi["patient"].nunique()),
        "n_normal_epi_dropped": int((cells["is_epi"] & cells["is_normal_tissue"]).sum()),
        "cells_per_patient": epi.groupby("patient").size().to_dict(),
    }
    print(json.dumps({"GSE205335_catalog": catalog}, indent=2), flush=True)
    epi = cap_barcodes(epi.set_index("barcode", drop=False), "patient", cap, SEED)
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
    adata.obs["dataset"] = "GSE205335"
    adata.obs["histology"] = adata.obs["cancer_subtype"].astype(str)
    adata.obs["patient_id"] = adata.obs["patient"].astype(str)
    adata.obs["author_lineage"] = adata.obs["lineage.total"].astype(str)
    adata.obs["author_subtype"] = adata.obs["lineage.sub"].astype(str)
    adata.obs["is_author_malignant"] = adata.obs["author_subtype"].eq("Malignant cells")
    adata.layers["counts"] = adata.X.copy()
    adata.obs = _sanitize_obs(adata.obs)
    print(f"GSE205335 written cells={adata.n_obs} genes={adata.n_vars}", flush=True)
    return adata, catalog


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=Path, default=Path("/tmp/gse205335_traj"))
    p.add_argument("--out", type=Path, default=Path("/tmp/gse205335_traj/epithelium.h5ad"))
    p.add_argument("--cap-per-patient", type=int, default=CAP_PER_PATIENT)
    args = p.parse_args()

    import anndata as ad  # noqa: F401

    adata, catalog = extract_gse205335(args.data, args.cap_per_patient)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.out)
    inv = {
        "accession": "GSE205335",
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_patients": int(adata.obs["patient_id"].nunique()),
        "by_subtype": adata.obs["author_subtype"].value_counts().to_dict(),
        "by_histology_cells": adata.obs["histology"].value_counts().to_dict(),
        "by_recist_arm_cells": adata.obs["recist_arm"].value_counts().to_dict(),
        "cap_per_patient": args.cap_per_patient,
        "gse148071_used": False,
        "gse131907_used": False,
        "dual_high": False,
        "catalog": catalog,
        "out": str(args.out),
    }
    args.out.with_suffix(".inventory.json").write_text(json.dumps(inv, indent=2))
    print(json.dumps({k: v for k, v in inv.items() if k != "catalog"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
