#!/usr/bin/env python3
"""Build a malignant-only AnnData for GSE205335 (Ahn/Lee ICI).

Author lineage.sub == Malignant cells on non-normal tissues.
Normal Lung / LN / Brain dropped. Patient is the unit. Cap after catalog.
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
CAP_PER_PATIENT = 400
SEED = 1
RESPONSE_MAP = {"PR": "R", "CR": "R", "SD": "NR", "PD": "NR", "NE": "NE"}


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


def gunzip_until_rds(src: Path, dest: Path) -> Path:
    """GEO ships a double-gzipped RDS. Peel gzip until the R XDR magic."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        with dest.open("rb") as handle:
            magic = handle.read(2)
        if magic == b"X\n":
            print(f"RDS ready {dest} ({dest.stat().st_size} bytes)", flush=True)
            return dest
    current = src
    tmp_dir = dest.parent
    for i in range(4):
        with current.open("rb") as handle:
            magic = handle.read(2)
        if magic == b"\x1f\x8b":
            nxt = tmp_dir / f"{dest.name}.peel{i}"
            print(f"gunzip peel {i}: {current}", flush=True)
            with gzip.open(current, "rb") as source, nxt.open("wb") as out:
                shutil.copyfileobj(source, out, 16 * 1024 * 1024)
            if current != src and current.exists():
                current.unlink()
            current = nxt
            continue
        break
    if current != dest:
        current.replace(dest)
    print(f"RDS ready {dest} ({dest.stat().st_size} bytes)", flush=True)
    return dest


def load_rds_matrix(path: Path):
    import rdata

    with tempfile.TemporaryDirectory(prefix="gse205335-") as tmp:
        rds_path = gunzip_until_rds(path, Path(tmp) / "matrix.rds")
        print("read GSE205335 RDS", flush=True)
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", message='Missing constructor for R class "dgCMatrix"'
            )
            obj = rdata.read_rds(rds_path)
    expected = {"i", "p", "Dim", "Dimnames", "x"}
    if not expected.issubset(vars(obj)):
        raise TypeError("RDS is not the expected Matrix::dgCMatrix object")
    genes = np.asarray(obj.Dimnames[0], dtype=str)
    barcodes = np.asarray(obj.Dimnames[1], dtype=str)
    print(f"build CSC {tuple(obj.Dim)}", flush=True)
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False)
    return matrix, genes, barcodes


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
        elif str(out[col].dtype) == "object" or pd.api.types.is_categorical_dtype(out[col]):
            out[col] = out[col].astype(str).replace({"nan": "NA", "None": "NA", "<NA>": "NA"})
    return out


def extract_malignant(data: Path, cap: int) -> "ad.AnnData":
    import anndata as ad

    ident_path = data / "GSE205335_Lung_IO_CellIdentity.txt.gz"
    soft_path = data / "GSE205335_family.soft.gz"
    rds_path = data / "GSE205335_Lung_IO_UMI_matrix.rds.gz"
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
                "platform",
            ]
        ],
        on="orig.ident",
        how="left",
        validate="many_to_one",
    )
    if cells["patient"].isna().any():
        missing = cells.loc[cells["patient"].isna(), "orig.ident"].unique().tolist()
        raise SystemExit(f"identity samples missing GEO metadata: {missing}")
    cells["is_normal_tissue"] = cells["tissue"].astype(str).str.startswith(NORMAL_TISSUE_PREFIX)
    cells["is_malignant"] = cells["lineage.sub"].eq("Malignant cells")
    mal = cells.loc[cells["is_malignant"] & ~cells["is_normal_tissue"]].copy()
    mal["response"] = mal["recist"].map(RESPONSE_MAP).fillna("NE")
    catalog = {
        "n_ident": int(len(ident)),
        "n_geo_patients": int(cells["patient"].nunique()),
        "n_malignant_non_normal": int(len(mal)),
        "n_malignant_normal_dropped": int(
            (cells["is_malignant"] & cells["is_normal_tissue"]).sum()
        ),
        "by_tissue": mal["tissue"].value_counts().to_dict(),
        "by_histology": mal["cancer_subtype"].value_counts().to_dict(),
        "by_recist": mal.groupby("patient")["recist"].first().value_counts().to_dict(),
        "n_patients_malignant": int(mal["patient"].nunique()),
        "per_patient_catalog": mal.groupby("patient").size().to_dict(),
        "patients_dropped_normal_only": sorted(
            set(cells["patient"].astype(str)) - set(mal["patient"].astype(str))
        ),
    }
    print(json.dumps({"GSE205335_malignant_catalog": catalog}, indent=2), flush=True)
    mal = cap_barcodes(mal.set_index("barcode", drop=False), "patient", cap, SEED)
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
    adata.obs["patient_id"] = adata.obs["patient"].astype(str)
    adata.obs["unit_id"] = adata.obs["patient"].astype(str)
    adata.obs["author_lineage"] = adata.obs["lineage.total"].astype(str)
    adata.obs["author_subtype"] = adata.obs["lineage.sub"].astype(str)
    adata.obs["recist"] = adata.obs["recist"].astype(str)
    adata.obs["response"] = adata.obs["response"].astype(str)
    adata.layers["counts"] = adata.X.copy()
    adata.obs = _sanitize_obs(adata.obs)
    print(f"GSE205335 malignant written cells={adata.n_obs} genes={adata.n_vars}", flush=True)
    return adata, catalog


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=Path, default=Path("/tmp/gse205335_slingshot_real"))
    p.add_argument(
        "--out",
        type=Path,
        default=Path("/tmp/gse205335_slingshot_real/malignant.h5ad"),
    )
    p.add_argument("--cap-per-patient", type=int, default=CAP_PER_PATIENT)
    args = p.parse_args()

    adata, catalog = extract_malignant(args.data, args.cap_per_patient)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.out)
    inv = {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_patients": int(adata.obs["patient_id"].nunique()),
        "cap_per_patient": args.cap_per_patient,
        "by_histology_cells": adata.obs["histology"].value_counts().to_dict(),
        "by_recist_patients": adata.obs.groupby("patient_id")["recist"]
        .first()
        .value_counts()
        .to_dict(),
        "by_response_patients": adata.obs.groupby("patient_id")["response"]
        .first()
        .value_counts()
        .to_dict(),
        "per_patient_analysis": adata.obs.groupby("patient_id").size().to_dict(),
        "catalog": catalog,
        "out": str(args.out),
        "compartment": "author lineage.sub == Malignant cells",
        "normal_tissues_dropped": True,
        "dual_high": False,
        "gse131907_added": False,
        "gse207422_added": False,
    }
    args.out.with_suffix(".inventory.json").write_text(json.dumps(inv, indent=2))
    print(json.dumps({k: v for k, v in inv.items() if k != "catalog"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
