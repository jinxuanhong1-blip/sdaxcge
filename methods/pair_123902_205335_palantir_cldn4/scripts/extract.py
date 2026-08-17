#!/usr/bin/env python3
"""Build a joint epithelial AnnData for GSE123902 + GSE205335.

ADDITIVE. CLDN4-only. No GSE148071. No dual-high gate.
GSE123902: marker epithelium (EPCAM|KRT8|KRT18|KRT19)>0 and PTPRC==0.
  Tumor/met: one library per donor (PRIMARY preferred).
  NORMAL libraries kept only as the Palantir root pool (not CLDN4-high).
GSE205335: author lineage.total == Epithelial cells; normal tissues dropped.
Donor (GSE123902) / patient (GSE205335) is the unit.
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

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE / "scripts"))
from gene_sets import EPI_POS  # noqa: E402

CAP_PER_UNIT = 250
SEED = 1
NORMAL_TISSUE_PREFIX = ("Normal ",)


def tissue_of(name: str) -> str:
    u = name.upper()
    if "NORMAL" in u:
        return "NORMAL"
    if "METASTASIS" in u:
        return "METASTASIS"
    if "PRIMARY" in u:
        return "PRIMARY"
    return "OTHER"


def donor_of(name: str) -> str:
    parts = Path(name).name.split("_")
    # GSM3516662_MSK_LX653_PRIMARY_TUMOUR_dense.csv.gz
    if len(parts) > 2:
        return parts[2]
    return Path(name).name


def marker_epi(df: pd.DataFrame) -> np.ndarray:
    cols = {str(c).upper(): c for c in df.columns}
    epi = np.zeros(len(df), dtype=bool)
    for g in EPI_POS:
        if g in cols:
            epi |= pd.to_numeric(df[cols[g]], errors="coerce").fillna(0).to_numpy() > 0
    ptprc = (
        pd.to_numeric(df[cols["PTPRC"]], errors="coerce").fillna(0).to_numpy()
        if "PTPRC" in cols
        else np.zeros(len(df))
    )
    return epi & (ptprc == 0)


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


def _to_sparse_cells_x_genes(df: pd.DataFrame) -> tuple[sparse.csr_matrix, list[str], list[str]]:
    """Laughney dense CSVs are cells × genes."""
    genes = [str(c).upper() for c in df.columns]
    cells = [str(i) for i in df.index]
    mat = sparse.csr_matrix(df.to_numpy(dtype=np.float32))
    return mat, cells, genes


def extract_gse123902(data: Path, cap: int) -> "ad.AnnData":
    import anndata as ad

    tar_path = data / "GSE123902" / "GSE123902_RAW.tar"
    units = pd.read_csv(HERE / "data" / "GSE123902_marker_units.tsv", sep="\t")
    units["patient"] = units["patient"].astype(str)
    tumor = units[units["tissue"].isin(["PRIMARY", "METASTASIS"])].copy()
    tumor = tumor.sort_values(["patient", "tissue"]).drop_duplicates("patient", keep="first")
    tumor = tumor[tumor["eligible"].astype(str).str.lower() == "true"]
    keep_tumor = set(tumor["patient"].astype(str))
    normals = units[units["tissue"].eq("NORMAL")].copy()
    print(
        json.dumps(
            {
                "GSE123902_locked": {
                    "n_tumor_donors": int(len(keep_tumor)),
                    "tumor_donors": sorted(keep_tumor),
                    "n_normal_libs": int(len(normals)),
                    "normal_donors": sorted(normals["patient"].astype(str).unique()),
                }
            },
            indent=2,
        ),
        flush=True,
    )

    blocks = []
    with tarfile.open(tar_path) as tf:
        members = []
        for m in tf.getmembers():
            if not m.name.endswith(".csv.gz"):
                continue
            name = Path(m.name).name
            donor = donor_of(name)
            tissue = tissue_of(name)
            if tissue == "NORMAL":
                members.append((donor, tissue, m, name, "root_pool"))
            elif tissue in {"PRIMARY", "METASTASIS"} and donor in keep_tumor:
                members.append((donor, tissue, m, name, "tumor"))
        members.sort(key=lambda x: (x[0], 0 if x[1] == "PRIMARY" else 1 if x[1] == "METASTASIS" else 2))

        chosen_tumor: dict[str, str] = {}
        for donor, tissue, m, name, role in members:
            if role == "tumor":
                if donor in chosen_tumor:
                    continue
                chosen_tumor[donor] = tissue
            print(f"  read {name} role={role}", flush=True)
            raw = gzip.GzipFile(fileobj=tf.extractfile(m))
            df = pd.read_csv(raw, index_col=0)
            epi = marker_epi(df)
            n_epi = int(epi.sum())
            if n_epi < 1:
                print(f"    SKIP {name}: 0 marker-epi", flush=True)
                continue
            sub = df.loc[epi].copy()
            sub.index = [f"{donor}:{tissue}:{i}" for i in sub.index.astype(str)]
            if len(sub) > cap:
                rng = np.random.default_rng(SEED + (0 if role == "tumor" else 17))
                pick = rng.choice(sub.index.to_numpy(), size=cap, replace=False)
                sub = sub.loc[pick]
            mat, cells, genes = _to_sparse_cells_x_genes(sub)
            obs = pd.DataFrame(index=cells)
            obs["dataset"] = "GSE123902"
            obs["donor"] = donor
            obs["patient"] = donor
            obs["tissue"] = tissue
            obs["Sample_Origin"] = tissue
            obs["is_normal_tissue"] = tissue == "NORMAL"
            obs["role"] = role
            obs["author_subtype"] = "NA"
            obs["author_lineage"] = "marker_epithelium"
            obs["malig_def"] = "marker_malig"
            obs["unit_id"] = np.where(
                obs["is_normal_tissue"],
                "GSE123902:" + donor + ":NORMAL",
                "GSE123902:" + donor,
            )
            obs["patient_id"] = donor
            obs["histology"] = "LUAD"
            adata_i = ad.AnnData(
                X=mat,
                obs=obs,
                var=pd.DataFrame(index=pd.Index(genes, name="gene")),
            )
            blocks.append(adata_i)
            print(f"    kept cells={adata_i.n_obs} genes={adata_i.n_vars}", flush=True)

    if not blocks:
        raise SystemExit("GSE123902: no epithelial cells extracted")
    import anndata as ad

    shared = blocks[0].var_names
    for b in blocks[1:]:
        shared = shared.intersection(b.var_names)
    blocks = [b[:, shared].copy() for b in blocks]
    out = ad.concat(blocks, axis=0, join="inner", merge="same")
    out.obs_names = "GSE123902:" + out.obs_names.astype(str)
    out.layers["counts"] = out.X.copy()
    print(f"GSE123902 written cells={out.n_obs} genes={out.n_vars}", flush=True)
    return out


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
    adata = ad.AnnData(
        X=sub.T.tocsr(),
        obs=obs,
        var=pd.DataFrame(index=pd.Index([str(g).upper() for g in genes], name="gene")),
    )
    adata.var_names_make_unique()
    adata.obs["dataset"] = "GSE205335"
    adata.obs["histology"] = adata.obs["cancer_subtype"].astype(str)
    adata.obs["donor"] = adata.obs["patient"].astype(str)
    adata.obs["Sample_Origin"] = adata.obs["tissue"].astype(str)
    adata.obs["is_normal_tissue"] = False
    adata.obs["role"] = "tumor"
    adata.obs["author_lineage"] = adata.obs["lineage.total"].astype(str)
    adata.obs["author_subtype"] = adata.obs["lineage.sub"].astype(str)
    adata.obs["malig_def"] = "author_malignant"
    adata.obs["unit_id"] = "GSE205335:" + adata.obs["patient"].astype(str)
    adata.obs["patient_id"] = adata.obs["patient"].astype(str)
    adata.obs_names = "GSE205335:" + adata.obs_names.astype(str)
    adata.layers["counts"] = adata.X.copy()
    print(f"GSE205335 written cells={adata.n_obs} genes={adata.n_vars}", flush=True)
    return adata


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
    p.add_argument("--data", type=Path, default=Path("/tmp/pair_123902_205335_palantir"))
    p.add_argument("--out", type=Path, default=Path("/tmp/pair_123902_205335_palantir/epithelium.h5ad"))
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
        "by_role": joint.obs["role"].astype(str).value_counts().to_dict(),
        "n_units": int(joint.obs["unit_id"].nunique()),
        "cap_per_unit": args.cap_per_unit,
        "gse148071_added": False,
        "dual_high": False,
        "out": str(args.out),
    }
    args.out.with_suffix(".inventory.json").write_text(json.dumps(inv, indent=2))
    print(json.dumps(inv, indent=2), flush=True)


if __name__ == "__main__":
    main()
