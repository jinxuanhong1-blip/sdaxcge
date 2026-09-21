#!/usr/bin/env python3
"""Concordant-4 malignant counts for CytoTRACE2.

Datasets ONLY: GSE123902 + GSE131907 + GSE205335 + GSE189357.
Not GSE148071 / GSE127465 / GSE154826 / GSE200563 / E-MTAB-13526.

Malignant rules (locked):
- GSE131907: author Cell_subtype == Malignant cells (Kim 2020). tLung has none.
- GSE205335: author lineage.sub == Malignant cells on non-normal tissue.
- GSE123902: marker gate (EPCAM|KRT8|KRT18|KRT19)>0 and PTPRC==0,
  PRIMARY_TUMOUR + METASTASIS only. Not CNV. NORMAL dropped.
- GSE189357: same marker gate. Not CNV.

QC: UMI >= 200 and n_genes >= 200. Then cap <= 200 cells / unit.
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
import shutil
import tarfile
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

CAP_PER_UNIT = 200
SEED = 1
MIN_UMI = 200
MIN_GENES = 200
EPI = ("EPCAM", "KRT8", "KRT18", "KRT19")
NAME_123902 = re.compile(
    r"(GSM\d+)_(MSK_LX[^_]+(?:B)?)_(PRIMARY_TUMOUR|METASTASIS|NORMAL)_dense\.csv\.gz"
)
# Locked pseudobulk malignant counts (pre this QC). Audit only.
EXPECTED_GATE = {
    "GSE123902": {
        "LX653": 145, "LX661": 117, "LX666": 786, "LX675": 618, "LX676": 147,
        "LX255B": 228, "LX679": 497, "LX680": 258, "LX681": 1207, "LX682": 421,
        "LX684": 173, "LX699": 46, "LX701": 90,
    },
    "GSE189357": {
        "TD1": 2186, "TD2": 1207, "TD3": 664, "TD4": 491, "TD5": 1736,
        "TD6": 1211, "TD7": 849, "TD8": 1019, "TD9": 4889,
    },
}
STAGE_189357 = {
    "TD1": "IAC", "TD2": "IAC", "TD3": "MIA", "TD4": "MIA", "TD5": "AIS",
    "TD6": "MIA", "TD7": "AIS", "TD8": "AIS", "TD9": "IAC",
}


def collapse_upper(X: sparse.spmatrix, genes: list[str]) -> tuple[sparse.csr_matrix, list[str]]:
    """Sum duplicate symbols after upper-casing. X is cells x genes."""
    genes_u = np.array([str(g).strip().upper() for g in genes], dtype=object)
    keep = np.array([(g != "") and (g != "NAN") and (not str(g).startswith("ENSG")) for g in genes_u])
    # If a row is a symbol that happens to be stored in col 0 as ENSG we already chose symbols.
    X = sparse.csr_matrix(X)[:, keep]
    genes_u = genes_u[keep]
    uniq, inv = np.unique(genes_u, return_inverse=True)
    coo = sparse.coo_matrix(X)
    mat = sparse.coo_matrix(
        (coo.data, (coo.row, inv[coo.col])),
        shape=(X.shape[0], len(uniq)),
    )
    mat.sum_duplicates()
    return mat.tocsr().astype(np.float32), list(uniq)


def gene_index(genes: list[str], name: str) -> int | None:
    try:
        return genes.index(name)
    except ValueError:
        return None


def marker_mask(X: sparse.csr_matrix, genes: list[str]) -> np.ndarray:
    epi_idx = [i for g in EPI if (i := gene_index(genes, g)) is not None]
    if not epi_idx:
        raise SystemExit(f"epithelial markers missing from {genes[:8]}...")
    epi = np.asarray(X[:, epi_idx].sum(axis=1)).ravel() > 0
    pt = gene_index(genes, "PTPRC")
    if pt is None:
        return epi
    ptprc = np.asarray(X[:, pt].todense()).ravel()
    return epi & (ptprc == 0)


def qc_mask(X: sparse.csr_matrix) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n_umi = np.asarray(X.sum(axis=1)).ravel()
    n_genes = np.asarray((X > 0).sum(axis=1)).ravel()
    return (n_umi >= MIN_UMI) & (n_genes >= MIN_GENES), n_umi, n_genes


def cap_rows(obs: pd.DataFrame, unit_col: str, cap: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    keep: list[int] = []
    for _, sub in obs.groupby(unit_col, sort=False):
        if len(sub) > cap:
            chosen = rng.choice(sub.index.to_numpy(), size=cap, replace=False)
            keep.extend(chosen.tolist())
        else:
            keep.extend(sub.index.tolist())
    return obs.loc[keep].copy()


def as_anndata(X: sparse.csr_matrix, obs: pd.DataFrame, genes: list[str]):
    import anndata as ad

    obs = obs.copy()
    for col in obs.columns:
        if pd.api.types.is_bool_dtype(obs[col]):
            obs[col] = obs[col].map({True: "True", False: "False"}).astype(str)
        elif pd.api.types.is_object_dtype(obs[col]):
            obs[col] = obs[col].astype(str)
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(genes, name="gene")))
    adata.layers["counts"] = adata.X.copy()
    return adata


def _donor_from_msk(raw_id: str) -> str:
    return raw_id[4:] if raw_id.startswith("MSK_") else raw_id


def extract_gse123902(data: Path, cap: int) -> tuple:
    import anndata as ad

    tar = data / "GSE123902" / "GSE123902_RAW.tar"
    pieces = []
    audit = []
    with tarfile.open(tar) as tf:
        names = [m.name for m in tf.getmembers() if m.isfile()]
        for name in names:
            m = NAME_123902.match(Path(name).name)
            if not m or m.group(3) == "NORMAL":
                continue
            donor = _donor_from_msk(m.group(2))
            tissue = "PRIMARY" if m.group(3) == "PRIMARY_TUMOUR" else "METASTASIS"
            print(f"  GSE123902 {Path(name).name}", flush=True)
            with gzip.GzipFile(fileobj=tf.extractfile(name)) as handle:
                header = handle.readline().decode().rstrip("\n").split(",")
                genes = [g.strip().strip('"') for g in header[1:]]
                barcodes = []
                rows = []
                for raw in handle:
                    line = raw.decode()
                    comma = line.find(",")
                    barcodes.append(line[:comma].strip())
                    rows.append(np.fromstring(line[comma + 1 :], sep=",", dtype=np.float32))
            X = sparse.csr_matrix(np.vstack(rows))
            del rows
            X, genes_u = collapse_upper(X, genes)
            gate = marker_mask(X, genes_u)
            Xg = X[gate]
            ok, n_umi, n_genes = qc_mask(Xg)
            bc = [f"{m.group(1)}_{barcodes[i]}" for i in np.flatnonzero(gate)]
            obs = pd.DataFrame(
                {
                    "barcode": bc,
                    "dataset": "GSE123902",
                    "unit_id": f"GSE123902:{donor}",
                    "patient_id": donor,
                    "tissue": tissue,
                    "histology": "LUAD",
                    "stage": "NA",
                    "malignant_rule": "marker_(EPCAM|KRT8|KRT18|KRT19)>0_PTPRC==0",
                    "n_umi": n_umi,
                    "n_genes_detected": n_genes,
                }
            )
            audit.append(
                {
                    "dataset": "GSE123902",
                    "patient_id": donor,
                    "file": Path(name).name,
                    "n_raw": int(X.shape[0]),
                    "n_gate": int(gate.sum()),
                    "n_qc": int(ok.sum()),
                    "expected_gate": EXPECTED_GATE["GSE123902"].get(donor),
                }
            )
            print(
                f"    donor={donor} raw={X.shape[0]} gate={int(gate.sum())} qc={int(ok.sum())}",
                flush=True,
            )
            if ok.sum() == 0:
                continue
            obs = obs.loc[ok].copy()
            obs.index = obs["barcode"].astype(str)
            pieces.append(as_anndata(Xg[ok], obs, genes_u))
            del X, Xg
    if not pieces:
        raise SystemExit("GSE123902: no marker-malignant cells")
    adata = ad.concat(pieces, join="outer", fill_value=0, merge="same")
    adata.obs_names_make_unique()
    # outer join may have introduced NA counts
    adata.X = sparse.csr_matrix(adata.X)
    adata.layers["counts"] = adata.X.copy()
    pre = adata.n_obs
    adata = _cap_anndata(adata, cap, SEED)
    catalog = {
        "dataset": "GSE123902",
        "malignant_rule": "marker epithelium, not CNV; NORMAL dropped",
        "n_units_files": len(audit),
        "n_cells_qc_before_cap": int(pre),
        "n_cells": int(adata.n_obs),
        "n_units": int(adata.obs["unit_id"].nunique()),
        "n_genes": int(adata.n_vars),
    }
    return adata, audit, catalog


def _cap_anndata(adata, cap: int, seed: int):
    obs = adata.obs.copy()
    obs["_row"] = np.arange(adata.n_obs)
    kept = cap_rows(obs.set_index("_row", drop=False), "unit_id", cap, seed)
    idx = kept["_row"].to_numpy()
    out = adata[idx].copy()
    out.layers["counts"] = sparse.csr_matrix(out.X)
    return out


def extract_gse189357(data: Path, cap: int) -> tuple:
    import anndata as ad

    tar_path = data / "GSE189357" / "GSE189357_RAW.tar"
    pieces = []
    audit = []
    with tarfile.open(tar_path) as tf:
        members = {Path(m.name).name: m for m in tf.getmembers() if m.isfile()}
        for sample in [f"TD{i}" for i in range(1, 10)]:
            feat = next(n for n in members if n.endswith(f"{sample}_features.tsv.gz"))
            bar = next(n for n in members if n.endswith(f"{sample}_barcodes.tsv.gz"))
            mtx = next(n for n in members if n.endswith(f"{sample}_matrix.mtx.gz"))
            print(f"  GSE189357 {sample}", flush=True)
            with gzip.GzipFile(fileobj=tf.extractfile(members[feat])) as handle:
                genes = []
                for raw in handle:
                    parts = raw.decode().rstrip("\n").split("\t")
                    if len(parts) >= 2 and not parts[1].startswith("ENSG"):
                        genes.append(parts[1])
                    else:
                        genes.append(parts[0])
            with gzip.GzipFile(fileobj=tf.extractfile(members[bar])) as handle:
                barcodes = [raw.decode().rstrip("\n").split("\t")[0] for raw in handle]
            with gzip.GzipFile(fileobj=tf.extractfile(members[mtx])) as handle:
                mat = mmread(handle).tocsr()
            if mat.shape[0] != len(genes):
                if mat.shape[1] == len(genes):
                    mat = mat.T.tocsr()
                else:
                    raise RuntimeError(f"{sample}: mtx {mat.shape} vs genes {len(genes)}")
            X = mat.T.tocsr()  # cells x genes
            del mat
            X, genes_u = collapse_upper(X, genes)
            gate = marker_mask(X, genes_u)
            Xg = X[gate]
            ok, n_umi, n_genes = qc_mask(Xg)
            bc = [f"{sample}_{barcodes[i]}" for i in np.flatnonzero(gate)]
            obs = pd.DataFrame(
                {
                    "barcode": bc,
                    "dataset": "GSE189357",
                    "unit_id": f"GSE189357:{sample}",
                    "patient_id": sample,
                    "tissue": "TUMOR",
                    "histology": "LUAD",
                    "stage": STAGE_189357[sample],
                    "malignant_rule": "marker_(EPCAM|KRT8|KRT18|KRT19)>0_PTPRC==0",
                    "n_umi": n_umi,
                    "n_genes_detected": n_genes,
                }
            )
            audit.append(
                {
                    "dataset": "GSE189357",
                    "patient_id": sample,
                    "file": mtx,
                    "n_raw": int(X.shape[0]),
                    "n_gate": int(gate.sum()),
                    "n_qc": int(ok.sum()),
                    "expected_gate": EXPECTED_GATE["GSE189357"].get(sample),
                }
            )
            print(
                f"    raw={X.shape[0]} gate={int(gate.sum())} qc={int(ok.sum())}",
                flush=True,
            )
            obs = obs.loc[ok].copy()
            obs.index = obs["barcode"].astype(str)
            pieces.append(as_anndata(Xg[ok], obs, genes_u))
            del X, Xg
    adata = ad.concat(pieces, join="outer", fill_value=0, merge="same")
    adata.obs_names_make_unique()
    adata.X = sparse.csr_matrix(adata.X)
    adata.layers["counts"] = adata.X.copy()
    pre = adata.n_obs
    adata = _cap_anndata(adata, cap, SEED + 3)
    catalog = {
        "dataset": "GSE189357",
        "malignant_rule": "marker epithelium, not CNV",
        "n_cells_qc_before_cap": int(pre),
        "n_cells": int(adata.n_obs),
        "n_units": int(adata.obs["unit_id"].nunique()),
        "n_genes": int(adata.n_vars),
    }
    return adata, audit, catalog


def stream_umi_subset(umi_path: Path, keep_ids: list[str]):
    keep_set = set(keep_ids)
    with gzip.open(umi_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        col_idx = [i for i, c in enumerate(cell_ids) if c in keep_set]
        if not col_idx:
            raise SystemExit("GSE131907: no requested cell IDs")
        ordered = [cell_ids[i] for i in col_idx]
        print(f"  UMI header cells={len(cell_ids)} keep={len(col_idx)}", flush=True)
        genes: list[str] = []
        data: list[np.ndarray] = []
        indices: list[np.ndarray] = []
        indptr = [0]
        nnz = 0
        for gi, line in enumerate(handle, start=1):
            raw = line.rstrip("\n")
            if not raw:
                continue
            tab0 = raw.find("\t")
            gene = raw[:tab0]
            vals = np.fromstring(raw[tab0 + 1 :], sep="\t", dtype=np.float32)
            if vals.size != len(cell_ids):
                raise SystemExit(f"GSE131907 row {gi} {gene}: {vals.size} != {len(cell_ids)}")
            sub = vals[col_idx]
            nz = np.flatnonzero(sub)
            if nz.size:
                data.append(sub[nz])
                indices.append(nz.astype(np.int32, copy=False))
                nnz += int(nz.size)
            indptr.append(nnz)
            genes.append(gene)
            if gi % 4000 == 0:
                print(f"    streamed {gi} genes nnz={nnz}", flush=True)
    data_a = np.concatenate(data) if data else np.array([], dtype=np.float32)
    idx_a = np.concatenate(indices) if indices else np.array([], dtype=np.int32)
    mat = sparse.csr_matrix(
        (data_a, idx_a, np.asarray(indptr, dtype=np.int64)),
        shape=(len(genes), len(ordered)),
        dtype=np.float32,
    )
    return ordered, genes, mat


def extract_gse131907(data: Path, cap: int) -> tuple:
    ann = pd.read_csv(
        data / "GSE131907" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz",
        sep="\t",
    )
    mal = ann.loc[ann["Cell_subtype"].astype(str).eq("Malignant cells")].copy()
    tlung_n = int(((ann["Sample_Origin"] == "tLung") & (ann["Cell_subtype"] == "Malignant cells")).sum())
    print(
        f"  GSE131907 author malignant={len(mal)} samples={mal['Sample'].nunique()} tLung_malignant={tlung_n}",
        flush=True,
    )
    mal["Index"] = mal["Index"].astype(str)
    mal = mal.set_index("Index", drop=False)
    mal["unit_id"] = "GSE131907:" + mal["Sample"].astype(str)
    capped = cap_rows(mal, "unit_id", cap, SEED + 1)
    ordered, genes, mat = stream_umi_subset(
        data / "GSE131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
        capped["Index"].astype(str).tolist(),
    )
    X = mat.T.tocsr()
    del mat
    X, genes_u = collapse_upper(X, genes)
    ok, n_umi, n_genes = qc_mask(X)
    obs = capped.set_index("Index").reindex(ordered)
    obs = obs.loc[np.asarray(ordered)[ok]].copy()
    X = X[ok]
    obs["barcode"] = obs.index.astype(str)
    obs["dataset"] = "GSE131907"
    obs["patient_id"] = obs["Sample"].astype(str)
    obs["unit_id"] = "GSE131907:" + obs["patient_id"]
    obs["tissue"] = obs["Sample_Origin"].astype(str)
    obs["histology"] = "LUAD"
    obs["stage"] = "NA"
    obs["malignant_rule"] = "author_Cell_subtype==Malignant cells"
    obs["n_umi"] = n_umi[ok]
    obs["n_genes_detected"] = n_genes[ok]
    keep_cols = [
        "barcode", "dataset", "unit_id", "patient_id", "tissue", "histology",
        "stage", "malignant_rule", "n_umi", "n_genes_detected",
    ]
    obs = obs[keep_cols]
    obs.index = obs["barcode"].astype(str)
    adata = as_anndata(X, obs, genes_u)
    catalog = {
        "dataset": "GSE131907",
        "malignant_rule": "author Cell_subtype==Malignant cells",
        "n_author_malignant": int(len(mal)),
        "n_author_samples": int(mal["Sample"].nunique()),
        "tLung_author_malignant": tlung_n,
        "n_cells_capped_before_qc": int(len(capped)),
        "n_cells": int(adata.n_obs),
        "n_units": int(adata.obs["unit_id"].nunique()),
        "n_genes": int(adata.n_vars),
    }
    by_sample = (
        mal.groupby(mal["Sample"].astype(str)).size().rename("n_author").to_frame()
    )
    by_sample["n_cap"] = capped.groupby(capped["Sample"].astype(str)).size()
    by_sample["n_qc"] = obs.groupby(obs["patient_id"]).size()
    audit = [
        {
            "dataset": "GSE131907",
            "patient_id": pid,
            "file": "raw_UMI",
            "n_raw": int(by_sample.loc[pid, "n_author"]) if pid in by_sample.index else 0,
            "n_gate": int(by_sample.loc[pid, "n_cap"]) if pid in by_sample.index else 0,
            "n_qc": int(by_sample.loc[pid, "n_qc"]) if pid in by_sample.index and pd.notna(by_sample.loc[pid, "n_qc"]) else 0,
            "expected_gate": "",
        }
        for pid in by_sample.index
    ]
    return adata, audit, catalog


def parse_geo_soft(path: Path) -> pd.DataFrame:
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    descriptions: list[str] = []
    titles: list[str] = []
    with gzip.open(path, "rt", errors="replace") as handle:
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


def gunzip_until_rds(src: Path, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    current = src
    tmp_dir = dest.parent
    for i in range(4):
        with current.open("rb") as handle:
            magic = handle.read(2)
        if magic == b"\x1f\x8b":
            nxt = tmp_dir / f"{dest.name}.peel{i}"
            print(f"  gunzip peel {i}", flush=True)
            with gzip.open(current, "rb") as source, nxt.open("wb") as out:
                shutil.copyfileobj(source, out, 16 * 1024 * 1024)
            if current != src and current.exists():
                current.unlink()
            current = nxt
            continue
        break
    if current != dest:
        current.replace(dest)
    return dest


def load_rds_matrix(path: Path):
    import rdata

    with tempfile.TemporaryDirectory(prefix="gse205335-") as tmp:
        rds_path = gunzip_until_rds(path, Path(tmp) / "matrix.rds")
        print(f"  read RDS ({rds_path.stat().st_size} bytes)", flush=True)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message='Missing constructor for R class "dgCMatrix"')
            obj = rdata.read_rds(rds_path)
    expected = {"i", "p", "Dim", "Dimnames", "x"}
    if not expected.issubset(vars(obj)):
        raise TypeError("RDS is not a dgCMatrix")
    genes = np.asarray(obj.Dimnames[0], dtype=str)
    barcodes = np.asarray(obj.Dimnames[1], dtype=str)
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False)
    return matrix, genes, barcodes


def extract_gse205335(data: Path, cap: int) -> tuple:
    ident = pd.read_csv(data / "GSE205335" / "GSE205335_Lung_IO_CellIdentity.txt.gz", sep="\t")
    meta = parse_geo_soft(data / "GSE205335" / "GSE205335_family.soft.gz")
    cells = ident.merge(
        meta[["orig.ident", "gsm", "patient", "tissue", "recist", "cancer_subtype", "tumor_stage"]],
        on="orig.ident",
        how="left",
        validate="many_to_one",
    )
    if cells["patient"].isna().any():
        missing = cells.loc[cells["patient"].isna(), "orig.ident"].unique().tolist()
        raise SystemExit(f"GSE205335 missing patient: {missing}")
    cells["is_normal"] = cells["tissue"].astype(str).str.startswith("Normal ")
    mal = cells.loc[
        cells["lineage.sub"].astype(str).eq("Malignant cells") & ~cells["is_normal"]
    ].copy()
    print(
        f"  GSE205335 author malignant non-normal={len(mal)} patients={mal['patient'].nunique()}",
        flush=True,
    )
    mal["unit_id"] = "GSE205335:" + mal["patient"].astype(str)
    mal = mal.set_index("barcode", drop=False)
    capped = cap_rows(mal, "unit_id", cap, SEED + 2)
    matrix, genes, barcodes = load_rds_matrix(
        data / "GSE205335" / "GSE205335_Lung_IO_UMI_matrix.rds.gz"
    )
    keep = set(capped["barcode"].astype(str))
    col_idx = np.array([i for i, b in enumerate(barcodes) if str(b) in keep], dtype=np.int64)
    print(f"  subset {len(col_idx)} / {len(barcodes)}", flush=True)
    sub = matrix[:, col_idx].tocsr()
    del matrix
    ordered = [str(barcodes[i]) for i in col_idx]
    X = sub.T.tocsr()
    del sub
    X, genes_u = collapse_upper(X, list(genes))
    ok, n_umi, n_genes = qc_mask(X)
    obs = capped.reindex(ordered)
    obs = obs.loc[np.asarray(ordered)[ok]].copy()
    X = X[ok]
    obs["barcode"] = obs.index.astype(str)
    obs["dataset"] = "GSE205335"
    obs["patient_id"] = obs["patient"].astype(str)
    obs["unit_id"] = "GSE205335:" + obs["patient_id"]
    obs["tissue"] = obs["tissue"].astype(str)
    obs["histology"] = obs["cancer_subtype"].astype(str)
    obs["stage"] = obs["tumor_stage"].astype(str)
    obs["malignant_rule"] = "author_lineage.sub==Malignant cells; non-normal"
    obs["n_umi"] = n_umi[ok]
    obs["n_genes_detected"] = n_genes[ok]
    keep_cols = [
        "barcode", "dataset", "unit_id", "patient_id", "tissue", "histology",
        "stage", "malignant_rule", "n_umi", "n_genes_detected",
    ]
    obs = obs[keep_cols]
    obs.index = obs["barcode"].astype(str)
    adata = as_anndata(X, obs, genes_u)
    by = mal.groupby(mal["patient"].astype(str)).size().rename("n_author").to_frame()
    by["n_cap"] = capped.groupby(capped["patient"].astype(str)).size()
    by["n_qc"] = obs.groupby(obs["patient_id"]).size()
    audit = []
    for pid in by.index:
        audit.append(
            {
                "dataset": "GSE205335",
                "patient_id": pid,
                "file": "UMI_matrix.rds",
                "n_raw": int(by.loc[pid, "n_author"]),
                "n_gate": int(by.loc[pid, "n_cap"]),
                "n_qc": int(by.loc[pid, "n_qc"]) if pd.notna(by.loc[pid, "n_qc"]) else 0,
                "expected_gate": "",
            }
        )
    catalog = {
        "dataset": "GSE205335",
        "malignant_rule": "author lineage.sub==Malignant cells; drop Normal *",
        "n_author_malignant_non_normal": int(len(mal)),
        "n_patients": int(mal["patient"].nunique()),
        "n_cells_capped_before_qc": int(len(capped)),
        "n_cells": int(adata.n_obs),
        "n_units": int(adata.obs["unit_id"].nunique()),
        "n_genes": int(adata.n_vars),
        "histology": adata.obs["histology"].value_counts().to_dict(),
    }
    return adata, audit, catalog


EXTRACTORS = {
    "GSE123902": extract_gse123902,
    "GSE189357": extract_gse189357,
    "GSE131907": extract_gse131907,
    "GSE205335": extract_gse205335,
}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=Path, default=Path("/tmp/concordant4_ct2"))
    p.add_argument("--out", type=Path, default=Path("/tmp/concordant4_ct2/h5ad"))
    p.add_argument("--cap-per-unit", type=int, default=CAP_PER_UNIT)
    p.add_argument("--dataset", action="append", default=None)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    which = args.dataset or list(EXTRACTORS)
    audits = []
    catalogs = {}
    for name in which:
        print(f"EXTRACT {name}", flush=True)
        adata, audit, catalog = EXTRACTORS[name](args.data, args.cap_per_unit)
        dest = args.out / f"{name}.h5ad"
        adata.write_h5ad(dest)
        print(f"wrote {dest} cells={adata.n_obs} genes={adata.n_vars} units={catalog['n_units']}", flush=True)
        audits.extend(audit)
        catalogs[name] = catalog
        del adata
    audit_path = args.out / "extract_audit.tsv"
    prev = pd.read_csv(audit_path, sep="\t") if audit_path.exists() else pd.DataFrame()
    new = pd.DataFrame(audits)
    if len(prev):
        prev = prev.loc[~prev["dataset"].isin(new["dataset"])]
        new = pd.concat([prev, new], ignore_index=True)
    new.to_csv(audit_path, sep="\t", index=False)
    cat_path = args.out / "catalog.json"
    old = json.loads(cat_path.read_text()) if cat_path.exists() else {}
    old.update(catalogs)
    cat_path.write_text(json.dumps(old, indent=2))
    print(json.dumps(catalogs, indent=2), flush=True)


if __name__ == "__main__":
    main()
