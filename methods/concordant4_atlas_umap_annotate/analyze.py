#!/usr/bin/env python3
"""Concordant-4 cell-level atlas: load, Harmony, Leiden, annotate.

Datasets ONLY: GSE123902 + GSE131907 + GSE205335 + GSE189357.
Public processed matrices. CLDN4-only. Methods/identity figure.

Does NOT re-run or re-audit PR #503 T/NK ρ (n=65) or IFN/MHC DE.
"""
from __future__ import annotations

import argparse
import gc
import gzip
import json
import re
import shutil
import tarfile
import tempfile
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

import anndata as ad
import scanpy as sc

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
TABLES = RESULTS / "tables"
FIGS = RESULTS / "figures"
DEFAULT_DATA = Path("/tmp/geo_atlas")

COHORTS = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
TUMOR_ORIGINS = {"tLung", "tL/B", "mLN", "PE", "mBrain"}
CAP_PER_UNIT = 350
MIN_GENES = 200
MIN_COUNTS = 500
MAX_MITO = 20.0
N_HVG = 2000
N_PCS = 30
N_NEIGHBORS = 15
LEIDEN_RES = 0.6
HARMONY_THETA = 2.0
PLOT_MAX = 80_000
SEED = 1

NAME_123902 = re.compile(
    r"(GSM\d+)_(MSK_LX[^_]+(?:B)?)_(PRIMARY_TUMOUR|METASTASIS|NORMAL)_dense\.csv\.gz"
)

MARKERS = {
    "malignant": ["EPCAM", "KRT8", "KRT18", "KRT19", "CLDN4"],
    "T": ["CD3D", "CD3E"],
    "NK": ["NKG7", "GNLY", "KLRD1"],
    "myeloid": ["CD68", "LYZ"],
    "B": ["MS4A1", "CD79A"],
    "endothelial": ["PECAM1", "VWF"],
    "fibroblast": ["COL1A1"],
}
DOT_GENES = [
    "EPCAM", "KRT8", "KRT18", "KRT19", "CLDN4",
    "CD3D", "CD3E", "NKG7", "GNLY", "KLRD1",
    "CD68", "LYZ", "MS4A1", "CD79A",
    "PECAM1", "VWF", "COL1A1",
]
COARSE_ORDER = ["malignant", "T", "NK", "myeloid", "B", "other"]
DATASET_COLORS = {
    "GSE123902": "#4c78a8",
    "GSE131907": "#f58518",
    "GSE205335": "#54a24b",
    "GSE189357": "#b279a2",
}
ANNO_COLORS = {
    "malignant": "#c0392b",
    "T": "#2471a3",
    "NK": "#16a085",
    "myeloid": "#e67e22",
    "B": "#8e44ad",
    "other": "#7f8c8d",
}

MM = 1.0 / 25.4


def nature_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans", "Helvetica"],
            "font.size": 7,
            "axes.titlesize": 7,
            "axes.labelsize": 7,
            "xtick.labelsize": 6,
            "ytick.labelsize": 6,
            "legend.fontsize": 6,
            "axes.linewidth": 0.5,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "xtick.major.size": 2.0,
            "ytick.major.size": 2.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.02,
        }
    )


def despine(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def panel_label(ax, letter: str, x: float = -0.12, y: float = 1.08) -> None:
    ax.text(
        x,
        y,
        letter,
        transform=ax.transAxes,
        fontsize=8,
        fontweight="bold",
        va="top",
        ha="left",
    )


def mito_mask(genes: list[str]) -> np.ndarray:
    return np.array(
        [g.upper().startswith("MT-") or g.upper().startswith("MT.") for g in genes],
        dtype=bool,
    )


def qc_from_counts(X: sparse.spmatrix, genes: list[str]) -> pd.DataFrame:
    if not sparse.isspmatrix_csr(X):
        X = X.tocsr()
    n_counts = np.asarray(X.sum(axis=1)).ravel()
    n_genes = np.asarray((X > 0).sum(axis=1)).ravel()
    mt = mito_mask(genes)
    if mt.any():
        mito = np.asarray(X[:, mt].sum(axis=1)).ravel()
        pct_mito = np.where(n_counts > 0, 100.0 * mito / n_counts, 0.0)
    else:
        pct_mito = np.zeros(X.shape[0], dtype=float)
    return pd.DataFrame({"n_counts": n_counts, "n_genes": n_genes, "pct_mito": pct_mito})


def pass_qc(qc: pd.DataFrame) -> np.ndarray:
    return (
        (qc["n_genes"].to_numpy() >= MIN_GENES)
        & (qc["n_counts"].to_numpy() >= MIN_COUNTS)
        & (qc["pct_mito"].to_numpy() < MAX_MITO)
    )


def cap_index(labels: np.ndarray, cap: int, rng: np.random.Generator) -> np.ndarray:
    keep: list[int] = []
    for lab in pd.unique(labels):
        idx = np.flatnonzero(labels == lab)
        if len(idx) > cap:
            idx = rng.choice(idx, size=cap, replace=False)
        keep.extend(idx.tolist())
    return np.asarray(sorted(keep), dtype=int)


def reservoir_push(store: list, n_seen: int, item, cap: int, rng: np.random.Generator) -> None:
    if len(store) < cap:
        store.append(item)
        return
    j = int(rng.integers(0, n_seen))
    if j < cap:
        store[j] = item


def genes_to_upper_unique(genes: list[str]) -> tuple[list[str], np.ndarray]:
    """Return unique uppercase gene names and a keep-index into the original list."""
    seen: dict[str, int] = {}
    keep: list[int] = []
    names: list[str] = []
    for i, g in enumerate(genes):
        u = str(g).strip().upper()
        if not u or u in seen:
            continue
        seen[u] = i
        keep.append(i)
        names.append(u)
    return names, np.asarray(keep, dtype=int)


def csr_from_dense_rows(rows: list[np.ndarray]) -> sparse.csr_matrix:
    if not rows:
        return sparse.csr_matrix((0, 0), dtype=np.float32)
    return sparse.csr_matrix(np.vstack(rows), dtype=np.float32)


# ---------------------------------------------------------------------------
# GSE123902 — Laughney 2020 dense UMI CSVs (tumor/met donors)
# ---------------------------------------------------------------------------
def load_gse123902(data: Path, cap: int, rng: np.random.Generator) -> tuple[ad.AnnData, dict]:
    tar = data / "GSE123902" / "GSE123902_RAW.tar"
    if not tar.is_file():
        raise SystemExit(f"missing {tar}")
    csv_dir = data / "GSE123902" / "csv"
    csv_dir.mkdir(parents=True, exist_ok=True)
    files = sorted(csv_dir.glob("GSM*_dense.csv.gz"))
    if not files:
        print(f"extract {tar}", flush=True)
        with tarfile.open(tar) as tf:
            tf.extractall(csv_dir)
        files = sorted(csv_dir.glob("GSM*_dense.csv.gz"))
    reservoirs: dict[str, list] = {}
    n_raw: dict[str, int] = {}
    n_qc: dict[str, int] = {}
    n_files = 0
    for path in files:
        m = NAME_123902.match(path.name)
        if not m:
            continue
        site = m.group(3)
        if site == "NORMAL":
            continue
        raw_id = m.group(2)
        donor = raw_id[4:] if raw_id.startswith("MSK_") else raw_id
        tissue = "PRIMARY" if site == "PRIMARY_TUMOUR" else "METASTASIS"
        n_files += 1
        reservoirs.setdefault(donor, [])
        n_raw.setdefault(donor, 0)
        n_qc.setdefault(donor, 0)
        with gzip.open(path, "rt") as handle:
            header = handle.readline().rstrip("\n").split(",")
            genes_u, keep_g = genes_to_upper_unique([g.strip() for g in header[1:]])
            mt = mito_mask(genes_u)
            for line in handle:
                parts = line.rstrip("\n").split(",")
                raw_vals = np.fromstring(",".join(parts[1:]), sep=",", dtype=np.float32)
                if raw_vals.size != len(header) - 1:
                    continue
                vals = raw_vals[keep_g]
                n_raw[donor] += 1
                n_counts = float(vals.sum())
                n_genes = int((vals > 0).sum())
                mito = float(vals[mt].sum()) if mt.any() else 0.0
                pct = 100.0 * mito / n_counts if n_counts > 0 else 0.0
                if n_genes < MIN_GENES or n_counts < MIN_COUNTS or pct >= MAX_MITO:
                    continue
                n_qc[donor] += 1
                item = (vals, genes_u, parts[0], m.group(1), tissue, path.name)
                reservoir_push(reservoirs[donor], n_qc[donor], item, cap, rng)
        print(
            f"  GSE123902 {path.name} donor={donor} raw_so_far={n_raw[donor]} "
            f"qc={n_qc[donor]} kept={len(reservoirs[donor])}",
            flush=True,
        )
    if not any(reservoirs.values()):
        raise SystemExit("GSE123902: no tumor CSVs")
    gene_sets = [set(item[1]) for items in reservoirs.values() for item in items]
    genes = sorted(set.intersection(*gene_sets)) if gene_sets else []
    print(f"  GSE123902 shared genes={len(genes)}", flush=True)
    gene_index = {g: i for i, g in enumerate(genes)}
    rows = []
    obs_rows = []
    for donor, items in reservoirs.items():
        for vals, genes_u, barcode, gsm, tissue, fname in items:
            aligned = np.zeros(len(genes), dtype=np.float32)
            lookup = np.fromiter((gene_index.get(g, -1) for g in genes_u), dtype=np.int32, count=len(genes_u))
            ok = lookup >= 0
            aligned[lookup[ok]] = vals[ok]
            rows.append(aligned)
            obs_rows.append(
                {
                    "barcode": f"{gsm}_{barcode}",
                    "dataset": "GSE123902",
                    "unit_id": donor,
                    "sample_id": gsm,
                    "patient": donor,
                    "tissue": tissue,
                    "author_label": "",
                    "author_malignant": False,
                }
            )
    X = csr_from_dense_rows(rows)
    obs = pd.DataFrame(obs_rows)
    obs.index = obs["barcode"].astype(str)
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))
    adata.var_names_make_unique()
    info = {
        "dataset": "GSE123902",
        "unit": "donor",
        "n_files_tumor": n_files,
        "n_units": int(obs["unit_id"].nunique()) if len(obs) else 0,
        "n_cells_raw": int(sum(n_raw.values())),
        "n_cells_qc": int(sum(n_qc.values())),
        "n_cells_used": int(adata.n_obs),
        "n_units_raw": int(len(n_raw)),
        "capped": True,
        "cap_per_unit": cap,
        "note": "Laughney 2020; PRIMARY+METASTASIS; NORMAL dropped; no author cell types",
        "per_unit_qc": n_qc,
        "per_unit_used": obs["unit_id"].value_counts().to_dict() if len(obs) else {},
    }
    print(f"GSE123902 used={adata.n_obs} units={info['n_units']} qc={info['n_cells_qc']}", flush=True)
    return adata, info


# ---------------------------------------------------------------------------
# GSE131907 — Kim 2020 UMI + author annotation (21 tumor-bearing units)
# ---------------------------------------------------------------------------
def stream_umi_subset(umi_path: Path, keep_ids: list[str]) -> tuple[list[str], list[str], sparse.csr_matrix]:
    keep_set = set(keep_ids)
    with gzip.open(umi_path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        col_idx = [i for i, c in enumerate(cell_ids) if c in keep_set]
        if not col_idx:
            raise SystemExit("GSE131907: no requested cell IDs in UMI header")
        ordered = [cell_ids[i] for i in col_idx]
        print(f"  UMI header cells={len(cell_ids)} keep={len(col_idx)}", flush=True)
        genes: list[str] = []
        data: list[np.ndarray] = []
        indices: list[np.ndarray] = []
        indptr = [0]
        nnz = 0
        for gi, line in enumerate(f, start=1):
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
                data.append(sub[nz].astype(np.float32, copy=False))
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


def load_gse131907(data: Path, cap: int, rng: np.random.Generator) -> tuple[ad.AnnData, dict]:
    ann_path = data / "GSE131907" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    umi_path = data / "GSE131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    ann = pd.read_csv(ann_path, sep="\t")
    tumor = ann[ann["Sample_Origin"].isin(TUMOR_ORIGINS)].copy()
    n_mal = tumor.groupby("Sample")["Cell_subtype"].apply(lambda s: int(s.eq("Malignant cells").sum()))
    keep_samples = set(n_mal[n_mal >= 20].index.astype(str))
    cells = tumor[tumor["Sample"].astype(str).isin(keep_samples)].copy()
    n_author = int(len(cells))
    n_units = int(cells["Sample"].nunique())
    # cap from author-QC cells, then stream
    cells = cells.set_index("Index", drop=False)
    keep_idx = cap_index(cells["Sample"].astype(str).to_numpy(), cap, rng)
    cells = cells.iloc[keep_idx].copy()
    ordered, genes_raw, mat = stream_umi_subset(umi_path, cells["Index"].astype(str).tolist())
    genes, keep_g = genes_to_upper_unique(genes_raw)
    X = mat[keep_g, :].T.tocsr()
    obs = cells.set_index("Index").reindex(ordered)
    qc = qc_from_counts(X, genes)
    ok = pass_qc(qc)
    X = X[ok]
    obs = obs.loc[np.asarray(ordered)[ok]].copy()
    obs["n_counts"] = qc.loc[ok, "n_counts"].to_numpy()
    obs["n_genes"] = qc.loc[ok, "n_genes"].to_numpy()
    obs["pct_mito"] = qc.loc[ok, "pct_mito"].to_numpy()
    obs["dataset"] = "GSE131907"
    obs["unit_id"] = obs["Sample"].astype(str)
    obs["sample_id"] = obs["Sample"].astype(str)
    obs["patient"] = obs["Sample"].astype(str)
    obs["tissue"] = obs["Sample_Origin"].astype(str)
    obs["author_label"] = obs["Cell_type"].astype(str)
    obs["author_malignant"] = obs["Cell_subtype"].astype(str).eq("Malignant cells")
    obs["barcode"] = obs.index.astype(str)
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))
    adata.var_names_make_unique()
    info = {
        "dataset": "GSE131907",
        "unit": "sample",
        "n_units": int(adata.obs["unit_id"].nunique()),
        "n_units_rule": n_units,
        "n_cells_author_in_units": n_author,
        "n_cells_capped": int(ok.size),
        "n_cells_qc": int(ok.sum()),
        "n_cells_used": int(adata.n_obs),
        "capped": True,
        "cap_per_unit": cap,
        "keep_samples": sorted(keep_samples),
        "note": "Kim 2020; tumor-bearing Sample_Origin; author Malignant cells n>=20 (21 samples); all compartments",
    }
    print(f"GSE131907 used={adata.n_obs} units={info['n_units']} author_in_units={n_author}", flush=True)
    return adata, info


# ---------------------------------------------------------------------------
# GSE205335 — Ahn/Lee ICI; author identity + dgCMatrix RDS
# ---------------------------------------------------------------------------
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


def gunzip_until_rds(src: Path, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        with dest.open("rb") as handle:
            if handle.read(2) == b"X\n":
                return dest
    current = src
    tmp_dir = dest.parent
    for i in range(4):
        with current.open("rb") as handle:
            magic = handle.read(2)
        if magic == b"\x1f\x8b":
            nxt = tmp_dir / f"{dest.name}.peel{i}"
            print(f"  gunzip peel {i}: {current}", flush=True)
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
        print(f"  read RDS {rds_path} ({rds_path.stat().st_size} bytes)", flush=True)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message='Missing constructor for R class "dgCMatrix"')
            obj = rdata.read_rds(rds_path)
    expected = {"i", "p", "Dim", "Dimnames", "x"}
    if not expected.issubset(vars(obj)):
        raise TypeError("RDS is not the expected Matrix::dgCMatrix object")
    genes = np.asarray(obj.Dimnames[0], dtype=str)
    barcodes = np.asarray(obj.Dimnames[1], dtype=str)
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False)
    return matrix, genes, barcodes


def load_gse205335(data: Path, cap: int, rng: np.random.Generator) -> tuple[ad.AnnData, dict]:
    ident = pd.read_csv(data / "GSE205335" / "GSE205335_Lung_IO_CellIdentity.txt.gz", sep="\t")
    meta = parse_geo_soft(data / "GSE205335" / "GSE205335_family.soft.gz")
    cells = ident.merge(
        meta[["orig.ident", "gsm", "patient", "tissue", "recist", "cancer_subtype"]],
        on="orig.ident",
        how="left",
        validate="many_to_one",
    )
    if cells["patient"].isna().any():
        missing = cells.loc[cells["patient"].isna(), "orig.ident"].unique().tolist()
        raise SystemExit(f"GSE205335 identity missing SOFT patient: {missing}")
    cells["is_normal_tissue"] = cells["tissue"].astype(str).str.startswith("Normal ")
    mal_patients = set(
        cells.loc[
            cells["lineage.sub"].eq("Malignant cells") & ~cells["is_normal_tissue"],
            "patient",
        ].astype(str)
    )
    use = cells.loc[~cells["is_normal_tissue"] & cells["patient"].astype(str).isin(mal_patients)].copy()
    n_author = int(len(use))
    n_units = int(use["patient"].nunique())
    use = use.set_index("barcode", drop=False)
    keep_idx = cap_index(use["patient"].astype(str).to_numpy(), cap, rng)
    use = use.iloc[keep_idx].copy()
    keep = set(use["barcode"].astype(str))
    matrix, genes_raw, barcodes = load_rds_matrix(data / "GSE205335" / "GSE205335_Lung_IO_UMI_matrix.rds.gz")
    col_idx = np.array([i for i, b in enumerate(barcodes) if b in keep], dtype=np.int64)
    print(f"  GSE205335 subset {len(col_idx)} / {len(barcodes)}", flush=True)
    sub = matrix[:, col_idx].tocsr()
    del matrix
    gc.collect()
    ordered = [str(barcodes[i]) for i in col_idx]
    genes, keep_g = genes_to_upper_unique(list(genes_raw))
    X = sub[keep_g, :].T.tocsr()
    obs = use.reindex(ordered)
    qc = qc_from_counts(X, genes)
    ok = pass_qc(qc)
    X = X[ok]
    obs = obs.loc[np.asarray(ordered)[ok]].copy()
    obs["dataset"] = "GSE205335"
    obs["unit_id"] = obs["patient"].astype(str)
    obs["sample_id"] = obs["orig.ident"].astype(str)
    obs["author_label"] = obs["lineage.sub"].astype(str)
    obs["author_malignant"] = obs["lineage.sub"].astype(str).eq("Malignant cells")
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))
    adata.var_names_make_unique()
    info = {
        "dataset": "GSE205335",
        "unit": "patient",
        "n_units": int(adata.obs["unit_id"].nunique()),
        "n_units_rule": n_units,
        "n_cells_author_in_units": n_author,
        "n_cells_capped": int(ok.size),
        "n_cells_qc": int(ok.sum()),
        "n_cells_used": int(adata.n_obs),
        "capped": True,
        "cap_per_unit": cap,
        "note": "Ahn/Lee eLife 2024; author lineage.sub; non-normal tissues; 22 patients with malignant; all compartments",
    }
    print(f"GSE205335 used={adata.n_obs} units={info['n_units']} author_in_units={n_author}", flush=True)
    return adata, info


# ---------------------------------------------------------------------------
# GSE189357 — Zhu 2022 10x MTX (9 patients, all compartments)
# ---------------------------------------------------------------------------
def _tar_members(tf: tarfile.TarFile) -> dict[str, tarfile.TarInfo]:
    return {m.name: m for m in tf.getmembers() if m.isfile()}


def _find_member(members: dict[str, tarfile.TarInfo], sample: str, tokens: tuple[str, ...]) -> str:
    hits = [
        n
        for n in members
        if (f"_{sample}_" in Path(n).name or f"_{sample}." in Path(n).name)
        and any(t in n.lower() for t in tokens)
    ]
    if not hits:
        raise FileNotFoundError(f"{sample}: no file matching {tokens}")
    hits.sort(key=len)
    return hits[0]


def _read_tsv_gz(tf: tarfile.TarFile, member: tarfile.TarInfo) -> list[list[str]]:
    rows = []
    with gzip.GzipFile(fileobj=tf.extractfile(member)) as handle:
        for raw in handle:
            rows.append(raw.decode().rstrip("\n").split("\t"))
    return rows


def load_gse189357(data: Path, cap: int, rng: np.random.Generator) -> tuple[ad.AnnData, dict]:
    tar_path = data / "GSE189357" / "GSE189357_RAW.tar"
    pieces = []
    per_qc = {}
    per_raw = {}
    with tarfile.open(tar_path) as tf:
        members = _tar_members(tf)
        print(f"  GSE189357 tar members={len(members)}", flush=True)
        for sample in [f"TD{i}" for i in range(1, 10)]:
            feat_n = _find_member(members, sample, ("features", "genes"))
            bar_n = _find_member(members, sample, ("barcodes",))
            mtx_n = _find_member(members, sample, ("matrix.mtx",))
            feat_rows = _read_tsv_gz(tf, members[feat_n])
            genes_raw = []
            for p in feat_rows:
                if len(p) >= 2 and not p[1].startswith("ENSG"):
                    genes_raw.append(p[1])
                else:
                    genes_raw.append(p[0])
            barcodes = [r[0] for r in _read_tsv_gz(tf, members[bar_n])]
            with gzip.GzipFile(fileobj=tf.extractfile(members[mtx_n])) as handle:
                mat = mmread(handle).tocsc()
            if mat.shape[0] != len(genes_raw):
                if mat.shape[1] == len(genes_raw):
                    mat = mat.T.tocsc()
                else:
                    raise RuntimeError(f"{sample}: mtx {mat.shape} vs genes {len(genes_raw)}")
            genes, keep_g = genes_to_upper_unique(genes_raw)
            X = mat[keep_g, :].T.tocsr()
            per_raw[sample] = int(X.shape[0])
            qc = qc_from_counts(X, genes)
            ok = pass_qc(qc)
            per_qc[sample] = int(ok.sum())
            X = X[ok]
            bc = np.asarray(barcodes)[ok]
            if X.shape[0] > cap:
                pick = rng.choice(X.shape[0], size=cap, replace=False)
                pick.sort()
                X = X[pick]
                bc = bc[pick]
            obs = pd.DataFrame(index=[f"{sample}_{b}" for b in bc])
            obs["barcode"] = obs.index
            obs["dataset"] = "GSE189357"
            obs["unit_id"] = sample
            obs["sample_id"] = sample
            obs["patient"] = sample
            obs["tissue"] = "TUMOR"
            obs["author_label"] = ""
            obs["author_malignant"] = False
            adata_s = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))
            adata_s.var_names_make_unique()
            pieces.append(adata_s)
            print(f"  {sample}: raw={per_raw[sample]} qc={per_qc[sample]} used={adata_s.n_obs}", flush=True)
            del mat
    # inner-join genes across patients
    common = set(pieces[0].var_names)
    for a in pieces[1:]:
        common &= set(a.var_names)
    common = sorted(common)
    aligned = [a[:, common].copy() for a in pieces]
    adata = ad.concat(aligned, axis=0, merge="same")
    adata.obs_names_make_unique()
    info = {
        "dataset": "GSE189357",
        "unit": "patient",
        "n_units": int(adata.obs["unit_id"].nunique()),
        "n_cells_raw": int(sum(per_raw.values())),
        "n_cells_qc": int(sum(per_qc.values())),
        "n_cells_used": int(adata.n_obs),
        "capped": True,
        "cap_per_unit": cap,
        "per_unit_qc": per_qc,
        "note": "Zhu/Wang 2022 AIS–IAC; 10x MTX; no author cell types; all compartments",
    }
    print(f"GSE189357 used={adata.n_obs} units={info['n_units']} qc={info['n_cells_qc']}", flush=True)
    return adata, info


def align_and_concat(adatas: list[ad.AnnData]) -> ad.AnnData:
    common = set(adatas[0].var_names)
    for a in adatas[1:]:
        common &= set(a.var_names)
    common = sorted(g for g in common if g and g != "NAN")
    print(f"inner-join genes={len(common)}", flush=True)
    pieces = []
    for a in adatas:
        b = a[:, common].copy()
        b.obs_names = [f"{r.dataset}_{r.barcode}" for r in b.obs.itertuples()]
        b.obs_names_make_unique()
        pieces.append(b)
    out = ad.concat(pieces, axis=0, merge="same")
    out.obs_names_make_unique()
    out.layers["counts"] = out.X.copy()
    return out


def integrate(adata: ad.AnnData) -> ad.AnnData:
    adata.obs["dataset"] = adata.obs["dataset"].astype("category")
    adata.obs["unit_id"] = adata.obs["unit_id"].astype(str)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    try:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=N_HVG,
            batch_key="dataset",
            flavor="seurat_v3",
            layer="counts",
        )
    except Exception as exc:
        print(f"seurat_v3 HVG failed ({exc}); falling back to seurat", flush=True)
        sc.pp.highly_variable_genes(
            adata, n_top_genes=N_HVG, batch_key="dataset", flavor="seurat"
        )
    adata.raw = adata
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata, n_comps=N_PCS, svd_solver="arpack")
    print("Harmony batch=dataset", flush=True)
    import harmonypy as hm

    ho = hm.run_harmony(
        np.asarray(adata.obsm["X_pca"]),
        adata.obs,
        ["dataset"],
        theta=HARMONY_THETA,
        max_iter_harmony=20,
        random_state=SEED,
    )
    z = np.asarray(ho.Z_corr)
    # harmonypy ≥1.1 returns cells × PCs; older builds returned PCs × cells.
    if z.ndim != 2:
        raise RuntimeError(f"Harmony Z_corr has unexpected shape {z.shape}")
    if z.shape[0] == adata.n_obs:
        adata.obsm["X_pca_harmony"] = z
    elif z.shape[1] == adata.n_obs:
        adata.obsm["X_pca_harmony"] = z.T
    else:
        raise RuntimeError(f"Harmony Z_corr {z.shape} does not match n_obs={adata.n_obs}")
    sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=N_NEIGHBORS, n_pcs=N_PCS)
    sc.tl.umap(adata, min_dist=0.4, spread=1.0, random_state=SEED)
    try:
        sc.tl.leiden(adata, resolution=LEIDEN_RES, flavor="igraph", n_iterations=2, random_state=SEED)
    except TypeError:
        sc.tl.leiden(adata, resolution=LEIDEN_RES, random_state=SEED)
    return adata


def _raw_vector(adata: ad.AnnData, gene: str) -> np.ndarray:
    if adata.raw is None or gene not in adata.raw.var_names:
        return np.zeros(adata.n_obs, dtype=float)
    x = adata.raw[:, gene].X
    if sparse.issparse(x):
        x = x.toarray()
    return np.asarray(x).ravel()


def marker_score(adata: ad.AnnData, genes: list[str]) -> np.ndarray:
    present = [g for g in genes if adata.raw is not None and g in adata.raw.var_names]
    if not present:
        return np.zeros(adata.n_obs, dtype=float)
    return np.mean(np.vstack([_raw_vector(adata, g) for g in present]), axis=0)


def annotate(adata: ad.AnnData) -> pd.DataFrame:
    scores = {name: marker_score(adata, genes) for name, genes in MARKERS.items()}
    for name, s in scores.items():
        adata.obs[f"score_{name}"] = s
    adata.obs["cldn4"] = _raw_vector(adata, "CLDN4")
    adata.obs["author_malignant"] = adata.obs["author_malignant"].astype(bool)

    rows = []
    labels = pd.Series(index=adata.obs_names, dtype=object)
    fine = pd.Series(index=adata.obs_names, dtype=object)
    for cl in sorted(adata.obs["leiden"].astype(str).unique(), key=lambda x: int(x) if x.isdigit() else x):
        mask = adata.obs["leiden"].astype(str).eq(cl)
        n = int(mask.sum())
        author_frac = float(adata.obs.loc[mask, "author_malignant"].mean())
        means = {k: float(v[mask.to_numpy()].mean()) for k, v in scores.items()}
        ranked = sorted(means.items(), key=lambda kv: kv[1], reverse=True)
        top_name, top_val = ranked[0]
        second_name, second_val = ranked[1]
        rule = []
        if author_frac >= 0.50:
            lab = "malignant"
            fine_lab = "malignant (author)"
            rule.append(f"author_malignant_frac={author_frac:.2f}≥0.50")
        elif top_val < 0.15:
            lab = "other"
            fine_lab = "other (low marker scores)"
            rule.append(f"top {top_name}={top_val:.3f}<0.15")
        elif (top_val - second_val) < 0.05 and second_val >= 0.15:
            lab = "other"
            fine_lab = f"mixed {top_name}/{second_name}"
            rule.append(
                f"mixed: {top_name}={top_val:.3f} vs {second_name}={second_val:.3f} (Δ<0.05)"
            )
        elif top_name == "malignant":
            lab = "malignant"
            fine_lab = "malignant (marker epithelium)"
            rule.append(f"EPCAM/KRT/CLDN4 score={top_val:.3f}")
        elif top_name in {"T", "NK", "myeloid", "B"}:
            lab = top_name
            fine_lab = top_name
            rule.append(f"{top_name} score={top_val:.3f}")
        elif top_name == "endothelial":
            lab = "other"
            fine_lab = "endothelial"
            rule.append(f"PECAM1/VWF score={top_val:.3f}")
        elif top_name == "fibroblast":
            lab = "other"
            fine_lab = "fibroblast"
            rule.append(f"COL1A1 score={top_val:.3f}")
        else:
            lab = "other"
            fine_lab = "other"
            rule.append("fallback")
        labels.loc[mask] = lab
        fine.loc[mask] = fine_lab
        # top genes from rank_genes_groups filled later
        rows.append(
            {
                "cluster": cl,
                "n_cells": n,
                "frac_cells": n / adata.n_obs,
                "label": lab,
                "label_fine": fine_lab,
                "rule": "; ".join(rule),
                "author_malignant_frac": author_frac,
                "n_datasets": int(adata.obs.loc[mask, "dataset"].nunique()),
                "datasets": ",".join(sorted(adata.obs.loc[mask, "dataset"].astype(str).unique())),
                **{f"score_{k}": v for k, v in means.items()},
                "top_markers": "",
            }
        )
    adata.obs["annotation"] = pd.Categorical(labels, categories=COARSE_ORDER)
    adata.obs["annotation_fine"] = fine.astype(str)
    table = pd.DataFrame(rows).sort_values("cluster")
    print("rank_genes_groups wilcoxon", flush=True)
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", n_genes=8, use_raw=True)
    names = adata.uns["rank_genes_groups"]["names"]
    top = []
    for cl in table["cluster"]:
        genes = [str(names[cl][i]) for i in range(min(5, len(names[cl])))]
        top.append(",".join(genes))
    table["top_markers"] = top
    return table


def subsample_plot(adata: ad.AnnData, rng: np.random.Generator) -> np.ndarray:
    n = adata.n_obs
    if n <= PLOT_MAX:
        return np.arange(n)
    return rng.choice(n, size=PLOT_MAX, replace=False)


def _umap_xy(adata: ad.AnnData, idx: np.ndarray) -> np.ndarray:
    return np.asarray(adata.obsm["X_umap"][idx])


def scatter_discrete(ax, xy, labels, colors, order, s=1.4, alpha=0.85) -> None:
    for lab in order:
        m = labels == lab
        if not m.any():
            continue
        ax.scatter(
            xy[m, 0],
            xy[m, 1],
            s=s,
            c=colors.get(lab, "#444444"),
            linewidths=0,
            alpha=alpha,
            rasterized=True,
            label=lab,
        )
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    despine(ax)
    ax.set_aspect("equal", adjustable="datalim")


def legend_compact(ax, handles=None, **kwargs) -> None:
    kw = dict(frameon=False, handletextpad=0.3, borderpad=0.1, labelspacing=0.15, markerscale=4.0)
    kw.update(kwargs)
    if handles is None:
        ax.legend(**kw)
    else:
        ax.legend(handles=handles, **kw)


def make_figure(adata: ad.AnnData, inventory: pd.DataFrame, ann: pd.DataFrame, idx: np.ndarray, path: Path) -> None:
    nature_style()
    xy = _umap_xy(adata, idx)
    obs = adata.obs.iloc[idx]
    fig = plt.figure(figsize=(180 * MM, 152 * MM), facecolor="white")
    gs = fig.add_gridspec(
        2,
        3,
        left=0.06,
        right=0.99,
        top=0.93,
        bottom=0.08,
        wspace=0.32,
        hspace=0.42,
        height_ratios=[1.0, 1.08],
    )

    # a — overview
    ax = fig.add_subplot(gs[0, 0])
    panel_label(ax, "a", x=-0.18, y=1.10)
    inv = inventory.set_index("dataset").loc[COHORTS]
    x = np.arange(len(COHORTS))
    w = 0.38
    ax.bar(x - w / 2, inv["n_units"], width=w, color="#4a4a4a", label="units")
    ax2 = ax.twinx()
    ax2.bar(x + w / 2, inv["n_cells_used"], width=w, color=[DATASET_COLORS[c] for c in COHORTS], label="cells used")
    ax.set_xticks(x)
    ax.set_xticklabels(COHORTS, rotation=25, ha="right")
    ax.set_ylabel("n units")
    ax2.set_ylabel("n cells (used)")
    ax.set_title("four public sets after QC + cap", pad=3)
    despine(ax)
    ax2.spines["top"].set_visible(False)
    n_u = int(inv["n_units"].sum())
    n_c = int(inv["n_cells_used"].sum())
    ax.text(
        0.0,
        1.02,
        f"N={n_u} units · {n_c} cells used",
        transform=ax.transAxes,
        fontsize=6,
        va="bottom",
    )
    # hide twin spine top already; keep left/bottom only on ax
    ax.legend(loc="upper left", frameon=False, fontsize=5.5)
    ax2.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{int(v):,}"))

    # b — dataset
    ax = fig.add_subplot(gs[0, 1])
    panel_label(ax, "b")
    scatter_discrete(ax, xy, obs["dataset"].astype(str).to_numpy(), DATASET_COLORS, COHORTS, s=1.2)
    ax.set_title("Harmony UMAP · dataset", pad=3)
    legend_compact(ax, loc="upper right", fontsize=5.5)

    # c — annotation
    ax = fig.add_subplot(gs[0, 2])
    panel_label(ax, "c")
    scatter_discrete(ax, xy, obs["annotation"].astype(str).to_numpy(), ANNO_COLORS, COARSE_ORDER, s=1.2)
    ax.set_title("annotation", pad=3)
    legend_compact(ax, loc="upper right", fontsize=5.5)

    # d — dotplot
    ax = fig.add_subplot(gs[1, 0])
    panel_label(ax, "d", x=-0.18, y=1.08)
    genes = [g for g in DOT_GENES if adata.raw is not None and g in adata.raw.var_names]
    groups = [g for g in COARSE_ORDER if (adata.obs["annotation"].astype(str) == g).any()]
    frac = np.zeros((len(groups), len(genes)))
    mean = np.zeros((len(groups), len(genes)))
    for i, g in enumerate(groups):
        m = adata.obs["annotation"].astype(str).eq(g).to_numpy()
        for j, gene in enumerate(genes):
            v = _raw_vector(adata, gene)[m]
            frac[i, j] = float(np.mean(v > 0)) if v.size else 0.0
            mean[i, j] = float(np.mean(v)) if v.size else 0.0
    if mean.max() > 0:
        mean_z = mean / mean.max()
    else:
        mean_z = mean
    yy, xx = np.meshgrid(np.arange(len(groups)), np.arange(len(genes)), indexing="ij")
    sizes = 8 + 62 * frac
    sca = ax.scatter(
        xx.ravel(),
        yy.ravel(),
        s=sizes.ravel(),
        c=mean_z.ravel(),
        cmap="YlOrRd",
        vmin=0,
        vmax=1,
        edgecolors="0.35",
        linewidths=0.2,
        rasterized=True,
    )
    ax.set_xticks(np.arange(len(genes)))
    ax.set_xticklabels(genes, rotation=90, ha="center", fontsize=5.5)
    ax.set_yticks(np.arange(len(groups)))
    ax.set_yticklabels(groups)
    ax.set_xlim(-0.6, len(genes) - 0.4)
    ax.set_ylim(-0.6, len(groups) - 0.4)
    ax.invert_yaxis()
    ax.set_title("canonical markers", pad=3)
    despine(ax)
    cb = fig.colorbar(sca, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("mean expr (rel.)", fontsize=5.5)
    cb.ax.tick_params(labelsize=5)
    # size legend
    for fr, lab in ((0.2, "20%"), (0.6, "60%"), (1.0, "100%")):
        ax.scatter([], [], s=8 + 62 * fr, c="0.6", edgecolors="0.35", linewidths=0.2, label=lab)
    ax.legend(title="% pos", frameon=False, fontsize=5, title_fontsize=5, loc="lower right")

    # e — CLDN4, malignant on top
    ax = fig.add_subplot(gs[1, 1])
    panel_label(ax, "e")
    cldn4 = obs["cldn4"].to_numpy(float)
    vmax = float(np.quantile(cldn4[cldn4 > 0], 0.95)) if (cldn4 > 0).any() else 1.0
    vmax = max(vmax, 0.5)
    mal = obs["annotation"].astype(str).eq("malignant").to_numpy()
    ax.scatter(
        xy[~mal, 0],
        xy[~mal, 1],
        s=1.0,
        c=cldn4[~mal],
        cmap="YlOrRd",
        vmin=0,
        vmax=vmax,
        linewidths=0,
        alpha=0.45,
        rasterized=True,
    )
    sca = ax.scatter(
        xy[mal, 0],
        xy[mal, 1],
        s=2.4,
        c=cldn4[mal],
        cmap="YlOrRd",
        vmin=0,
        vmax=vmax,
        linewidths=0.12,
        edgecolors="0.15",
        alpha=0.95,
        rasterized=True,
    )
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    despine(ax)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_title("CLDN4  ·  malignant outlined", pad=3)
    cb = fig.colorbar(sca, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label("CLDN4", fontsize=5.5)
    cb.ax.tick_params(labelsize=5)

    # f — Leiden
    ax = fig.add_subplot(gs[1, 2])
    panel_label(ax, "f")
    cl = obs["leiden"].astype(str).to_numpy()
    clusters = sorted(pd.unique(cl), key=lambda x: int(x) if str(x).isdigit() else x)
    cmap = plt.get_cmap("tab20", max(len(clusters), 1))
    cols = {c: cmap(i % cmap.N) for i, c in enumerate(clusters)}
    scatter_discrete(ax, xy, cl, cols, clusters, s=1.2, alpha=0.85)
    ax.set_title(f"Leiden res={LEIDEN_RES}  ({len(clusters)} clusters)", pad=3)
    if len(clusters) <= 16:
        legend_compact(ax, loc="upper right", fontsize=4.5, markerscale=3.5, ncol=2)
    else:
        ax.text(0.02, 0.98, f"{len(clusters)} clusters", transform=ax.transAxes, va="top", fontsize=6)

    fig.savefig(path.with_suffix(".png"), dpi=600, facecolor="white")
    fig.savefig(path.with_suffix(".pdf"), dpi=600, facecolor="white")
    fig.savefig(path.with_suffix(".svg"), dpi=600, facecolor="white")
    plt.close(fig)
    print(f"wrote {path}.png/pdf/svg", flush=True)


def write_finding(inventory: pd.DataFrame, ann: pd.DataFrame, adata: ad.AnnData, infos: list[dict], plot_n: int, path: Path) -> None:
    inv = inventory.set_index("dataset")
    n_cells = int(inv["n_cells_used"].sum())
    n_units = int(inv["n_units"].sum())
    n_clust = int(adata.obs["leiden"].nunique())
    by = adata.obs["annotation"].astype(str).value_counts().to_dict()
    lines = [
        "# Concordant-4 cell-level atlas (methods / identity)",
        "",
        "ADDITIVE. **CLDN4-only.** This is the missing cell-level atlas figure",
        "for how the four public scRNA sets were taken in, integrated, clustered,",
        "and annotated. It is **methods/identity, not a new claim test.**",
        "",
        "Datasets ONLY: **GSE123902 + GSE131907 + GSE205335 + GSE189357**.",
        "Not GSE148071 / GSE127465 / GSE154826 / GSE207422. No dual-high.",
        "",
        "**PR #503 T/NK ρ and IFN/MHC numbers were not re-run and were not re-audited.**",
        "The locked four-set T/NK result remains n=65, %pos ρ=−0.531. Do not replace",
        "that patient/sample/donor number with a cell count from this atlas.",
        "",
        "## Honest n",
        "",
        f"- **n_units = {n_units}** (13 donors + 21 samples + 22 patients + 9 patients).",
        f"- **n_cells used (after QC + cap) = {n_cells}**.",
        f"- **n_clusters (Leiden {LEIDEN_RES}) = {n_clust}**.",
        f"- Plot used **{plot_n}** cells"
        + (
            f" (subsampled from {n_cells}; plot cap {PLOT_MAX})."
            if plot_n < n_cells
            else f" (all used cells; below the {PLOT_MAX} plot cap; no subsample)."
        ),
        "",
        "| dataset | unit | n_units | n_cells (author/raw in units) | n_cells QC | n_cells used | cap |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for c in COHORTS:
        r = inv.loc[c]
        lines.append(
            f"| {c} | {r['unit']} | {int(r['n_units'])} | {int(r['n_cells_input'])} | "
            f"{int(r['n_cells_qc'])} | {int(r['n_cells_used'])} | {int(r['cap_per_unit'])} |"
        )
    lines += [
        "",
        f"Annotation (used cells): {by}.",
        "",
        "## QC / Harmony / Leiden",
        "",
        f"- QC: n_genes ≥ {MIN_GENES}, n_counts ≥ {MIN_COUNTS}, mitochondrial % < {MAX_MITO}.",
        f"- Cap: ≤{CAP_PER_UNIT} cells / unit when the QC-pass set is larger (memory).",
        f"- Concatenate (inner gene join) → HVG {N_HVG} (Seurat flavor / seurat_v3 if `skmisc` is present, `batch_key=dataset`) → PCA {N_PCS}.",
        f"- Harmony: `batch=dataset`, theta={HARMONY_THETA}, max_iter=20. Sample was not a second Harmony key.",
        f"- Neighbors k={N_NEIGHBORS} on `X_pca_harmony`; UMAP; Leiden resolution **{LEIDEN_RES}**.",
        "- Public processed matrices only. Author 36.5 GB GSE123902 H5 and GSE131907 log2TPM text were not used.",
        "",
        "## Annotation rules",
        "",
        "Cluster → label. Do not invent cell types. Mixed clusters are labelled mixed and mapped to **other** on the UMAP.",
        "",
        "1. If ≥50% of cells in the cluster carry an **author malignant** label",
        "   (GSE131907 `Cell_subtype==Malignant cells`; GSE205335 `lineage.sub==Malignant cells`) → malignant (author).",
        "2. Else highest mean marker score on log-normalized expression:",
        "   EPCAM/KRT8/KRT18/KRT19/CLDN4 epithelium-malignant; CD3D/CD3E T; NKG7/GNLY/KLRD1 NK;",
        "   CD68/LYZ myeloid; MS4A1/CD79A B; PECAM1/VWF endothelial; COL1A1 fibroblast.",
        "3. If the top score is <0.15 → other. If the top two scores differ by <0.05 and both ≥0.15 → mixed (other).",
        "4. Endothelial and fibroblast clusters are **other** on the coarse UMAP; the fine label is in `annotation.tsv`.",
        "5. GSE123902 and GSE189357 have no author cell-type column; they use markers only.",
        "",
        "## Cluster table",
        "",
        "| cluster | n | label | fine | top markers | rule |",
        "|---|---:|---|---|---|---|",
    ]
    for r in ann.itertuples():
        lines.append(
            f"| {r.cluster} | {r.n_cells} | {r.label} | {r.label_fine} | {r.top_markers} | {r.rule} |"
        )
    lines += [
        "",
        "## What this is not",
        "",
        "- Not a re-audit of PR #503 T/NK ρ (n=65) or malignant IFN/MHC DE.",
        "- Not a dual-high TACSTD2∩CLDN4 object.",
        "- Not GSE148071 / GSE127465 / GSE154826 / GSE207422.",
        "- Not evidence that CLDN4 *causes* T/NK exclusion. Identity figure only.",
        "- Cell counts are not the inferential n. Inferential n for the thesis remains the PR #503 units.",
        "",
        "## Files",
        "",
        "- `results/tables/annotation.tsv` — cluster → label + top markers",
        "- `results/tables/inventory.tsv` — honest n per dataset",
        "- `results/figures/fig_atlas_umap.png` / `.svg` / `.pdf` — Figure 0 / Fig. atlas",
        "",
        "Reproduce:",
        "",
        "```bash",
        "python3 methods/concordant4_atlas_umap_annotate/download.py",
        "python3 methods/concordant4_atlas_umap_annotate/analyze.py",
        "```",
        "",
    ]
    path.write_text("\n".join(lines))
    print(f"wrote {path}", flush=True)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=Path, default=DEFAULT_DATA)
    p.add_argument("--cap", type=int, default=CAP_PER_UNIT)
    args = p.parse_args()
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    sc.settings.verbosity = 2
    sc.settings.figdir = str(FIGS)

    loaders = {
        "GSE123902": load_gse123902,
        "GSE131907": load_gse131907,
        "GSE205335": load_gse205335,
        "GSE189357": load_gse189357,
    }
    cache = args.data / "cache"
    cache.mkdir(parents=True, exist_ok=True)
    adatas = []
    infos = []
    for name, fn in loaders.items():
        h5 = cache / f"{name}.h5ad"
        js = cache / f"{name}.json"
        if h5.exists() and js.exists():
            print(f"==== load {name} (cache) ====", flush=True)
            a = ad.read_h5ad(h5)
            info = json.loads(js.read_text())
        else:
            print(f"==== load {name} ====", flush=True)
            a, info = fn(args.data, args.cap, rng)
            a.write_h5ad(h5)
            js.write_text(json.dumps(info, indent=2, default=str))
        adatas.append(a)
        infos.append(info)
        gc.collect()

    print("==== concat / Harmony / Leiden ====", flush=True)
    adata = align_and_concat(adatas)
    del adatas
    gc.collect()
    adata = integrate(adata)
    print("==== annotate ====", flush=True)
    ann = annotate(adata)
    TABLES.joinpath("annotation.tsv").write_text(ann.to_csv(sep="\t", index=False))

    inv_rows = []
    for info in infos:
        inv_rows.append(
            {
                "dataset": info["dataset"],
                "unit": info["unit"],
                "n_units": info["n_units"],
                "n_cells_input": info.get("n_cells_raw", info.get("n_cells_author_in_units", np.nan)),
                "n_cells_qc": info.get("n_cells_qc", info.get("n_cells_used")),
                "n_cells_used": info["n_cells_used"],
                "cap_per_unit": info.get("cap_per_unit", args.cap),
                "note": info.get("note", ""),
            }
        )
    inventory = pd.DataFrame(inv_rows)
    inventory.to_csv(TABLES / "inventory.tsv", sep="\t", index=False)
    adata.obs.to_csv(TABLES / "cell_obs.tsv", sep="\t")

    rng2 = np.random.default_rng(SEED)
    idx = subsample_plot(adata, rng2)
    plot_n = int(len(idx))
    print(f"==== figure (plot {plot_n} / {adata.n_obs}) ====", flush=True)
    make_figure(adata, inventory, ann, idx, FIGS / "fig_atlas_umap")

    summary = {
        "n_cells": int(adata.n_obs),
        "n_units": int(inventory["n_units"].sum()),
        "n_clusters": int(adata.obs["leiden"].nunique()),
        "n_plot": plot_n,
        "qc": {"min_genes": MIN_GENES, "min_counts": MIN_COUNTS, "max_mito": MAX_MITO},
        "harmony": {"batch": "dataset", "theta": HARMONY_THETA, "n_pcs": N_PCS, "n_hvg": N_HVG},
        "leiden_resolution": LEIDEN_RES,
        "pr503_tnk_not_rerun": True,
        "dual_high": False,
        "datasets": COHORTS,
        "annotation_counts": adata.obs["annotation"].astype(str).value_counts().to_dict(),
        "infos": infos,
    }
    (RESULTS / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    write_finding(inventory, ann, adata, infos, plot_n, HERE / "FINDING.md")
    print(json.dumps({k: summary[k] for k in ("n_cells", "n_units", "n_clusters", "n_plot")}, indent=2), flush=True)


if __name__ == "__main__":
    main()
