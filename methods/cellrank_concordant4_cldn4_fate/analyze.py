#!/usr/bin/env python3
"""CellRank 2 fate toward a barrier CLDN4-high state in concordant-4 epithelium.

Datasets ONLY: GSE123902 + GSE131907 + GSE205335 + GSE189357.
Not GSE148071, GSE127465, GSE154826, GSE200563, E-MTAB-13526, and not
public mouse scRNA. Does not re-estimate the locked patient-level
malignant CLDN4 %pos vs T/NK association (n=65, ρ=−0.531).

Population: tumor epithelial / author-malignant cells. Normal lung is
not mixed in as a root. Each dataset is embedded on its own
(Harmony on the unit key). No cross-dataset graph.

Terminals are pre-specified, not GPCCA-discovered:
  barrier_cldn4_high — CLDN4 at or above the dataset 80th percentile,
      barrier score (CLDN4 held out) at or above the median, and
      airway score below the 75th percentile.
  cycle_cldn4_low — competing sink: cycle score at or above the 80th
      percentile and CLDN4 at or below the median.
Root, only to orient diffusion pseudotime: AT2-like score at or above
the 75th percentile and CLDN4 at or below the median. If that root is
missing, or barrier cells are not downstream of it, the primary kernel
is CellRank ConnectivityKernel. Otherwise the primary kernel is the
pre-specified mixture 0.8 * PseudotimeKernel + 0.2 * ConnectivityKernel.

The inferential unit is the patient / donor / sample. The reported fate
is the mean absorption probability of non-terminal cells in that unit.
Terminal cells are excluded from the mean (their fate is 1 by construction).

A discrete kNN contact matrix and a multinomial logistic of neighbor
state are secondary, undirected descriptions of the same graph.
"""
from __future__ import annotations

import argparse
import faulthandler
import gc
import gzip
import json
import logging
import re
import shutil
import tarfile
import tempfile
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.io import mmread

import anndata as ad
import scanpy as sc

faulthandler.enable()
logging.getLogger("harmonypy").setLevel(logging.WARNING)

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
TABLES = RESULTS / "tables"
FIGS = RESULTS / "figures"
CACHE = Path("/tmp/epi_cache")
DEFAULT_DATA = Path("/tmp/geo_atlas")

COHORTS = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
TUMOR_ORIGINS = {"tLung", "tL/B", "mLN", "PE", "mBrain"}
UNIT_KIND = {
    "GSE123902": "donor",
    "GSE131907": "sample",
    "GSE205335": "patient",
    "GSE189357": "patient",
}
CAP_PER_UNIT = 400
MIN_GENES = 200
MIN_COUNTS = 500
MAX_MITO = 20.0
N_HVG = 2000
N_PCS = 30
N_NEIGHBORS = 15
SEED = 1
MIN_TERMINAL = 25
MIN_ROOT = 30
MIN_NONTERMINAL_UNIT = 15
MIN_WITHIN_UNIT = 80

AT2_GENES = ["SFTPC", "SFTPB", "SFTPA1", "SFTPA2", "NAPSA", "SLC34A2"]
BARRIER_WO = ["CLDN3", "CLDN7", "CDH1", "OCLN", "TJP1", "EPCAM", "TACSTD2", "KRT8", "KRT19"]
CYCLE_GENES = ["MKI67", "TOP2A", "STMN1", "PCNA", "HMGB2"]
AIRWAY_GENES = ["FOXJ1", "CAPS", "PIFO", "TPPP3", "KRT5"]
EPI_POS = ["EPCAM", "KRT8", "KRT18", "KRT19"]

NAME_123902 = re.compile(
    r"(GSM\d+)_(MSK_LX[^_]+(?:B)?)_(PRIMARY_TUMOUR|METASTASIS|NORMAL)_dense\.csv\.gz"
)
DATASET_COLORS = {
    "GSE123902": "#4c78a8",
    "GSE131907": "#f58518",
    "GSE205335": "#54a24b",
    "GSE189357": "#b279a2",
}
STATE_COLORS = {
    "at2_like": "#4c78a8",
    "intermediate": "#bdbdbd",
    "barrier_cldn4_high": "#c0392b",
    "cycle_cldn4_low": "#e67e22",
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
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.02,
        }
    )


def despine(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def mito_mask(genes: list[str]) -> np.ndarray:
    return np.array([g.startswith("MT-") or g.startswith("MT.") for g in genes], dtype=bool)


def genes_to_upper_unique(genes: list[str]) -> tuple[list[str], np.ndarray]:
    seen: dict[str, int] = {}
    keep: list[int] = []
    names: list[str] = []
    for i, g in enumerate(genes):
        u = str(g).strip().strip('"').upper()
        if not u or u in seen:
            continue
        seen[u] = i
        keep.append(i)
        names.append(u)
    return names, np.asarray(keep, dtype=int)


def qc_from_counts(X: sparse.spmatrix, genes: list[str]) -> pd.DataFrame:
    if not sparse.isspmatrix_csr(X):
        X = X.tocsr()
    n_counts = np.asarray(X.sum(axis=1)).ravel()
    n_genes = np.asarray((X > 0).sum(axis=1)).ravel()
    mt = mito_mask(genes)
    if mt.any():
        mito = np.asarray(X[:, mt].sum(axis=1)).ravel()
        pct = np.where(n_counts > 0, 100.0 * mito / n_counts, 0.0)
    else:
        pct = np.zeros(X.shape[0], dtype=float)
    return pd.DataFrame({"n_counts": n_counts, "n_genes": n_genes, "pct_mito": pct})


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
        keep.extend(int(i) for i in idx)
    return np.asarray(sorted(keep), dtype=int)


def reservoir_push(store: list, n_seen: int, item, cap: int, rng: np.random.Generator) -> None:
    if len(store) < cap:
        store.append(item)
        return
    j = int(rng.integers(0, n_seen))
    if j < cap:
        store[j] = item


def csr_from_dense_rows(rows: list[np.ndarray]) -> sparse.csr_matrix:
    if not rows:
        return sparse.csr_matrix((0, 0), dtype=np.float32)
    return sparse.csr_matrix(np.vstack(rows), dtype=np.float32)


def make_adata(X: sparse.spmatrix, obs: pd.DataFrame, genes: list[str]) -> ad.AnnData:
    obs = obs.copy()
    obs.index = pd.Index(obs["barcode"].astype(str))
    if not obs.index.is_unique:
        obs.index = pd.Index([f"{b}_{i}" for i, b in enumerate(obs.index)])
        obs["barcode"] = obs.index
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))
    adata.var_names_make_unique()
    return adata


def epi_mask_from_X(X: sparse.spmatrix, genes: list[str]) -> np.ndarray:
    gmap = {g: i for i, g in enumerate(genes)}
    pos_idx = [gmap[g] for g in EPI_POS if g in gmap]
    if not pos_idx:
        return np.zeros(X.shape[0], dtype=bool)
    pos = np.asarray(X[:, pos_idx].sum(axis=1)).ravel() > 0
    if "PTPRC" in gmap:
        neg = np.asarray(X[:, [gmap["PTPRC"]]].sum(axis=1)).ravel() > 0
    else:
        neg = np.zeros(X.shape[0], dtype=bool)
    return pos & ~neg


def dense_epi_ok(vals: np.ndarray, gmap: dict[str, int]) -> bool:
    if not any(vals[gmap[g]] > 0 for g in EPI_POS if g in gmap):
        return False
    if "PTPRC" in gmap and vals[gmap["PTPRC"]] > 0:
        return False
    return True


# ---------------------------------------------------------------------------
# Loaders — epithelial / malignant only
# ---------------------------------------------------------------------------
def load_gse123902(data: Path, cap: int, rng: np.random.Generator) -> tuple[ad.AnnData, dict]:
    tar = data / "GSE123902" / "GSE123902_RAW.tar"
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
    n_epi: dict[str, int] = {}
    n_files = 0
    for path in files:
        m = NAME_123902.match(path.name)
        if not m:
            print(f"  skip name {path.name}", flush=True)
            continue
        if m.group(3) == "NORMAL":
            continue
        raw_id = m.group(2)
        donor = raw_id[4:] if raw_id.startswith("MSK_") else raw_id
        tissue = "PRIMARY" if m.group(3) == "PRIMARY_TUMOUR" else "METASTASIS"
        n_files += 1
        reservoirs.setdefault(donor, [])
        n_raw.setdefault(donor, 0)
        n_qc.setdefault(donor, 0)
        n_epi.setdefault(donor, 0)
        with gzip.open(path, "rt") as handle:
            header = handle.readline().rstrip("\n").split(",")
            genes_u, keep_g = genes_to_upper_unique([g.strip() for g in header[1:]])
            gmap = {g: i for i, g in enumerate(genes_u)}
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
                if not dense_epi_ok(vals, gmap):
                    continue
                n_epi[donor] += 1
                item = (vals, genes_u, parts[0], m.group(1), tissue)
                reservoir_push(reservoirs[donor], n_epi[donor], item, cap, rng)
        print(
            f"  GSE123902 {path.name} donor={donor} qc={n_qc[donor]} epi={n_epi[donor]} kept={len(reservoirs[donor])}",
            flush=True,
        )
    gene_sets = [set(item[1]) for items in reservoirs.values() for item in items]
    if not gene_sets:
        raise RuntimeError("GSE123902: no epithelial tumor cells")
    genes = sorted(set.intersection(*gene_sets))
    gene_index = {g: i for i, g in enumerate(genes)}
    rows, obs_rows = [], []
    for donor, items in reservoirs.items():
        for vals, genes_u, barcode, gsm, tissue in items:
            aligned = np.zeros(len(genes), dtype=np.float32)
            lookup = np.fromiter((gene_index.get(g, -1) for g in genes_u), dtype=np.int32, count=len(genes_u))
            ok = lookup >= 0
            aligned[lookup[ok]] = vals[ok]
            rows.append(aligned)
            obs_rows.append(
                {
                    "barcode": f"GSE123902_{gsm}_{barcode}",
                    "dataset": "GSE123902",
                    "unit_id": donor,
                    "tissue": tissue,
                    "author_malignant": False,
                    "gate": "marker_EPCAM_KRT_PTPRC0",
                }
            )
    adata = make_adata(csr_from_dense_rows(rows), pd.DataFrame(obs_rows), genes)
    info = {
        "dataset": "GSE123902",
        "unit": "donor",
        "n_files_tumor": n_files,
        "n_units": int(adata.obs["unit_id"].nunique()),
        "n_cells_raw": int(sum(n_raw.values())),
        "n_cells_qc": int(sum(n_qc.values())),
        "n_cells_epi_qc": int(sum(n_epi.values())),
        "n_cells_used": int(adata.n_obs),
        "cap_per_unit": cap,
        "note": "Laughney 2020 tumor/met donors; NORMAL dropped; marker epithelium (EPCAM|KRT8|18|19)>0 and PTPRC==0",
    }
    print(f"GSE123902 used={adata.n_obs} units={info['n_units']} epi_qc={info['n_cells_epi_qc']}", flush=True)
    return adata, info


def stream_umi_subset(umi_path: Path, keep_ids: list[str]) -> tuple[list[str], list[str], sparse.csr_matrix]:
    keep_set = set(keep_ids)
    with gzip.open(umi_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        col_idx = np.flatnonzero(np.fromiter((c in keep_set for c in cell_ids), dtype=bool, count=len(cell_ids)))
        if col_idx.size == 0:
            raise RuntimeError("GSE131907: no requested cell IDs in UMI header")
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
                raise RuntimeError(f"GSE131907 row {gi} {gene}: {vals.size} != {len(cell_ids)}")
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


def load_gse131907(data: Path, cap: int, rng: np.random.Generator) -> tuple[ad.AnnData, dict]:
    ann = pd.read_csv(data / "GSE131907" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz", sep="\t")
    tumor = ann[ann["Sample_Origin"].isin(TUMOR_ORIGINS)].copy()
    mal = tumor[tumor["Cell_subtype"].astype(str).eq("Malignant cells")].copy()
    n_per = mal.groupby("Sample").size()
    keep_samples = set(n_per[n_per >= 20].index.astype(str))
    cells = mal[mal["Sample"].astype(str).isin(keep_samples)].copy()
    n_author = int(len(cells))
    ordered, genes_raw, mat = stream_umi_subset(
        data / "GSE131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
        cells["Index"].astype(str).tolist(),
    )
    genes, keep_g = genes_to_upper_unique(genes_raw)
    X = mat[keep_g, :].T.tocsr()
    del mat
    gc.collect()
    obs = cells.set_index("Index").reindex(ordered)
    qc = qc_from_counts(X, genes)
    ok = pass_qc(qc)
    X = X[ok]
    obs = obs.loc[np.asarray(ordered)[ok]].copy()
    keep_idx = cap_index(obs["Sample"].astype(str).to_numpy(), cap, rng)
    X = X[keep_idx]
    obs = obs.iloc[keep_idx].copy()
    obs["dataset"] = "GSE131907"
    obs["unit_id"] = obs["Sample"].astype(str)
    obs["tissue"] = obs["Sample_Origin"].astype(str)
    obs["author_malignant"] = True
    obs["gate"] = "author_Malignant_cells"
    obs["barcode"] = ["GSE131907_" + str(i) for i in obs.index]
    adata = make_adata(X, obs, genes)
    info = {
        "dataset": "GSE131907",
        "unit": "sample",
        "n_units": int(adata.obs["unit_id"].nunique()),
        "n_cells_raw": n_author,
        "n_cells_qc": int(ok.sum()),
        "n_cells_epi_qc": int(ok.sum()),
        "n_cells_used": int(adata.n_obs),
        "cap_per_unit": cap,
        "note": "Kim 2020 author Malignant cells in tumor-bearing samples with n>=20; nLung AT2 not included",
    }
    print(f"GSE131907 used={adata.n_obs} units={info['n_units']} author={n_author}", flush=True)
    return adata, info


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
    metadata["orig.ident"] = metadata["description"].str.replace("_", "-", regex=False) + "-" + read_end + "P"
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
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(int(v) for v in obj.Dim))
    return matrix, genes, barcodes


def load_gse205335(data: Path, cap: int, rng: np.random.Generator) -> tuple[ad.AnnData, dict]:
    ident = pd.read_csv(data / "GSE205335" / "GSE205335_Lung_IO_CellIdentity.txt.gz", sep="\t")
    meta = parse_geo_soft(data / "GSE205335" / "GSE205335_family.soft.gz")
    cells = ident.merge(
        meta[["orig.ident", "patient", "tissue", "recist", "cancer_subtype"]],
        on="orig.ident",
        how="left",
        validate="many_to_one",
    )
    if cells["patient"].isna().any():
        missing = cells.loc[cells["patient"].isna(), "orig.ident"].unique().tolist()[:8]
        raise RuntimeError(f"GSE205335 identity missing SOFT patient: {missing}")
    cells["is_normal_tissue"] = cells["tissue"].astype(str).str.startswith("Normal ")
    use = cells.loc[
        cells["lineage.sub"].astype(str).eq("Malignant cells") & ~cells["is_normal_tissue"]
    ].copy()
    n_author = int(len(use))
    n_per = use.groupby("patient").size()
    keep_patients = set(n_per[n_per >= 20].index.astype(str))
    use = use[use["patient"].astype(str).isin(keep_patients)].copy()
    keep = set(use["barcode"].astype(str))
    matrix, genes_raw, barcodes = load_rds_matrix(data / "GSE205335" / "GSE205335_Lung_IO_UMI_matrix.rds.gz")
    col_idx = np.flatnonzero(np.fromiter((b in keep for b in barcodes), dtype=bool, count=len(barcodes)))
    print(f"  GSE205335 subset {len(col_idx)} / {len(barcodes)}", flush=True)
    sub = matrix[:, col_idx].tocsr()
    del matrix
    gc.collect()
    ordered = [str(barcodes[i]) for i in col_idx]
    genes, keep_g = genes_to_upper_unique(list(genes_raw))
    X = sub[keep_g, :].T.tocsr()
    del sub
    obs = use.set_index("barcode").reindex(ordered)
    qc = qc_from_counts(X, genes)
    ok = pass_qc(qc)
    X = X[ok]
    obs = obs.loc[np.asarray(ordered)[ok]].copy()
    keep_idx = cap_index(obs["patient"].astype(str).to_numpy(), cap, rng)
    X = X[keep_idx]
    obs = obs.iloc[keep_idx].copy()
    obs["dataset"] = "GSE205335"
    obs["unit_id"] = obs["patient"].astype(str)
    obs["author_malignant"] = True
    obs["gate"] = "author_Malignant_cells"
    obs["barcode"] = ["GSE205335_" + str(i) for i in obs.index]
    adata = make_adata(X, obs, genes)
    info = {
        "dataset": "GSE205335",
        "unit": "patient",
        "n_units": int(adata.obs["unit_id"].nunique()),
        "n_cells_raw": n_author,
        "n_cells_qc": int(ok.sum()),
        "n_cells_epi_qc": int(ok.sum()),
        "n_cells_used": int(adata.n_obs),
        "cap_per_unit": cap,
        "note": "Ahn/Lee; author Malignant cells; normal tissue dropped; patients with n>=20",
    }
    print(f"GSE205335 used={adata.n_obs} units={info['n_units']} author={n_author}", flush=True)
    return adata, info


def _tar_members(tf: tarfile.TarFile) -> dict[str, tarfile.TarInfo]:
    return {m.name: m for m in tf.getmembers() if m.isfile()}


def _find_member(members: dict[str, tarfile.TarInfo], sample: str, token: str) -> str:
    hits = [n for n in members if f"_{sample}_" in Path(n).name and token in n.lower()]
    if not hits:
        raise FileNotFoundError(f"{sample}: no file matching {token}")
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
    per_qc, per_epi, per_raw = {}, {}, {}
    with tarfile.open(tar_path) as tf:
        members = _tar_members(tf)
        for sample in [f"TD{i}" for i in range(1, 10)]:
            feat_rows = _read_tsv_gz(tf, members[_find_member(members, sample, "features")])
            genes_raw = []
            for parts in feat_rows:
                if len(parts) >= 2 and not parts[1].startswith("ENSG"):
                    genes_raw.append(parts[1])
                else:
                    genes_raw.append(parts[0])
            barcodes = [r[0] for r in _read_tsv_gz(tf, members[_find_member(members, sample, "barcodes")])]
            with gzip.GzipFile(fileobj=tf.extractfile(members[_find_member(members, sample, "matrix.mtx")])) as handle:
                mat = mmread(handle).tocsc()
            if mat.shape[0] != len(genes_raw):
                if mat.shape[1] == len(genes_raw):
                    mat = mat.T.tocsc()
                else:
                    raise RuntimeError(f"{sample}: mtx {mat.shape} vs genes {len(genes_raw)}")
            genes, keep_g = genes_to_upper_unique(genes_raw)
            X = mat[keep_g, :].T.tocsr()
            del mat
            per_raw[sample] = int(X.shape[0])
            qc = qc_from_counts(X, genes)
            ok = pass_qc(qc)
            per_qc[sample] = int(ok.sum())
            X = X[ok]
            bc = np.asarray(barcodes)[ok]
            epi = epi_mask_from_X(X, genes)
            per_epi[sample] = int(epi.sum())
            X = X[epi]
            bc = bc[epi]
            if X.shape[0] > cap:
                pick = rng.choice(X.shape[0], size=cap, replace=False)
                pick.sort()
                X = X[pick]
                bc = bc[pick]
            obs = pd.DataFrame(
                {
                    "barcode": [f"GSE189357_{sample}_{b}" for b in bc],
                    "dataset": "GSE189357",
                    "unit_id": sample,
                    "tissue": "TUMOR",
                    "author_malignant": False,
                    "gate": "marker_EPCAM_KRT_PTPRC0",
                }
            )
            pieces.append(ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes)))
            print(f"  {sample}: raw={per_raw[sample]} qc={per_qc[sample]} epi={per_epi[sample]} used={X.shape[0]}", flush=True)
    if not pieces:
        raise RuntimeError("GSE189357: no samples")
    common = set(pieces[0].var_names)
    for a in pieces[1:]:
        common &= set(a.var_names)
    common_l = sorted(common)
    aligned = [a[:, common_l].copy() for a in pieces]
    adata = ad.concat(aligned, join="inner", merge="same")
    adata.obs_names_make_unique()
    info = {
        "dataset": "GSE189357",
        "unit": "patient",
        "n_units": int(adata.obs["unit_id"].nunique()),
        "n_cells_raw": int(sum(per_raw.values())),
        "n_cells_qc": int(sum(per_qc.values())),
        "n_cells_epi_qc": int(sum(per_epi.values())),
        "n_cells_used": int(adata.n_obs),
        "cap_per_unit": cap,
        "note": "Zhu 2022 TD1–TD9; no author cell type; marker epithelium; normal not in this matrix",
    }
    print(f"GSE189357 used={adata.n_obs} units={info['n_units']} epi_qc={info['n_cells_epi_qc']}", flush=True)
    return adata, info


LOADERS = {
    "GSE123902": load_gse123902,
    "GSE131907": load_gse131907,
    "GSE205335": load_gse205335,
    "GSE189357": load_gse189357,
}


# ---------------------------------------------------------------------------
# States, kernels, contact tracing
# ---------------------------------------------------------------------------
def gene_vector(adata: ad.AnnData, gene: str) -> np.ndarray:
    if gene not in adata.var_names:
        return np.zeros(adata.n_obs, dtype=np.float64)
    x = adata[:, gene].X
    if sparse.issparse(x):
        return np.asarray(x.todense()).ravel().astype(np.float64)
    return np.asarray(x).ravel().astype(np.float64)


def score_or_zero(adata: ad.AnnData, genes: list[str], name: str) -> list[str]:
    present = [g for g in genes if g in adata.var_names]
    if len(present) < 2:
        adata.obs[name] = 0.0
        return present
    sc.tl.score_genes(adata, present, score_name=name, random_state=SEED)
    return present


def assign_states(
    cldn4: np.ndarray,
    score_at2: np.ndarray,
    score_barrier: np.ndarray,
    score_cycle: np.ndarray,
    score_airway: np.ndarray,
) -> tuple[np.ndarray, dict]:
    n = len(cldn4)
    q80 = float(np.quantile(cldn4, 0.80))
    q50 = float(np.quantile(cldn4, 0.50))
    bmed = float(np.median(score_barrier))
    c80 = float(np.quantile(score_cycle, 0.80))
    a75 = float(np.quantile(score_at2, 0.75))
    air75 = float(np.quantile(score_airway, 0.75))
    info = {
        "cldn4_q80": q80,
        "cldn4_q50": q50,
        "barrier_median": bmed,
        "cycle_q80": c80,
        "at2_q75": a75,
        "airway_q75": air75,
        "cycle_relaxed": False,
        "airway_filter_relaxed": False,
        "estimable": True,
        "reason": "",
    }
    if not np.isfinite(q80) or np.sum(cldn4 > 0) < MIN_TERMINAL:
        info["estimable"] = False
        info["reason"] = "fewer than 25 cells express CLDN4"
        return np.array(["intermediate"] * n, dtype=object), info
    hi = cldn4 >= q80 if q80 > 0 else cldn4 > 0
    barrier = hi & (score_barrier >= bmed) & (score_airway < air75)
    if int(barrier.sum()) < MIN_TERMINAL:
        barrier = hi & (score_barrier >= bmed)
        info["airway_filter_relaxed"] = True
    if int(barrier.sum()) < MIN_TERMINAL:
        info["estimable"] = False
        info["reason"] = f"barrier terminal n={int(barrier.sum())} < {MIN_TERMINAL}"
        return np.array(["intermediate"] * n, dtype=object), info
    cycle = (score_cycle >= c80) & (cldn4 <= q50) & ~barrier
    if int(cycle.sum()) < MIN_TERMINAL:
        c70 = float(np.quantile(score_cycle, 0.70))
        cycle = (score_cycle >= c70) & (cldn4 <= q50) & ~barrier
        info["cycle_relaxed"] = True
        info["cycle_q70"] = c70
    if int(cycle.sum()) < MIN_TERMINAL:
        info["estimable"] = False
        info["reason"] = f"competing cycle terminal n={int(cycle.sum())} < {MIN_TERMINAL}"
        return np.array(["intermediate"] * n, dtype=object), info
    at2 = (score_at2 >= a75) & (cldn4 <= q50) & ~barrier & ~cycle
    state = np.array(["intermediate"] * n, dtype=object)
    state[at2] = "at2_like"
    state[cycle] = "cycle_cldn4_low"
    state[barrier] = "barrier_cldn4_high"
    info["n_barrier"] = int(barrier.sum())
    info["n_cycle"] = int(cycle.sum())
    info["n_at2"] = int(at2.sum())
    info["n_intermediate"] = int((state == "intermediate").sum())
    info["root_ok"] = int(at2.sum()) >= MIN_ROOT
    return state, info


def embed(adata: ad.AnnData) -> None:
    n_comps = int(min(N_PCS, adata.n_obs - 1, max(2, adata.n_vars - 1)))
    n_top = int(min(N_HVG, max(50, adata.n_vars // 2)))
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor="seurat")
    hvg = adata.var_names[adata.var["highly_variable"].to_numpy()]
    sub = adata[:, hvg].copy()
    sc.pp.scale(sub, max_value=10)
    sc.pp.pca(sub, n_comps=n_comps, random_state=SEED, zero_center=True)
    adata.obsm["X_pca"] = np.asarray(sub.obsm["X_pca"])
    del sub
    gc.collect()
    use = "X_pca"
    n_units = int(adata.obs["unit_id"].nunique())
    if n_units >= 2 and adata.n_obs >= 50:
        import harmonypy

        ho = harmonypy.run_harmony(
            adata.obsm["X_pca"],
            adata.obs,
            ["unit_id"],
            theta=2.0,
            max_iter_harmony=20,
            random_state=SEED,
            verbose=False,
            ncores=1,
        )
        Z = np.asarray(ho.Z_corr)
        if Z.shape == adata.obsm["X_pca"].shape:
            adata.obsm["X_pca_harmony"] = Z
        elif Z.T.shape == adata.obsm["X_pca"].shape:
            adata.obsm["X_pca_harmony"] = Z.T
        else:
            raise RuntimeError(f"Harmony Z_corr shape {Z.shape} vs pca {adata.obsm['X_pca'].shape}")
        use = "X_pca_harmony"
    k = int(min(N_NEIGHBORS, max(5, adata.n_obs // 10)))
    k = min(k, adata.n_obs - 1)
    sc.pp.neighbors(adata, use_rep=use, n_neighbors=k, random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    adata.uns["embed_use"] = use
    adata.uns["n_neighbors"] = k


def orient_pseudotime(adata: ad.AnnData) -> dict:
    out = {"orientation_ok": False, "root_mode": "none", "frac_dpt_finite": 0.0}
    state = adata.obs["state"].astype(str).to_numpy()
    if int((state == "at2_like").sum()) < MIN_ROOT:
        out["reason"] = "AT2-like root below 30 cells; ConnectivityKernel only"
        return out
    at2_idx = np.flatnonzero(state == "at2_like")
    scores = adata.obs["score_at2"].to_numpy()
    iroot = int(at2_idx[np.argmax(scores[at2_idx])])
    adata.uns["iroot"] = iroot
    sc.tl.diffmap(adata)
    sc.tl.dpt(adata)
    dpt = adata.obs["dpt_pseudotime"].to_numpy(dtype=float)
    finite = np.isfinite(dpt)
    out["frac_dpt_finite"] = float(finite.mean())
    if finite.mean() < 0.95:
        out["reason"] = f"DPT finite fraction {finite.mean():.3f} < 0.95"
        return out
    def _mean(lab: str) -> float:
        m = (state == lab) & finite
        return float(dpt[m].mean()) if m.any() else float("nan")
    out["dpt_at2"] = _mean("at2_like")
    out["dpt_intermediate"] = _mean("intermediate")
    out["dpt_barrier"] = _mean("barrier_cldn4_high")
    out["dpt_cycle"] = _mean("cycle_cldn4_low")
    if not (out["dpt_barrier"] > out["dpt_at2"] + 0.02):
        out["reason"] = "barrier terminal is not downstream of the AT2-like root"
        return out
    out["orientation_ok"] = True
    out["root_mode"] = "at2_like"
    out["reason"] = "DPT from highest AT2-like cell; barrier mean DPT exceeds AT2"
    return out


def _fate_barrier(kernel, adata: ad.AnnData) -> np.ndarray:
    import cellrank as cr

    est = cr.estimators.GPCCA(kernel)
    terminals = {
        "barrier_cldn4_high": list(adata.obs_names[adata.obs["state"].astype(str).eq("barrier_cldn4_high")]),
        "cycle_cldn4_low": list(adata.obs_names[adata.obs["state"].astype(str).eq("cycle_cldn4_low")]),
    }
    est.set_terminal_states(terminals)
    errors = []
    for solver, tol, check in (("gmres", 1e-8, 1e-3), ("direct", 1e-8, 1e-2)):
        try:
            est.compute_fate_probabilities(
                use_petsc=False,
                n_jobs=1,
                backend="threading",
                show_progress_bar=False,
                solver=solver,
                tol=tol,
                check_sum_tol=check,
            )
            arr = np.asarray(est.fate_probabilities["barrier_cldn4_high"]).ravel().astype(np.float64)
            return arr
        except Exception as exc:  # solver failure is a recorded fallback, not a silent pass
            errors.append(f"{solver}: {type(exc).__name__}: {exc}")
            print(f"  fate solver failed ({solver}): {exc}", flush=True)
    raise RuntimeError(" | ".join(errors))


def run_kernels(adata: ad.AnnData, orient: dict) -> dict:
    import cellrank as cr

    ck = cr.kernels.ConnectivityKernel(adata)
    ck.compute_transition_matrix()
    fate_c = _fate_barrier(ck, adata)
    result = {
        "fate_connectivity": fate_c,
        "fate_combined": np.full(adata.n_obs, np.nan),
        "primary": "connectivity",
        "kernel_note": orient.get("reason", ""),
    }
    if orient.get("orientation_ok"):
        pk = cr.kernels.PseudotimeKernel(adata, time_key="dpt_pseudotime")
        pk.compute_transition_matrix(show_progress_bar=False, n_jobs=1, backend="threading")
        combined = 0.8 * pk + 0.2 * ck
        result["fate_combined"] = _fate_barrier(combined, adata)
        result["primary"] = "pseudotime_0.8_connectivity_0.2"
        result["kernel_note"] = "0.8 PseudotimeKernel + 0.2 ConnectivityKernel; " + orient.get("reason", "")
    return result


def gpcca_sensitivity(adata: ad.AnnData) -> dict:
    """Macrostates from the connectivity kernel. CLDN4 names the basin after the fact."""
    import cellrank as cr

    out: dict = {"ok": False}
    try:
        ck = cr.kernels.ConnectivityKernel(adata)
        ck.compute_transition_matrix()
        g = cr.estimators.GPCCA(ck)
        n_states = 4 if adata.n_obs >= 400 else 3
        g.compute_macrostates(n_states=n_states, n_cells=min(30, max(10, adata.n_obs // 50)))
        mem = g.macrostates_memberships
        cldn = adata.obs["cldn4"].to_numpy()
        barrier = adata.obs["score_barrier"].to_numpy()
        rows = []
        best_name, best_cldn = None, -np.inf
        for name in list(mem.names):
            w = np.asarray(mem[name]).ravel().astype(np.float64)
            w = np.clip(w, 0, None)
            if w.sum() <= 0:
                continue
            cmean = float(np.average(cldn, weights=w))
            bmean = float(np.average(barrier, weights=w))
            rows.append({"macrostate": str(name), "wmean_cldn4": cmean, "wmean_barrier_wo": bmean})
            if cmean > best_cldn:
                best_cldn = cmean
                best_name = str(name)
        g.set_terminal_states(states=None)
        g.compute_fate_probabilities(
            use_petsc=False, n_jobs=1, backend="threading", show_progress_bar=False, solver="direct", check_sum_tol=1e-2
        )
        fp = g.fate_probabilities
        fate = np.asarray(fp[best_name]).ravel().astype(np.float64)
        hard = g.macrostates.astype(str)
        non = hard.to_numpy() != best_name
        out.update(
            {
                "ok": True,
                "n_macrostates": n_states,
                "top_cldn4_macrostate": best_name,
                "macrostates": rows,
                "mean_fate_nonmember": float(np.nanmean(fate[non])) if non.any() else float("nan"),
                "fate": fate,
            }
        )
    except Exception as exc:
        out["reason"] = f"{type(exc).__name__}: {exc}"
        print(f"  GPCCA sensitivity skipped: {out['reason']}", flush=True)
    return out


def contact_matrix(states: np.ndarray, connectivities: sparse.spmatrix) -> tuple[np.ndarray, list[str]]:
    labs = ["at2_like", "intermediate", "barrier_cldn4_high", "cycle_cldn4_low"]
    C = connectivities.tocsr().astype(np.float64).copy()
    C.setdiag(0)
    C.eliminate_zeros()
    deg = np.asarray(C.sum(axis=1)).ravel()
    inv = np.zeros_like(deg)
    pos = deg > 0
    inv[pos] = 1.0 / deg[pos]
    C = sparse.diags(inv) @ C
    index = {s: np.flatnonzero(states == s) for s in labs}
    T = np.zeros((len(labs), len(labs)))
    for i, s in enumerate(labs):
        rows = index[s]
        rows = rows[deg[rows] > 0] if len(rows) else rows
        if len(rows) == 0:
            continue
        block = C[rows]
        for j, tname in enumerate(labs):
            cols = index[tname]
            if len(cols) == 0:
                continue
            T[i, j] = float(np.asarray(block[:, cols].sum()) / len(rows))
        ssum = T[i].sum()
        if ssum > 0:
            T[i] /= ssum
    return T, labs


def absorption_from_T(T: np.ndarray, labs: list[str]) -> dict[str, float]:
    transient = [s for s in ("at2_like", "intermediate") if s in labs]
    absorbing = ["barrier_cldn4_high", "cycle_cldn4_low"]
    # rebuild in transient-then-absorbing order
    order = transient + absorbing
    ix = [labs.index(s) for s in order]
    M = T[np.ix_(ix, ix)]
    nt = len(transient)
    if nt == 0:
        return {}
    Q = M[:nt, :nt]
    R = M[:nt, nt:]
    try:
        B = np.linalg.solve(np.eye(nt) - Q, R)
    except np.linalg.LinAlgError:
        return {}
    target = absorbing.index("barrier_cldn4_high")
    return {transient[i]: float(B[i, target]) for i in range(nt)}


def multinomial_contact(adata: ad.AnnData, rng: np.random.Generator) -> pd.DataFrame:
    from sklearn.linear_model import LogisticRegression

    C = adata.obsp["connectivities"].tocsr().astype(np.float64).copy()
    C.setdiag(0)
    C.eliminate_zeros()
    states = adata.obs["state"].astype(str).to_numpy()
    Xcols = np.column_stack(
        [
            adata.obs["score_at2"].to_numpy(),
            adata.obs["score_barrier"].to_numpy(),
            adata.obs["score_cycle"].to_numpy(),
        ]
    )
    # CLDN4 is not a predictor.
    ys = []
    xs = []
    indptr, indices, data = C.indptr, C.indices, C.data
    for i in range(adata.n_obs):
        a, b = indptr[i], indptr[i + 1]
        if a == b:
            continue
        w = data[a:b].astype(np.float64)
        wsum = w.sum()
        if wsum <= 0:
            continue
        j = int(indices[a + int(rng.choice(b - a, p=w / wsum))])
        ys.append(states[j])
        xs.append(Xcols[i])
    if len(ys) < 50:
        return pd.DataFrame()
    y = np.asarray(ys)
    X = np.asarray(xs)
    if len(np.unique(y)) < 2:
        return pd.DataFrame()
    clf = LogisticRegression(solver="lbfgs", max_iter=500)
    clf.fit(X, y)
    rows = []
    for ci, cls in enumerate(clf.classes_):
        for fi, feat in enumerate(["score_at2", "score_barrier_wo_cldn4", "score_cycle"]):
            rows.append({"class": str(cls), "feature": feat, "coef": float(clf.coef_[ci, fi]), "n_contacts": int(len(y))})
    return pd.DataFrame(rows)


def prepare_scores(adata: ad.AnnData) -> dict[str, list[str]]:
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    used = {
        "at2": score_or_zero(adata, AT2_GENES, "score_at2"),
        "barrier_wo": score_or_zero(adata, BARRIER_WO, "score_barrier"),
        "cycle": score_or_zero(adata, CYCLE_GENES, "score_cycle"),
        "airway": score_or_zero(adata, AIRWAY_GENES, "score_airway"),
    }
    adata.obs["cldn4"] = gene_vector(adata, "CLDN4")
    return used


def run_dataset(adata: ad.AnnData, rng: np.random.Generator) -> tuple[pd.DataFrame, dict, pd.DataFrame, pd.DataFrame]:
    genes_used = prepare_scores(adata)
    state, st_info = assign_states(
        adata.obs["cldn4"].to_numpy(),
        adata.obs["score_at2"].to_numpy(),
        adata.obs["score_barrier"].to_numpy(),
        adata.obs["score_cycle"].to_numpy(),
        adata.obs["score_airway"].to_numpy(),
    )
    adata.obs["state"] = pd.Categorical(state)
    summary = {
        "dataset": str(adata.obs["dataset"].iloc[0]),
        "n_cells": int(adata.n_obs),
        "n_units": int(adata.obs["unit_id"].nunique()),
        "genes_used": genes_used,
        "states": st_info,
        "estimable": bool(st_info.get("estimable")),
    }
    empty_cells = pd.DataFrame()
    empty_coef = pd.DataFrame()
    empty_T = pd.DataFrame()
    if not st_info.get("estimable"):
        summary["reason"] = st_info.get("reason", "not estimable")
        print(f"  {summary['dataset']} not estimable: {summary['reason']}", flush=True)
        return empty_cells, summary, empty_coef, empty_T
    embed(adata)
    orient = orient_pseudotime(adata)
    summary["orientation"] = orient
    print(
        f"  {summary['dataset']} states barrier={st_info['n_barrier']} cycle={st_info['n_cycle']} "
        f"at2={st_info['n_at2']} intermediate={st_info['n_intermediate']} "
        f"orient={orient.get('orientation_ok')} ({orient.get('reason')})",
        flush=True,
    )
    kernels = run_kernels(adata, orient)
    summary["primary_kernel"] = kernels["primary"]
    summary["kernel_note"] = kernels["kernel_note"]
    primary = kernels["fate_combined"] if kernels["primary"].startswith("pseudotime") else kernels["fate_connectivity"]
    adata.obs["fate_primary"] = primary
    adata.obs["fate_connectivity"] = kernels["fate_connectivity"]
    adata.obs["fate_combined"] = kernels["fate_combined"]
    gpcca = gpcca_sensitivity(adata)
    summary["gpcca"] = {k: v for k, v in gpcca.items() if k != "fate"}
    if gpcca.get("ok"):
        adata.obs["fate_gpcca_topcldn4"] = gpcca["fate"]
    else:
        adata.obs["fate_gpcca_topcldn4"] = np.nan
    states_arr = adata.obs["state"].astype(str).to_numpy()
    T, labs = contact_matrix(states_arr, adata.obsp["connectivities"])
    absorb = absorption_from_T(T, labs)
    summary["contact_absorption"] = absorb
    T_rows = []
    for i, src in enumerate(labs):
        for j, dst in enumerate(labs):
            T_rows.append({"dataset": summary["dataset"], "from_state": src, "to_state": dst, "contact_prob": float(T[i, j])})
    coef = multinomial_contact(adata, rng)
    if len(coef):
        coef.insert(0, "dataset", summary["dataset"])
    within = within_unit_table(adata, rng)
    summary["within_unit_n"] = int(within["within_unit_cr_fate"].notna().sum()) if len(within) else 0
    umap = np.asarray(adata.obsm["X_umap"])
    non = ~adata.obs["state"].astype(str).isin(["barrier_cldn4_high", "cycle_cldn4_low"])
    cells = pd.DataFrame(
        {
            "dataset": adata.obs["dataset"].astype(str).to_numpy(),
            "unit_id": adata.obs["unit_id"].astype(str).to_numpy(),
            "state": adata.obs["state"].astype(str).to_numpy(),
            "nonterminal": non.to_numpy(),
            "cldn4": adata.obs["cldn4"].to_numpy(),
            "score_at2": adata.obs["score_at2"].to_numpy(),
            "score_barrier": adata.obs["score_barrier"].to_numpy(),
            "score_cycle": adata.obs["score_cycle"].to_numpy(),
            "score_airway": adata.obs["score_airway"].to_numpy(),
            "dpt_pseudotime": adata.obs["dpt_pseudotime"].to_numpy() if "dpt_pseudotime" in adata.obs else np.nan,
            "fate_primary": adata.obs["fate_primary"].to_numpy(),
            "fate_connectivity": adata.obs["fate_connectivity"].to_numpy(),
            "fate_combined": adata.obs["fate_combined"].to_numpy(),
            "fate_gpcca_topcldn4": adata.obs["fate_gpcca_topcldn4"].to_numpy(),
            "umap1": umap[:, 0],
            "umap2": umap[:, 1],
        }
    )
    cells = cells.merge(within, on="unit_id", how="left")
    # composition checks
    barrier = cells["state"].eq("barrier_cldn4_high")
    summary["mean_cldn4_barrier"] = float(cells.loc[barrier, "cldn4"].mean())
    summary["mean_cldn4_cycle"] = float(cells.loc[cells["state"].eq("cycle_cldn4_low"), "cldn4"].mean())
    summary["mean_barrier_score_barrier"] = float(cells.loc[barrier, "score_barrier"].mean())
    summary["mean_barrier_score_rest"] = float(cells.loc[~barrier, "score_barrier"].mean())
    air_cut = float(st_info["airway_q75"])
    summary["airway_frac_in_barrier"] = float((cells.loc[barrier, "score_airway"] >= air_cut).mean()) if barrier.any() else float("nan")
    return cells, summary, coef, pd.DataFrame(T_rows)


def within_unit_table(adata: ad.AnnData, rng: np.random.Generator) -> pd.DataFrame:
    """CellRank connectivity fate inside one unit. No cross-unit edges."""
    rows = []
    for unit, idx in adata.obs.groupby("unit_id", observed=True).groups.items():
        rec = {
            "unit_id": str(unit),
            "within_n": int(len(idx)),
            "within_unit_cr_fate": np.nan,
            "within_unit_contact_absorb": np.nan,
            "within_status": "skipped_n<80",
        }
        if len(idx) < MIN_WITHIN_UNIT:
            rows.append(rec)
            continue
        sub = adata[list(idx)].copy()
        try:
            state, info = assign_states(
                sub.obs["cldn4"].to_numpy(),
                sub.obs["score_at2"].to_numpy(),
                sub.obs["score_barrier"].to_numpy(),
                sub.obs["score_cycle"].to_numpy(),
                sub.obs["score_airway"].to_numpy(),
            )
            if not info.get("estimable"):
                rec["within_status"] = info.get("reason", "not estimable")
                rows.append(rec)
                continue
            sub.obs["state"] = pd.Categorical(state)
            n_comps = int(min(15, sub.n_obs // 5, sub.n_vars - 1))
            n_comps = max(n_comps, 2)
            n_top = int(min(800, max(50, sub.n_vars // 3)))
            sc.pp.highly_variable_genes(sub, n_top_genes=n_top, flavor="seurat")
            hvg = sub.var_names[sub.var["highly_variable"].to_numpy()]
            piece = sub[:, hvg].copy()
            sc.pp.scale(piece, max_value=10)
            sc.pp.pca(piece, n_comps=n_comps, random_state=SEED)
            sub.obsm["X_pca"] = np.asarray(piece.obsm["X_pca"])
            del piece
            k = int(min(10, max(5, sub.n_obs // 15)))
            k = min(k, sub.n_obs - 1)
            sc.pp.neighbors(sub, use_rep="X_pca", n_neighbors=k, random_state=SEED)
            st = sub.obs["state"].astype(str).to_numpy()
            import cellrank as cr

            ck = cr.kernels.ConnectivityKernel(sub)
            ck.compute_transition_matrix()
            fate = _fate_barrier(ck, sub)
            non = ~np.isin(st, ["barrier_cldn4_high", "cycle_cldn4_low"])
            if int(non.sum()) >= 10:
                rec["within_unit_cr_fate"] = float(np.nanmean(fate[non]))
            T, labs = contact_matrix(st, sub.obsp["connectivities"])
            absorb = absorption_from_T(T, labs)
            if "intermediate" in absorb:
                rec["within_unit_contact_absorb"] = absorb["intermediate"]
            elif "at2_like" in absorb:
                rec["within_unit_contact_absorb"] = absorb["at2_like"]
            rec["within_status"] = "ok"
        except Exception as exc:
            rec["within_status"] = f"{type(exc).__name__}: {exc}"
        rows.append(rec)
        del sub
    return pd.DataFrame(rows)


def unit_summary(cells: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (dataset, unit), g in cells.groupby(["dataset", "unit_id"], observed=True):
        non = g[g["nonterminal"]]
        rec = {
            "dataset": dataset,
            "unit_id": unit,
            "n_cells": int(len(g)),
            "n_nonterminal": int(len(non)),
            "n_barrier": int(g["state"].eq("barrier_cldn4_high").sum()),
            "n_cycle": int(g["state"].eq("cycle_cldn4_low").sum()),
            "n_at2": int(g["state"].eq("at2_like").sum()),
            "mean_fate_primary": float(non["fate_primary"].mean()) if len(non) else np.nan,
            "mean_fate_connectivity": float(non["fate_connectivity"].mean()) if len(non) else np.nan,
            "mean_fate_combined": float(non["fate_combined"].mean()) if len(non) and non["fate_combined"].notna().any() else np.nan,
            "mean_fate_gpcca": float(non["fate_gpcca_topcldn4"].mean()) if len(non) and non["fate_gpcca_topcldn4"].notna().any() else np.nan,
            "within_unit_cr_fate": float(g["within_unit_cr_fate"].iloc[0]) if "within_unit_cr_fate" in g and pd.notna(g["within_unit_cr_fate"].iloc[0]) else np.nan,
            "within_unit_contact_absorb": float(g["within_unit_contact_absorb"].iloc[0]) if "within_unit_contact_absorb" in g and pd.notna(g["within_unit_contact_absorb"].iloc[0]) else np.nan,
        }
        rec["included"] = bool(rec["n_nonterminal"] >= MIN_NONTERMINAL_UNIT and np.isfinite(rec["mean_fate_primary"]))
        rows.append(rec)
    return pd.DataFrame(rows)


def wilcoxon_vs_half(x: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    if len(x) < 5 or np.allclose(x, x[0]):
        return float("nan")
    try:
        return float(stats.wilcoxon(x - 0.5, alternative="two-sided").pvalue)
    except ValueError:
        return float("nan")


def dl_meta(effects: np.ndarray, variances: np.ndarray) -> dict:
    effects = np.asarray(effects, dtype=float)
    variances = np.asarray(variances, dtype=float)
    ok = np.isfinite(effects) & np.isfinite(variances) & (variances > 0)
    effects, variances = effects[ok], variances[ok]
    k = len(effects)
    if k == 0:
        return {"k": 0}
    w = 1.0 / variances
    mu = float(np.sum(w * effects) / np.sum(w))
    Q = float(np.sum(w * (effects - mu) ** 2))
    df = k - 1
    C = float(np.sum(w) - np.sum(w ** 2) / np.sum(w)) if k > 1 else np.nan
    tau2 = float(max(0.0, (Q - df) / C)) if k > 1 and C > 0 else 0.0
    w2 = 1.0 / (variances + tau2)
    mu2 = float(np.sum(w2 * effects) / np.sum(w2))
    se = float(np.sqrt(1.0 / np.sum(w2)))
    z = mu2 / se if se > 0 else float("nan")
    p = float(2 * stats.norm.sf(abs(z))) if np.isfinite(z) else float("nan")
    I2 = float(max(0.0, (Q - df) / Q)) if Q > 0 and k > 1 else 0.0
    return {"k": k, "mu": mu2, "se": se, "z": z, "p": p, "I2": I2, "tau2": tau2, "Q": Q}


def cohort_rows(units: pd.DataFrame, summaries: list[dict]) -> tuple[pd.DataFrame, dict]:
    by_ds = {s["dataset"]: s for s in summaries}
    rows = []
    effects, variances, labels = [], [], []
    for ds in COHORTS:
        g = units[(units["dataset"] == ds) & units["included"]] if len(units) else pd.DataFrame()
        s = by_ds.get(ds, {})
        if len(g) == 0:
            rows.append({"dataset": ds, "estimable": False, "reason": s.get("reason", s.get("states", {}).get("reason", "no included units"))})
            continue
        x = g["mean_fate_primary"].to_numpy()
        mean = float(np.mean(x))
        sd = float(np.std(x, ddof=1)) if len(x) > 1 else float("nan")
        se = sd / np.sqrt(len(x)) if np.isfinite(sd) else float("nan")
        rec = {
            "dataset": ds,
            "estimable": True,
            "unit": UNIT_KIND[ds],
            "primary_kernel": s.get("primary_kernel"),
            "kernel_note": s.get("kernel_note"),
            "n_units_included": int(len(g)),
            "n_units_loaded": int(s.get("n_units", g["unit_id"].nunique())),
            "n_cells": int(s.get("n_cells", 0)),
            "median_unit_fate": float(np.median(x)),
            "mean_unit_fate": mean,
            "sd_unit_fate": sd,
            "se_unit_fate": se,
            "mean_minus_0.5": mean - 0.5,
            "frac_units_gt_0.5": float(np.mean(x > 0.5)),
            "wilcoxon_p_vs_0.5": wilcoxon_vs_half(x),
            "mean_fate_connectivity_units": float(g["mean_fate_connectivity"].mean()),
            "mean_within_unit_cr": float(g["within_unit_cr_fate"].mean()) if g["within_unit_cr_fate"].notna().any() else float("nan"),
            "n_within_unit": int(g["within_unit_cr_fate"].notna().sum()),
            "mean_within_contact": float(g["within_unit_contact_absorb"].mean()) if g["within_unit_contact_absorb"].notna().any() else float("nan"),
            "contact_absorb_intermediate": (s.get("contact_absorption") or {}).get("intermediate", float("nan")),
            "contact_absorb_at2": (s.get("contact_absorption") or {}).get("at2_like", float("nan")),
            "gpcca_mean_fate_nonmember": (s.get("gpcca") or {}).get("mean_fate_nonmember", float("nan")),
            "airway_frac_in_barrier": s.get("airway_frac_in_barrier"),
            "mean_cldn4_barrier": s.get("mean_cldn4_barrier"),
            "mean_cldn4_cycle": s.get("mean_cldn4_cycle"),
            "mean_barrier_score_barrier": s.get("mean_barrier_score_barrier"),
            "mean_barrier_score_rest": s.get("mean_barrier_score_rest"),
            "n_barrier": (s.get("states") or {}).get("n_barrier"),
            "n_cycle": (s.get("states") or {}).get("n_cycle"),
            "n_at2": (s.get("states") or {}).get("n_at2"),
            "n_intermediate": (s.get("states") or {}).get("n_intermediate"),
            "orientation_ok": (s.get("orientation") or {}).get("orientation_ok"),
            "dpt_at2": (s.get("orientation") or {}).get("dpt_at2"),
            "dpt_barrier": (s.get("orientation") or {}).get("dpt_barrier"),
            "dpt_cycle": (s.get("orientation") or {}).get("dpt_cycle"),
        }
        rows.append(rec)
        if len(x) >= 4 and np.isfinite(se) and se > 0:
            effects.append(mean - 0.5)
            variances.append(se ** 2)
            labels.append(ds)
    meta = dl_meta(np.asarray(effects), np.asarray(variances))
    meta["datasets"] = labels
    return pd.DataFrame(rows), meta


def fmt_p(p) -> str:
    try:
        p = float(p)
    except (TypeError, ValueError):
        return "NA"
    if not np.isfinite(p):
        return "NA"
    if p < 1e-3:
        return f"{p:.2e}"
    return f"{p:.4f}"


def fmt(x, nd=3) -> str:
    try:
        x = float(x)
    except (TypeError, ValueError):
        return "NA"
    if not np.isfinite(x):
        return "NA"
    return f"{x:.{nd}f}"


def write_finding(units: pd.DataFrame, cohort: pd.DataFrame, meta: dict, inventory: pd.DataFrame, versions: dict) -> None:
    included = units[units["included"]] if len(units) else units
    lines = []
    lines.append("# CellRank 2 fate toward barrier CLDN4-high (concordant-4)")
    lines.append("")
    lines.append("ADDITIVE. **CLDN4-only.** This does not replace the locked concordant-4 patient result")
    lines.append("(malignant CLDN4 %pos vs T/NK, n=65, ρ=−0.531, P=1.65×10⁻⁵, I²=0).")
    lines.append("Cell counts below are not that n.")
    lines.append("")
    lines.append("Question: in tumor epithelium from the same four public scRNA cohorts,")
    lines.append("what is the CellRank 2 absorption probability of cells that are not already")
    lines.append("in a terminal state, toward a pre-specified **barrier CLDN4-high** sink,")
    lines.append("when the competing sink is cycle-high and CLDN4-low?")
    lines.append("")
    lines.append("## Population")
    lines.append("")
    lines.append("- GSE123902 (Laughney): tumor/metastasis donors; marker epithelium `(EPCAM|KRT8|KRT18|KRT19)>0` and `PTPRC==0`. NORMAL files dropped.")
    lines.append("- GSE131907 (Kim): author `Cell_subtype==Malignant cells` in tumor-bearing samples with ≥20 malignant cells. nLung AT2/AT1/ciliated cells are not in the graph.")
    lines.append("- GSE205335 (Ahn/Lee): author `lineage.sub==Malignant cells`, normal tissue dropped, patients with ≥20 malignant cells.")
    lines.append("- GSE189357 (Zhu): TD1–TD9, same marker epithelium gate. No author cell-type column.")
    lines.append("- Not included: GSE148071, GSE127465, GSE154826, GSE200563, E-MTAB-13526, public mouse, private 8KL.")
    lines.append("- Cap: ≤400 epithelial cells per unit after QC (n_genes≥200, n_UMI≥500, mito%<20), seed=1.")
    lines.append("- Each dataset is its own graph. Harmony corrects the unit key inside the dataset. Datasets are not concatenated.")
    lines.append("")
    lines.append("## State definitions (pre-specified)")
    lines.append("")
    lines.append("- Barrier score **holds CLDN4 out**: CLDN3, CLDN7, CDH1, OCLN, TJP1, EPCAM, TACSTD2, KRT8, KRT19.")
    lines.append("- `barrier_cldn4_high`: CLDN4 ≥ dataset 80th percentile (or CLDN4>0 if that percentile is 0) AND barrier score ≥ median AND airway score (FOXJ1, CAPS, PIFO, TPPP3, KRT5) < 75th percentile. Airway filter relaxes only if the terminal would otherwise have <25 cells.")
    lines.append("- `cycle_cldn4_low`: cycle score (MKI67, TOP2A, STMN1, PCNA, HMGB2) ≥ 80th percentile AND CLDN4 ≤ median. Quantile relaxes to 70th once if n<25.")
    lines.append("- `at2_like` root: AT2 score (SFTPC, SFTPB, SFTPA1, SFTPA2, NAPSA, SLC34A2) ≥ 75th percentile AND CLDN4 ≤ median, only if ≥30 cells. This root orients diffusion pseudotime. It is not a normal-lung atlas root.")
    lines.append("- Primary kernel if that root exists **and** mean DPT of the barrier terminal exceeds mean DPT of AT2-like cells, with ≥95% finite DPT: **0.8 × PseudotimeKernel + 0.2 × ConnectivityKernel** (CellRank 2 default mixture). Otherwise primary kernel = **ConnectivityKernel** alone. Velocity is off (these GEO matrices have no spliced/unspliced counts).")
    lines.append("- Fate of a unit = mean `P(barrier_cldn4_high)` over its non-terminal cells. Units with <15 non-terminal cells are listed and not tested. Terminal cells are left out of the mean.")
    lines.append("")
    lines.append("## Honest n")
    lines.append("")
    if len(inventory):
        lines.append("| dataset | unit | cells QC / epithelial QC | cells used | units used |")
        lines.append("|---|---|---:|---:|---:|")
        for _, r in inventory.iterrows():
            lines.append(
                f"| {r['dataset']} | {r.get('unit','')} | {int(r['n_cells_qc'])} / {int(r['n_cells_epi_qc'])} | {int(r['n_cells_used'])} | {int(r['n_units'])} |"
            )
    lines.append("")
    n_units = int(included["unit_id"].nunique()) if len(included) else 0
    n_by = included.groupby("dataset").size().to_dict() if len(included) else {}
    lines.append(
        f"Units entering the fate test (non-terminal cells ≥15): **n={n_units}** "
        f"({', '.join(f'{ds} {n_by.get(ds, 0)}' for ds in COHORTS)})."
    )
    lines.append("Do not quote cell counts as n.")
    lines.append("")
    lines.append("## Primary result")
    lines.append("")
    if meta.get("k", 0) >= 2:
        lines.append(
            f"DerSimonian–Laird meta-analysis of cohort mean unit fate minus 0.5 "
            f"(k={meta['k']}: {', '.join(meta.get('datasets', []))}): "
            f"**Δ = {fmt(meta['mu'])}** (SE {fmt(meta['se'])}, descriptive z={fmt(meta.get('z'), 2)}, "
            f"p={fmt_p(meta.get('p'))}, I²={fmt(100*meta.get('I2', float('nan')), 1)}%)."
        )
        lines.append("Positive Δ means the unit-average fate sits above the two-sink coin-flip of 0.5.")
    elif meta.get("k", 0) == 1:
        lines.append(
            f"Only one cohort entered the meta-analysis ({', '.join(meta.get('datasets', []))}): "
            f"Δ = {fmt(meta['mu'])} (SE {fmt(meta['se'])})."
        )
    else:
        lines.append("No cohort produced an estimable two-sink fate. See reasons in the cohort table.")
    lines.append("")
    lines.append("| dataset | kernel | units | median fate | mean fate | fraction >0.5 | Wilcoxon p vs 0.5 | contact absorption (intermediate) | within-unit CellRank mean |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|---:|")
    if len(cohort):
        for _, r in cohort.iterrows():
            if not r.get("estimable", False):
                lines.append(f"| {r['dataset']} | — | — | — | — | — | — | — | {r.get('reason','not estimable')} |")
                continue
            lines.append(
                f"| {r['dataset']} | {r.get('primary_kernel')} | {int(r['n_units_included'])} | "
                f"{fmt(r['median_unit_fate'])} | {fmt(r['mean_unit_fate'])} | {fmt(r['frac_units_gt_0.5'], 2)} | "
                f"{fmt_p(r['wilcoxon_p_vs_0.5'])} | {fmt(r['contact_absorb_intermediate'])} | {fmt(r['mean_within_unit_cr'])} |"
            )
    lines.append("")
    lines.append("Wilcoxon p-values are descriptive (units inside one cohort). They are not a substitute for the meta-analytic Δ.")
    lines.append("")
    lines.append("## What the terminals actually are")
    lines.append("")
    lines.append("| dataset | n barrier / cycle / AT2 / intermediate | mean CLDN4 barrier vs cycle | mean barrier-score (CLDN4 held out) barrier vs rest | airway score ≥q75 inside barrier | DPT barrier vs AT2 | orientation ok |")
    lines.append("|---|---|---|---|---:|---|---|")
    if len(cohort):
        for _, r in cohort.iterrows():
            if not r.get("estimable", False):
                continue
            lines.append(
                f"| {r['dataset']} | {r.get('n_barrier')} / {r.get('n_cycle')} / {r.get('n_at2')} / {r.get('n_intermediate')} | "
                f"{fmt(r.get('mean_cldn4_barrier'))} vs {fmt(r.get('mean_cldn4_cycle'))} | "
                f"{fmt(r.get('mean_barrier_score_barrier'))} vs {fmt(r.get('mean_barrier_score_rest'))} | "
                f"{fmt(r.get('airway_frac_in_barrier'), 2)} | "
                f"{fmt(r.get('dpt_barrier'))} vs {fmt(r.get('dpt_at2'))} | {r.get('orientation_ok')} |"
            )
    lines.append("")
    lines.append("## Secondary descriptions (not the primary n)")
    lines.append("")
    lines.append("- **Connectivity-only CellRank** is always computed. When orientation fails it is the primary kernel; otherwise it is a sensitivity that does not use pseudotime. Unit means are in `results/tables/unit_fate.tsv` (`mean_fate_connectivity`).")
    lines.append("- **Within-unit CellRank** (ConnectivityKernel on that unit alone, states re-defined inside the unit, ≥80 epithelial cells) does not let Harmony stitch patients. Cohort means of those unit fates are in the table above. Missing units are too small or lack both sinks.")
    lines.append("- **Contact matrix:** row-normalized kNN weights aggregated by the four states, then absorption of the intermediate (and AT2-like) state into the barrier sink. This is an undirected snapshot Markov model. It is the contact-tracing estimator. Table: `results/tables/contact_transition.tsv`.")
    lines.append("- **Multinomial logistic:** one random graph neighbor per cell; class = neighbor state; predictors = AT2 score, barrier score with CLDN4 held out, and cycle score. CLDN4 itself is not a predictor. Coefficients are cell-level and dependent, so they are not a p-value claim. Table: `results/tables/multinomial_coefs.tsv`.")
    lines.append("- **GPCCA macrostates** on the connectivity kernel are fit without using CLDN4 to choose them. The macrostate with the highest membership-weighted CLDN4 is named afterwards. Its fate among non-members is a sensitivity, not a test against 0.5 (there are 3–4 macrostates). Details are in `results/summary.json`.")
    lines.append("")
    lines.append("## What this is not")
    lines.append("")
    lines.append("- Not evidence that CLDN4 causes T/NK exclusion. The locked exclusion result is a different analysis.")
    lines.append("- Not a timed trajectory and not RNA velocity. Public matrices here have no spliced/unspliced layer.")
    lines.append("- Not a statement that AT2 becomes tumor. The root is an AT2-like program inside the tumor epithelial gate, and only when those cells exist and sit upstream.")
    lines.append("- Not a cross-dataset Harmony object, and not a merge with mouse or with the private 8KL matrices.")
    lines.append("- A fate near 0.5 means the two pre-specified sinks split the non-terminal mass. That is a result, not a failed audit of the n=65 correlation.")
    lines.append("")
    lines.append("## Software")
    lines.append("")
    ver = ", ".join(f"{k} {v}" for k, v in versions.items())
    lines.append(f"Python packages: {ver}.")
    lines.append("")
    lines.append("## Reproduce")
    lines.append("")
    lines.append("```bash")
    lines.append("python3 methods/cellrank_concordant4_cldn4_fate/download.py")
    lines.append("python3 methods/cellrank_concordant4_cldn4_fate/analyze.py")
    lines.append("```")
    lines.append("")
    lines.append("Figures: `results/figures/fig_cellrank_fate.png`.")
    text = "\n".join(lines) + "\n"
    (HERE / "FINDING.md").write_text(text)
    (RESULTS / "FINDING.md").write_text(text)


def plot_figure(cells: pd.DataFrame, units: pd.DataFrame) -> None:
    nature_style()
    FIGS.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(183 * MM, 118 * MM))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.15, 1.05], hspace=0.55, wspace=0.38)
    included = units[units["included"]] if len(units) else units
    vmin, vmax = 0.0, 1.0
    for i, ds in enumerate(COHORTS):
        ax = fig.add_subplot(gs[0, i])
        sub = cells[cells["dataset"] == ds] if len(cells) else pd.DataFrame()
        ax.set_title(ds, loc="left", fontweight="medium")
        if len(sub) == 0 or "umap1" not in sub:
            ax.text(0.5, 0.5, "not estimable", ha="center", va="center", transform=ax.transAxes, fontsize=6)
            ax.set_axis_off()
            continue
        sca = ax.scatter(
            sub["umap1"], sub["umap2"], c=sub["fate_primary"], s=1.2, cmap="viridis",
            vmin=vmin, vmax=vmax, linewidths=0, rasterized=True,
        )
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_linewidth(0.4)
        if i == 0:
            ax.set_ylabel("UMAP2")
        ax.set_xlabel("UMAP1")
        if i == 3:
            cbar = fig.colorbar(sca, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_label("P(barrier CLDN4-high)", fontsize=6)
            cbar.ax.tick_params(labelsize=5, width=0.4, length=2)
            cbar.outline.set_linewidth(0.4)
    ax = fig.add_subplot(gs[1, :])
    rng = np.random.default_rng(SEED)
    if len(included):
        for i, ds in enumerate(COHORTS):
            y = included.loc[included["dataset"] == ds, "mean_fate_primary"].to_numpy()
            if len(y) == 0:
                continue
            x = i + rng.uniform(-0.12, 0.12, size=len(y))
            ax.scatter(x, y, s=10, color=DATASET_COLORS[ds], linewidths=0, zorder=3, label=f"{ds} (n={len(y)})")
            ax.hlines(np.median(y), i - 0.22, i + 0.22, colors="black", linewidths=0.8, zorder=4)
    ax.axhline(0.5, color="#666666", linewidth=0.6, linestyle="--", zorder=1)
    ax.set_xlim(-0.6, 3.6)
    ax.set_xticks(range(4))
    ax.set_xticklabels(COHORTS, rotation=15, ha="right")
    ax.set_ylabel("Unit mean fate of non-terminal cells")
    ax.set_ylim(-0.02, 1.02)
    despine(ax)
    ax.legend(frameon=False, loc="lower right", ncol=2)
    fig.savefig(FIGS / "fig_cellrank_fate.png", dpi=200)
    fig.savefig(FIGS / "fig_cellrank_fate.pdf")
    plt.close(fig)


def package_versions() -> dict:
    import importlib.metadata as md

    out = {}
    for name in ["cellrank", "scanpy", "anndata", "numpy", "pandas", "scipy", "scikit-learn", "harmonypy"]:
        try:
            out[name] = md.version(name)
        except md.PackageNotFoundError:
            out[name] = "missing"
    return out


def json_safe(obj):
    if isinstance(obj, dict):
        return {str(k): json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [json_safe(v) for v in obj]
    if isinstance(obj, (np.floating, float)):
        x = float(obj)
        return x if np.isfinite(x) else None
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    if obj is None or isinstance(obj, str):
        return obj
    return str(obj)


def run_smoke() -> None:
    rng = np.random.default_rng(SEED)
    n = 360
    trunk = np.column_stack([np.linspace(0, 1, 140), np.zeros(140)])
    arm0 = np.column_stack([np.linspace(1, 2, 110), np.linspace(0, 0.8, 110)])
    arm1 = np.column_stack([np.linspace(1, 2, 110), np.linspace(0, -0.8, 110)])
    xy = np.vstack([trunk, arm0, arm1])
    xy += rng.normal(0, 0.01, xy.shape)
    # expression that matches the state rules along the arms
    cldn4 = np.concatenate([
        np.linspace(0, 0.2, 140),
        np.linspace(0.2, 3.0, 110),
        np.linspace(0.2, 0.1, 110),
    ])
    at2 = np.concatenate([np.linspace(3, 0.2, 140), np.full(110, 0.1), np.full(110, 0.1)])
    cycle = np.concatenate([np.full(140, 0.1), np.full(110, 0.1), np.linspace(0.2, 3.0, 110)])
    barrier = cldn4.copy()
    X = np.column_stack([
        cldn4, at2, barrier, cycle,
        rng.random(n), rng.random(n),
    ]).astype(np.float32)
    # store as counts-like positive
    X = np.clip(X, 0, None)
    adata = ad.AnnData(X)
    adata.var_names = ["CLDN4", "SFTPC", "KRT8", "MKI67", "EPCAM", "CDH1"]
    adata.obs["dataset"] = "smoke"
    adata.obs["unit_id"] = np.where(np.arange(n) < n // 2, "U1", "U2")
    adata.obs["barcode"] = [f"c{i}" for i in range(n)]
    adata.obs_names = adata.obs["barcode"]
    # manually plant scores and bypass score_genes by writing obs after a fake log
    adata.obs["cldn4"] = cldn4
    adata.obs["score_at2"] = at2
    adata.obs["score_barrier"] = barrier
    adata.obs["score_cycle"] = cycle
    adata.obs["score_airway"] = 0.0
    state, info = assign_states(cldn4, at2, barrier, cycle, np.zeros(n))
    assert info["estimable"], info
    assert info["n_barrier"] >= MIN_TERMINAL
    assert info["n_cycle"] >= MIN_TERMINAL
    adata.obs["state"] = pd.Categorical(state)
    adata.obsm["X_pca"] = xy
    sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=10, random_state=SEED)
    T, labs = contact_matrix(state, adata.obsp["connectivities"])
    absorb = absorption_from_T(T, labs)
    assert "intermediate" in absorb or "at2_like" in absorb
    assert 0.0 <= list(absorb.values())[0] <= 1.0
    orient = {"orientation_ok": False, "reason": "smoke connectivity only"}
    kernels = run_kernels(adata, orient)
    non = ~np.isin(state, ["barrier_cldn4_high", "cycle_cldn4_low"])
    mu = float(np.nanmean(kernels["fate_connectivity"][non]))
    assert 0.0 <= mu <= 1.0, mu
    term = kernels["fate_connectivity"][state == "barrier_cldn4_high"]
    assert np.allclose(term, 1.0, atol=1e-5)
    print(f"SMOKE_OK absorb={absorb} mean_nonterminal_fate={mu:.3f} n_barrier={info['n_barrier']}", flush=True)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=Path, default=DEFAULT_DATA)
    p.add_argument("--datasets", nargs="*", default=COHORTS)
    p.add_argument("--cap", type=int, default=CAP_PER_UNIT)
    p.add_argument("--smoke", action="store_true")
    p.add_argument("--rebuild", action="store_true")
    args = p.parse_args()
    if args.smoke:
        run_smoke()
        return
    sc.settings.verbosity = 1
    sc.settings.n_jobs = 1
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    np.random.seed(SEED)
    inventories = []
    cell_frames = []
    summaries = []
    coef_frames = []
    trans_frames = []
    for ds in args.datasets:
        if ds not in LOADERS:
            raise SystemExit(f"refusing dataset outside concordant-4: {ds}")
        print(f"=== {ds} ===", flush=True)
        cache = CACHE / f"{ds}.h5ad"
        try:
            if cache.exists() and not args.rebuild:
                print(f"  cache {cache}", flush=True)
                adata = sc.read_h5ad(cache)
                info = {
                    "dataset": ds,
                    "unit": UNIT_KIND[ds],
                    "n_units": int(adata.obs["unit_id"].nunique()),
                    "n_cells_raw": int(adata.n_obs),
                    "n_cells_qc": int(adata.uns.get("inventory", {}).get("n_cells_qc", adata.n_obs)),
                    "n_cells_epi_qc": int(adata.uns.get("inventory", {}).get("n_cells_epi_qc", adata.n_obs)),
                    "n_cells_used": int(adata.n_obs),
                    "cap_per_unit": args.cap,
                    "note": "loaded from cache",
                }
                if "inventory" in adata.uns:
                    info.update(dict(adata.uns["inventory"]))
            else:
                adata, info = LOADERS[ds](args.data, args.cap, rng)
                adata.uns["inventory"] = {
                    k: (int(v) if isinstance(v, (np.integer, int)) else v) for k, v in info.items()
                }
                adata.write_h5ad(cache)
        except Exception as exc:
            print(f"  LOAD FAILED {ds}: {type(exc).__name__}: {exc}", flush=True)
            inventories.append(
                {
                    "dataset": ds,
                    "unit": UNIT_KIND[ds],
                    "n_units": 0,
                    "n_cells_raw": 0,
                    "n_cells_qc": 0,
                    "n_cells_epi_qc": 0,
                    "n_cells_used": 0,
                    "cap_per_unit": args.cap,
                    "note": f"load failed: {type(exc).__name__}: {exc}",
                }
            )
            summaries.append(
                {
                    "dataset": ds,
                    "estimable": False,
                    "reason": f"load failed: {type(exc).__name__}: {exc}",
                    "n_units": 0,
                    "n_cells": 0,
                }
            )
            continue
        inventories.append(info)
        try:
            cells, summary, coef, trans = run_dataset(adata, rng)
        except Exception as exc:
            print(f"  FAILED {ds}: {type(exc).__name__}: {exc}", flush=True)
            summary = {"dataset": ds, "estimable": False, "reason": f"{type(exc).__name__}: {exc}", "n_units": info["n_units"], "n_cells": info["n_cells_used"]}
            cells, coef, trans = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
        summaries.append(summary)
        if len(cells):
            cell_frames.append(cells)
        if len(coef):
            coef_frames.append(coef)
        if len(trans):
            trans_frames.append(trans)
        del adata
        gc.collect()
    inventory = pd.DataFrame(inventories)
    cells = pd.concat(cell_frames, ignore_index=True) if cell_frames else pd.DataFrame()
    units = unit_summary(cells) if len(cells) else pd.DataFrame()
    cohort, meta = cohort_rows(units, summaries)
    versions = package_versions()
    inventory.to_csv(TABLES / "inventory.tsv", sep="\t", index=False)
    if len(cells):
        cells.drop(columns=[c for c in cells.columns if c.startswith("within_status")], errors="ignore")
        # within status is not on cells; keep the fate columns. Drop bulky duplicate within fields? keep them.
        cells.to_csv(TABLES / "cell_fate.tsv.gz", sep="\t", index=False)
    if len(units):
        units.to_csv(TABLES / "unit_fate.tsv", sep="\t", index=False)
    cohort.to_csv(TABLES / "cohort_summary.tsv", sep="\t", index=False)
    if coef_frames:
        pd.concat(coef_frames, ignore_index=True).to_csv(TABLES / "multinomial_coefs.tsv", sep="\t", index=False)
    if trans_frames:
        pd.concat(trans_frames, ignore_index=True).to_csv(TABLES / "contact_transition.tsv", sep="\t", index=False)
    payload = {"versions": versions, "meta": meta, "summaries": summaries, "cap": args.cap, "seed": SEED}
    (RESULTS / "summary.json").write_text(json.dumps(json_safe(payload), indent=2))
    write_finding(units, cohort, meta, inventory, versions)
    if len(cells) and len(units):
        plot_figure(cells, units)
    print("DONE", flush=True)
    print(json.dumps(json_safe({"meta": meta, "cohort_estimable": cohort["estimable"].tolist() if len(cohort) and "estimable" in cohort else []}), indent=2))


if __name__ == "__main__":
    main()
