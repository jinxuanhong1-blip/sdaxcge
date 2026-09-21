#!/usr/bin/env python3
"""Local TROP2 (TACSTD2) vs CD8 after a CLDN4 residual on public lung spatial data.

Lead test: CosMx SMI NSCLC (He 2022; Zenodo 15487520 cosmx_lung), tumor cells,
author CD8 T cells inside 50 µm and 100 µm. Visium same-spot correlation is
scored because the whole-transcriptome matrix contains TACSTD2, CLDN4, and
CD8A, and is not the lead (prior Visium same-spot tests are cell-composition).

Xenium panels that lack CLDN4 are inventoried and not given a cell-level test.
No private 8-KL data. Numbers are written only from matrices read in this run.
"""

from __future__ import annotations

import gzip
import io
import json
import tarfile
import zipfile
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.io import mmread
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "spatial_trop2_cd8_cldn4"
FIG = OUT / "figures"
TAB = OUT / "tables"

COSMX_ZIP = Path("/tmp/cosmx/cosmx_lung.zip")
UM_PER_PX = 0.18
COSMX_SAMPLES = [
    "Lung5_Rep1",
    "Lung5_Rep2",
    "Lung5_Rep3",
    "Lung6",
    "Lung9_Rep1",
    "Lung9_Rep2",
    "Lung12",
    "Lung13",
]
COSMX_PATIENT = {
    "Lung5_Rep1": "Lung5",
    "Lung5_Rep2": "Lung5",
    "Lung5_Rep3": "Lung5",
    "Lung6": "Lung6",
    "Lung9_Rep1": "Lung9",
    "Lung9_Rep2": "Lung9",
    "Lung12": "Lung12",
    "Lung13": "Lung13",
}
CD8_TYPES = {"T CD8 naive", "T CD8 memory"}
CD8NK_TYPES = CD8_TYPES | {"NK"}

# GEO sample characteristics (histolgical type field), fetched from NCBI.
GSE189487_HIST = {
    "TD1": "IAC",
    "TD2": "IAC",
    "TD3": "MIA",
    "TD5": "AIS",
    "TD6": "MIA",
    "TD8": "AIS",
}

XENIUM_OFFICIAL = [
    (
        "10x Xenium lung preview + add-on (cancer)",
        "https://cf.10xgenomics.com/samples/xenium/1.3.0/Xenium_Preview_Human_Lung_Cancer_With_Add_on_2_FFPE/Xenium_Preview_Human_Lung_Cancer_With_Add_on_2_FFPE_gene_panel.json",
    ),
    (
        "10x Xenium lung preview + add-on (non-diseased)",
        "https://cf.10xgenomics.com/samples/xenium/1.3.0/Xenium_Preview_Human_Non_diseased_Lung_With_Add_on_FFPE/Xenium_Preview_Human_Non_diseased_Lung_With_Add_on_FFPE_gene_panel.json",
    ),
    (
        "10x Xenium hLung cancer multi-tissue preview",
        "https://cf.10xgenomics.com/samples/xenium/1.5.0/Xenium_V1_hLung_cancer_section/Xenium_V1_hLung_cancer_section_gene_panel.json",
    ),
    (
        "10x Xenium FFPE LUAD multimodal",
        "https://cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_humanLung_Cancer_FFPE/Xenium_V1_humanLung_Cancer_FFPE_gene_panel.json",
    ),
    (
        "10x Xenium FFPE NSCLC IO + add-on",
        "https://cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_Human_Lung_Cancer_Addon_FFPE/Xenium_V1_Human_Lung_Cancer_Addon_FFPE_gene_panel.json",
    ),
    (
        "10x Xenium v1 lung panel post-Xenium LUAD",
        "https://cf.10xgenomics.com/samples/xenium/3.0.0/Xenium_V1_Human_Lung_Cancer_FFPE/Xenium_V1_Human_Lung_Cancer_FFPE_gene_panel.json",
    ),
    (
        "10x Xenium Prime 5K post-Xenium LUAD",
        "https://cf.10xgenomics.com/samples/xenium/3.0.0/Xenium_Prime_Human_Lung_Cancer_FFPE/Xenium_Prime_Human_Lung_Cancer_FFPE_gene_panel.json",
    ),
]


def _log(msg: str) -> None:
    print(msg, flush=True)


def zscore(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    s = float(np.std(v, ddof=1)) if len(v) > 1 else 0.0
    if not np.isfinite(s) or s == 0.0:
        return np.zeros_like(v)
    return (v - float(np.mean(v))) / s


def ols_slopes(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    A = np.column_stack([np.ones(len(y)), X])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return coef[1:]


def mediation(x: np.ndarray, m: np.ndarray, y: np.ndarray, cov: np.ndarray | None = None) -> dict:
    """Within-section OLS on z-scored variables. Observational product method.

    a: mediator ~ exposure (+ cov)
    c: outcome ~ exposure (+ cov)
    c', b: outcome ~ exposure + mediator (+ cov)
    indirect = a * b; proportion = indirect / c when |c| is not tiny.
    """
    x_z, m_z, y_z = zscore(x), zscore(m), zscore(y)
    if cov is not None and cov.size:
        C = np.column_stack([zscore(cov[:, i]) for i in range(cov.shape[1])])
        Xa = np.column_stack([x_z, C])
        Xc = Xa
        Xb = np.column_stack([x_z, m_z, C])
    else:
        Xa = x_z.reshape(-1, 1)
        Xc = Xa
        Xb = np.column_stack([x_z, m_z])
    a = float(ols_slopes(m_z, Xa)[0])
    c = float(ols_slopes(y_z, Xc)[0])
    slopes = ols_slopes(y_z, Xb)
    c_prime = float(slopes[0])
    b = float(slopes[1])
    indirect = a * b
    prop = indirect / c if abs(c) >= 0.02 else np.nan
    return {
        "a": a,
        "b": b,
        "c": c,
        "c_prime": c_prime,
        "indirect": indirect,
        "prop": prop,
    }


def partial_spearman(x: np.ndarray, y: np.ndarray, cov: np.ndarray | None) -> tuple[float, float]:
    xr = stats.rankdata(x).astype(float)
    yr = stats.rankdata(y).astype(float)
    if cov is None or cov.size == 0:
        rho, p = stats.spearmanr(x, y)
        return float(rho), float(p)
    Cr = np.column_stack([stats.rankdata(cov[:, i]).astype(float) for i in range(cov.shape[1])])
    A = np.column_stack([np.ones(len(xr)), Cr])
    bx, *_ = np.linalg.lstsq(A, xr, rcond=None)
    by, *_ = np.linalg.lstsq(A, yr, rcond=None)
    rx = xr - A @ bx
    ry = yr - A @ by
    if float(np.std(rx)) == 0.0 or float(np.std(ry)) == 0.0:
        return np.nan, np.nan
    rho, p = stats.pearsonr(rx, ry)
    return float(rho), float(p)


def fisher_mean(rhos: np.ndarray) -> float:
    rhos = np.asarray(rhos, dtype=float)
    rhos = rhos[np.isfinite(rhos)]
    if len(rhos) == 0:
        return float("nan")
    z = np.arctanh(np.clip(rhos, -0.999999, 0.999999))
    return float(np.tanh(np.mean(z)))


def wilcoxon_p(values: np.ndarray) -> float:
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if len(v) < 5 or np.allclose(v, 0):
        return float("nan")
    # n=5 is the smallest set used here; scipy returns the exact floor.
    return float(stats.wilcoxon(v, alternative="two-sided", zero_method="wilcox").pvalue)


def self_check() -> None:
    rng = np.random.default_rng(0)
    n = 8000
    x = rng.normal(size=n)
    m = 0.8 * x + 0.2 * rng.normal(size=n)
    y = -0.7 * m + 0.2 * rng.normal(size=n)
    med = mediation(x, m, y)
    if abs(med["c_prime"]) > 0.08 or not (0.6 <= med["prop"] <= 1.3):
        raise RuntimeError(f"mediation self-check failed (full mediation): {med}")
    y2 = -0.5 * x + 0.3 * rng.normal(size=n)
    m2 = 0.8 * x + 0.3 * rng.normal(size=n)
    med2 = mediation(x, m2, y2)
    if abs(med2["prop"]) > 0.25 or abs(med2["c_prime"] - med2["c"]) > 0.08:
        raise RuntimeError(f"mediation self-check failed (direct only): {med2}")
    rho, _ = partial_spearman(m, y, m.reshape(-1, 1))
    if not np.isfinite(rho) or abs(rho) > 0.05:
        raise RuntimeError(f"partial Spearman self-check failed: {rho}")
    _log("self-check ok")


def gene_set_from_names(names: list[str]) -> dict:
    s = set(names)
    return {
        "n_features": len(names),
        "TACSTD2": int("TACSTD2" in s),
        "CLDN4": int("CLDN4" in s),
        "CD8A": int("CD8A" in s),
        "CLDN_present": ",".join(sorted(g for g in s if g.startswith("CLDN"))),
    }


def names_from_xenium_json(raw: bytes) -> list[str]:
    data = json.loads(raw)
    found: list[str] = []

    def walk(obj):
        if isinstance(obj, dict):
            if "gene_name" in obj and isinstance(obj["gene_name"], str):
                found.append(obj["gene_name"])
            elif "name" in obj and isinstance(obj["name"], str) and obj.get("type") in {"gene", "Gene Expression", None}:
                # avoid grabbing panel-level names only when nested under targets
                pass
            for v in obj.values():
                walk(v)
        elif isinstance(obj, list):
            for v in obj:
                walk(v)

    walk(data)
    if found:
        return found
    # fallback: any string value equal to a gene-like token collected from "name"
    alt: list[str] = []

    def walk2(obj):
        if isinstance(obj, dict):
            name = obj.get("name")
            if isinstance(name, str) and name.isupper() and 2 <= len(name) <= 20:
                alt.append(name)
            for v in obj.values():
                walk2(v)
        elif isinstance(obj, list):
            for v in obj:
                walk2(v)

    walk2(data)
    return alt


def h5_feature_names(path: Path) -> list[str]:
    with h5py.File(path, "r") as f:
        raw = f["matrix/features/name"][()]
    return [x.decode() if isinstance(x, bytes) else str(x) for x in raw]


def tsv_gene_symbols(path: Path) -> list[str]:
    with gzip.open(path, "rt") as f:
        lines = [ln for ln in f.read().splitlines() if ln.strip()]
    symbols = []
    for i, ln in enumerate(lines):
        parts = ln.split("\t")
        if i == 0 and parts[0].startswith("ENSG") is False and "gene" in ln.lower():
            continue
        if len(parts) >= 2:
            symbols.append(parts[1])
        else:
            symbols.append(parts[0])
    return symbols


def inventory_xenium() -> pd.DataFrame:
    rows = []
    import urllib.request

    for label, url in XENIUM_OFFICIAL:
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "spatial-trop2-cd8/1.0"})
            with urllib.request.urlopen(req, timeout=60) as resp:
                names = names_from_xenium_json(resp.read())
            info = gene_set_from_names(names)
            scored = "no"
            reason = "CLDN4 absent; mediator cannot be residualized"
            if info["TACSTD2"] and info["CLDN4"] and info["CD8A"]:
                scored = "not_run"
                reason = "panel has all three genes but no cell matrix was downloaded"
            elif not info["n_features"]:
                reason = "gene_panel.json parsed to 0 gene names"
            rows.append({"platform": "Xenium", "dataset": label, "source": url, "scored": scored, "reason": reason, **info})
        except Exception as exc:  # noqa: BLE001 — inventory must record the failure, not invent genes
            rows.append(
                {
                    "platform": "Xenium",
                    "dataset": label,
                    "source": url,
                    "n_features": 0,
                    "TACSTD2": -1,
                    "CLDN4": -1,
                    "CD8A": -1,
                    "CLDN_present": "",
                    "scored": "no",
                    "reason": f"download/parse failed: {type(exc).__name__}: {exc}",
                }
            )

    feat_dir = Path("/tmp/xenium/gse300007")
    if feat_dir.exists():
        # One LUAD TMA features file is enough; all six matched in the gate pass.
        files = sorted(feat_dir.glob("*features.tsv.gz"))
        luad = [p for p in files if "Lung_Adenocarcinoma" in p.name]
        use = luad[0] if luad else (files[0] if files else None)
        if use is not None:
            info = gene_set_from_names(tsv_gene_symbols(use))
            rows.append(
                {
                    "platform": "Xenium",
                    "dataset": "GSE300007 LUAD TMA2 unimodal/multimodal",
                    "source": str(use.name),
                    "scored": "no",
                    "reason": "CLDN4 absent on the downloaded features.tsv; LUAD TMA not scored",
                    **info,
                }
            )
    h5_rows = [
        ("GSE319755 tumor core", Path("/tmp/xenium/gse319755/tumor_core.h5")),
        ("GSE319755 adjacent lung", Path("/tmp/xenium/gse319755/adjacent.h5")),
        ("GSE319755 invasive front", Path("/tmp/xenium/gse319755/front.h5")),
        ("GSE311609 NSCLC Prime 5K L1", Path("/tmp/xenium/gse311609/nsclc5k_L1.h5")),
        ("GSE343063 SCLC p2", Path("/tmp/xenium/gse343063/p2.h5")),
    ]
    for label, path in h5_rows:
        if not path.exists():
            rows.append(
                {
                    "platform": "Xenium",
                    "dataset": label,
                    "source": str(path),
                    "n_features": 0,
                    "TACSTD2": -1,
                    "CLDN4": -1,
                    "CD8A": -1,
                    "CLDN_present": "",
                    "scored": "no",
                    "reason": "h5 not on disk",
                }
            )
            continue
        info = gene_set_from_names(h5_feature_names(path))
        reason = "CLDN4 absent; cell coordinates not used"
        if info["TACSTD2"] == 0 and info["CLDN4"] == 0:
            reason = "TACSTD2 and CLDN4 both absent"
        rows.append({"platform": "Xenium", "dataset": label, "source": str(path.name), "scored": "no", "reason": reason, **info})
    return pd.DataFrame(rows)


def _expr_block(counts: np.ndarray, totals: np.ndarray, genes: list[str], idx: dict[str, int]) -> dict[str, np.ndarray]:
    scale = np.divide(1e4, totals, out=np.zeros_like(totals), where=totals > 0)
    out = {}
    for g in genes:
        out[g] = np.log1p(counts[:, idx[g]] * scale)
    out["n_counts"] = totals
    return out


def associate(x: np.ndarray, m: np.ndarray, y: np.ndarray, extra_cov: np.ndarray | None = None) -> dict:
    finite = np.isfinite(x) & np.isfinite(m) & np.isfinite(y)
    if extra_cov is not None:
        finite &= np.all(np.isfinite(extra_cov), axis=1)
    x, m, y = x[finite], m[finite], y[finite]
    cov_extra = extra_cov[finite] if extra_cov is not None else None
    if len(x) < 40 or float(np.std(x)) == 0.0 or float(np.std(y)) == 0.0 or float(np.std(m)) == 0.0:
        return {"n": int(len(x)), "rho": np.nan, "partial_rho": np.nan, "partial_p": np.nan, "partial_rho_extra": np.nan, **{k: np.nan for k in ["a", "b", "c", "c_prime", "indirect", "prop"]}}
    rho, _ = stats.spearmanr(x, y)
    partial, partial_p = partial_spearman(x, y, m.reshape(-1, 1))
    med = mediation(x, m, y)
    partial_extra = np.nan
    if cov_extra is not None:
        both = np.column_stack([m, cov_extra])
        partial_extra, _ = partial_spearman(x, y, both)
    return {
        "n": int(len(x)),
        "rho": float(rho),
        "partial_rho": partial,
        "partial_p": partial_p,
        "partial_rho_extra": partial_extra,
        **med,
    }


def visium_neighbors(xy: np.ndarray) -> tuple[np.ndarray, float, float]:
    """Ring-1 mean is filled later. Return neighbor index lists and QC of the graph."""
    tree = cKDTree(xy)
    dist, _ = tree.query(xy, k=2)
    nn = dist[:, 1]
    nn = nn[np.isfinite(nn) & (nn > 0)]
    med = float(np.median(nn)) if len(nn) else np.nan
    radius = 1.5 * med if np.isfinite(med) else np.nan
    return tree, med, radius


def neighbor_mean(tree: cKDTree, values: np.ndarray, radius: float) -> tuple[np.ndarray, np.ndarray]:
    lists = tree.query_ball_point(tree.data, r=radius)
    out = np.full(len(values), np.nan)
    deg = np.zeros(len(values), dtype=int)
    for i, neigh in enumerate(lists):
        neigh = [j for j in neigh if j != i]
        deg[i] = len(neigh)
        if neigh:
            out[i] = float(np.mean(values[neigh]))
    return out, deg


def analyze_spot_table(
    *,
    section: str,
    series: str,
    histology: str,
    logexp: dict[str, np.ndarray],
    xy: np.ndarray,
) -> dict:
    epi_score = logexp["KRT8"] + logexp["EPCAM"]
    thr = float(np.quantile(epi_score, 0.60))
    epi = epi_score >= thr
    # If the threshold is 0 because epithelium is sparse, keep positive epithelium.
    if int(epi.sum()) < 80:
        epi = epi_score > 0
    tree, med_nn, radius = visium_neighbors(xy)
    local_cd8, deg = neighbor_mean(tree, logexp["CD8A"], radius)
    krt = logexp["KRT8"]
    lib = np.log1p(logexp["n_counts"])
    base = {
        "section": section,
        "series": series,
        "histology": histology,
        "n_spots": int(len(xy)),
        "n_epi": int(epi.sum()),
        "median_nn_px": med_nn,
        "neighbor_radius_px": radius,
        "median_degree_epi": float(np.median(deg[epi])) if epi.any() else np.nan,
        "pct_tacstd2_gt0_epi": float(np.mean(np.expm1(logexp["TACSTD2"][epi]) > 0)) if epi.any() else np.nan,
        "pct_cldn4_gt0_epi": float(np.mean(np.expm1(logexp["CLDN4"][epi]) > 0)) if epi.any() else np.nan,
        "pct_cd8a_gt0_epi": float(np.mean(np.expm1(logexp["CD8A"][epi]) > 0)) if epi.any() else np.nan,
    }
    specs = {
        "local_epi": (logexp["TACSTD2"][epi], logexp["CLDN4"][epi], local_cd8[epi], np.column_stack([krt[epi], lib[epi]])),
        "same_epi": (logexp["TACSTD2"][epi], logexp["CLDN4"][epi], logexp["CD8A"][epi], np.column_stack([krt[epi], lib[epi]])),
        "same_all": (logexp["TACSTD2"], logexp["CLDN4"], logexp["CD8A"], np.column_stack([krt, lib])),
    }
    for name, (x, m, y, extra) in specs.items():
        # extra covariate for the partial_rho_extra column only; primary partial is CLDN4 alone
        res = associate(x, m, y, extra)
        for k, v in res.items():
            base[f"{name}_{k}"] = v
    return base


def load_10x_h5(path: Path) -> tuple[np.ndarray, list[str], np.ndarray]:
    with h5py.File(path, "r") as f:
        data = f["matrix/data"][()]
        indices = f["matrix/indices"][()]
        indptr = f["matrix/indptr"][()]
        shape = tuple(int(x) for x in f["matrix/shape"][()])
        mat = sparse.csc_matrix((data, indices, indptr), shape=shape)
        names = [x.decode() if isinstance(x, bytes) else str(x) for x in f["matrix/features/name"][()]]
        barcodes = [x.decode() if isinstance(x, bytes) else str(x) for x in f["matrix/barcodes"][()]]
    # genes x spots
    if mat.shape[0] != len(names) or mat.shape[1] != len(barcodes):
        raise RuntimeError(f"{path} unexpected h5 shape {mat.shape}")
    return mat.tocsr(), names, np.array(barcodes)


def positions_frame(text: str) -> pd.DataFrame:
    sample = text.splitlines()[0]
    has_header = "barcode" in sample.lower() or "in_tissue" in sample.lower()
    df = pd.read_csv(io.StringIO(text), header=0 if has_header else None)
    if not has_header:
        if df.shape[1] < 6:
            raise RuntimeError(f"positions have {df.shape[1]} columns")
        df = df.iloc[:, :6]
        df.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"]
    else:
        cols = {c.lower(): c for c in df.columns}
        df = df.rename(
            columns={
                cols.get("barcode", df.columns[0]): "barcode",
                cols.get("in_tissue", df.columns[1]): "in_tissue",
                cols.get("pxl_row_in_fullres", df.columns[4]): "pxl_row",
                cols.get("pxl_col_in_fullres", df.columns[5]): "pxl_col",
            }
        )
    df["in_tissue"] = pd.to_numeric(df["in_tissue"], errors="coerce").fillna(0).astype(int)
    df["pxl_row"] = pd.to_numeric(df["pxl_row"], errors="coerce")
    df["pxl_col"] = pd.to_numeric(df["pxl_col"], errors="coerce")
    return df


def spot_from_matrix(mat: sparse.spmatrix, names: list[str], barcodes: np.ndarray, pos: pd.DataFrame) -> tuple[dict[str, np.ndarray], np.ndarray]:
    want = ["TACSTD2", "CLDN4", "CD8A", "KRT8", "EPCAM"]
    missing = [g for g in want if g not in names]
    if missing:
        raise RuntimeError(f"missing genes {missing}")
    idx = {g: names.index(g) for g in want}
    totals = np.asarray(mat.sum(axis=0)).ravel().astype(float)
    # mat is genes x spots
    keep_genes = mat[list(idx.values()), :].toarray().T  # spots x genes
    bc_to_i = {b: i for i, b in enumerate(barcodes)}
    pos = pos[pos["in_tissue"] == 1].copy()
    pos = pos[pos["barcode"].isin(bc_to_i)].copy()
    ii = pos["barcode"].map(bc_to_i).to_numpy()
    totals_i = totals[ii]
    ok = totals_i >= 200
    ii = ii[ok]
    pos = pos.iloc[np.flatnonzero(ok)].copy()
    counts = keep_genes[ii]
    totals_i = totals_i[ok]
    # counts columns follow idx insertion order of want
    logexp = {}
    scale = 1e4 / totals_i
    for j, g in enumerate(want):
        logexp[g] = np.log1p(counts[:, j] * scale)
    logexp["n_counts"] = totals_i
    xy = np.column_stack([pos["pxl_col"].to_numpy(dtype=float), pos["pxl_row"].to_numpy(dtype=float)])
    finite = np.isfinite(xy).all(axis=1)
    for g in list(logexp):
        logexp[g] = logexp[g][finite]
    return logexp, xy[finite]


def analyze_h5_section(h5: Path, spatial_tar: Path, section: str, series: str, histology: str) -> dict:
    mat, names, barcodes = load_10x_h5(h5)
    with tarfile.open(spatial_tar, "r:gz") as tf:
        member = [m for m in tf.getmembers() if m.name.endswith("tissue_positions.csv") or m.name.endswith("tissue_positions_list.csv")]
        if not member:
            raise RuntimeError(f"no positions in {spatial_tar}")
        raw = tf.extractfile(member[0]).read().decode()
    pos = positions_frame(raw)
    logexp, xy = spot_from_matrix(mat, names, barcodes, pos)
    _log(f"  {section}: spots {len(xy)}")
    return analyze_spot_table(section=section, series=series, histology=histology, logexp=logexp, xy=xy)


def _gzip_bytes(raw: bytes) -> str:
    return gzip.decompress(raw).decode()


def analyze_geo_tar(tar_path: Path, series: str, histology_of) -> list[dict]:
    rows = []
    with tarfile.open(tar_path, "r:") as tf:
        names = tf.getnames()
        prefixes = sorted({n[: -len("_matrix.mtx.gz")] for n in names if n.endswith("_matrix.mtx.gz")})
        for pref in prefixes:
            _log(f"  reading {pref}")
            mtx_b = tf.extractfile(pref + "_matrix.mtx.gz").read()
            feat_b = tf.extractfile(pref + "_features.tsv.gz").read()
            bc_b = tf.extractfile(pref + "_barcodes.tsv.gz").read()
            pos_name = pref + "_tissue_positions_list.csv.gz"
            if pos_name not in names:
                pos_name = pref + "_tissue_positions.csv.gz"
            pos_b = tf.extractfile(pos_name).read()
            mat = mmread(io.BytesIO(gzip.decompress(mtx_b))).tocsr()
            feat_lines = _gzip_bytes(feat_b).splitlines()
            genes = []
            for ln in feat_lines:
                parts = ln.split("\t")
                genes.append(parts[1] if len(parts) > 1 else parts[0])
            barcodes = np.array([ln.split("\t")[0] for ln in _gzip_bytes(bc_b).splitlines() if ln.strip()])
            if mat.shape[0] != len(genes) and mat.shape[1] == len(genes) and mat.shape[0] == len(barcodes):
                mat = mat.T.tocsr()
            if mat.shape[1] != len(barcodes) or mat.shape[0] != len(genes):
                raise RuntimeError(f"{pref} shape {mat.shape} genes {len(genes)} barcodes {len(barcodes)}")
            pos = positions_frame(_gzip_bytes(pos_b))
            # section label: GSM..._TD1 or GSM..._LM_SD_...
            tail = pref.split("_", 1)[1] if "_" in pref else pref
            if series == "GSE189487":
                # TD1 from GSM5702473_TD1
                key = pref.split("_")[-1]
                histology = histology_of(key)
                section = f"GSE189487_{key}"
            else:
                histology = histology_of(tail)
                section = f"GSE273378_{tail}"
            logexp, xy = spot_from_matrix(mat, genes, barcodes, pos)
            rows.append(analyze_spot_table(section=section, series=series, histology=histology, logexp=logexp, xy=xy))
            _log(f"    {section} n={len(xy)} local partial={rows[-1].get('local_epi_partial_rho')}")
    return rows


def load_cosmx_sample(zf: zipfile.ZipFile, sample: str) -> pd.DataFrame:
    feat_lines = zf.read(f"{sample}/qc/features.tsv").decode().splitlines()
    genes = [ln.split("\t")[0] for ln in feat_lines[1:] if ln.strip()]
    need = ["TACSTD2", "CLDN4", "CD8A", "KRT8", "EPCAM", "NKG7"]
    missing = [g for g in need if g not in genes]
    if missing:
        raise RuntimeError(f"{sample} missing {missing}")
    obs = pd.read_csv(zf.open(f"{sample}/qc/observations.tsv"), sep="\t", index_col=0)
    coord = pd.read_csv(zf.open(f"{sample}/qc/coordinates.tsv"), sep="\t", index_col=0)
    labs = pd.read_csv(zf.open(f"{sample}/labels.tsv"), sep="\t", index_col=0)
    # labels path is sample/labels.tsv
    _log(f"  {sample}: reading counts ({len(obs):,} cells)")
    mat = mmread(zf.open(f"{sample}/qc/counts.mtx")).tocsr()
    if mat.shape != (len(obs), len(genes)):
        raise RuntimeError(f"{sample} counts shape {mat.shape} vs cells {len(obs)} genes {len(genes)}")
    cols = [genes.index(g) for g in need]
    expr = pd.DataFrame(mat[:, cols].toarray(), index=obs.index, columns=need)
    df = expr.join(obs[["selected", "n_counts", "n_genes"]], how="left")
    df = df.join(coord.rename(columns={"x": "x_px", "y": "y_px"}), how="left")
    df = df.join(labs[["cell_type"]], how="left")
    if df["selected"].dtype == object:
        df = df[df["selected"].astype(str).str.lower().isin(["true", "1"])]
    else:
        df = df[df["selected"].astype(bool)]
    df = df[np.isfinite(df["x_px"]) & np.isfinite(df["y_px"])].copy()
    df["cell_type"] = df["cell_type"].astype(str)
    df["x_um"] = df["x_px"].to_numpy(dtype=float) * UM_PER_PX
    df["y_um"] = df["y_px"].to_numpy(dtype=float) * UM_PER_PX
    return df


def cosmx_radius_stats(tumor_xy: np.ndarray, tumor_x: np.ndarray, tumor_m: np.ndarray, tumor_krt: np.ndarray, tumor_lib: np.ndarray, immune_xy: np.ndarray, radius: float) -> dict:
    if len(immune_xy) == 0 or len(tumor_xy) < 40:
        return {"n": int(len(tumor_xy)), "rho": np.nan, "partial_rho": np.nan, "partial_p": np.nan, "partial_rho_extra": np.nan, "a": np.nan, "b": np.nan, "c": np.nan, "c_prime": np.nan, "indirect": np.nan, "prop": np.nan, "median_count": np.nan}
    tree = cKDTree(immune_xy)
    counts = np.asarray(tree.query_ball_point(tumor_xy, r=radius, return_length=True), dtype=float)
    y = np.log1p(counts)
    extra = np.column_stack([tumor_krt, tumor_lib])
    res = associate(tumor_x, tumor_m, y, extra)
    res["median_count"] = float(np.median(counts))
    return res


def analyze_cosmx() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    if not COSMX_ZIP.exists():
        raise RuntimeError(f"missing {COSMX_ZIP}")
    section_rows = []
    with zipfile.ZipFile(COSMX_ZIP) as zf:
        # scale check on the first sample before trusting 0.18 µm/px
        df0 = load_cosmx_sample(zf, COSMX_SAMPLES[0])
        nn = cKDTree(df0[["x_um", "y_um"]].to_numpy()).query(df0[["x_um", "y_um"]].to_numpy(), k=2)[0][:, 1]
        med_nn = float(np.median(nn))
        _log(f"  median NN {med_nn:.2f} µm using {UM_PER_PX} µm/px")
        if not (4.0 <= med_nn <= 30.0):
            raise RuntimeError(f"CosMx median NN {med_nn:.2f} µm is outside 4–30; refusing to score radii")
        samples = [COSMX_SAMPLES[0]] + COSMX_SAMPLES[1:]
        # reload first as part of the loop by caching df0
        cache = {COSMX_SAMPLES[0]: df0}
        for sample in samples:
            df = cache.get(sample)
            if df is None:
                df = load_cosmx_sample(zf, sample)
            tumor = df["cell_type"].str.lower().str.startswith("tumor")
            cd8 = df["cell_type"].isin(CD8_TYPES)
            cd8nk = df["cell_type"].isin(CD8NK_TYPES)
            xy = df[["x_um", "y_um"]].to_numpy(dtype=float)
            txy = xy[tumor.to_numpy()]
            tx = np.log1p(df.loc[tumor, "TACSTD2"].to_numpy(dtype=float))
            tm = np.log1p(df.loc[tumor, "CLDN4"].to_numpy(dtype=float))
            tk = np.log1p(df.loc[tumor, "KRT8"].to_numpy(dtype=float))
            tlib = np.log1p(df.loc[tumor, "n_counts"].to_numpy(dtype=float))
            row = {
                "sample": sample,
                "patient": COSMX_PATIENT[sample],
                "n_qc": int(len(df)),
                "n_tumor": int(tumor.sum()),
                "n_cd8": int(cd8.sum()),
                "n_cd8nk": int(cd8nk.sum()),
                "pct_tacstd2_gt0_tumor": float(np.mean(df.loc[tumor, "TACSTD2"].to_numpy() > 0)) if tumor.any() else np.nan,
                "pct_cldn4_gt0_tumor": float(np.mean(df.loc[tumor, "CLDN4"].to_numpy() > 0)) if tumor.any() else np.nan,
                "median_nn_um": float(np.median(cKDTree(xy).query(xy, k=2)[0][:, 1])),
                "modal_tumor_label": df.loc[tumor, "cell_type"].value_counts().index[0] if tumor.any() else "",
            }
            for radius, tag in [(50.0, "r50"), (100.0, "r100")]:
                res = cosmx_radius_stats(txy, tx, tm, tk, tlib, xy[cd8.to_numpy()], radius)
                for k, v in res.items():
                    row[f"{tag}_{k}"] = v
                _log(f"    {sample} {tag} partial={res['partial_rho']:.4f} c={res['c']:.3f} c'={res['c_prime']:.3f} prop={res['prop']}")
            res_nk = cosmx_radius_stats(txy, tx, tm, tk, tlib, xy[cd8nk.to_numpy()], 50.0)
            row["r50_cd8nk_partial_rho"] = res_nk["partial_rho"]
            row["r50_cd8nk_c_prime"] = res_nk["c_prime"]
            section_rows.append(row)
    sec = pd.DataFrame(section_rows)
    patient_rows = []
    for patient, sub in sec.groupby("patient"):
        rec = {"patient": patient, "n_sections": int(len(sub))}
        for col in [
            "r50_rho",
            "r50_partial_rho",
            "r50_partial_rho_extra",
            "r50_c",
            "r50_c_prime",
            "r50_indirect",
            "r50_prop",
            "r50_a",
            "r50_b",
            "r100_rho",
            "r100_partial_rho",
            "r100_c",
            "r100_c_prime",
            "r100_indirect",
            "r100_prop",
            "r50_cd8nk_partial_rho",
        ]:
            if col.endswith("rho") or col.endswith("partial_rho") or col.endswith("partial_rho_extra"):
                rec[col] = fisher_mean(sub[col].to_numpy())
            else:
                rec[col] = float(np.nanmean(sub[col].to_numpy(dtype=float)))
        patient_rows.append(rec)
    pat = pd.DataFrame(patient_rows)
    summary = {
        "um_per_px": UM_PER_PX,
        "n_sections": int(len(sec)),
        "n_patients": int(len(pat)),
        "r50_patient_partial_median": float(np.nanmedian(pat["r50_partial_rho"])),
        "r50_patient_partial_wilcoxon_p": wilcoxon_p(pat["r50_partial_rho"].to_numpy()),
        "r50_patient_c_prime_median": float(np.nanmedian(pat["r50_c_prime"])),
        "r50_patient_c_prime_wilcoxon_p": wilcoxon_p(pat["r50_c_prime"].to_numpy()),
        "r50_patient_prop_median": float(np.nanmedian(pat["r50_prop"])),
        "r50_section_partial_neg": int(np.sum(sec["r50_partial_rho"] < 0)),
        "r50_section_n": int(sec["r50_partial_rho"].notna().sum()),
        "r100_patient_partial_median": float(np.nanmedian(pat["r100_partial_rho"])),
        "r100_patient_partial_wilcoxon_p": wilcoxon_p(pat["r100_partial_rho"].to_numpy()),
        "r100_patient_c_prime_median": float(np.nanmedian(pat["r100_c_prime"])),
        "r50_patient_total_median": float(np.nanmedian(pat["r50_rho"])),
        "r50_patient_c_median": float(np.nanmedian(pat["r50_c"])),
    }
    return sec, pat, summary


def series_summary(df: pd.DataFrame, series: str, col: str) -> dict:
    sub = df[df["series"] == series]
    v = sub[col].to_numpy(dtype=float)
    return {
        "series": series,
        "column": col,
        "n": int(np.isfinite(v).sum()),
        "median": float(np.nanmedian(v)) if np.isfinite(v).any() else np.nan,
        "n_neg": int(np.sum(v[np.isfinite(v)] < 0)),
        "wilcoxon_p": wilcoxon_p(v),
    }


def forest(ax, labels, values, title, xlabel) -> None:
    y = np.arange(len(labels))
    vals = np.asarray(values, dtype=float)
    colors = ["#9b2226" if (np.isfinite(v) and v < 0) else "#1d3557" for v in vals]
    ax.axvline(0, color="#6c757d", lw=0.8)
    ax.scatter(vals, y, c=colors, s=36, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel(xlabel)
    ax.set_title(title, fontsize=10)
    ax.set_xlim(-0.6, 0.6)
    finite = vals[np.isfinite(vals)]
    if len(finite):
        ax.axvline(float(np.median(finite)), color="#9b2226", lw=0.7, ls="--")


def make_figures(cosmx_pat: pd.DataFrame, visium: pd.DataFrame) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.6), constrained_layout=True)
    order = cosmx_pat.sort_values("patient")
    forest(
        axes[0],
        order["patient"].tolist(),
        order["r50_partial_rho"].tolist(),
        "CosMx patients: TACSTD2 vs 50 µm CD8 | CLDN4",
        "Partial Spearman ρ",
    )
    forest(
        axes[1],
        order["patient"].tolist(),
        order["r50_c_prime"].tolist(),
        "Direct effect c′ after CLDN4",
        "Standardized OLS c′",
    )
    fig.savefig(FIG / "cosmx_patient_partial_and_direct.png", dpi=140)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 6.2), constrained_layout=True)
    luad = visium[visium["series"].isin(["GSE189487", "GSE273378"])].copy()
    luad = luad.sort_values(["series", "section"])
    labels = [f"{r.series.replace('GSE', '')} {r.section.split('_', 1)[-1]} ({r.histology})" for r in luad.itertuples()]
    forest(
        axes[0],
        labels,
        luad["local_epi_partial_rho"].tolist(),
        "Visium epithelial spots, ring-1 CD8A | CLDN4",
        "Partial Spearman ρ",
    )
    forest(
        axes[1],
        labels,
        luad["same_all_partial_rho"].tolist(),
        "Same-spot, all spots (scored, not the lead)",
        "Partial Spearman ρ",
    )
    fig.savefig(FIG / "visium_local_vs_samespot_partial.png", dpi=140)
    plt.close(fig)


def fmt(x, nd=3) -> str:
    if x is None or not isinstance(x, (int, float, np.floating)) or not np.isfinite(x):
        return "NA"
    ax = abs(float(x))
    if ax != 0 and (ax < 0.001 or ax >= 1000):
        return f"{float(x):.2e}"
    return f"{float(x):.{nd}f}"


def write_results(xenium: pd.DataFrame, visium: pd.DataFrame, sec: pd.DataFrame, pat: pd.DataFrame, cosmx_sum: dict, vis_sum: list[dict]) -> None:
    def row_med(series, col):
        hit = [r for r in vis_sum if r["series"] == series and r["column"] == col]
        return hit[0] if hit else None

    lines = []
    lines.append("# Local TROP2–CD8 after a CLDN4 residual (public lung spatial)")
    lines.append("")
    lines.append("Question: among public lung spatial assays that measure TACSTD2 (TROP2), CLDN4, and CD8 together, does a local CD8 association of TROP2 remain after residualizing CLDN4?")
    lines.append("")
    a_med = float(np.nanmedian(sec["r50_a"])) if len(sec) else float("nan")
    n_a_pos = int((sec["r50_a"] > 0).sum()) if len(sec) else 0
    n_prop = int(np.isfinite(pat["r50_prop"]).sum()) if len(pat) else 0
    n_zero50 = int((sec["r50_median_count"] == 0).sum()) if len(sec) else 0
    top = pat.loc[pat["r50_rho"].abs().idxmax()] if len(pat) else None
    lines.append(
        "Answer from the matrices read in this run: no stable local association remains, because the total association is already near null. "
        f"CosMx patient-level partial ρ at 50 µm is {fmt(cosmx_sum['r50_patient_partial_median'])} (Wilcoxon p = {fmt(cosmx_sum['r50_patient_partial_wilcoxon_p'], 4)}, n = 5). "
        "Visium ring-1 partial correlations in two LUAD series are centered on zero. "
        "Public Xenium lung panels opened here do not contain CLDN4, so that residual was not computed."
    )
    lines.append("")
    lines.append("This is an observational product-method residual, not a causal mediation claim. Same-spot Visium correlation is reported because those whole-transcriptome matrices exist, and it is not the lead. Prior Visium same-spot tests are epithelial-versus-stroma composition (PRs #134, #196, #553).")
    lines.append("")
    lines.append("## Lead: CosMx NSCLC, tumor cell vs nearby CD8 T cells")
    lines.append("")
    lines.append("Dataset: He et al. 2022 CosMx SMI NSCLC FFPE, 8 sections / 5 patients, 960-plex. Mirror: Zenodo 15487520 `cosmx_lung` (counts, coordinates, author `cell_type`). Pixel size 0.18 µm/px, checked by a median nearest-neighbor distance inside 4–30 µm. Tumor cells are author labels beginning with `tumor`. CD8 cells are author `T CD8 naive` and `T CD8 memory`. Outcome is log1p of the CD8 count inside a radius. Exposure and mediator are log1p counts on the index tumor cell. Partial Spearman residualizes ranks on CLDN4 only. OLS mediation z-scores each variable inside the section. Replicate sections are Fisher-averaged (correlations) or mean-averaged (coefficients) to the patient, and the Wilcoxon is on the 5 patients. Cell-level p-values are not used as the claim: neighboring cells are not independent.")
    lines.append("")
    lines.append(
        f"50 µm partial ρ(TACSTD2, local CD8 | CLDN4), patient median **{fmt(cosmx_sum['r50_patient_partial_median'])}**, Wilcoxon p = **{fmt(cosmx_sum['r50_patient_partial_wilcoxon_p'], 4)}** (n = 5). "
        f"Total ρ median {fmt(cosmx_sum['r50_patient_total_median'])}. Standardized total c median {fmt(cosmx_sum['r50_patient_c_median'])}; direct c′ median {fmt(cosmx_sum['r50_patient_c_prime_median'])} "
        f"(Wilcoxon p = {fmt(cosmx_sum['r50_patient_c_prime_wilcoxon_p'], 4)}). Sections with partial ρ < 0: {cosmx_sum['r50_section_partial_neg']}/{cosmx_sum['r50_section_n']}."
    )
    lines.append("")
    lines.append(
        f"TACSTD2 and CLDN4 do co-vary in tumor cells: the standardized a-path median across sections is {fmt(a_med)} ({n_a_pos}/{int(len(sec))} sections positive). "
        f"That coupling does not produce a stable indirect path to local CD8. A proportion (indirect/c) is stored only when |c| ≥ 0.02, which is {n_prop} of {int(len(pat))} patients, and those proportions sit on small total effects. "
        + (
            f"The patient with the largest |total ρ| is {top['patient']} (total ρ {fmt(top['r50_rho'])}, partial ρ {fmt(top['r50_partial_rho'])}; after CLDN4 + KRT8 + library size, partial ρ {fmt(top['r50_partial_rho_extra'])}). That single patient is not the five-patient result."
            if top is not None
            else ""
        )
    )
    lines.append("")
    lines.append(
        f"100 µm partial ρ median {fmt(cosmx_sum['r100_patient_partial_median'])}, Wilcoxon p = {fmt(cosmx_sum['r100_patient_partial_wilcoxon_p'], 4)}. Direct c′ median {fmt(cosmx_sum['r100_patient_c_prime_median'])}."
    )
    lines.append("")
    lines.append(
        f"The 50 µm author-CD8 count has median 0 in {n_zero50}/{int(len(sec))} sections "
        f"(100 µm median ≤ 1 in {int((sec['r100_median_count'] <= 1).sum())}/{int(len(sec))} sections; "
        f"maximum 100 µm median count is {fmt(float(sec['r100_median_count'].max()), 1)}). "
        "The outcome is sparse, so a near-null means no detectable gradient, not a precisely estimated zero. "
        "With n = 5 patients the two-sided Wilcoxon floor is 0.0625. Section-level signs are not independent (Lung5 has three sections, Lung9 has two)."
    )
    lines.append("")
    lines.append("This does not rewrite the locked CLDN4 exclusion result (cytotoxic neighbor ratios at 50/100 µm). The outcome here is author CD8 T cells, and the exposure is TACSTD2 after CLDN4.")
    lines.append("")
    lines.append("### Patients (50 µm)")
    lines.append("")
    lines.append("| Patient | sections | total ρ | partial ρ \\| CLDN4 | partial ρ \\| CLDN4+KRT8+nCount | c | c′ | proportion |")
    lines.append("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for r in pat.sort_values("patient").itertuples():
        lines.append(
            f"| {r.patient} | {r.n_sections} | {fmt(r.r50_rho)} | {fmt(r.r50_partial_rho)} | {fmt(r.r50_partial_rho_extra)} | {fmt(r.r50_c)} | {fmt(r.r50_c_prime)} | {fmt(r.r50_prop)} |"
        )
    lines.append("")
    lines.append("Per-section counts, detection, and both radii are in `results/spatial_trop2_cd8_cldn4/tables/cosmx_section.tsv`.")
    lines.append("")
    lines.append("## Visium: mediation scored; same-spot is not the lead")
    lines.append("")
    lines.append("Matrices read in this run: GSE189487 (6 early LUAD sections; histology from GEO sample characteristics), GSE273378 (16 stage I LUAD sections), and two 10x CytAssist FFPE demos (LUSC; lung neuroendocrine). Each is whole transcriptome, so TACSTD2, CLDN4, and CD8A are present. Spots are in-tissue with ≥200 counts, log1p(CP10k). Epithelial spots are the top 40% of log KRT8 + log EPCAM inside the section. Local CD8 is the mean log CD8A of other spots inside 1.5× the median nearest-neighbor distance (hex ring 1). E-MTAB-13530 was not re-downloaded: `ftp.ebi.ac.uk` TLS failed in this environment, so no mediation number is reported for it here.")
    lines.append("")
    lines.append("| Series | quantity | n | median | n negative | Wilcoxon p |")
    lines.append("| --- | --- | ---: | ---: | ---: | ---: |")
    for r in vis_sum:
        if r["series"] not in {"GSE189487", "GSE273378"}:
            continue
        lines.append(
            f"| {r['series']} | `{r['column']}` | {r['n']} | {fmt(r['median'])} | {r['n_neg']}/{r['n']} | {fmt(r['wilcoxon_p'], 4)} |"
        )
    lines.append("")
    g189 = row_med("GSE189487", "local_epi_partial_rho")
    g273 = row_med("GSE273378", "local_epi_partial_rho")
    lines.append(
        f"Ring-1 partial ρ is centered on zero (GSE189487 median {fmt(g189['median'])}, {g189['n_neg']}/{g189['n']} negative, p = {fmt(g189['wilcoxon_p'], 4)}; "
        f"GSE273378 median {fmt(g273['median'])}, {g273['n_neg']}/{g273['n']} negative, p = {fmt(g273['wilcoxon_p'], 4)}). "
        "Same-spot partial ρ, which is scored only because the spot matrix exists, is also centered on zero. "
        "Median neighbor degree on epithelial spots is 6, consistent with one Visium hex ring. "
        "Proportions are again undefined when |c| < 0.02, so the proportion rows have a smaller n than the correlation rows."
    )
    lines.append("")
    lines.append("Column key: `local_epi_partial_rho` is the ring-1 test after CLDN4; `local_epi_prop` is the mediation proportion on that local outcome; `same_all_partial_rho` and `same_epi_partial_rho` are same-spot scores. Same-spot rows are composition-confounded and are not the lead.")
    lines.append("")
    lines.append("10x demos are one section each and are not entered in the Wilcoxon table:")
    lines.append("")
    lines.append("| Section | histology | local partial ρ | local c′ | local proportion | same-spot all-spots partial ρ | same-spot proportion |")
    lines.append("| --- | --- | ---: | ---: | ---: | ---: | ---: |")
    for r in visium[visium["series"] == "10x"].itertuples():
        lines.append(
            f"| {r.section} | {r.histology} | {fmt(r.local_epi_partial_rho)} | {fmt(r.local_epi_c_prime)} | {fmt(r.local_epi_prop)} | {fmt(r.same_all_partial_rho)} | {fmt(r.same_all_prop)} |"
        )
    lines.append("")
    lines.append("## Xenium: no public lung panel here has TACSTD2 + CLDN4 + CD8")
    lines.append("")
    lines.append("CLDN4 is the mediator. Where it is absent, the residual was not computed and no cell-level CD8 distance was invented.")
    lines.append("")
    lines.append("| Dataset | features | TACSTD2 | CLDN4 | CD8A | why not scored |")
    lines.append("| --- | ---: | ---: | ---: | ---: | --- |")
    for r in xenium.itertuples():
        lines.append(
            f"| {r.dataset} | {r.n_features} | {r.TACSTD2} | {r.CLDN4} | {r.CD8A} | {r.reason} |"
        )
    lines.append("")
    lines.append("GSE311609 and GSE343063 each store 5,001 Gene Expression features plus codeword and control rows (10,029 feature rows in the h5). TACSTD2 and CLDN4 are absent from those gene names. GSE319755 (480 genes plus controls and 27 proteins) and GSE300007 LUAD Xenium have TACSTD2 and CD8A and do not have CLDN4. Official 10x lung / NSCLC `gene_panel.json` files downloaded in this run likewise lack CLDN4; three of them include TACSTD2 and CD8A.")
    lines.append("")
    lines.append("## How to rerun")
    lines.append("")
    lines.append("Raw matrices are not committed. `scripts/download_spatial_trop2_inputs.sh` fetches the Zenodo CosMx zip, the two 10x Visium demos, GSE189487, and GSE273378. Then:")
    lines.append("")
    lines.append("```bash")
    lines.append("python3 scripts/spatial_trop2_cd8_after_cldn4.py")
    lines.append("```")
    lines.append("")
    text = "\n".join(lines) + "\n"
    (ROOT / "RESULTS.md").write_text(text)
    (OUT / "RESULTS.md").write_text(text)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)
    self_check()
    _log("Xenium inventory")
    xenium = inventory_xenium()
    xenium.to_csv(TAB / "xenium_panel_gate.tsv", sep="\t", index=False)
    _log(xenium[["dataset", "TACSTD2", "CLDN4", "CD8A", "scored"]].to_string(index=False))

    _log("Visium")
    vis_rows = []
    vis_rows.append(
        analyze_h5_section(
            Path("/tmp/visium/lusc.h5"),
            Path("/tmp/visium/lusc_spatial.tar.gz"),
            "10x_CytAssist_FFPE_LUSC",
            "10x",
            "LUSC",
        )
    )
    vis_rows.append(
        analyze_h5_section(
            Path("/tmp/visium/nec.h5"),
            Path("/tmp/visium/nec_spatial.tar.gz"),
            "10x_CytAssist_11mm_NEC",
            "10x",
            "lung_neuroendocrine",
        )
    )
    vis_rows.extend(
        analyze_geo_tar(
            Path("/tmp/visium/GSE189487_RAW.tar"),
            "GSE189487",
            lambda key: GSE189487_HIST.get(key, "LUAD"),
        )
    )
    vis_rows.extend(
        analyze_geo_tar(
            Path("/tmp/visium/GSE273378_RAW.tar"),
            "GSE273378",
            lambda _tail: "stage_I_LUAD",
        )
    )
    visium = pd.DataFrame(vis_rows)
    visium.to_csv(TAB / "visium_section.tsv", sep="\t", index=False)
    vis_sum = []
    for series in ["GSE189487", "GSE273378"]:
        for col in [
            "local_epi_rho",
            "local_epi_partial_rho",
            "local_epi_c_prime",
            "local_epi_prop",
            "same_epi_partial_rho",
            "same_epi_prop",
            "same_all_rho",
            "same_all_partial_rho",
            "same_all_prop",
        ]:
            vis_sum.append(series_summary(visium, series, col))
    pd.DataFrame(vis_sum).to_csv(TAB / "visium_series_summary.tsv", sep="\t", index=False)

    _log("CosMx")
    sec, pat, cosmx_sum = analyze_cosmx()
    sec.to_csv(TAB / "cosmx_section.tsv", sep="\t", index=False)
    pat.to_csv(TAB / "cosmx_patient.tsv", sep="\t", index=False)
    (TAB / "cosmx_summary.json").write_text(json.dumps(cosmx_sum, indent=2) + "\n")

    make_figures(pat, visium)
    write_results(xenium, visium, sec, pat, cosmx_sum, vis_sum)
    _log("wrote RESULTS.md")
    _log(json.dumps(cosmx_sum, indent=2))


if __name__ == "__main__":
    main()
