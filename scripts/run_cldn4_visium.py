#!/usr/bin/env python3
"""CLDN4-only Visium LUAD histology analysis.

Per section: Spearman(CLDN4, CD8A); among epithelial spots, Q4 vs Q1 CLDN4
for nearest CD8A-high distance, neighbor CD8A, and KRT8-residualized neighbor CD8A.
Stratify by pathologist histology (lepidic vs acinar; STAS vs non-STAS) when labels exist.
"""
from __future__ import annotations

import gzip
import json
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats
from sklearn.linear_model import LinearRegression
from sklearn.neighbors import NearestNeighbors

warnings.filterwarnings("ignore", category=UserWarning)

ROOT = Path("/workspace")
DATA = ROOT / "data"
MAPS = ROOT / "maps"
RES = ROOT / "results"
MAPS.mkdir(exist_ok=True)
RES.mkdir(exist_ok=True)

EPI_GENES = ["EPCAM", "KRT8", "KRT19", "KRT7", "NAPSA", "SFTPB"]
LEPIDIC_SIG = ["SFTPB", "SFTPC", "NAPSA", "PGC", "SLC34A2", "LAMP3", "PEBP4", "SFTPD"]
INVASIVE_SIG = ["CEACAM5", "CEACAM6", "MUC21", "CRABP2", "AGR2", "HMGA1", "S100P"]
TUMOR_LABELS = {
    "Lepidic",
    "Acinar",
    "Papillary",
    "Solid",
    "Micropapillary",
    "Cribriform",
    "STAS",
}
LEPIDIC_LABELS = {"Lepidic"}
ACINAR_LABELS = {"Acinar"}
STAS_LABELS = {"STAS"}


def _open(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def gene_index(names: np.ndarray) -> dict[str, int]:
    idx = {}
    for i, n in enumerate(names):
        idx.setdefault(str(n), i)
    return idx


def get_vec(mat: sparse.spmatrix, names: np.ndarray, symbol: str) -> np.ndarray:
    idx = gene_index(names)
    if symbol in idx:
        return np.asarray(mat.getrow(idx[symbol]).todense()).ravel().astype(np.float64)
    # case-insensitive / alias
    up = {k.upper(): v for k, v in idx.items()}
    if symbol.upper() in up:
        return np.asarray(mat.getrow(up[symbol.upper()]).todense()).ravel().astype(np.float64)
    return np.zeros(mat.shape[1], dtype=np.float64)


def sig_score(mat, names, genes) -> np.ndarray:
    vecs = []
    for g in genes:
        v = get_vec(mat, names, g)
        if v.sum() > 0:
            vecs.append(v)
    if not vecs:
        return np.zeros(mat.shape[1])
    X = np.vstack(vecs)
    # z-score genes then mean
    mu = X.mean(axis=1, keepdims=True)
    sd = X.std(axis=1, keepdims=True)
    sd[sd == 0] = 1
    return ((X - mu) / sd).mean(axis=0)


def load_h5(path: Path):
    import h5py

    with h5py.File(path, "r") as f:
        grp = f["matrix"]
        data = grp["data"][:]
        indices = grp["indices"][:]
        indptr = grp["indptr"][:]
        shape = tuple(grp["shape"][:])
        mat = sparse.csc_matrix((data, indices, indptr), shape=shape)
        barcodes = [b.decode() if isinstance(b, bytes) else str(b) for b in grp["barcodes"][:]]
        feat = grp["features"]
        if "name" in feat:
            names = np.array(
                [x.decode() if isinstance(x, bytes) else str(x) for x in feat["name"][:]]
            )
        else:
            names = np.array(
                [x.decode() if isinstance(x, bytes) else str(x) for x in feat["id"][:]]
            )
    return mat.tocsr(), np.array(barcodes), names


def load_mtx(mtx, barcodes, features):
    mat = sparse.load_npz(mtx) if str(mtx).endswith(".npz") else None
    if mat is None:
        from scipy.io import mmread

        mat = mmread(str(mtx)).tocsr()
    with _open(barcodes) as fh:
        bcs = np.array([ln.strip() for ln in fh if ln.strip()])
    with _open(features) as fh:
        names = []
        for ln in fh:
            parts = ln.rstrip("\n").split("\t")
            names.append(parts[1] if len(parts) > 1 else parts[0])
    return mat.tocsr(), bcs, np.array(names)


def load_positions(path: Path) -> pd.DataFrame:
    if path.suffix == ".parquet" or path.name.endswith(".parquet"):
        df = pd.read_parquet(path)
        # Visium HD parquet
        colmap = {c.lower(): c for c in df.columns}
        bc = colmap.get("barcode") or list(df.columns)[0]
        out = pd.DataFrame({"barcode": df[bc].astype(str)})
        for src, dst in [
            ("in_tissue", "in_tissue"),
            ("array_row", "array_row"),
            ("array_col", "array_col"),
            ("pxl_row_in_fullres", "pxl_row"),
            ("pxl_col_in_fullres", "pxl_col"),
        ]:
            if src in colmap:
                out[dst] = df[colmap[src]].values
        return out
    # csv / csv.gz, with or without header
    raw = pd.read_csv(path, header=None)
    if raw.shape[1] >= 6 and str(raw.iloc[0, 0]).lower() in {"barcode", "barcodes"}:
        raw = pd.read_csv(path)
        cols = list(raw.columns)
        out = pd.DataFrame({"barcode": raw[cols[0]].astype(str)})
        # try named
        def pick(*cands):
            for c in cands:
                if c in raw.columns:
                    return raw[c].values
            return None

        out["in_tissue"] = pick("in_tissue") 
        out["array_row"] = pick("array_row")
        out["array_col"] = pick("array_col")
        out["pxl_row"] = pick("pxl_row_in_fullres", "pxl_row")
        out["pxl_col"] = pick("pxl_col_in_fullres", "pxl_col")
        if out["in_tissue"].isna().all():
            out["in_tissue"] = 1
        return out
    # no header: barcode, in_tissue, array_row, array_col, pxl_row, pxl_col
    out = pd.DataFrame(
        {
            "barcode": raw[0].astype(str),
            "in_tissue": raw[1].astype(int) if raw.shape[1] > 1 else 1,
            "array_row": raw[2] if raw.shape[1] > 2 else 0,
            "array_col": raw[3] if raw.shape[1] > 3 else 0,
            "pxl_row": raw[4] if raw.shape[1] > 4 else raw[2],
            "pxl_col": raw[5] if raw.shape[1] > 5 else raw[3],
        }
    )
    return out


def align_spots(mat, barcodes, pos: pd.DataFrame):
    pos = pos.copy()
    pos["barcode"] = pos["barcode"].astype(str)
    # strip possible suffixes
    bmap = {b: i for i, b in enumerate(barcodes)}
    keep_i = []
    keep_b = []
    for b in pos["barcode"]:
        if b in bmap:
            keep_i.append(bmap[b])
            keep_b.append(b)
    if not keep_i:
        # try without suffix
        bmap2 = {b.split("_")[0]: i for i, b in enumerate(barcodes)}
        for b in pos["barcode"]:
            key = b.split("_")[0]
            if key in bmap2:
                keep_i.append(bmap2[key])
                keep_b.append(b)
    pos = pos.set_index("barcode").loc[keep_b].reset_index()
    mat = mat[:, np.array(keep_i)]
    barcodes = np.array(keep_b)
    if "in_tissue" in pos.columns:
        mask = pos["in_tissue"].fillna(1).astype(float) > 0
        if mask.sum() >= 50:
            mat = mat[:, np.where(mask.values)[0]]
            barcodes = barcodes[mask.values]
            pos = pos.loc[mask.values].reset_index(drop=True)
    return mat, barcodes, pos


def lognorm(mat: sparse.spmatrix) -> sparse.csr_matrix:
    ncount = np.asarray(mat.sum(axis=0)).ravel()
    ncount[ncount == 0] = 1
    # CPM-like 1e4
    scale = 1e4 / ncount
    # scale columns of csr (genes x spots)
    mat = mat.tocsc(copy=True)
    mat.data = mat.data * np.repeat(scale, np.diff(mat.indptr))
    mat.data = np.log1p(mat.data)
    return mat.tocsr()


def ols_residual(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    ok = np.isfinite(y) & np.isfinite(x)
    resid = np.full_like(y, np.nan, dtype=np.float64)
    if ok.sum() < 10:
        return resid
    mdl = LinearRegression()
    mdl.fit(x[ok].reshape(-1, 1), y[ok])
    resid[ok] = y[ok] - mdl.predict(x[ok].reshape(-1, 1))
    return resid


def mw_effect(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 5 or len(b) < 5:
        return dict(n1=len(a), n2=len(b), median1=np.nan, median2=np.nan, delta=np.nan, p=np.nan, u=np.nan)
    try:
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    except ValueError:
        u, p = np.nan, np.nan
    return dict(
        n1=int(len(a)),
        n2=int(len(b)),
        median1=float(np.median(a)),
        median2=float(np.median(b)),
        delta=float(np.median(a) - np.median(b)),
        p=float(p) if np.isfinite(p) else np.nan,
        u=float(u) if np.isfinite(u) else np.nan,
    )


def spearman_safe(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 20:
        return np.nan, np.nan, int(m.sum())
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p), int(m.sum())


def analyze_section(
    name: str,
    dataset: str,
    mat,
    names,
    barcodes,
    pos: pd.DataFrame,
    pathology: pd.Series | None,
    section_histology: str,
    pitch_um: float,
    make_map: bool,
) -> list[dict]:
    nfeat, nspot = mat.shape
    ncount = np.asarray(mat.sum(axis=0)).ravel()
    ngene = np.asarray((mat > 0).sum(axis=0)).ravel()
    qc = (ncount >= 100) & (ngene >= 50)
    if qc.sum() < 80:
        qc = ncount >= 50
    mat_q = mat[:, np.where(qc)[0]]
    barcodes = barcodes[qc]
    pos = pos.loc[qc.values if hasattr(qc, "values") else qc].reset_index(drop=True)
    if pathology is not None:
        pathology = pd.Series(pathology.values, index=list(pathology.index)).reindex(list(barcodes))
    ln = lognorm(mat_q)
    cldn4 = get_vec(ln, names, "CLDN4")
    cd8a = get_vec(ln, names, "CD8A")
    krt8 = get_vec(ln, names, "KRT8")
    epcam = get_vec(ln, names, "EPCAM")
    krt19 = get_vec(ln, names, "KRT19")
    epi_score = (krt8 + epcam + krt19) / 3.0
    lep = sig_score(ln, names, LEPIDIC_SIG)
    inv = sig_score(ln, names, INVASIVE_SIG)

    # coordinates in microns
    if "pxl_col" in pos.columns and pos["pxl_col"].notna().all():
        xy = np.column_stack([pos["pxl_col"].astype(float), pos["pxl_row"].astype(float)])
        # convert pixels to um approximately via array pitch
        if "array_col" in pos.columns:
            ac = pos["array_col"].astype(float).values
            ar = pos["array_row"].astype(float).values
            dpx = np.nanmedian(np.abs(np.diff(np.sort(np.unique(ac)))))
            # use array coords * pitch (Visium hex: col*100, row*100*sqrt(3)/2 ~ 86.6)
            xy_um = np.column_stack([ac * pitch_um, ar * pitch_um * np.sqrt(3) / 2.0])
        else:
            xy_um = xy
    else:
        ac = pos["array_col"].astype(float).values
        ar = pos["array_row"].astype(float).values
        xy_um = np.column_stack([ac * pitch_um, ar * pitch_um * np.sqrt(3) / 2.0])

    # epithelial: pathology tumor labels, else top 40% epi score among QC
    if pathology is not None and pathology.notna().any():
        lab = pathology.fillna("").astype(str)
        epi_mask = lab.isin(TUMOR_LABELS).values
        if epi_mask.sum() < 40:
            epi_mask = epi_score >= np.quantile(epi_score, 0.60)
    else:
        epi_mask = epi_score >= np.quantile(epi_score, 0.60)
        lab = pd.Series([""] * len(cldn4), index=range(len(cldn4)))
        # computational histology
        hist_comp = np.where(lep - inv > 0.15, "Lepidic_like", np.where(inv - lep > 0.15, "Acinar_like", "Intermediate"))
        lab = pd.Series(hist_comp)

    # CD8A-high: top quartile among QC spots with CD8A>0, fallback top 15%
    pos_cd8 = cd8a[cd8a > 0]
    if len(pos_cd8) >= 20:
        thr = np.quantile(cd8a, 0.75)
        cd8_high = cd8a >= max(thr, np.quantile(pos_cd8, 0.50))
    else:
        cd8_high = cd8a >= np.quantile(cd8a, 0.90)
    if cd8_high.sum() < 10:
        cd8_high = cd8a >= np.quantile(cd8a, 0.85)

    # nearest CD8A-high distance for each spot
    dist_nn = np.full(len(cldn4), np.nan)
    if cd8_high.sum() >= 3:
        nn = NearestNeighbors(n_neighbors=1, algorithm="kd_tree")
        nn.fit(xy_um[cd8_high])
        d, _ = nn.kneighbors(xy_um)
        dist_nn = d.ravel()
        # zero out self (a CD8-high spot distance to itself)
        dist_nn[cd8_high] = 0.0

    # neighbor CD8A: mean of spots within 1.6 * pitch (first hex ring ~100um visium)
    neigh_cd8 = np.full(len(cldn4), np.nan)
    radius = 1.65 * pitch_um
    nnr = NearestNeighbors(radius=radius, algorithm="kd_tree")
    nnr.fit(xy_um)
    neighs = nnr.radius_neighbors(xy_um, return_distance=False)
    for i, ix in enumerate(neighs):
        ix = ix[ix != i]
        if len(ix):
            neigh_cd8[i] = float(cd8a[ix].mean())
        else:
            neigh_cd8[i] = 0.0

    krt8_resid_cldn4 = ols_residual(cldn4, krt8)
    neigh_cd8_krt8_resid = ols_residual(neigh_cd8, krt8)

    rows = []

    def quartiles(mask, subset_name):
        idx = np.where(mask)[0]
        if len(idx) < 20:
            return None
        vals = cldn4[idx]
        q1t, q3t = np.quantile(vals, [0.25, 0.75])
        q1 = idx[vals <= q1t]
        q4 = idx[vals >= q3t]
        if len(q1) < 8 or len(q4) < 8:
            return None
        r_all, p_all, n_all = spearman_safe(cldn4[idx], cd8a[idx])
        r_k, p_k, n_k = spearman_safe(krt8_resid_cldn4[idx], cd8a[idx])
        d_nn = mw_effect(dist_nn[q4], dist_nn[q1])
        d_nb = mw_effect(neigh_cd8[q4], neigh_cd8[q1])
        d_kr = mw_effect(neigh_cd8_krt8_resid[q4], neigh_cd8_krt8_resid[q1])
        d_krt8 = mw_effect(krt8[q4], krt8[q1])
        return dict(
            dataset=dataset,
            section=name,
            subset=subset_name,
            section_histology=section_histology,
            n_spots=int(len(cldn4)),
            n_subset=int(len(idx)),
            n_q1=int(len(q1)),
            n_q4=int(len(q4)),
            n_cd8_high=int(cd8_high.sum()),
            cldn4_cd8a_spearman_all=spearman_safe(cldn4, cd8a)[0],
            cldn4_cd8a_spearman_all_p=spearman_safe(cldn4, cd8a)[1],
            cldn4_cd8a_spearman_epi=r_all,
            cldn4_cd8a_spearman_epi_p=p_all,
            cldn4_krt8_residual_vs_cd8a_spearman=r_k,
            cldn4_krt8_residual_vs_cd8a_p=p_k,
            q4_minus_q1_nn_cd8high_um=d_nn["delta"],
            q4_vs_q1_nn_cd8high_p=d_nn["p"],
            q4_nn_median_um=d_nn["median1"],
            q1_nn_median_um=d_nn["median2"],
            q4_minus_q1_neighbor_cd8a=d_nb["delta"],
            q4_vs_q1_neighbor_cd8a_p=d_nb["p"],
            q4_neighbor_cd8a_median=d_nb["median1"],
            q1_neighbor_cd8a_median=d_nb["median2"],
            q4_minus_q1_krt8resid_neighbor_cd8a=d_kr["delta"],
            q4_vs_q1_krt8resid_neighbor_cd8a_p=d_kr["p"],
            q4_minus_q1_krt8=d_krt8["delta"],
            q4_vs_q1_krt8_p=d_krt8["p"],
            mean_cldn4_subset=float(np.mean(cldn4[idx])),
            mean_cd8a_subset=float(np.mean(cd8a[idx])),
            mean_krt8_subset=float(np.mean(krt8[idx])),
        )

    rows.append(quartiles(epi_mask, "epithelial"))
    # histology strata
    if pathology is not None and pathology.notna().any():
        labv = pathology.fillna("").astype(str).values
        rows.append(quartiles(epi_mask & np.isin(labv, list(LEPIDIC_LABELS)), "lepidic"))
        rows.append(quartiles(epi_mask & np.isin(labv, list(ACINAR_LABELS)), "acinar"))
        rows.append(quartiles(epi_mask & np.isin(labv, list(STAS_LABELS)), "STAS"))
        rows.append(quartiles(epi_mask & ~np.isin(labv, list(STAS_LABELS)), "nonSTAS_tumor"))
        for extra in ["Papillary", "Solid", "Micropapillary"]:
            rows.append(quartiles(epi_mask & (labv == extra), extra.lower()))
    else:
        labv = lab.values if hasattr(lab, "values") else np.array(lab)
        rows.append(quartiles(epi_mask & (labv == "Lepidic_like"), "lepidic_like"))
        rows.append(quartiles(epi_mask & (labv == "Acinar_like"), "acinar_like"))

    rows = [r for r in rows if r is not None]

    if make_map:
        fig, axes = plt.subplots(2, 2, figsize=(10, 9), dpi=140)
        x = xy_um[:, 0]
        y = -xy_um[:, 1]
        def sc(ax, val, title, cmap="viridis"):
            ax.scatter(x, y, c=val, s=4, cmap=cmap, linewidths=0)
            ax.set_title(title, fontsize=10)
            ax.set_aspect("equal")
            ax.axis("off")
        sc(axes[0, 0], cldn4, "CLDN4 (log1p)", "YlOrRd")
        sc(axes[0, 1], cd8a, "CD8A (log1p)", "YlGnBu")
        qcat = np.full(len(cldn4), 0.0)
        if epi_mask.sum() >= 20:
            ev = cldn4[epi_mask]
            q1t, q3t = np.quantile(ev, [0.25, 0.75])
            qcat[epi_mask] = 2
            qcat[epi_mask & (cldn4 <= q1t)] = 1
            qcat[epi_mask & (cldn4 >= q3t)] = 3
        sc(axes[1, 0], qcat, "Epithelial CLDN4 Q1(1)/mid(2)/Q4(3)", "coolwarm")
        if pathology is not None and pathology.notna().any():
            labv = pathology.fillna("Unlabeled").astype(str).values
            codes, uniques = pd.factorize(labv)
            axes[1, 1].scatter(x, y, c=codes, s=4, cmap="tab20", linewidths=0)
            axes[1, 1].set_title("Pathology labels", fontsize=10)
            axes[1, 1].set_aspect("equal")
            axes[1, 1].axis("off")
        else:
            sc(axes[1, 1], lep - inv, "Lepidic minus invasive signature", "PiYG")
        fig.suptitle(f"{dataset} | {name} | {section_histology}", fontsize=11)
        fig.tight_layout()
        fig.savefig(MAPS / f"{dataset}__{name}.png")
        plt.close(fig)

    return rows


def load_geo273378():
    d = DATA / "GSE273378"
    # parse SOFT-derived section histology from filenames + a sidecar
    meta = {
        "GSM8427428_LM_SD_1216_1": ("VI", True),
        "GSM8427429_LM_SD_16": ("VI", True),
        "GSM8427430_LM_SD_11": ("NST", False),
        "GSM8427431_LM_SD_2": ("NST", False),
        "GSM8427432_LM_SD_3": ("NST", False),
        "GSM8427433_LM_SD_4": ("VI", True),
        "GSM8427434_LM_SD_5": ("NST", False),
        "GSM8427435_LM_SD_6": ("VI", True),
        "GSM8427436_LM_SD_7": ("LMP", False),
        "GSM8427437_LM_SD_1216_8": ("VI", False),
        "GSM8427438_LM_SD_9": ("VI", True),
        "GSM8427439_LM_SD_10": ("NST", False),
        "GSM8427440_LM_SD_1216_12": ("NST", False),
        "GSM8427441_LM_SD_13": ("NST", False),
        "GSM8427442_LM_SD_1216_14": ("NST", False),
        "GSM8427443_LM_SD_15": ("NST", False),
    }
    out = []
    for mtx in sorted(d.glob("*_matrix.mtx.gz")):
        stem = mtx.name.replace("_matrix.mtx.gz", "")
        bc = d / f"{stem}_barcodes.tsv.gz"
        ft = d / f"{stem}_features.tsv.gz"
        posf = d / f"{stem}_tissue_positions_list.csv.gz"
        pathf = d / f"{stem}_pathology.csv.gz"
        mat, barcodes, names = load_mtx(mtx, bc, ft)
        pos = load_positions(posf)
        mat, barcodes, pos = align_spots(mat, barcodes, pos)
        path = pd.read_csv(pathf)
        pmap = dict(zip(path.iloc[:, 0].astype(str), path.iloc[:, 1].astype(str)))
        pathology = pd.Series([pmap.get(b, "") for b in barcodes], index=barcodes)
        grade, has_vi_frame = meta.get(stem, ("unk", False))
        labs = pathology.value_counts()
        has_lep = int(labs.get("Lepidic", 0))
        has_aci = int(labs.get("Acinar", 0))
        has_stas = int(labs.get("STAS", 0))
        hist = f"grade={grade};Lepidic={has_lep};Acinar={has_aci};STAS={has_stas}"
        out.append(
            dict(
                name=stem,
                dataset="GSE273378",
                mat=mat,
                names=names,
                barcodes=barcodes,
                pos=pos,
                pathology=pathology,
                section_histology=hist,
                pitch_um=100.0,
            )
        )
    return out


def load_geo300676():
    d = DATA / "GSE300676"
    out = []
    for h5 in sorted(d.glob("*_filtered_feature_bc_matrix.h5")):
        stem = h5.name.replace("_filtered_feature_bc_matrix.h5", "")
        spat = d / f"{stem}_spatial" / "spatial"
        posf = spat / "tissue_positions_list.csv"
        if not posf.exists():
            posf = spat / "tissue_positions.csv"
        mat, barcodes, names = load_h5(h5)
        pos = load_positions(posf)
        mat, barcodes, pos = align_spots(mat, barcodes, pos)
        case = stem.split("_")[-2] if "_" in stem else stem
        out.append(
            dict(
                name=stem,
                dataset="GSE300676",
                mat=mat,
                names=names,
                barcodes=barcodes,
                pos=pos,
                pathology=None,
                section_histology=f"mixed_lepidic_filigree_mPAP:{case}",
                pitch_um=100.0,
            )
        )
    return out


def load_geo189487():
    d = DATA / "GSE189487"
    hist = {
        "GSM5702473_TD1": "IAC",
        "GSM5702474_TD2": "IAC",
        "GSM5702475_TD3": "MIA",
        "GSM5702476_TD5": "AIS",
        "GSM5702477_TD6": "MIA",
        "GSM5702478_TD8": "AIS",
    }
    out = []
    for mtx in sorted(d.glob("*_matrix.mtx.gz")):
        stem = mtx.name.replace("_matrix.mtx.gz", "")
        mat, barcodes, names = load_mtx(
            mtx, d / f"{stem}_barcodes.tsv.gz", d / f"{stem}_features.tsv.gz"
        )
        pos = load_positions(d / f"{stem}_tissue_positions_list.csv.gz")
        mat, barcodes, pos = align_spots(mat, barcodes, pos)
        h = hist.get(stem, "unk")
        out.append(
            dict(
                name=stem,
                dataset="GSE189487",
                mat=mat,
                names=names,
                barcodes=barcodes,
                pos=pos,
                pathology=None,
                section_histology=h,
                pitch_um=100.0,
            )
        )
    return out


def load_kero():
    mapping = {
        "TSU-23": ("AIS_MIA_lepidic-like", 100.0),
        "TSU-36": ("AIS_MIA_lepidic-like", 100.0),
        "LUAD_No_17": ("IA_invasive", 100.0),
        "FFPE_LUAD_No3_A": ("IA_invasive_FFPE", 100.0),
    }
    out = []
    for folder, (hist, pitch) in mapping.items():
        p = DATA / "KERO" / folder
        h5 = p / "filtered_feature_bc_matrix.h5"
        posf = p / "spatial" / "tissue_positions_list.csv"
        if not posf.exists():
            posf = p / "spatial" / "tissue_positions.csv"
        mat, barcodes, names = load_h5(h5)
        pos = load_positions(posf)
        mat, barcodes, pos = align_spots(mat, barcodes, pos)
        out.append(
            dict(
                name=folder,
                dataset="KERO_AdSpatial2024",
                mat=mat,
                names=names,
                barcodes=barcodes,
                pos=pos,
                pathology=None,
                section_histology=hist,
                pitch_um=pitch,
            )
        )
    return out


def load_visium_hd():
    p = DATA / "VisiumHD_STAS" / "binned_outputs" / "square_016um"
    h5 = p / "filtered_feature_bc_matrix.h5"
    posf = p / "spatial" / "tissue_positions.parquet"
    mat, barcodes, names = load_h5(h5)
    pos = load_positions(posf)
    mat, barcodes, pos = align_spots(mat, barcodes, pos)
    return [
        dict(
            name="postXenium_Prime5K_Exp2_16um",
            dataset="VisiumHD_10x_LUAD_STAS",
            mat=mat,
            names=names,
            barcodes=barcodes,
            pos=pos,
            pathology=None,
            section_histology="LUAD_with_STAS_focus_G2G3",
            pitch_um=16.0,
        )
    ]


def main():
    loaders = [
        ("GSE273378", load_geo273378),
        ("GSE300676", load_geo300676),
        ("GSE189487", load_geo189487),
        ("KERO", load_kero),
        ("VisiumHD", load_visium_hd),
    ]
    all_rows = []
    map_priority = {
        "GSM8427428_LM_SD_1216_1",
        "GSM8427443_LM_SD_15",
        "GSM8427430_LM_SD_11",
        "GSM9066292_CRC_mPAP1_A",
        "GSM9066288_CRC_mPAP3_A",
        "GSM5702476_TD5",
        "GSM5702473_TD1",
        "TSU-23",
        "LUAD_No_17",
        "FFPE_LUAD_No3_A",
        "postXenium_Prime5K_Exp2_16um",
    }
    for label, fn in loaders:
        print(f"== loading {label}", flush=True)
        try:
            secs = fn()
        except Exception as e:
            print(f"FAILED load {label}: {e}", flush=True)
            continue
        print(f"  {len(secs)} sections", flush=True)
        for sec in secs:
            print(f"  analyze {sec['name']} spots={sec['mat'].shape[1]} genes={sec['mat'].shape[0]}", flush=True)
            try:
                rows = analyze_section(
                    name=sec["name"],
                    dataset=sec["dataset"],
                    mat=sec["mat"],
                    names=sec["names"],
                    barcodes=sec["barcodes"],
                    pos=sec["pos"],
                    pathology=sec["pathology"],
                    section_histology=sec["section_histology"],
                    pitch_um=sec["pitch_um"],
                    make_map=sec["name"] in map_priority or sec["dataset"] == "GSE273378",
                )
                all_rows.extend(rows)
                print(f"    subsets={ [r['subset'] for r in rows] }", flush=True)
            except Exception as e:
                print(f"  FAILED {sec['name']}: {e}", flush=True)
                import traceback

                traceback.print_exc()
    df = pd.DataFrame(all_rows)
    df.to_csv(RES / "section_metrics.csv", index=False)
    df.to_json(RES / "section_metrics.json", orient="records", indent=2)
    print(df.head(20).to_string())
    print("wrote", RES / "section_metrics.csv", "n=", len(df))


if __name__ == "__main__":
    main()
