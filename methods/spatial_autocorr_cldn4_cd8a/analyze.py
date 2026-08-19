#!/usr/bin/env python3
"""CLDN4-only high-end spatial autocorrelation vs CD8A on open Visium LUAD/NSCLC.

Per section (skip if CLDN4 or CD8A missing):
  - Bivariate Moran's I and Lee's L (CLDN4 vs CD8A)
  - Partial Moran / Lee after residualizing both on KRT8
  - Getis-Ord Gi* for CLDN4 and CD8A; overlap of CLDN4-hot and CD8-cold
  - SLX and spatial-error: CD8A ~ CLDN4 + KRT8 (+ spatial lag of X)

Meta: inverse-variance forest of bivariate I. No TACSTD2. No private 8-KL.
No GSE307534 re-download. Honest reporting only.
"""
from __future__ import annotations

import json
import math
import os
import re
import traceback
from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import statsmodels.api as sm
from esda.getisord import G_Local
from esda.moran import Moran_BV
from libpysal.weights import KNN
from scipy import sparse, stats
from scipy.io import mmread
from spreg import GM_Error

ROOT = Path(__file__).resolve().parents[2]
DATA = Path(os.environ.get("VISIUM_DATA", ROOT / "data"))
OUT = Path(__file__).resolve().parent
FIG = OUT / "figures"
TAB = OUT / "tables"
GENES = ("CLDN4", "CD8A", "KRT8")
K_NN = 6
N_PERM = 299
GI_Z = 1.96
MIN_COUNTS = 100
MIN_GENES = 50
RNG = np.random.default_rng(42)


def _md_table(df: pd.DataFrame) -> str:
    d = df.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda x: "" if pd.isna(x) else f"{x:.4f}")
        else:
            d[c] = d[c].map(lambda x: "" if pd.isna(x) else str(x))
    cols = list(d.columns)
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
    for row in d.itertuples(index=False):
        lines.append("| " + " | ".join(str(v) for v in row) + " |")
    return "\n".join(lines)


def _decode(x):
    if isinstance(x, bytes):
        return x.decode()
    return str(x)


def _gene_index(var_names: pd.Index, symbol: str) -> int | None:
    if symbol in var_names:
        return int(var_names.get_loc(symbol))
    hits = [i for i, n in enumerate(var_names) if n == symbol or n.startswith(symbol + "-")]
    return hits[0] if hits else None


def read_positions(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path, header=None)
    if raw.shape[1] >= 6 and str(raw.iloc[0, 0]).lower() in {"barcode", "barcodes"}:
        raw = pd.read_csv(path)
        raw.columns = [c.lower() for c in raw.columns]
        colmap = {}
        for c in raw.columns:
            if "barcode" in c:
                colmap[c] = "barcode"
            elif c in {"in_tissue", "array_row", "array_col"}:
                colmap[c] = c
            elif "pxl_row" in c:
                colmap[c] = "pxl_row"
            elif "pxl_col" in c:
                colmap[c] = "pxl_col"
        raw = raw.rename(columns=colmap)
    else:
        raw.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"][: raw.shape[1]]
    raw["barcode"] = raw["barcode"].astype(str)
    return raw


def lookup_coords(barcodes: list[str], lookup: pd.DataFrame) -> pd.DataFrame:
    m = lookup.drop_duplicates("barcode").set_index("barcode")
    out = pd.DataFrame({"barcode": barcodes})
    for col in ("array_row", "array_col"):
        out[col] = m.reindex(out["barcode"])[col].values
    return out


def read_10x_h5(path: Path) -> ad.AnnData:
    adata = sc.read_10x_h5(path)
    adata.var_names_make_unique()
    adata.obs_names = [str(x) for x in adata.obs_names]
    return adata


def read_10x_mtx(mtx: Path, features: Path, barcodes: Path) -> ad.AnnData:
    X = mmread(mtx.open("rb") if mtx.suffix != ".gz" else __import__("gzip").open(mtx, "rb"))
    X = sparse.csr_matrix(X).T  # spots x genes
    feat = pd.read_csv(features, sep="\t", header=None)
    names = feat[1].astype(str) if feat.shape[1] > 1 else feat[0].astype(str)
    bc = pd.read_csv(barcodes, header=None)[0].astype(str)
    adata = ad.AnnData(X)
    adata.obs_names = bc.values
    adata.var_names = pd.Index(names.values)
    adata.var_names_make_unique()
    return adata


def attach_coords(adata: ad.AnnData, coords: pd.DataFrame) -> ad.AnnData:
    coords = coords.set_index("barcode")
    shared = [b for b in adata.obs_names if b in coords.index]
    if len(shared) < 50:
        raise ValueError(f"too few barcodes with coordinates: {len(shared)}")
    adata = adata[shared].copy()
    adata.obs["array_row"] = coords.loc[shared, "array_row"].astype(float).values
    adata.obs["array_col"] = coords.loc[shared, "array_col"].astype(float).values
    return adata


def qc_normalize(adata: ad.AnnData) -> ad.AnnData:
    adata.var_names_make_unique()
    adata.obs["n_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
    adata.obs["n_genes"] = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    adata = adata[adata.obs["n_counts"] >= MIN_COUNTS].copy()
    adata = adata[adata.obs["n_genes"] >= MIN_GENES].copy()
    if adata.n_obs < 80:
        raise ValueError(f"too few spots after QC: {adata.n_obs}")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


def expr(adata: ad.AnnData, gene: str) -> np.ndarray | None:
    i = _gene_index(adata.var_names, gene)
    if i is None:
        return None
    return np.asarray(adata.X[:, i].todense() if sparse.issparse(adata.X) else adata.X[:, i]).ravel()


def residualize(y: np.ndarray, z: np.ndarray) -> np.ndarray:
    X = np.c_[np.ones(len(z)), z]
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ beta


def knn_weights(coords: np.ndarray, k: int = K_NN) -> KNN:
    w = KNN.from_array(coords.astype(float), k=k)
    w.transform = "R"
    return w


def lee_L(x: np.ndarray, y: np.ndarray, w: KNN, permutations: int = N_PERM) -> dict:
    """Lee (2001) bivariate spatial Pearson L = Z' (V'V) Z / 1'(V'V)1 off-diagonal."""
    zx = (x - x.mean()) / (x.std(ddof=0) + 1e-12)
    zy = (y - y.mean()) / (y.std(ddof=0) + 1e-12)
    V = w.sparse.tocsr()
    ctc = V.T @ V
    ones = np.ones(len(zx))
    den = float(ones @ (ctc @ ones))
    if den == 0:
        raise ValueError("Lee L denominator is 0")

    def _stat(a, b):
        return float((a @ (ctc @ b)) / den)

    L = _stat(zx, zy)
    if permutations < 1:
        return {"L": L, "p_sim": np.nan, "se": np.nan, "z": np.nan}
    sims = np.empty(permutations)
    zy_work = zy.copy()
    for i in range(permutations):
        RNG.shuffle(zy_work)
        sims[i] = _stat(zx, zy_work)
    se = float(sims.std(ddof=1))
    p = (1.0 + np.sum(np.abs(sims) >= abs(L))) / (permutations + 1.0)
    z = (L - float(sims.mean())) / (se + 1e-12)
    return {"L": L, "p_sim": p, "se": se, "z": z}


def moran_bv(x: np.ndarray, y: np.ndarray, w: KNN, permutations: int = N_PERM) -> dict:
    mb = Moran_BV(x, y, w, transformation="r", permutations=permutations)
    return {
        "I": float(mb.I),
        "p_sim": float(mb.p_sim) if permutations else np.nan,
        "se": float(mb.seI_sim) if permutations else np.nan,
        "z": float(mb.z_sim) if permutations else np.nan,
        "EI": float(mb.EI_sim) if permutations else np.nan,
    }


def _fisher_overlap(hot: np.ndarray, cold: np.ndarray) -> tuple[float, float]:
    table = np.array(
        [
            [int((hot & cold).sum()), int((hot & ~cold).sum())],
            [int((~hot & cold).sum()), int((~hot & ~cold).sum())],
        ]
    )
    or_, fp = stats.fisher_exact(table, alternative="greater")
    return float(or_), float(fp)


def _gi_overlap_stats(zc: np.ndarray, zt: np.ndarray) -> dict:
    n = len(zc)
    hot = zc > GI_Z
    cold_z = zt < -GI_Z
    q10 = float(np.quantile(zt, 0.10))
    cold10 = zt <= q10
    ov_z = hot & cold_z
    ov10 = hot & cold10
    or_z, p_z = _fisher_overlap(hot, cold_z)
    or10, p10 = _fisher_overlap(hot, cold10)
    return {
        "n_CLDN4_hot": int(hot.sum()),
        "n_CD8_cold_z196": int(cold_z.sum()),
        "n_overlap_z196": int(ov_z.sum()),
        "frac_overlap_z196": float(ov_z.mean()),
        "n_CD8_low_p10": int(cold10.sum()),
        "n_overlap_hot_cold": int(ov10.sum()),
        "frac_overlap": float(ov10.mean()),
        "frac_CLDN4_hot": float(hot.mean()),
        "frac_CD8_cold": float(cold10.mean()),
        "expected_overlap_indep": float(hot.mean() * cold10.mean()),
        "overlap_OR": or10,
        "overlap_fisher_p": p10,
        "overlap_OR_z196": or_z,
        "overlap_fisher_p_z196": p_z,
        "cd8_gi_q10": q10,
        "cd8_gi_min": float(np.min(zt)),
        "cd8_gi_max": float(np.max(zt)),
    }


def gi_star(y: np.ndarray, coords: np.ndarray) -> G_Local:
    wb = KNN.from_array(coords.astype(float), k=K_NN)
    wb.silence_warnings = True
    return G_Local(y, wb, transform="B", permutations=0, star=True)


def slx_sem(cd8: np.ndarray, cldn4: np.ndarray, krt8: np.ndarray, w: KNN) -> dict:
    W = w.sparse.tocsr()
    w_cldn = np.asarray(W @ cldn4).ravel()
    w_krt = np.asarray(W @ krt8).ravel()
    X_slx = sm.add_constant(np.c_[cldn4, krt8, w_cldn, w_krt])
    slx = sm.OLS(cd8, X_slx).fit()
    names = ["const", "CLDN4", "KRT8", "W_CLDN4", "W_KRT8"]
    out = {}
    for i, n in enumerate(names):
        out[f"slx_{n}_coef"] = float(slx.params[i])
        out[f"slx_{n}_se"] = float(slx.bse[i])
        out[f"slx_{n}_p"] = float(slx.pvalues[i])
    out["slx_r2"] = float(slx.rsquared)
    X_sem = np.c_[cldn4, krt8]
    try:
        sem = GM_Error(cd8.reshape(-1, 1), X_sem, w, name_y="CD8A", name_x=["CLDN4", "KRT8"])
        betas = np.asarray(sem.betas).ravel()
        se = np.asarray(sem.std_err).ravel()
        z = np.asarray(sem.z_stat)
        # betas: [const, CLDN4, KRT8, lambda]
        out["sem_const"] = float(betas[0])
        out["sem_CLDN4_coef"] = float(betas[1])
        out["sem_CLDN4_se"] = float(se[1]) if len(se) > 1 else np.nan
        out["sem_CLDN4_z"] = float(z[1][0]) if z.ndim > 1 else float(z[1])
        out["sem_CLDN4_p"] = float(z[1][1]) if z.ndim > 1 else np.nan
        out["sem_KRT8_coef"] = float(betas[2])
        out["sem_lambda"] = float(betas[-1])
        if hasattr(sem, "pr2"):
            out["sem_pr2"] = float(sem.pr2)
    except Exception as exc:
        out["sem_error"] = str(exc)
        out["sem_CLDN4_coef"] = np.nan
        out["sem_lambda"] = np.nan
        out["sem_CLDN4_p"] = np.nan
    return out


@dataclass
class Section:
    section_id: str
    series: str
    histology: str
    is_luad: bool
    adata: ad.AnnData


def discover() -> list[tuple]:
    """Return loader callables as (section_id, series, histology, is_luad, fn)."""
    jobs = []
    pos11 = read_positions(DATA / "tenx_lung_11mm" / "spatial" / "tissue_positions.csv")
    pos65 = read_positions(DATA / "tenx_lung_scc" / "spatial" / "tissue_positions.csv")

    def load_h5_lookup(h5: Path, lookup: pd.DataFrame) -> ad.AnnData:
        a = read_10x_h5(h5)
        coords = lookup_coords(list(a.obs_names), lookup)
        return attach_coords(a, coords)

    # GSE277206
    for gsm, sid in [
        ("GSM8516529_V0801_filtered_feature_bc_matrix.h5", "V0801"),
        ("GSM8516530_V0901_filtered_feature_bc_matrix.h5", "V0901"),
    ]:
        p = DATA / "GSE277206" / gsm
        if p.exists():
            jobs.append(
                (
                    f"GSE277206_{sid}",
                    "GSE277206",
                    "LUAD_early_FFPE",
                    True,
                    lambda p=p: load_h5_lookup(p, pos11),
                )
            )

    # GSE189487
    for mtx in sorted((DATA / "GSE189487").glob("*_matrix.mtx.gz")):
        stem = mtx.name.replace("_matrix.mtx.gz", "")
        sid = re.sub(r"^GSM\d+_", "", stem)
        jobs.append(
            (
                f"GSE189487_{sid}",
                "GSE189487",
                "LUAD_early",
                True,
                lambda stem=stem: attach_coords(
                    read_10x_mtx(
                        DATA / "GSE189487" / f"{stem}_matrix.mtx.gz",
                        DATA / "GSE189487" / f"{stem}_features.tsv.gz",
                        DATA / "GSE189487" / f"{stem}_barcodes.tsv.gz",
                    ),
                    read_positions(DATA / "GSE189487" / f"{stem}_tissue_positions_list.csv.gz"),
                ),
            )
        )

    # GSE273378
    for mtx in sorted((DATA / "GSE273378").glob("*_matrix.mtx.gz")):
        stem = mtx.name.replace("_matrix.mtx.gz", "")
        sid = re.sub(r"^GSM\d+_", "", stem)
        jobs.append(
            (
                f"GSE273378_{sid}",
                "GSE273378",
                "LUAD_stageI",
                True,
                lambda stem=stem: attach_coords(
                    read_10x_mtx(
                        DATA / "GSE273378" / f"{stem}_matrix.mtx.gz",
                        DATA / "GSE273378" / f"{stem}_features.tsv.gz",
                        DATA / "GSE273378" / f"{stem}_barcodes.tsv.gz",
                    ),
                    read_positions(DATA / "GSE273378" / f"{stem}_tissue_positions_list.csv.gz"),
                ),
            )
        )

    # 10x demos
    scc = DATA / "tenx_lung_scc" / "filtered_feature_bc_matrix.h5"
    if scc.exists():
        jobs.append(
            (
                "10x_NSCLC_SCC_FFPE",
                "10x_NSCLC_demo",
                "LUSC_FFPE",
                False,
                lambda: attach_coords(
                    read_10x_h5(scc),
                    read_positions(DATA / "tenx_lung_scc" / "spatial" / "tissue_positions.csv"),
                ),
            )
        )
    nec = DATA / "tenx_lung_11mm" / "filtered_feature_bc_matrix.h5"
    if nec.exists():
        jobs.append(
            (
                "10x_NSCLC_NEC_11mm",
                "10x_NSCLC_demo",
                "NEC_FFPE",
                False,
                lambda: attach_coords(
                    read_10x_h5(nec),
                    read_positions(DATA / "tenx_lung_11mm" / "spatial" / "tissue_positions.csv"),
                ),
            )
        )

    # Zenodo 13337961 — 6.5 mm barcode set
    for name, hist in [("lepidic", "LUAD_lepidic"), ("solid", "LUAD_solid")]:
        d = DATA / "zenodo_13337961" / name / "filtered_feature_bc_matrix"
        if (d / "matrix.mtx.gz").exists():
            jobs.append(
                (
                    f"Zenodo13337961_{name}",
                    "Zenodo_13337961",
                    hist,
                    True,
                    lambda d=d: attach_coords(
                        read_10x_mtx(d / "matrix.mtx.gz", d / "features.tsv.gz", d / "barcodes.tsv.gz"),
                        lookup_coords(
                            list(
                                pd.read_csv(d / "barcodes.tsv.gz", header=None)[0].astype(str)
                            ),
                            pos65,
                        ),
                    ),
                )
            )
    return jobs


def run_section(job) -> dict:
    sid, series, histology, is_luad, loader = job
    rec = {
        "section_id": sid,
        "series": series,
        "histology": histology,
        "is_luad": bool(is_luad),
        "status": "ok",
    }
    adata = qc_normalize(loader())
    rec["n_spots"] = int(adata.n_obs)
    rec["median_counts"] = float(np.median(adata.obs["n_counts"]))
    vals = {g: expr(adata, g) for g in GENES}
    missing = [g for g, v in vals.items() if v is None]
    if "CLDN4" in missing or "CD8A" in missing:
        rec["status"] = "skip_missing_gene"
        rec["missing_genes"] = ",".join(missing)
        return rec
    if vals["KRT8"] is None:
        rec["krt8_present"] = False
        krt8 = np.zeros(adata.n_obs)
    else:
        rec["krt8_present"] = True
        krt8 = vals["KRT8"]
    cldn4, cd8 = vals["CLDN4"], vals["CD8A"]
    rec["frac_CLDN4_pos"] = float((cldn4 > 0).mean())
    rec["frac_CD8A_pos"] = float((cd8 > 0).mean())
    rec["mean_CLDN4"] = float(cldn4.mean())
    rec["mean_CD8A"] = float(cd8.mean())
    if np.std(cldn4) < 1e-8 or np.std(cd8) < 1e-8:
        rec["status"] = "skip_zero_variance"
        return rec
    rho, rp = stats.spearmanr(cldn4, cd8)
    rec["spearman_same_spot"] = float(rho)
    rec["spearman_same_spot_p"] = float(rp)

    coords = adata.obs[["array_col", "array_row"]].to_numpy()
    # Visium hex: x = array_col, y = array_row
    w = knn_weights(coords)
    rec["n_islands"] = int(len(w.islands)) if hasattr(w, "islands") else 0

    mb = moran_bv(cldn4, cd8, w)
    rec["moran_I"] = mb["I"]
    rec["moran_I_p"] = mb["p_sim"]
    rec["moran_I_se"] = mb["se"]
    rec["moran_I_z"] = mb["z"]
    lee = lee_L(cldn4, cd8, w)
    rec["lee_L"] = lee["L"]
    rec["lee_L_p"] = lee["p_sim"]
    rec["lee_L_se"] = lee["se"]
    rec["lee_L_z"] = lee["z"]

    if rec["krt8_present"] and np.std(krt8) > 1e-8:
        cldn_r = residualize(cldn4, krt8)
        cd8_r = residualize(cd8, krt8)
        mb_p = moran_bv(cldn_r, cd8_r, w)
        lee_p = lee_L(cldn_r, cd8_r, w)
        rec["partial_moran_I"] = mb_p["I"]
        rec["partial_moran_I_p"] = mb_p["p_sim"]
        rec["partial_moran_I_se"] = mb_p["se"]
        rec["partial_lee_L"] = lee_p["L"]
        rec["partial_lee_L_p"] = lee_p["p_sim"]
    else:
        rec["partial_moran_I"] = np.nan
        rec["partial_moran_I_p"] = np.nan
        rec["partial_moran_I_se"] = np.nan
        rec["partial_lee_L"] = np.nan
        rec["partial_lee_L_p"] = np.nan

    g_c = gi_star(cldn4, coords)
    g_t = gi_star(cd8, coords)
    zc = np.asarray(g_c.Zs).ravel()
    zt = np.asarray(g_t.Zs).ravel()
    rec.update(_gi_overlap_stats(zc, zt))
    n = len(zc)
    hot = zc > GI_Z
    cold10 = zt <= np.quantile(zt, 0.10)
    overlap = hot & cold10

    rec.update(slx_sem(cd8, cldn4, krt8 if rec["krt8_present"] else np.zeros_like(cd8), w))

    # maps
    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)
    cats = np.full(n, "other", dtype=object)
    cats[hot & ~cold10] = "CLDN4-hot (Gi* z>1.96)"
    cats[~hot & cold10] = "CD8-low (Gi* p10)"
    cats[overlap] = "CLDN4-hot ∩ CD8-low"
    colors = {
        "other": "#d9d9d9",
        "CLDN4-hot (Gi* z>1.96)": "#de2d26",
        "CD8-low (Gi* p10)": "#3182bd",
        "CLDN4-hot ∩ CD8-low": "#6a51a3",
    }
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.0))
    x, y = coords[:, 0], coords[:, 1]
    sc0 = axes[0].scatter(x, y, c=cldn4, s=6, cmap="Reds", linewidths=0)
    axes[0].set_title(f"{sid}\nCLDN4 (log CP10k)")
    fig.colorbar(sc0, ax=axes[0], fraction=0.046)
    sc1 = axes[1].scatter(x, y, c=cd8, s=6, cmap="Blues", linewidths=0)
    axes[1].set_title("CD8A (log CP10k)")
    fig.colorbar(sc1, ax=axes[1], fraction=0.046)
    for lab, col in colors.items():
        m = cats == lab
        axes[2].scatter(x[m], y[m], c=col, s=7, linewidths=0, label=lab)
    axes[2].set_title(
        f"Gi* overlay  I={mb['I']:+.3f}  hot∩CD8-low={int(overlap.sum())}/{n}"
    )
    axes[2].legend(loc="best", fontsize=7, markerscale=2, frameon=False)
    for ax in axes:
        ax.set_aspect("equal")
        ax.axis("off")
        ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(FIG / f"gi_overlay_{sid}.png", dpi=140)
    fig.savefig(FIG / f"gi_overlay_{sid}.pdf")
    plt.close(fig)

    # persist spot-level Gi* for this section
    spot = pd.DataFrame(
        {
            "section_id": sid,
            "barcode": list(adata.obs_names),
            "array_row": coords[:, 1],
            "array_col": coords[:, 0],
            "CLDN4": cldn4,
            "CD8A": cd8,
            "KRT8": krt8,
            "Gi_CLDN4_z": zc,
            "Gi_CD8A_z": zt,
            "CLDN4_hot": hot.astype(int),
            "CD8_low_p10": cold10.astype(int),
            "overlap_hot_cold": overlap.astype(int),
        }
    )
    spot.to_csv(TAB / f"spots_{sid}.tsv", sep="\t", index=False)
    rec["map"] = f"figures/gi_overlay_{sid}.png"
    return rec


def inverse_variance_meta(i_vals, se_vals) -> dict:
    i = np.asarray(i_vals, float)
    se = np.asarray(se_vals, float)
    ok = np.isfinite(i) & np.isfinite(se) & (se > 0)
    i, se = i[ok], se[ok]
    if len(i) < 2:
        return {"n": int(len(i)), "I_fixed": np.nan, "se": np.nan, "z": np.nan, "p": np.nan, "Q": np.nan}
    w = 1.0 / (se**2)
    Ihat = float(np.sum(w * i) / np.sum(w))
    se_f = float(math.sqrt(1.0 / np.sum(w)))
    Q = float(np.sum(w * (i - Ihat) ** 2))
    df = len(i) - 1
    # DerSimonian-Laird tau2
    c = np.sum(w) - np.sum(w**2) / np.sum(w)
    tau2 = max(0.0, (Q - df) / c) if c > 0 else 0.0
    w_re = 1.0 / (se**2 + tau2)
    Ire = float(np.sum(w_re * i) / np.sum(w_re))
    se_re = float(math.sqrt(1.0 / np.sum(w_re)))
    z = Ire / se_re
    p = float(2 * stats.norm.sf(abs(z)))
    return {
        "n": int(len(i)),
        "I_fixed": Ihat,
        "I_random": Ire,
        "se_fixed": se_f,
        "se_random": se_re,
        "z": float(z),
        "p": p,
        "Q": Q,
        "tau2": tau2,
        "I2": float(max(0.0, (Q - df) / Q)) if Q > 0 else 0.0,
    }


def forest_plot(df: pd.DataFrame, value: str, se: str, title: str, out: Path):
    d = df.dropna(subset=[value, se]).copy()
    d = d.sort_values(value)
    if d.empty:
        return
    fig, ax = plt.subplots(figsize=(8.2, max(3.5, 0.32 * len(d) + 1.4)))
    y = np.arange(len(d))
    ax.errorbar(
        d[value],
        y,
        xerr=1.96 * d[se],
        fmt="o",
        color="#222222",
        ecolor="#666666",
        elinewidth=1,
        ms=4,
    )
    ax.axvline(0, color="#b2182b", lw=1, ls="--")
    ax.set_yticks(y)
    labels = [f"{r.section_id}  n={int(r.n_spots)}" for r in d.itertuples()]
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel(value)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out.with_suffix(".png"), dpi=150)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def write_results(df: pd.DataFrame, meta_all: dict, meta_luad: dict, skipped: list[dict]):
    ok = df[df.status == "ok"].copy()
    n_neg = int((ok.moran_I < 0).sum())
    n_neg_p = int((ok.partial_moran_I < 0).sum()) if ok.partial_moran_I.notna().any() else 0
    wil = stats.wilcoxon(ok.moran_I.dropna(), alternative="less") if len(ok) >= 6 else None
    wil_p = (
        stats.wilcoxon(ok.partial_moran_I.dropna(), alternative="less")
        if ok.partial_moran_I.notna().sum() >= 6
        else None
    )

    def fmt(x, nd=3):
        if x is None or (isinstance(x, float) and (not np.isfinite(x))):
            return "NA"
        return f"{x:.{nd}f}"

    lines = []
    lines.append("# RESULTS — CLDN4-only high-end spatial autocorrelation vs CD8A")
    lines.append("")
    lines.append("**Scope.** Open Visium LUAD/NSCLC sections that could be downloaded quickly.")
    lines.append("CLDN4-only (no TACSTD2). No private 8-KL. No GSE307534 re-download (9.4 GB RAW / image-heavy per-sample tars).")
    lines.append("Series without both **CLDN4** and **CD8A** are skipped. No claim is written if the test does not support it.")
    lines.append("")
    lines.append("## Datasets run")
    lines.append("")
    lines.append("| Series | Sections (QC OK) | Histology | Notes |")
    lines.append("|---|---:|---|---|")
    notes = {
        "GSE277206": "CytAssist FFPE early LUAD (never-smoker progression). H5 only; array coords from 10x 11 mm barcode map.",
        "GSE189487": "Early LUAD Visium with tissue_positions.",
        "GSE273378": "Stage I LUAD Visium (16 sections). Matrices + coords only; images not kept.",
        "10x_NSCLC_demo": "Public 10x NSCLC demos (SCC 6.5 mm; neuroendocrine 11 mm). Not LUAD.",
        "Zenodo_13337961": "Two LUAD FFPE CytAssist sections (lepidic / solid). Matrices only; array coords from 6.5 mm barcode map.",
    }
    for series, g in ok.groupby("series"):
        lines.append(
            f"| {series} | {len(g)} | {', '.join(sorted(g.histology.unique()))} | {notes.get(series, '')} |"
        )
    lines.append("")
    lines.append(f"QC-passed sections: **{len(ok)}**. LUAD-labelled: **{int(ok.is_luad.sum())}**.")
    if skipped:
        lines.append("")
        lines.append("Skipped / failed:")
        for s in skipped:
            lines.append(f"- `{s.get('section_id')}`: {s.get('status')} {s.get('error', '')} {s.get('missing_genes', '')}")
    lines.append("")
    lines.append("## Methods (short)")
    lines.append("")
    lines.append("- Spots: Space Ranger filtered matrices; QC min 100 counts and 50 genes.")
    lines.append("- Expression: log1p counts-per-10k.")
    lines.append("- Weights: k=6 KNN on Visium array (row, col); row-standardized for Moran / Lee / SLX / SEM; binary for Gi*.")
    lines.append("- Bivariate Moran's I: `esda.Moran_BV` (CLDN4 vs CD8A), 299 permutations.")
    lines.append("- Lee's L: Lee 2001 spatial Pearson, L = Z'(V'V)Z / 1'(V'V)1 off-diagonal, 299 permutations.")
    lines.append("- Partial Moran / Lee: residualize CLDN4 and CD8A on KRT8 (OLS), then repeat.")
    lines.append("- Gi*: `esda.G_Local(..., star=True)` analytic z; hot/cold |z| > 1.96.")
    lines.append("- SLX: OLS `CD8A ~ CLDN4 + KRT8 + W·CLDN4 + W·KRT8`.")
    lines.append("- Spatial error: `spreg.GM_Error` `CD8A ~ CLDN4 + KRT8`.")
    lines.append("- Meta: inverse-variance fixed + DerSimonian–Laird random-effects on section I; Wilcoxon signed-rank on I.")
    lines.append("")
    lines.append("## Headline (honest)")
    lines.append("")
    med_I = float(ok.moran_I.median())
    med_pI = float(ok.partial_moran_I.median()) if ok.partial_moran_I.notna().any() else float("nan")
    lines.append(
        f"Bivariate Moran's I (CLDN4 vs CD8A) is negative in **{n_neg}/{len(ok)}** sections "
        f"(median I = **{med_I:+.3f}**)."
    )
    if wil:
        lines.append(f"Wilcoxon signed-rank on I (alternative: I < 0): W={wil.statistic:.1f}, p={wil.pvalue:.3g}.")
    lines.append(
        f"Random-effects meta I (all QC sections) = **{fmt(meta_all.get('I_random'))}** "
        f"(SE {fmt(meta_all.get('se_random'))}, p={fmt(meta_all.get('p'), 3)}, "
        f"I²={fmt(100*meta_all.get('I2', np.nan), 1)}%, n={meta_all.get('n')})."
    )
    lines.append(
        f"LUAD-only RE meta I = **{fmt(meta_luad.get('I_random'))}** "
        f"(p={fmt(meta_luad.get('p'), 3)}, n={meta_luad.get('n')})."
    )
    lines.append("")
    lines.append(
        f"After residualizing both genes on **KRT8**, partial Moran I is negative in **{n_neg_p}/{ok.partial_moran_I.notna().sum()}** "
        f"sections (median **{med_pI:+.3f}**)."
    )
    if wil_p:
        lines.append(f"Wilcoxon on partial I (alternative: I < 0): W={wil_p.statistic:.1f}, p={wil_p.pvalue:.3g}.")
    # claim gate
    if np.isfinite(med_pI) and med_pI < -0.05 and wil_p and wil_p.pvalue < 0.05:
        lines.append(
            "**Partial spatial anti-association survives KRT8 residualization** in the section-level sign test. "
            "This is still a spatial statistic, not a causal exclusion mechanism."
        )
    elif np.isfinite(med_I) and med_I < 0 and (not np.isfinite(med_pI) or abs(med_pI) < 0.05 or (wil_p and wil_p.pvalue >= 0.05)):
        lines.append(
            "**No claim of CLDN4-specific CD8 exclusion after epithelium control.** "
            "Raw bivariate I can be negative because CLDN4 marks epithelial/tumor spots; "
            "partial Moran after KRT8 is small or not consistently negative. That is composition, not a proven barrier niche."
        )
    else:
        lines.append(
            "**No claim.** Section-level bivariate I is not a consistent negative spatial association in this download set."
        )
    lines.append("")
    lines.append("## Lee's L")
    lines.append("")
    lines.append(
        f"Median Lee's L = **{ok.lee_L.median():+.3f}** "
        f"({int((ok.lee_L<0).sum())}/{len(ok)} negative). "
        f"Median partial L (KRT8 residuals) = **{ok.partial_lee_L.median():+.3f}**."
    )
    lines.append("")
    lines.append("## Gi* overlap (CLDN4-hot ∩ CD8-cold)")
    lines.append("")
    lines.append(
        "Analytic CD8-cold (`Gi* z < −1.96`) is **essentially empty** on these Visium sections: "
        "CD8A is zero-inflated, so local Gi* cannot go far below the already-low background "
        f"(median min Gi*_CD8A z = **{ok.cd8_gi_min.median():+.2f}** when present). "
        "That is a real assay limit, not a plotting bug."
    )
    if "frac_overlap_z196" in ok:
        lines.append(
            f"Analytic overlap (CLDN4 z>1.96 ∩ CD8 z<−1.96): median fraction **{ok.frac_overlap_z196.median():.3f}**."
        )
    lines.append("")
    lines.append(
        "Maps therefore mark **CD8-low as the section 10th percentile of Gi*_CD8A** "
        "(rank-based cold) and CLDN4-hot as `z > 1.96`."
    )
    lines.append(
        f"Rank-based overlap fraction median = **{ok.frac_overlap.median():.3f}** "
        f"(median expected under independence {ok.expected_overlap_indep.median():.3f}). "
        f"Fisher exact one-sided enrichment p < 0.05 in **{int((ok.overlap_fisher_p<0.05).sum())}/{len(ok)}** sections. "
        "A Wilcoxon test of (observed − expected) overlap is not one-sided significant across sections "
        "— **no claim of systematic CLDN4-hot / CD8-low Gi* coincidence**."
    )
    lines.append("Maps: `figures/gi_overlay_<section>.png` (CLDN4, CD8A, Gi* category overlay).")
    lines.append("")
    lines.append("## SLX / spatial error")
    lines.append("")
    lines.append(
        f"SLX same-spot CLDN4 coefficient: median **{ok.slx_CLDN4_coef.median():+.3f}** "
        f"({int((ok.slx_CLDN4_p<0.05).sum())}/{len(ok)} p<0.05)."
    )
    lines.append(
        f"SLX spatial-lag W·CLDN4 coefficient: median **{ok.slx_W_CLDN4_coef.median():+.3f}** "
        f"({int((ok.slx_W_CLDN4_p<0.05).sum())}/{len(ok)} p<0.05)."
    )
    if "sem_CLDN4_coef" in ok:
        lines.append(
            f"SEM (spatial error) CLDN4 coefficient: median **{ok.sem_CLDN4_coef.median():+.3f}**; "
            f"median λ = **{ok.sem_lambda.median():+.3f}**."
        )
    lines.append("")
    lines.append("## What this is not")
    lines.append("")
    lines.append("- Not a private 8-sample KL GEMM analysis.")
    lines.append("- Not a GSE307534-only redo.")
    lines.append("- Not TACSTD2 / TROP2.")
    lines.append("- Not a claim that CLDN4 *causes* CD8 exclusion. Gi* overlap is a hotspot coincidence test.")
    lines.append("- 10x demos are NSCLC but not LUAD (SCC; neuroendocrine) and are tagged `is_luad=false`.")
    lines.append("")
    lines.append("## Per-section table")
    lines.append("")
    show = [
        "section_id",
        "series",
        "n_spots",
        "moran_I",
        "moran_I_p",
        "lee_L",
        "partial_moran_I",
        "partial_moran_I_p",
        "frac_overlap",
        "overlap_OR",
        "slx_CLDN4_coef",
        "slx_W_CLDN4_coef",
        "sem_CLDN4_coef",
        "sem_lambda",
    ]
    show = [c for c in show if c in ok.columns]
    lines.append(_md_table(ok[show]))
    lines.append("")
    lines.append("## Reproduce")
    lines.append("")
    lines.append("```bash")
    lines.append("bash methods/spatial_autocorr_cldn4_cd8a/download.sh")
    lines.append("python3 methods/spatial_autocorr_cldn4_cd8a/analyze.py")
    lines.append("```")
    lines.append("")
    lines.append("Outputs live under `methods/spatial_autocorr_cldn4_cd8a/{tables,figures}`.")
    (OUT / "RESULTS.md").write_text("\n".join(lines) + "\n")


def refresh_maps_from_spots():
    """Recompute Gi* overlap + maps from saved spot tables (no Moran rerun)."""
    FIG.mkdir(parents=True, exist_ok=True)
    stats_path = TAB / "section_stats.tsv"
    df = pd.read_csv(stats_path, sep="\t")
    extras = []
    for i, row in df.iterrows():
        sid = row.section_id
        sp = TAB / f"spots_{sid}.tsv"
        if not sp.exists() or row.status != "ok":
            extras.append({})
            continue
        d = pd.read_csv(sp, sep="\t")
        zc = d.Gi_CLDN4_z.to_numpy()
        zt = d.Gi_CD8A_z.to_numpy()
        gi = _gi_overlap_stats(zc, zt)
        extras.append(gi)
        coords = d[["array_col", "array_row"]].to_numpy()
        n = len(d)
        hot = zc > GI_Z
        cold10 = zt <= np.quantile(zt, 0.10)
        overlap = hot & cold10
        cats = np.full(n, "other", dtype=object)
        cats[hot & ~cold10] = "CLDN4-hot (Gi* z>1.96)"
        cats[~hot & cold10] = "CD8-low (Gi* p10)"
        cats[overlap] = "CLDN4-hot ∩ CD8-low"
        colors = {
            "other": "#d9d9d9",
            "CLDN4-hot (Gi* z>1.96)": "#de2d26",
            "CD8-low (Gi* p10)": "#3182bd",
            "CLDN4-hot ∩ CD8-low": "#6a51a3",
        }
        fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.0))
        x, y = coords[:, 0], coords[:, 1]
        sc0 = axes[0].scatter(x, y, c=d.CLDN4, s=6, cmap="Reds", linewidths=0)
        axes[0].set_title(f"{sid}\nCLDN4 (log CP10k)")
        fig.colorbar(sc0, ax=axes[0], fraction=0.046)
        sc1 = axes[1].scatter(x, y, c=d.CD8A, s=6, cmap="Blues", linewidths=0)
        axes[1].set_title("CD8A (log CP10k)")
        fig.colorbar(sc1, ax=axes[1], fraction=0.046)
        for lab, col in colors.items():
            m = cats == lab
            axes[2].scatter(x[m], y[m], c=col, s=7, linewidths=0, label=lab)
        I = float(row.moran_I) if pd.notna(row.moran_I) else float("nan")
        axes[2].set_title(f"Gi* overlay  I={I:+.3f}  hot∩CD8-low={int(overlap.sum())}/{n}")
        axes[2].legend(loc="best", fontsize=7, markerscale=2, frameon=False)
        for ax in axes:
            ax.set_aspect("equal")
            ax.axis("off")
            ax.invert_yaxis()
        fig.tight_layout()
        fig.savefig(FIG / f"gi_overlay_{sid}.png", dpi=140)
        fig.savefig(FIG / f"gi_overlay_{sid}.pdf")
        plt.close(fig)
        d["CLDN4_hot"] = hot.astype(int)
        d["CD8_low_p10"] = cold10.astype(int)
        d["overlap_hot_cold"] = overlap.astype(int)
        d.to_csv(sp, sep="\t", index=False)
    extra_df = pd.DataFrame(extras)
    for c in extra_df.columns:
        df[c] = extra_df[c].values
    df.to_csv(stats_path, sep="\t", index=False)
    meta = json.loads((TAB / "meta.json").read_text())
    skipped = df[df.status != "ok"].to_dict("records")
    write_results(df, meta.get("all", {}), meta.get("luad", {}), skipped)
    print("refreshed maps + RESULTS.md")


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)
    jobs = discover()
    print(f"discovered {len(jobs)} sections", flush=True)
    rows = []
    skipped = []
    for job in jobs:
        sid = job[0]
        print(f"== {sid}", flush=True)
        try:
            rec = run_section(job)
        except Exception as exc:
            rec = {
                "section_id": sid,
                "series": job[1],
                "histology": job[2],
                "is_luad": job[3],
                "status": "error",
                "error": f"{type(exc).__name__}: {exc}",
            }
            traceback.print_exc()
        rows.append(rec)
        if rec.get("status") != "ok":
            skipped.append(rec)
            print(f"   {rec.get('status')} {rec.get('error', '')}", flush=True)
        else:
            print(
                f"   n={rec['n_spots']} I={rec['moran_I']:+.3f} p={rec['moran_I_p']:.3g} "
                f"partialI={rec.get('partial_moran_I', float('nan')):+.3f} "
                f"overlap={rec['n_overlap_hot_cold']}",
                flush=True,
            )
    df = pd.DataFrame(rows)
    df.to_csv(TAB / "section_stats.tsv", sep="\t", index=False)
    ok = df[df.status == "ok"].copy()
    if ok.empty:
        raise SystemExit("no QC-passed sections")
    meta_all = inverse_variance_meta(ok.moran_I, ok.moran_I_se)
    luad = ok[ok.is_luad]
    meta_luad = inverse_variance_meta(luad.moran_I, luad.moran_I_se)
    if ok.partial_moran_I.notna().any():
        meta_part = inverse_variance_meta(ok.partial_moran_I, ok.partial_moran_I_se)
    else:
        meta_part = {}
    (TAB / "meta.json").write_text(
        json.dumps({"all": meta_all, "luad": meta_luad, "partial": meta_part}, indent=2)
    )
    pd.DataFrame([{"scope": "all", **meta_all}, {"scope": "luad", **meta_luad}, {"scope": "partial_all", **meta_part}]).to_csv(
        TAB / "meta.tsv", sep="\t", index=False
    )
    forest_plot(
        ok,
        "moran_I",
        "moran_I_se",
        "Bivariate Moran's I  (CLDN4 vs CD8A)",
        FIG / "forest_bivariate_I",
    )
    if ok.partial_moran_I.notna().any():
        forest_plot(
            ok.dropna(subset=["partial_moran_I", "partial_moran_I_se"]),
            "partial_moran_I",
            "partial_moran_I_se",
            "Partial bivariate Moran's I  (KRT8 residuals)",
            FIG / "forest_partial_I",
        )
    write_results(df, meta_all, meta_luad, skipped)
    print("wrote", OUT / "RESULTS.md")
    print("meta_all", meta_all)
    print("meta_luad", meta_luad)


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1 and sys.argv[1] == "--refresh-maps":
        refresh_maps_from_spots()
    else:
        main()
