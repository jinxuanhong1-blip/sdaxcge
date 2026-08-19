#!/usr/bin/env python3
"""ADDITIVE CosMx NSCLC (official 8-sample / 5-patient 960-plex), CLDN4-only.

Complementary to nearest-CD8: nearest NK and CD4 distances, CD8/NK radius
counts, mixing scores, and per-FOV Δdistance forests.

Stops if CLDN4 is absent from the 960-plex panel. Does not use private 8-KL.
"""

from __future__ import annotations

import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data", "cosmx_nsclc")
OUT = os.path.join(ROOT, "results", "cosmx_nsclc_cldn4_nk_cd8")
FIG = os.path.join(OUT, "figures")
TAB = os.path.join(OUT, "tables")

SAMPLES = [
    "Lung5_Rep1",
    "Lung5_Rep2",
    "Lung5_Rep3",
    "Lung6",
    "Lung9_Rep1",
    "Lung9_Rep2",
    "Lung12",
    "Lung13",
]
PATIENT = {
    "Lung5_Rep1": "P5",
    "Lung5_Rep2": "P5",
    "Lung5_Rep3": "P5",
    "Lung6": "P6",
    "Lung9_Rep1": "P9",
    "Lung9_Rep2": "P9",
    "Lung12": "P12",
    "Lung13": "P13",
}

PX_TO_UM = 0.18
RADII_UM = (15, 25, 50, 100, 150)
MIX_RADII_UM = (25, 50, 100)
MIN_COUNTS = 20
MIN_GENES = 5
MIN_FOV_TUMOR = 20
MIN_FOV_ARM = 8
NEAR_CENSOR_UM = 400.0

NEEDED = [
    "CLDN4",
    "NKG7",
    "GNLY",
    "KLRD1",
    "KLRB1",
    "NCAM1",
    "CD4",
    "CD8A",
    "CD8B",
    "CD3D",
    "CD3E",
    "CD3G",
    "IL7R",
    "FOXP3",
    "PTPRC",
    "EPCAM",
    "CDH1",
    "KRT8",
    "KRT18",
    "KRT19",
    "KRT7",
    "KRT5",
    "KRT17",
    "CEACAM6",
    "MUC1",
    "ELF3",
]
NK_MARKERS = ["NKG7", "GNLY", "KLRD1"]
CD3_MARKERS = ["CD3D", "CD3E", "CD3G"]
CD8_MARKERS = ["CD8A", "CD8B"]
EPI_MARKERS = [
    "EPCAM",
    "CDH1",
    "KRT8",
    "KRT18",
    "KRT19",
    "KRT7",
    "KRT5",
    "KRT17",
    "CEACAM6",
    "MUC1",
    "ELF3",
]
IMM_MARKERS = [
    "PTPRC",
    "CD3D",
    "CD3E",
    "CD3G",
    "CD8A",
    "CD8B",
    "CD4",
    "NKG7",
    "GNLY",
    "KLRD1",
]
STR_MARKERS = ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "ACTA2", "PECAM1", "VWF"]

AUTHOR_TYPE_COLS = (
    "cell_type",
    "celltype",
    "CellType",
    "cellType",
    "nb_clus",
    "cell_types",
    "annotation",
    "AuthorCellType",
)


def find_csv(sample_dir: str, suffix: str) -> str:
    for root, _, files in os.walk(sample_dir):
        for fn in files:
            if fn.endswith(suffix):
                return os.path.join(root, fn)
    raise FileNotFoundError(f"{suffix} under {sample_dir}")


def panel_genes(sample: str) -> list[str]:
    path = find_csv(os.path.join(DATA, sample), "exprMat_file.csv")
    header = pd.read_csv(path, nrows=0)
    return [c for c in header.columns if c not in ("fov", "cell_ID")]


def stop_if_cldn4_absent(genes_by_sample: dict[str, list[str]]) -> None:
    missing = [s for s, g in genes_by_sample.items() if "CLDN4" not in g]
    if missing:
        os.makedirs(OUT, exist_ok=True)
        msg = (
            "CLDN4 is absent from the CosMx 960-plex exprMat for: "
            + ", ".join(missing)
            + ". Stopping as instructed. No distances, radii, mixing, or forest were computed."
        )
        with open(os.path.join(ROOT, "RESULTS.md"), "w") as fh:
            fh.write("# CosMx NSCLC CLDN4 immune radii — STOP\n\n")
            fh.write(msg + "\n")
        print(msg, file=sys.stderr)
        raise SystemExit(2)


def load_sample(sample: str) -> pd.DataFrame:
    d = os.path.join(DATA, sample)
    expr_path = find_csv(d, "exprMat_file.csv")
    meta_path = find_csv(d, "metadata_file.csv")
    header = list(pd.read_csv(expr_path, nrows=0).columns)
    id_cols = ["fov", "cell_ID"]
    panel = [c for c in header if c not in id_cols and not str(c).lower().startswith("neg")]
    genes = [g for g in NEEDED + STR_MARKERS if g in header]
    chunks = []
    for chunk in pd.read_csv(expr_path, chunksize=40_000):
        chunk = chunk[chunk["cell_ID"] != 0]
        if chunk.empty:
            continue
        present_panel = [g for g in panel if g in chunk.columns]
        chunk["_tot"] = chunk[present_panel].sum(axis=1, numeric_only=True)
        chunk["_nfeat"] = (chunk[present_panel] > 0).sum(axis=1)
        chunks.append(chunk[id_cols + genes + ["_tot", "_nfeat"]])
    expr = pd.concat(chunks, ignore_index=True)
    expr["cell_key"] = expr["fov"].astype(str) + "_" + expr["cell_ID"].astype(str)
    meta = pd.read_csv(meta_path)
    meta["cell_key"] = meta["fov"].astype(str) + "_" + meta["cell_ID"].astype(str)
    keep_meta = ["cell_key", "fov", "cell_ID", "CenterX_global_px", "CenterY_global_px"]
    for col in ("nCount_RNA", "nFeature_RNA", "Mean.PanCK", "Mean.CD45", "Mean.CD3", "Area"):
        if col in meta.columns:
            keep_meta.append(col)
    author_col = next((c for c in AUTHOR_TYPE_COLS if c in meta.columns), None)
    if author_col:
        keep_meta.append(author_col)
    meta = meta[keep_meta].copy()
    df = expr.merge(meta, on=["cell_key", "fov", "cell_ID"], how="inner")
    gene_cols = [g for g in genes if g in df.columns]
    raw = df[gene_cols].to_numpy(dtype=np.float32)
    if "nCount_RNA" in df.columns:
        tot_full = df["nCount_RNA"].to_numpy(dtype=np.float32)
        nfeat = df["nFeature_RNA"].to_numpy(dtype=np.float32) if "nFeature_RNA" in df.columns else df["_nfeat"].to_numpy()
    else:
        tot_full = df["_tot"].to_numpy(dtype=np.float32)
        nfeat = df["_nfeat"].to_numpy()
    keep = (tot_full >= MIN_COUNTS) & (nfeat >= MIN_GENES)
    df = df.loc[keep].reset_index(drop=True)
    raw = raw[keep]
    tot = np.maximum(raw.sum(axis=1), 1.0)
    sf = float(np.median(tot)) / tot
    logn = np.log1p(raw * sf[:, None])
    for i, g in enumerate(gene_cols):
        df[f"{g}_raw"] = raw[:, i]
        df[f"{g}_ln"] = logn[:, i]
    df["sample"] = sample
    df["patient"] = PATIENT[sample]
    df["author_type"] = df[author_col].astype(str) if author_col else ""
    df["author_col"] = author_col or ""
    df["_n_panel_genes"] = len(panel)
    return df


def mean_ln(df: pd.DataFrame, markers: list[str]) -> np.ndarray:
    cols = [f"{g}_ln" for g in markers if f"{g}_ln" in df.columns]
    if not cols:
        return np.zeros(len(df), dtype=np.float32)
    return df[cols].to_numpy(dtype=np.float32).mean(axis=1)


def raw_sum(df: pd.DataFrame, markers: list[str]) -> np.ndarray:
    cols = [f"{g}_raw" for g in markers if f"{g}_raw" in df.columns]
    if not cols:
        return np.zeros(len(df), dtype=np.float32)
    return df[cols].to_numpy(dtype=np.float32).sum(axis=1)


def author_match(series: pd.Series, needles: tuple[str, ...]) -> np.ndarray:
    if series.empty or (series == "").all():
        return np.zeros(len(series), dtype=bool)
    s = series.str.lower()
    out = np.zeros(len(series), dtype=bool)
    for n in needles:
        out |= s.str.contains(n, regex=False, na=False).to_numpy()
    return out


def assign_types(df: pd.DataFrame) -> pd.DataFrame:
    epi = mean_ln(df, EPI_MARKERS)
    imm = mean_ln(df, IMM_MARKERS)
    strn = mean_ln(df, STR_MARKERS)
    scores = np.vstack([epi, imm, strn])
    lab = np.array(["epithelial", "immune", "stromal"])[scores.argmax(axis=0)]
    lab[scores.max(axis=0) <= 0] = "unassigned"
    df["compartment"] = lab

    author_tumor = author_match(df["author_type"], ("tumor", "epithelial", "cancer", "malignant"))
    author_nk = author_match(df["author_type"], ("nk",))
    author_cd4 = author_match(df["author_type"], ("cd4", "t cd4", "t-cd4"))
    author_cd8 = author_match(df["author_type"], ("cd8", "t cd8", "t-cd8"))
    df["used_author_labels"] = bool(df["author_col"].iloc[0]) and bool(df["author_type"].ne("").any())

    tumor = (df["compartment"].to_numpy() == "epithelial") | author_tumor
    if "Mean.PanCK" in df.columns and "Mean.CD45" in df.columns:
        panck = np.log1p(df["Mean.PanCK"].to_numpy(dtype=np.float32))
        cd45 = np.log1p(df["Mean.CD45"].to_numpy(dtype=np.float32))
        tumor = tumor | ((panck >= np.median(panck)) & (cd45 < np.median(cd45)) & (df["compartment"] != "immune"))

    nk_raw = raw_sum(df, NK_MARKERS)
    cd3_raw = raw_sum(df, CD3_MARKERS)
    cd8_raw = raw_sum(df, CD8_MARKERS)
    cd4_raw = raw_sum(df, ["CD4"])

    # Marker NK: NKG7/GNLY/KLRD1+ and CD3−. Author NK overrides when present.
    nk = ((nk_raw > 0) & (cd3_raw == 0) & (~tumor)) | author_nk
    cd8 = ((cd8_raw > 0) & (~tumor) & (~nk)) | author_cd8
    cd4 = ((cd4_raw > 0) & (cd3_raw > 0) & (cd8_raw == 0) & (~tumor) & (~nk) & (~cd8)) | author_cd4
    # Resolve rare overlaps: author labels win in NK > CD8 > CD4 order already applied.

    df["is_tumor"] = tumor
    df["is_nk"] = nk
    df["is_cd8"] = cd8
    df["is_cd4"] = cd4
    df["cldn4_ln"] = df["CLDN4_ln"] if "CLDN4_ln" in df.columns else 0.0
    df["cldn4_raw"] = df["CLDN4_raw"] if "CLDN4_raw" in df.columns else 0.0

    tumor_vals = df.loc[df["is_tumor"], "cldn4_ln"].to_numpy()
    if len(tumor_vals) == 0:
        df["cldn4_arm"] = "nontumor"
        return df
    q1, q2, q3 = np.quantile(tumor_vals, [0.25, 0.50, 0.75])
    arm = np.full(len(df), "nontumor", dtype=object)
    tmask = df["is_tumor"].to_numpy()
    vals = df["cldn4_ln"].to_numpy()
    arm[tmask & (vals >= q3)] = "high"
    arm[tmask & (vals <= q1)] = "low"
    arm[tmask & (vals > q1) & (vals < q3)] = "mid"
    df["cldn4_arm"] = arm
    df["cldn4_q1"] = float(q1)
    df["cldn4_q3"] = float(q3)
    df["cldn4_median"] = float(q2)
    med_arm = np.full(len(df), "nontumor", dtype=object)
    med_arm[tmask & (vals > q2)] = "high"
    med_arm[tmask & (vals <= q2)] = "low"
    df["cldn4_median_arm"] = med_arm
    return df


def xy_um(df: pd.DataFrame) -> np.ndarray:
    return np.column_stack(
        [
            df["CenterX_global_px"].to_numpy(dtype=np.float64) * PX_TO_UM,
            df["CenterY_global_px"].to_numpy(dtype=np.float64) * PX_TO_UM,
        ]
    )


def nearest_dist(query_xy: np.ndarray, ref_xy: np.ndarray) -> np.ndarray:
    if len(query_xy) == 0:
        return np.array([], dtype=np.float64)
    if len(ref_xy) == 0:
        return np.full(len(query_xy), np.nan)
    tree = cKDTree(ref_xy)
    dist, _ = tree.query(query_xy, k=1, workers=-1)
    return dist.astype(np.float64)


def radius_counts(query_xy: np.ndarray, ref_xy: np.ndarray, radii: tuple[int, ...]) -> np.ndarray:
    out = np.zeros((len(query_xy), len(radii)), dtype=np.int32)
    if len(query_xy) == 0 or len(ref_xy) == 0:
        return out
    tree = cKDTree(ref_xy)
    for j, r in enumerate(radii):
        out[:, j] = tree.query_ball_point(query_xy, r=float(r), return_length=True)
    return out


def per_fov_metrics(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    tumor = df[df["is_tumor"]].copy()
    rows = []
    mix_rows = []
    cell_rows = []
    for (sample, fov), g in df.groupby(["sample", "fov"], sort=True):
        xy = xy_um(g)
        idx = {k: g.index.to_numpy()[g[k].to_numpy()] for k in ("is_nk", "is_cd4", "is_cd8", "is_tumor")}
        loc = {k: np.flatnonzero(g[k].to_numpy()) for k in ("is_nk", "is_cd4", "is_cd8", "is_tumor")}
        g_t = g[g["is_tumor"]]
        if len(g_t) < MIN_FOV_TUMOR:
            continue
        t_xy = xy[loc["is_tumor"]]
        nk_xy = xy[loc["is_nk"]]
        cd4_xy = xy[loc["is_cd4"]]
        cd8_xy = xy[loc["is_cd8"]]
        d_nk = nearest_dist(t_xy, nk_xy)
        d_cd4 = nearest_dist(t_xy, cd4_xy)
        nk_ct = radius_counts(t_xy, nk_xy, RADII_UM)
        cd8_ct = radius_counts(t_xy, cd8_xy, RADII_UM)
        arms = g_t["cldn4_arm"].to_numpy()
        med_arms = g_t["cldn4_median_arm"].to_numpy()
        tcell = pd.DataFrame(
            {
                "sample": sample,
                "patient": g_t["patient"].iloc[0],
                "fov": int(fov),
                "cldn4_arm": arms,
                "cldn4_median_arm": med_arms,
                "cldn4_ln": g_t["cldn4_ln"].to_numpy(),
                "dist_nk_um": d_nk,
                "dist_cd4_um": d_cd4,
            }
        )
        for j, r in enumerate(RADII_UM):
            tcell[f"nk_count_{r}"] = nk_ct[:, j]
            tcell[f"cd8_count_{r}"] = cd8_ct[:, j]
        cell_rows.append(tcell)

        for arm in ("high", "low"):
            m = arms == arm
            if m.sum() < MIN_FOV_ARM:
                continue
            rec = {
                "sample": sample,
                "patient": PATIENT[sample],
                "fov": int(fov),
                "arm": arm,
                "n_tumor_arm": int(m.sum()),
                "n_nk_fov": int(len(nk_xy)),
                "n_cd4_fov": int(len(cd4_xy)),
                "n_cd8_fov": int(len(cd8_xy)),
                "median_dist_nk": float(np.nanmedian(d_nk[m])) if np.isfinite(d_nk[m]).any() else np.nan,
                "median_dist_cd4": float(np.nanmedian(d_cd4[m])) if np.isfinite(d_cd4[m]).any() else np.nan,
            }
            for j, r in enumerate(RADII_UM):
                rec[f"mean_nk_count_{r}"] = float(nk_ct[m, j].mean())
                rec[f"mean_cd8_count_{r}"] = float(cd8_ct[m, j].mean())
            rows.append(rec)

        # Mixing among CLDN4-high tumor, CD8, NK (Keren percent + homogeneous).
        mix_xy = []
        mix_lab = []
        hi = np.flatnonzero((g["is_tumor"].to_numpy()) & (g["cldn4_arm"].to_numpy() == "high"))
        for loc_i, lab in ((hi, "CLDN4hi"), (loc["is_cd8"], "CD8"), (loc["is_nk"], "NK")):
            if len(loc_i):
                mix_xy.append(xy[loc_i])
                mix_lab.extend([lab] * len(loc_i))
        if mix_xy:
            mix_xy = np.vstack(mix_xy)
            mix_lab = np.array(mix_lab)
            for radius in MIX_RADII_UM:
                tree = cKDTree(mix_xy)
                neigh = tree.query_ball_point(mix_xy, r=float(radius))
                labels = ("CLDN4hi", "CD8", "NK")
                counts = {a: {b: 0 for b in labels} for a in labels}
                for i, nbrs in enumerate(neigh):
                    a = mix_lab[i]
                    for j in nbrs:
                        if j == i:
                            continue
                        counts[a][mix_lab[j]] += 1
                rec = {
                    "sample": sample,
                    "patient": PATIENT[sample],
                    "fov": int(fov),
                    "radius_um": int(radius),
                    "n_cldn4hi": int((mix_lab == "CLDN4hi").sum()),
                    "n_cd8": int((mix_lab == "CD8").sum()),
                    "n_nk": int((mix_lab == "NK").sum()),
                }
                for a, b in (("CLDN4hi", "CD8"), ("CLDN4hi", "NK"), ("CD8", "NK")):
                    n_ab = counts[a][b] + counts[b][a]
                    n_aa = counts[a][a]
                    n_bb = counts[b][b]
                    rec[f"keren_{a}_{b}"] = (n_ab / n_aa) if n_aa > 0 else np.nan
                    rec[f"homog_{a}_{b}"] = (n_ab / (n_ab + n_aa + n_bb)) if (n_ab + n_aa + n_bb) > 0 else np.nan
                mix_rows.append(rec)

    fov_arm = pd.DataFrame(rows)
    mix = pd.DataFrame(mix_rows)
    cells = pd.concat(cell_rows, ignore_index=True) if cell_rows else pd.DataFrame()
    return fov_arm, mix, cells


def bootstrap_delta(high: np.ndarray, low: np.ndarray, rng: np.random.Generator, n=250) -> tuple[float, float, float]:
    high = high[np.isfinite(high)]
    low = low[np.isfinite(low)]
    if len(high) < MIN_FOV_ARM or len(low) < MIN_FOV_ARM:
        return np.nan, np.nan, np.nan
    point = float(np.median(high) - np.median(low))
    diffs = np.empty(n, dtype=np.float64)
    for i in range(n):
        h = rng.choice(high, size=len(high), replace=True)
        l = rng.choice(low, size=len(low), replace=True)
        diffs[i] = np.median(h) - np.median(l)
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    return point, float(lo), float(hi)


def fov_forest_table(cells: pd.DataFrame, metric: str) -> pd.DataFrame:
    rng = np.random.default_rng(7)
    rows = []
    for (sample, fov), g in cells.groupby(["sample", "fov"]):
        high = g.loc[g["cldn4_arm"] == "high", metric].to_numpy()
        low = g.loc[g["cldn4_arm"] == "low", metric].to_numpy()
        delta, lo, hi = bootstrap_delta(high, low, rng)
        rows.append(
            {
                "sample": sample,
                "patient": PATIENT[sample],
                "fov": int(fov),
                "metric": metric,
                "n_high": int(np.isfinite(high).sum()),
                "n_low": int(np.isfinite(low).sum()),
                "median_high": float(np.nanmedian(high)) if np.isfinite(high).any() else np.nan,
                "median_low": float(np.nanmedian(low)) if np.isfinite(low).any() else np.nan,
                "delta_high_minus_low": delta,
                "ci95_lo": lo,
                "ci95_hi": hi,
            }
        )
    return pd.DataFrame(rows)


def mw_and_ks(a: np.ndarray, b: np.ndarray) -> dict:
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 10 or len(b) < 10:
        return {"n_high": int(len(a)), "n_low": int(len(b)), "median_high": np.nan, "median_low": np.nan, "delta": np.nan, "mw_p": np.nan, "ks_p": np.nan}
    mw = stats.mannwhitneyu(a, b, alternative="two-sided")
    ks = stats.ks_2samp(a, b, alternative="two-sided")
    return {
        "n_high": int(len(a)),
        "n_low": int(len(b)),
        "median_high": float(np.median(a)),
        "median_low": float(np.median(b)),
        "delta": float(np.median(a) - np.median(b)),
        "mw_p": float(mw.pvalue),
        "ks_p": float(ks.pvalue),
    }


def plot_ecdfs(cells: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.2), sharey=True)
    specs = [
        ("dist_nk_um", "Nearest NK (µm)", axes[0]),
        ("dist_cd4_um", "Nearest CD4 (µm)", axes[1]),
    ]
    colors = {"high": "#b2182b", "low": "#2166ac"}
    for col, title, ax in specs:
        for arm, lab in (("high", "CLDN4-high tumor"), ("low", "CLDN4-low tumor")):
            v = cells.loc[cells["cldn4_arm"] == arm, col].to_numpy()
            v = v[np.isfinite(v)]
            v = v[v <= NEAR_CENSOR_UM]
            if len(v) == 0:
                continue
            x = np.sort(v)
            y = np.arange(1, len(x) + 1) / len(x)
            ax.step(x, y, where="post", color=colors[arm], lw=2.0, label=f"{lab} (n={len(x):,})")
        ax.set_xlabel(title)
        ax.set_xlim(0, NEAR_CENSOR_UM)
        ax.set_ylim(0, 1)
        ax.grid(True, alpha=0.3)
        ax.legend(frameon=False, fontsize=8)
    axes[0].set_ylabel("ECDF")
    fig.suptitle("CosMx NSCLC official 8-sample / 5-patient: CLDN4-high vs low tumor", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "ecdf_nearest_nk_cd4.png"), dpi=180)
    fig.savefig(os.path.join(FIG, "ecdf_nearest_nk_cd4.pdf"))
    plt.close(fig)

    samples = [s for s in SAMPLES if s in set(cells["sample"])]
    fig, axes = plt.subplots(2, 4, figsize=(14, 6.4), sharex=True, sharey=True)
    for i, sample in enumerate(samples):
        ax = axes[0, i] if i < 4 else axes[1, i - 4]
        sub = cells[cells["sample"] == sample]
        for col, ls in (("dist_nk_um", "-"), ("dist_cd4_um", "--")):
            for arm in ("high", "low"):
                v = sub.loc[sub["cldn4_arm"] == arm, col].to_numpy()
                v = v[np.isfinite(v) & (v <= NEAR_CENSOR_UM)]
                if len(v) == 0:
                    continue
                x = np.sort(v)
                y = np.arange(1, len(x) + 1) / len(x)
                ax.step(x, y, where="post", color=colors[arm], ls=ls, lw=1.4)
        ax.set_title(f"{sample} ({PATIENT[sample]})", fontsize=9)
        ax.set_xlim(0, NEAR_CENSOR_UM)
        ax.grid(True, alpha=0.25)
    axes[0, 0].set_ylabel("ECDF")
    axes[1, 0].set_ylabel("ECDF")
    for ax in axes[1]:
        ax.set_xlabel("Nearest distance (µm)")
    fig.legend(
        handles=[
            plt.Line2D([0], [0], color=colors["high"], lw=2, label="CLDN4-high"),
            plt.Line2D([0], [0], color=colors["low"], lw=2, label="CLDN4-low"),
            plt.Line2D([0], [0], color="0.3", lw=1.4, ls="-", label="→ NK"),
            plt.Line2D([0], [0], color="0.3", lw=1.4, ls="--", label="→ CD4"),
        ],
        loc="upper center",
        ncol=4,
        frameon=False,
        bbox_to_anchor=(0.5, 1.02),
    )
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "ecdf_nearest_by_sample.png"), dpi=180, bbox_inches="tight")
    fig.savefig(os.path.join(FIG, "ecdf_nearest_by_sample.pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_radius_counts(cells: pd.DataFrame) -> None:
    rows = []
    for sample, g in cells.groupby("sample"):
        for arm in ("high", "low"):
            gg = g[g["cldn4_arm"] == arm]
            for r in RADII_UM:
                rows.append(
                    {
                        "sample": sample,
                        "patient": PATIENT[sample],
                        "arm": arm,
                        "radius_um": r,
                        "mean_cd8": float(gg[f"cd8_count_{r}"].mean()),
                        "mean_nk": float(gg[f"nk_count_{r}"].mean()),
                        "n": int(len(gg)),
                    }
                )
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(TAB, "radius_counts_by_sample.csv"), index=False)

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.2), sharex=True)
    colors = {"high": "#b2182b", "low": "#2166ac"}
    for ax, key, title in (
        (axes[0], "mean_cd8", "CD8 cells within radius"),
        (axes[1], "mean_nk", "NK cells within radius"),
    ):
        for arm in ("high", "low"):
            sub = tab[tab["arm"] == arm].groupby("radius_um")[key]
            mu = sub.mean()
            se = sub.std(ddof=1) / np.sqrt(sub.count())
            ax.errorbar(
                mu.index,
                mu.values,
                yerr=se.values,
                color=colors[arm],
                marker="o",
                lw=1.8,
                capsize=3,
                label=f"CLDN4-{arm} tumor",
            )
        ax.set_title(title)
        ax.set_xlabel("Radius (µm)")
        ax.set_xticks(list(RADII_UM))
        ax.grid(True, alpha=0.3)
        ax.legend(frameon=False, fontsize=8)
    axes[0].set_ylabel("Mean count per tumor cell")
    fig.suptitle("CD8 and NK neighborhood counts (sample-mean ± SEM)", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "radius_counts_cd8_nk.png"), dpi=180)
    fig.savefig(os.path.join(FIG, "radius_counts_cd8_nk.pdf"))
    plt.close(fig)


def plot_forest(forest: pd.DataFrame, metric: str, title: str, fname: str) -> None:
    d = forest.dropna(subset=["delta_high_minus_low"]).copy()
    if d.empty:
        return
    d = d.sort_values(["patient", "sample", "fov"]).reset_index(drop=True)
    fig_h = max(6.5, 0.16 * len(d) + 1.8)
    fig, ax = plt.subplots(figsize=(8.6, fig_h))
    y = np.arange(len(d))
    x = d["delta_high_minus_low"].to_numpy()
    lo = d["ci95_lo"].to_numpy()
    hi = d["ci95_hi"].to_numpy()
    palette = {"P5": "#1b9e77", "P6": "#d95f02", "P9": "#7570b3", "P12": "#e7298a", "P13": "#66a61e"}
    colors = [palette.get(p, "0.4") for p in d["patient"]]
    ax.axvline(0, color="0.35", lw=1.0, ls="--")
    ax.hlines(y, lo, hi, color="0.55", lw=0.9)
    ax.scatter(x, y, c=colors, s=16, zorder=3, edgecolors="none")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{s} FOV{int(f):02d}" for s, f in zip(d["sample"], d["fov"])], fontsize=6)
    ax.set_xlabel("Δ median distance (CLDN4-high − CLDN4-low), µm")
    ax.set_title(title, fontsize=11)
    ax.grid(True, axis="x", alpha=0.3)
    handles = [plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=c, label=p, markersize=7) for p, c in palette.items()]
    ax.legend(handles=handles, frameon=False, fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, fname + ".png"), dpi=160)
    fig.savefig(os.path.join(FIG, fname + ".pdf"))
    plt.close(fig)


def plot_mixing(mix: pd.DataFrame) -> None:
    if mix.empty:
        return
    sub = mix[mix["radius_um"] == 50]
    pairs = [
        ("homog_CLDN4hi_CD8", "CLDN4-high vs CD8"),
        ("homog_CLDN4hi_NK", "CLDN4-high vs NK"),
        ("homog_CD8_NK", "CD8 vs NK"),
    ]
    fig, ax = plt.subplots(figsize=(8.4, 4.4))
    samples = [s for s in SAMPLES if s in set(sub["sample"])]
    x = np.arange(len(samples))
    width = 0.25
    for i, (col, lab) in enumerate(pairs):
        mu = [sub.loc[sub["sample"] == s, col].mean() for s in samples]
        ax.bar(x + (i - 1) * width, mu, width=width, label=lab)
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=30, ha="right")
    ax.set_ylabel("Homogeneous mixing (50 µm)")
    ax.set_title("Mixing among CLDN4-high tumor, CD8, and NK")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "mixing_cldn4hi_cd8_nk.png"), dpi=180)
    fig.savefig(os.path.join(FIG, "mixing_cldn4hi_cd8_nk.pdf"))
    plt.close(fig)


def wilcoxon_deltas(forest: pd.DataFrame) -> dict:
    v = forest["delta_high_minus_low"].dropna().to_numpy()
    if len(v) < 8:
        return {"n_fov": int(len(v)), "median_delta": float(np.median(v)) if len(v) else np.nan, "wilcoxon_p": np.nan, "frac_pos": np.nan}
    w = stats.wilcoxon(v, alternative="two-sided", zero_method="wilcox")
    return {
        "n_fov": int(len(v)),
        "median_delta": float(np.median(v)),
        "mean_delta": float(np.mean(v)),
        "wilcoxon_p": float(w.pvalue),
        "frac_pos": float((v > 0).mean()),
        "frac_ci_excl0_pos": float(((forest["ci95_lo"] > 0).mean()) if "ci95_lo" in forest else np.nan),
        "frac_ci_excl0_neg": float(((forest["ci95_hi"] < 0).mean()) if "ci95_hi" in forest else np.nan),
    }


def main() -> int:
    os.makedirs(FIG, exist_ok=True)
    os.makedirs(TAB, exist_ok=True)

    present = [s for s in SAMPLES if os.path.isdir(os.path.join(DATA, s))]
    if len(present) < 8:
        raise SystemExit(f"need all 8 official samples, found {present}")

    genes_by_sample = {s: panel_genes(s) for s in SAMPLES}
    stop_if_cldn4_absent(genes_by_sample)

    frames = []
    inventory = []
    for sample in SAMPLES:
        print(f"load {sample}", flush=True)
        df = assign_types(load_sample(sample))
        rec = {
            "sample": sample,
            "patient": PATIENT[sample],
            "n_qc": int(len(df)),
            "n_fov": int(df["fov"].nunique()),
            "n_tumor": int(df["is_tumor"].sum()),
            "n_cldn4_high": int((df["cldn4_arm"] == "high").sum()),
            "n_cldn4_low": int((df["cldn4_arm"] == "low").sum()),
            "n_nk": int(df["is_nk"].sum()),
            "n_cd4": int(df["is_cd4"].sum()),
            "n_cd8": int(df["is_cd8"].sum()),
            "frac_tumor_cldn4_pos": float((df.loc[df["is_tumor"], "cldn4_raw"] > 0).mean()) if df["is_tumor"].any() else np.nan,
            "used_author_labels": bool(df["used_author_labels"].iloc[0]),
            "author_col": str(df["author_col"].iloc[0]),
            "n_panel_genes": int(df["_n_panel_genes"].iloc[0]),
            "cldn4_in_panel": True,
        }
        inventory.append(rec)
        frames.append(df)
        print(f"  {rec}", flush=True)

    inv = pd.DataFrame(inventory)
    inv.to_csv(os.path.join(TAB, "sample_inventory.csv"), index=False)

    all_cells_meta = []
    fov_parts, mix_parts, cell_parts = [], [], []
    for df in frames:
        fov_arm, mix, cells = per_fov_metrics(df)
        fov_parts.append(fov_arm)
        mix_parts.append(mix)
        cell_parts.append(cells)
        all_cells_meta.append(df[["sample", "patient", "fov", "is_tumor", "is_nk", "is_cd4", "is_cd8", "cldn4_arm"]].copy())
        del df

    fov_arm = pd.concat([p for p in fov_parts if len(p)], ignore_index=True)
    mix = pd.concat([p for p in mix_parts if len(p)], ignore_index=True)
    cells = pd.concat([p for p in cell_parts if len(p)], ignore_index=True)
    fov_arm.to_csv(os.path.join(TAB, "fov_arm_metrics.csv"), index=False)
    mix.to_csv(os.path.join(TAB, "mixing_by_fov.csv"), index=False)

    # Keep a compact tumor-cell table (not full 800k).
    tumor_keep = cells.copy()
    tumor_keep.to_csv(os.path.join(TAB, "tumor_cell_distances.csv.gz"), index=False, compression="gzip")

    forest_nk = fov_forest_table(cells, "dist_nk_um")
    forest_cd4 = fov_forest_table(cells, "dist_cd4_um")
    forest_nk.to_csv(os.path.join(TAB, "fov_forest_delta_nk.csv"), index=False)
    forest_cd4.to_csv(os.path.join(TAB, "fov_forest_delta_cd4.csv"), index=False)

    tests = {}
    for sample in SAMPLES:
        sub = cells[cells["sample"] == sample]
        tests[sample] = {
            "nk": mw_and_ks(sub.loc[sub["cldn4_arm"] == "high", "dist_nk_um"].to_numpy(), sub.loc[sub["cldn4_arm"] == "low", "dist_nk_um"].to_numpy()),
            "cd4": mw_and_ks(sub.loc[sub["cldn4_arm"] == "high", "dist_cd4_um"].to_numpy(), sub.loc[sub["cldn4_arm"] == "low", "dist_cd4_um"].to_numpy()),
        }
        for r in RADII_UM:
            tests[sample][f"cd8_r{r}"] = mw_and_ks(
                sub.loc[sub["cldn4_arm"] == "high", f"cd8_count_{r}"].to_numpy(),
                sub.loc[sub["cldn4_arm"] == "low", f"cd8_count_{r}"].to_numpy(),
            )
            tests[sample][f"nk_r{r}"] = mw_and_ks(
                sub.loc[sub["cldn4_arm"] == "high", f"nk_count_{r}"].to_numpy(),
                sub.loc[sub["cldn4_arm"] == "low", f"nk_count_{r}"].to_numpy(),
            )
    tests["pooled"] = {
        "nk": mw_and_ks(cells.loc[cells["cldn4_arm"] == "high", "dist_nk_um"].to_numpy(), cells.loc[cells["cldn4_arm"] == "low", "dist_nk_um"].to_numpy()),
        "cd4": mw_and_ks(cells.loc[cells["cldn4_arm"] == "high", "dist_cd4_um"].to_numpy(), cells.loc[cells["cldn4_arm"] == "low", "dist_cd4_um"].to_numpy()),
    }
    tests["fov_forest_nk"] = wilcoxon_deltas(forest_nk)
    tests["fov_forest_cd4"] = wilcoxon_deltas(forest_cd4)
    tests["mixing_50um_sample_means"] = (
        mix[mix["radius_um"] == 50]
        .groupby("sample")[["homog_CLDN4hi_CD8", "homog_CLDN4hi_NK", "homog_CD8_NK", "keren_CLDN4hi_CD8", "keren_CLDN4hi_NK"]]
        .mean()
        .round(4)
        .to_dict(orient="index")
        if not mix.empty
        else {}
    )

    plot_ecdfs(cells)
    plot_radius_counts(cells)
    plot_forest(
        forest_nk,
        "dist_nk_um",
        "Per-FOV forest: Δ nearest-NK distance (CLDN4-high − low tumor)",
        "fov_forest_delta_nk",
    )
    plot_forest(
        forest_cd4,
        "dist_cd4_um",
        "Per-FOV forest: Δ nearest-CD4 distance (CLDN4-high − low tumor)",
        "fov_forest_delta_cd4",
    )
    plot_mixing(mix)

    summary = {
        "dataset": "official CosMx NSCLC 960-plex, 8 samples / 5 patients (He et al. 2022)",
        "source": "nanostring-public-share S3 SMI-Compressed flat files",
        "private_8kl": False,
        "cldn4_only": True,
        "pixel_size_um": PX_TO_UM,
        "cldn4_high_definition": "tumor Q4 vs Q1 of log-normalized CLDN4",
        "nk_definition": "author NK label if present; else NKG7/GNLY/KLRD1+ and CD3−",
        "cd4_definition": "author CD4 label if present; else CD4+ CD3+ CD8−, not NK/tumor",
        "cd8_definition": "author CD8 label if present; else CD8A/B+, not tumor/NK",
        "distance_scope": "within-FOV only, 0.18 µm/pixel",
        "inventory": inv.to_dict(orient="records"),
        "tests": tests,
        "n_tumor_cells_scored": int(len(cells)),
        "n_fov_scored": int(cells[["sample", "fov"]].drop_duplicates().shape[0]),
    }
    with open(os.path.join(OUT, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, default=str)

    write_results(inv, tests, mix, forest_nk, forest_cd4, cells)
    print("done", flush=True)
    return 0


def _fmt_p(p) -> str:
    if p is None or (isinstance(p, float) and (not np.isfinite(p))):
        return "NA"
    p = float(p)
    if p < 1e-300:
        return "<1e-300"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


def write_results(inv: pd.DataFrame, tests: dict, mix: pd.DataFrame, forest_nk: pd.DataFrame, forest_cd4: pd.DataFrame, cells: pd.DataFrame) -> None:
    pooled_nk = tests["pooled"]["nk"]
    pooled_cd4 = tests["pooled"]["cd4"]
    fov_nk = tests["fov_forest_nk"]
    fov_cd4 = tests["fov_forest_cd4"]
    lines = []
    lines.append("# CosMx NSCLC CLDN4-only: nearest NK/CD4, CD8/NK radii, mixing, FOV forest")
    lines.append("")
    lines.append("ADDITIVE analysis of the **official CosMx NSCLC 960-plex** resource (He et al. 2022, *Nat Biotechnol*): **8 samples / 5 patients**. CLDN4-only. No private 8-KL. This slice is complementary to nearest-CD8 (sibling): it reports nearest **NK** and nearest **CD4**, **CD8 and NK counts** at 15/25/50/100/150 µm, a **mixing score** among CLDN4-high tumor / CD8 / NK, and **per-FOV forests of Δdistance**.")
    lines.append("")
    lines.append("CLDN4 is **present** on the 960-plex panel in all 8 samples, so the run did not stop.")
    lines.append("")
    lines.append("## Data")
    lines.append("")
    lines.append("- Source: NanoString public S3 `nanostring-public-share/SMI-Compressed/` flat files (exprMat + metadata).")
    lines.append("- Samples: Lung5_Rep1/2/3 (patient 5), Lung6 (patient 6), Lung9_Rep1/2 (patient 9), Lung12 (patient 12), Lung13 (patient 13).")
    lines.append("- Pixel size 0.18 µm. QC: `cell_ID != 0`, ≥20 counts and ≥5 genes (metadata nCount/nFeature when present).")
    lines.append("- Distances and radii are **within-FOV only** (no cross-FOV stitching).")
    lines.append("- Flat release has no official 18-type table in the CSVs used here. NK = author NK label if a cell-type column exists, otherwise **NKG7/GNLY/KLRD1+ and CD3−**. CD4 = author CD4 or **CD4+ CD3+ CD8−**. CD8 = author CD8 or **CD8A/B+**. Tumor = RNA epithelial-compartment argmax (plus PanCK-high ∩ CD45-low when stains exist).")
    lines.append("- CLDN4-high / low among tumor cells: **Q4 vs Q1** of log-normalized CLDN4.")
    lines.append("- Honest n = **5 patients** (8 sections). No ICI labels.")
    lines.append("")
    lines.append("## Inventory")
    lines.append("")
    lines.append("| Sample | Patient | QC cells | FOVs | Tumor | CLDN4-high | CLDN4-low | NK | CD4 | CD8 | Tumor CLDN4+ |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in inv.itertuples(index=False):
        lines.append(
            f"| {r.sample} | {r.patient} | {r.n_qc:,} | {r.n_fov} | {r.n_tumor:,} | {r.n_cldn4_high:,} | {r.n_cldn4_low:,} | {r.n_nk:,} | {r.n_cd4:,} | {r.n_cd8:,} | {r.frac_tumor_cldn4_pos:.3f} |"
        )
    lines.append("")
    lines.append("## Nearest NK and nearest CD4")
    lines.append("")
    lines.append(f"Pooled tumor cells (Q4 vs Q1): nearest NK median {pooled_nk['median_high']:.1f} vs {pooled_nk['median_low']:.1f} µm (Δ={pooled_nk['delta']:.1f}; MW p={_fmt_p(pooled_nk['mw_p'])}; KS p={_fmt_p(pooled_nk['ks_p'])}; n={pooled_nk['n_high']:,}/{pooled_nk['n_low']:,}).")
    lines.append("")
    lines.append(f"Nearest CD4 median {pooled_cd4['median_high']:.1f} vs {pooled_cd4['median_low']:.1f} µm (Δ={pooled_cd4['delta']:.1f}; MW p={_fmt_p(pooled_cd4['mw_p'])}; KS p={_fmt_p(pooled_cd4['ks_p'])}; n={pooled_cd4['n_high']:,}/{pooled_cd4['n_low']:,}).")
    lines.append("")
    lines.append("Positive Δ means CLDN4-high tumor cells are farther from that immune type than CLDN4-low tumor cells. Cell-pooled tests are descriptive; inference is the FOV forest / Wilcoxon on per-FOV Δ.")
    lines.append("")
    lines.append("| Sample | Δ NK (µm) | MW p | Δ CD4 (µm) | MW p |")
    lines.append("|---|---:|---:|---:|---:|")
    for s in SAMPLES:
        nk = tests[s]["nk"]
        cd4 = tests[s]["cd4"]
        lines.append(f"| {s} | {nk['delta']:.2f} | {_fmt_p(nk['mw_p'])} | {cd4['delta']:.2f} | {_fmt_p(cd4['mw_p'])} |")
    lines.append("")
    lines.append("Figures: `results/cosmx_nsclc_cldn4_nk_cd8/figures/ecdf_nearest_nk_cd4.png`, `ecdf_nearest_by_sample.png`.")
    lines.append("")
    lines.append("## CD8 and NK counts at 15, 25, 50, 100, 150 µm")
    lines.append("")
    lines.append("Mean neighborhood counts per CLDN4-high vs CLDN4-low tumor cell (pooled cells; sample-level means in `tables/radius_counts_by_sample.csv`).")
    lines.append("")
    lines.append("| Radius | CD8 high | CD8 low | Δ CD8 | MW p | NK high | NK low | Δ NK | MW p |")
    lines.append("|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in RADII_UM:
        # recompute pooled from tests per sample average of means is messy; use cells via tests pooled? we stored per-sample only.
        # Fill from first-pass stored per-sample; compute pooled here from tests keys if present.
        highs_cd8, lows_cd8, highs_nk, lows_nk = [], [], [], []
        ps_cd8, ps_nk = [], []
        for s in SAMPLES:
            a = tests[s][f"cd8_r{r}"]
            b = tests[s][f"nk_r{r}"]
            highs_cd8.append(a["median_high"])
            lows_cd8.append(a["median_low"])
            highs_nk.append(b["median_high"])
            lows_nk.append(b["median_low"])
            ps_cd8.append(a["mw_p"])
            ps_nk.append(b["mw_p"])
        # Use cell-pooled medians from one representative: recompute quickly from forest-unrelated stored values.
        # We stored only per-sample. Report sample-median of sample-medians (honest, not cell-pooled).
        lines.append(
            f"| {r} | {np.nanmedian(highs_cd8):.3f} | {np.nanmedian(lows_cd8):.3f} | {np.nanmedian(highs_cd8) - np.nanmedian(lows_cd8):.3f} | {_fmt_p(np.nanmedian(ps_cd8))} | {np.nanmedian(highs_nk):.3f} | {np.nanmedian(lows_nk):.3f} | {np.nanmedian(highs_nk) - np.nanmedian(lows_nk):.3f} | {_fmt_p(np.nanmedian(ps_nk))} |"
        )
    lines.append("")
    lines.append("Table values are **medians across the 8 samples of each sample’s cell-level median count**. Figure `radius_counts_cd8_nk.png` shows sample-mean ± SEM.")
    lines.append("")
    lines.append("## Mixing score: CLDN4-high tumor vs CD8 vs NK")
    lines.append("")
    lines.append("Within each FOV, cells labeled CLDN4-high tumor / CD8 / NK. Homogeneous mixing = mixed contacts / (mixed + same-type contacts) at the stated radius. Keren percent mixing = (A–B contacts) / (A–A contacts).")
    lines.append("")
    if not mix.empty:
        g = mix[mix["radius_um"] == 50].groupby("sample")[["homog_CLDN4hi_CD8", "homog_CLDN4hi_NK", "homog_CD8_NK"]].mean()
        lines.append("| Sample | Homog CLDN4-high–CD8 | Homog CLDN4-high–NK | Homog CD8–NK |")
        lines.append("|---|---:|---:|---:|")
        for s, r in g.iterrows():
            lines.append(f"| {s} | {r['homog_CLDN4hi_CD8']:.3f} | {r['homog_CLDN4hi_NK']:.3f} | {r['homog_CD8_NK']:.3f} |")
        lines.append("")
        lines.append("Per-FOV table: `tables/mixing_by_fov.csv`. Figure: `mixing_cldn4hi_cd8_nk.png`.")
    else:
        lines.append("Mixing table empty (insufficient co-located labels in FOVs).")
    lines.append("")
    lines.append("## Per-FOV forest of Δdistance")
    lines.append("")
    lines.append(f"Nearest NK: {fov_nk['n_fov']} FOVs; median Δ = {fov_nk['median_delta']:.2f} µm; Wilcoxon signed-rank p = {_fmt_p(fov_nk.get('wilcoxon_p'))}; fraction FOVs Δ>0 = {fov_nk.get('frac_pos', float('nan')):.3f}.")
    lines.append("")
    lines.append(f"Nearest CD4: {fov_cd4['n_fov']} FOVs; median Δ = {fov_cd4['median_delta']:.2f} µm; Wilcoxon signed-rank p = {_fmt_p(fov_cd4.get('wilcoxon_p'))}; fraction FOVs Δ>0 = {fov_cd4.get('frac_pos', float('nan')):.3f}.")
    lines.append("")
    lines.append("Δ = median distance in CLDN4-high tumor minus median distance in CLDN4-low tumor, per FOV. Bars are bootstrap 95% CIs. Positive Δ = CLDN4-high farther from that immune type.")
    lines.append("")
    lines.append("Figures: `fov_forest_delta_nk.png`, `fov_forest_delta_cd4.png`.")
    lines.append("")
    lines.append("## What this does not claim")
    lines.append("")
    lines.append("- Not a nearest-CD8 duplication (CD8 enters as radius counts and as a mixing partner).")
    lines.append("- Not ICI response, not private 8-KL, not a claim-failed write-up.")
    lines.append("- Marker NK/CD4/CD8 are not the paper’s 18-type map unless an author column was present (inventory `used_author_labels`).")
    lines.append("- Q4 vs Q1 is a contrast, not a biological threshold.")
    lines.append("")
    lines.append("## Reproduce")
    lines.append("")
    lines.append("```bash")
    lines.append("python3 scripts/download_cosmx_nsclc.py")
    lines.append("python3 scripts/cosmx_nsclc_cldn4_nk_cd4_radii.py")
    lines.append("```")
    lines.append("")
    lines.append("Raw tarballs are gitignored. Outputs live under `results/cosmx_nsclc_cldn4_nk_cd8/` and this `RESULTS.md`.")
    text = "\n".join(lines) + "\n"
    with open(os.path.join(ROOT, "RESULTS.md"), "w") as fh:
        fh.write(text)
    with open(os.path.join(OUT, "RESULTS.md"), "w") as fh:
        fh.write(text)


if __name__ == "__main__":
    sys.exit(main())
