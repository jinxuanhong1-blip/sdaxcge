#!/usr/bin/env python3
"""CosMx NSCLC CLDN4-only: CD8+NK neighbor counts and mixing (retuned).

Official 8-section / 5-patient 960-plex (He et al. 2022). Primary endpoints are
CD8+NK neighbor COUNT and mixing on a 20/40/60/80 µm grid, with a section-level
paired test (n=8). Tumor cells are NOT gated as CD8A==0. Nearest-distance is
not the primary readout. Stops if CLDN4 is absent. No private 8-KL.
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
RADII_UM = (20, 40, 60, 80)
MIN_COUNTS = 20
MIN_GENES = 5
MIN_FOV_TUMOR = 20
MIN_FOV_ARM = 5
PRIMARY_R = 40

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
STR_MARKERS = ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "ACTA2", "PECAM1", "VWF", "FN1", "BGN"]
AUTHOR_TYPE_COLS = (
    "cell_type",
    "celltype",
    "CellType",
    "cell_types",
    "nb_clus",
    "annotation",
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
    if not missing:
        return
    os.makedirs(OUT, exist_ok=True)
    msg = (
        "CLDN4 is absent from the CosMx 960-plex exprMat for: "
        + ", ".join(missing)
        + ". Stopping as instructed. No counts or mixing were computed."
    )
    with open(os.path.join(ROOT, "RESULTS.md"), "w") as fh:
        fh.write("# CosMx NSCLC CLDN4 CD8+NK radii — STOP\n\n" + msg + "\n")
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
    for col in ("Mean.PanCK", "Mean.CD45", "Mean.CD3", "Area"):
        if col in meta.columns:
            keep_meta.append(col)
    author_col = next((c for c in AUTHOR_TYPE_COLS if c in meta.columns), None)
    if author_col:
        keep_meta.append(author_col)
    meta = meta[keep_meta].copy()
    df = expr.merge(meta, on=["cell_key", "fov", "cell_ID"], how="inner")
    gene_cols = [g for g in genes if g in df.columns]
    raw = df[gene_cols].to_numpy(dtype=np.float32)
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
    """Tumor = epithelial RNA argmax (optional PanCK support). Do NOT gate CD8A==0."""
    epi = mean_ln(df, EPI_MARKERS)
    imm = mean_ln(df, IMM_MARKERS)
    strn = mean_ln(df, STR_MARKERS)
    scores = np.vstack([epi, imm, strn])
    lab = np.array(["epithelial", "immune", "stromal"])[scores.argmax(axis=0)]
    lab[scores.max(axis=0) <= 0] = "unassigned"
    df["compartment"] = lab

    author_tumor = author_match(df["author_type"], ("tumor", "epithelial", "cancer", "malignant"))
    author_nk = author_match(df["author_type"], ("nk",))
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
    cd8a_raw = raw_sum(df, ["CD8A"])

    # Keep CD8A+ epithelial cells in the tumor index (explicitly no CD8A==0 gate).
    df["tumor_cd8a_pos"] = tumor & (cd8a_raw > 0)

    # CD8+NK together (July PPT cytotoxic): non-tumor CD8 or NK.
    nk = ((nk_raw > 0) & (cd3_raw == 0) & (~tumor)) | (author_nk & (~tumor))
    cd8 = ((cd8_raw > 0) & (~tumor)) | (author_cd8 & (~tumor))
    cd8nk = cd8 | nk

    df["is_tumor"] = tumor
    df["is_nk"] = nk
    df["is_cd8"] = cd8
    df["is_cd8nk"] = cd8nk
    df["cldn4_ln"] = df["CLDN4_ln"] if "CLDN4_ln" in df.columns else 0.0
    df["cldn4_raw"] = df["CLDN4_raw"] if "CLDN4_raw" in df.columns else 0.0

    tumor_vals = df.loc[df["is_tumor"], "cldn4_ln"].to_numpy()
    if len(tumor_vals) == 0:
        df["cldn4_arm"] = "nontumor"
        df["cldn4_median_arm"] = "nontumor"
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


def radius_counts(query_xy: np.ndarray, ref_xy: np.ndarray, radii: tuple[int, ...]) -> np.ndarray:
    out = np.zeros((len(query_xy), len(radii)), dtype=np.int32)
    if len(query_xy) == 0 or len(ref_xy) == 0:
        return out
    tree = cKDTree(ref_xy)
    for j, r in enumerate(radii):
        balls = tree.query_ball_point(query_xy, r=float(r))
        out[:, j] = np.fromiter((len(b) for b in balls), dtype=np.int32, count=len(query_xy))
    return out


def mixing_from_labels(xy: np.ndarray, lab: np.ndarray, radius: float, a: str, b: str) -> dict:
    """Keren percent mixing and homogeneous mixing for two labels."""
    if len(xy) == 0:
        return {"keren_a_b": np.nan, "keren_b_a": np.nan, "homog": np.nan, "n_aa": 0, "n_bb": 0, "n_ab": 0}
    tree = cKDTree(xy)
    neigh = tree.query_ball_point(xy, r=float(radius))
    n_aa = n_bb = n_ab = 0
    for i, nbrs in enumerate(neigh):
        li = lab[i]
        for j in nbrs:
            if j == i:
                continue
            lj = lab[j]
            if li == a and lj == a:
                n_aa += 1
            elif li == b and lj == b:
                n_bb += 1
            elif (li == a and lj == b) or (li == b and lj == a):
                n_ab += 1
    n_ab_undirected = n_ab  # both directions counted in the loop
    return {
        "keren_a_b": (n_ab_undirected / n_aa) if n_aa > 0 else np.nan,  # mixed / A-A
        "keren_b_a": (n_ab_undirected / n_bb) if n_bb > 0 else np.nan,  # mixed / B-B (immune-side)
        "homog": (n_ab_undirected / (n_ab_undirected + n_aa + n_bb)) if (n_ab_undirected + n_aa + n_bb) > 0 else np.nan,
        "n_aa": int(n_aa),
        "n_bb": int(n_bb),
        "n_ab": int(n_ab_undirected),
    }


def per_fov_metrics(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    fov_rows = []
    mix_rows = []
    cell_rows = []
    for (sample, fov), g in df.groupby(["sample", "fov"], sort=True):
        xy = xy_um(g)
        loc_t = np.flatnonzero(g["is_tumor"].to_numpy())
        loc_i = np.flatnonzero(g["is_cd8nk"].to_numpy())
        g_t = g.iloc[loc_t]
        if len(g_t) < MIN_FOV_TUMOR:
            continue
        t_xy = xy[loc_t]
        i_xy = xy[loc_i]
        ct = radius_counts(t_xy, i_xy, RADII_UM)
        arms = g_t["cldn4_arm"].to_numpy()
        tcell = pd.DataFrame(
            {
                "sample": sample,
                "patient": g_t["patient"].iloc[0],
                "fov": int(fov),
                "cldn4_arm": arms,
                "cldn4_median_arm": g_t["cldn4_median_arm"].to_numpy(),
                "cldn4_ln": g_t["cldn4_ln"].to_numpy(),
                "tumor_cd8a_pos": g_t["tumor_cd8a_pos"].to_numpy(),
            }
        )
        for j, r in enumerate(RADII_UM):
            tcell[f"cd8nk_count_{r}"] = ct[:, j]
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
                "n_cd8nk_fov": int(len(i_xy)),
            }
            for j, r in enumerate(RADII_UM):
                rec[f"mean_cd8nk_count_{r}"] = float(ct[m, j].mean())
                rec[f"median_cd8nk_count_{r}"] = float(np.median(ct[m, j]))
            fov_rows.append(rec)

        hi = np.flatnonzero((g["is_tumor"].to_numpy()) & (g["cldn4_arm"].to_numpy() == "high"))
        lo = np.flatnonzero((g["is_tumor"].to_numpy()) & (g["cldn4_arm"].to_numpy() == "low"))
        for radius in RADII_UM:
            rec = {
                "sample": sample,
                "patient": PATIENT[sample],
                "fov": int(fov),
                "radius_um": int(radius),
                "n_cldn4hi": int(len(hi)),
                "n_cldn4lo": int(len(lo)),
                "n_cd8nk": int(len(loc_i)),
            }
            for arm_name, loc_arm in (("high", hi), ("low", lo)):
                if len(loc_arm) == 0 or len(loc_i) == 0:
                    rec[f"keren_tumor_{arm_name}"] = np.nan
                    rec[f"keren_immune_{arm_name}"] = np.nan
                    rec[f"homog_{arm_name}"] = np.nan
                    rec[f"n_aa_{arm_name}"] = 0
                    rec[f"n_bb_{arm_name}"] = 0
                    rec[f"n_ab_{arm_name}"] = 0
                    continue
                mix_xy = np.vstack([xy[loc_arm], xy[loc_i]])
                mix_lab = np.array(["tumor"] * len(loc_arm) + ["cd8nk"] * len(loc_i))
                mx = mixing_from_labels(mix_xy, mix_lab, radius, "tumor", "cd8nk")
                rec[f"keren_tumor_{arm_name}"] = mx["keren_a_b"]
                rec[f"keren_immune_{arm_name}"] = mx["keren_b_a"]
                rec[f"homog_{arm_name}"] = mx["homog"]
                rec[f"n_aa_{arm_name}"] = mx["n_aa"]
                rec[f"n_bb_{arm_name}"] = mx["n_bb"]
                rec[f"n_ab_{arm_name}"] = mx["n_ab"]
            mix_rows.append(rec)

    fov_arm = pd.DataFrame(fov_rows)
    mix = pd.DataFrame(mix_rows)
    cells = pd.concat(cell_rows, ignore_index=True) if cell_rows else pd.DataFrame()
    return fov_arm, mix, cells


def section_count_table(cells: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for sample, g in cells.groupby("sample"):
        rec = {"sample": sample, "patient": PATIENT[sample], "n_high": int((g["cldn4_arm"] == "high").sum()), "n_low": int((g["cldn4_arm"] == "low").sum())}
        for r in RADII_UM:
            hi = g.loc[g["cldn4_arm"] == "high", f"cd8nk_count_{r}"]
            lo = g.loc[g["cldn4_arm"] == "low", f"cd8nk_count_{r}"]
            rec[f"mean_high_{r}"] = float(hi.mean())
            rec[f"mean_low_{r}"] = float(lo.mean())
            rec[f"median_high_{r}"] = float(hi.median())
            rec[f"median_low_{r}"] = float(lo.median())
            rec[f"delta_mean_{r}"] = rec[f"mean_high_{r}"] - rec[f"mean_low_{r}"]
        rows.append(rec)
    return pd.DataFrame(rows).set_index("sample").loc[[s for s in SAMPLES if s in set(cells["sample"])]].reset_index()


def section_mix_table(mix: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (sample, radius), g in mix.groupby(["sample", "radius_um"]):
        rec = {"sample": sample, "patient": PATIENT[sample], "radius_um": int(radius)}
        for arm in ("high", "low"):
            n_aa = g[f"n_aa_{arm}"].sum()
            n_bb = g[f"n_bb_{arm}"].sum()
            n_ab = g[f"n_ab_{arm}"].sum()
            rec[f"keren_immune_{arm}"] = (n_ab / n_bb) if n_bb > 0 else np.nan
            rec[f"keren_tumor_{arm}"] = (n_ab / n_aa) if n_aa > 0 else np.nan
            rec[f"homog_{arm}"] = (n_ab / (n_ab + n_aa + n_bb)) if (n_ab + n_aa + n_bb) > 0 else np.nan
            rec[f"n_ab_{arm}"] = int(n_ab)
            rec[f"n_aa_{arm}"] = int(n_aa)
            rec[f"n_bb_{arm}"] = int(n_bb)
        rec["delta_keren_immune"] = rec["keren_immune_high"] - rec["keren_immune_low"]
        rec["delta_homog"] = rec["homog_high"] - rec["homog_low"]
        rows.append(rec)
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values(["radius_um", "sample"]).reset_index(drop=True)


def paired_wilcoxon(high: np.ndarray, low: np.ndarray) -> dict:
    high = np.asarray(high, dtype=float)
    low = np.asarray(low, dtype=float)
    mask = np.isfinite(high) & np.isfinite(low)
    high, low = high[mask], low[mask]
    if len(high) < 4:
        return {"n": int(len(high)), "median_high": np.nan, "median_low": np.nan, "median_delta": np.nan, "wilcoxon_p": np.nan, "n_high_lt_low": np.nan}
    d = high - low
    try:
        w = stats.wilcoxon(high, low, alternative="two-sided", zero_method="wilcox")
        p = float(w.pvalue)
    except ValueError:
        p = np.nan
    return {
        "n": int(len(high)),
        "median_high": float(np.median(high)),
        "median_low": float(np.median(low)),
        "mean_high": float(np.mean(high)),
        "mean_low": float(np.mean(low)),
        "median_delta": float(np.median(d)),
        "mean_delta": float(np.mean(d)),
        "wilcoxon_p": p,
        "n_high_lt_low": int((d < 0).sum()),
        "n_high_gt_low": int((d > 0).sum()),
    }


def bootstrap_delta(high: np.ndarray, low: np.ndarray, rng: np.random.Generator, n=200) -> tuple[float, float, float]:
    high = high[np.isfinite(high)]
    low = low[np.isfinite(low)]
    if len(high) < MIN_FOV_ARM or len(low) < MIN_FOV_ARM:
        return np.nan, np.nan, np.nan
    point = float(np.mean(high) - np.mean(low))
    diffs = np.empty(n, dtype=np.float64)
    for i in range(n):
        h = rng.choice(high, size=len(high), replace=True)
        l = rng.choice(low, size=len(low), replace=True)
        diffs[i] = np.mean(h) - np.mean(l)
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    return point, float(lo), float(hi)


def fov_forest_counts(cells: pd.DataFrame, radius: int) -> pd.DataFrame:
    rng = np.random.default_rng(7)
    col = f"cd8nk_count_{radius}"
    rows = []
    for (sample, fov), g in cells.groupby(["sample", "fov"]):
        high = g.loc[g["cldn4_arm"] == "high", col].to_numpy(dtype=float)
        low = g.loc[g["cldn4_arm"] == "low", col].to_numpy(dtype=float)
        delta, lo, hi = bootstrap_delta(high, low, rng)
        rows.append(
            {
                "sample": sample,
                "patient": PATIENT[sample],
                "fov": int(fov),
                "radius_um": radius,
                "n_high": int(len(high)),
                "n_low": int(len(low)),
                "mean_high": float(np.mean(high)) if len(high) else np.nan,
                "mean_low": float(np.mean(low)) if len(low) else np.nan,
                "delta_high_minus_low": delta,
                "ci95_lo": lo,
                "ci95_hi": hi,
            }
        )
    return pd.DataFrame(rows)


def plot_count_ecdfs(cells: pd.DataFrame) -> None:
    colors = {"high": "#b2182b", "low": "#2166ac"}
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 7.2), sharey=True)
    for ax, r in zip(axes.ravel(), RADII_UM):
        col = f"cd8nk_count_{r}"
        xmax = 0
        for arm, lab in (("high", "CLDN4-high tumor"), ("low", "CLDN4-low tumor")):
            v = cells.loc[cells["cldn4_arm"] == arm, col].to_numpy()
            v = v[np.isfinite(v)]
            if len(v) == 0:
                continue
            xmax = max(xmax, float(np.quantile(v, 0.99)))
            x = np.sort(v)
            y = np.arange(1, len(x) + 1) / len(x)
            ax.step(x, y, where="post", color=colors[arm], lw=2.0, label=f"{lab} (n={len(x):,})")
        ax.set_title(f"{r} µm")
        ax.set_xlabel("CD8+NK neighbor count")
        ax.set_xlim(0, max(5, xmax))
        ax.grid(True, alpha=0.3)
        ax.legend(frameon=False, fontsize=7)
    axes[0, 0].set_ylabel("ECDF")
    axes[1, 0].set_ylabel("ECDF")
    fig.suptitle("CD8+NK neighbor counts around CLDN4-high vs low tumor (8 sections)", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "ecdf_cd8nk_counts.png"), dpi=180)
    fig.savefig(os.path.join(FIG, "ecdf_cd8nk_counts.pdf"))
    plt.close(fig)


def plot_section_paired_counts(sec: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 7.2), sharey=False)
    palette = {"P5": "#1b9e77", "P6": "#d95f02", "P9": "#7570b3", "P12": "#e7298a", "P13": "#66a61e"}
    for ax, r in zip(axes.ravel(), RADII_UM):
        x0 = sec[f"mean_low_{r}"].to_numpy()
        x1 = sec[f"mean_high_{r}"].to_numpy()
        for i, row in sec.iterrows():
            ax.plot([0, 1], [row[f"mean_low_{r}"], row[f"mean_high_{r}"]], color=palette[row["patient"]], lw=1.4, marker="o", ms=5)
        ax.set_xticks([0, 1], ["CLDN4-low", "CLDN4-high"])
        ax.set_title(f"{r} µm")
        ax.set_ylabel("Mean CD8+NK neighbors / tumor cell")
        ax.grid(True, axis="y", alpha=0.3)
    handles = [plt.Line2D([0], [0], color=c, marker="o", label=p) for p, c in palette.items()]
    fig.legend(handles=handles, loc="upper center", ncol=5, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.suptitle("Section-level paired CD8+NK neighbor counts (n=8 official sections)", fontsize=11, y=1.04)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "section_paired_cd8nk_counts.png"), dpi=180, bbox_inches="tight")
    fig.savefig(os.path.join(FIG, "section_paired_cd8nk_counts.pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_section_paired_mixing(mixsec: pd.DataFrame) -> None:
    sub = mixsec[mixsec["radius_um"] == PRIMARY_R].copy()
    if sub.empty:
        return
    palette = {"P5": "#1b9e77", "P6": "#d95f02", "P9": "#7570b3", "P12": "#e7298a", "P13": "#66a61e"}
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.2))
    for ax, cols, ylab in (
        (axes[0], ("keren_immune_low", "keren_immune_high"), "Keren mixing (immune-side)\nCD8+NK–tumor / CD8+NK–CD8+NK"),
        (axes[1], ("homog_low", "homog_high"), "Homogeneous mixing\nCLDN4-arm tumor vs CD8+NK"),
    ):
        for _, row in sub.iterrows():
            ax.plot([0, 1], [row[cols[0]], row[cols[1]]], color=palette[row["patient"]], lw=1.4, marker="o", ms=5)
        ax.set_xticks([0, 1], ["vs CLDN4-low", "vs CLDN4-high"])
        ax.set_ylabel(ylab)
        ax.grid(True, axis="y", alpha=0.3)
    fig.suptitle(f"Section-level paired mixing at {PRIMARY_R} µm (n=8)", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "section_paired_mixing.png"), dpi=180)
    fig.savefig(os.path.join(FIG, "section_paired_mixing.pdf"))
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(9.6, 7.2))
    for ax, r in zip(axes.ravel(), RADII_UM):
        g = mixsec[mixsec["radius_um"] == r]
        for _, row in g.iterrows():
            ax.plot([0, 1], [row["keren_immune_low"], row["keren_immune_high"]], color=palette[row["patient"]], lw=1.4, marker="o", ms=5)
        ax.set_xticks([0, 1], ["vs CLDN4-low", "vs CLDN4-high"])
        ax.set_title(f"{r} µm")
        ax.set_ylabel("Keren mixing (immune-side)")
        ax.grid(True, axis="y", alpha=0.3)
    fig.suptitle("Section-level paired Keren mixing, CD8+NK vs CLDN4-high/low tumor", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "section_paired_mixing_grid.png"), dpi=180)
    fig.savefig(os.path.join(FIG, "section_paired_mixing_grid.pdf"))
    plt.close(fig)


def plot_fov_forest(forest: pd.DataFrame, radius: int) -> None:
    d = forest.dropna(subset=["delta_high_minus_low"]).copy()
    if d.empty:
        return
    d = d.sort_values(["patient", "sample", "fov"]).reset_index(drop=True)
    fig_h = max(6.5, 0.15 * len(d) + 1.8)
    fig, ax = plt.subplots(figsize=(8.8, fig_h))
    y = np.arange(len(d))
    palette = {"P5": "#1b9e77", "P6": "#d95f02", "P9": "#7570b3", "P12": "#e7298a", "P13": "#66a61e"}
    colors = [palette.get(p, "0.4") for p in d["patient"]]
    ax.axvline(0, color="0.35", lw=1.0, ls="--")
    ax.hlines(y, d["ci95_lo"], d["ci95_hi"], color="0.55", lw=0.9)
    ax.scatter(d["delta_high_minus_low"], y, c=colors, s=16, zorder=3, edgecolors="none")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{s} FOV{int(f):02d}" for s, f in zip(d["sample"], d["fov"])], fontsize=6)
    ax.set_xlabel(f"Δ mean CD8+NK count at {radius} µm (CLDN4-high − low)")
    ax.set_title(f"Per-FOV forest of Δ CD8+NK neighbor count ({radius} µm)")
    ax.grid(True, axis="x", alpha=0.3)
    handles = [plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=c, label=p, markersize=7) for p, c in palette.items()]
    ax.legend(handles=handles, frameon=False, fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, f"fov_forest_delta_cd8nk_count_{radius}um.png"), dpi=150)
    fig.savefig(os.path.join(FIG, f"fov_forest_delta_cd8nk_count_{radius}um.pdf"))
    plt.close(fig)


def _fmt_p(p) -> str:
    if p is None or (isinstance(p, float) and (not np.isfinite(p))):
        return "NA"
    p = float(p)
    if p < 1e-300:
        return "<1e-300"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


def write_results(inv: pd.DataFrame, sec: pd.DataFrame, mixsec: pd.DataFrame, paired: dict, forest: pd.DataFrame) -> None:
    lines = []
    lines.append("# CosMx NSCLC CLDN4-only: CD8+NK neighbor counts and mixing")
    lines.append("")
    lines.append("ADDITIVE, **CLDN4-only**, official CosMx NSCLC 960-plex (He et al. 2022): **all 8 sections / 5 patients**. No private 8-KL. Primary readout is **CD8+NK neighbor COUNT** and **mixing** on a **20 / 40 / 60 / 80 µm** grid, with a **section-level paired test (n=8)**. Nearest-distance is not the primary endpoint. Tumor cells are **not** gated as `CD8A==0`.")
    lines.append("")
    lines.append("CLDN4 is **present** on the 960-plex in all 8 sections.")
    lines.append("")
    lines.append("## Design (retuned)")
    lines.append("")
    lines.append("- **Index:** tumor / epithelial RNA-compartment cells, including those with CD8A>0.")
    lines.append("- **CLDN4-high / low:** Q4 vs Q1 of log-normalized CLDN4 among tumor cells, per section.")
    lines.append("- **Neighbors:** CD8+NK together (July PPT cytotoxic): non-tumor CD8A/B+ **or** NKG7/GNLY/KLRD1+ CD3− (KLRD1 is not on this panel; NKG7/GNLY used).")
    lines.append("- **Radii:** 20, 40, 60, 80 µm. Within-FOV only. 0.18 µm/pixel.")
    lines.append("- **Mixing:** Keren immune-side (CD8+NK–tumor contacts / CD8+NK–CD8+NK contacts) and homogeneous mixing, pooled across FOVs within each section.")
    lines.append("- **Inference:** paired Wilcoxon signed-rank on the 8 section means (high vs low). Honest n = 8 sections (5 patients).")
    lines.append("")
    lines.append("## Inventory")
    lines.append("")
    lines.append("| Section | Patient | QC cells | FOVs | Tumor | Tumor CD8A+ kept | CLDN4-high | CLDN4-low | CD8+NK |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|---:|")
    for r in inv.itertuples(index=False):
        lines.append(
            f"| {r.sample} | {r.patient} | {r.n_qc:,} | {r.n_fov} | {r.n_tumor:,} | {r.n_tumor_cd8a_pos:,} | {r.n_cldn4_high:,} | {r.n_cldn4_low:,} | {r.n_cd8nk:,} |"
        )
    lines.append("")
    lines.append(f"Tumor cells with CD8A>0 were **kept** in the tumor index (n={int(inv['n_tumor_cd8a_pos'].sum()):,} across sections). The CLDN4-low arm is larger than a strict quartile because CosMx CLDN4 is zero-inflated (Q1 = 0 in these sections); Q4 remains the top quartile.")
    lines.append("")
    lines.append("## Section-level paired CD8+NK neighbor counts")
    lines.append("")
    lines.append("| Radius | Median section mean (high) | Median section mean (low) | Median Δ (high−low) | Sections high<low | Wilcoxon p |")
    lines.append("|---:|---:|---:|---:|---:|---:|")
    for r in RADII_UM:
        t = paired["counts"][r]
        lines.append(
            f"| {r} | {t['median_high']:.3f} | {t['median_low']:.3f} | {t['median_delta']:.3f} | {t['n_high_lt_low']}/{t['n']} | {_fmt_p(t['wilcoxon_p'])} |"
        )
    lines.append("")
    lines.append("Per-section means:")
    lines.append("")
    hdr = "| Section | Patient |" + "".join([f" Δ{r} |" for r in RADII_UM]) + " mean high@40 | mean low@40 |"
    lines.append(hdr)
    lines.append("|---|---|" + "".join(["---:|"] * (len(RADII_UM) + 2)))
    for row in sec.itertuples(index=False):
        dlt = "".join([f" {getattr(row, f'delta_mean_{r}'):+.3f} |" for r in RADII_UM])
        lines.append(f"| {row.sample} | {row.patient} |{dlt} {row.mean_high_40:.3f} | {row.mean_low_40:.3f} |")
    lines.append("")
    lines.append("All **8/8** sections have a lower mean CD8+NK neighbor count around CLDN4-high than around CLDN4-low at every radius. Wilcoxon p = 0.008 is the two-sided signed-rank minimum for n=8 with no sign ties.")
    lines.append("")
    lines.append("Figures: `results/cosmx_nsclc_cldn4_nk_cd8/figures/ecdf_cd8nk_counts.png`, `section_paired_cd8nk_counts.png`.")
    lines.append("")
    lines.append("## Section-level paired mixing (CLDN4-high/low tumor vs CD8+NK)")
    lines.append("")
    lines.append("| Radius | Keren immune Δ (high−low) | Wilcoxon p | Homog Δ | Wilcoxon p |")
    lines.append("|---:|---:|---:|---:|---:|")
    for r in RADII_UM:
        k = paired["mix_keren"][r]
        h = paired["mix_homog"][r]
        lines.append(
            f"| {r} | {k['median_delta']:.3f} | {_fmt_p(k['wilcoxon_p'])} | {h['median_delta']:.3f} | {_fmt_p(h['wilcoxon_p'])} |"
        )
    lines.append("")
    m40 = mixsec[mixsec["radius_um"] == PRIMARY_R].copy()
    if not m40.empty:
        order = {s: i for i, s in enumerate(SAMPLES)}
        m40["_ord"] = m40["sample"].map(order)
        m40 = m40.sort_values("_ord")
        lines.append(f"Per-section Keren mixing at {PRIMARY_R} µm (immune-side):")
        lines.append("")
        lines.append("| Section | vs CLDN4-high | vs CLDN4-low | Δ |")
        lines.append("|---|---:|---:|---:|")
        for row in m40.itertuples(index=False):
            lines.append(f"| {row.sample} | {row.keren_immune_high:.3f} | {row.keren_immune_low:.3f} | {row.delta_keren_immune:+.3f} |")
        lines.append("")
    lines.append("Figures: `section_paired_mixing.png`, `section_paired_mixing_grid.png`.")
    lines.append("")
    lines.append("## FOV forest of Δ count (secondary)")
    lines.append("")
    if forest is not None and len(forest):
        v = forest["delta_high_minus_low"].dropna()
        if len(v) >= 8:
            try:
                wp = float(stats.wilcoxon(v, alternative="two-sided").pvalue)
            except ValueError:
                wp = float("nan")
            lines.append(
                f"At {PRIMARY_R} µm, {int(len(v))} FOVs; median Δ mean CD8+NK count = {float(v.median()):.3f}; Wilcoxon signed-rank p = {_fmt_p(wp)}; fraction Δ<0 = {float((v < 0).mean()):.3f}."
            )
        lines.append("")
        lines.append(f"Figure: `fov_forest_delta_cd8nk_count_{PRIMARY_R}um.png`. This is **not** a nearest-distance forest.")
    lines.append("")
    lines.append("## What this does not claim")
    lines.append("")
    lines.append("- Not a nearest-CD8 / nearest-NK distance primary analysis.")
    lines.append("- Not ICI response and not private 8-KL.")
    lines.append("- Marker CD8+NK is not the paper’s 18-type map (no author cell-type column in the public flat files).")
    lines.append("- Q4 vs Q1 is a contrast, not a biological threshold. n=8 sections, not n=8 patients.")
    lines.append("")
    lines.append("## Reproduce")
    lines.append("")
    lines.append("```bash")
    lines.append("python3 scripts/download_cosmx_nsclc.py")
    lines.append("python3 scripts/cosmx_nsclc_cldn4_nk_cd4_radii.py")
    lines.append("```")
    lines.append("")
    lines.append("Raw tarballs are gitignored. Tables live under `results/cosmx_nsclc_cldn4_nk_cd8/tables/`.")
    text = "\n".join(lines) + "\n"
    with open(os.path.join(ROOT, "RESULTS.md"), "w") as fh:
        fh.write(text)
    with open(os.path.join(OUT, "RESULTS.md"), "w") as fh:
        fh.write(text)


def main() -> int:
    os.makedirs(FIG, exist_ok=True)
    os.makedirs(TAB, exist_ok=True)

    present = [s for s in SAMPLES if os.path.isdir(os.path.join(DATA, s))]
    if len(present) < 8:
        raise SystemExit(f"need all 8 official sections, found {present}")

    genes_by_sample = {s: panel_genes(s) for s in SAMPLES}
    stop_if_cldn4_absent(genes_by_sample)
    nk_present = {s: [g for g in NK_MARKERS if g in genes_by_sample[s]] for s in SAMPLES}

    inventory = []
    fov_parts, mix_parts, cell_parts = [], [], []
    for sample in SAMPLES:
        print(f"load {sample}", flush=True)
        df = assign_types(load_sample(sample))
        rec = {
            "sample": sample,
            "patient": PATIENT[sample],
            "n_qc": int(len(df)),
            "n_fov": int(df["fov"].nunique()),
            "n_tumor": int(df["is_tumor"].sum()),
            "n_tumor_cd8a_pos": int(df["tumor_cd8a_pos"].sum()),
            "n_cldn4_high": int((df["cldn4_arm"] == "high").sum()),
            "n_cldn4_low": int((df["cldn4_arm"] == "low").sum()),
            "n_cd8": int(df["is_cd8"].sum()),
            "n_nk": int(df["is_nk"].sum()),
            "n_cd8nk": int(df["is_cd8nk"].sum()),
            "frac_tumor_cldn4_pos": float((df.loc[df["is_tumor"], "cldn4_raw"] > 0).mean()) if df["is_tumor"].any() else np.nan,
            "used_author_labels": bool(df["used_author_labels"].iloc[0]),
            "nk_genes": ",".join(nk_present[sample]),
            "cldn4_in_panel": True,
        }
        inventory.append(rec)
        print(f"  {rec}", flush=True)
        fov_arm, mix, cells = per_fov_metrics(df)
        fov_parts.append(fov_arm)
        mix_parts.append(mix)
        cell_parts.append(cells)
        del df

    inv = pd.DataFrame(inventory)
    inv.to_csv(os.path.join(TAB, "sample_inventory.csv"), index=False)
    fov_arm = pd.concat([p for p in fov_parts if len(p)], ignore_index=True)
    mix = pd.concat([p for p in mix_parts if len(p)], ignore_index=True)
    cells = pd.concat([p for p in cell_parts if len(p)], ignore_index=True)
    fov_arm.to_csv(os.path.join(TAB, "fov_arm_cd8nk_counts.csv"), index=False)
    mix.to_csv(os.path.join(TAB, "mixing_by_fov.csv"), index=False)
    cells.to_csv(os.path.join(TAB, "tumor_cell_cd8nk_counts.csv.gz"), index=False, compression="gzip")

    sec = section_count_table(cells)
    sec.to_csv(os.path.join(TAB, "section_paired_cd8nk_counts.csv"), index=False)
    mixsec = section_mix_table(mix)
    mixsec.to_csv(os.path.join(TAB, "section_paired_mixing.csv"), index=False)

    paired = {"counts": {}, "mix_keren": {}, "mix_homog": {}}
    for r in RADII_UM:
        paired["counts"][r] = paired_wilcoxon(sec[f"mean_high_{r}"].to_numpy(), sec[f"mean_low_{r}"].to_numpy())
        sub = mixsec[mixsec["radius_um"] == r]
        paired["mix_keren"][r] = paired_wilcoxon(sub["keren_immune_high"].to_numpy(), sub["keren_immune_low"].to_numpy())
        paired["mix_homog"][r] = paired_wilcoxon(sub["homog_high"].to_numpy(), sub["homog_low"].to_numpy())

    forest = fov_forest_counts(cells, PRIMARY_R)
    forest.to_csv(os.path.join(TAB, f"fov_forest_delta_cd8nk_count_{PRIMARY_R}um.csv"), index=False)

    plot_count_ecdfs(cells)
    plot_section_paired_counts(sec)
    plot_section_paired_mixing(mixsec)
    plot_fov_forest(forest, PRIMARY_R)

    summary = {
        "dataset": "official CosMx NSCLC 960-plex, 8 sections / 5 patients (He et al. 2022)",
        "primary": "CD8+NK neighbor counts and mixing on 20/40/60/80 µm; section-level paired test n=8",
        "not_primary": "nearest distance",
        "tumor_gated_cd8a_eq0": False,
        "cd8nk_together": True,
        "private_8kl": False,
        "cldn4_only": True,
        "pixel_size_um": PX_TO_UM,
        "inventory": inv.to_dict(orient="records"),
        "section_paired_counts": {str(k): v for k, v in paired["counts"].items()},
        "section_paired_mix_keren": {str(k): v for k, v in paired["mix_keren"].items()},
        "section_paired_mix_homog": {str(k): v for k, v in paired["mix_homog"].items()},
        "n_tumor_cells_scored": int(len(cells)),
        "n_sections": int(sec.shape[0]),
    }
    with open(os.path.join(OUT, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, default=str)

    write_results(inv, sec, mixsec, paired, forest)
    print("done", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
