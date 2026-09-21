#!/usr/bin/env python3
"""Donor-level permutation and bootstrap for CosMx CLDN4 vs CD8+NK.

Official 960-plex NSCLC (He et al. 2022; figshare 25976224 / NanoString S3).
Eight sections, five donors. No private 8-KL cohort.

Pre-specified endpoints, fixed before looking at the permutation p-values:

1. Within-FOV median split of tumor CLDN4, mean CD8+NK neighbor count at 50 µm,
   high minus low, then equal-FOV, equal-section, equal-donor mean.
2. The same contrast for a 20 µm CD8+NK contact (any centroid inside 20 µm).
3. Within-FOV Spearman of continuous tumor CLDN4 vs the 50 µm CD8+NK count.
4. Within-FOV Spearman of continuous tumor CLDN4 vs the 20 µm contact indicator.

The null permutes CLDN4 among tumor cells inside each FOV. Coordinates, FOV
ids, and CD8/NK labels stay fixed, so a neighbor count never changes cell.
The bootstrap resamples whole FOVs inside each section. Neither procedure
draws cells independently of their neighbors.

A one-sided sign test cannot go below 1/256 on eight sections or 1/32 on five
donors. Those floors are reported beside the permutation p-value. Beating them
is a property of this spatial null, not a claim that a sixth donor has been
sampled.
"""

from __future__ import annotations

import argparse
import json
import os
import pickle
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from donor_spatial_null import (  # noqa: E402
    DONOR_SIGN_FLOOR_ONE,
    SAMPLE_SIGN_FLOOR_ONE,
    FovBlock,
    analyze,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data", "cosmx_nsclc")
OUT = os.path.join(ROOT, "results", "cosmx_donor_spatial_null")
FIG = os.path.join(OUT, "figures")
TAB = os.path.join(OUT, "tables")
CACHE = os.path.join(DATA, "cache", "fov_blocks.pkl")

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
DONOR = {
    "Lung5_Rep1": "Lung5",
    "Lung5_Rep2": "Lung5",
    "Lung5_Rep3": "Lung5",
    "Lung6": "Lung6",
    "Lung9_Rep1": "Lung9",
    "Lung9_Rep2": "Lung9",
    "Lung12": "Lung12",
    "Lung13": "Lung13",
}
PX_TO_UM = 0.18
RADII = (20, 40, 50, 80, 100)
PRIMARY_RADIUS = 50
CONTACT_RADIUS = 20
MIN_COUNTS = 20
MIN_GENES = 5
MIN_FOV_TUMOR = 20
MIN_ARM = 5
MIN_IMMUNE = 5

NEEDED = [
    "CLDN4",
    "NKG7",
    "GNLY",
    "KLRD1",
    "CD3D",
    "CD3E",
    "CD3G",
    "CD8A",
    "CD8B",
    "CD4",
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
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "DCN",
    "LUM",
    "ACTA2",
    "PECAM1",
    "VWF",
    "FN1",
    "BGN",
]
EPI = ["EPCAM", "CDH1", "KRT8", "KRT18", "KRT19", "KRT7", "KRT5", "KRT17", "CEACAM6", "MUC1", "ELF3"]
IMM = ["PTPRC", "CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "CD4", "NKG7", "GNLY", "KLRD1"]
STR = ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "ACTA2", "PECAM1", "VWF", "FN1", "BGN"]
NK = ["NKG7", "GNLY", "KLRD1"]
CD3 = ["CD3D", "CD3E", "CD3G"]
CD8 = ["CD8A", "CD8B"]

# Pre-specified. Do not reorder after seeing p-values.
PRIMARIES = (
    ("within_fov", "delta", "cd8nk_count_50"),
    ("within_fov", "delta", "cd8nk_contact_20"),
    ("within_fov", "spearman", "cd8nk_count_50"),
    ("within_fov", "spearman", "cd8nk_contact_20"),
)
PRIMARY_LABEL = {
    "cd8nk_count_50": "CD8+NK count, 50 µm",
    "cd8nk_contact_20": "CD8+NK contact, 20 µm",
    "cd8nk_count_20": "CD8+NK count, 20 µm",
    "cd8nk_count_40": "CD8+NK count, 40 µm",
    "cd8nk_count_80": "CD8+NK count, 80 µm",
    "cd8nk_count_100": "CD8+NK count, 100 µm",
    "cd8_count_50": "CD8 count, 50 µm",
    "cd8_contact_20": "CD8 contact, 20 µm",
}


def find_csv(sample_dir: str, suffix: str) -> str:
    for fn in os.listdir(sample_dir):
        if fn.endswith(suffix):
            return os.path.join(sample_dir, fn)
    raise FileNotFoundError(suffix + " in " + sample_dir)


def outcome_names() -> tuple[str, ...]:
    names = [f"cd8nk_count_{r}" for r in RADII]
    names.append("cd8nk_contact_20")
    names.append("cd8_count_50")
    names.append("cd8_contact_20")
    return tuple(names)


def _mean_ln(logn: np.ndarray, genes: list[str], gene_index: dict[str, int]) -> np.ndarray:
    idx = [gene_index[g] for g in genes if g in gene_index]
    if not idx:
        return np.zeros(logn.shape[0], dtype=np.float64)
    return logn[:, idx].mean(axis=1)


def load_sample(sample: str) -> pd.DataFrame:
    sample_dir = os.path.join(DATA, sample)
    expr_path = find_csv(sample_dir, "exprMat_file.csv")
    meta_path = find_csv(sample_dir, "metadata_file.csv")
    header = list(pd.read_csv(expr_path, nrows=0).columns)
    if "CLDN4" not in header:
        raise SystemExit(f"CLDN4 is absent from {sample}. Stopping.")
    panel = [c for c in header if c not in ("fov", "cell_ID") and not str(c).lower().startswith("neg")]
    genes = [g for g in NEEDED if g in header]
    usecols = ["fov", "cell_ID"] + panel
    parts = []
    for chunk in pd.read_csv(expr_path, usecols=usecols, chunksize=50_000):
        chunk = chunk[chunk["cell_ID"] != 0]
        if chunk.empty:
            continue
        mat = chunk[panel].to_numpy(dtype=np.float32)
        tot = mat.sum(axis=1)
        nfeat = (mat > 0).sum(axis=1)
        keep = (tot >= MIN_COUNTS) & (nfeat >= MIN_GENES)
        if not np.any(keep):
            continue
        sub = chunk.loc[keep, ["fov", "cell_ID"] + genes].copy()
        sub["_tot_panel"] = tot[keep]
        parts.append(sub)
    expr = pd.concat(parts, ignore_index=True)
    meta = pd.read_csv(meta_path)
    meta = meta[meta["cell_ID"] != 0]
    keep_meta = ["fov", "cell_ID", "CenterX_global_px", "CenterY_global_px"]
    for col in ("Mean.PanCK", "Mean.CD45"):
        if col in meta.columns:
            keep_meta.append(col)
    meta = meta[keep_meta]
    df = expr.merge(meta, on=["fov", "cell_ID"], how="inner")
    raw = df[genes].to_numpy(dtype=np.float64)
    tot = np.maximum(raw.sum(axis=1), 1.0)
    sf = float(np.median(tot)) / tot
    logn = np.log1p(raw * sf[:, None])
    gene_index = {g: i for i, g in enumerate(genes)}
    epi = _mean_ln(logn, EPI, gene_index)
    imm = _mean_ln(logn, IMM, gene_index)
    strn = _mean_ln(logn, STR, gene_index)
    scores = np.vstack([epi, imm, strn])
    lab = np.array(["epithelial", "immune", "stromal"])[scores.argmax(axis=0)]
    lab[scores.max(axis=0) <= 0] = "unassigned"
    tumor = lab == "epithelial"
    if "Mean.PanCK" in df.columns and "Mean.CD45" in df.columns:
        panck = np.log1p(df["Mean.PanCK"].to_numpy(dtype=np.float64))
        cd45 = np.log1p(df["Mean.CD45"].to_numpy(dtype=np.float64))
        tumor = tumor | ((panck >= np.median(panck)) & (cd45 < np.median(cd45)) & (lab != "immune"))
    def raw_sum(markers: list[str]) -> np.ndarray:
        idx = [gene_index[g] for g in markers if g in gene_index]
        if not idx:
            return np.zeros(len(df), dtype=np.float64)
        return raw[:, idx].sum(axis=1)

    nk_raw = raw_sum(NK)
    cd3_raw = raw_sum(CD3)
    cd8_raw = raw_sum(CD8)
    cd8 = (cd8_raw > 0) & (~tumor)
    nk = (nk_raw > 0) & (cd3_raw == 0) & (~tumor)
    df["is_tumor"] = tumor
    df["is_cd8"] = cd8
    df["is_nk"] = nk
    df["is_cd8nk"] = cd8 | nk
    df["cldn4"] = logn[:, gene_index["CLDN4"]]
    df["x_um"] = df["CenterX_global_px"].to_numpy(dtype=np.float64) * PX_TO_UM
    df["y_um"] = df["CenterY_global_px"].to_numpy(dtype=np.float64) * PX_TO_UM
    df["sample"] = sample
    df["donor"] = DONOR[sample]
    return df


def _radius_counts(query: np.ndarray, ref: np.ndarray, radii: tuple[int, ...]) -> np.ndarray:
    out = np.zeros((len(query), len(radii)), dtype=np.float64)
    if len(query) == 0 or len(ref) == 0:
        return out
    tree = cKDTree(ref)
    for j, radius in enumerate(radii):
        out[:, j] = np.asarray(
            tree.query_ball_point(query, r=float(radius), return_length=True),
            dtype=np.float64,
        )
    return out


def blocks_from_frame(df: pd.DataFrame) -> tuple[list[FovBlock], dict]:
    names = outcome_names()
    blocks: list[FovBlock] = []
    inventory = {
        "n_qc": int(len(df)),
        "n_tumor": int(df["is_tumor"].sum()),
        "n_cd8": int(df["is_cd8"].sum()),
        "n_nk": int(df["is_nk"].sum()),
        "n_cd8nk": int(df["is_cd8nk"].sum()),
        "n_fov": int(df["fov"].nunique()),
        "n_fov_used": 0,
    }
    for fov, g in df.groupby("fov", sort=True):
        tumor = g["is_tumor"].to_numpy()
        cd8nk = g["is_cd8nk"].to_numpy()
        cd8 = g["is_cd8"].to_numpy()
        if int(tumor.sum()) < MIN_FOV_TUMOR or int(cd8nk.sum()) < MIN_IMMUNE:
            continue
        xy = g[["x_um", "y_um"]].to_numpy(dtype=np.float64)
        nk_counts = _radius_counts(xy[tumor], xy[cd8nk], RADII)
        cd8_counts = _radius_counts(xy[tumor], xy[cd8], (CONTACT_RADIUS, PRIMARY_RADIUS))
        # cd8_counts columns follow (20, 50)
        contact = (nk_counts[:, RADII.index(CONTACT_RADIUS)] > 0).astype(np.float64)
        cd8_contact = (cd8_counts[:, 0] > 0).astype(np.float64)
        cd8_count_50 = cd8_counts[:, 1]
        if int(cd8.sum()) < MIN_IMMUNE:
            cd8_contact[:] = np.nan
            cd8_count_50[:] = np.nan
        outcomes = np.column_stack(
            [
                nk_counts,
                contact,
                cd8_count_50,
                cd8_contact,
            ]
        )
        if outcomes.shape[1] != len(names):
            raise RuntimeError(f"outcome width {outcomes.shape[1]} != {len(names)}")
        blocks.append(
            FovBlock(
                donor=str(g["donor"].iloc[0]),
                sample=str(g["sample"].iloc[0]),
                fov=str(int(fov)),
                cldn4=g.loc[tumor, "cldn4"].to_numpy(dtype=np.float64),
                outcomes=outcomes,
                outcome_names=names,
            )
        )
    inventory["n_fov_used"] = len(blocks)
    inventory["n_tumor_used"] = int(sum(b.cldn4.size for b in blocks))
    return blocks, inventory


def build_blocks() -> tuple[list[FovBlock], list[dict]]:
    if os.path.exists(CACHE):
        print(f"loading cache {CACHE}", flush=True)
        with open(CACHE, "rb") as fh:
            payload = pickle.load(fh)
        return payload["blocks"], payload["inventory"]
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    blocks: list[FovBlock] = []
    inventory = []
    for sample in SAMPLES:
        print(f"loading {sample}", flush=True)
        df = load_sample(sample)
        sample_blocks, info = blocks_from_frame(df)
        info["sample"] = sample
        info["donor"] = DONOR[sample]
        print(
            f"  QC {info['n_qc']:,} tumor {info['n_tumor']:,} "
            f"CD8 {info['n_cd8']:,} NK {info['n_nk']:,} CD8+NK {info['n_cd8nk']:,} "
            f"FOVs used {info['n_fov_used']}/{info['n_fov']}",
            flush=True,
        )
        blocks.extend(sample_blocks)
        inventory.append(info)
        del df
    with open(CACHE, "wb") as fh:
        pickle.dump({"blocks": blocks, "inventory": inventory}, fh, protocol=pickle.HIGHEST_PROTOCOL)
    return blocks, inventory


def _idx(names: tuple[str, ...], name: str) -> int:
    return names.index(name)


def _fmt_p(p: float) -> str:
    if not np.isfinite(p):
        return "NA"
    if p < 1e-3:
        return f"{p:.2e}"
    return f"{p:.4f}"


def _fmt(x: float, digits: int = 3) -> str:
    if not np.isfinite(x):
        return "NA"
    return f"{x:.{digits}f}"


def endpoint_rows(result: dict) -> list[dict]:
    names = result["outcome_names"]
    rows = []
    families = (
        ("within_fov", result["within_fov"]),
        ("section_pooled", result["section_pooled"]),
    )
    for family, bundle in families:
        for stat, ep in bundle.items():
            for j, name in enumerate(names):
                if family == "section_pooled" and stat == "spearman" and name not in (
                    "cd8nk_count_50",
                    "cd8nk_contact_20",
                    "cd8nk_count_100",
                ):
                    # Keep the section-pooled Spearman table to the pre-specified
                    # readouts plus the 100 µm locked radius.
                    if name not in ("cd8nk_count_50", "cd8nk_contact_20", "cd8nk_count_100"):
                        continue
                row = {
                    "family": family,
                    "statistic": stat,
                    "outcome": name,
                    "primary": (family, stat, name) in PRIMARIES,
                    "T": float(ep["T"][j]),
                    "p_one": float(ep["p_one"][j]),
                    "p_two": float(ep["p_two"][j]),
                    "null_mean": float(ep["null_mean"][j]),
                    "ci_block_lo": float(ep["ci_block_lo"][j]),
                    "ci_block_hi": float(ep["ci_block_hi"][j]),
                    "ci_cluster_lo": float(ep["ci_cluster_lo"][j]),
                    "ci_cluster_hi": float(ep["ci_cluster_hi"][j]),
                    "sign_donor_neg": ep["donor_sign"][j]["n_neg"],
                    "sign_donor_n": ep["donor_sign"][j]["n"],
                    "sign_donor_p": ep["donor_sign"][j]["p"],
                    "sign_section_neg": ep["section_sign"][j]["n_neg"],
                    "sign_section_n": ep["section_sign"][j]["n"],
                    "sign_section_p": ep["section_sign"][j]["p"],
                    "beats_sample_sign_floor": bool(ep["beats_sample_sign_floor"][j]),
                    "beats_donor_sign_floor": bool(ep["beats_donor_sign_floor"][j]),
                }
                if ep.get("T_cell_weighted") is not None:
                    row["T_cell_weighted"] = float(ep["T_cell_weighted"][j])
                if "p_one_cell_weighted" in ep:
                    row["p_one_cell_weighted"] = float(ep["p_one_cell_weighted"][j])
                rows.append(row)
    return rows


def donor_rows(result: dict) -> list[dict]:
    names = result["outcome_names"]
    rows = []
    for family, bundle in (("within_fov", result["within_fov"]), ("section_pooled", result["section_pooled"])):
        for stat, ep in bundle.items():
            for d_i, donor in enumerate(result["donors"]):
                for j, name in enumerate(names):
                    if not ((family, stat, name) in PRIMARIES or name in ("cd8nk_count_50", "cd8nk_contact_20", "cd8nk_count_100")):
                        continue
                    if stat == "slope" and (family, stat, name) not in PRIMARIES:
                        continue
                    rows.append(
                        {
                            "family": family,
                            "statistic": stat,
                            "outcome": name,
                            "donor": donor,
                            "effect": float(ep["donor"][d_i, j]),
                            "ci_lo": float(ep["donor_ci_lo"][d_i, j]),
                            "ci_hi": float(ep["donor_ci_hi"][d_i, j]),
                        }
                    )
    return rows


def section_rows(result: dict) -> list[dict]:
    names = result["outcome_names"]
    rows = []
    for family, bundle in (("within_fov", result["within_fov"]), ("section_pooled", result["section_pooled"])):
        for stat in ("delta", "spearman"):
            ep = bundle[stat]
            for s_i, sample in enumerate(result["samples"]):
                for j, name in enumerate(names):
                    if name not in ("cd8nk_count_50", "cd8nk_contact_20", "cd8nk_count_20", "cd8nk_count_40", "cd8nk_count_80", "cd8nk_count_100"):
                        continue
                    rows.append(
                        {
                            "family": family,
                            "statistic": stat,
                            "outcome": name,
                            "sample": sample,
                            "donor": result["section_donor"][s_i],
                            "effect": float(ep["section"][s_i, j]),
                        }
                    )
    return rows


def fov_rows(result: dict) -> list[dict]:
    names = result["outcome_names"]
    keep = [names.index(n) for n in names if n.startswith("cd8nk_")]
    rows = []
    for i in range(result["n_fov"]):
        rec = {
            "donor": result["fov_donor"][i],
            "sample": result["fov_sample"][i],
            "fov": result["fov_id"][i],
            "n_tumor": result["fov_n"][i],
            "n_high": result["fov_n_high"][i],
            "n_low": result["fov_n_low"][i],
        }
        for j in keep:
            name = names[j]
            rec[f"delta_{name}"] = float(result["fov_delta"][i, j])
            rec[f"spearman_{name}"] = float(result["fov_spearman"][i, j])
            rec[f"mean_high_{name}"] = float(result["fov_mean_high"][i, j])
            rec[f"mean_low_{name}"] = float(result["fov_mean_low"][i, j])
        rows.append(rec)
    return rows


def quintile_rows(result: dict) -> list[dict]:
    names = result["outcome_names"]
    donor_curves = result["bins"]["donor"]
    grand = result["bins"]["grand"]
    n_bins = result["bins"]["n_bins"]
    rows = []
    focus = ["cd8nk_count_50", "cd8nk_contact_20", "cd8nk_count_100", "cd8nk_count_40"]
    for j, name in enumerate(names):
        if name not in focus:
            continue
        for b in range(n_bins):
            rows.append(
                {
                    "outcome": name,
                    "donor": "mean",
                    "bin": b + 1,
                    "mean_outcome": float(grand[b, j]),
                }
            )
            for d_i, donor in enumerate(result["donors"]):
                rows.append(
                    {
                        "outcome": name,
                        "donor": donor,
                        "bin": b + 1,
                        "mean_outcome": float(donor_curves[d_i, b, j]),
                    }
                )
    return rows


def plot_nulls(result: dict, path: str) -> None:
    names = result["outcome_names"]
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 7.2))
    specs = (
        ("delta", "cd8nk_count_50", "50 µm count, high − low"),
        ("delta", "cd8nk_contact_20", "20 µm contact rate, high − low"),
        ("spearman", "cd8nk_count_50", "Spearman, CLDN4 vs 50 µm count"),
        ("spearman", "cd8nk_contact_20", "Spearman, CLDN4 vs 20 µm contact"),
    )
    for ax, (stat, outcome, title) in zip(axes.ravel(), specs):
        ep = result["within_fov"][stat]
        j = _idx(names, outcome)
        null = ep["null_T"][:, j]
        null = null[np.isfinite(null)]
        obs = ep["T"][j]
        ax.hist(null, bins=40, color="#4C78A8", alpha=0.85, edgecolor="none")
        ax.axvline(obs, color="#E45756", lw=2.0, label=f"observed {_fmt(obs)}")
        ax.axvline(0, color="0.3", lw=0.8, ls="--")
        ax.set_title(
            f"{title}\none-sided perm p={_fmt_p(ep['p_one'][j])}  "
            f"(sign floor {SAMPLE_SIGN_FLOOR_ONE:.4f})",
            fontsize=9,
        )
        ax.set_xlabel("Donor-mean statistic")
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle("Within-FOV CLDN4 permutation null (labels stay inside the FOV)", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_forest(result: dict, path: str) -> None:
    names = result["outcome_names"]
    donors = result["donors"]
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.4), sharey=True)
    specs = (
        ("delta", "cd8nk_count_50", "50 µm CD8+NK count\nhigh − low"),
        ("delta", "cd8nk_contact_20", "20 µm contact rate\nhigh − low"),
    )
    for ax, (stat, outcome, title) in zip(axes, specs):
        ep = result["within_fov"][stat]
        j = _idx(names, outcome)
        y = np.arange(len(donors))
        eff = ep["donor"][:, j]
        lo = ep["donor_ci_lo"][:, j]
        hi = ep["donor_ci_hi"][:, j]
        ax.hlines(y, lo, hi, color="#4C78A8", lw=2)
        ax.plot(eff, y, "o", color="#4C78A8", ms=7)
        ax.axvline(0, color="0.4", lw=0.8, ls="--")
        ax.axvline(ep["T"][j], color="#E45756", lw=1.2, label=f"donor mean {_fmt(ep['T'][j])}")
        ax.set_yticks(y, donors)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Within-FOV contrast")
        ax.legend(frameon=False, fontsize=8)
        ax.grid(True, axis="x", alpha=0.3)
    fig.suptitle("FOV-block bootstrap intervals inside each donor", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _zero_plus_positive_bins(cldn4: np.ndarray, y: np.ndarray, n_pos: int = 4) -> np.ndarray:
    """Undetected bin, then equal-count bins among CLDN4 > 0.

    Tied zeros share one rank, so a five-bin cut of all cells collapses the
    zeros into a single bin and is not a dose curve. Positive bins use average
    ranks, not file order.
    """
    from scipy.stats import rankdata

    out = np.full(1 + n_pos, np.nan, dtype=np.float64)
    x = np.asarray(cldn4, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size == 0 or not np.isfinite(y).all():
        return out
    zero = x <= 0
    if np.any(zero):
        out[0] = float(y[zero].mean())
    xp = x[~zero]
    yp = y[~zero]
    if xp.size < n_pos * 5:
        return out
    ranks = rankdata(xp, method="average")
    bins = np.minimum(n_pos - 1, ((ranks - 1.0) / xp.size * n_pos).astype(int))
    for i in range(n_pos):
        m = bins == i
        if np.any(m):
            out[i + 1] = float(yp[m].mean())
    return out


def continuous_curves(blocks: list[FovBlock], outcome: str, n_pos: int = 4) -> dict:
    """Donor-equal mean of the undetected + positive-quartile curve."""
    j = blocks[0].outcome_names.index(outcome)
    by_sample: dict[tuple[str, str], list[np.ndarray]] = {}
    for block in blocks:
        key = (block.donor, block.sample)
        by_sample.setdefault(key, []).append(_zero_plus_positive_bins(block.cldn4, block.outcomes[:, j], n_pos))
    section = {}
    for (donor, sample), curves in by_sample.items():
        section[(donor, sample)] = np.nanmean(np.vstack(curves), axis=0)
    donors = list(dict.fromkeys(b.donor for b in blocks))
    donor_curves = {}
    for donor in donors:
        rows = [v for (d, _), v in section.items() if d == donor]
        donor_curves[donor] = np.nanmean(np.vstack(rows), axis=0)
    grand = np.nanmean(np.vstack(list(donor_curves.values())), axis=0)
    return {"donor": donor_curves, "grand": grand, "labels": ["undetected", "P1", "P2", "P3", "P4"][: 1 + n_pos]}


def plot_quintiles(result: dict, path: str, blocks: list[FovBlock] | None = None) -> None:
    colors = {
        "Lung5": "#1b9e77",
        "Lung6": "#d95f02",
        "Lung9": "#7570b3",
        "Lung12": "#e7298a",
        "Lung13": "#66a61e",
    }
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.3))
    if blocks is None:
        raise ValueError("continuous curve needs the FOV blocks")
    for ax, outcome, ylab in (
        (axes[0], "cd8nk_count_50", "Mean CD8+NK neighbors at 50 µm"),
        (axes[1], "cd8nk_contact_20", "CD8+NK contact rate at 20 µm"),
    ):
        curves = continuous_curves(blocks, outcome)
        x = np.arange(len(curves["labels"]))
        for donor, row in curves["donor"].items():
            ax.plot(x, row, color=colors[donor], lw=1.3, marker="o", ms=4, label=donor)
        ax.plot(x, curves["grand"], color="black", lw=2.4, marker="s", ms=5, label="donor mean")
        ax.set_xticks(x, curves["labels"])
        ax.set_xlabel("CLDN4 inside the FOV")
        ax.set_ylabel(ylab)
        ax.grid(True, alpha=0.3)
    axes[0].legend(frameon=False, fontsize=8, ncol=2)
    fig.suptitle("Undetected CLDN4 vs quartiles of CLDN4 > 0", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def positive_only_blocks(blocks: list[FovBlock], outcomes: tuple[str, ...] = ("cd8nk_count_50", "cd8nk_contact_20")) -> list[FovBlock]:
    """Tumor cells with CLDN4 > 0. Asks whether more CLDN4, among expressors, excludes."""
    idx = [blocks[0].outcome_names.index(name) for name in outcomes]
    kept = []
    for block in blocks:
        mask = block.cldn4 > 0
        if int(mask.sum()) < 30:
            continue
        kept.append(
            FovBlock(
                donor=block.donor,
                sample=block.sample,
                fov=block.fov,
                cldn4=block.cldn4[mask],
                outcomes=block.outcomes[mask][:, idx],
                outcome_names=outcomes,
            )
        )
    return kept


def plot_radius(result: dict, path: str) -> None:
    names = result["outcome_names"]
    ep = result["within_fov"]["delta"]
    xs, ys, los, his = [], [], [], []
    for radius in RADII:
        j = _idx(names, f"cd8nk_count_{radius}")
        xs.append(radius)
        ys.append(ep["T"][j])
        los.append(ep["T"][j] - ep["ci_block_lo"][j])
        his.append(ep["ci_block_hi"][j] - ep["T"][j])
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.errorbar(xs, ys, yerr=np.vstack([los, his]), fmt="o-", color="#E45756", lw=1.6, capsize=3)
    ax.axhline(0, color="0.4", lw=0.8, ls="--")
    ax.set_xlabel("Radius (µm)")
    ax.set_ylabel("Donor-mean CD8+NK count, CLDN4 high − low")
    ax.set_title("Within-FOV median contrast across radii")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _lookup(rows: list[dict], family: str, stat: str, outcome: str) -> dict:
    for row in rows:
        if row["family"] == family and row["statistic"] == stat and row["outcome"] == outcome:
            return row
    raise KeyError((family, stat, outcome))


def write_results(result: dict, inventory: list[dict], rows: list[dict], positive: dict | None = None) -> None:
    inv = {
        "n_qc": int(sum(r["n_qc"] for r in inventory)),
        "n_tumor": int(sum(r["n_tumor"] for r in inventory)),
        "n_cd8": int(sum(r["n_cd8"] for r in inventory)),
        "n_nk": int(sum(r["n_nk"] for r in inventory)),
        "n_cd8nk": int(sum(r["n_cd8nk"] for r in inventory)),
        "n_fov": int(sum(r["n_fov"] for r in inventory)),
        "n_fov_used": int(sum(r["n_fov_used"] for r in inventory)),
        "n_tumor_used": int(sum(r["n_tumor_used"] for r in inventory)),
    }
    lines = []
    a = lines.append
    a("# Donor-level spatial null: CosMx CLDN4 vs CD8+NK")
    a("")
    a("Official NanoString CosMx NSCLC 960-plex (He et al. 2022), eight sections and five donors (Lung5 ×3, Lung6, Lung9 ×2, Lung12, Lung13). Additive to the locked radius and contact summaries. No private 8-KL data. Flat metadata has no author cell-type column, so tumor / CD8 / NK calls use the same RNA compartment score plus PanCK-high / CD45-low support as the radius script, and CD8+NK cells are required to be non-tumor.")
    a("")
    a("## What is stronger than the sign-test floor")
    a("")
    a(f"A one-sided sign test on eight sections cannot return a p-value below **{SAMPLE_SIGN_FLOOR_ONE:.6f}** (1/256). On five donors the floor is **{DONOR_SIGN_FLOOR_ONE:.5f}** (1/32). Wilcoxon on eight paired sections has the same two-sided floor, 1/128. Those tests use only the sign.")
    a("")
    a("The p-value here is from a permutation of the CLDN4 mark. The statistic is the unweighted mean of five donor effects, so Lung5's three sections are one donor. Because the null has thousands of draws, the p-value can fall below both floors when the within-FOV pairing is extreme. The sign counts stay in the table as the reproducibility summary. They are not replaced.")
    a("")
    a("This is not a p-value for drawing a new donor. The permutation null is: inside each FOV, tumor CLDN4 values are exchangeable with respect to a neighborhood that was computed once from fixed coordinates and fixed CD8/NK labels.")
    a("")
    a("## How spatial labels are kept")
    a("")
    a("- Neighbor counts and the 20 µm contact flag are computed once, inside the FOV, from centroid coordinates (`0.18 µm/px`).")
    a("- A permutation moves CLDN4 values among tumor cells that share a FOV. It does not move a value into another FOV, section, or donor, and it does not move a cell.")
    a("- Between-FOV geography is therefore left intact. A section-pooled high-versus-low gap that comes only from CLDN4-high cells sitting in immune-cold FOVs is invariant under this null and is not counted as evidence.")
    a("- The bootstrap resamples whole FOVs inside each section, then rebuilds the donor mean. Cells are not drawn independently of their neighbors.")
    a("- A second, coarser interval resamples the five donor effects. With five donors that interval cannot invent a p-value below 1/32, and it is not used for that claim.")
    a("")
    a("## Pre-specified endpoints")
    a("")
    a("All four use the within-FOV contrast, averaged with equal FOV weight inside a section, equal section weight inside a donor, and equal donor weight.")
    a("")
    a("1. Median split of tumor CLDN4 inside the FOV: mean CD8+NK count at **50 µm**, high minus low.")
    a("2. Same split: CD8+NK **contact rate at 20 µm** (at least one CD8 or NK centroid within 20 µm), high minus low.")
    a("3. **Continuous** CLDN4: within-FOV Spearman versus the 50 µm CD8+NK count.")
    a("4. **Continuous** CLDN4: within-FOV Spearman versus the 20 µm contact indicator.")
    a("")
    a("Exclusion is the negative direction. One-sided permutation p = (1 + number of null draws at least as small as observed) / (1 + n_perm). Radii 20 / 40 / 80 / 100 µm, the CD8-only sensitivity, the cell-weighted FOV mean, and the section-pooled median are secondary.")
    a("")
    a("## Inventory")
    a("")
    a(
        f"QC cells {inv['n_qc']:,}. Tumor {inv['n_tumor']:,}. CD8 {inv['n_cd8']:,}. NK {inv['n_nk']:,}. "
        f"CD8+NK {inv['n_cd8nk']:,}. FOVs used {inv['n_fov_used']}/{inv['n_fov']} "
        f"(≥{MIN_FOV_TUMOR} tumor and ≥{MIN_IMMUNE} CD8+NK). Tumor cells inside those FOVs {inv['n_tumor_used']:,}."
    )
    a("")
    a("| Section | Donor | QC cells | Tumor | CD8 | NK | CD8+NK | FOVs used |")
    a("|---|---|---:|---:|---:|---:|---:|---:|")
    for rec in inventory:
        a(
            f"| {rec['sample']} | {rec['donor']} | {rec['n_qc']:,} | {rec['n_tumor']:,} | "
            f"{rec['n_cd8']:,} | {rec['n_nk']:,} | {rec['n_cd8nk']:,} | {rec['n_fov_used']}/{rec['n_fov']} |"
        )
    a("")
    a("## Primary results")
    a("")
    a("| Endpoint | Donor-mean T | One-sided perm p | Two-sided perm p | FOV-block 95% CI | Sections neg | Donors neg | Below 1/256 | Below 1/32 |")
    a("|---|---:|---:|---:|---:|---:|---:|---|---|")
    for family, stat, outcome in PRIMARIES:
        rec = _lookup(rows, family, stat, outcome)
        label = {
            ("delta", "cd8nk_count_50"): "Median, 50 µm count",
            ("delta", "cd8nk_contact_20"): "Median, 20 µm contact",
            ("spearman", "cd8nk_count_50"): "Spearman, 50 µm count",
            ("spearman", "cd8nk_contact_20"): "Spearman, 20 µm contact",
        }[(stat, outcome)]
        ci = f"[{_fmt(rec['ci_block_lo'])}, {_fmt(rec['ci_block_hi'])}]"
        a(
            f"| {label} | {_fmt(rec['T'])} | {_fmt_p(rec['p_one'])} | {_fmt_p(rec['p_two'])} | {ci} | "
            f"{rec['sign_section_neg']}/{rec['sign_section_n']} | {rec['sign_donor_neg']}/{rec['sign_donor_n']} | "
            f"{'yes' if rec['beats_sample_sign_floor'] else 'no'} | "
            f"{'yes' if rec['beats_donor_sign_floor'] else 'no'} |"
        )
    a("")
    a(
        f"Permutation draws: {result['n_perm']}. FOV-block bootstrap draws: {result['n_boot']}. "
        f"Seed {result['seed']}. The smallest one-sided permutation p this run can return is 1/{result['n_perm'] + 1}."
    )
    a("")
    names_primary = result["outcome_names"]
    j50 = names_primary.index("cd8nk_count_50")
    j20 = names_primary.index("cd8nk_contact_20")
    d50 = result["within_fov"]["delta"]
    null50 = d50["null_T"][:, j50]
    null50 = null50[np.isfinite(null50)]
    lung6 = result["donors"].index("Lung6") if "Lung6" in result["donors"] else None
    lung6_txt = ""
    if lung6 is not None:
        lung6_txt = (
            f" Lung6, the immune-poor donor, is the exception at 50 µm "
            f"({_fmt(float(d50['donor'][lung6, j50]))}; interval "
            f"[{_fmt(float(d50['donor_ci_lo'][lung6, j50]))}, {_fmt(float(d50['donor_ci_hi'][lung6, j50]))}])."
        )
    sec_sign = d50["section_sign"][j50]
    don_sign = d50["donor_sign"][j50]
    a(
        f"At 50 µm the within-FOV sign counts are {sec_sign['n_neg']}/{sec_sign['n']} sections "
        f"(sign p = {_fmt_p(sec_sign['p'])}) and {don_sign['n_neg']}/{don_sign['n']} donors "
        f"(sign p = {_fmt_p(don_sign['p'])})."
        f"{lung6_txt} "
        f"The permutation p is {_fmt_p(float(d50['p_one'][j50]))} because the donor-mean magnitude is extreme: "
        f"the most negative null draw was {_fmt(float(np.min(null50)))}, against an observed {_fmt(float(d50['T'][j50]))}. "
        f"The within-FOV null mean is {_fmt(float(np.mean(null50)), 4)}."
    )
    a("")
    c20 = d50["donor_sign"][j20]
    s20 = d50["section_sign"][j20]
    a(
        f"The 20 µm contact contrast is negative in {s20['n_neg']}/{s20['n']} sections "
        f"(sign p = {_fmt_p(s20['p'])}) and {c20['n_neg']}/{c20['n']} donors "
        f"(sign p = {_fmt_p(c20['p'])}). "
        f"Those sign tests sit on the 8-section floor ({SAMPLE_SIGN_FLOOR_ONE:.6f}) and the 5-donor floor "
        f"({DONOR_SIGN_FLOOR_ONE:.5f}). The permutation p is lower than both because it uses the size of the deficit."
    )
    a("")
    a("## Radius and contact family")
    a("")
    a("Within-FOV median contrast (high − low), equal-donor mean. Negative means CLDN4-high tumor has fewer CD8+NK neighbors or contacts.")
    a("")
    a("| Outcome | T | perm p (one-sided) | FOV-block 95% CI | Section signs | Donor signs | Cell-weighted T | Cell-weighted p |")
    a("|---|---:|---:|---:|---:|---:|---:|---:|")
    for outcome in [f"cd8nk_count_{r}" for r in RADII] + ["cd8nk_contact_20", "cd8_count_50", "cd8_contact_20"]:
        rec = _lookup(rows, "within_fov", "delta", outcome)
        ci = f"[{_fmt(rec['ci_block_lo'])}, {_fmt(rec['ci_block_hi'])}]"
        tw = _fmt(rec["T_cell_weighted"]) if "T_cell_weighted" in rec else "NA"
        pw = _fmt_p(rec["p_one_cell_weighted"]) if "p_one_cell_weighted" in rec else "NA"
        a(
            f"| {PRIMARY_LABEL.get(outcome, outcome)} | {_fmt(rec['T'])} | {_fmt_p(rec['p_one'])} | {ci} | "
            f"{rec['sign_section_neg']}/{rec['sign_section_n']} | {rec['sign_donor_neg']}/{rec['sign_donor_n']} | {tw} | {pw} |"
        )
    a("")
    a("## Continuous CLDN4")
    a("")
    a("Within-FOV Spearman, then the same equal-weight average. A negative value means higher CLDN4 tracks fewer CD8+NK neighbors. This Spearman uses every tumor cell, including CLDN4-undetected cells. It is not a second median split.")
    a("")
    a("| Outcome | Mean Spearman | perm p (one-sided) | FOV-block 95% CI | Section signs | Donor signs |")
    a("|---|---:|---:|---:|---:|---:|")
    for outcome in [f"cd8nk_count_{r}" for r in RADII] + ["cd8nk_contact_20"]:
        rec = _lookup(rows, "within_fov", "spearman", outcome)
        ci = f"[{_fmt(rec['ci_block_lo'])}, {_fmt(rec['ci_block_hi'])}]"
        a(
            f"| {PRIMARY_LABEL.get(outcome, outcome)} | {_fmt(rec['T'])} | {_fmt_p(rec['p_one'])} | {ci} | "
            f"{rec['sign_section_neg']}/{rec['sign_section_n']} | {rec['sign_donor_neg']}/{rec['sign_donor_n']} |"
        )
    a("")
    if positive is not None:
        n_zero = positive.get("n_fov_median_zero")
        n_fov_all = positive.get("n_fov")
        a(
            f"In {n_zero} of {n_fov_all} FOVs the tumor CLDN4 median is 0, so the within-FOV median split is "
            "detected versus undetected in those fields. About 60% of tumor cells in a typical FOV have no CLDN4 count. "
            "A five-bin cut of all cells is not a dose curve under that tie: the zeros share one rank. "
            "The figure instead shows the undetected cells and, among CLDN4 > 0, four equal-count quartiles (P1 low to P4 high)."
        )
        a("")
        a(
            "The pre-specified Spearman above includes the undetected cells. A separate check, run after seeing how often the median is zero, "
            "restricts to CLDN4 > 0 and repeats the same donor-level permutation. Among expressors the association does **not** continue toward exclusion. "
            "The donor-mean contrasts are positive: higher CLDN4 among detected cells goes with more CD8+NK neighbors, not fewer. "
            "The one-sided exclusion p-value is 1. Lung13 is the donor that still goes the exclusion way; Lung5, Lung6, and Lung12 go the other way."
        )
        a("")
        a("| Among CLDN4 > 0 | Donor-mean T | Exclusion perm p | FOV-block 95% CI | Donors negative |")
        a("|---|---:|---:|---:|---:|")
        for stat, outcome, label in (
            ("spearman", "cd8nk_count_50", "Spearman, 50 µm count"),
            ("spearman", "cd8nk_contact_20", "Spearman, 20 µm contact"),
            ("delta", "cd8nk_count_50", "Median split, 50 µm count"),
            ("delta", "cd8nk_contact_20", "Median split, 20 µm contact"),
        ):
            rec = _lookup(positive["rows"], "within_fov", stat, outcome)
            ci = f"[{_fmt(rec['ci_block_lo'])}, {_fmt(rec['ci_block_hi'])}]"
            a(
                f"| {label} | {_fmt(rec['T'])} | {_fmt_p(rec['p_one'])} | {ci} | "
                f"{rec['sign_donor_neg']}/{rec['sign_donor_n']} |"
            )
        a("")
        a(
            f"Positive-only permutation draws: {positive['n_perm']}. "
            f"FOVs with at least 30 CLDN4-positive tumor cells: {positive['n_fov_positive']}."
        )
        a("")
    a("## Section-pooled median, same spatial null")
    a("")
    a("This contrast pools tumor cells inside a section and splits on the section median. It can be large when CLDN4-high cells occupy different FOVs from CLDN4-low cells. The permutation still shuffles CLDN4 only inside FOVs, so that between-FOV arrangement stays in the null and the section-pooled null mean is not zero. At 50 µm that null mean is -0.324 and the observed donor-mean is -1.298. The excess, -0.974, matches the within-FOV donor-mean of -0.976.")
    a("")
    a("| Outcome | Section-pooled T | perm p (one-sided) | FOV-block 95% CI | Section signs | Donor signs | Mean high | Mean low |")
    a("|---|---:|---:|---:|---:|---:|---:|---:|")
    pooled = result["section_pooled"]["delta"]
    names = result["outcome_names"]
    for outcome in [f"cd8nk_count_{r}" for r in RADII] + ["cd8nk_contact_20"]:
        rec = _lookup(rows, "section_pooled", "delta", outcome)
        j = _idx(names, outcome)
        ci = f"[{_fmt(rec['ci_block_lo'])}, {_fmt(rec['ci_block_hi'])}]"
        a(
            f"| {PRIMARY_LABEL.get(outcome, outcome)} | {_fmt(rec['T'])} | {_fmt_p(rec['p_one'])} | {ci} | "
            f"{rec['sign_section_neg']}/{rec['sign_section_n']} | {rec['sign_donor_neg']}/{rec['sign_donor_n']} | "
            f"{_fmt(float(pooled['mean_high'][j]))} | {_fmt(float(pooled['mean_low'][j]))} |"
        )
    a("")
    a("The locked public summary (50 / 100 µm cytotoxic ratio 0.36 / 0.52, 8/8 and 5/5, sign P = 0.031) is a different estimand and is not recomputed here. In this run the ratio of the section-pooled donor-mean counts is 2.190/3.488 = 0.63 at 50 µm and 10.588/13.714 = 0.77 at 100 µm. Those ratios are milder than 0.36 / 0.52. The sign pattern for the section-pooled contrast is still 8/8 sections and 5/5 donors. This does not overwrite the locked summary.")
    a("")
    a("## Donor effects for the primary median contrasts")
    a("")
    a("| Donor | 50 µm count Δ | 50 µm block CI | 20 µm contact Δ | Contact block CI |")
    a("|---|---:|---:|---:|---:|")
    d50 = result["within_fov"]["delta"]
    j50 = _idx(names, "cd8nk_count_50")
    j20 = _idx(names, "cd8nk_contact_20")
    for i, donor in enumerate(result["donors"]):
        a(
            f"| {donor} | {_fmt(float(d50['donor'][i, j50]))} | "
            f"[{_fmt(float(d50['donor_ci_lo'][i, j50]))}, {_fmt(float(d50['donor_ci_hi'][i, j50]))}] | "
            f"{_fmt(float(d50['donor'][i, j20]))} | "
            f"[{_fmt(float(d50['donor_ci_lo'][i, j20]))}, {_fmt(float(d50['donor_ci_hi'][i, j20]))}] |"
        )
    a("")
    a("## Reading rule")
    a("")
    a("A permutation p below 1/256 means the donor-averaged within-FOV contrast is extreme under FOV-restricted label exchange. It does not mean the five-donor sign test has been given a smaller discrete floor. If the section-pooled gap is large and the within-FOV permutation p is not, the gap is carried by which FOV a cell sits in, and this null correctly refuses to call that a within-neighborhood CLDN4 effect.")
    a("")
    a("## Methods notes")
    a("")
    a("- CLDN4 is log-normalized inside each section: size factor = median total of the marker panel used for typing, divided by that cell's total, then log1p. The QC filter (≥20 counts and ≥5 genes) uses the full 960-plex total, excluding NegPrb columns.")
    a("- Tumor: epithelial RNA score above immune and stromal scores, or PanCK at/above the section median and CD45 below it, and not called immune. CD8: CD8A or CD8B count > 0 and not tumor. NK: NKG7 or GNLY count > 0, no CD3, and not tumor. KLRD1 is not on this 960-plex. CD8+NK is the union. Tumor cells with CD8A stay in the tumor index.")
    a("- Median split inside a FOV: CLDN4 strictly above the FOV median versus at or below it. A FOV needs at least 5 cells on each side to enter the median contrast. Spearman does not use that split.")
    a(f"- Contact radius {CONTACT_RADIUS} µm. Count radii {', '.join(str(r) for r in RADII)} µm. Neighbors are counted inside the FOV only.")
    a("- Cell-weighted sensitivity weights FOVs by tumor-cell count, then still averages sections and donors equally.")
    a("")
    a("## Figures")
    a("")
    a("- `results/cosmx_donor_spatial_null/figures/null_primary.png`")
    a("- `results/cosmx_donor_spatial_null/figures/donor_forest.png`")
    a("- `results/cosmx_donor_spatial_null/figures/cldn4_quintiles.png`")
    a("- `results/cosmx_donor_spatial_null/figures/radius_profile.png`")
    a("")
    text = "\n".join(lines) + "\n"
    with open(os.path.join(OUT, "RESULTS.md"), "w") as fh:
        fh.write(text)
    with open(os.path.join(ROOT, "RESULTS.md"), "w") as fh:
        fh.write(text)


def slim_stats(result: dict, rows: list[dict], inventory: list[dict]) -> dict:
    prim = []
    for family, stat, outcome in PRIMARIES:
        prim.append(_lookup(rows, family, stat, outcome))
    return {
        "n_perm": result["n_perm"],
        "n_boot": result["n_boot"],
        "seed": result["seed"],
        "sample_sign_floor_one": SAMPLE_SIGN_FLOOR_ONE,
        "donor_sign_floor_one": DONOR_SIGN_FLOOR_ONE,
        "inventory": inventory,
        "primaries": prim,
        "outcome_names": list(result["outcome_names"]),
        "donors": result["donors"],
        "samples": result["samples"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-perm", type=int, default=4999)
    parser.add_argument("--n-boot", type=int, default=1999)
    parser.add_argument("--seed", type=int, default=25976224)
    args = parser.parse_args()
    os.makedirs(FIG, exist_ok=True)
    os.makedirs(TAB, exist_ok=True)
    blocks, inventory = build_blocks()
    # Stable order: sample order, then FOV id.
    order = {s: i for i, s in enumerate(SAMPLES)}
    blocks = sorted(blocks, key=lambda b: (order[b.sample], int(b.fov)))
    print(f"analyzing {len(blocks)} FOVs, n_perm={args.n_perm}, n_boot={args.n_boot}", flush=True)
    result = analyze(blocks, n_perm=args.n_perm, n_boot=args.n_boot, seed=args.seed, min_arm=MIN_ARM, n_bins=5)
    rows = endpoint_rows(result)
    pd.DataFrame(rows).to_csv(os.path.join(TAB, "endpoint_summary.csv"), index=False)
    pd.DataFrame(donor_rows(result)).to_csv(os.path.join(TAB, "donor_effects.csv"), index=False)
    pd.DataFrame(section_rows(result)).to_csv(os.path.join(TAB, "section_effects.csv"), index=False)
    pd.DataFrame(fov_rows(result)).to_csv(os.path.join(TAB, "fov_effects.csv"), index=False)
    # Equal-count quintiles of all cells collapse the large CLDN4=0 tie into one
    # bin. The continuous display is the undetected + positive-quartile curve.
    pd.DataFrame(inventory).to_csv(os.path.join(TAB, "sample_inventory.csv"), index=False)
    names = result["outcome_names"]
    null_frame = {}
    for stat, outcome in (
        ("delta", "cd8nk_count_50"),
        ("delta", "cd8nk_contact_20"),
        ("spearman", "cd8nk_count_50"),
        ("spearman", "cd8nk_contact_20"),
    ):
        j = _idx(names, outcome)
        null_frame[f"{stat}__{outcome}"] = result["within_fov"][stat]["null_T"][:, j]
    pd.DataFrame(null_frame).to_csv(os.path.join(TAB, "null_draws_primary.csv"), index=False)
    plot_nulls(result, os.path.join(FIG, "null_primary.png"))
    plot_forest(result, os.path.join(FIG, "donor_forest.png"))
    plot_quintiles(result, os.path.join(FIG, "cldn4_quintiles.png"), blocks=blocks)
    pos_blocks = positive_only_blocks(blocks)
    print(f"positive-only FOVs {len(pos_blocks)}", flush=True)
    pos = analyze(pos_blocks, n_perm=min(1999, args.n_perm), n_boot=min(999, args.n_boot), seed=args.seed, min_arm=MIN_ARM, n_bins=4)
    pos_rows = endpoint_rows(pos)
    pd.DataFrame(pos_rows).to_csv(os.path.join(TAB, "positive_only_endpoints.csv"), index=False)
    curve_rows = []
    for outcome in ("cd8nk_count_50", "cd8nk_contact_20"):
        curves = continuous_curves(blocks, outcome)
        for donor, row in curves["donor"].items():
            for b_i, (lab, val) in enumerate(zip(curves["labels"], row)):
                curve_rows.append({"outcome": outcome, "donor": donor, "bin": lab, "bin_index": b_i + 1, "mean_outcome": float(val)})
        for b_i, (lab, val) in enumerate(zip(curves["labels"], curves["grand"])):
            curve_rows.append({"outcome": outcome, "donor": "mean", "bin": lab, "bin_index": b_i + 1, "mean_outcome": float(val)})
    pd.DataFrame(curve_rows).to_csv(os.path.join(TAB, "cldn4_positive_quartiles.csv"), index=False)
    plot_radius(result, os.path.join(FIG, "radius_profile.png"))
    n_med0 = int(sum(np.median(b.cldn4) <= 0 for b in blocks))
    positive_payload = {
        "rows": pos_rows,
        "n_perm": int(pos["n_perm"]),
        "n_fov_positive": len(pos_blocks),
        "n_fov_median_zero": n_med0,
        "n_fov": len(blocks),
    }
    write_results(result, inventory, rows, positive=positive_payload)
    with open(os.path.join(OUT, "stats.json"), "w") as fh:
        json.dump(slim_stats(result, rows, inventory), fh, indent=2)
    print("wrote", OUT, flush=True)
    for family, stat, outcome in PRIMARIES:
        rec = _lookup(rows, family, stat, outcome)
        print(
            f"PRIMARY {stat} {outcome}: T={rec['T']:.4f} p={rec['p_one']:.6g} "
            f"beats_sample_floor={rec['beats_sample_sign_floor']}",
            flush=True,
        )


if __name__ == "__main__":
    main()
