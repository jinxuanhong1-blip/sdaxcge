#!/usr/bin/env python3
"""GSE221322 GeoMx protein DSP: barrier (EpCAM/PanCk) vs immune proteins.

CLDN4 and TROP2/TACSTD2 are not on the 68-plex nCounter protein panel.
EpCAM (TACSTD1) and PanCk are the closest deposited barrier/epithelial proteins.

Public inputs only (GEO GSE221322 QC + norm CSVs + series matrix).
"""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

CONTROLS = ["Ms IgG1", "Rb IgG", "Ms IgG2a"]
HOUSEKEEPERS = ["Histone H3", "GAPDH", "S6"]
BARRIER = ["EpCAM", "PanCk"]
PRIMARY_IMMUNE = ["CD8", "CD3", "PD-1"]
NEIGHBORHOOD = [
    "CD8",
    "CD3",
    "PD-1",
    "CD45",
    "CD4",
    "GZMB",
    "CD20",
    "FOXP3",
    "CD68",
    "HLA-DR",
    "PD-L1",
    "CD56",
]
ABSENT_TARGETS = ["CLDN4", "TROP2", "TACSTD2", "Claudin-4", "Trop-2"]

PROTEIN_CLASS = {
    "Ms IgG1": "isotype_control",
    "Rb IgG": "isotype_control",
    "Ms IgG2a": "isotype_control",
    "Histone H3": "housekeeper",
    "GAPDH": "housekeeper",
    "S6": "housekeeper",
    "EpCAM": "barrier_epithelial",
    "PanCk": "barrier_epithelial",
    "CD8": "immune_T",
    "CD3": "immune_T",
    "PD-1": "immune_checkpoint",
    "CD45": "immune_leukocyte",
    "CD4": "immune_T",
    "GZMB": "immune_cytotoxic",
    "GZMA": "immune_cytotoxic",
    "CD20": "immune_B",
    "FOXP3": "immune_Treg",
    "CD68": "immune_myeloid",
    "CD163": "immune_myeloid",
    "CD14": "immune_myeloid",
    "CD11c": "immune_myeloid",
    "CD66b": "immune_myeloid",
    "HLA-DR": "immune_antigen_presentation",
    "PD-L1": "immune_checkpoint",
    "PD-L2": "immune_checkpoint",
    "CTLA4": "immune_checkpoint",
    "LAG3": "immune_checkpoint",
    "Tim-3": "immune_checkpoint",
    "VISTA": "immune_checkpoint",
    "IDO1": "immune_checkpoint",
    "CD56": "immune_NK",
    "CD45RO": "immune_memory",
}


def bh_fdr(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    adj = np.empty(n, dtype=float)
    prev = 1.0
    for rank_from_end, idx in enumerate(order[::-1]):
        rank = n - rank_from_end
        val = min(prev, p[idx] * n / rank)
        adj[idx] = val
        prev = val
    return adj.tolist()


def spearman_ci(rho: float, n: int, alpha: float = 0.05) -> tuple[float, float]:
    if n <= 3 or not np.isfinite(rho) or abs(rho) >= 1:
        return (float("nan"), float("nan"))
    z = np.arctanh(rho)
    se = 1.0 / math.sqrt(n - 3)
    zcrit = stats.norm.ppf(1 - alpha / 2)
    return float(np.tanh(z - zcrit * se)), float(np.tanh(z + zcrit * se))


def spearman_row(x: pd.Series, y: pd.Series, stratum: str, a: str, b: str) -> dict:
    m = x.notna() & y.notna()
    n = int(m.sum())
    if n < 5:
        rho = p = lo = hi = float("nan")
    else:
        rho, p = stats.spearmanr(x[m], y[m])
        lo, hi = spearman_ci(float(rho), n)
    return {
        "stratum": stratum,
        "x": a,
        "y": b,
        "n": n,
        "rho": float(rho) if np.isfinite(rho) else np.nan,
        "p": float(p) if np.isfinite(p) else np.nan,
        "rho_ci95_lo": lo,
        "rho_ci95_hi": hi,
    }


def partial_spearman(x, y, z) -> tuple[float, float]:
    rx, ry, rz = stats.rankdata(x), stats.rankdata(y), stats.rankdata(z)

    def resid(a, b):
        sl, inter = stats.linregress(b, a)[:2]
        return a - (sl * b + inter)

    r, p = stats.pearsonr(resid(rx, rz), resid(ry, rz))
    return float(r), float(p)


def parse_series_matrix(path: Path) -> pd.DataFrame:
    fields: dict[str, list[list[str]]] = defaultdict(list)
    titles = None
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith("!Sample_"):
            continue
        parts = [p.strip('"') for p in line.split("\t")]
        key = parts[0].replace("!Sample_", "")
        if key == "title":
            titles = parts[1:]
        fields[key].append(parts[1:])
    if titles is None:
        raise SystemExit(f"no !Sample_title in {path}")
    smeta = pd.DataFrame({"title": titles})
    char_names = [
        "roi",
        "segment",
        "area",
        "aoinucleicount",
        "patient_id",
        "response",
        "followup",
        "status",
    ]
    for name, row in zip(char_names, fields["characteristics_ch1"]):
        smeta[name] = [x.split(":", 1)[1].strip() if ":" in x else x for x in row]
    smeta = smeta.set_index("title")
    smeta["area"] = pd.to_numeric(smeta["area"], errors="coerce")
    smeta["aoinucleicount"] = pd.to_numeric(smeta["aoinucleicount"], errors="coerce")
    smeta["followup"] = pd.to_numeric(smeta["followup"], errors="coerce")
    smeta["core"] = [str(x).split("|")[1].strip() for x in smeta.index]
    return smeta


def load_inputs(data_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    qc_path = data_dir / "GSE221322_4301_protein_QC.csv"
    norm_path = data_dir / "GSE221322_4301_protein_norm.csv"
    mtx_path = data_dir / "GSE221322_series_matrix.txt"
    if not qc_path.exists() or not norm_path.exists() or not mtx_path.exists():
        raise SystemExit(
            f"Missing GEO files in {data_dir}. Run scripts/gse221322_dsp/00_fetch.sh"
        )
    qc = pd.read_csv(qc_path, index_col=0)
    norm = pd.read_csv(norm_path, index_col=0)
    if list(qc.index) != list(norm.index) or list(qc.columns) != list(norm.columns):
        raise SystemExit("QC and norm matrices are not aligned")
    smeta = parse_series_matrix(mtx_path)
    if set(smeta.index) != set(norm.index):
        raise SystemExit("series matrix titles do not match protein matrix rows")
    smeta = smeta.loc[norm.index]
    return qc, norm, smeta


def panel_table(proteins: list[str]) -> pd.DataFrame:
    rows = []
    lower = {p.lower(): p for p in proteins}
    for name in ABSENT_TARGETS:
        rows.append(
            {
                "protein": name,
                "on_deposited_panel": name.lower() in lower,
                "class": "requested_absent",
                "role": "CLDN4/TROP2 not on GPL29263 68-plex used here",
            }
        )
    for p in proteins:
        cls = PROTEIN_CLASS.get(p, "other")
        if p in BARRIER:
            role = "closest barrier / epithelial proxy (CLDN4/TROP2 absent)"
        elif p in PRIMARY_IMMUNE:
            role = "prespecified T / PD-1"
        elif p in NEIGHBORHOOD:
            role = "immune neighborhood"
        elif cls == "isotype_control":
            role = "isotype control"
        elif cls == "housekeeper":
            role = "housekeeper / loading"
        else:
            role = cls
        rows.append(
            {
                "protein": p,
                "on_deposited_panel": True,
                "class": cls,
                "role": role,
            }
        )
    return pd.DataFrame(rows)


def qc_signal_table(qc: pd.DataFrame, smeta: pd.DataFrame) -> pd.DataFrame:
    igg = qc[CONTROLS].mean(axis=1)
    rows = []
    for p in BARRIER + NEIGHBORHOOD + CONTROLS + HOUSEKEEPERS:
        if p not in qc.columns:
            continue
        ratio = qc[p] / igg.replace(0, np.nan)
        rows.append(
            {
                "protein": p,
                "n_aoi": int(qc[p].notna().sum()),
                "qc_median": float(qc[p].median()),
                "qc_min": float(qc[p].min()),
                "qc_max": float(qc[p].max()),
                "median_sn_vs_mean_igg": float(ratio.median()),
                "frac_aoi_sn_gt_2": float((ratio > 2).mean()),
                "frac_aoi_sn_gt_3": float((ratio > 3).mean()),
            }
        )
    return pd.DataFrame(rows)


def correlate_block(df: pd.DataFrame, xs: list[str], ys: list[str], stratum: str) -> pd.DataFrame:
    rows = []
    for a in xs:
        for b in ys:
            if a == b:
                continue
            rows.append(spearman_row(df[a], df[b], stratum, a, b))
    out = pd.DataFrame(rows)
    if not out.empty:
        # BH within each barrier protein (12 neighborhood tests), matching methods.md
        q = np.full(len(out), np.nan)
        for x in out["x"].unique():
            idx = out.index[out["x"] == x]
            q[idx] = bh_fdr(out.loc[idx, "p"].fillna(1).tolist())
        out["q_bh_within_barrier"] = q
    return out


def mannwhitney_highlow(df: pd.DataFrame, split_col: str, ys: list[str], stratum: str) -> pd.DataFrame:
    med = float(df[split_col].median())
    hi = df[df[split_col] >= med]
    lo = df[df[split_col] < med]
    rows = []
    for y in ys:
        u, p = stats.mannwhitneyu(hi[y], lo[y], alternative="two-sided")
        rows.append(
            {
                "stratum": stratum,
                "split": split_col,
                "split_median": med,
                "y": y,
                "n_high": int(len(hi)),
                "n_low": int(len(lo)),
                "high_median": float(hi[y].median()),
                "low_median": float(lo[y].median()),
                "delta_high_minus_low": float(hi[y].median() - lo[y].median()),
                "U": float(u),
                "p": float(p),
            }
        )
    out = pd.DataFrame(rows)
    out["q_bh"] = bh_fdr(out["p"].tolist())
    return out


def style_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=8)
    ax.yaxis.label.set_size(9)
    ax.xaxis.label.set_size(9)


def savefig(fig, out_dir: Path, stem: str):
    fig.tight_layout()
    fig.savefig(out_dir / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(out_dir / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def fig_forest(corr: pd.DataFrame, title: str, out_dir: Path, stem: str):
    sub = corr.copy()
    sub["label"] = sub["x"] + " vs " + sub["y"]
    sub = sub.sort_values(["x", "rho"], ascending=[True, True])
    fig, ax = plt.subplots(figsize=(7.2, 6.4))
    y = np.arange(len(sub))
    colors = ["#1f4e79" if x == "EpCAM" else "#8a4b08" for x in sub["x"]]
    ax.axvline(0, color="#888888", lw=0.8)
    ax.hlines(y, sub["rho_ci95_lo"], sub["rho_ci95_hi"], color="#b0b0b0", lw=1.2)
    ax.scatter(sub["rho"], y, c=colors, s=28, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(sub["label"], fontsize=8)
    ax.set_xlabel("Spearman ρ (Fisher-z 95% CI)")
    ax.set_title(title, fontsize=10)
    for i, (_, r) in enumerate(sub.iterrows()):
        star = " *" if r["p"] < 0.05 else ""
        ax.text(
            0.98,
            i,
            f"ρ={r['rho']:.2f} p={r['p']:.3g} n={int(r['n'])}{star}",
            va="center",
            ha="right",
            fontsize=6.5,
            transform=ax.get_yaxis_transform(),
            color="#333333",
        )
    style_axes(ax)
    ax.set_xlim(-0.85, 0.85)
    savefig(fig, out_dir, stem)


def fig_scatter(tum: pd.DataFrame, out_dir: Path):
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.6))
    pairs = [("CD45", "#1b4f72"), ("CD8", "#196f3d")]
    for ax, (y, color) in zip(axes, pairs):
        ax.scatter(tum["EpCAM"], tum[y], s=22, c=color, alpha=0.85, edgecolors="none")
        rho, p = stats.spearmanr(tum["EpCAM"], tum[y])
        ax.set_xlabel("Tumour AOI EpCAM (norm)")
        ax.set_ylabel(f"Tumour AOI {y} (norm)")
        ax.set_title(f"ρ = {rho:.2f}, p = {p:.3g}, n = {len(tum)}", fontsize=9)
        style_axes(ax)
    fig.suptitle("GSE221322 tumour segments — EpCAM vs leukocyte / CD8 protein", fontsize=10)
    savefig(fig, out_dir, "scatter_tumour_epcam_cd45_cd8")


def fig_boxes(tum: pd.DataFrame, out_dir: Path):
    med = float(tum["EpCAM"].median())
    hi = tum[tum["EpCAM"] >= med]
    lo = tum[tum["EpCAM"] < med]
    fig, axes = plt.subplots(1, 3, figsize=(8.6, 3.5))
    for ax, y in zip(axes, ["CD45", "CD8", "CD3"]):
        data = [lo[y].values, hi[y].values]
        bp = ax.boxplot(data, tick_labels=["EpCAM-low", "EpCAM-high"], widths=0.55, patch_artist=True)
        for patch, c in zip(bp["boxes"], ["#aed6f1", "#1f4e79"]):
            patch.set_facecolor(c)
            patch.set_alpha(0.85)
        for medn in bp["medians"]:
            medn.set_color("black")
        u, p = stats.mannwhitneyu(hi[y], lo[y], alternative="two-sided")
        ax.set_ylabel(f"{y} (norm)")
        ax.set_title(f"MW p = {p:.3g}\nn_lo={len(lo)} n_hi={len(hi)}", fontsize=8)
        style_axes(ax)
    fig.suptitle("Tumour AOIs split at median EpCAM", fontsize=10)
    savefig(fig, out_dir, "box_tumour_epcam_highlow")


def fig_heatmap(tum_corr: pd.DataFrame, st_corr: pd.DataFrame, out_dir: Path):
    def pivot(df, xs):
        sub = df[df["x"].isin(xs)]
        mat = sub.pivot(index="y", columns="x", values="rho")
        return mat.reindex(NEIGHBORHOOD)

    m1 = pivot(tum_corr, BARRIER)
    m2 = pivot(st_corr, BARRIER)
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 5.6))
    vmin, vmax = -0.55, 0.55
    for ax, mat, title in zip(
        axes,
        [m1, m2],
        ["Tumour AOI Spearman ρ", "Stroma AOI Spearman ρ"],
    ):
        im = ax.imshow(mat.values, cmap="RdBu_r", vmin=vmin, vmax=vmax, aspect="auto")
        ax.set_xticks(range(mat.shape[1]))
        ax.set_xticklabels(mat.columns, fontsize=8)
        ax.set_yticks(range(mat.shape[0]))
        ax.set_yticklabels(mat.index, fontsize=8)
        ax.set_title(title, fontsize=9)
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                val = mat.values[i, j]
                if np.isfinite(val):
                    ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7)
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.7, label="Spearman ρ")
    fig.suptitle("Barrier vs immune neighborhood proteins (author-normalized)", fontsize=10)
    fig.savefig(out_dir / "heatmap_rho_compartments.png", dpi=300, bbox_inches="tight")
    fig.savefig(out_dir / "heatmap_rho_compartments.pdf", bbox_inches="tight")
    plt.close(fig)


def fig_cross(cross: pd.DataFrame, out_dir: Path):
    sub = cross[cross["x_tumour"] == "EpCAM"].copy()
    sub = sub.set_index("y_stroma").loc[NEIGHBORHOOD].reset_index()
    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    y = np.arange(len(sub))
    ax.axvline(0, color="#888888", lw=0.8)
    ax.hlines(y, sub["rho_ci95_lo"], sub["rho_ci95_hi"], color="#b0b0b0", lw=1.2)
    ax.scatter(sub["rho"], y, c="#6c3483", s=28, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(sub["y_stroma"], fontsize=8)
    ax.set_xlabel("Spearman ρ (tumour EpCAM vs same-core stroma protein)")
    ax.set_title(f"Cross-compartment, n = {int(sub['n'].iloc[0])} paired cores", fontsize=10)
    style_axes(ax)
    ax.set_xlim(-0.7, 0.7)
    savefig(fig, out_dir, "forest_cross_tumour_epcam_vs_stroma")


def write_stats(path: Path, d: dict):
    lines = []
    for k, v in d.items():
        if isinstance(v, float):
            lines.append(f"{k}\t{v:.6g}")
        else:
            lines.append(f"{k}\t{v}")
    path.write_text("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", type=Path, default=Path("data/gse221322"))
    ap.add_argument("--out-dir", type=Path, default=Path("results/gse221322_dsp"))
    args = ap.parse_args()
    data_dir = args.data_dir
    out = args.out_dir
    fig_dir = out / "figures"
    tab_dir = out / "tables"
    fig_dir.mkdir(parents=True, exist_ok=True)
    tab_dir.mkdir(parents=True, exist_ok=True)

    qc, norm, smeta = load_inputs(data_dir)
    proteins = list(norm.columns)
    smeta = smeta.copy()
    smeta["fail_nuclei"] = smeta["aoinucleicount"] < 20

    panel = panel_table(proteins)
    panel.to_csv(tab_dir / "panel_inventory.csv", index=False)

    sn = qc_signal_table(qc, smeta)
    sn.to_csv(tab_dir / "qc_signal_vs_igg.csv", index=False)

    aoi = smeta.join(norm)
    aoi["igg_mean"] = qc.loc[aoi.index, CONTROLS].mean(axis=1)
    aoi_keep = aoi[~aoi["fail_nuclei"]].copy()
    tum = aoi_keep[aoi_keep["segment"] == "Tumour"].copy()
    st = aoi_keep[aoi_keep["segment"] == "Stroma"].copy()

    aoi_keep[["core", "segment", "patient_id", "response", "status", "followup", "area", "aoinucleicount"] + BARRIER + NEIGHBORHOOD].to_csv(
        tab_dir / "aoi_compact.csv"
    )

    tum_corr = correlate_block(tum, BARRIER, NEIGHBORHOOD, "tumour_aoi")
    st_corr = correlate_block(st, BARRIER, NEIGHBORHOOD, "stroma_aoi")

    def patient_mean(df):
        cols = BARRIER + NEIGHBORHOOD
        g = df.groupby("patient_id")
        outp = g[cols].mean()
        outp["response"] = g["response"].first()
        outp["n_cores"] = g["core"].nunique()
        return outp

    tp = patient_mean(tum)
    sp = patient_mean(st)
    tp_corr = correlate_block(tp, BARRIER, NEIGHBORHOOD, "tumour_patient")
    sp_corr = correlate_block(sp, BARRIER, NEIGHBORHOOD, "stroma_patient")

    paired_cores = sorted(set(tum["core"]) & set(st["core"]))
    t_c = tum.set_index("core")
    s_c = st.set_index("core")
    cross_rows = []
    for a in BARRIER:
        for b in NEIGHBORHOOD:
            row = spearman_row(
                t_c.loc[paired_cores, a],
                s_c.loc[paired_cores, b],
                "cross_core_tumourX_stromaY",
                a,
                b,
            )
            row["x_tumour"] = a
            row["y_stroma"] = b
            cross_rows.append(row)
    cross = pd.DataFrame(cross_rows)
    cross["q_bh"] = bh_fdr(cross["p"].fillna(1).tolist())

    common_pts = sorted(set(tp.index) & set(sp.index))
    cross_p_rows = []
    for a in BARRIER:
        for b in NEIGHBORHOOD:
            row = spearman_row(
                tp.loc[common_pts, a],
                sp.loc[common_pts, b],
                "cross_patient_tumourX_stromaY",
                a,
                b,
            )
            row["x_tumour"] = a
            row["y_stroma"] = b
            cross_p_rows.append(row)
    cross_p = pd.DataFrame(cross_p_rows)

    all_corr = pd.concat([tum_corr, st_corr, tp_corr, sp_corr], ignore_index=True)
    all_corr.to_csv(tab_dir / "spearman_barrier_vs_immune.csv", index=False)
    cross.to_csv(tab_dir / "spearman_cross_core.csv", index=False)
    cross_p.to_csv(tab_dir / "spearman_cross_patient.csv", index=False)

    mw_epcam = mannwhitney_highlow(tum, "EpCAM", NEIGHBORHOOD, "tumour_aoi")
    mw_panck = mannwhitney_highlow(tum, "PanCk", PRIMARY_IMMUNE + ["CD45"], "tumour_aoi")
    pd.concat([mw_epcam, mw_panck], ignore_index=True).to_csv(
        tab_dir / "mannwhitney_highlow.csv", index=False
    )

    # Partial Spearman (tumour EpCAM vs immune | PanCk or nuclei)
    partial_rows = []
    for y in PRIMARY_IMMUNE + ["CD45", "CD4"]:
        for zname, z in [("PanCk", tum["PanCk"]), ("nuclei", tum["aoinucleicount"]), ("Histone H3", tum["Histone H3"])]:
            r, p = partial_spearman(tum["EpCAM"], tum[y], z)
            partial_rows.append(
                {"x": "EpCAM", "y": y, "covariate": zname, "n": len(tum), "partial_rho": r, "p": p}
            )
    pd.DataFrame(partial_rows).to_csv(tab_dir / "partial_spearman_tumour.csv", index=False)

    # Area-normalized QC sensitivity
    qc_tum = smeta.join(qc)
    qc_tum = qc_tum[(qc_tum["segment"] == "Tumour") & (~qc_tum["fail_nuclei"])].copy()
    sens_rows = []
    for method, denom in [("area_1e5", qc_tum["area"] / 1e5), ("per_nucleus", qc_tum["aoinucleicount"])]:
        tmp = {}
        for p in BARRIER + NEIGHBORHOOD:
            tmp[p] = np.log2(qc_tum[p] / denom)
        tmp = pd.DataFrame(tmp, index=qc_tum.index)
        block = correlate_block(tmp, BARRIER, PRIMARY_IMMUNE + ["CD45", "CD4"], f"tumour_qc_{method}")
        block["normalization"] = method
        sens_rows.append(block)
    pd.concat(sens_rows, ignore_index=True).to_csv(tab_dir / "spearman_qc_area_nuclei_sensitivity.csv", index=False)

    # ICI response (sanity check vs published stromal EpCAM; not the B6 test)
    def response_mw(df, stratum):
        sub = df[df["response"].isin(["Responder", "Non-responder"])]
        rows = []
        for y in BARRIER + PRIMARY_IMMUNE + ["CD45"]:
            a = sub.loc[sub["response"] == "Responder", y]
            b = sub.loc[sub["response"] == "Non-responder", y]
            u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            rows.append(
                {
                    "stratum": stratum,
                    "protein": y,
                    "n_responder": int(len(a)),
                    "n_nonresponder": int(len(b)),
                    "responder_median": float(a.median()),
                    "nonresponder_median": float(b.median()),
                    "U": float(u),
                    "p": float(p),
                }
            )
        return pd.DataFrame(rows)

    resp = pd.concat(
        [response_mw(tp, "tumour_patient"), response_mw(sp, "stroma_patient")],
        ignore_index=True,
    )
    resp.to_csv(tab_dir / "response_mannwhitney.csv", index=False)

    # Key numbers
    def grab(df, stratum, x, y):
        hit = df[(df["stratum"] == stratum) & (df["x"] == x) & (df["y"] == y)]
        return hit.iloc[0]

    e_cd45 = grab(tum_corr, "tumour_aoi", "EpCAM", "CD45")
    e_cd8 = grab(tum_corr, "tumour_aoi", "EpCAM", "CD8")
    e_cd3 = grab(tum_corr, "tumour_aoi", "EpCAM", "CD3")
    e_pd1 = grab(tum_corr, "tumour_aoi", "EpCAM", "PD-1")
    p_cd45 = grab(tp_corr, "tumour_patient", "EpCAM", "CD45")
    s_cd45 = grab(st_corr, "stroma_aoi", "EpCAM", "CD45")
    x_pd1 = cross[(cross["x_tumour"] == "EpCAM") & (cross["y_stroma"] == "PD-1")].iloc[0]
    epcam_panck_t = stats.spearmanr(tum["EpCAM"], tum["PanCk"])
    cd8_cd3_t = stats.spearmanr(tum["CD8"], tum["CD3"])
    pd1_igg = stats.spearmanr(tum["PD-1"], tum["Ms IgG1"])
    part_cd45_panck = partial_spearman(tum["EpCAM"], tum["CD45"], tum["PanCk"])
    mw_cd45 = mw_epcam[mw_epcam["y"] == "CD45"].iloc[0]
    st_epcam_resp = resp[(resp["stratum"] == "stroma_patient") & (resp["protein"] == "EpCAM")].iloc[0]
    pd1_sn = sn[sn["protein"] == "PD-1"].iloc[0]

    # BH among the 12 tumour EpCAM vs neighborhood tests
    epcam_tum = tum_corr[tum_corr["x"] == "EpCAM"].copy()
    epcam_tum = epcam_tum[epcam_tum["y"].isin(NEIGHBORHOOD)]
    q_map = dict(zip(epcam_tum["y"], bh_fdr(epcam_tum["p"].tolist())))

    key = {
        "n_proteins_deposited": len(proteins),
        "cldn4_on_panel": False,
        "trop2_on_panel": False,
        "n_aoi_deposited": int(len(aoi)),
        "n_cores_deposited": int(aoi["core"].nunique()),
        "n_patients": int(aoi["patient_id"].nunique()),
        "n_tumour_aoi_qc": int(len(tum)),
        "n_stroma_aoi_qc": int(len(st)),
        "n_tumour_patients": int(len(tp)),
        "n_stroma_patients": int(len(sp)),
        "n_paired_cores": int(len(paired_cores)),
        "dropped_tumour_aoi": "4301 Protein | 026 | Tumour (0 nuclei, area 54.1)",
        "tumour_EpCAM_CD45_rho": float(e_cd45["rho"]),
        "tumour_EpCAM_CD45_p": float(e_cd45["p"]),
        "tumour_EpCAM_CD45_q": float(q_map["CD45"]),
        "tumour_EpCAM_CD8_rho": float(e_cd8["rho"]),
        "tumour_EpCAM_CD8_p": float(e_cd8["p"]),
        "tumour_EpCAM_CD3_rho": float(e_cd3["rho"]),
        "tumour_EpCAM_CD3_p": float(e_cd3["p"]),
        "tumour_EpCAM_PD1_rho": float(e_pd1["rho"]),
        "tumour_EpCAM_PD1_p": float(e_pd1["p"]),
        "patient_EpCAM_CD45_rho": float(p_cd45["rho"]),
        "patient_EpCAM_CD45_p": float(p_cd45["p"]),
        "stroma_EpCAM_CD45_rho": float(s_cd45["rho"]),
        "stroma_EpCAM_CD45_p": float(s_cd45["p"]),
        "cross_EpCAM_stroma_PD1_rho": float(x_pd1["rho"]),
        "cross_EpCAM_stroma_PD1_p": float(x_pd1["p"]),
        "tumour_EpCAM_PanCk_rho": float(epcam_panck_t.statistic),
        "tumour_EpCAM_PanCk_p": float(epcam_panck_t.pvalue),
        "tumour_CD8_CD3_rho": float(cd8_cd3_t.statistic),
        "tumour_PD1_IgG1_rho": float(pd1_igg.statistic),
        "tumour_PD1_IgG1_p": float(pd1_igg.pvalue),
        "pd1_median_sn_vs_igg": float(pd1_sn["median_sn_vs_mean_igg"]),
        "pd1_frac_sn_gt_2": float(pd1_sn["frac_aoi_sn_gt_2"]),
        "partial_EpCAM_CD45_PanCk_rho": part_cd45_panck[0],
        "partial_EpCAM_CD45_PanCk_p": part_cd45_panck[1],
        "mw_EpCAM_CD45_p": float(mw_cd45["p"]),
        "mw_EpCAM_CD45_delta": float(mw_cd45["delta_high_minus_low"]),
        "stroma_EpCAM_response_p": float(st_epcam_resp["p"]),
        "stroma_EpCAM_R_median": float(st_epcam_resp["responder_median"]),
        "stroma_EpCAM_NR_median": float(st_epcam_resp["nonresponder_median"]),
    }
    write_stats(out / "stats.txt", key)
    tp.to_csv(tab_dir / "patient_tumour_means.csv")
    sp.to_csv(tab_dir / "patient_stroma_means.csv")

    fig_forest(
        tum_corr,
        f"Tumour AOIs (n = {len(tum)}) — barrier vs immune protein",
        fig_dir,
        "forest_tumour_barrier_vs_immune",
    )
    fig_forest(
        st_corr,
        f"Stroma AOIs (n = {len(st)}) — barrier vs immune protein",
        fig_dir,
        "forest_stroma_barrier_vs_immune",
    )
    fig_scatter(tum, fig_dir)
    fig_boxes(tum, fig_dir)
    fig_heatmap(tum_corr, st_corr, fig_dir)
    fig_cross(cross, fig_dir)

    provenance = {
        "dataset": "GSE221322",
        "platform": "GPL29263 NanoString GeoMx Human Protein for nCounter 2020",
        "paper": "Monkman et al., Immunology 2023, PMID 37022147, DOI 10.1111/imm.13646",
        "files": [
            "GSE221322_4301_protein_QC.csv",
            "GSE221322_4301_protein_norm.csv",
            "GSE221322_series_matrix.txt",
        ],
        "urls": {
            "geo": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221322",
            "qc": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE221322&format=file&file=GSE221322_4301_protein_QC.csv.gz",
            "norm": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE221322&format=file&file=GSE221322_4301_protein_norm.csv.gz",
            "matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221322/matrix/GSE221322_series_matrix.txt.gz",
        },
        "primary_matrix": "author-deposited protein_norm.csv",
        "absent_on_panel": ["CLDN4", "TROP2", "TACSTD2"],
        "closest_barrier": ["EpCAM", "PanCk"],
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "key": key,
    }
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps(key, indent=2))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
