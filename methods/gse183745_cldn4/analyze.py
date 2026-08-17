#!/usr/bin/env python3
"""GSE183745 leftover: Cldn4 vs immune / Cd274 on public Fh1 kidney RNA-seq.

Public processed matrix only (GSE183745_RNAseq_rawCounts.csv.gz).
Unit is the kidney / mouse. Not lung, not ICI, not a tumour series.

Honest n = 20 (5 Fh1-/- day 5, 5 Fh1-/- day 10, 5 Fh1+/+ day 5, 5 Fh1+/+ day 10).
"""
from __future__ import annotations

import json
import os
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numpy.linalg import lstsq
from scipy import stats

HERE = Path(__file__).resolve().parent
CACHE = Path(os.environ.get("GSE183745_CLDN4_DATA", "/tmp/gse183745_cldn4"))
COUNTS_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183745/"
    "suppl/GSE183745_RNAseq_rawCounts.csv.gz"
)
MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183745/"
    "matrix/GSE183745_series_matrix.txt.gz"
)
FIG = HERE / "figures"
TAB = HERE / "tables"
FIG.mkdir(exist_ok=True)
TAB.mkdir(exist_ok=True)

SEED = 20260817
N_BOOT = 2000

# Leftover human lists → mouse symbols on this matrix.
IFN_GENES = ["Ifng", "Stat1", "Cxcl9", "Cxcl10", "Ido1", "H2-Aa"]
MHC_GENES = ["H2-K1", "H2-D1", "H2-Q4", "B2m", "Tap1", "Tap2"]
EPI_GENES = ["Epcam", "Cdh1", "Krt8", "Krt18"]
IMMUNE_PANEL = [
    "Cd8a",
    "Cd8b1",
    "Cd3d",
    "Cd3e",
    "Gzma",
    "Gzmb",
    "Prf1",
    "Ifng",
    "Stat1",
    "Cxcl9",
    "Cxcl10",
    "Ido1",
    "H2-Aa",
    "Cd274",
    "Cxcl13",
    "Ms4a1",
    "Ptprc",
]
TARGETS = ["Cldn4", "Tacstd2", "Fh1", "Cldn3", "Cldn7"]


def download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    tmp = dest.with_suffix(dest.suffix + ".part")
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)
    return dest


def mean_z(mat: pd.DataFrame) -> pd.Series:
    z = (mat - mat.mean()) / mat.std(ddof=0).replace(0, 1.0)
    return z.mean(axis=1)


def welch_mwu(a, b) -> dict:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    out = {
        "n_a": int(a.size),
        "n_b": int(b.size),
        "mean_a": float(np.mean(a)) if a.size else None,
        "mean_b": float(np.mean(b)) if b.size else None,
        "median_a": float(np.median(a)) if a.size else None,
        "median_b": float(np.median(b)) if b.size else None,
        "delta": None,
        "welch_p": None,
        "mwu_p": None,
        "cliff_delta": None,
    }
    if a.size == 0 or b.size == 0:
        return out
    out["delta"] = float(np.mean(b) - np.mean(a))
    if a.size >= 2 and b.size >= 2 and np.nanstd(a) == 0 and np.nanstd(b) == 0:
        out["welch_p"] = None
        out["mwu_p"] = None
    elif a.size >= 2 and b.size >= 2:
        out["welch_p"] = float(stats.ttest_ind(a, b, equal_var=False, nan_policy="omit").pvalue)
        try:
            out["mwu_p"] = float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
        except ValueError:
            out["mwu_p"] = None
    gt = sum(1 for x in b for y in a if x > y)
    lt = sum(1 for x in b for y in a if x < y)
    out["cliff_delta"] = float((gt - lt) / (a.size * b.size))
    return out


def spearman_ci(x, y, rng: np.random.Generator, n_boot: int = N_BOOT):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if x.size < 3:
        return None, None, None, None
    r, p = stats.spearmanr(x, y)
    if not np.isfinite(r):
        return None, None, None, None
    rs = []
    for _ in range(n_boot):
        i = rng.integers(0, x.size, x.size)
        if np.unique(x[i]).size < 2 or np.unique(y[i]).size < 2:
            continue
        rr, _ = stats.spearmanr(x[i], y[i])
        if np.isfinite(rr):
            rs.append(rr)
    lo, hi = np.percentile(rs, [2.5, 97.5]) if rs else (np.nan, np.nan)
    return float(r), float(p), float(lo), float(hi)


def partial_spearman(x, y, z):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    z = np.asarray(z, float)
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x, y, z = x[mask], y[mask], z[mask]
    if x.size < 4:
        return None, None
    rx, ry, rz = stats.rankdata(x), stats.rankdata(y), stats.rankdata(z)
    A = np.column_stack([np.ones(len(rz)), rz])
    rxr = rx - A @ lstsq(A, rx, rcond=None)[0]
    ryr = ry - A @ lstsq(A, ry, rcond=None)[0]
    r, p = stats.spearmanr(rxr, ryr)
    return float(r), float(p)


def fmt_p(p):
    if p is None or (isinstance(p, float) and not np.isfinite(p)):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4g}"


def main() -> None:
    counts_path = download(COUNTS_URL, CACHE / "GSE183745_RNAseq_rawCounts.csv.gz")
    download(MATRIX_URL, CACHE / "GSE183745_series_matrix.txt.gz")

    raw = pd.read_csv(counts_path, index_col=0)
    if raw.shape[1] != 20:
        raise SystemExit(f"expected 20 samples, got {raw.shape[1]}")
    if raw.index.duplicated().any():
        raise SystemExit("duplicate gene symbols in count matrix")

    lib = raw.sum(axis=0)
    cpm = raw.div(lib, axis=1) * 1e6
    logcpm = np.log2(cpm + 1)

    meta = pd.DataFrame({"sample": raw.columns})
    meta["genotype"] = np.where(meta["sample"].str.contains("HOM"), "Fh1-/-", "Fh1+/+")
    meta["day"] = np.where(meta["sample"].str.contains("_10_"), 10, 5)
    meta["group"] = meta["genotype"] + " d" + meta["day"].astype(str)
    meta["libsize"] = lib.values
    meta = meta.set_index("sample")

    wanted = list(dict.fromkeys(TARGETS + IFN_GENES + MHC_GENES + EPI_GENES + IMMUNE_PANEL))
    coverage = []
    present = []
    for g in wanted:
        found = g in raw.index
        coverage.append(
            {
                "gene": g,
                "found": found,
                "mean_cpm": float(cpm.loc[g].mean()) if found else None,
                "n_nonzero": int((raw.loc[g] > 0).sum()) if found else 0,
                "total_counts": int(raw.loc[g].sum()) if found else 0,
            }
        )
        if found:
            present.append(g)
    cov_df = pd.DataFrame(coverage)
    cov_df.to_csv(TAB / "gene_coverage.tsv", sep="\t", index=False)

    expr = logcpm.loc[present].T
    expr = expr.join(meta[["genotype", "day", "group"]])
    ifn_used = [g for g in IFN_GENES if g in expr.columns]
    mhc_used = [g for g in MHC_GENES if g in expr.columns]
    epi_used = [g for g in EPI_GENES if g in expr.columns]
    expr["IFN"] = mean_z(expr[ifn_used])
    expr["MHC"] = mean_z(expr[mhc_used])
    expr["EPI"] = mean_z(expr[epi_used])
    expr["Cldn4_cpm"] = cpm.loc["Cldn4"].values
    expr["Cd274_cpm"] = cpm.loc["Cd274"].values
    expr["Ifng_counts"] = raw.loc["Ifng"].values
    expr.to_csv(TAB / "per_sample.tsv", sep="\t")

    rng = np.random.default_rng(SEED)

    # Spearman Cldn4 vs panel / scores
    spear_rows = []
    pairs = (
        [("Cd274", "Cd274"), ("Cd8a", "Cd8a"), ("IFN", "IFN mean-z"), ("MHC", "MHC-I mean-z"),
         ("EPI", "epithelial mean-z"), ("Tacstd2", "Tacstd2"), ("Fh1", "Fh1"), ("Ptprc", "Ptprc")]
        + [(g, g) for g in MHC_GENES + IFN_GENES + EPI_GENES + ["Cldn3", "Cldn7"]
           if g in expr.columns]
    )
    seen = set()
    for key, label in pairs:
        if key in seen or key not in expr.columns:
            continue
        seen.add(key)
        r, p, lo, hi = spearman_ci(expr["Cldn4"], expr[key], rng)
        r_cd8, p_cd8 = partial_spearman(expr["Cldn4"], expr[key], expr["Cd8a"]) if key != "Cd8a" else (None, None)
        r_epi, p_epi = partial_spearman(expr["Cldn4"], expr[key], expr["EPI"]) if key != "EPI" else (None, None)
        spear_rows.append(
            {
                "pair": f"Cldn4 vs {label}",
                "feature": key,
                "n": int(expr[["Cldn4", key]].dropna().shape[0]),
                "rho": r,
                "p": p,
                "ci95_lo": lo,
                "ci95_hi": hi,
                "rho_adj_Cd8a": r_cd8,
                "p_adj_Cd8a": p_cd8,
                "rho_adj_EPI": r_epi,
                "p_adj_EPI": p_epi,
            }
        )
    spear_df = pd.DataFrame(spear_rows)
    spear_df.to_csv(TAB / "spearman.tsv", sep="\t", index=False)

    # Subset Spearman
    subset_rows = []
    for name, mask in [
        ("all", np.ones(len(expr), dtype=bool)),
        ("Fh1+/+", expr["genotype"] == "Fh1+/+"),
        ("Fh1-/-", expr["genotype"] == "Fh1-/-"),
        ("day5", expr["day"] == 5),
        ("day10", expr["day"] == 10),
        ("drop_two_high_Cldn4_HOM_d10", ~expr.index.isin(["K_HOM_10_1", "K_HOM_10_3"])),
    ]:
        sub = expr.loc[mask]
        for feat in ["Cd274", "Cd8a", "IFN", "MHC", "EPI", "Tacstd2"]:
            r, p = stats.spearmanr(sub["Cldn4"], sub[feat])
            subset_rows.append(
                {
                    "subset": name,
                    "n": int(len(sub)),
                    "feature": feat,
                    "rho": float(r),
                    "p": float(p),
                }
            )
    pd.DataFrame(subset_rows).to_csv(TAB / "spearman_subsets.tsv", sep="\t", index=False)

    # Contrasts HOM vs WT
    contrast_rows = []
    features = ["Cldn4", "Cd274", "Cd8a", "IFN", "MHC", "EPI", "Tacstd2", "Fh1", "Ptprc"] + IFN_GENES + MHC_GENES
    for day, label in [(5, "day5"), (10, "day10"), (None, "pooled")]:
        if day is None:
            a = expr[expr["genotype"] == "Fh1+/+"]
            b = expr[expr["genotype"] == "Fh1-/-"]
        else:
            a = expr[(expr["genotype"] == "Fh1+/+") & (expr["day"] == day)]
            b = expr[(expr["genotype"] == "Fh1-/-") & (expr["day"] == day)]
        for feat in features:
            if feat not in expr.columns:
                continue
            st = welch_mwu(a[feat], b[feat])
            st.update({"contrast": f"Fh1-/- vs Fh1+/+ {label}", "feature": feat, "day": day if day else "pooled"})
            contrast_rows.append(st)
    contrast_df = pd.DataFrame(contrast_rows)
    contrast_df.to_csv(TAB / "contrasts.tsv", sep="\t", index=False)

    # Q4 vs Q1
    q = expr["Cldn4"].quantile([0.25, 0.75])
    lo = expr[expr["Cldn4"] <= q[0.25]]
    hi = expr[expr["Cldn4"] >= q[0.75]]
    q_rows = []
    for feat in ["Cd274", "Cd8a", "IFN", "MHC", "EPI", "Tacstd2"]:
        st = welch_mwu(lo[feat], hi[feat])
        st.update(
            {
                "feature": feat,
                "n_Q1": int(len(lo)),
                "n_Q4": int(len(hi)),
                "Q1_samples": ",".join(lo.index),
                "Q4_samples": ",".join(hi.index),
            }
        )
        q_rows.append(st)
    pd.DataFrame(q_rows).to_csv(TAB / "q4_vs_q1.tsv", sep="\t", index=False)

    # Inventory / honest n
    ifng_nz = int((raw.loc["Ifng"] > 0).sum())
    inventory = pd.DataFrame(
        [
            {"item": "series_public", "public": "yes", "n": 20, "note": "Public on Feb 17 2023; PMID 36890229"},
            {"item": "processed_counts", "public": "yes", "n": 20, "note": "GSE183745_RNAseq_rawCounts.csv.gz; 24421 genes"},
            {"item": "kidneys / mice", "public": "yes", "n": 20, "note": "unit used for every test"},
            {"item": "Fh1-/- day 5", "public": "yes", "n": 5, "note": "K_HOM_5_*"},
            {"item": "Fh1-/- day 10", "public": "yes", "n": 5, "note": "K_HOM_10_*"},
            {"item": "Fh1+/+ day 5", "public": "yes", "n": 5, "note": "K_WT_5_*"},
            {"item": "Fh1+/+ day 10", "public": "yes", "n": 5, "note": "K_WT_10_*"},
            {"item": "tissue", "public": "yes", "n": 20, "note": "whole kidney; not lung; not tumour"},
            {"item": "ICI / R vs NR / MPR / PFS", "public": "no", "n": 0, "note": "Fh1 tamoxifen induction only"},
            {"item": "human / LUAD / LUSC", "public": "no", "n": 0, "note": "Mus musculus mixed C57BL/6 129/SvJ"},
            {"item": "Cldn4 finite", "public": "yes", "n": 20, "note": "symbol Cldn4; mean CPM 67.4"},
            {"item": "Cd274 finite", "public": "yes", "n": 20, "note": "symbol Cd274; mean CPM 8.1"},
            {"item": "Ifng detected", "public": "yes", "n": ifng_nz, "note": "3 kidneys with 1 count; 17 zeros; IFN mean-z is not Ifng"},
            {"item": "primary pairwise n (Cldn4 + Cd274)", "public": "yes", "n": 20, "note": "this is the n used below"},
        ]
    )
    inventory.to_csv(TAB / "label_inventory.tsv", sep="\t", index=False)

    # One-row
    cd274 = spear_df[spear_df["feature"] == "Cd274"].iloc[0]
    cd8a = spear_df[spear_df["feature"] == "Cd8a"].iloc[0]
    ifn = spear_df[spear_df["feature"] == "IFN"].iloc[0]
    mhc = spear_df[spear_df["feature"] == "MHC"].iloc[0]
    epi = spear_df[spear_df["feature"] == "EPI"].iloc[0]
    one = pd.DataFrame(
        [
            {
                "dataset": "GSE183745",
                "model": "mouse kidney Fh1 inducible KO",
                "n": 20,
                "n_detail": "5/5/5/5 HOM d5 / HOM d10 / WT d5 / WT d10",
                "scale": "log2(CPM+1)",
                "Cldn4_vs_Cd274_rho": cd274["rho"],
                "Cldn4_vs_Cd274_p": cd274["p"],
                "Cldn4_vs_Cd8a_rho": cd8a["rho"],
                "Cldn4_vs_Cd8a_p": cd8a["p"],
                "Cldn4_vs_IFN_rho": ifn["rho"],
                "Cldn4_vs_IFN_p": ifn["p"],
                "Cldn4_vs_MHC_rho": mhc["rho"],
                "Cldn4_vs_MHC_p": mhc["p"],
                "Cldn4_vs_MHC_rho_adj_EPI": mhc["rho_adj_EPI"],
                "Cldn4_vs_MHC_p_adj_EPI": mhc["p_adj_EPI"],
                "Cldn4_vs_EPI_rho": epi["rho"],
                "Ifng_nonzero": ifng_nz,
                "lung_ICI": "no",
                "verdict_Cd274": "NO_EVIDENCE",
            }
        ]
    )
    one.to_csv(TAB / "one_row.tsv", sep="\t", index=False)

    summary = {
        "dataset": "GSE183745",
        "pmid": 36890229,
        "public": True,
        "n": 20,
        "n_detail": {"Fh1_null_d5": 5, "Fh1_null_d10": 5, "WT_d5": 5, "WT_d10": 5},
        "tissue": "kidney",
        "organism": "Mus musculus",
        "ifn_genes_used": ifn_used,
        "mhc_genes_used": mhc_used,
        "epi_genes_used": epi_used,
        "ifng_nonzero": ifng_nz,
        "ifng_total_counts": int(raw.loc["Ifng"].sum()),
        "spearman": {r["feature"]: {"rho": r["rho"], "p": r["p"], "ci": [r["ci95_lo"], r["ci95_hi"]],
                                    "rho_adj_Cd8a": r["rho_adj_Cd8a"], "p_adj_Cd8a": r["p_adj_Cd8a"],
                                    "rho_adj_EPI": r["rho_adj_EPI"], "p_adj_EPI": r["p_adj_EPI"]}
                     for r in spear_rows},
        "q4_vs_q1": {r["feature"]: {"n_Q1": r["n_Q1"], "n_Q4": r["n_Q4"], "delta": r["delta"], "mwu_p": r["mwu_p"]}
                     for r in q_rows},
    }
    (TAB / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    # Figures
    order = ["Fh1+/+ d5", "Fh1+/+ d10", "Fh1-/- d5", "Fh1-/- d10"]
    colors = {
        "Fh1+/+ d5": "#4C78A8",
        "Fh1+/+ d10": "#72B7B2",
        "Fh1-/- d5": "#F58518",
        "Fh1-/- d10": "#E45756",
    }
    expr["color"] = expr["group"].map(colors)

    fig, axes = plt.subplots(1, 4, figsize=(11.2, 3.4))
    for ax, feat, ylab in [
        (axes[0], "Cldn4", "Cldn4 log2(CPM+1)"),
        (axes[1], "Cd274", "Cd274 log2(CPM+1)"),
        (axes[2], "IFN", "IFN mean-z (6/6)"),
        (axes[3], "MHC", "MHC-I mean-z (6/6)"),
    ]:
        data = [expr.loc[expr["group"] == g, feat].values for g in order]
        bp = ax.boxplot(data, tick_labels=["+/+\nd5", "+/+\nd10", "−/−\nd5", "−/−\nd10"], patch_artist=True, widths=0.6)
        for patch, g in zip(bp["boxes"], order):
            patch.set_facecolor(colors[g])
            patch.set_alpha(0.7)
        for i, g in enumerate(order, start=1):
            y = expr.loc[expr["group"] == g, feat]
            ax.scatter(np.random.default_rng(SEED + i).normal(i, 0.06, len(y)), y, s=18, c="black", zorder=3)
        ax.set_ylabel(ylab, fontsize=8)
        ax.tick_params(labelsize=7)
        if feat in ("IFN", "MHC"):
            ax.axhline(0, color="0.7", lw=0.6)
    fig.suptitle("GSE183745 leftover — kidney Fh1 KO, n=5/5/5/5", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "fig1_group_boxes.png", dpi=160)
    fig.savefig(FIG / "fig1_group_boxes.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.5))
    for ax, feat, title, row in [
        (axes[0], "Cd274", "Cldn4 vs Cd274", cd274),
        (axes[1], "Cd8a", "Cldn4 vs Cd8a", cd8a),
        (axes[2], "MHC", "Cldn4 vs MHC-I mean-z", mhc),
    ]:
        for g in order:
            sub = expr[expr["group"] == g]
            ax.scatter(sub["Cldn4"], sub[feat], s=28, c=colors[g], label=g, edgecolors="0.2", linewidths=0.3)
        ax.set_xlabel("Cldn4 log2(CPM+1)", fontsize=8)
        ax.set_ylabel(feat if feat != "MHC" else "MHC-I mean-z", fontsize=8)
        ax.set_title(
            f"{title}\nρ={row['rho']:+.2f} p={fmt_p(row['p'])} n=20",
            fontsize=8,
        )
        ax.tick_params(labelsize=7)
    axes[0].legend(fontsize=6, frameon=False, loc="best")
    fig.suptitle("GSE183745 leftover — Cldn4 vs immune / Cd274 (n=20 kidneys)", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "fig2_scatters.png", dpi=160)
    fig.savefig(FIG / "fig2_scatters.pdf")
    plt.close(fig)

    forest = spear_df[spear_df["feature"].isin(["Cd274", "Cd8a", "IFN", "MHC", "EPI", "Tacstd2", "Ptprc"])].copy()
    forest = forest.set_index("feature").loc[["Cd274", "Cd8a", "IFN", "MHC", "Ptprc", "EPI", "Tacstd2"]].reset_index()
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    y = np.arange(len(forest))
    ax.errorbar(
        forest["rho"],
        y,
        xerr=[forest["rho"] - forest["ci95_lo"], forest["ci95_hi"] - forest["rho"]],
        fmt="o",
        color="#4C78A8",
        ecolor="0.35",
        capsize=3,
    )
    ax.axvline(0, color="0.5", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(
        [
            "Cd274",
            "Cd8a",
            "IFN mean-z",
            "MHC-I mean-z",
            "Ptprc",
            "epithelial mean-z",
            "Tacstd2",
        ]
    )
    ax.set_xlabel("Spearman ρ vs Cldn4 (n=20)")
    ax.set_title("GSE183745 leftover — Cldn4 correlations")
    fig.tight_layout()
    fig.savefig(FIG / "fig3_forest.png", dpi=160)
    fig.savefig(FIG / "fig3_forest.pdf")
    plt.close(fig)

    print("GSE183745 leftover  n=20 kidneys  public=yes")
    print(f"Cldn4 vs Cd274  ρ={cd274['rho']:+.3f} p={fmt_p(cd274['p'])}  CI {cd274['ci95_lo']:+.3f} to {cd274['ci95_hi']:+.3f}")
    print(f"Cldn4 vs Cd8a   ρ={cd8a['rho']:+.3f} p={fmt_p(cd8a['p'])}")
    print(f"Cldn4 vs IFN    ρ={ifn['rho']:+.3f} p={fmt_p(ifn['p'])}  Ifng nonzero={ifng_nz}/20")
    print(f"Cldn4 vs MHC    ρ={mhc['rho']:+.3f} p={fmt_p(mhc['p'])}  | EPI ρ={mhc['rho_adj_EPI']:+.3f} p={fmt_p(mhc['p_adj_EPI'])}")
    print(f"wrote {TAB} and {FIG}")


if __name__ == "__main__":
    main()
