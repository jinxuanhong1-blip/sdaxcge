#!/usr/bin/env python3
"""GSE207422 malignant CLDN4-high vs low: TJ/keratin and IFN in the same cells.

Patient is the inferential unit. Within each post-treatment patient, malignant
cells are split on CLDN4 log1p(CP10k) Q4 vs Q1. Module means are compared with
a paired Wilcoxon. Cell-level p-values are stored and labeled exploratory
(pseudoreplication).

CLDN4 is the splitter and is excluded from the TJ module.
TACSTD2 is a companion gene and is never a gate.
Author CopyKAT barcodes are not on GEO. A3-malignant-like = epithelial AND
zero UMI for SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3 (same rule as the dual-high
CLDN4 slice).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gene_sets import (
    A3_NORMAL_LUNG,
    BROAD_NORMAL_LUNG,
    LINEAGES,
    MODULES,
)

HERE = Path(__file__).resolve().parents[1]

PAPER_GROUP = {
    "BD_immune01": "TN",
    "BD_immune02": "NMPR",
    "BD_immune03": "MPR",
    "BD_immune04": "NMPR",
    "BD_immune05": "TN",
    "BD_immune06": "MPR",  # pCR
    "BD_immune07": "NMPR",
    "BD_immune08": "TN",
    "BD_immune09": "NMPR",
    "BD_immune10": "NMPR",
    "BD_immune11": "MPR",
    "BD_immune12": "NMPR",
    "BD_immune13": "NMPR",
    "BD_immune14": "MPR",
    "BD_immune15": "NMPR",
}

COLOR = {"MPR": "#d1495b", "NMPR": "#2c6eaf", "TN": "#6b6b6b"}

# Primary: ≥20 malignant cells and ≥8 cells in each CLDN4 tail.
MIN_MAL_PRIMARY = 20
MIN_TAIL = 8
MIN_MAL_SENS = 10
MIN_TAIL_SENS = 5


def score_log1p(expr: dict[str, np.ndarray], genes: list[str], n: int) -> np.ndarray:
    mats = [np.log1p(expr[g].astype(np.float32)) for g in genes if g in expr]
    if not mats:
        return np.zeros(n, dtype=np.float32)
    return np.mean(np.vstack(mats), axis=0)


def assign_lineage(expr: dict[str, np.ndarray], n: int) -> np.ndarray:
    names = list(LINEAGES)
    scores = np.vstack([score_log1p(expr, LINEAGES[k], n) for k in names])
    best = scores.argmax(axis=0)
    best_val = scores.max(axis=0)
    second = np.partition(scores, -2, axis=0)[-2]
    assigned = np.array(names, dtype=object)[best]
    assigned[(best_val < 0.15) | ((best_val - second) < 0.05)] = "other"
    return assigned


def log1p_cp10k(umi: np.ndarray, lib: np.ndarray) -> np.ndarray:
    umi = np.asarray(umi, dtype=np.float64)
    lib = np.asarray(lib, dtype=np.float64)
    out = np.full(umi.shape, np.nan, dtype=np.float64)
    ok = lib > 0
    out[ok] = np.log1p(umi[ok] / lib[ok] * 1e4)
    return out


def module_cp10k(expr: dict[str, np.ndarray], genes: list[str], lib: np.ndarray) -> tuple[np.ndarray, list[str]]:
    present = [g for g in genes if g in expr]
    if not present:
        return np.full(len(lib), np.nan, dtype=np.float64), present
    stacked = np.vstack([log1p_cp10k(expr[g], lib) for g in present])
    return stacked.mean(axis=0), present


def load_sample_meta(path: Path) -> pd.DataFrame:
    raw = pd.read_excel(path)
    raw = raw.dropna(subset=["Sample"]).copy()
    raw = raw[~raw["Sample"].astype(str).str.contains("RECIST|MPR:|NMPR:|pCR:", regex=True)]
    raw["Sample"] = raw["Sample"].astype(str)
    raw["paper_group"] = raw["Sample"].map(PAPER_GROUP)
    raw["path_response"] = raw["Pathologic Response"].replace({"pCR": "MPR"})
    raw["timing"] = np.where(
        raw["Resource"].astype(str).str.contains("Pre", case=False, na=False),
        "pre",
        "post",
    )
    return raw


def spearman_row(x, y, contrast: str) -> dict:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = int(len(x))
    if n < 4:
        return {
            "contrast": contrast,
            "n": n,
            "spearman_rho": np.nan,
            "spearman_p": np.nan,
            "note": "too_few_samples",
        }
    rho, p = stats.spearmanr(x, y)
    return {
        "contrast": contrast,
        "n": n,
        "spearman_rho": float(rho),
        "spearman_p": float(p),
        "note": "",
    }


def wilcoxon_paired(a, b, contrast: str) -> dict:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m], b[m]
    n = int(len(a))
    if n < 3:
        return {
            "contrast": contrast,
            "n": n,
            "mean_high": float(np.mean(a)) if n else np.nan,
            "mean_low": float(np.mean(b)) if n else np.nan,
            "median_delta": np.nan,
            "wilcoxon_w": np.nan,
            "p_value": np.nan,
            "n_high_gt_low": int(np.sum(a > b)) if n else 0,
            "note": "too_few_samples",
        }
    d = a - b
    if np.allclose(d, 0):
        w, p = 0.0, 1.0
    else:
        try:
            w, p = stats.wilcoxon(a, b, alternative="two-sided", zero_method="wilcox")
            w, p = float(w), float(p)
        except ValueError:
            w, p = np.nan, np.nan
    return {
        "contrast": contrast,
        "n": n,
        "mean_high": float(np.mean(a)),
        "mean_low": float(np.mean(b)),
        "median_delta": float(np.median(d)),
        "mean_delta": float(np.mean(d)),
        "wilcoxon_w": w,
        "p_value": p,
        "n_high_gt_low": int(np.sum(a > b)),
        "note": "",
    }


def mwu_row(a, b, contrast: str) -> dict:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    n_a, n_b = int(len(a)), int(len(b))
    if n_a < 2 or n_b < 2:
        return {
            "contrast": contrast,
            "n_a": n_a,
            "n_b": n_b,
            "mean_a": np.nan,
            "mean_b": np.nan,
            "mwu_u": np.nan,
            "p_value": np.nan,
            "note": "too_few_samples",
        }
    method = "exact" if (n_a + n_b) <= 20 else "asymptotic"
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided", method=method)
    return {
        "contrast": contrast,
        "n_a": n_a,
        "n_b": n_b,
        "mean_a": float(np.mean(a)),
        "mean_b": float(np.mean(b)),
        "delta_mean": float(np.mean(a) - np.mean(b)),
        "mwu_u": float(u),
        "p_value": float(p),
        "note": "",
    }


def fmt_p(p) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    return f"{p:.2e}" if p < 0.001 else f"{p:.4f}"


def fmt_r(r) -> str:
    if r is None or not np.isfinite(r):
        return "NA"
    return f"{r:.3f}"


def fmt_paired(row: dict) -> str:
    if row.get("note") == "too_few_samples":
        return f"n={row['n']} (too few)"
    return (
        f"high {row['mean_high']:.3f} vs low {row['mean_low']:.3f} "
        f"(Δmed={row['median_delta']:+.3f}); W={row['wilcoxon_w']:.1f}; "
        f"p={fmt_p(row['p_value'])}; n={row['n']} paired"
    )


def fmt_rho(row: dict) -> str:
    if row.get("note") == "too_few_samples":
        return f"n={row['n']} (too few)"
    return f"ρ={row['spearman_rho']:.3f}, p={fmt_p(row['spearman_p'])}, n={row['n']}"


def savefig(fig, path: Path) -> None:
    fig.savefig(path.with_suffix(".png"), dpi=160, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", type=Path, default=Path("data/GSE207422"))
    ap.add_argument("--outdir", type=Path, default=HERE)
    args = ap.parse_args()
    outdir = args.outdir
    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    extracted = np.load(args.workdir / "extracted_cldn4_tj_ifn.npz", allow_pickle=True)
    cells = extracted["cells"]
    genes = list(extracted["genes"])
    mat = extracted["mat"]
    total = extracted["total"]
    n_nonzero = extracted["n_nonzero_genes"]
    n_genes_in_matrix = int(extracted["n_genes_in_matrix"][0])
    missing = [str(x) for x in extracted["missing"]] if "missing" in extracted.files else []
    expr = {g: mat[i] for i, g in enumerate(genes)}
    n = len(cells)
    if "CLDN4" not in expr:
        raise SystemExit("CLDN4 missing from extracted matrix")

    sample_ids = np.array([str(c).rsplit("_", 1)[0] for c in cells])
    lineage = assign_lineage(expr, n)
    is_epi = lineage == "epithelial"

    a3_normal = np.zeros(n, dtype=np.int32)
    for g in A3_NORMAL_LUNG:
        if g in expr:
            a3_normal += expr[g]
    is_malig = is_epi & (a3_normal == 0)

    broad_score = score_log1p(expr, [g for g in BROAD_NORMAL_LUNG if g in expr], n)
    epi_cut = float(np.quantile(broad_score[is_epi], 0.60)) if is_epi.any() else np.inf
    is_malig_broad = is_epi & (broad_score < epi_cut)

    cldn4_cp = log1p_cp10k(expr["CLDN4"], total)
    tac_cp = log1p_cp10k(expr["TACSTD2"], total) if "TACSTD2" in expr else np.full(n, np.nan)

    module_scores: dict[str, np.ndarray] = {}
    module_present: dict[str, list[str]] = {}
    for name, glist in MODULES.items():
        sc, present = module_cp10k(expr, glist, total)
        module_scores[name] = sc
        module_present[name] = present

    meta = load_sample_meta(args.workdir / "GSE207422_NSCLC_scRNAseq_metadata.xlsx")
    meta_map = meta.set_index("Sample")

    # Per-cell table (malignant A3 only written compressed later if needed)
    per_cell = pd.DataFrame(
        {
            "cell": cells,
            "Sample": sample_ids,
            "lineage": lineage,
            "is_epi": is_epi,
            "is_malig_a3": is_malig,
            "is_malig_broad": is_malig_broad,
            "total_umi": total,
            "n_genes": n_nonzero,
            "CLDN4_log1p_cp10k": cldn4_cp,
            "TACSTD2_log1p_cp10k": tac_cp,
        }
    )
    for name, sc in module_scores.items():
        per_cell[f"mod_{name}"] = sc
    per_cell["paper_group"] = per_cell["Sample"].map(PAPER_GROUP)
    per_cell["timing"] = per_cell["Sample"].map(
        lambda s: meta_map.loc[s, "timing"] if s in meta_map.index else np.nan
    )

    # Per-patient occupancy and means
    patient_rows = []
    paired_rows = []
    for sid, g in per_cell.groupby("Sample", sort=True):
        group = PAPER_GROUP.get(sid, "NA")
        timing = str(g["timing"].iloc[0]) if len(g) else "NA"
        rec = {
            "Sample": sid,
            "paper_group": group,
            "timing": timing,
            "n_cells": int(len(g)),
            "n_epithelial": int(g["is_epi"].sum()),
            "n_malig_a3": int(g["is_malig_a3"].sum()),
            "n_malig_broad": int(g["is_malig_broad"].sum()),
        }
        for flag, tag in (("is_malig_a3", "a3"), ("is_epi", "epi"), ("is_malig_broad", "broad")):
            sub = g.loc[g[flag]]
            rec[f"mean_CLDN4_{tag}"] = float(sub["CLDN4_log1p_cp10k"].mean()) if len(sub) else np.nan
            rec[f"mean_TACSTD2_{tag}"] = float(sub["TACSTD2_log1p_cp10k"].mean()) if len(sub) else np.nan
            rec[f"pct_CLDN4pos_{tag}"] = float((sub["CLDN4_log1p_cp10k"] > 0).mean() * 100) if len(sub) else np.nan
            for name in MODULES:
                rec[f"mean_{name}_{tag}"] = float(sub[f"mod_{name}"].mean()) if len(sub) else np.nan
        patient_rows.append(rec)

        # Within-patient CLDN4 Q4 vs Q1 on A3-malignant (and epithelial sensitivity)
        for flag, tag, min_n, min_tail in (
            ("is_malig_a3", "a3", MIN_MAL_PRIMARY, MIN_TAIL),
            ("is_malig_a3", "a3_sens", MIN_MAL_SENS, MIN_TAIL_SENS),
            ("is_epi", "epi", MIN_MAL_PRIMARY, MIN_TAIL),
            ("is_malig_broad", "broad", MIN_MAL_PRIMARY, MIN_TAIL),
        ):
            sub = g.loc[g[flag]].copy()
            n_sub = int(len(sub))
            if n_sub < min_n:
                paired_rows.append(
                    {
                        "Sample": sid,
                        "paper_group": group,
                        "timing": timing,
                        "compartment": tag,
                        "n_cells": n_sub,
                        "n_high": 0,
                        "n_low": 0,
                        "eligible": False,
                        "reason": f"n_cells<{min_n}",
                    }
                )
                continue
            q = sub["CLDN4_log1p_cp10k"]
            q75 = float(q.quantile(0.75))
            q25 = float(q.quantile(0.25))
            hi = sub.loc[q >= q75]
            lo = sub.loc[q <= q25]
            if len(hi) < min_tail or len(lo) < min_tail:
                paired_rows.append(
                    {
                        "Sample": sid,
                        "paper_group": group,
                        "timing": timing,
                        "compartment": tag,
                        "n_cells": n_sub,
                        "n_high": int(len(hi)),
                        "n_low": int(len(lo)),
                        "eligible": False,
                        "reason": f"tail<{min_tail}",
                    }
                )
                continue
            prow = {
                "Sample": sid,
                "paper_group": group,
                "timing": timing,
                "compartment": tag,
                "n_cells": n_sub,
                "n_high": int(len(hi)),
                "n_low": int(len(lo)),
                "eligible": True,
                "reason": "",
                "cldn4_high": float(hi["CLDN4_log1p_cp10k"].mean()),
                "cldn4_low": float(lo["CLDN4_log1p_cp10k"].mean()),
                "tacstd2_high": float(hi["TACSTD2_log1p_cp10k"].mean()),
                "tacstd2_low": float(lo["TACSTD2_log1p_cp10k"].mean()),
                "umi_high": float(hi["total_umi"].mean()),
                "umi_low": float(lo["total_umi"].mean()),
                "ngenes_high": float(hi["n_genes"].mean()),
                "ngenes_low": float(lo["n_genes"].mean()),
            }
            for name in MODULES:
                prow[f"{name}_high"] = float(hi[f"mod_{name}"].mean())
                prow[f"{name}_low"] = float(lo[f"mod_{name}"].mean())
                prow[f"{name}_delta"] = prow[f"{name}_high"] - prow[f"{name}_low"]
            paired_rows.append(prow)

    patients = pd.DataFrame(patient_rows)
    paired = pd.DataFrame(paired_rows)

    # Tests: primary = post A3-malignant, eligible under MIN_MAL_PRIMARY / MIN_TAIL
    tests = []
    post = patients[patients["timing"] == "post"].copy()

    def add_test(kind, row):
        rec = {"kind": kind, **row}
        tests.append(rec)

    # Sample-level Spearman on post patients with ≥20 A3-malignant
    post_a3 = post[post["n_malig_a3"] >= MIN_MAL_PRIMARY].copy()
    for mod in MODULES:
        add_test(
            "sample_spearman",
            spearman_row(
                post_a3["mean_CLDN4_a3"],
                post_a3[f"mean_{mod}_a3"],
                f"post_a3_CLDN4_vs_{mod}_min{MIN_MAL_PRIMARY}",
            ),
        )
    add_test(
        "sample_spearman",
        spearman_row(post_a3["mean_CLDN4_a3"], post_a3["mean_TACSTD2_a3"], "post_a3_CLDN4_vs_TACSTD2"),
    )

    # Complete-case epithelial (all 12 post)
    post_epi = post[post["n_epithelial"] >= MIN_MAL_PRIMARY].copy()
    for mod in ("tj_keratin", "ifn_isg"):
        add_test(
            "sample_spearman",
            spearman_row(
                post_epi["mean_CLDN4_epi"],
                post_epi[f"mean_{mod}_epi"],
                f"post_epi_CLDN4_vs_{mod}",
            ),
        )

    # Paired within-patient high vs low
    for compartment in ("a3", "a3_sens", "epi", "broad"):
        sub = paired[(paired["compartment"] == compartment) & (paired["eligible"]) & (paired["timing"] == "post")]
        for mod in MODULES:
            add_test(
                "paired_high_vs_low",
                {
                    **wilcoxon_paired(
                        sub[f"{mod}_high"],
                        sub[f"{mod}_low"],
                        f"post_{compartment}_{mod}_CLDN4high_vs_low",
                    ),
                    "compartment": compartment,
                    "module": mod,
                    "n_patients_eligible": int(len(sub)),
                    "patients": ",".join(sub["Sample"].tolist()),
                },
            )
        add_test(
            "paired_high_vs_low",
            {
                **wilcoxon_paired(sub["tacstd2_high"], sub["tacstd2_low"], f"post_{compartment}_TACSTD2_CLDN4high_vs_low"),
                "compartment": compartment,
                "module": "TACSTD2",
                "n_patients_eligible": int(len(sub)),
            },
        )
        add_test(
            "paired_high_vs_low",
            {
                **wilcoxon_paired(sub["umi_high"], sub["umi_low"], f"post_{compartment}_UMI_CLDN4high_vs_low"),
                "compartment": compartment,
                "module": "total_UMI",
                "n_patients_eligible": int(len(sub)),
            },
        )
        # Do the same patients that gain TJ/keratin lose IFN?
        if len(sub) >= 4:
            add_test(
                "delta_spearman",
                {
                    **spearman_row(
                        sub["tj_keratin_delta"],
                        sub["ifn_isg_delta"],
                        f"post_{compartment}_delta_tj_keratin_vs_ifn_isg",
                    ),
                    "compartment": compartment,
                },
            )

    # NMPR vs MPR on sample-level module means (honest 8 vs 4 or fewer)
    for tag, n_col in (("a3", "n_malig_a3"), ("epi", "n_epithelial")):
        use = post[post[n_col] >= MIN_MAL_PRIMARY]
        nmpr = use[use["paper_group"] == "NMPR"]
        mpr = use[use["paper_group"] == "MPR"]
        for col in [f"mean_CLDN4_{tag}", f"mean_tj_keratin_{tag}", f"mean_ifn_isg_{tag}"]:
            add_test(
                "nmpr_vs_mpr",
                {
                    **mwu_row(nmpr[col], mpr[col], f"post_{col}_NMPR_vs_MPR"),
                    "n_nmpr": int(len(nmpr)),
                    "n_mpr": int(len(mpr)),
                    "patients_nmpr": ",".join(nmpr["Sample"].tolist()),
                    "patients_mpr": ",".join(mpr["Sample"].tolist()),
                },
            )

    # Cell-level exploratory (A3 malignant, post only)
    post_cells = per_cell[(per_cell["timing"] == "post") & (per_cell["is_malig_a3"])]
    for mod in ("tj_keratin", "ifn_isg", "tj_no_cldn4", "keratin", "ifn_core6", "ctrl_oxphos"):
        add_test(
            "cell_spearman_exploratory",
            {
                **spearman_row(
                    post_cells["CLDN4_log1p_cp10k"],
                    post_cells[f"mod_{mod}"],
                    f"post_a3_cell_CLDN4_vs_{mod}",
                ),
                "note2": "pseudoreplication; not a claim",
            },
        )

    tests_df = pd.DataFrame(tests)

    # Lineage counts
    lineage_counts = (
        per_cell.groupby(["Sample", "lineage"]).size().unstack(fill_value=0).reset_index()
    )

    # Figures
    prim = paired[(paired["compartment"] == "a3") & (paired["eligible"]) & (paired["timing"] == "post")].copy()
    prim = prim.sort_values("Sample")

    # Extra figure 1: paired high vs low for TJ/keratin and IFN (same patients, same cells)
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.2))
    for ax, mod, title in (
        (axes[0], "tj_keratin", "TJ/keratin module"),
        (axes[1], "ifn_isg", "IFN (ISG core) module"),
    ):
        xs = np.array([0, 1])
        for _, row in prim.iterrows():
            col = COLOR.get(row["paper_group"], "#444")
            ax.plot(
                xs,
                [row[f"{mod}_low"], row[f"{mod}_high"]],
                color=col,
                alpha=0.75,
                lw=1.4,
                marker="o",
                ms=5,
            )
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["CLDN4-low\n(Q1)", "CLDN4-high\n(Q4)"])
        ax.set_ylabel("mean log1p(CP10k)")
        ax.set_title(title)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    handles = [
        plt.Line2D([0], [0], color=COLOR["NMPR"], marker="o", label="NMPR"),
        plt.Line2D([0], [0], color=COLOR["MPR"], marker="o", label="MPR"),
    ]
    axes[1].legend(handles=handles, frameon=False, loc="best")
    fig.suptitle(
        f"GSE207422 malignant CLDN4-high vs low, same cells\n"
        f"paired patients n={len(prim)} (post; A3-malignant ≥{MIN_MAL_PRIMARY}, tail ≥{MIN_TAIL})",
        fontsize=11,
    )
    fig.tight_layout()
    savefig(fig, figdir / "fig_extra_paired_tj_ifn")

    # Extra figure 2: ΔTJ/keratin vs ΔIFN (same patients)
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    for grp, sub in prim.groupby("paper_group"):
        ax.scatter(
            sub["tj_keratin_delta"],
            sub["ifn_isg_delta"],
            c=COLOR[grp],
            s=48,
            label=grp,
            zorder=3,
        )
        for _, r in sub.iterrows():
            ax.annotate(r["Sample"].replace("BD_immune", "P"), (r["tj_keratin_delta"], r["ifn_isg_delta"]),
                        textcoords="offset points", xytext=(4, 3), fontsize=7, color=COLOR[grp])
    ax.axhline(0, color="#888", lw=0.8)
    ax.axvline(0, color="#888", lw=0.8)
    ax.set_xlabel("Δ TJ/keratin (CLDN4-high − low)")
    ax.set_ylabel("Δ IFN ISG (CLDN4-high − low)")
    ax.set_title(f"Same-cell module deltas, n={len(prim)} paired patients")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    savefig(fig, figdir / "fig_extra_delta_scatter")

    # Extra figure 3: honest n / compartment occupancy
    fig, ax = plt.subplots(figsize=(8.2, 3.8))
    post_ord = post.sort_values("Sample")
    x = np.arange(len(post_ord))
    ax.bar(x, post_ord["n_malig_a3"], color=[COLOR[g] for g in post_ord["paper_group"]], alpha=0.9)
    ax.axhline(MIN_MAL_PRIMARY, color="#333", ls="--", lw=0.9, label=f"primary floor n={MIN_MAL_PRIMARY}")
    ax.set_xticks(x)
    ax.set_xticklabels([s.replace("BD_immune", "P") for s in post_ord["Sample"]], rotation=0)
    ax.set_ylabel("A3-malignant cells")
    ax.set_title("Honest n: A3-malignant occupancy (12 post patients)")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    savefig(fig, figdir / "fig_extra_honest_n")

    # Sample-level CLDN4 vs modules
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.0))
    for ax, mod, title in (
        (axes[0], "tj_keratin", "TJ/keratin"),
        (axes[1], "ifn_isg", "IFN ISG"),
    ):
        for grp, sub in post_a3.groupby("paper_group"):
            ax.scatter(sub["mean_CLDN4_a3"], sub[f"mean_{mod}_a3"], c=COLOR[grp], s=48, label=grp)
            for _, r in sub.iterrows():
                ax.annotate(
                    r["Sample"].replace("BD_immune", "P"),
                    (r["mean_CLDN4_a3"], r[f"mean_{mod}_a3"]),
                    textcoords="offset points",
                    xytext=(4, 3),
                    fontsize=7,
                    color=COLOR[grp],
                )
        ax.set_xlabel("malignant mean CLDN4 log1p(CP10k)")
        ax.set_ylabel(f"malignant mean {title}")
        ax.set_title(title)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[1].legend(frameon=False)
    fig.suptitle(f"Sample-level (post, A3-malignant ≥{MIN_MAL_PRIMARY}, n={len(post_a3)})")
    fig.tight_layout()
    savefig(fig, figdir / "fig_sample_cldn4_vs_modules")

    # Compact primary table
    prim_tests = []
    for t in tests:
        if t["kind"] == "paired_high_vs_low" and t.get("compartment") == "a3" and t.get("module") in (
            "tj_keratin",
            "ifn_isg",
            "tj_no_cldn4",
            "keratin",
            "ifn_core6",
            "mhc1_apm",
            "ctrl_oxphos",
            "TACSTD2",
            "total_UMI",
        ):
            prim_tests.append(t)
        if t["kind"] == "sample_spearman" and "post_a3_CLDN4_vs_" in t.get("contrast", ""):
            if any(x in t["contrast"] for x in ("tj_keratin", "ifn_isg", "tj_no_cldn4", "keratin", "ifn_core6", "ctrl_oxphos", "TACSTD2")):
                prim_tests.append(t)
        if t["kind"] == "delta_spearman" and t.get("compartment") == "a3":
            prim_tests.append(t)

    table = pd.DataFrame(prim_tests)
    # Write a human-facing compact TSV
    compact_rows = []
    a3_paired = {t["module"]: t for t in tests if t["kind"] == "paired_high_vs_low" and t.get("compartment") == "a3"}
    for mod, label in (
        ("tj_keratin", "TJ/keratin (CLDN4 excluded)"),
        ("tj_no_cldn4", "TJ only (CLDN4 excluded)"),
        ("keratin", "keratin (simple+basal)"),
        ("ifn_isg", "IFN ISG core"),
        ("ifn_core6", "IFN user core-6"),
        ("mhc1_apm", "MHC-I APM"),
        ("ctrl_oxphos", "control OXPHOS"),
        ("TACSTD2", "companion TACSTD2"),
        ("total_UMI", "library size (total UMI)"),
    ):
        t = a3_paired.get(mod, {})
        compact_rows.append(
            {
                "test": f"paired CLDN4-high vs low, {label}",
                "unit": "patient (within-tumor Q4 vs Q1 malignant cells)",
                "n": t.get("n", 0),
                "n_high_gt_low": t.get("n_high_gt_low", ""),
                "mean_high": t.get("mean_high", np.nan),
                "mean_low": t.get("mean_low", np.nan),
                "median_delta": t.get("median_delta", np.nan),
                "stat": t.get("wilcoxon_w", np.nan),
                "p": t.get("p_value", np.nan),
                "patients": t.get("patients", ""),
            }
        )
    for t in tests:
        if t["kind"] == "sample_spearman" and t["contrast"] in (
            f"post_a3_CLDN4_vs_tj_keratin_min{MIN_MAL_PRIMARY}",
            f"post_a3_CLDN4_vs_ifn_isg_min{MIN_MAL_PRIMARY}",
            "post_a3_CLDN4_vs_TACSTD2",
        ):
            compact_rows.append(
                {
                    "test": t["contrast"],
                    "unit": "patient mean (A3-malignant)",
                    "n": t.get("n", 0),
                    "n_high_gt_low": "",
                    "mean_high": np.nan,
                    "mean_low": np.nan,
                    "median_delta": t.get("spearman_rho", np.nan),
                    "stat": t.get("spearman_rho", np.nan),
                    "p": t.get("spearman_p", np.nan),
                    "patients": "",
                }
            )
        if t["kind"] == "delta_spearman" and t.get("compartment") == "a3":
            compact_rows.append(
                {
                    "test": "paired-delta Spearman TJ/keratin vs IFN ISG",
                    "unit": "patient (same-cell deltas)",
                    "n": t.get("n", 0),
                    "n_high_gt_low": "",
                    "mean_high": np.nan,
                    "mean_low": np.nan,
                    "median_delta": t.get("spearman_rho", np.nan),
                    "stat": t.get("spearman_rho", np.nan),
                    "p": t.get("spearman_p", np.nan),
                    "patients": "",
                }
            )
    compact = pd.DataFrame(compact_rows)

    # Eligibility audit
    elig = paired[paired["compartment"] == "a3"][["Sample", "paper_group", "timing", "n_cells", "n_high", "n_low", "eligible", "reason"]]

    # Summary JSON
    def pick(kind, contrast=None, module=None, compartment="a3"):
        for t in tests:
            if t["kind"] != kind:
                continue
            if contrast and t.get("contrast") != contrast:
                continue
            if module and t.get("module") != module:
                continue
            if kind in ("paired_high_vs_low", "delta_spearman") and t.get("compartment") != compartment:
                continue
            return t
        return {}

    paired_tj = pick("paired_high_vs_low", module="tj_keratin")
    paired_ifn = pick("paired_high_vs_low", module="ifn_isg")
    dropped_post = elig[(elig["timing"] == "post") & (~elig["eligible"])]
    summary = {
        "dataset": "GSE207422",
        "n_cells_matrix": int(n),
        "n_genes_matrix": n_genes_in_matrix,
        "n_epithelial": int(is_epi.sum()),
        "n_malig_a3": int(is_malig.sum()),
        "n_post_patients": 12,
        "n_post_a3_eligible_primary": int(len(prim)),
        "eligible_patients": prim["Sample"].tolist(),
        "dropped_post_a3": dropped_post.to_dict(orient="records"),
        "missing_genes": missing,
        "module_genes_present": {k: v for k, v in module_present.items()},
        "primary_paired_tj_keratin": paired_tj,
        "primary_paired_ifn_isg": paired_ifn,
        "min_mal_primary": MIN_MAL_PRIMARY,
        "min_tail": MIN_TAIL,
        "note": (
            "Unit is the patient. Cell-level Spearman is exploratory. "
            "CLDN4 is excluded from the TJ module. Dual-high was not run."
        ),
    }

    patients.to_csv(outdir / "per_patient.tsv", sep="\t", index=False)
    paired.to_csv(outdir / "paired_high_low.tsv", sep="\t", index=False)
    tests_df.to_csv(outdir / "tests.tsv", sep="\t", index=False)
    compact.to_csv(outdir / "primary_table.tsv", sep="\t", index=False)
    elig.to_csv(outdir / "eligibility.tsv", sep="\t", index=False)
    lineage_counts.to_csv(outdir / "lineage_counts.tsv", sep="\t", index=False)
    meta.to_csv(outdir / "sample_metadata.tsv", sep="\t", index=False)
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")

    # Do not write 92k-row per-cell table into git; write a malignant-only slim table
    slim_cols = [
        "cell", "Sample", "paper_group", "timing", "lineage",
        "is_malig_a3", "is_epi", "total_umi", "n_genes",
        "CLDN4_log1p_cp10k", "TACSTD2_log1p_cp10k",
        "mod_tj_keratin", "mod_ifn_isg", "mod_tj_no_cldn4", "mod_keratin",
        "mod_ifn_core6", "mod_mhc1_apm", "mod_ctrl_oxphos",
    ]
    mal_slim = per_cell.loc[per_cell["is_malig_a3"] | per_cell["is_epi"], slim_cols]
    mal_slim.to_csv(outdir / "malignant_epithelial_scores.tsv.gz", sep="\t", index=False, compression="gzip")

    print("=== primary paired A3 (post) ===", flush=True)
    print(f"eligible n={len(prim)} patients: {prim['Sample'].tolist()}", flush=True)
    print("TJ/keratin:", fmt_paired(paired_tj), flush=True)
    print("IFN ISG:   ", fmt_paired(paired_ifn), flush=True)
    print(f"wrote tables under {outdir}", flush=True)


if __name__ == "__main__":
    main()
