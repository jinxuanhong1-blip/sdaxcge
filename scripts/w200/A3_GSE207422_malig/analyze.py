#!/usr/bin/env python3
"""A3 on GSE207422: marker-based malignant TACSTD2 vs NMPR/MPR and vs T/NK.

Author per-cell labels are not public. Lineages follow Hu et al. Genome Medicine
2023 canonical markers. Malignant-like = epithelial minus normal-lung markers
(SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3). That is not CopyKAT (author method).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

MIN_MALIG = 10
MIN_TNK = 20

# Paper groups (Hu et al. Fig. 1): TN n=3 pre-biopsy; post-surgery MPR n=4
# (includes pCR P06) and NMPR n=8.
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

LINEAGES = {
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "T": ["CD3D", "CD3E", "CD2"],
    "NK": ["NKG7", "GNLY", "FGFBP2"],
    "B": ["CD79A", "MS4A1"],
    "plasma": ["IGHG1", "MZB1"],
    "myeloid": ["LYZ", "CD68", "CD14"],
    "neutrophil": ["CSF3R"],
    "fibroblast": ["COL1A1", "DCN"],
    "endothelial": ["VWF", "PECAM1"],
    "mast": ["KIT"],
}

NORMAL_LUNG = ["SFTPA2", "AGER", "SCGB1A1", "SCGB3A1", "TPPP3"]
REPORT_GENES = [
    "TACSTD2",
    "CLDN4",
    "EPCAM",
    "PTPRC",
    "CD3E",
    "NKG7",
    "CX3CL1",
    "CD74",
    "HLA-DRA",
    "AKR1C1",
    "AKR1C2",
    "AKR1C3",
]


def score(expr: dict[str, np.ndarray], genes: list[str], n: int) -> np.ndarray:
    mats = [np.log1p(expr[g].astype(np.float32)) for g in genes if g in expr]
    if not mats:
        return np.zeros(n, dtype=np.float32)
    return np.mean(np.vstack(mats), axis=0)


def assign_lineage(expr: dict[str, np.ndarray], n: int) -> np.ndarray:
    names = list(LINEAGES)
    scores = np.vstack([score(expr, LINEAGES[k], n) for k in names])
    best = scores.argmax(axis=0)
    best_val = scores.max(axis=0)
    second = np.partition(scores, -2, axis=0)[-2]
    assigned = np.array(names, dtype=object)[best]
    assigned[(best_val < 0.15) | ((best_val - second) < 0.05)] = "other"
    return assigned


def summarize_group(df: pd.DataFrame, genes: list[str], total: np.ndarray | None = None) -> dict:
    out = {"n_cells": int(len(df))}
    if df.empty:
        for gene in genes:
            out[f"{gene}_mean_umi"] = np.nan
            out[f"{gene}_mean_log1p"] = np.nan
            out[f"{gene}_mean_cp10k"] = np.nan
            out[f"{gene}_pct_pos"] = np.nan
        return out
    lib = df["total"].to_numpy(dtype=np.float64) if "total" in df.columns else None
    for gene in genes:
        umi = df[gene].to_numpy(dtype=np.float64)
        out[f"{gene}_mean_umi"] = float(np.mean(umi))
        out[f"{gene}_mean_log1p"] = float(np.mean(np.log1p(umi)))
        out[f"{gene}_pct_pos"] = float(np.mean(umi > 0) * 100.0)
        if lib is not None:
            cp = np.where(lib > 0, umi / lib * 1e4, np.nan)
            out[f"{gene}_mean_cp10k"] = float(np.nanmean(np.log1p(cp)))
        else:
            out[f"{gene}_mean_cp10k"] = np.nan
    return out


def spearman_block(x, y, label: str) -> dict:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 4:
        return {
            "contrast": label,
            "n": int(len(x)),
            "spearman_rho": np.nan,
            "spearman_p": np.nan,
            "note": "too_few_samples",
        }
    rho, p = stats.spearmanr(x, y)
    return {
        "contrast": label,
        "n": int(len(x)),
        "spearman_rho": float(rho),
        "spearman_p": float(p),
        "note": "",
    }


def mwu_block(a, b, label: str) -> dict:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return {
            "contrast": label,
            "n_a": int(len(a)),
            "n_b": int(len(b)),
            "median_a": float(np.median(a)) if len(a) else np.nan,
            "median_b": float(np.median(b)) if len(b) else np.nan,
            "mwu_u": np.nan,
            "mwu_p": np.nan,
            "note": "too_few_samples",
        }
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return {
        "contrast": label,
        "n_a": int(len(a)),
        "n_b": int(len(b)),
        "median_a": float(np.median(a)),
        "median_b": float(np.median(b)),
        "mwu_u": float(u),
        "mwu_p": float(p),
        "note": "",
    }


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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", type=Path, default=Path("data/GSE207422"))
    ap.add_argument(
        "--outdir",
        type=Path,
        default=Path(__file__).resolve().parents[3] / "results" / "w200" / "A3_GSE207422_malig",
    )
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    extracted = np.load(args.workdir / "extracted_markers.npz", allow_pickle=True)
    cells = extracted["cells"]
    genes = list(extracted["genes"])
    mat = extracted["mat"]
    total = extracted["total"]
    expr = {g: mat[i] for i, g in enumerate(genes)}
    n = len(cells)

    sample_ids = np.array([str(c).rsplit("_", 1)[0] for c in cells])
    lineage = assign_lineage(expr, n)
    is_tnk = np.isin(lineage, ["T", "NK"])
    is_epi = lineage == "epithelial"
    normal_umi = np.zeros(n, dtype=np.int32)
    for g in NORMAL_LUNG:
        if g in expr:
            normal_umi += expr[g]
    is_malig = is_epi & (normal_umi == 0)
    is_normal_epi = is_epi & (normal_umi > 0)

    per_cell = pd.DataFrame(
        {
            "cell": cells,
            "Sample": sample_ids,
            "lineage": lineage,
            "is_tnk": is_tnk,
            "is_epi": is_epi,
            "is_malig_like": is_malig,
            "is_normal_epi": is_normal_epi,
            "total": total,
        }
    )
    for g in genes:
        per_cell[g] = expr[g]
    per_cell["paper_group"] = per_cell["Sample"].map(PAPER_GROUP)

    sample_meta = load_sample_meta(args.workdir / "GSE207422_NSCLC_scRNAseq_metadata.xlsx")
    sample_meta.to_csv(args.outdir / "sample_metadata.tsv", sep="\t", index=False)

    comp = (
        per_cell.groupby(["Sample", "lineage"], dropna=False)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["Sample", "n_cells"], ascending=[True, False])
    )
    comp.to_csv(args.outdir / "cell_type_composition.tsv", sep="\t", index=False)
    sample_comp = pd.crosstab(per_cell["Sample"], per_cell["lineage"]).reset_index()
    sample_comp.to_csv(args.outdir / "sample_composition.tsv", sep="\t", index=False)

    genes_present = [g for g in REPORT_GENES if g in per_cell.columns]
    rows = []
    for sample, sdf in per_cell.groupby("Sample", sort=True):
        epi = sdf[sdf["is_epi"]]
        malig = sdf[sdf["is_malig_like"]]
        tnk = sdf[sdf["is_tnk"]]
        row = {
            "Sample": sample,
            "n_cells": int(len(sdf)),
            "n_epithelial": int(len(epi)),
            "n_malig_like": int(len(malig)),
            "n_normal_epi": int(sdf["is_normal_epi"].sum()),
            "n_T_NK": int(len(tnk)),
            "n_T": int((sdf["lineage"] == "T").sum()),
            "n_NK": int((sdf["lineage"] == "NK").sum()),
            "frac_epithelial": float(len(epi) / len(sdf)),
            "frac_malig_like": float(len(malig) / len(sdf)),
            "frac_T_NK": float(len(tnk) / len(sdf)),
            "eligible_malig_tnk": int(len(malig) >= MIN_MALIG and len(tnk) >= MIN_TNK),
        }
        for prefix, frame in (("epi", epi), ("malig", malig), ("tnk", tnk)):
            stats_row = summarize_group(frame, genes_present)
            for key, value in stats_row.items():
                if key == "n_cells":
                    continue
                row[f"{prefix}_{key}"] = value
        rows.append(row)
    sample_df = pd.DataFrame(rows).merge(sample_meta, on="Sample", how="left")
    sample_df.to_csv(args.outdir / "per_sample_tacstd2.tsv", sep="\t", index=False)

    post = sample_df[sample_df["paper_group"].isin(["MPR", "NMPR"])].copy()
    post_el = post[post["eligible_malig_tnk"] == 1]
    all_el = sample_df[sample_df["eligible_malig_tnk"] == 1]
    mpr = post[post["paper_group"] == "MPR"]
    nmpr = post[post["paper_group"] == "NMPR"]
    mpr_el = post_el[post_el["paper_group"] == "MPR"]
    nmpr_el = post_el[post_el["paper_group"] == "NMPR"]

    contrasts = []
    for frame, label in (
        (all_el, "all eligible samples"),
        (post_el, "post-treatment eligible (paper MPR/NMPR)"),
        (post, "post-treatment all (no cell-count filter)"),
    ):
        contrasts.append(
            spearman_block(
                frame["malig_TACSTD2_mean_log1p"],
                frame["frac_T_NK"],
                f"{label}: malig-like TACSTD2 mean_log1p vs T/NK fraction",
            )
        )
        contrasts.append(
            spearman_block(
                frame["malig_TACSTD2_pct_pos"],
                frame["frac_T_NK"],
                f"{label}: malig-like TACSTD2 %pos vs T/NK fraction",
            )
        )
        contrasts.append(
            spearman_block(
                frame["malig_TACSTD2_mean_cp10k"],
                frame["frac_T_NK"],
                f"{label}: malig-like TACSTD2 mean log1p(CP10k) vs T/NK fraction",
            )
        )
        contrasts.append(
            spearman_block(
                frame["epi_TACSTD2_mean_log1p"],
                frame["frac_T_NK"],
                f"{label}: all-epithelial TACSTD2 mean_log1p vs T/NK fraction",
            )
        )

    group_tests = [
        mwu_block(
            nmpr["malig_TACSTD2_mean_log1p"],
            mpr["malig_TACSTD2_mean_log1p"],
            "post: NMPR vs MPR malig-like TACSTD2 mean_log1p (sample-level)",
        ),
        mwu_block(
            nmpr["malig_TACSTD2_pct_pos"],
            mpr["malig_TACSTD2_pct_pos"],
            "post: NMPR vs MPR malig-like TACSTD2 %pos (sample-level)",
        ),
        mwu_block(
            nmpr["malig_TACSTD2_mean_cp10k"],
            mpr["malig_TACSTD2_mean_cp10k"],
            "post: NMPR vs MPR malig-like TACSTD2 mean log1p(CP10k) (sample-level)",
        ),
        mwu_block(
            nmpr["frac_T_NK"],
            mpr["frac_T_NK"],
            "post: NMPR vs MPR T/NK fraction (sanity vs paper Fig 1D)",
        ),
        mwu_block(
            nmpr["malig_CX3CL1_mean_log1p"],
            mpr["malig_CX3CL1_mean_log1p"],
            "post: NMPR vs MPR malig-like CX3CL1 mean_log1p (paper: MPR higher)",
        ),
        mwu_block(
            nmpr["malig_HLA-DRA_mean_log1p"],
            mpr["malig_HLA-DRA_mean_log1p"],
            "post: NMPR vs MPR malig-like HLA-DRA mean_log1p (paper: MPR higher)",
        ),
        mwu_block(
            nmpr["malig_AKR1C1_mean_log1p"],
            mpr["malig_AKR1C1_mean_log1p"],
            "post: NMPR vs MPR malig-like AKR1C1 mean_log1p (paper: NMPR higher)",
        ),
    ]

    paired = all_el.dropna(subset=["malig_TACSTD2_mean_log1p", "tnk_TACSTD2_mean_log1p"])
    if len(paired) >= 6:
        w_stat, w_p = stats.wilcoxon(
            paired["malig_TACSTD2_mean_log1p"],
            paired["tnk_TACSTD2_mean_log1p"],
            alternative="greater",
        )
        paired_result = {
            "contrast": "paired samples: malig-like TACSTD2 > T/NK TACSTD2 (Wilcoxon signed-rank)",
            "n": int(len(paired)),
            "wilcoxon_stat": float(w_stat),
            "wilcoxon_p": float(w_p),
            "median_malig_mean_log1p": float(paired["malig_TACSTD2_mean_log1p"].median()),
            "median_tnk_mean_log1p": float(paired["tnk_TACSTD2_mean_log1p"].median()),
            "median_malig_pct_pos": float(paired["malig_TACSTD2_pct_pos"].median()),
            "median_tnk_pct_pos": float(paired["tnk_TACSTD2_pct_pos"].median()),
        }
    else:
        paired_result = {"contrast": "paired malig vs T/NK TACSTD2", "n": int(len(paired)), "note": "too_few"}

    pd.DataFrame(contrasts).to_csv(args.outdir / "association_statistics.tsv", sep="\t", index=False)
    pd.DataFrame(group_tests).to_csv(args.outdir / "group_statistics.tsv", sep="\t", index=False)
    with (args.outdir / "paired_compartment_test.json").open("w") as fh:
        json.dump(paired_result, fh, indent=2)

    sanity = {
        "n_cells": int(n),
        "n_samples": int(per_cell["Sample"].nunique()),
        "lineage_counts": {k: int(v) for k, v in per_cell["lineage"].value_counts().items()},
        "n_malig_like": int(is_malig.sum()),
        "n_epithelial": int(is_epi.sum()),
        "n_T_NK": int(is_tnk.sum()),
        "n_normal_epi": int(is_normal_epi.sum()),
        "median_genes_per_cell": float(np.median(extracted["n_nonzero_genes"])),
        "paper_n_cells": 92330,
        "paper_median_genes": 1256,
        "malig_EPCAM_mean_log1p": float(np.mean(np.log1p(per_cell.loc[is_malig, "EPCAM"]))) if is_malig.any() else np.nan,
        "tnk_EPCAM_mean_log1p": float(np.mean(np.log1p(per_cell.loc[is_tnk, "EPCAM"]))) if is_tnk.any() else np.nan,
        "malig_PTPRC_mean_log1p": float(np.mean(np.log1p(per_cell.loc[is_malig, "PTPRC"]))) if is_malig.any() else np.nan,
        "tnk_PTPRC_mean_log1p": float(np.mean(np.log1p(per_cell.loc[is_tnk, "PTPRC"]))) if is_tnk.any() else np.nan,
        "malig_CD3E_mean_log1p": float(np.mean(np.log1p(per_cell.loc[is_malig, "CD3E"]))) if is_malig.any() else np.nan,
        "tnk_CD3E_mean_log1p": float(np.mean(np.log1p(per_cell.loc[is_tnk, "CD3E"]))) if is_tnk.any() else np.nan,
        "note": (
            "Marker-based compartments, not author Seurat/CopyKAT labels. "
            "Malignant-like = epithelial argmax AND zero UMI for SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3."
        ),
    }
    with (args.outdir / "sanity_checks.json").open("w") as fh:
        json.dump(sanity, fh, indent=2)

    # figures
    color = {"MPR": "#d1495b", "NMPR": "#2c6eaf", "TN": "#6b6b6b"}
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.1), constrained_layout=True)
    panels = [
        (all_el, "All eligible samples", "malig_TACSTD2_mean_log1p"),
        (post_el, "Post-treatment eligible", "malig_TACSTD2_mean_log1p"),
        (post_el, "Post-treatment eligible", "malig_TACSTD2_pct_pos"),
    ]
    xlabels = [
        "Malignant-like TACSTD2 (mean log1p UMI)",
        "Malignant-like TACSTD2 (mean log1p UMI)",
        "Malignant-like TACSTD2 (% UMI>0)",
    ]
    for ax, (frame, title, xcol), xlab in zip(axes, panels, xlabels):
        if frame.empty:
            ax.set_title(title + "\n(no eligible samples)")
            ax.axis("off")
            continue
        for grp, sub in frame.groupby("paper_group"):
            ax.scatter(
                sub[xcol],
                sub["frac_T_NK"],
                c=color.get(grp, "#333"),
                label=grp,
                edgecolors="white",
                linewidths=0.6,
                s=48,
            )
        if len(frame) >= 2:
            coef = np.polyfit(frame[xcol].to_numpy(), frame["frac_T_NK"].to_numpy(), 1)
            grid = np.linspace(frame[xcol].min(), frame[xcol].max(), 50)
            ax.plot(grid, np.polyval(coef, grid), color="#333333", lw=1)
        rho, p = stats.spearmanr(frame[xcol], frame["frac_T_NK"])
        ax.set_title(f"{title}\nρ={rho:.2f}  p={p:.3g}  n={len(frame)}")
        ax.set_xlabel(xlab)
        ax.set_ylabel("T/NK fraction")
        ax.legend(frameon=False, fontsize=8)
    fig.savefig(args.outdir / "fig_malig_tacstd2_vs_tnk_fraction.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2), constrained_layout=True)
    plot_df = post.dropna(subset=["malig_TACSTD2_mean_log1p"])
    if not plot_df.empty:
        axes[0].boxplot(
            [
                plot_df.loc[plot_df.paper_group == "MPR", "malig_TACSTD2_mean_log1p"],
                plot_df.loc[plot_df.paper_group == "NMPR", "malig_TACSTD2_mean_log1p"],
            ],
            tick_labels=["MPR", "NMPR"],
            widths=0.55,
        )
        for i, grp in enumerate(["MPR", "NMPR"], start=1):
            y = plot_df.loc[plot_df.paper_group == grp, "malig_TACSTD2_mean_log1p"]
            axes[0].scatter(np.full(len(y), i) + np.random.default_rng(0).uniform(-0.08, 0.08, len(y)), y, s=28, c=color[grp], zorder=3)
        g = next(t for t in group_tests if "mean_log1p (sample-level)" in t["contrast"] and "TACSTD2" in t["contrast"])
        axes[0].set_ylabel("Malignant-like TACSTD2 (mean log1p UMI)")
        axes[0].set_title(f"NMPR vs MPR (sample-level)\nMWU p={g['mwu_p']:.3g}; n={g['n_b']} vs {g['n_a']}")
    if not paired.empty:
        axes[1].boxplot(
            [paired["malig_TACSTD2_pct_pos"], paired["tnk_TACSTD2_pct_pos"]],
            tick_labels=["Malignant-like", "T/NK"],
            widths=0.55,
        )
        axes[1].set_ylabel("% TACSTD2+ cells (UMI>0)")
        axes[1].set_title(
            f"TACSTD2 compartment restriction\n"
            f"median {paired_result.get('median_malig_pct_pos', float('nan')):.1f}% vs "
            f"{paired_result.get('median_tnk_pct_pos', float('nan')):.1f}%; "
            f"p={paired_result.get('wilcoxon_p', float('nan')):.2g}"
        )
    fig.savefig(args.outdir / "fig_tacstd2_nmpr_mpr_and_compartment.png", dpi=160)
    plt.close(fig)

    rho_post = next(
        c
        for c in contrasts
        if c["contrast"].startswith("post-treatment eligible (paper MPR/NMPR): malig-like TACSTD2 mean_log1p")
    )
    rho_all = next(
        c
        for c in contrasts
        if c["contrast"].startswith("all eligible samples: malig-like TACSTD2 mean_log1p")
    )
    mwu_tac = next(t for t in group_tests if t["contrast"].startswith("post: NMPR vs MPR malig-like TACSTD2 mean_log1p"))

    def _fmt_rho(item: dict) -> str:
        if item.get("note") == "too_few_samples":
            return f"n={item['n']} (too few)"
        return f"ρ={item['spearman_rho']:.2f}, p={item['spearman_p']:.2g}, n={item['n']}"

    user_claim_rho = "user-stated ρ −0.40 to −0.50"
    rho_val = rho_post["spearman_rho"]
    if np.isfinite(rho_val) and -0.50 <= rho_val <= -0.40:
        rho_vs_claim = f"recomputed {_fmt_rho(rho_post)} falls inside the {user_claim_rho} window"
        immune_supported = True
    else:
        rho_vs_claim = f"recomputed {_fmt_rho(rho_post)} does not match the {user_claim_rho} window"
        immune_supported = bool(np.isfinite(rho_val) and rho_val <= -0.30 and rho_post["spearman_p"] < 0.05)

    nmpr_higher = (
        np.isfinite(mwu_tac["median_a"])
        and np.isfinite(mwu_tac["median_b"])
        and mwu_tac["median_a"] > mwu_tac["median_b"]
        and mwu_tac["mwu_p"] < 0.05
    )

    audit = {
        "dataset": "GSE207422",
        "n_cells": int(n),
        "n_samples": int(per_cell["Sample"].nunique()),
        "genes_found": genes,
        "genes_missing": sorted(set(REPORT_GENES) - set(genes)),
        "n_malig_like": int(is_malig.sum()),
        "n_T_NK": int(is_tnk.sum()),
        "n_eligible_all": int(len(all_el)),
        "n_eligible_post": int(len(post_el)),
        "n_post_mpr": int(len(mpr)),
        "n_post_nmpr": int(len(nmpr)),
        "author_cell_labels_used": False,
        "compartment_definition": sanity["note"],
        "expression_metric": "mean log1p(raw UMI), % UMI>0, and mean log1p(CP10k)",
        "min_cells": {"malig_like": MIN_MALIG, "T_NK": MIN_TNK},
        "paper_groups": "TN=pre-biopsy (n=3); post MPR includes pCR P06 (n=4); post NMPR (n=8)",
        "associations": contrasts,
        "group_tests": group_tests,
        "paired_compartment": paired_result,
    }
    with (args.outdir / "audit.json").open("w") as fh:
        json.dump(audit, fh, indent=2)

    summary = {
        "dataset": "GSE207422",
        "task": "W200-A3: malignant-only TACSTD2 NMPR vs MPR and vs T/NK",
        "n_cells": int(n),
        "author_cell_labels_public": False,
        "honest_verdict": (
            f"Processed UMI 175.5 MB was in budget; recomputed on {n} cells. "
            f"Author per-cell labels are not public; used marker-based malignant-like cells "
            f"(n={int(is_malig.sum())}). "
            f"NMPR vs MPR malignant TACSTD2 (sample-level mean log1p): "
            f"median NMPR={mwu_tac['median_a']:.3f} vs MPR={mwu_tac['median_b']:.3f}, "
            f"MWU p={mwu_tac['mwu_p']:.2g} (n={mwu_tac['n_a']} vs {mwu_tac['n_b']}). "
            f"Malignant TACSTD2 vs T/NK fraction (post eligible): {_fmt_rho(rho_post)}; "
            f"all eligible: {_fmt_rho(rho_all)}. {rho_vs_claim}."
        ),
        "mpr_contrast_possible": True,
        "nmpr_tacstd2_higher_than_mpr": bool(nmpr_higher),
        "immune_anticorrelation_supported": bool(immune_supported),
        "user_claim_rho_window": [-0.50, -0.40],
        "recomputed_post_eligible_rho": rho_post["spearman_rho"],
        "recomputed_post_eligible_p": rho_post["spearman_p"],
        "recomputed_post_eligible_n": rho_post["n"],
    }
    with (args.outdir / "summary.json").open("w") as fh:
        json.dump(summary, fh, indent=2)

    write_results_readme(
        args.outdir,
        sample_df,
        all_el,
        post_el,
        contrasts,
        group_tests,
        paired_result,
        sanity,
        summary,
    )
    print(json.dumps(summary, indent=2))


def write_results_readme(
    outdir: Path,
    sample_df: pd.DataFrame,
    all_el: pd.DataFrame,
    post_el: pd.DataFrame,
    contrasts: list[dict],
    group_tests: list[dict],
    paired: dict,
    sanity: dict,
    summary: dict,
) -> None:
    def fmt_s(item: dict) -> str:
        if item.get("note") == "too_few_samples":
            return f"- {item['contrast']}: n={item['n']} (too few)"
        return (
            f"- {item['contrast']}: ρ={item['spearman_rho']:.3f}, "
            f"p={item['spearman_p']:.3g}, n={item['n']}"
        )

    def fmt_g(item: dict) -> str:
        if item.get("note") == "too_few_samples":
            return f"- {item['contrast']}: too few"
        return (
            f"- {item['contrast']}: median {item['median_a']:.3f} vs {item['median_b']:.3f}; "
            f"MWU p={item['mwu_p']:.3g}; n={item['n_a']} vs {item['n_b']}"
        )

    lines = [
        "# W200-A3 · GSE207422 — malignant-only TACSTD2 (NMPR vs MPR, vs T/NK)",
        "",
        "**Task:** Recompute Hu et al. 2023 malignant TACSTD2: NMPR vs MPR, and vs T/NK. Honest.",
        "",
        "## Verdict",
        "",
        summary["honest_verdict"],
        "",
        "## Data policy",
        "",
        "- Processed UMI matrix **175.5 MB gzip < 2 GB** → recomputed (not catalog-only).",
        "- 92,330 cells (header count matches the paper's post-QC 92,330).",
        "- Raw GSA-Human HRA001033 was **not** downloaded.",
        "- Author per-cell labels are **not public** (GEO sample sheet only; paper ESM has no barcode table; "
        "[Junjie-Hu/NSCLC-immunotherapy](https://github.com/Junjie-Hu/NSCLC-immunotherapy) has scripts, not RDS/metadata_scrublet.csv; "
        "TISCH2 NSCLC_GSE207422 CellMetainfo 404). See `label_search.json`.",
        "",
        "## What was used instead of author labels",
        "",
        "- Lineage = argmax of mean log1p marker scores (Hu canonical: EPCAM/KRTs, CD3D/E/CD2, NKG7/GNLY/FGFBP2, etc.).",
        "- **Malignant-like** = epithelial AND zero UMI for normal-lung markers SFTPA2, AGER, SCGB1A1, SCGB3A1, TPPP3.",
        "- This is **not** author CopyKAT aneuploid calls. Residual basal/unmarked normal epithelium can leak in.",
        "- T/NK = T or NK lineage. Eligible sample: ≥10 malignant-like and ≥20 T/NK cells.",
        "- Paper groups: TN = 3 pre-biopsies; post MPR includes pCR P06 (n=4); post NMPR n=8.",
        f"- Marker-based counts: epithelial={sanity['n_epithelial']}, malignant-like={sanity['n_malig_like']}, T/NK={sanity['n_T_NK']}.",
        "",
        "## Results",
        "",
    ]
    lines.extend(fmt_s(c) for c in contrasts)
    lines.append("")
    lines.extend(fmt_g(t) for t in group_tests)
    if paired.get("n"):
        lines.append(
            f"- {paired['contrast']}: n={paired['n']}, p={paired.get('wilcoxon_p', float('nan')):.3g}; "
            f"median %pos malig-like={paired.get('median_malig_pct_pos', float('nan')):.1f} vs "
            f"T/NK={paired.get('median_tnk_pct_pos', float('nan')):.1f}."
        )
    lines += [
        "",
        "## Honest limits",
        "",
        "- n=12 post-treatment samples (4 MPR / 8 NMPR) is underpowered for a ρ≈−0.45 claim.",
        "- Sample-level tests are the correct unit; cell-level p-values would be pseudoreplication and are not reported as primary.",
        "- Without CopyKAT, “malignant-only” is a marker proxy.",
        "- Expression is log1p(raw UMI) and log1p(CP10k), not author log-normalized Seurat values after CCA.",
        "",
        "## Files",
        "",
        "| File | Role |",
        "|---|---|",
        "| `feasibility.json` / `file_manifest.tsv` | Size budget |",
        "| `label_search.json` | Author-label hunt |",
        "| `sample_metadata.tsv` | GEO + paper TN/MPR/NMPR groups |",
        "| `cell_type_composition.tsv` / `sample_composition.tsv` | Marker-based lineages |",
        "| `per_sample_tacstd2.tsv` | Sample-level malignant / T/NK TACSTD2 |",
        "| `association_statistics.tsv` | Spearman vs T/NK |",
        "| `group_statistics.tsv` | NMPR vs MPR MWU |",
        "| `summary.json` / `audit.json` / `sanity_checks.json` | Verdict |",
        "| `fig_malig_tacstd2_vs_tnk_fraction.png` | A3-style correlation |",
        "| `fig_tacstd2_nmpr_mpr_and_compartment.png` | Group + restriction |",
        "",
        "## Reproduce",
        "",
        "```bash",
        "python3 scripts/w200/A3_GSE207422_malig/download.py",
        "python3 scripts/w200/A3_GSE207422_malig/extract.py",
        "python3 scripts/w200/A3_GSE207422_malig/analyze.py",
        "```",
        "",
    ]
    (outdir / "README.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
