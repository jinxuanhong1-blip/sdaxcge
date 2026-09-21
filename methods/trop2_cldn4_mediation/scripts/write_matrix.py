#!/usr/bin/env python3
"""Write the one-page evidence matrix from the analysis tables."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"


def fmt_p(p) -> str:
    p = float(p)
    if not np.isfinite(p):
        return "NA"
    if p <= 0:
        return "p<0.001"
    if p < 1e-3:
        return f"{p:.1e}"
    return f"{p:.3f}"


def fmt_rho(r, p=None) -> str:
    if not np.isfinite(r):
        return "—"
    text = f"{float(r):+.3f}"
    if p is not None and np.isfinite(float(p)):
        ps = fmt_p(p)
        text += f" ({ps})" if ps.startswith("p<") else f" (p={ps})"
    return text


def fmt_att(a) -> str:
    try:
        val = float(a)
    except (TypeError, ValueError):
        return "—"
    if not np.isfinite(val):
        return "—"
    return f"{val:+.0f}%"


def fmt_acme(a, lo, hi) -> str:
    if not np.isfinite(a):
        return "—"
    if np.isfinite(lo) and np.isfinite(hi):
        return f"{float(a):+.3f} ({float(lo):+.3f}, {float(hi):+.3f})"
    return f"{float(a):+.3f}"


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    cohort = pd.read_csv(TAB / "concordant4_cohort.tsv", sep="\t", keep_default_na=False)
    meta = pd.read_csv(TAB / "concordant4_meta.tsv", sep="\t", keep_default_na=False)
    cos = pd.read_csv(TAB / "cosmx_summary.tsv", sep="\t", keep_default_na=False)
    between = pd.read_csv(TAB / "cosmx_between_section.tsv", sep="\t", keep_default_na=False)

    rows = []
    order = [("donor", "cd8nk_50um"), ("donor", "cd8nk_100um")]
    labels = {
        "cd8nk_50um": "CD8+NK neighbors, 50 µm",
        "cd8nk_100um": "CD8+NK neighbors, 100 µm",
        "pct_pos": "malignant % detected",
        "mean_log1p_umi": "malignant mean log1p(UMI)",
    }
    for unit, readout in order:
        r = cos[(cos["unit"] == unit) & (cos.readout == readout)].iloc[0]
        rows.append({
            "dataset": "He2022 CosMx",
            "readout": labels[readout],
            "n": "5 donors",
            "tac_raw": fmt_rho(r.rho_tacstd2, r.rho_tacstd2_p),
            "cldn_raw": fmt_rho(r.rho_cldn4, r.rho_cldn4_p),
            "tac_partial": fmt_rho(r.rho_tacstd2_given_cldn4, r.rho_tacstd2_given_cldn4_p),
            "cldn_partial": fmt_rho(r.rho_cldn4_given_tacstd2, r.rho_cldn4_given_tacstd2_p),
            "attenuation": fmt_att(r.attenuation_of_mean_rho),
            "acme": fmt_acme(r.acme, r.acme_lo, r.acme_hi),
            "verdict": r.verdict,
        })
    cohort_order = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
    readout_order = ["pct_pos", "mean_log1p_umi"]
    for ds in cohort_order:
        for readout in readout_order:
            r = cohort[(cohort.dataset == ds) & (cohort.readout == readout)].iloc[0]
            rows.append({
                "dataset": ds,
                "readout": labels[readout],
                "n": str(int(r.n)),
                "tac_raw": fmt_rho(r.rho_tacstd2, r.p_tacstd2),
                "cldn_raw": fmt_rho(r.rho_cldn4, r.p_cldn4),
                "tac_partial": fmt_rho(r.rho_tacstd2_given_cldn4, r.p_tacstd2_given_cldn4),
                "cldn_partial": fmt_rho(r.rho_cldn4_given_tacstd2, r.p_cldn4_given_tacstd2),
                "attenuation": fmt_att(r.attenuation_pct),
                "acme": fmt_acme(r.acme, r.acme_lo, r.acme_hi),
                "verdict": r.verdict,
            })
    for readout in readout_order:
        r = meta[meta.readout == readout].iloc[0]
        rows.append({
            "dataset": "concordant-4",
            "readout": labels[readout],
            "n": "65 (4 cohorts)",
            "tac_raw": fmt_rho(r.rho_tacstd2, r.p_tacstd2),
            "cldn_raw": fmt_rho(r.rho_cldn4, r.p_cldn4),
            "tac_partial": fmt_rho(r.rho_tacstd2_given_cldn4, r.p_tacstd2_given_cldn4),
            "cldn_partial": fmt_rho(r.rho_cldn4_given_tacstd2, r.p_cldn4_given_tacstd2),
            "attenuation": fmt_att(r.attenuation_pct),
            "acme": fmt_acme(r.acme, r.acme_lo, r.acme_hi),
            "verdict": r.verdict,
        })

    matrix = pd.DataFrame(rows)
    matrix.to_csv(TAB / "evidence_matrix.tsv", sep="\t", index=False)

    header = [
        "Dataset", "Immune readout", "n",
        "TACSTD2–immune", "CLDN4–immune",
        "TACSTD2 given CLDN4", "CLDN4 given TACSTD2",
        "% attenuation", "ACME (95% interval)", "Call",
    ]
    keys = ["dataset", "readout", "n", "tac_raw", "cldn_raw", "tac_partial", "cldn_partial", "attenuation", "acme", "verdict"]
    lines = [
        "| " + " | ".join(header) + " |",
        "| " + " | ".join(["---"] * len(header)) + " |",
    ]
    for rec in rows:
        lines.append("| " + " | ".join(str(rec[k]) for k in keys) + " |")
    table_md = "\n".join(lines) + "\n"
    (TAB / "evidence_matrix.md").write_text(table_md)

    # Figure
    colors = {"null": "#f4f4f4", "partial": "#fff4d6", "opposite": "#fde8e8", "support": "#e5f6e8"}
    fig_h = 1.6 + 0.42 * len(rows)
    fig, ax = plt.subplots(figsize=(18.2, fig_h))
    ax.axis("off")
    cell = [[rec[k] for k in keys] for rec in rows]
    table = ax.table(cellText=cell, colLabels=header, loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.55)
    widths = [0.10, 0.145, 0.075, 0.115, 0.115, 0.125, 0.13, 0.075, 0.145, 0.055]
    for (r, c), cell_obj in table.get_celld().items():
        cell_obj.set_width(widths[c])
    for (r, c), cell_obj in table.get_celld().items():
        cell_obj.set_edgecolor("#cccccc")
        cell_obj.set_linewidth(0.4)
        if r == 0:
            cell_obj.set_facecolor("#1f2933")
            cell_obj.set_text_props(color="white", fontweight="bold")
            continue
        verdict = rows[r - 1]["verdict"]
        if c == 9:
            cell_obj.set_facecolor(colors.get(verdict, "white"))
            cell_obj.set_text_props(fontweight="bold")
        elif r % 2 == 0:
            cell_obj.set_facecolor("#fafafa")
    ax.set_title(
        "TROP2–immune depends partly on CLDN4: observational evidence matrix",
        loc="left", fontsize=12, pad=8,
    )
    fig.tight_layout()
    fig.savefig(FIG / "evidence_matrix.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIG / "evidence_matrix.pdf", bbox_inches="tight")
    plt.close()

    # Between-section sensitivity, for the write-up, not a second claim.
    b50 = between[between.readout == "cd8nk_50um"].iloc[0]
    b100 = between[between.readout == "cd8nk_100um"].iloc[0]
    sens = pd.DataFrame([
        {
            "readout": "50 µm",
            "rho_tacstd2": b50.rho_tacstd2,
            "p_tacstd2": b50.p_tacstd2,
            "rho_cldn4": b50.rho_cldn4,
            "p_cldn4": b50.p_cldn4,
            "partial_tac": b50.rho_tacstd2_given_cldn4,
            "partial_cldn": b50.rho_cldn4_given_tacstd2,
            "attenuation_pct": b50.attenuation_pct,
            "acme": b50.acme,
        },
        {
            "readout": "100 µm",
            "rho_tacstd2": b100.rho_tacstd2,
            "p_tacstd2": b100.p_tacstd2,
            "rho_cldn4": b100.rho_cldn4,
            "p_cldn4": b100.p_cldn4,
            "partial_tac": b100.rho_tacstd2_given_cldn4,
            "partial_cldn": b100.rho_cldn4_given_tacstd2,
            "attenuation_pct": b100.attenuation_pct,
            "acme": b100.acme,
        },
    ])
    sens.to_csv(TAB / "cosmx_between_section_display.tsv", sep="\t", index=False)
    from scipy import stats
    inv = pd.read_csv(TAB / "cosmx_section_inventory.tsv", sep="\t")
    loo = []
    for y, label in (("mean_cd8nk_50", "50 µm"), ("mean_cd8nk_100", "100 µm")):
        for i, row in inv.iterrows():
            sub = inv.drop(index=i)
            rc, pc = stats.spearmanr(sub.mean_cldn4_log1p_cp10k, sub[y])
            rt, pt = stats.spearmanr(sub.mean_tacstd2_log1p_cp10k, sub[y])
            loo.append({
                "readout": label,
                "left_out": row["sample"],
                "rho_cldn4": rc,
                "p_cldn4": pc,
                "rho_tacstd2": rt,
                "p_tacstd2": pt,
            })
    pd.DataFrame(loo).to_csv(TAB / "cosmx_between_section_loo.tsv", sep="\t", index=False)
    print(table_md)
    print("calls", matrix.verdict.value_counts().to_dict())


if __name__ == "__main__":
    sys.exit(main())
