#!/usr/bin/env python3
"""Write Part2 verdict matrix + three-layer rank figure from committed tables."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def main():
    c4 = json.loads((TAB / "c4_sweep_summary.json").read_text())
    tcga = json.loads((TAB / "tcga_sweep_summary.json").read_text())
    cosmx = json.loads((TAB / "cosmx_sweep_summary.json").read_text())

    c4_max_rho = c4["baseline_unadj_pct_full7"]
    c4_max_margin = c4["primary_max_margin"]
    c4_ranks = read_tsv(TAB / "c4_primary_rank.tsv")
    # also load locked-style unadj ranks from grid
    grid = read_tsv(TAB / "c4_sweep_grid.tsv")
    unadj = next(
        r
        for r in grid
        if r["filter"] == "all4" and r["panel"] == "full7" and r["score"] == "pct" and r["adjustment"] == "none"
    )

    tcga_anchor = tcga["anchor_pr748_like"]
    cosmx_anchor = cosmx["anchor_q4q1_50"]
    cosmx_primary = cosmx["primary"]

    verdict = [
        {
            "layer": "concordant4",
            "objective": "max_|rho|_among_CLDN4_lead",
            "spec": "all4|full7|pct|none",
            "cldn4_metric": "rho",
            "cldn4_value": c4_max_rho["cldn4_rho"],
            "lead_gene": "CLDN4",
            "runner_gene": c4_max_rho["runner_gene"],
            "margin": c4_max_rho["margin_runner_minus_cldn4"],
            "cldn4_clearly_wins": "point_estimate_yes;perm_vs_CLDN7_no",
            "part1_ok": False,
            "part2_ok": True,
            "note": "Locked-scale |ρ|=0.531; I2=0%; n=65. Separable from CLDN1/OCLN/F11R/CDH1, not from CLDN7/CLDN3.",
        },
        {
            "layer": "concordant4",
            "objective": "max_margin_among_CLDN4_lead",
            "spec": f"{c4_max_margin['filter']}|{c4_max_margin['panel']}|{c4_max_margin['score']}|{c4_max_margin['adjustment']}",
            "cldn4_metric": "rho",
            "cldn4_value": c4_max_margin["cldn4_rho"],
            "lead_gene": "CLDN4",
            "runner_gene": c4_max_margin["runner_gene"],
            "margin": c4_max_margin["margin_runner_minus_cldn4"],
            "cldn4_clearly_wins": "point_estimate_yes;perm_vs_CLDN3_no",
            "part1_ok": False,
            "part2_ok": True,
            "note": "Largest margin among eligible; |ρ| smaller than locked unadj. Searched maximum.",
        },
        {
            "layer": "tcga",
            "objective": "max_margin_CLDN4_lead",
            "spec": "none_eligible",
            "cldn4_metric": "rho",
            "cldn4_value": tcga_anchor["cldn4_rho"],
            "lead_gene": tcga_anchor["lead_gene"],
            "runner_gene": "",
            "margin": "",
            "cldn4_clearly_wins": False,
            "part1_ok": False,
            "part2_ok": True,
            "note": f"0/60 specs have CLDN4 as lead. Anchor CD8/KRT8-19/lung+funnel: lead={tcga_anchor['lead_gene']} {float(str(tcga_anchor['rank_table']).split(';')[0].split('=')[1]):.3f}; CLDN4 ρ={tcga_anchor['cldn4_rho']:.3f} (rank 5/6).",
        },
        {
            "layer": "cosmx",
            "objective": "anchor_q4q1_50",
            "spec": "q4_q1@50um",
            "cldn4_metric": "ratio",
            "cldn4_value": cosmx_anchor["cldn4_ratio"],
            "lead_gene": cosmx_anchor["lead_gene"],
            "runner_gene": "",
            "margin": "",
            "cldn4_clearly_wins": False,
            "part1_ok": False,
            "part2_ok": True,
            "note": "Only CLDN4+CDH1 on 960 panel; CLDN3/7/OCLN/F11R OFF. CDH1 colder. Locked 0.36/0.52 not replaced.",
        },
        {
            "layer": "cosmx",
            "objective": "max_margin_CLDN4_coldest_vs_CDH1",
            "spec": f"{cosmx_primary['mode']}@{cosmx_primary['radius_um']}um",
            "cldn4_metric": "ratio",
            "cldn4_value": cosmx_primary["cldn4_ratio"],
            "lead_gene": "CLDN4",
            "runner_gene": cosmx_primary["runner_gene"],
            "margin": cosmx_primary["margin_runner_minus_cldn4"],
            "cldn4_clearly_wins": "vs_CDH1_point_estimate_only;not_8of8;sister_claudins_OFF",
            "part1_ok": False,
            "part2_ok": True,
            "note": "5/16 specs CLDN4 colder than CDH1; none 8/8 & 5/5 for CLDN4 on those cuts.",
        },
    ]
    write_tsv(TAB / "ppt_verdict_matrix.tsv", verdict)

    # Three-panel figure
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.2))

    # C4: max |rho| unadj ranks from rank_table
    genes, rhos = [], []
    for part in unadj["rank_table"].split(";"):
        g, v = part.split("=")
        genes.append(g)
        rhos.append(float(v))
    ax = axes[0]
    colors = ["#b45309" if g == "CLDN4" else "#4b5563" for g in genes]
    ax.barh(range(len(genes))[::-1], rhos, color=colors)
    ax.set_yticks(range(len(genes))[::-1])
    ax.set_yticklabels(genes)
    ax.axvline(0, color="#9ca3af", lw=1)
    ax.set_xlabel("DL ρ vs T/NK")
    ax.set_title("Concordant-4\nmax |ρ| (CLDN4 leads)")

    # TCGA anchor ranks
    genes, rhos = [], []
    for part in tcga_anchor["rank_table"].split(";"):
        g, v = part.split("=")
        genes.append(g)
        rhos.append(float(v))
    ax = axes[1]
    colors = ["#b45309" if g == "CLDN4" else "#4b5563" for g in genes]
    ax.barh(range(len(genes))[::-1], rhos, color=colors)
    ax.set_yticks(range(len(genes))[::-1])
    ax.set_yticklabels(genes)
    ax.axvline(0, color="#9ca3af", lw=1)
    ax.set_xlabel("KRT-partial ρ vs CD8")
    ax.set_title("TCGA\nCLDN4 does not lead (0/60)")

    # CosMx: show both anchor and primary
    ax = axes[2]
    labels = [
        "anchor q4/q1@50\nCLDN4",
        "anchor q4/q1@50\nCDH1",
        "max-margin det@25\nCLDN4",
        "max-margin det@25\nCDH1",
    ]
    vals = [
        float(cosmx_anchor["cldn4_ratio"]),
        float(cosmx_anchor["rank_table"].split(";")[0].split("=")[1])
        if cosmx_anchor["lead_gene"] == "CDH1"
        else float(cosmx_anchor["rank_table"].split(";")[1].split("=")[1]),
        float(cosmx_primary["cldn4_ratio"]),
        float(cosmx_primary["runner_ratio"]),
    ]
    # fix anchor CDH1 value from rank_table
    for part in cosmx_anchor["rank_table"].split(";"):
        g, v = part.split("=")
        if g == "CDH1":
            vals[1] = float(v)
        if g == "CLDN4":
            vals[0] = float(v)
    colors = ["#b45309", "#4b5563", "#b45309", "#4b5563"]
    ax.barh(range(len(labels))[::-1], vals, color=colors)
    ax.set_yticks(range(len(labels))[::-1])
    ax.set_yticklabels(labels, fontsize=8)
    ax.axvline(1.0, color="#9ca3af", lw=1)
    ax.set_xlabel("CD8+NK neighbor ratio")
    ax.set_title("CosMx (CLDN3/7/OCLN/F11R OFF)\nCLDN4 wins only vs CDH1 at short r")

    fig.suptitle(
        "Part2 TJ gene screen — maximize |ρ| / margin so CLDN4 wins\n"
        "Do not force into Part1: TCGA fails; CosMx panel-limited; C4 not perm-separable from CLDN3/7",
        fontsize=10,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "ppt_three_layer_max_effect.png", dpi=150)
    fig.savefig(FIG / "ppt_three_layer_max_effect.pdf")
    plt.close(fig)
    print("wrote verdict + figure")


if __name__ == "__main__":
    main()
