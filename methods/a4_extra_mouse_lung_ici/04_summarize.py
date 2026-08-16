#!/usr/bin/env python3
"""Build clean primary tables + a small figure from scored contrasts."""
from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

OUT = Path("results/a4_extra_mouse_lung_ici")
LOG_SPACE = {
    "GSE239485",
    "GSE297630",
    "GSE169194",
    "GSE137396",
    "GSE246922",
    "GSE338923",
    "GSE256071",
}


def main() -> None:
    contrasts = json.loads((OUT / "contrasts.json").read_text())
    samples = list(csv.DictReader((OUT / "per_sample.tsv").open(), delimiter="\t"))

    # recompute log2FC from per-sample for honesty
    by = defaultdict(list)
    for r in samples:
        if not r.get("value"):
            continue
        by[(r["accession"], r["gene"], r["group"])].append(float(r["value"]))

    clean = []
    for c in contrasts:
        if c["gene"] not in ("Tacstd2", "Cldn4"):
            continue
        a = by.get((c["accession"], c["gene"], c["group_a"]), [])
        b = by.get((c["accession"], c["gene"], c["group_b"]), [])
        rec = {
            "accession": c["accession"],
            "model": c["model"],
            "design": c["design"],
            "gene": c["gene"],
            "group_a": c["group_a"],
            "group_b": c["group_b"],
            "n_a": c.get("n_a"),
            "n_b": c.get("n_b"),
            "mean_a": c.get("mean_a"),
            "mean_b": c.get("mean_b"),
            "welch_p": c.get("welch_p"),
            "mwu_p": c.get("mwu_p"),
            "icb": c.get("icb_contrast"),
            "response": c.get("response_contrast"),
            "genotype": c.get("genotype_contrast"),
            "note": c.get("note"),
        }
        if a and b:
            rec["n_a"] = len(a)
            rec["n_b"] = len(b)
            rec["mean_a"] = float(np.mean(a))
            rec["mean_b"] = float(np.mean(b))
            if c["accession"] in LOG_SPACE:
                rec["log2FC_b_vs_a"] = float(np.mean(b) - np.mean(a))
                rec["value_space"] = "already_log_mean_diff"
            else:
                rec["log2FC_b_vs_a"] = float(np.log2((np.mean(b) + 1.0) / (np.mean(a) + 1.0)))
                rec["value_space"] = "log2_(mean+1)"
        else:
            # author DEG rows (GSE244452)
            rec["log2FC_b_vs_a"] = c.get("log2fc_mean")
            rec["value_space"] = "author_DEG_KP_vs_KL" if c["accession"] == "GSE244452" else c.get("space")
            if c["accession"] == "GSE244452" and c.get("log2fc_mean") is not None:
                rec["log2FC_KL_vs_KP"] = -float(c["log2fc_mean"])
                rec["author_p"] = c.get("welch_p")
        clean.append(rec)

    keys = [
        "accession",
        "model",
        "design",
        "gene",
        "group_a",
        "group_b",
        "n_a",
        "n_b",
        "mean_a",
        "mean_b",
        "log2FC_b_vs_a",
        "log2FC_KL_vs_KP",
        "welch_p",
        "mwu_p",
        "author_p",
        "value_space",
        "icb",
        "response",
        "genotype",
        "note",
    ]
    with (OUT / "primary_contrasts.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore", delimiter="\t")
        w.writeheader()
        for r in clean:
            w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in keys})

    # figure: Tacstd2 ICB contrasts with n>=2 both arms
    fig_rows = [
        r
        for r in clean
        if r["gene"] == "Tacstd2"
        and r.get("icb")
        and r.get("n_a", 0) >= 2
        and r.get("n_b", 0) >= 2
        and r.get("log2FC_b_vs_a") is not None
        and r["accession"] != "GSE244452"
    ]
    fig_rows.sort(key=lambda r: r["log2FC_b_vs_a"])
    if fig_rows:
        labels = [f"{r['accession']} {r['design'][:40]}" for r in fig_rows]
        xs = [r["log2FC_b_vs_a"] for r in fig_rows]
        cols = ["#b85c38" if (r.get("welch_p") or 1) < 0.05 else "#4a6fa5" for r in fig_rows]
        fig, ax = plt.subplots(figsize=(9, 5.2))
        ax.barh(range(len(xs)), xs, color=cols)
        ax.set_yticks(range(len(xs)))
        ax.set_yticklabels(labels, fontsize=7)
        ax.axvline(0, color="black", lw=0.8)
        ax.set_xlabel("Tacstd2 log2FC (ICB / later arm vs control)")
        ax.set_title("Extra mouse lung ICI RNA (not TISMO): Tacstd2")
        fig.tight_layout()
        fig.savefig(OUT / "fig_tacstd2_icb_log2fc.png", dpi=140)
        plt.close()

    (OUT / "summary_counts.json").write_text(
        json.dumps(
            {
                "n_contrasts_tacstd2_cldn4": len(clean),
                "n_icb_tacstd2_n_ge_2": len(fig_rows),
                "tismo_taken_as_given": "49/64 Tacstd2 up after ICB, p=5.8e-5",
                "response_labels": "none with per-mouse R vs NR",
            },
            indent=2,
        )
    )
    print("wrote", OUT / "primary_contrasts.tsv", "n=", len(clean))


if __name__ == "__main__":
    main()
