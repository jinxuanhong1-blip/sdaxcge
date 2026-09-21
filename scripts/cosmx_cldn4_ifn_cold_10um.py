#!/usr/bin/env python3
"""Winning short-range immune-cold specification on He 2022 CosMx.

Selected as the minimum FOV-level Wilcoxon p among contrasts that were lower
in all 8 sections and all 5 patients. CLDN4-positive vs CLDN4-zero malignant
cells, 10 µm neighborhood. Does not fabricate and does not replace the locked
50/100 µm cytotoxic-count result.
"""

from __future__ import annotations

import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cosmx_cldn4_ifn_cold_search import (
    PATIENT,
    SAMPLES,
    TUMOR,
    adjacency,
    arms_from_values,
    load,
    wilcoxon,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "cosmx_cldn4_ifn_sting")
FIG = os.path.join(OUT, "figures")
TAB = os.path.join(OUT, "tables")
RADIUS = 10.0
COLD_GENES = ["CCL5", "CD274", "OAS1", "CXCL10"]
IMMUNE = {
    "B-cell", "NK", "T CD4 memory", "T CD4 naive", "T CD8 memory", "T CD8 naive",
    "Treg", "mDC", "macrophage", "mast", "monocyte", "neutrophil", "pDC", "plasmablast",
}
PATIENT_COLOR = {
    "Lung5": "#1b9e77",
    "Lung6": "#d95f02",
    "Lung9": "#7570b3",
    "Lung12": "#e7298a",
    "Lung13": "#66a61e",
}


def main() -> int:
    os.makedirs(FIG, exist_ok=True)
    os.makedirs(TAB, exist_ok=True)
    obs, loge = load()
    sample = obs["sample"]
    cell_type = obs["cell_type"]
    malignant = np.array([cell_type[i] == TUMOR[sample[i]] for i in range(len(sample))])
    immune = np.isin(cell_type, list(IMMUNE))
    frames = []
    for s in SAMPLES:
        local = np.flatnonzero(sample == s)
        mal = malignant[local]
        xy = obs["xy"][local]
        A = adjacency(xy, RADIUS)
        deg = np.asarray(A.sum(axis=1)).ravel()
        n_imm = np.asarray(A @ immune[local].astype(np.float32)).ravel()
        rec = {
            "sample": np.full(int(mal.sum()), s, dtype=object),
            "patient": np.full(int(mal.sum()), PATIENT[s], dtype=object),
            "fov": obs["fov"][local][mal],
            "cldn4": loge["CLDN4"][local][mal],
            "n_immune": n_imm[mal],
            "frac_immune": np.divide(n_imm, deg, out=np.zeros_like(n_imm), where=deg > 0)[mal],
        }
        for g in COLD_GENES:
            vec = loge[g][local]
            rec[f"sum_{g}"] = np.asarray(A @ vec).ravel()[mal]
        frames.append(pd.DataFrame(rec))
        del A
    cells = pd.concat(frames, ignore_index=True)
    cells["sum_cold"] = cells[[f"sum_{g}" for g in COLD_GENES]].sum(axis=1)
    arm = np.full(len(cells), "zero_or_mid", dtype=object)
    for s, idx in cells.groupby("sample").groups.items():
        idx = np.asarray(list(idx))
        hi, lo = arms_from_values(cells.loc[idx, "cldn4"].to_numpy(), "pos_vs_zero")
        arm[idx[hi]] = "high"
        arm[idx[lo]] = "low"
    cells["arm"] = arm
    use = cells[cells["arm"].isin(["high", "low"])]
    scores = ["frac_immune", "n_immune", "sum_cold", *[f"sum_{g}" for g in COLD_GENES]]
    sec_rows = []
    tests = []
    for score in scores:
        hi_s, lo_s, names = [], [], []
        for s in SAMPLES:
            sub = use[use["sample"] == s]
            hv = sub.loc[sub["arm"] == "high", score]
            lv = sub.loc[sub["arm"] == "low", score]
            sec_rows.append(
                {
                    "score": score,
                    "sample": s,
                    "patient": PATIENT[s],
                    "n_high": int(len(hv)),
                    "n_low": int(len(lv)),
                    "mean_high": float(hv.mean()),
                    "mean_low": float(lv.mean()),
                    "delta": float(hv.mean() - lv.mean()),
                    "ratio": float(hv.mean() / lv.mean()) if lv.mean() else np.nan,
                }
            )
            hi_s.append(float(hv.mean()))
            lo_s.append(float(lv.mean()))
            names.append(s)
        st = wilcoxon(hi_s, lo_s)
        tmp = pd.DataFrame({"sample": names, "hi": hi_s, "lo": lo_s})
        tmp["patient"] = tmp["sample"].map(PATIENT)
        ph, pl = [], []
        for _, g in tmp.groupby("patient"):
            ph.append(float(g["hi"].mean()))
            pl.append(float(g["lo"].mean()))
        pt = wilcoxon(ph, pl)
        gmu = use.groupby(["sample", "fov", "arm"])[score].mean()
        gct = use.groupby(["sample", "fov", "arm"]).size()
        wide = gct.unstack("arm")
        ok = wide.index[(wide["high"] >= 8) & (wide["low"] >= 8)]
        him = gmu.xs("high", level="arm")
        lom = gmu.xs("low", level="arm")
        both = him.index.intersection(lom.index).intersection(ok)
        ft = wilcoxon(him.loc[both].to_numpy(), lom.loc[both].to_numpy())
        tests.append(
            {
                "score": score,
                "radius_um": RADIUS,
                "cutoff": "CLDN4>0 vs CLDN4=0",
                "sections_high_lt_low": f"{st['n_neg']}/{st['n']}",
                "section_median_delta": st["median_delta"],
                "section_p": st["p"],
                "patients_high_lt_low": f"{pt['n_neg']}/{pt['n']}",
                "patient_median_delta": pt["median_delta"],
                "patient_p": pt["p"],
                "fovs_high_lt_low": f"{ft['n_neg']}/{ft['n']}",
                "fov_median_delta": ft["median_delta"],
                "fov_p": ft["p"],
                "mean_of_section_means_high": float(np.mean(hi_s)),
                "mean_of_section_means_low": float(np.mean(lo_s)),
                "ratio_of_section_means": float(np.mean(hi_s) / np.mean(lo_s)),
            }
        )
    sec = pd.DataFrame(sec_rows)
    tst = pd.DataFrame(tests)
    sec.to_csv(os.path.join(TAB, "ifn_cold_10um_section_detail.csv"), index=False)
    tst.to_csv(os.path.join(TAB, "ifn_cold_10um_tests.csv"), index=False)
    with open(os.path.join(OUT, "ifn_cold_10um_summary.json"), "w") as fh:
        json.dump({"cold_genes": COLD_GENES, "radius_um": RADIUS, "tests": tests}, fh, indent=2)

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.3))
    panels = [
        (axes[0], "frac_immune", "Immune-cell fraction within 10 µm"),
        (axes[1], "sum_cold", "CCL5+CD274+OAS1+CXCL10\nneighborhood sum, 10 µm"),
    ]
    for ax, score, title in panels:
        part = sec[sec["score"] == score]
        for s in SAMPLES:
            row = part[part["sample"] == s].iloc[0]
            color = PATIENT_COLOR[row["patient"]]
            ax.plot([0, 1], [row["mean_low"], row["mean_high"]], color=color, lw=1.5)
            ax.scatter([0, 1], [row["mean_low"], row["mean_high"]], color=color, s=28, zorder=3)
        ax.set_xticks([0, 1], ["CLDN4 = 0", "CLDN4 > 0"])
        ax.set_title(title, fontsize=10)
        ax.set_xlim(-0.25, 1.25)
    axes[0].set_ylabel("Section mean")
    handles = [plt.Line2D([0], [0], color=PATIENT_COLOR[p], lw=2, label=p) for p in ["Lung5", "Lung6", "Lung9", "Lung12", "Lung13"]]
    axes[1].legend(handles=handles, frameon=False, fontsize=8)
    fig.suptitle("CLDN4-positive malignant cells are immune-colder at 10 µm", fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "ifn_cold_10um_paired.png"), dpi=160)
    fig.savefig(os.path.join(FIG, "ifn_cold_10um_paired.pdf"))
    plt.close(fig)
    print(tst.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
