"""Within-patient PRE -> ON-treatment dynamics of TACSTD2/CLDN4 (GSE91061, Riaz).

This is the only cohort here with true same-patient paired pre / early-on-treatment
biopsies during PD-1 blockade (nivolumab). It is melanoma, not lung, so it is used
as an ORTHOGONAL sensitivity check on the *direction* of on-treatment change that
the lung neoadjuvant data (pre vs post, cross-patient) cannot resolve.

Response groups follow the deposited annotation: PRCR (PR/CR) = responder,
PD = non-responder, SD = stable (reported separately).
"""
from __future__ import annotations

import gzip

import numpy as np
import pandas as pd

import config as C
from common import compare_groups, paired_test, savefig, PALETTE
import matplotlib.pyplot as plt


def _load_meta() -> pd.DataFrame:
    titles, resp, visit = [], [], []
    with gzip.open(C.raw_path("GSE91061_matrix"), "rt") as fh:
        for line in fh:
            cells = [x.strip('"') for x in line.rstrip("\n").split("\t")]
            if line.startswith("!Sample_title"):
                titles = cells[1:]
            elif line.startswith("!Sample_characteristics_ch1"):
                if any("visit" in c for c in cells):
                    visit = [c.split(":", 1)[1].strip() for c in cells[1:]]
                elif any(c.startswith("response") for c in cells):
                    resp = [c.split(":", 1)[1].strip() for c in cells[1:]]
    meta = pd.DataFrame({"title": titles, "visit": visit, "response_raw": resp})
    meta["patient"] = meta["title"].str.extract(r"(Pt\d+)")
    meta["timepoint"] = meta["visit"].str.lower().map({"pre": "pre", "on": "on"})
    meta["response"] = meta["response_raw"].map(
        {"PRCR": "responder", "PD": "non_responder", "SD": "stable", "UNK": "unknown"})
    return meta.set_index("title")


def main() -> None:
    fpkm = pd.read_csv(C.raw_path("GSE91061_fpkm"), index_col=0)
    fpkm.index = fpkm.index.astype(str)
    meta = _load_meta()

    log = np.log2(fpkm + 1.0)
    recs = []
    for gene in C.TARGET_GENES:
        eid = C.ENTREZ_IDS[gene]
        vals = log.loc[eid]
        for title, v in vals.items():
            if title in meta.index:
                recs.append({"title": title, "gene": gene, "expr": float(v),
                             **meta.loc[title][["patient", "timepoint", "response"]].to_dict()})
    long = pd.DataFrame(recs).dropna(subset=["timepoint"])
    long.to_csv(C.TABLES_DIR / "gse91061_pre_on_long.csv", index=False)

    rows = []
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax, gene in zip(axes, C.TARGET_GENES):
        sub = long[long.gene == gene]
        wide = sub.pivot_table(index=["patient", "response"], columns="timepoint",
                               values="expr").dropna(subset=["pre", "on"]).reset_index()
        # paired test across all patients with both timepoints
        stat, p, med_delta = paired_test(wide["pre"].values, wide["on"].values)
        rows.append({"dataset": "GSE91061", "gene": gene,
                     "comparison": "paired on_vs_pre (all patients)",
                     "n_pairs": int(len(wide)),
                     "median_delta_on_minus_pre": med_delta,
                     "test": "Wilcoxon signed-rank", "statistic": stat, "p_value": p})

        # delta by response group
        wide["delta"] = wide["on"] - wide["pre"]
        for grp in ["responder", "non_responder"]:
            d = wide.loc[wide.response == grp, "delta"].values
            rows.append({"dataset": "GSE91061", "gene": gene,
                         "comparison": f"on-pre delta [{grp}]",
                         "n_pairs": int(len(d)),
                         "median_delta_on_minus_pre": float(np.median(d)) if len(d) else np.nan,
                         "test": "descriptive", "statistic": np.nan, "p_value": np.nan})
        # responder vs non-responder delta
        dr = wide.loc[wide.response == "responder", "delta"].values
        dnr = wide.loc[wide.response == "non_responder", "delta"].values
        cmp = compare_groups("GSE91061", gene, "on-pre delta: responder_vs_nonresponder",
                             dr, dnr, "responder", "non_responder", log_scale_input=True)
        rows.append(cmp.as_row())

        # plot paired lines coloured by response
        cmap = {"responder": PALETTE["responder"], "non_responder": PALETTE["non_responder"],
                "stable": "#999999", "unknown": "#cccccc"}
        for _, r in wide.iterrows():
            ax.plot([0, 1], [r["pre"], r["on"]], "-o", color=cmap.get(r["response"], "#ccc"),
                    alpha=0.7, markersize=4)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["Pre", "On"])
        ax.set_ylabel(f"{gene} (log2 FPKM+1)")
        ax.set_title(f"GSE91061 {gene}: paired pre -> on (p={p:.3f}, n={len(wide)})", fontsize=10)
        ax.spines[["top", "right"]].set_visible(False)
    handles = [plt.Line2D([0], [0], color=PALETTE["responder"], marker="o", label="responder (PR/CR)"),
               plt.Line2D([0], [0], color=PALETTE["non_responder"], marker="o", label="non-responder (PD)"),
               plt.Line2D([0], [0], color="#999999", marker="o", label="stable (SD)")]
    axes[1].legend(handles=handles, fontsize=8, frameon=False)
    fig.suptitle("Riaz melanoma (GSE91061) - within-patient pre->on nivolumab dynamics", fontsize=11)
    savefig(fig, C.FIGURES_DIR / "fig3_gse91061_paired_pre_on.png")

    out = pd.DataFrame(rows)
    out.to_csv(C.TABLES_DIR / "gse91061_paired_stats.csv", index=False)
    print(out.to_string())


if __name__ == "__main__":
    main()
