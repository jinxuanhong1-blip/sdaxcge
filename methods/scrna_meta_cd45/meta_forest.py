#!/usr/bin/env python3
"""Pool patient-level TACSTD2/CLDN4 detection fraction vs MPR. Forest + extra table."""

from __future__ import annotations

import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from config import RESULTS
from stats_lib import iv_meta, two_group


def _study_effect(patient: pd.DataFrame, dataset: str, gene: str, note: str) -> dict:
    a = patient.loc[patient["mpr_any"] == "yes", f"{gene}_frac_pos"].to_numpy()
    b = patient.loc[patient["mpr_any"] == "no", f"{gene}_frac_pos"].to_numpy()
    t = two_group(a, b, "MPR", "non-MPR")
    t.update({
        "dataset": dataset,
        "gene": gene,
        "metric": f"{gene}_frac_pos",
        "note": note,
        "n_cells_mpr": int(patient.loc[patient["mpr_any"] == "yes", "n_cells"].sum()),
        "n_cells_nonmpr": int(patient.loc[patient["mpr_any"] == "no", "n_cells"].sum()),
    })
    return t


def _forest(effects: list[dict], meta: dict, title: str, out_png: str, out_pdf: str) -> None:
    rows = list(effects)
    labels = []
    ys = []
    los = []
    his = []
    annot = []
    for e in rows:
        labels.append(f"{e['dataset']}  {e['n_a']} vs {e['n_b']}")
        ys.append(e["g_hedges_g"])
        los.append(e["g_ci_lo"])
        his.append(e["g_ci_hi"])
        annot.append(f"g={e['g_hedges_g']:+.2f}  p={e['mw_p']:.3g}")
    if meta.get("k", 0) >= 1:
        labels.append(f"RE meta  {meta['n_total_mpr']} vs {meta['n_total_nonmpr']}")
        ys.append(meta["random_g"])
        los.append(meta["random_ci_lo"])
        his.append(meta["random_ci_hi"])
        annot.append(f"g={meta['random_g']:+.2f}  p={meta['random_p']:.3g}")

    n = len(labels)
    fig, ax = plt.subplots(figsize=(8.2, 1.3 + 0.55 * n))
    y_pos = np.arange(n)[::-1]
    for i, (y, lo, hi, lab_i) in enumerate(zip(ys, los, his, y_pos)):
        is_meta = i == n - 1 and meta.get("k", 0) >= 1
        ax.plot([lo, hi], [lab_i, lab_i], color="#222" if is_meta else "#4C72B0", lw=2)
        ax.plot(y, lab_i, "D" if is_meta else "o", color="#222" if is_meta else "#4C72B0", ms=7)
        ax.text(0.99, lab_i, annot[i], transform=ax.get_yaxis_transform(),
                va="center", ha="right", fontsize=8, color="#333")
    ax.axvline(0, color="#888", lw=1, ls="--")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Hedges' g  (MPR − non-MPR)\nnegative = lower detection fraction in MPR")
    ax.set_title(title, loc="left", fontsize=11)
    ax.set_xlim(min(los + [-1.2]) - 0.15, max(his + [0.4]) + 0.15)
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    fig.savefig(out_pdf)
    plt.close(fig)


def main() -> None:
    g243 = pd.read_csv(RESULTS / "gse243013_patient_level.tsv", sep="\t")
    g229 = pd.read_csv(RESULTS / "gse229353_patient_level.tsv", sep="\t")
    g243 = g243[g243["mpr_any"].isin(["yes", "no"])].copy()
    if "mtx_complete" in g229.columns:
        g229 = g229[g229["mtx_complete"].astype(str).isin(["True", "true", "1"])].copy()

    table_rows = []
    forest_by_gene = {}
    for gene in ("TACSTD2", "CLDN4"):
        e243 = _study_effect(
            g243, "GSE243013", gene,
            "CD45+ deposited immune atlas; GEO pCR+MPR vs non-MPR",
        )
        e229 = _study_effect(
            g229, "GSE229353", gene,
            "CD45+ bead-sorted; paper Table S1; P03 MTX truncated on GEO and dropped",
        )
        meta = iv_meta([
            {"hedges_g": e243["g_hedges_g"], "se": e243["g_se"], "n1": e243["n_a"], "n2": e243["n_b"]},
            {"hedges_g": e229["g_hedges_g"], "se": e229["g_se"], "n1": e229["n_a"], "n2": e229["n_b"]},
        ])
        forest_by_gene[gene] = {"studies": [e243, e229], "meta": meta}
        for e in (e243, e229):
            table_rows.append({
                "gene": gene,
                "dataset": e["dataset"],
                "compartment": "immune library (not malignant RNA)",
                "n_mpr": e["n_a"],
                "n_nonmpr": e["n_b"],
                "median_frac_mpr": e["g_median_a"],
                "median_frac_nonmpr": e["g_median_b"],
                "mean_frac_mpr": e["g_mean_a"],
                "mean_frac_nonmpr": e["g_mean_b"],
                "mannwhitney_p": e["mw_p"],
                "cliffs_delta_mpr_minus_nonmpr": e["mw_cliffs_delta"],
                "hedges_g": e["g_hedges_g"],
                "hedges_g_ci_lo": e["g_ci_lo"],
                "hedges_g_ci_hi": e["g_ci_hi"],
                "note": e["note"],
            })
        table_rows.append({
            "gene": gene,
            "dataset": "RE_meta",
            "compartment": "immune library (not malignant RNA)",
            "n_mpr": meta["n_total_mpr"],
            "n_nonmpr": meta["n_total_nonmpr"],
            "median_frac_mpr": np.nan,
            "median_frac_nonmpr": np.nan,
            "mean_frac_mpr": np.nan,
            "mean_frac_nonmpr": np.nan,
            "mannwhitney_p": meta["random_p"],
            "cliffs_delta_mpr_minus_nonmpr": np.nan,
            "hedges_g": meta["random_g"],
            "hedges_g_ci_lo": meta["random_ci_lo"],
            "hedges_g_ci_hi": meta["random_ci_hi"],
            "note": f"DerSimonian-Laird RE; I2={meta['I2']:.2f}; Q={meta['Q']:.2f}; k={meta['k']}",
        })
        _forest(
            [e243, e229],
            meta,
            f"{gene} detection fraction in CD45+ / immune libraries vs MPR",
            str(RESULTS / f"forest_{gene.lower()}_frac_pos.png"),
            str(RESULTS / f"forest_{gene.lower()}_frac_pos.pdf"),
        )

    extra = pd.DataFrame(table_rows)
    extra.to_csv(RESULTS / "extra_table_immune_leak_vs_mpr.tsv", sep="\t", index=False)
    (RESULTS / "meta_summary.json").write_text(json.dumps(forest_by_gene, indent=2, default=float))

    # Combined patient file for the extra
    keep = [
        "dataset", "patient_id", "mpr_any", "n_cells",
        "TACSTD2_frac_pos", "CLDN4_frac_pos", "compartment", "not_malignant_rna",
    ]
    a = g243.rename(columns={"sampleID": "patient_id"})
    b = g229.rename(columns={"patient": "patient_id"})
    cols = [c for c in keep if c in a.columns]
    pd.concat([a[cols], b[[c for c in keep if c in b.columns]]], ignore_index=True).to_csv(
        RESULTS / "all_patients_immune_leak.tsv", sep="\t", index=False
    )
    print(extra.to_string(index=False))


if __name__ == "__main__":
    main()
