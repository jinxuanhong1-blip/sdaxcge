#!/usr/bin/env python3
"""Analyze the anti-PD1 NSCLC Imaging Mass Cytometry dataset (Zenodo 8041882).

The RDS holds a SingleCellExperiment with a 41-marker IMC panel x 99,659 single
cells, annotated with PatientID, ImmuneStatus (Tcell_low/high), TherapyStatus
(pre/post anti-PD1), SampleID and a cluster label. We parse it with `rdata`
(no R required), then:
  1. Recover the marker panel and per-cell metadata.
  2. Report the panel and dataset shape.
  3. Compare mean marker expression pre vs post anti-PD1 therapy.
  4. Save tidy CSVs and a summary figure.

Note: TACSTD2 / CLDN4 are transcriptomic markers and are NOT part of this
targeted protein panel, so they cannot be quantified here; this dataset instead
characterises the immune/tumour micro-environment of ICI-treated NSCLC.
"""
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import rdata

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[2]
RDS = ROOT / "results" / "fable_zenodo" / "downloads" / "8041882" / "sce-clustering-allsamples-original_final.RDS"
OUT = ROOT / "results" / "fable_zenodo"


def to_str_list(arr):
    out = []
    for x in list(arr):
        if hasattr(x, "value"):
            x = x.value
        out.append(str(x))
    return out


def cat_to_series(cat):
    # rdata converts R factors to pandas Categorical
    try:
        return pd.Series(np.asarray(cat).astype(str))
    except Exception:
        return pd.Series([str(v) for v in list(cat)])


def main():
    conv = rdata.conversion.convert(rdata.parser.parse_file(str(RDS)))

    # Marker panel names (row order)
    markers = to_str_list(conv.rowRanges.partitioning.NAMES)
    n_markers = len(markers)

    # Assay matrices: shape (markers, cells)
    assays = conv.assays.data.listData
    exprs = np.asarray(assays["exprs"])      # arcsinh-transformed counts
    counts = np.asarray(assays["counts"])
    logc = np.asarray(assays["logcounts"])
    n_cells = exprs.shape[1]
    assert exprs.shape[0] == n_markers, (exprs.shape, n_markers)

    # Cell metadata
    cd = conv.colData.listData
    meta = pd.DataFrame({
        "PatientID": cat_to_series(cd["PatientID"]),
        "ImmuneStatus": cat_to_series(cd["ImmuneStatus"]),
        "TherapyStatus": cat_to_series(cd["TherapyStatus"]),
        "SampleID": pd.Series(to_str_list(cd["SampleID"])),
        "Cluster": cat_to_series(cd["Cluster_15.sub"]),
    })
    print(f"IMC dataset: {n_markers} markers x {n_cells} cells")
    print("Markers:", ", ".join(markers))
    print("\nTherapyStatus counts:\n", meta["TherapyStatus"].value_counts().to_string())
    print("\nImmuneStatus counts:\n", meta["ImmuneStatus"].value_counts().to_string())
    print("\nPatients:", sorted(meta["PatientID"].unique()))

    # Save panel + metadata summary
    pd.DataFrame({"marker": markers}).to_csv(OUT / "ici_imc_panel.csv", index=False)
    meta_summary = {
        "n_markers": n_markers,
        "n_cells": n_cells,
        "therapy_status": meta["TherapyStatus"].value_counts().to_dict(),
        "immune_status": meta["ImmuneStatus"].value_counts().to_dict(),
        "n_patients": int(meta["PatientID"].nunique()),
        "n_samples": int(meta["SampleID"].nunique()),
    }
    import json
    (OUT / "ici_imc_summary.json").write_text(json.dumps(meta_summary, indent=2, default=str))

    # Pre vs post anti-PD1: mean arcsinh expression per marker
    expr_df = pd.DataFrame(exprs.T, columns=markers)
    expr_df["TherapyStatus"] = meta["TherapyStatus"].values
    grp = expr_df.groupby("TherapyStatus")[markers].mean().T
    # Order by pre/post difference if both present
    have = [c for c in ("pre", "post") if c in grp.columns]
    if len(have) == 2:
        grp["delta_post_minus_pre"] = grp["post"] - grp["pre"]
        grp = grp.sort_values("delta_post_minus_pre", ascending=False)
    grp.to_csv(OUT / "ici_imc_marker_pre_post.csv")
    print("\n=== Mean arcsinh marker expression: pre vs post anti-PD1 ===")
    print(grp.round(3).to_string())

    # Figure: pre vs post grouped bars
    if len(have) == 2:
        order = grp.index.tolist()
        y = np.arange(len(order))
        fig, ax = plt.subplots(figsize=(7, max(4, 0.32 * len(order))))
        ax.barh(y - 0.2, grp.loc[order, "pre"], height=0.4, label="pre anti-PD1", color="#4C72B0")
        ax.barh(y + 0.2, grp.loc[order, "post"], height=0.4, label="post anti-PD1", color="#C44E52")
        ax.set_yticks(y)
        ax.set_yticklabels(order, fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel("mean arcsinh expression")
        ax.set_title("Anti-PD1 NSCLC IMC (Zenodo 8041882)\nmarker expression pre vs post therapy")
        ax.legend()
        fig.tight_layout()
        fig.savefig(OUT / "ici_imc_pre_post_markers.png", dpi=150)
        print(f"\nWrote figure {OUT / 'ici_imc_pre_post_markers.png'}")
    print(f"Wrote {OUT / 'ici_imc_marker_pre_post.csv'}, ici_imc_panel.csv, ici_imc_summary.json")


if __name__ == "__main__":
    main()
