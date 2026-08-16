#!/usr/bin/env python3
"""Analyse TACSTD2 / CLDN4 detectability and ICI-response association in the
GSE285888 baseline-PBMC scRNA-seq dataset (33 NSCLC patients, 222,144 cells).

Design
------
* Expression values are raw UMI counts. We report a detection rate
  (fraction of cells with count > 0), which is normalisation-independent, and a
  CP10K-normalised mean (count / nCount_RNA * 1e4) for magnitude.
* Response labels (SubType): CR (complete response), DR (durable response),
  PD (progressive disease); five patients are labelled only by irAE and are
  excluded from the efficacy contrast. irAE severity: non / mild / severe.
* Statistics use patient-level pseudobulk values (one number per patient) so the
  unit of analysis is the patient, not the cell. Groups are compared with the
  Mann-Whitney U test (two-sided). Effector genes PRF1/GZMB act as a positive
  control: the source study reports them as response-associated, so a working
  pipeline should recover signal there even if the epithelial targets cannot.

Outputs (results/fable_blood_ici/):
    gse285888_gene_detectability.csv
    gse285888_detection_by_response.csv
    gse285888_patient_pseudobulk.csv
    gse285888_response_stats.csv
    gse285888_irae_stats.csv
    gse285888_detectability.png
    gse285888_summary.json
"""
import os
import json
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

GEO_DIR = os.environ.get("GEO_DIR", "/tmp/geo")
RESULTS = os.path.join(os.path.dirname(__file__), "..", "..",
                       "results", "fable_blood_ici")
RESULTS = os.path.abspath(RESULTS)
IN = os.path.join(GEO_DIR, "gse285888_panel_cells.parquet")

TARGETS = ["TACSTD2", "CLDN4"]
EPITHELIAL = ["EPCAM", "KRT8", "KRT18"]
IMMUNE = ["PTPRC", "CD3D", "CD8A", "MS4A1", "CD14", "NKG7"]
EFFECTOR = ["GZMB", "PRF1", "IL1B", "CXCL8"]
PANEL = TARGETS + EPITHELIAL + IMMUNE + EFFECTOR


def cp10k(df, genes):
    """CP10K normalise raw counts by per-cell total UMI (nCount_RNA)."""
    factor = 1e4 / df["nCount_RNA"].to_numpy()[:, None]
    return pd.DataFrame(df[genes].to_numpy() * factor, index=df.index,
                        columns=genes)


def main():
    os.makedirs(RESULTS, exist_ok=True)
    df = pd.read_parquet(IN)
    n_cells = len(df)
    norm = cp10k(df, PANEL)

    # ---- 1. Global detectability -------------------------------------------
    rows = []
    for g in PANEL:
        detected = (df[g] > 0)
        rows.append({
            "gene": g,
            "class": ("target" if g in TARGETS else
                      "epithelial" if g in EPITHELIAL else
                      "immune" if g in IMMUNE else "effector"),
            "n_cells": n_cells,
            "n_detected": int(detected.sum()),
            "pct_cells_detected": round(100 * detected.mean(), 4),
            "mean_raw_count": round(float(df[g].mean()), 5),
            "mean_cp10k": round(float(norm[g].mean()), 4),
        })
    det = pd.DataFrame(rows).sort_values("pct_cells_detected", ascending=False)
    det.to_csv(os.path.join(RESULTS, "gse285888_gene_detectability.csv"),
               index=False)

    # ---- 2. Detection rate by response group -------------------------------
    by_resp = []
    for g in TARGETS + ["EPCAM", "PRF1", "GZMB"]:
        for grp, sub in df.groupby("SubType"):
            by_resp.append({
                "gene": g, "SubType": grp, "n_cells": len(sub),
                "pct_cells_detected": round(100 * (sub[g] > 0).mean(), 4),
                "mean_cp10k": round(float(cp10k(sub, [g])[g].mean()), 4),
            })
    by_resp = pd.DataFrame(by_resp)
    by_resp.to_csv(os.path.join(RESULTS, "gse285888_detection_by_response.csv"),
                   index=False)

    # ---- 3. Patient-level pseudobulk ---------------------------------------
    pb_rows = []
    for pid, sub in df.groupby("orig.ident"):
        subn = cp10k(sub, PANEL)
        rec = {"patient": pid, "SubType": sub["SubType"].iloc[0],
               "irAE": sub["irAE"].iloc[0], "n_cells": len(sub)}
        for g in PANEL:
            rec[f"{g}_detrate"] = float((sub[g] > 0).mean())
            rec[f"{g}_cp10k"] = float(subn[g].mean())
        pb_rows.append(rec)
    pb = pd.DataFrame(pb_rows)
    pb.to_csv(os.path.join(RESULTS, "gse285888_patient_pseudobulk.csv"),
              index=False)

    # ---- 4. Response association (Mann-Whitney on patient pseudobulk) ------
    def mw(a, b):
        a, b = np.asarray(a, float), np.asarray(b, float)
        if len(a) < 2 or len(b) < 2 or (np.all(a == a[0]) and np.all(b == b[0])
                                        and a[0] == b[0]):
            return np.nan, np.nan
        u, p = mannwhitneyu(a, b, alternative="two-sided")
        return float(u), float(p)

    contrasts = {
        "CR_vs_PD": (pb[pb.SubType == "CR"], pb[pb.SubType == "PD"]),
        "responder(CR+DR)_vs_PD": (pb[pb.SubType.isin(["CR", "DR"])],
                                   pb[pb.SubType == "PD"]),
    }
    resp_rows = []
    for g in PANEL:
        for cname, (ga, gb) in contrasts.items():
            for metric in ["detrate", "cp10k"]:
                col = f"{g}_{metric}"
                u, p = mw(ga[col], gb[col])
                resp_rows.append({
                    "gene": g, "contrast": cname, "metric": metric,
                    "group_a_mean": round(float(ga[col].mean()), 6),
                    "group_b_mean": round(float(gb[col].mean()), 6),
                    "n_a": len(ga), "n_b": len(gb),
                    "U": u, "p_value": p,
                })
    resp = pd.DataFrame(resp_rows)
    resp.to_csv(os.path.join(RESULTS, "gse285888_response_stats.csv"),
                index=False)

    # ---- 5. irAE association (severe vs non) -------------------------------
    irae_rows = []
    ga = pb[pb.irAE == "severe_irAE"]
    gb = pb[pb.irAE == "non"]
    for g in PANEL:
        for metric in ["detrate", "cp10k"]:
            col = f"{g}_{metric}"
            u, p = mw(ga[col], gb[col])
            irae_rows.append({
                "gene": g, "contrast": "severe_irAE_vs_non", "metric": metric,
                "severe_mean": round(float(ga[col].mean()), 6),
                "non_mean": round(float(gb[col].mean()), 6),
                "n_severe": len(ga), "n_non": len(gb), "U": u, "p_value": p,
            })
    irae = pd.DataFrame(irae_rows)
    irae.to_csv(os.path.join(RESULTS, "gse285888_irae_stats.csv"), index=False)

    # ---- 6. Figure: detection rate per gene --------------------------------
    fig, ax = plt.subplots(figsize=(9, 5))
    colors = {"target": "#d62728", "epithelial": "#ff7f0e",
              "immune": "#1f77b4", "effector": "#2ca02c"}
    d = det.sort_values("pct_cells_detected")
    ax.barh(d["gene"], d["pct_cells_detected"],
            color=[colors[c] for c in d["class"]])
    ax.set_xlabel("% of PBMC cells with detected expression (count > 0)")
    ax.set_title("GSE285888 baseline PBMC (n=222,144 cells, 33 NSCLC patients)\n"
                 "TACSTD2 / CLDN4 detectability vs reference markers")
    handles = [plt.Rectangle((0, 0), 1, 1, color=colors[c]) for c in colors]
    ax.legend(handles, list(colors.keys()), title="gene class",
              loc="lower right")
    for i, (g, v) in enumerate(zip(d["gene"], d["pct_cells_detected"])):
        ax.text(v + 0.5, i, f"{v:.2f}%", va="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS, "gse285888_detectability.png"), dpi=140)

    # ---- 7. Machine-readable summary ---------------------------------------
    summary = {
        "dataset": "GSE285888",
        "n_cells": int(n_cells),
        "n_patients": int(df["orig.ident"].nunique()),
        "response_groups": df.groupby("orig.ident")["SubType"].first()
                             .value_counts().to_dict(),
        "targets": {},
    }
    for g in TARGETS:
        row = det[det.gene == g].iloc[0]
        summary["targets"][g] = {
            "pct_cells_detected": float(row["pct_cells_detected"]),
            "mean_cp10k": float(row["mean_cp10k"]),
        }
    with open(os.path.join(RESULTS, "gse285888_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    print(det.to_string(index=False))
    print("\nResponse stats (targets + effector controls):")
    print(resp[resp.gene.isin(TARGETS + ["PRF1", "GZMB"])].to_string(index=False))


if __name__ == "__main__":
    main()
