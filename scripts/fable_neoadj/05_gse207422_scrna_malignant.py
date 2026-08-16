"""GSE207422 scRNA -- malignant/epithelial TACSTD2 vs MPR (leftover).

UMI matrix is genes x ~92k cells (BD_immune01..15). No published per-cell
type labels are on GEO, so epithelial/malignant cells are gated as
EPCAM>0 and PTPRC==0. Most samples are POST-treatment resections
(3 pre-treatment biopsies only, of which 1 is NE), so this is NOT a
baseline-biomarker test -- residual malignant TACSTD2 after therapy.
"""
from __future__ import annotations

import gzip
import json
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import GENES, DATA, FIG, TAB, mwu_stats

MATRIX = os.path.join(DATA, "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz")
META = os.path.join(DATA, "GSE207422_scRNA_metadata.xlsx")
TARGETS = GENES + ["EPCAM", "PTPRC", "KRT19", "KRT18"]


def load_meta():
    df = pd.read_excel(META, sheet_name="sheet1")
    df = df[df["Sample"].notna()].copy()
    df["Sample"] = df["Sample"].astype(str).str.strip()
    resp = df["Pathologic Response"].astype(str).str.strip()
    df["MPR"] = np.where(resp.str.upper().str.startswith("MPR") | (resp.str.lower() == "pcr"), 1,
                         np.where(resp.str.upper() == "NMPR", 0, np.nan))
    df["timing"] = np.where(df["Resource"].astype(str).str.contains("Pre", case=False),
                            "pre", "post")
    return df.set_index("Sample")


def extract_targets():
    with gzip.open(MATRIX, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cells = [c.strip('"') for c in header[1:]]
        want = set(TARGETS)
        rows = {}
        for line in fh:
            g = line.split("\t", 1)[0].strip('"')
            if g in want:
                rows[g] = np.array(line.rstrip("\n").split("\t")[1:], float)
                if len(rows) == len(want):
                    break
    return cells, rows


def main():
    meta = load_meta()
    cells, rows = extract_targets()
    sample = np.array(["_".join(c.split("_")[:2]) for c in cells])  # BD_immuneXX
    epcam = rows["EPCAM"]
    ptprc = rows["PTPRC"]
    malig = (epcam > 0) & (ptprc == 0)

    recs = []
    for sid in sorted(set(sample)):
        if sid not in meta.index:
            continue
        m = (sample == sid) & malig
        n_mal = int(m.sum())
        n_all = int((sample == sid).sum())
        tot = 0.0
        for g in TARGETS:
            tot += rows[g][m].sum() if n_mal else 0.0
        # use sum of extracted genes as a poor-man's size factor is biased;
        # instead report mean UMI of the gene among malignant cells and %pos
        row = {
            "sample": sid,
            "patient": meta.loc[sid, "Patient"],
            "timing": meta.loc[sid, "timing"],
            "response_raw": meta.loc[sid, "Pathologic Response"],
            "MPR": meta.loc[sid, "MPR"],
            "n_cells": n_all,
            "n_malignant": n_mal,
        }
        for g in GENES:
            v = rows[g][m]
            row[f"{g}_mean_umi"] = float(np.mean(v)) if n_mal else np.nan
            row[f"{g}_pct_pos"] = float(100.0 * np.mean(v > 0)) if n_mal else np.nan
            # library-size proxy: EPCAM+PTPRC- cells, CP10K using sum of marker genes is bad;
            # use mean log1p UMI (comparable within this matrix)
            row[f"{g}_mean_log1p"] = float(np.mean(np.log1p(v))) if n_mal else np.nan
        recs.append(row)
    df = pd.DataFrame(recs)
    df["pass_qc"] = df["n_malignant"] >= 20
    df.to_csv(os.path.join(TAB, "GSE207422_scrna_malignant_per_sample.csv"), index=False)
    print(df.to_string())

    # primary contrast: POST-treatment samples with MPR label + QC
    qc = df[(df["pass_qc"]) & (df["MPR"].notna()) & (df["timing"] == "post")].copy()
    results = []
    for g in GENES:
        col = f"{g}_mean_log1p"
        pos = qc.loc[qc["MPR"] == 1, col].values
        neg = qc.loc[qc["MPR"] == 0, col].values
        st = mwu_stats(pos, neg)
        row = {"dataset": "GSE207422_scRNA_malig_post", "gene": g,
               "comparison": "MPR_vs_nonMPR", "unit": "mean_log1p_UMI_EPCAM+_PTPRC-"}
        row.update(st)
        results.append(row)
    res = pd.DataFrame(results)
    res.to_csv(os.path.join(TAB, "GSE207422_scrna_malignant_stats.csv"), index=False)
    print(res.to_string())

    _boxplots(qc)

    summary = {
        "dataset": "GSE207422_scRNA_malig_post",
        "n_samples_total": int(len(df)),
        "n_post_qc": int(len(qc)),
        "n_MPR": int((qc["MPR"] == 1).sum()),
        "n_nonMPR": int((qc["MPR"] == 0).sum()),
        "gate": "EPCAM>0 & PTPRC==0 (no published malignant labels on GEO)",
        "timing": "mostly post-treatment residual tumor",
    }
    with open(os.path.join(TAB, "GSE207422_scrna_malignant_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))


def _boxplots(qc):
    if len(qc) == 0:
        return
    fig, axes = plt.subplots(1, len(GENES), figsize=(4.2 * len(GENES), 4.2))
    for ax, g in zip(axes, GENES):
        col = f"{g}_mean_log1p"
        groups = [qc.loc[qc["MPR"] == 0, col].values, qc.loc[qc["MPR"] == 1, col].values]
        bp = ax.boxplot(groups, widths=0.6, patch_artist=True, showfliers=False)
        for patch, c in zip(bp["boxes"], ["#8aa1b1", "#d1495b"]):
            patch.set_facecolor(c)
            patch.set_alpha(0.55)
        rng = np.random.default_rng(2)
        for i, vals in enumerate(groups, 1):
            if len(vals):
                ax.scatter(rng.normal(i, 0.06, len(vals)), vals, s=22, color="#22303c", zorder=3)
        st = mwu_stats(groups[1], groups[0]) if len(groups[0]) and len(groups[1]) else {"p": np.nan, "auc": np.nan}
        ax.set_xticks([1, 2])
        ax.set_xticklabels([f"non-MPR\n(n={len(groups[0])})", f"MPR/pCR\n(n={len(groups[1])})"])
        ax.set_ylabel("malignant mean log1p UMI")
        ax.set_title(f"{g}\nMWU p={st['p']:.3f}, AUC={st['auc']:.2f}")
    fig.suptitle("GSE207422 post-tx EPCAM+/PTPRC- TACSTD2/CLDN4 vs MPR", y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "GSE207422_scrna_malignant_boxplots.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
