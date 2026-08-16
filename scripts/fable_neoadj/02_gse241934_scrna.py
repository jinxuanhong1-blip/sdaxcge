"""GSE241934 -- scRNA validation (NEOTIDE/CTONG2104 + real-world cohort).

FACS-sorted living cells (7AAD-/CD235a-) from resected NSCLC tumors after
neoadjuvant anti-PD-1 +/- chemotherapy, annotated by cell type incl. an
epithelial/tumor compartment ("Epi"), with per-patient pathologic response
(MPR / pCR / non-MPR) (Zhang et al., PMID from series).

We build a per-patient EPITHELIAL pseudobulk of TACSTD2 / CLDN4 (CP10K,
log1p) and test it against MPR. Epithelial restriction isolates the
tumor-intrinsic ADC-target signal from immune/stromal contamination.
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

COHORTS = {
    "IIT": ("GSE241934_IIT_Matrix.mtx.gz", "GSE241934_IIT_barcodes.tsv.gz",
            "GSE241934_IIT_features.tsv.gz", "GSE241934_IIT_Meta.txt.gz"),
    "Real": ("GSE241934_Real_Matrix.mtx.gz", "GSE241934_RWC_barcodes.tsv.gz",
             "GSE241934_Real_features.tsv.gz", "GSE241934_Real_Meta.txt.gz"),
}
EPI_LABEL = "Epi"
TARGETS = GENES + ["EPCAM"]


def load_features(path):
    with gzip.open(os.path.join(DATA, path), "rt") as fh:
        return [l.split("\t")[0] for l in fh]


def load_meta(path):
    cols = {}
    with gzip.open(os.path.join(DATA, path), "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(hdr)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            cols[p[ix["cellID"]]] = (
                p[ix["major_cell_type"]], p[ix["sampleID"]],
                p[ix["Pathological Response"]], float(p[ix["nCount_RNA"]]))
    return cols


def process_cohort(name):
    mtx, bc, feat, meta = COHORTS[name]
    features = load_features(feat)
    gene_row = {}  # 1-based row index -> gene
    for i, g in enumerate(features, 1):
        if g in TARGETS:
            gene_row[i] = g
    metad = load_meta(meta)

    # column (1-based) -> cellID
    with gzip.open(os.path.join(DATA, bc), "rt") as fh:
        barcodes = [l.strip() for l in fh]
    ncol = len(barcodes)

    is_epi = np.zeros(ncol + 1, bool)
    sample_of = np.empty(ncol + 1, object)
    resp_of = np.empty(ncol + 1, object)
    for j, cid in enumerate(barcodes, 1):
        m = metad.get(cid)
        if m is None:
            continue
        mct, sid, resp, _ = m
        sample_of[j] = sid
        resp_of[j] = resp
        is_epi[j] = (mct == EPI_LABEL)

    # per-sample epithelial cell count + total counts (denominator)
    epi_cells = {}
    epi_total = {}
    for j in range(1, ncol + 1):
        if is_epi[j]:
            sid = sample_of[j]
            epi_cells[sid] = epi_cells.get(sid, 0) + 1
            epi_total[sid] = epi_total.get(sid, 0.0) + metad[barcodes[j - 1]][3]

    # gene count sums over epithelial cells, per sample; and detection counts
    gsum = {g: {} for g in TARGETS}
    gdet = {g: {} for g in TARGETS}
    with gzip.open(os.path.join(DATA, mtx), "rt") as fh:
        for line in fh:
            if line.startswith("%"):
                continue
            # first non-comment line is the size header
            dims = line.split()
            break
        for line in fh:
            r_s, c_s, v_s = line.split()
            r = int(r_s)
            g = gene_row.get(r)
            if g is None:
                continue
            c = int(c_s)
            if not is_epi[c]:
                continue
            v = float(v_s)
            sid = sample_of[c]
            gsum[g][sid] = gsum[g].get(sid, 0.0) + v
            gdet[g][sid] = gdet[g].get(sid, 0) + 1

    rows = []
    for sid in sorted(epi_cells):
        n_epi = epi_cells[sid]
        tot = epi_total[sid]
        # any epi cell of this sample carries the response label
        resp = None
        for j in range(1, ncol + 1):
            if is_epi[j] and sample_of[j] == sid:
                resp = resp_of[j]
                break
        row = {"cohort": name, "patient": sid, "n_epi_cells": n_epi, "response_raw": resp}
        for g in TARGETS:
            s = gsum[g].get(sid, 0.0)
            row[f"{g}_cp10k_log1p"] = float(np.log1p(1e4 * s / tot)) if tot > 0 else np.nan
            row[f"{g}_pct_pos"] = 100.0 * gdet[g].get(sid, 0) / n_epi if n_epi else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def responder(label):
    lab = str(label).strip().lower()
    if lab in ("mpr", "pcr"):
        return 1
    if lab in ("non-mpr", "nonmpr", "non_mpr"):
        return 0
    return np.nan


def main():
    parts = [process_cohort(c) for c in COHORTS]
    df = pd.concat(parts, ignore_index=True)
    df["MPR"] = df["response_raw"].map(responder)
    # keep patients with >=20 epithelial cells for a stable pseudobulk
    df["pass_qc"] = df["n_epi_cells"] >= 20
    df.to_csv(os.path.join(TAB, "GSE241934_scrna_per_patient.csv"), index=False)

    qc = df[df["pass_qc"] & df["MPR"].notna()].copy()
    results = []
    for g in GENES:
        col = f"{g}_cp10k_log1p"
        pos = qc.loc[qc["MPR"] == 1, col].values
        neg = qc.loc[qc["MPR"] == 0, col].values
        st = mwu_stats(pos, neg)
        row = {"dataset": "GSE241934_scRNA_epi", "gene": g,
               "comparison": "MPR_vs_nonMPR", "unit": "epi_pseudobulk_cp10k_log1p"}
        row.update(st)
        results.append(row)
    res = pd.DataFrame(results)
    res.to_csv(os.path.join(TAB, "GSE241934_scrna_stats.csv"), index=False)
    print(res.to_string())

    _boxplots(qc)

    summary = {
        "dataset": "GSE241934_scRNA_epi",
        "n_patients_total": int(len(df)),
        "n_patients_qc": int(len(qc)),
        "n_MPR": int((qc["MPR"] == 1).sum()),
        "n_nonMPR": int((qc["MPR"] == 0).sum()),
        "compartment": "epithelial (Epi) pseudobulk",
        "qc": ">=20 epithelial cells/patient",
    }
    with open(os.path.join(TAB, "GSE241934_scrna_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))


def _boxplots(qc):
    fig, axes = plt.subplots(1, len(GENES), figsize=(4.2 * len(GENES), 4.2))
    for ax, g in zip(axes, GENES):
        col = f"{g}_cp10k_log1p"
        groups = [qc.loc[qc["MPR"] == 0, col].values, qc.loc[qc["MPR"] == 1, col].values]
        bp = ax.boxplot(groups, widths=0.6, patch_artist=True, showfliers=False)
        for patch, c in zip(bp["boxes"], ["#8aa1b1", "#d1495b"]):
            patch.set_facecolor(c)
            patch.set_alpha(0.55)
        rng = np.random.default_rng(1)
        for i, vals in enumerate(groups, 1):
            ax.scatter(rng.normal(i, 0.06, len(vals)), vals, s=22, color="#22303c", alpha=0.8, zorder=3)
        st = mwu_stats(groups[1], groups[0])
        ax.set_xticks([1, 2])
        ax.set_xticklabels([f"non-MPR\n(n={len(groups[0])})", f"MPR/pCR\n(n={len(groups[1])})"])
        ax.set_ylabel("epithelial pseudobulk (log1p CP10K)")
        ax.set_title(f"{g}\nMWU p={st['p']:.3f}, AUC={st['auc']:.2f}")
    fig.suptitle("GSE241934 tumor-epithelial ADC targets vs MPR (n=45 patients)", y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "GSE241934_scrna_boxplots.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
