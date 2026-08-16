"""GSE205335 leftover -- malignant TACSTD2 vs ICI response (RECIST, not MPR).

HONEST MISMATCH: this is advanced/metastatic lung cancer treated with ICI,
endpoint is RECIST (PR / SD / PD), NOT neoadjuvant MPR/pCR. Included because
it is the only public lung-ICI scRNA series with published malignant-cell
labels and TACSTD2 on GEO under 2GB (Kim/Ahn/Lee, GSE205335).

Malignant cells: lineage.sub == 'Malignant cells' (28,512 cells).
Per-patient malignant pseudobulk (log1p CP10K) is tested PR vs SD+PD.
NE and normal-tissue-only samples are dropped.
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

IDENT = os.path.join(DATA, "GSE205335_Lung_IO_CellIdentity.txt.gz")
COUNTS = os.path.join(TAB, "GSE205335_target_counts.csv")

# GEO GSM characteristics (patient, tissue, stage, subtype, RECIST)
# orig.ident in the matrix is SITE-NN-3P/5P matching the title suffix.
SAMPLE_MAP = [
    # orig.ident, patient, tissue, stage, subtype, recist
    ("EBUS-06-3P", "P1006", "Lung", "IV", "ADC", "PR"),
    ("EFFUSION-06-3P", "P1006", "Effusion", "IV", "ADC", "PR"),
    ("EBUS-27-5P", "P1027", "LN", "IV", "ADC", "PR"),
    ("EBUS-37-3P", "P1037", "LN", "IVA", "SQ", "PR"),
    ("EBUS-90-5P", "P1090", "LN", "IIIA", "SQ", "PR"),
    ("EBUS-17-5P", "P1017", "LN", "IV", "SQ", "SD"),
    ("LM-17-5P", "P1017", "Liver", "IV", "SQ", "SD"),
    ("PCNB-01-5P", "P4001", "Lung", "IV", "ADC", "SD"),
    ("EBUS-30-5P", "P1030", "LN", "IV", "ADC", "PD"),
    ("EBUS-62-5P", "P1062", "LN", "IV", "ADC", "PD"),
    ("EBUS-76-3P", "P1076", "LN", "IV", "ADC", "PD"),
    ("NECK-05-3P", "P1076", "LN", "IV", "ADC", "PD"),
    ("EBUS-89-3P", "P1089", "LN", "IV", "ADC", "PD"),
    ("EBUS-119-3P", "P1119", "LN", "IV", "ADC", "PD"),
    ("EBUS-16-5P", "P1016", "LN", "ED", "SCLC", "PR"),
    ("LM-16-5P", "P1016", "Liver", "ED", "SCLC", "PR"),
    ("LM-115-3P", "P1115", "Liver", "ED", "SCLC", "PR"),
    ("EBUS-56-3P", "P1056", "LN", "IV", "NUT", "SD"),
    ("LM-56-3P", "P1056", "Liver", "IV", "NUT", "SD"),
    ("EBUS-15-3P", "P1015", "LN", "IIIA", "ADC", "PD"),
    ("EBUS-25-3P", "P1025", "LN", "IV", "SCLC", "PD"),
    ("EBUS-84-5P", "P1084", "LN", "IV", "ADC", "NE"),
    ("EBUS-18-5P", "P1018", "Bronchus", "IV", "ADC", "NE"),
    ("BRONCHO-63-5P", "P1063", "Bronchus", "IV", "ADC", "NE"),
    ("EBUS-63-5P", "P1063", "LN", "IV", "ADC", "NE"),
    ("EBUS-79-5P", "P1079", "LN", "IV", "ADC", "NE"),
    ("EBUS-72-5P", "P1072", "LN", "ED", "SCLC", "NE"),
    ("LUNG-T31-3P", "P0031", "Lung", "IIIA", "ADC", "NE"),
    ("LUNG-N31-3P", "P0031", "Normal Lung", "IIIA", "ADC", "NE"),
    ("LN-16-5P", "P2016", "Normal LN", "IIB", "SQ", "PD"),
    ("LN-01-3P", "P2001", "Normal LN", "IIB", "ADC", "PR"),
    ("LN-09-5P", "P2009", "Normal LN", "IB", "ADC", "PR"),
    ("NS-32-5P", "P3032", "Normal Brain", "IV", "ADC", "PR"),
]
NORMAL_TISSUES = {"Normal Lung", "Normal LN", "Normal Brain"}


def load_identity():
    return pd.read_csv(IDENT, sep="\t")


def main():
    ident = load_identity()
    counts = pd.read_csv(COUNTS)
    smap = pd.DataFrame(SAMPLE_MAP, columns=["orig.ident", "patient", "tissue",
                                             "stage", "subtype", "recist"])
    df = ident.merge(counts, on="barcode", how="inner")
    df = df.merge(smap, on="orig.ident", how="left")
    print("barcode overlap", len(df), "of identity", len(ident), "counts", len(counts))
    print("lineage.sub:", df["lineage.sub"].value_counts().to_dict())

    mal = df[(df["lineage.sub"] == "Malignant cells") & (~df["tissue"].isin(NORMAL_TISSUES))].copy()
    # per-patient malignant pseudobulk
    rows = []
    for pid, sub in mal.groupby("patient"):
        tot = sub["nCount"].sum()
        recist = sub["recist"].iloc[0]
        subtype = sub["subtype"].iloc[0]
        tissues = ",".join(sorted(sub["tissue"].unique()))
        row = {"patient": pid, "recist": recist, "subtype": subtype, "tissues": tissues,
               "n_malignant": int(len(sub))}
        for g in GENES:
            s = sub[g].sum()
            row[f"{g}_cp10k_log1p"] = float(np.log1p(1e4 * s / tot)) if tot > 0 else np.nan
            row[f"{g}_pct_pos"] = float(100.0 * (sub[g] > 0).mean())
        rows.append(row)
    per = pd.DataFrame(rows)
    per["pass_qc"] = per["n_malignant"] >= 20
    per["responder"] = per["recist"].map({"PR": 1, "SD": 0, "PD": 0})
    per.to_csv(os.path.join(TAB, "GSE205335_malignant_per_patient.csv"), index=False)
    print(per.to_string())

    qc = per[per["pass_qc"] & per["responder"].notna()].copy()
    # also NSCLC-only (drop SCLC/NUT) as a sensitivity
    nsclc = qc[qc["subtype"].isin(["ADC", "SQ"])].copy()

    results = []
    for label, sub in [("all_with_RECIST", qc), ("NSCLC_ADC_SQ", nsclc)]:
        for g in GENES:
            col = f"{g}_cp10k_log1p"
            pos = sub.loc[sub["responder"] == 1, col].values
            neg = sub.loc[sub["responder"] == 0, col].values
            st = mwu_stats(pos, neg)
            row = {"dataset": f"GSE205335_malig_{label}", "gene": g,
                   "comparison": "RECIST_PR_vs_SDPD", "unit": "malig_pseudobulk_cp10k_log1p"}
            row.update(st)
            results.append(row)
    res = pd.DataFrame(results)
    res.to_csv(os.path.join(TAB, "GSE205335_malignant_stats.csv"), index=False)
    print(res.to_string())

    _boxplots(qc, nsclc)

    summary = {
        "dataset": "GSE205335_malignant",
        "mismatch": "advanced/metastatic ICI, RECIST not neoadjuvant MPR",
        "n_patients_qc": int(len(qc)),
        "n_PR": int((qc["responder"] == 1).sum()),
        "n_SDPD": int((qc["responder"] == 0).sum()),
        "n_NSCLC_qc": int(len(nsclc)),
        "compartment": "lineage.sub == Malignant cells",
    }
    with open(os.path.join(TAB, "GSE205335_malignant_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))


def _boxplots(qc, nsclc):
    fig, axes = plt.subplots(1, len(GENES), figsize=(4.2 * len(GENES), 4.2))
    for ax, g in zip(axes, GENES):
        col = f"{g}_cp10k_log1p"
        groups = [qc.loc[qc["responder"] == 0, col].values,
                  qc.loc[qc["responder"] == 1, col].values]
        bp = ax.boxplot(groups, widths=0.6, patch_artist=True, showfliers=False)
        for patch, c in zip(bp["boxes"], ["#8aa1b1", "#d1495b"]):
            patch.set_facecolor(c)
            patch.set_alpha(0.55)
        rng = np.random.default_rng(3)
        for i, vals in enumerate(groups, 1):
            if len(vals):
                ax.scatter(rng.normal(i, 0.06, len(vals)), vals, s=22, color="#22303c", zorder=3)
        st = mwu_stats(groups[1], groups[0])
        ax.set_xticks([1, 2])
        ax.set_xticklabels([f"SD/PD\n(n={len(groups[0])})", f"PR\n(n={len(groups[1])})"])
        ax.set_ylabel("malignant pseudobulk (log1p CP10K)")
        ax.set_title(f"{g}\nMWU p={st['p']:.3f}, AUC={st['auc']:.2f}")
    fig.suptitle("GSE205335 leftover: malignant TACSTD2/CLDN4 vs RECIST (NOT neoadjuvant MPR)", y=1.03)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "GSE205335_malignant_boxplots.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # NSCLC-only
    fig, axes = plt.subplots(1, len(GENES), figsize=(4.2 * len(GENES), 4.2))
    for ax, g in zip(axes, GENES):
        col = f"{g}_cp10k_log1p"
        groups = [nsclc.loc[nsclc["responder"] == 0, col].values,
                  nsclc.loc[nsclc["responder"] == 1, col].values]
        bp = ax.boxplot(groups, widths=0.6, patch_artist=True, showfliers=False)
        for patch, c in zip(bp["boxes"], ["#8aa1b1", "#d1495b"]):
            patch.set_facecolor(c)
            patch.set_alpha(0.55)
        rng = np.random.default_rng(4)
        for i, vals in enumerate(groups, 1):
            if len(vals):
                ax.scatter(rng.normal(i, 0.06, len(vals)), vals, s=22, color="#22303c", zorder=3)
        st = mwu_stats(groups[1], groups[0])
        ax.set_xticks([1, 2])
        ax.set_xticklabels([f"SD/PD\n(n={len(groups[0])})", f"PR\n(n={len(groups[1])})"])
        ax.set_ylabel("malignant pseudobulk (log1p CP10K)")
        ax.set_title(f"{g} (NSCLC only)\nMWU p={st['p']:.3f}, AUC={st['auc']:.2f}")
    fig.suptitle("GSE205335 NSCLC (ADC/SQ) malignant vs RECIST", y=1.03)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "GSE205335_malignant_NSCLC_boxplots.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
