#!/usr/bin/env python3
"""Cross-dataset synthesis focused on IFN and MHC-I after true CLDN4 / TACSTD2 KD/KO.

Reads only the per-dataset tables already written by scripts 04-09. Does not
re-download anything. Produces:
  results/kdko/cross_ifn_mhci_set_summary.tsv
  results/kdko/cross_ifn_mhci_gene_log2FC.tsv
  results/kdko/cross_reciprocal_CLDN4_TACSTD2.tsv
  results/kdko/fig_cross_ifn_mhci.png
"""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import lib_kdko as L

# Contrasts that are true CLDN4 or TACSTD2 perturbations (not comparators).
# `set_file` / `set_filter` point at the already-written set-level table.
# `kind` is either 'sample' (replicate-based) or 'gene' (competitive ranking).
CONTRASTS = [
    dict(id="GSE334497_4T1_Trop2KO",
         accession="GSE334497", gene="TACSTD2", species="mouse",
         system="4T1 mammary tumour, BALB/c immunocompetent",
         n="5 KO vs 5 WT", kind="sample",
         set_file="GSE334497_set_stats_KOvsWT.tsv", set_filter=None,
         gene_file="GSE334497_focus_genes_KOvsWT.tsv",
         fc_col="log2FC", p_col="p_welch", q_col="q_welch_BH"),
    dict(id="GSE289287_T47D_TACSTD2KO_xg",
         accession="GSE289287", gene="TACSTD2", species="human",
         system="T-47D xenograft, NRG immunodeficient",
         n="4 KO vs 3 WT", kind="sample",
         set_file="GSE289287_set_stats_Trop2KO_xenograft.tsv", set_filter=None,
         gene_file="GSE289287_focus_genes_Trop2KO_xenograft.tsv",
         fc_col="log2FC_KO_vs_WT", p_col="p_deseq2", q_col="q_deseq2_BH"),
    dict(id="GSE15212_SW480_TACSTD2_siRNA",
         accession="GSE15212", gene="TACSTD2", species="human",
         system="SW480 colon adenocarcinoma, 72 h siRNA",
         n="6 KD vs 9 neg-ctrl", kind="sample",
         set_file="GSE15212_set_stats_KDpooled_vs_negctrl.tsv", set_filter=None,
         gene_file="GSE15212_focus_genes_pooledKD.tsv",
         fc_col="log2FC", p_col="p_welch", q_col="q_welch_BH"),
    dict(id="GSE207704_T47D_MCF7_CLDN4KO",
         accession="GSE207704", gene="CLDN4", species="human",
         system="T47D + MCF7 CRISPR KO (group-mean FPKM)",
         n="2 lines, replicates collapsed", kind="gene",
         set_file="GSE207704_set_stats_CLDN4KO.tsv",
         set_filter=("contrast", "mean_of_both_lines"),
         gene_file="GSE207704_focus_genes_CLDN4KO.tsv",
         fc_col="log2FC_mean_bothlines", p_col=None, q_col=None),
    dict(id="GSE50927_lung_Cldn4KO_naive",
         accession="GSE50927", gene="CLDN4", species="mouse",
         system="whole lung, naive (no VILI)",
         n="1 KO vs 1 WT", kind="gene",
         set_file="GSE50927_set_stats_all_contrasts.tsv",
         set_filter=("contrast", "Cldn4KO_vs_WT_naive_lung"),
         gene_file="GSE50927_focus_genes_all_contrasts.tsv",
         fc_col=None, p_col=None, q_col=None),  # multi-index, handled below
    dict(id="GSE22493_SKOV3_CLDN4_siRNA",
         accession="GSE22493", gene="CLDN4", species="human",
         system="SKOV-3 vs CLDN4-OE reference, 2-colour array",
         n="3 arrays", kind="gene",
         set_file="GSE22493_set_stats_CLDN4KD.tsv", set_filter=None,
         gene_file="GSE22493_focus_genes_CLDN4KD.tsv",
         fc_col="mean_log2ratio_KD_over_control",
         p_col="p_onesample_t", q_col="q_onesample_BH"),
]

IFN_SETS = ["IFN_ALPHA_TYPE1", "IFN_GAMMA_TYPE2", "ANTIGEN_PRESENTATION_MHC1"]
KEY_GENES_H = ["TACSTD2", "CLDN4", "ISG15", "STAT1", "STAT2", "IRF7",
               "IFIT1", "IFIT3", "MX1", "OAS1", "B2M", "HLA-A", "TAP1",
               "PSMB9", "NLRC5", "CXCL9", "CXCL10", "CD274"]
KEY_GENES_M = ["Tacstd2", "Cldn4", "Isg15", "Stat1", "Stat2", "Irf7",
               "Ifit1", "Ifit3", "Mx1", "Oas1", "B2m", "H2-K1", "Tap1",
               "Psmb9", "Nlrc5", "Cxcl9", "Cxcl10", "Cd274"]
DISPLAY = ["TACSTD2", "CLDN4", "ISG15", "STAT1", "STAT2", "IRF7",
           "IFIT1", "IFIT3", "MX1", "OAS1", "B2M", "HLA-A/H2-K1",
           "TAP1", "PSMB9", "NLRC5", "CXCL9", "CXCL10", "CD274"]


def _read_sets(c: dict) -> pd.DataFrame:
    df = pd.read_csv(os.path.join(L.RES, c["set_file"]), sep="\t")
    if c["set_filter"]:
        k, v = c["set_filter"]
        df = df[df[k] == v].copy()
    return df.set_index("set_name")


def _read_genes(c: dict) -> pd.DataFrame:
    path = os.path.join(L.RES, c["gene_file"])
    if c["id"].startswith("GSE50927"):
        df = pd.read_csv(path, sep="\t", header=[0, 1], index_col=0)
        sub = df["Cldn4KO_vs_WT_naive_lung"].copy()
        sub.columns = ["log2FC", "p", "q"]
        return sub
    df = pd.read_csv(path, sep="\t", index_col=0)
    return df


def extract_set_row(c: dict, set_name: str) -> dict:
    s = _read_sets(c)
    if set_name not in s.index:
        return {"contrast": c["id"], "set_name": set_name, "n_genes_measured": 0}
    r = s.loc[set_name]
    out = {"contrast": c["id"], "accession": c["accession"],
           "target_gene": c["gene"], "species": c["species"],
           "system": c["system"], "n": c["n"], "set_name": set_name,
           "inference": c["kind"],
           "n_genes_measured": int(r.get("n_genes_measured", np.nan)
                                   if pd.notna(r.get("n_genes_measured", np.nan))
                                   else 0)}
    if c["kind"] == "sample":
        out.update({
            "effect": float(r["sample_delta"]) if pd.notna(r.get("sample_delta")) else np.nan,
            "effect_label": "Δ mean-z set score (pert − control)",
            "p": float(r["sample_p"]) if pd.notna(r.get("sample_p")) else np.nan,
            "q": float(r["sample_q_BH"]) if pd.notna(r.get("sample_q_BH")) else np.nan,
            "cohens_d": float(r["sample_cohens_d"]) if pd.notna(r.get("sample_cohens_d")) else np.nan,
            "comp_median_log2FC": float(r["comp_median_log2FC_in_set"])
            if pd.notna(r.get("comp_median_log2FC_in_set")) else np.nan,
            "comp_p": float(r["comp_mannwhitney_p"])
            if pd.notna(r.get("comp_mannwhitney_p")) else np.nan,
        })
    else:
        # gene-level competitive ranking
        med = r.get("median_stat_in_set", r.get("median_log2FC_in_set"))
        p = r.get("mannwhitney_p", r.get("competitive_mannwhitney_p"))
        q = r.get("q_BH", r.get("competitive_q_BH"))
        out.update({
            "effect": float(med) if pd.notna(med) else np.nan,
            "effect_label": "median log2FC of set genes (competitive)",
            "p": float(p) if pd.notna(p) else np.nan,
            "q": float(q) if pd.notna(q) else np.nan,
            "cohens_d": np.nan,
            "comp_median_log2FC": float(med) if pd.notna(med) else np.nan,
            "comp_p": float(p) if pd.notna(p) else np.nan,
        })
    return out


def extract_gene_fc(c: dict) -> dict[str, dict]:
    g = _read_genes(c)
    mapping = dict(zip(KEY_GENES_H, KEY_GENES_M if c["species"] == "mouse"
                       else KEY_GENES_H))
    # HLA-A display uses H2-K1 in mouse
    out = {}
    for h, m in zip(KEY_GENES_H, KEY_GENES_M if c["species"] == "mouse"
                    else KEY_GENES_H):
        key = m if c["species"] == "mouse" else h
        if key not in g.index:
            out[h] = {"log2FC": np.nan, "p": np.nan, "q": np.nan, "symbol_used": key,
                      "present": False}
            continue
        row = g.loc[key]
        if c["fc_col"] and c["fc_col"] in row.index:
            fc = row[c["fc_col"]]
        elif "log2FC" in row.index:
            fc = row["log2FC"]
        else:
            fc = np.nan
        p = row[c["p_col"]] if c["p_col"] and c["p_col"] in row.index else (
            row["p"] if "p" in row.index else np.nan)
        q = row[c["q_col"]] if c["q_col"] and c["q_col"] in row.index else (
            row["q"] if "q" in row.index else np.nan)
        out[h] = {"log2FC": float(fc) if pd.notna(fc) else np.nan,
                  "p": float(p) if pd.notna(p) else np.nan,
                  "q": float(q) if pd.notna(q) else np.nan,
                  "symbol_used": key, "present": True}
    return out


def main() -> None:
    os.makedirs(L.RES, exist_ok=True)
    set_rows = [extract_set_row(c, s) for c in CONTRASTS for s in IFN_SETS]
    setdf = pd.DataFrame(set_rows)
    L.write_tsv(setdf, "cross_ifn_mhci_set_summary.tsv", index=False)

    gene_rows = []
    gene_mat = {}
    gene_p = {}
    for c in CONTRASTS:
        d = extract_gene_fc(c)
        gene_mat[c["id"]] = {h: d[h]["log2FC"] for h in KEY_GENES_H}
        gene_p[c["id"]] = {h: d[h]["p"] for h in KEY_GENES_H}
        for h, rec in d.items():
            gene_rows.append({"contrast": c["id"], "accession": c["accession"],
                              "target_gene": c["gene"], "query_symbol": h,
                              **rec})
    gdf = pd.DataFrame(gene_rows)
    L.write_tsv(gdf, "cross_ifn_mhci_gene_log2FC.tsv", index=False)

    rec_rows = []
    for c in CONTRASTS:
        d = extract_gene_fc(c)
        rec_rows.append({
            "contrast": c["id"], "accession": c["accession"],
            "perturbed_gene": c["gene"], "n": c["n"],
            "CLDN4_log2FC": d["CLDN4"]["log2FC"],
            "CLDN4_p": d["CLDN4"]["p"],
            "TACSTD2_log2FC": d["TACSTD2"]["log2FC"],
            "TACSTD2_p": d["TACSTD2"]["p"],
            "other_gene_tested": (
                "TACSTD2" if c["gene"] == "CLDN4" else "CLDN4"),
            "other_present": (
                d["TACSTD2"]["present"] if c["gene"] == "CLDN4"
                else d["CLDN4"]["present"]),
        })
    L.write_tsv(pd.DataFrame(rec_rows), "cross_reciprocal_CLDN4_TACSTD2.tsv",
                index=False)

    make_figure(setdf, gene_mat, gene_p)
    print("done")


def make_figure(setdf: pd.DataFrame, gene_mat: dict, gene_p: dict) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(14, 12),
                             gridspec_kw={"height_ratios": [1.05, 1.35]})

    # ---- panel A: IFN / MHC-I set effects ---------------------------------
    ax = axes[0]
    ids = [c["id"] for c in CONTRASTS]
    short = {
        "GSE334497_4T1_Trop2KO": "GSE334497\n4T1 Trop2 KO\nn=5 vs 5  *sample*",
        "GSE289287_T47D_TACSTD2KO_xg": "GSE289287\nT-47D Trop2 KO xg\nn=4 vs 3  *sample*",
        "GSE15212_SW480_TACSTD2_siRNA": "GSE15212\nSW480 TACSTD2 siRNA\nn=6 vs 9  *sample*",
        "GSE207704_T47D_MCF7_CLDN4KO": "GSE207704\nT47D/MCF7 CLDN4 KO\ngroup means  *gene*",
        "GSE50927_lung_Cldn4KO_naive": "GSE50927\nlung Cldn4 KO\nn=1 vs 1  *gene*",
        "GSE22493_SKOV3_CLDN4_siRNA": "GSE22493\nSKOV-3 CLDN4 siRNA\n3 arrays  *gene*",
    }
    colours = {"IFN_ALPHA_TYPE1": "#cc3311",
               "IFN_GAMMA_TYPE2": "#ee7733",
               "ANTIGEN_PRESENTATION_MHC1": "#4477aa"}
    labels = {"IFN_ALPHA_TYPE1": "IFN-I (type I ISGs)",
              "IFN_GAMMA_TYPE2": "IFN-II (IFN-γ / GBP / CXCL9-11)",
              "ANTIGEN_PRESENTATION_MHC1": "MHC-I antigen-processing"}
    x = np.arange(len(ids))
    w = 0.26
    for i, sname in enumerate(IFN_SETS):
        sub = setdf[setdf["set_name"] == sname].set_index("contrast").reindex(ids)
        ax.bar(x + (i - 1) * w, sub["effect"].fillna(0), width=w,
               color=colours[sname], label=labels[sname])
        for j, (eff, p, kind) in enumerate(
                zip(sub["effect"], sub["p"], sub["inference"])):
            if pd.notna(p) and p < 0.05:
                ax.text(x[j] + (i - 1) * w, eff + (0.04 if eff >= 0 else -0.04),
                        "*", ha="center",
                        va="bottom" if eff >= 0 else "top", fontsize=12)
    ax.axhline(0, c="k", lw=.7)
    ax.set_xticks(x)
    ax.set_xticklabels([short[i] for i in ids], fontsize=8)
    ax.set_ylabel("set effect\n(sample Δz  |  gene-level median log2FC)")
    ax.set_title("A. IFN and MHC-I programmes after true CLDN4 / TACSTD2 "
                 "perturbation\n* = p<0.05 on the test that matches the design "
                 "(sample-level Welch for n≥3 replicates; competitive "
                 "Mann-Whitney of gene fold-changes otherwise). "
                 "Left three = TACSTD2 loss; right three = CLDN4 loss.")
    ax.legend(frameon=False, fontsize=9, loc="upper right")

    # ---- panel B: per-gene heatmap ----------------------------------------
    ax = axes[1]
    mat = pd.DataFrame(gene_mat).T[KEY_GENES_H]
    mat.columns = DISPLAY
    mat.index = [short[i].split("\n")[0] + "  " +
                 short[i].split("\n")[1] for i in ids]
    vmax = 2.0
    im = ax.imshow(mat.to_numpy(dtype=float), aspect="auto", cmap="RdBu_r",
                   vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(DISPLAY)))
    ax.set_xticklabels(DISPLAY, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(len(mat)))
    ax.set_yticklabels(mat.index, fontsize=9)
    pmat = pd.DataFrame(gene_p).T[KEY_GENES_H]
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.iloc[i, j]
            if pd.isna(v):
                ax.text(j, i, "·", ha="center", va="center", color="#666666",
                        fontsize=9)
                continue
            p = pmat.iloc[i, j]
            mark = ""
            if pd.notna(p) and p < 0.05:
                mark = "*"
            ax.text(j, i, f"{v:+.2f}{mark}", ha="center", va="center",
                    fontsize=7,
                    color="white" if abs(v) > 1.1 else "black")
    fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01,
                 label="log2FC (perturbation / control)")
    ax.set_title("B. Per-gene log2FC. * = gene-level p<0.05 where a replicate- "
                 "or DESeq2/edgeR p-value exists. '·' = gene not measured. "
                 "GSE207704 has no per-gene p (replicates collapsed). "
                 "GSE50927 p-values are unreplicated edgeR (shown as * if "
                 "submitter FDR-unadjusted p<0.05).")

    fig.suptitle("CLDN4 / TACSTD2 true KD-KO series — IFN and MHC-I synthesis",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = os.path.join(L.RES, "fig_cross_ifn_mhci.png")
    fig.savefig(out, dpi=160)
    print(f"  wrote {os.path.basename(out)}")
    plt.close(fig)


if __name__ == "__main__":
    main()
