#!/usr/bin/env python3
"""GSE50927 -- mouse whole-lung Cldn4 knockout vs wild type (+/- ventilator injury).

This is the only public LUNG CLDN4 perturbation transcriptome found (see
notes/kdko_catalog.tsv). Samples (GEO SOFT record):
    GSM1232580 WT no VILI      GSM1232581 Cldn4 KO no VILI
    GSM1232582 WT VILI         GSM1232583 Cldn4 KO VILIlow
                               GSM1232584 Cldn4 KO VILIhigh
Whole lung, mixed 129S6/C57BL/6/BALB/c background, ventilator-induced lung
injury (VILI) 40 cmH2O for 2 h.

*** THE CENTRAL LIMITATION OF THIS SERIES ***
n = 1 mouse per group. GEO holds no count matrix, only four submitter-supplied
edgeR tables. An edgeR fit on an unreplicated design must assume a dispersion
rather than estimate it, so the PValue/FDR columns in those tables are NOT
replicate-based evidence. They are reported here verbatim, always labelled
`submitter_edgeR_*`, and no claim in the write-up rests on them alone. The
fold-change ranking is still usable for descriptive gene-set enrichment, which
is what the competitive tests below do (their n is a number of genes).

Sign convention verified empirically: Cldn4 logFC = -6.06 in the WT-vs-KO table,
so logFC = KO - WT in that table, and KO - WT within VILI in the two VILI tables.
"""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import lib_kdko as L

GSE = "GSE50927"
TABLES = {
    # file -> (contrast label, what logFC means, n per side)
    "GSE50927_Cldn4lungWTvsKOgenes.csv.gz": (
        "Cldn4KO_vs_WT_naive_lung", "Cldn4 KO minus WT, no ventilator injury", 1, 1),
    "GSE50927_VILIwtkoloGenes.csv.gz": (
        "Cldn4KO_VILIlow_vs_WT_VILI", "Cldn4 KO (low-injury) minus WT, both VILI", 1, 1),
    "GSE50927_VILIwtkohiGenes.csv.gz": (
        "Cldn4KO_VILIhigh_vs_WT_VILI", "Cldn4 KO (high-injury) minus WT, both VILI", 1, 1),
    "GSE50927_VILIwtGenes.csv.gz": (
        "WT_VILI_vs_WT_naive", "VILI minus no-VILI within WT (injury reference, "
        "not a Cldn4 perturbation)", 1, 1),
}
FOCUS = ["Cldn4", "Tacstd2", "Epcam", "Cldn1", "Cldn3", "Cldn5", "Cldn7",
         "Cldn18", "Ocln", "Tjp1", "Tjp2", "Cdh1", "Vim", "Cd274", "B2m",
         "H2-K1", "Tap1", "Psmb9", "Nlrc5", "Stat1", "Stat2", "Irf7", "Isg15",
         "Ifit1", "Ifit3", "Mx1", "Oasl2", "Cxcl9", "Cxcl10", "Ccl5", "Ido1",
         "Gbp2", "Cd8a", "Cd3e", "Gzmb", "Prf1", "Ifng", "Ptprc", "Nkg7",
         "Il6", "Tnf", "Cxcl1", "Cxcl2"]


def load_table(fname: str) -> pd.DataFrame:
    df = pd.read_csv(os.path.join(L.DATA, GSE, fname))
    symcol = "Marker.Symbol" if "Marker.Symbol" in df.columns else "GeneSymbol"
    df = df[df[symcol].notna()]
    df = (df.assign(_absfc=df["logFC"].abs())
            .sort_values("logCPM", ascending=False)
            .drop_duplicates(symcol)
            .set_index(symcol))
    out = df[["logFC", "logCPM", "PValue", "FDR"]].copy()
    out.columns = ["logFC", "logCPM", "submitter_edgeR_PValue",
                   "submitter_edgeR_FDR"]
    return out


def main() -> None:
    os.makedirs(L.RES, exist_ok=True)
    sets = L.load_gene_sets("mouse")
    all_focus, all_sets, validity = {}, [], []

    for fname, (tag, meaning, n_a, n_b) in TABLES.items():
        t = load_table(fname)
        t.insert(0, "contrast", tag)
        t.insert(1, "logFC_meaning", meaning)
        t.insert(2, "n_per_group", f"{n_a} vs {n_b} (UNREPLICATED)")
        L.write_tsv(t.sort_values("submitter_edgeR_PValue"),
                    f"{GSE}_gene_stats_{tag}.tsv")
        all_focus[tag] = t.reindex([g for g in FOCUS if g in t.index])

        comp = L.competitive_from_ranking(t["logFC"], sets)
        comp.insert(0, "contrast", tag)
        comp["inference_unit"] = "genes (NOT biological replicates; n=1 per group)"
        all_sets.append(comp)

        if "Cldn4" in t.index:
            r = t.loc["Cldn4"]
            validity.append({
                "contrast": tag, "gene": "Cldn4", "logFC": float(r["logFC"]),
                "pct_of_reference": float(100 * 2 ** float(r["logFC"])),
                "logCPM": float(r["logCPM"]),
                "submitter_edgeR_PValue": float(r["submitter_edgeR_PValue"]),
                "submitter_edgeR_FDR": float(r["submitter_edgeR_FDR"]),
                "n_per_group": f"{n_a} vs {n_b}",
                "caveat": "unreplicated design; edgeR p-value assumes a "
                          "dispersion rather than estimating it"})
            print(f"  {tag}: Cldn4 logFC {r['logFC']:+.2f} "
                  f"({100 * 2 ** float(r['logFC']):.2f}% of reference)")

    L.write_tsv(pd.DataFrame(validity), f"{GSE}_ko_validity.tsv", index=False)
    L.write_tsv(pd.concat(all_sets, ignore_index=True),
                f"{GSE}_set_stats_all_contrasts.tsv", index=False)

    focus_wide = pd.concat(
        {k: v[["logFC", "submitter_edgeR_PValue", "submitter_edgeR_FDR"]]
         for k, v in all_focus.items()}, axis=1)
    L.write_tsv(focus_wide, f"{GSE}_focus_genes_all_contrasts.tsv")

    make_figures(all_focus, pd.concat(all_sets, ignore_index=True), sets)
    print("done")


def make_figures(all_focus, setdf, sets) -> None:
    naive = all_focus["Cldn4KO_vs_WT_naive_lung"]
    full_naive = load_table("GSE50927_Cldn4lungWTvsKOgenes.csv.gz")
    fig, axes = plt.subplots(2, 3, figsize=(17, 10))

    ax = axes[0, 0]
    fdr = full_naive["submitter_edgeR_FDR"].fillna(1.0)
    ax.scatter(full_naive["logFC"],
               -np.log10(full_naive["submitter_edgeR_PValue"].clip(lower=1e-300)),
               s=5, c=np.where(fdr < 0.05, "#cc3311", "#bbbbbb"), alpha=.55)
    for g in ["Cldn4", "Cldn18", "Tacstd2", "Isg15", "Ifit1", "Cxcl9", "Cd8a"]:
        if g in full_naive.index:
            ax.annotate(g, (full_naive.loc[g, "logFC"],
                            -np.log10(max(full_naive.loc[
                                g, "submitter_edgeR_PValue"], 1e-300))), fontsize=8)
    ax.axvline(0, lw=.5, c="k")
    ax.set_xlabel("logFC (Cldn4 KO − WT), naive lung")
    ax.set_ylabel("-log10 submitter edgeR p (UNREPLICATED, n=1 vs 1)")
    ax.set_title(f"GSE50927 naive lung\n{int((fdr<0.05).sum())} of "
                 f"{len(full_naive)} genes at submitter FDR<0.05")

    ax = axes[0, 1]
    genes = [g for g in ["Cldn4", "Cldn1", "Cldn3", "Cldn5", "Cldn7", "Cldn18",
                         "Tacstd2", "Epcam", "Ocln", "Tjp1"] if g in naive.index]
    fc = naive.loc[genes, "logFC"]
    ax.barh(np.arange(len(genes)), fc,
            color=["#cc3311" if v > 0 else "#4477aa" for v in fc])
    ax.set_yticks(np.arange(len(genes)))
    ax.set_yticklabels(genes, fontsize=9)
    ax.invert_yaxis()
    ax.axvline(0, c="k", lw=.6)
    ax.set_xlabel("logFC (Cldn4 KO − WT)")
    ax.set_title("Claudins / junctions / TACSTD2\nnaive Cldn4 KO lung")
    for i, g in enumerate(genes):
        if naive.loc[g, "submitter_edgeR_FDR"] < 0.05:
            ax.text(fc[g] + (.2 if fc[g] > 0 else -.2), i, "†", va="center",
                    ha="left" if fc[g] > 0 else "right", fontsize=11)

    ax = axes[0, 2]
    genes = [g for g in ["Isg15", "Ifit1", "Ifit3", "Mx1", "Oasl2", "Stat1",
                         "Irf7", "Cxcl9", "Cxcl10", "Gbp2", "B2m", "H2-K1",
                         "Tap1", "Psmb9", "Cd274", "Cd8a", "Cd3e", "Gzmb",
                         "Ifng", "Ptprc"] if g in naive.index]
    fc = naive.loc[genes, "logFC"]
    ax.barh(np.arange(len(genes)), fc,
            color=["#cc3311" if v > 0 else "#4477aa" for v in fc])
    ax.set_yticks(np.arange(len(genes)))
    ax.set_yticklabels(genes, fontsize=8)
    ax.invert_yaxis()
    ax.axvline(0, c="k", lw=.6)
    ax.set_xlabel("logFC (Cldn4 KO − WT)")
    ax.set_title("IFN / APM / immune genes\nnaive Cldn4 KO lung († submitter FDR<0.05)")
    for i, g in enumerate(genes):
        if naive.loc[g, "submitter_edgeR_FDR"] < 0.05:
            ax.text(fc[g] + (.08 if fc[g] > 0 else -.08), i, "†", va="center",
                    ha="left" if fc[g] > 0 else "right", fontsize=11)

    ax = axes[1, 0]
    order = ["TARGET_GENES", "CLAUDIN_FAMILY", "TIGHT_JUNCTION_CORE",
             "EPITHELIAL_IDENTITY", "IFN_ALPHA_TYPE1", "IFN_GAMMA_TYPE2",
             "ANTIGEN_PRESENTATION_MHC1", "T_CELL_CYTOTOXICITY",
             "T_CELL_INFILTRATION_BROAD", "CHEMOKINE_T_RECRUITMENT",
             "ICI_CHECKPOINT_GENES", "NFKB_INFLAMMATORY"]
    s = setdf[setdf["contrast"] == "Cldn4KO_vs_WT_naive_lung"].set_index(
        "set_name").reindex(order)
    y = np.arange(len(order))
    ax.barh(y, s["median_stat_in_set"].fillna(0),
            color=["#cc3311" if v > 0 else "#4477aa"
                   for v in s["median_stat_in_set"].fillna(0)])
    ax.set_yticks(y)
    ax.set_yticklabels([n.replace("_", " ").title() for n in order], fontsize=8)
    ax.invert_yaxis()
    ax.axvline(0, c="k", lw=.6)
    ax.set_xlabel("median logFC of set genes (KO − WT)")
    ax.set_title("Naive lung: set shift\n* = competitive q<0.05 (gene-level)")
    for i, (dd, q) in enumerate(zip(s["median_stat_in_set"], s["q_BH"])):
        if pd.notna(q) and q < 0.05:
            ax.text(dd + (.02 if dd > 0 else -.02), i, "*", va="center",
                    ha="left" if dd > 0 else "right", fontsize=13)

    ax = axes[1, 1]
    tags = ["Cldn4KO_vs_WT_naive_lung", "Cldn4KO_VILIlow_vs_WT_VILI",
            "Cldn4KO_VILIhigh_vs_WT_VILI"]
    show = ["IFN_ALPHA_TYPE1", "IFN_GAMMA_TYPE2", "CLAUDIN_FAMILY",
            "TIGHT_JUNCTION_CORE", "T_CELL_CYTOTOXICITY", "NFKB_INFLAMMATORY"]
    w = 0.26
    for i, tag in enumerate(tags):
        s = setdf[setdf["contrast"] == tag].set_index("set_name").reindex(show)
        ax.bar(np.arange(len(show)) + (i - 1) * w, s["median_stat_in_set"].fillna(0),
               width=w, label=tag.replace("Cldn4KO_", "").replace("_", " "))
    ax.set_xticks(np.arange(len(show)))
    ax.set_xticklabels([n.replace("_", "\n") for n in show], fontsize=7)
    ax.axhline(0, c="k", lw=.6)
    ax.set_ylabel("median logFC of set genes")
    ax.set_title("Cldn4 KO effect across injury states")
    ax.legend(frameon=False, fontsize=7)

    ax = axes[1, 2]
    a = load_table("GSE50927_VILIwtkoloGenes.csv.gz")["logFC"]
    b = load_table("GSE50927_VILIwtkohiGenes.csv.gz")["logFC"]
    common = a.index.intersection(b.index)
    from scipy import stats as st
    r, pr = st.pearsonr(a[common], b[common])
    ax.scatter(a[common], b[common], s=4, alpha=.2, c="#999999")
    for g in ["Cldn4", "Isg15", "Cxcl9", "Tacstd2", "Cldn18"]:
        if g in common:
            ax.scatter(a[g], b[g], s=70, c="#cc3311", zorder=5)
            ax.annotate(g, (a[g], b[g]), fontsize=9)
    ax.axhline(0, lw=.5, c="k")
    ax.axvline(0, lw=.5, c="k")
    ax.set_xlabel("logFC KO VILIlow − WT VILI")
    ax.set_ylabel("logFC KO VILIhigh − WT VILI")
    ax.set_title(f"Two Cldn4 KO VILI mice vs the same WT VILI mouse\n"
                 f"Pearson r={r:.3f}, n={len(common)} genes "
                 "(shared reference inflates r)")

    fig.suptitle("GSE50927 — mouse whole-lung Cldn4 knockout, the only public "
                 "lung CLDN4 perturbation; n=1 animal per group (UNREPLICATED)",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    out = os.path.join(L.RES, f"fig_{GSE}_cldn4ko_lung.png")
    fig.savefig(out, dpi=160)
    print(f"  wrote {os.path.basename(out)}")
    plt.close(fig)


if __name__ == "__main__":
    main()
