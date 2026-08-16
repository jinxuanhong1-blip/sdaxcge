#!/usr/bin/env python3
"""GSE15212 -- human SW480 colon adenocarcinoma, TACSTD2 siRNA knockdown.

Design (from the GEO series matrix; the series is an RNAi signature panel in which
TACSTD2 is one of five targeted genes):
    TACSTD2 siRNA #1, 72 h : GSM379987-89  (TACSTD2.1.72h.4/5/6)   n=3
    TACSTD2 siRNA #4, 72 h : GSM379990-92  (TACSTD2.4.72h.4/5/6)   n=3
    negative-control siRNA #2, 72 h : GSM379960-68 (Neg.2.72h.1-9)  n=9
    mock (no siRNA), 72 h           : GSM379957-59 (Neg.0.72h.1-3)  n=3

Two independent siRNAs against the same gene is the strongest available design
for on-target attribution in a knockdown experiment, so both the pooled contrast
and the per-siRNA contrasts are reported. Platform GPL4133 (Agilent 4x44K),
log2 quantile-normalised by the submitters.
"""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import lib_kdko as L

GSE = "GSE15212"
SI1 = ["GSM379987", "GSM379988", "GSM379989"]
SI4 = ["GSM379990", "GSM379991", "GSM379992"]
KD = SI1 + SI4
NEG2 = [f"GSM3799{n}" for n in range(60, 69)]          # Neg.2.72h.1-9
NEG2_MATCHED = ["GSM379963", "GSM379964", "GSM379965"]  # Neg.2.72h.4/5/6
MOCK = ["GSM379957", "GSM379958", "GSM379959"]         # Neg.0.72h.1-3


def load() -> pd.DataFrame:
    expr, _ = L.read_series_matrix(
        os.path.join(L.DATA, GSE, f"{GSE}_series_matrix.txt.gz"))
    plat = L.read_platform_table(os.path.join(L.DATA, GSE, f"{GSE}_family.soft.gz"))
    plat = plat[plat["CONTROL_TYPE"].astype(str).str.lower().isin(["false", ""])]
    sym = plat["GENE_SYMBOL"].replace("", np.nan)
    expr.index = expr.index.astype(str)
    sym.index = sym.index.astype(str)
    need = KD + NEG2 + MOCK
    missing = [c for c in need if c not in expr.columns]
    if missing:
        raise SystemExit(f"samples missing: {missing}")
    mat = L.collapse_to_symbol(expr[need], sym)
    # Agilent log2 signal: drop probes near the detection floor in all samples
    return mat[mat.max(axis=1) > mat.stack().quantile(0.15)]


def contrast(mat, kd, ctrl, tag, sets, out_prefix) -> dict:
    gs = L.per_gene_tests(mat, kd, ctrl, "KD", "CTRL")
    gs.insert(0, "contrast", tag)
    L.write_tsv(gs.sort_values("p_welch"), f"{out_prefix}_gene_stats_{tag}.tsv")
    ss = L.set_level_tests(mat, gs, sets, kd, ctrl, "TACSTD2_KD", "control")
    ss.insert(0, "contrast", tag)
    L.write_tsv(ss, f"{out_prefix}_set_stats_{tag}.tsv", index=False)
    v = mat.loc["TACSTD2"]
    t, p = stats.ttest_ind(v[kd], v[ctrl], equal_var=False)
    rec = {"contrast": tag, "n_KD": len(kd), "n_control": len(ctrl),
           "TACSTD2_mean_log2_KD": float(v[kd].mean()),
           "TACSTD2_mean_log2_control": float(v[ctrl].mean()),
           "TACSTD2_log2FC": float(v[kd].mean() - v[ctrl].mean()),
           "TACSTD2_pct_of_control": float(
               100 * 2 ** (v[kd].mean() - v[ctrl].mean())),
           "TACSTD2_welch_t": float(t), "TACSTD2_welch_p": float(p),
           "n_genes_q_lt_0.05": int((gs["q_welch_BH"] < 0.05).sum()),
           "n_genes_tested": int(len(gs))}
    return rec, gs, ss


def main() -> None:
    os.makedirs(L.RES, exist_ok=True)
    mat = load()
    sets = L.load_gene_sets("human")
    print(f"{GSE}: {mat.shape[0]} genes x {mat.shape[1]} samples")

    validity, gstats, sstats = [], {}, {}
    for tag, kd, ctrl in [
        ("KDpooled_vs_negctrl", KD, NEG2),
        ("siRNA1_vs_negctrl", SI1, NEG2),
        ("siRNA4_vs_negctrl", SI4, NEG2),
        ("KDpooled_vs_matchedreps", KD, NEG2_MATCHED),
        ("KDpooled_vs_mock", KD, MOCK),
    ]:
        rec, gs, ss = contrast(mat, kd, ctrl, tag, sets, GSE)
        validity.append(rec)
        gstats[tag], sstats[tag] = gs, ss
        print(f"  {tag}: TACSTD2 {rec['TACSTD2_log2FC']:+.2f} log2 "
              f"({rec['TACSTD2_pct_of_control']:.0f}% of control), "
              f"p={L.fmt_p(rec['TACSTD2_welch_p'])}, "
              f"{rec['n_genes_q_lt_0.05']} genes at q<0.05")
    L.write_tsv(pd.DataFrame(validity), f"{GSE}_kd_validity.tsv", index=False)

    # on-target reproducibility across the two independent siRNAs
    a = gstats["siRNA1_vs_negctrl"]["log2FC"]
    b = gstats["siRNA4_vs_negctrl"]["log2FC"]
    common = a.index.intersection(b.index)
    r, pr = stats.pearsonr(a[common], b[common])
    rho, prho = stats.spearmanr(a[common], b[common])
    pd.DataFrame([{"comparison": "siRNA1_vs_siRNA4_log2FC",
                   "pearson_r": r, "pearson_p": pr, "spearman_rho": rho,
                   "spearman_p": prho, "n_genes": len(common)}]).to_csv(
        os.path.join(L.RES, f"{GSE}_sirna_concordance.tsv"), sep="\t", index=False)
    print(f"  siRNA1 vs siRNA4 log2FC concordance: r={r:.3f} "
          f"(p={L.fmt_p(pr)}), n={len(common)} genes")

    focus = ["TACSTD2", "CLDN4", "CLDN1", "CLDN3", "CLDN7", "EPCAM", "OCLN",
             "TJP1", "TJP2", "CDH1", "VIM", "DSG2", "CD274", "B2M", "HLA-A",
             "HLA-B", "TAP1", "PSMB9", "STAT1", "STAT2", "IRF7", "IRF9",
             "ISG15", "IFIT1", "IFIT3", "MX1", "MX2", "OAS1", "OAS2", "OASL",
             "CXCL9", "CXCL10", "IDO1", "GBP1", "GBP2", "NFKBIA", "CXCL8", "IL6"]
    L.write_tsv(gstats["KDpooled_vs_negctrl"].reindex(
        [g for g in focus if g in mat.index]), f"{GSE}_focus_genes_pooledKD.tsv")

    make_figures(mat, gstats, sstats, sets, r)
    print("done")


def make_figures(mat, gstats, sstats, sets, sirna_r) -> None:
    gs = gstats["KDpooled_vs_negctrl"]
    ss = sstats["KDpooled_vs_negctrl"].set_index("set_name")
    fig, axes = plt.subplots(2, 3, figsize=(17, 10))

    ax = axes[0, 0]
    groups = [("neg-ctrl siRNA", NEG2, "#4477aa"), ("mock", MOCK, "#88ccee"),
              ("TACSTD2 siRNA#1", SI1, "#cc3311"),
              ("TACSTD2 siRNA#4", SI4, "#ee7733")]
    for i, (lab, cols, c) in enumerate(groups):
        ax.scatter(np.full(len(cols), i) + np.random.uniform(-.08, .08, len(cols)),
                   mat.loc["TACSTD2", cols], color=c, s=48)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([g[0] for g in groups], rotation=25, ha="right", fontsize=8)
    ax.set_ylabel("TACSTD2 log2 signal (Agilent, quantile-norm.)")
    ax.set_title("GSE15212 SW480: knockdown validity")

    ax = axes[0, 1]
    q = gs["q_welch_BH"].fillna(1.0)
    ax.scatter(gs["log2FC"], -np.log10(gs["p_welch"].clip(lower=1e-300)), s=5,
               c=np.where(q < 0.05, "#cc3311", "#bbbbbb"), alpha=.5)
    for g in ["TACSTD2", "CLDN4", "CLDN1", "CLDN7", "ISG15", "IFIT1", "STAT1",
              "CD274", "B2M"]:
        if g in gs.index and pd.notna(gs.loc[g, "p_welch"]):
            ax.annotate(g, (gs.loc[g, "log2FC"],
                            -np.log10(max(gs.loc[g, "p_welch"], 1e-300))), fontsize=8)
    ax.axvline(0, lw=.5, c="k")
    ax.set_xlabel("log2FC (TACSTD2 KD / neg-ctrl siRNA)")
    ax.set_ylabel("-log10 Welch p")
    ax.set_title(f"Pooled KD (n=6) vs neg-ctrl (n=9)\n{int((q<0.05).sum())} of "
                 f"{len(gs)} genes at BH q<0.05")

    ax = axes[0, 2]
    order = ["TARGET_GENES", "CLAUDIN_FAMILY", "TIGHT_JUNCTION_CORE",
             "DESMOSOME_ADHERENS", "EPITHELIAL_IDENTITY", "EMT_CORE",
             "IFN_ALPHA_TYPE1", "IFN_GAMMA_TYPE2", "ANTIGEN_PRESENTATION_MHC1",
             "CHEMOKINE_T_RECRUITMENT", "ICI_CHECKPOINT_GENES",
             "TIS_18GENE_AYERS", "NFKB_INFLAMMATORY"]
    s = ss.reindex(order)
    y = np.arange(len(order))
    ax.barh(y, s["sample_delta"].fillna(0),
            color=["#cc3311" if v > 0 else "#4477aa"
                   for v in s["sample_delta"].fillna(0)])
    ax.set_yticks(y)
    ax.set_yticklabels([n.replace("_", " ").title() for n in order], fontsize=8)
    ax.invert_yaxis()
    ax.axvline(0, c="k", lw=.6)
    ax.set_xlabel("Δ mean-z set score (KD − control)")
    ax.set_title("Set score shift; * = sample-level p<0.05")
    for i, (dd, p) in enumerate(zip(s["sample_delta"], s["sample_p"])):
        if pd.notna(p) and p < 0.05:
            ax.text(dd + (.02 if dd > 0 else -.02), i, "*", va="center",
                    ha="left" if dd > 0 else "right", fontsize=13)

    ax = axes[1, 0]
    a = gstats["siRNA1_vs_negctrl"]["log2FC"]
    b = gstats["siRNA4_vs_negctrl"]["log2FC"]
    common = a.index.intersection(b.index)
    ax.scatter(a[common], b[common], s=4, alpha=.2, c="#888888")
    for g in ["TACSTD2", "CLDN4", "CLDN1", "CLDN7", "EPCAM"]:
        if g in common:
            ax.scatter(a[g], b[g], s=70, c="#cc3311", zorder=5)
            ax.annotate(g, (a[g], b[g]), fontsize=9)
    ax.axhline(0, lw=.5, c="k")
    ax.axvline(0, lw=.5, c="k")
    ax.set_xlabel("log2FC siRNA #1 vs neg-ctrl")
    ax.set_ylabel("log2FC siRNA #4 vs neg-ctrl")
    ax.set_title(f"Two independent siRNAs\nPearson r={sirna_r:.3f}, "
                 f"n={len(common)} genes")

    ax = axes[1, 1]
    show = ["CLAUDIN_FAMILY", "TIGHT_JUNCTION_CORE", "IFN_ALPHA_TYPE1",
            "IFN_GAMMA_TYPE2", "ANTIGEN_PRESENTATION_MHC1", "EPITHELIAL_IDENTITY"]
    for i, name in enumerate(show):
        sc = L.set_score_matrix(mat, sets[name])
        if sc is None:
            continue
        ax.scatter(np.full(len(NEG2), i) - .18, sc[NEG2], color="#4477aa", s=26,
                   label="neg-ctrl" if i == 0 else None)
        ax.scatter(np.full(len(KD), i) + .18, sc[KD], color="#cc3311", s=26,
                   label="TACSTD2 KD" if i == 0 else None)
    ax.set_xticks(range(len(show)))
    ax.set_xticklabels([n.replace("_", "\n") for n in show], fontsize=7)
    ax.axhline(0, lw=.5, c="k")
    ax.set_ylabel("per-sample mean-z score")
    ax.set_title("Per-sample set scores")
    ax.legend(frameon=False, fontsize=9)

    ax = axes[1, 2]
    genes = [g for g in ["TACSTD2", "CLDN1", "CLDN3", "CLDN4", "CLDN7", "EPCAM",
                         "OCLN", "TJP1", "TJP2", "CDH1", "DSG2", "ISG15",
                         "IFIT1", "MX1", "STAT1", "B2M", "HLA-A", "CD274"]
             if g in gs.index]
    fc = gs.loc[genes, "log2FC"]
    ax.barh(np.arange(len(genes)), fc,
            color=["#cc3311" if v > 0 else "#4477aa" for v in fc])
    ax.set_yticks(np.arange(len(genes)))
    ax.set_yticklabels(genes, fontsize=8)
    ax.invert_yaxis()
    ax.axvline(0, c="k", lw=.6)
    ax.set_xlabel("log2FC (TACSTD2 KD / neg-ctrl)")
    ax.set_title("Junction / IFN / APM genes\n* = Welch p<0.05")
    for i, g in enumerate(genes):
        if pd.notna(gs.loc[g, "p_welch"]) and gs.loc[g, "p_welch"] < 0.05:
            ax.text(fc[g] + (.02 if fc[g] > 0 else -.02), i, "*", va="center",
                    ha="left" if fc[g] > 0 else "right", fontsize=12)

    fig.suptitle("GSE15212 — human SW480 TACSTD2 siRNA knockdown, 72 h, "
                 "two independent siRNAs (n=6) vs negative-control siRNA (n=9)",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    out = os.path.join(L.RES, f"fig_{GSE}_tacstd2_kd_SW480.png")
    fig.savefig(out, dpi=160)
    print(f"  wrote {os.path.basename(out)}")
    plt.close(fig)


if __name__ == "__main__":
    np.random.seed(0)
    main()
