#!/usr/bin/env python3
"""GSE289287 -- human T-47D TACSTD2 (Trop-2) knockout xenografts vs WT xenografts.

Design (from the GEO SOFT record):
    WT  xenografts : GSM8788419-21  animals 2812, 2808, 2810              n=3
    Trop-2 KO xg   : GSM8788422-25  animals 2815, 2817, 2818, 2807        n=4
    Host: NRG (immunodeficient) female mice, mammary fat pad.

The submitters' own DESeq2 table is used for gene-level inference (it is a
proper negative-binomial fit, better than re-deriving a t-test from normalised
counts), and the per-sample normalised counts embedded in the same file are used
for sample-level set scores.

IMPORTANT interpretive limit: NRG mice lack T, B and NK cells, so this contrast
can only report the TUMOUR-INTRINSIC arm of the mechanism (claudin/tight-junction
and cell-intrinsic IFN wiring). It cannot report immune infiltration.
GEO holds no in vitro Trop-2 KO arm for this series (only WT and DSG2 KO cell
lines), so the DSG2 KO cell-line contrast is analysed alongside as a
junction-perturbation comparator, not as a TACSTD2 perturbation.
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

GSE = "GSE289287"
WT_COLS = ["OV.2812.RNA_normCounts", "OV.2808.RNA_normCounts", "OV.2810.RNA_normCounts"]
KO_COLS = ["OV.2815.RNA_normCounts", "OV.2817.RNA_normCounts",
           "OV.2818.RNA_normCounts", "OV.2807.RNA_normCounts"]
GSM = {"OV.2812": "GSM8788419", "OV.2808": "GSM8788420", "OV.2810": "GSM8788421",
       "OV.2815": "GSM8788422", "OV.2817": "GSM8788423", "OV.2818": "GSM8788424",
       "OV.2807": "GSM8788425"}


def load(fname: str, wt_cols: list[str], ko_cols: list[str]):
    df = pd.read_csv(os.path.join(L.DATA, GSE, fname), sep="\t")
    df = df[df["Feature_name"].notna()]
    df = df[df["biotype"] == "protein_coding"]
    # DESeq2 stats, collapsed to one row per symbol (keep highest baseMean)
    de = (df.sort_values("baseMean", ascending=False)
            .drop_duplicates("Feature_name")
            .set_index("Feature_name"))
    mat = np.log2(de[wt_cols + ko_cols].astype(float) + 1.0)
    mat.columns = [c.replace(".RNA_normCounts", "") for c in mat.columns]
    keep = mat.mean(axis=1) > 1
    return de[keep], mat[keep]


def main() -> None:
    os.makedirs(L.RES, exist_ok=True)
    sets = L.load_gene_sets("human")

    de, mat = load("GSE289287_DESeq2-Trop2KO_tumors_vs_WT.tsv.gz", WT_COLS, KO_COLS)
    wt = [c.replace(".RNA_normCounts", "") for c in WT_COLS]
    ko = [c.replace(".RNA_normCounts", "") for c in KO_COLS]
    print(f"{GSE} Trop2KO xenografts: {mat.shape[0]} protein-coding genes, "
          f"KO n={len(ko)} vs WT n={len(wt)}")

    # ---- KO validity -------------------------------------------------------
    v = mat.loc["TACSTD2"]
    t, p = stats.ttest_ind(v[ko], v[wt], equal_var=False)
    validity = {"gene": "TACSTD2", "mean_log2_KO": float(v[ko].mean()),
                "mean_log2_WT": float(v[wt].mean()),
                "log2FC_KO_vs_WT_myttest": float(v[ko].mean() - v[wt].mean()),
                "log2FC_submitter_DESeq2": float(de.loc["TACSTD2", "log2FoldChange"]),
                "padj_submitter_DESeq2": float(de.loc["TACSTD2", "padj"]),
                "pct_of_WT": float(100 * (2 ** v[ko].mean()) / (2 ** v[wt].mean())),
                "welch_t": float(t), "welch_p": float(p),
                "n_KO": len(ko), "n_WT": len(wt)}
    print("  KO validity:", {k: (round(x, 4) if isinstance(x, float) else x)
                             for k, x in validity.items()})
    pd.DataFrame([validity]).to_csv(
        os.path.join(L.RES, f"{GSE}_ko_validity.tsv"), sep="\t", index=False)

    # ---- gene level: submitter DESeq2 stats --------------------------------
    gs = de[["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]].copy()
    gs.columns = ["baseMean", "log2FC_KO_vs_WT", "lfcSE", "deseq2_stat",
                  "p_deseq2", "q_deseq2_BH"]
    gs["n_KO"], gs["n_WT"] = len(ko), len(wt)
    L.write_tsv(gs.sort_values("p_deseq2"), f"{GSE}_gene_stats_Trop2KO_xenograft.tsv")

    # ---- set level ---------------------------------------------------------
    gs_for_sets = gs.rename(columns={"log2FC_KO_vs_WT": "log2FC",
                                     "deseq2_stat": "t_stat"})
    sstats = L.set_level_tests(mat, gs_for_sets, sets, ko, wt, "TACSTD2_KO", "WT")
    L.write_tsv(sstats, f"{GSE}_set_stats_Trop2KO_xenograft.tsv", index=False)

    focus = ["TACSTD2", "CLDN4", "CLDN3", "CLDN7", "CLDN1", "EPCAM", "OCLN",
             "TJP1", "TJP2", "DSG2", "DSP", "JUP", "CDH1", "VIM", "CD274",
             "B2M", "HLA-A", "HLA-B", "TAP1", "PSMB9", "NLRC5", "STAT1",
             "STAT2", "IRF7", "IRF9", "ISG15", "IFIT1", "IFIT3", "MX1",
             "OAS1", "CXCL9", "CXCL10", "CXCL11", "IDO1", "GBP1", "GBP2",
             "NFKBIA", "CXCL8", "IL6"]
    L.write_tsv(gs.reindex([g for g in focus if g in gs.index]),
                f"{GSE}_focus_genes_Trop2KO_xenograft.tsv")

    # ---- DSG2 KO cell-line comparator (junction perturbation, not TACSTD2) --
    d2 = pd.read_csv(
        os.path.join(L.DATA, GSE, "GSE289287_DESeq2-DSG2KO_cells_vs_DSG2WT.tsv.gz"),
        sep="\t")
    d2 = d2[(d2["Feature_name"].notna()) & (d2["biotype"] == "protein_coding")]
    d2 = (d2.sort_values("baseMean", ascending=False)
            .drop_duplicates("Feature_name").set_index("Feature_name"))
    comp = L.competitive_from_ranking(d2["stat"].astype(float), sets)
    comp.insert(0, "contrast", "DSG2_KO_vs_WT_cells_COMPARATOR_not_TACSTD2")
    L.write_tsv(comp, f"{GSE}_set_stats_DSG2KO_cells_comparator.tsv", index=False)
    dsg_focus = d2.reindex([g for g in focus if g in d2.index])[
        ["baseMean", "log2FoldChange", "stat", "pvalue", "padj"]]
    L.write_tsv(dsg_focus, f"{GSE}_focus_genes_DSG2KO_cells_comparator.tsv")

    make_figures(mat, gs, sstats, sets, ko, wt, d2)
    print("done")


def make_figures(mat, gs, sstats, sets, ko, wt, d2) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(17, 10))

    ax = axes[0, 0]
    genes = ["TACSTD2", "CLDN4", "CLDN3", "CLDN7", "CLDN1", "EPCAM"]
    genes = [g for g in genes if g in mat.index]
    for i, g in enumerate(genes):
        ax.scatter(np.full(len(wt), i) - .16 + np.random.uniform(-.05, .05, len(wt)),
                   mat.loc[g, wt], color="#4477aa", s=38,
                   label="WT" if i == 0 else None)
        ax.scatter(np.full(len(ko), i) + .16 + np.random.uniform(-.05, .05, len(ko)),
                   mat.loc[g, ko], color="#cc3311", s=38,
                   label="TACSTD2 KO" if i == 0 else None)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("log2(DESeq2 normalised count + 1)")
    ax.set_title("GSE289287 T-47D xenografts:\ntarget + claudins (n=4 KO vs 3 WT)")
    ax.legend(frameon=False, fontsize=9)

    ax = axes[0, 1]
    q = gs["q_deseq2_BH"].fillna(1.0)
    ax.scatter(gs["log2FC_KO_vs_WT"],
               -np.log10(gs["p_deseq2"].clip(lower=1e-300)), s=5,
               c=np.where(q < 0.05, "#cc3311", "#bbbbbb"), alpha=.55)
    for g in ["TACSTD2", "CLDN4", "CLDN3", "CLDN7", "CDH1", "CXCL10", "B2M", "CD274"]:
        if g in gs.index and pd.notna(gs.loc[g, "p_deseq2"]):
            ax.annotate(g, (gs.loc[g, "log2FC_KO_vs_WT"],
                            -np.log10(max(gs.loc[g, "p_deseq2"], 1e-300))), fontsize=8)
    ax.axvline(0, lw=.5, c="k")
    ax.set_xlabel("log2FC (TACSTD2 KO / WT)")
    ax.set_ylabel("-log10 DESeq2 p")
    ax.set_title(f"Submitter DESeq2\n{int((q<0.05).sum())} of {len(gs)} "
                 "protein-coding genes at q<0.05")

    ax = axes[0, 2]
    order = ["TARGET_GENES", "CLAUDIN_FAMILY", "TIGHT_JUNCTION_CORE",
             "DESMOSOME_ADHERENS", "EPITHELIAL_IDENTITY", "EMT_CORE",
             "IFN_ALPHA_TYPE1", "IFN_GAMMA_TYPE2", "ANTIGEN_PRESENTATION_MHC1",
             "CHEMOKINE_T_RECRUITMENT", "ICI_CHECKPOINT_GENES",
             "NFKB_INFLAMMATORY"]
    s = sstats.set_index("set_name").reindex(order)
    y = np.arange(len(order))
    ax.barh(y, s["sample_delta"].fillna(0),
            color=["#cc3311" if v > 0 else "#4477aa"
                   for v in s["sample_delta"].fillna(0)])
    ax.set_yticks(y)
    ax.set_yticklabels([n.replace("_", " ").title() for n in order], fontsize=8)
    ax.invert_yaxis()
    ax.axvline(0, c="k", lw=.6)
    ax.set_xlabel("Δ mean-z set score (KO − WT)")
    ax.set_title("Set score shift; * = sample-level p<0.05")
    for i, (dd, p) in enumerate(zip(s["sample_delta"], s["sample_p"])):
        if pd.notna(p) and p < 0.05:
            ax.text(dd + (.03 if dd > 0 else -.03), i, "*", va="center",
                    ha="left" if dd > 0 else "right", fontsize=13)

    ax = axes[1, 0]
    show = ["CLAUDIN_FAMILY", "TIGHT_JUNCTION_CORE", "IFN_ALPHA_TYPE1",
            "IFN_GAMMA_TYPE2", "ANTIGEN_PRESENTATION_MHC1", "EPITHELIAL_IDENTITY"]
    for i, name in enumerate(show):
        sc = L.set_score_matrix(mat, sets[name])
        if sc is None:
            continue
        ax.scatter(np.full(len(wt), i) - .16, sc[wt], color="#4477aa", s=34,
                   label="WT" if i == 0 else None)
        ax.scatter(np.full(len(ko), i) + .16, sc[ko], color="#cc3311", s=34,
                   label="TACSTD2 KO" if i == 0 else None)
    ax.set_xticks(range(len(show)))
    ax.set_xticklabels([n.replace("_", "\n") for n in show], fontsize=7)
    ax.axhline(0, lw=.5, c="k")
    ax.set_ylabel("per-sample mean-z score")
    ax.set_title("Per-sample set scores")
    ax.legend(frameon=False, fontsize=9)

    ax = axes[1, 1]
    cl = [g for g in sorted(sets["CLAUDIN_FAMILY"]) if g in gs.index]
    fc = gs.loc[cl, "log2FC_KO_vs_WT"]
    ax.barh(np.arange(len(cl)), fc,
            color=["#cc3311" if v > 0 else "#4477aa" for v in fc])
    ax.set_yticks(np.arange(len(cl)))
    ax.set_yticklabels(cl, fontsize=8)
    ax.invert_yaxis()
    ax.axvline(0, c="k", lw=.6)
    ax.set_xlabel("log2FC (TACSTD2 KO / WT)")
    ax.set_title("Claudin family, TACSTD2 KO xenografts")
    for i, g in enumerate(cl):
        if pd.notna(gs.loc[g, "q_deseq2_BH"]) and gs.loc[g, "q_deseq2_BH"] < 0.05:
            ax.text(fc[g] + (.1 if fc[g] > 0 else -.1), i, "*", va="center",
                    ha="left" if fc[g] > 0 else "right", fontsize=12)

    ax = axes[1, 2]
    common = [g for g in gs.index if g in d2.index]
    x = gs.loc[common, "log2FC_KO_vs_WT"].astype(float)
    y2 = d2.loc[common, "log2FoldChange"].astype(float)
    ok = x.notna() & y2.notna()
    r, pr = stats.pearsonr(x[ok], y2[ok])
    rho, prho = stats.spearmanr(x[ok], y2[ok])
    ax.scatter(x[ok], y2[ok], s=4, alpha=.25, c="#888888")
    for g in ["TACSTD2", "DSG2", "CLDN4", "CLDN3", "CLDN7"]:
        if g in common:
            ax.scatter(x[g], y2[g], s=60, c="#cc3311", zorder=5)
            ax.annotate(g, (x[g], y2[g]), fontsize=9)
    ax.axhline(0, lw=.5, c="k")
    ax.axvline(0, lw=.5, c="k")
    ax.set_xlabel("log2FC TACSTD2 KO vs WT (xenograft)")
    ax.set_ylabel("log2FC DSG2 KO vs WT (cells)")
    ax.set_title(f"TACSTD2 KO vs DSG2 KO comparator\nPearson r={r:.3f} "
                 f"(p={L.fmt_p(pr)}), Spearman rho={rho:.3f}, n={int(ok.sum())} genes")

    fig.suptitle("GSE289287 — human T-47D TACSTD2 (Trop-2) knockout xenografts, "
                 "NRG immunodeficient host (tumour-intrinsic arm only)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    out = os.path.join(L.RES, f"fig_{GSE}_tacstd2ko_T47D.png")
    fig.savefig(out, dpi=160)
    print(f"  wrote {os.path.basename(out)}")
    plt.close(fig)

    pd.DataFrame([{"comparison": "TACSTD2_KO_xenograft_vs_DSG2_KO_cells",
                   "pearson_r": r, "pearson_p": pr, "spearman_rho": rho,
                   "spearman_p": prho, "n_genes": int(ok.sum())}]).to_csv(
        os.path.join(L.RES, f"{GSE}_trop2_vs_dsg2_concordance.tsv"),
        sep="\t", index=False)


if __name__ == "__main__":
    np.random.seed(0)
    main()
