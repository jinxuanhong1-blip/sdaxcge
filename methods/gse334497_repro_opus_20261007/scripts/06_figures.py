"""Figures for RESULTS.md (all read from tables/ written by scripts 01-05)."""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW, TAB, FIG = (os.path.join(HERE, d) for d in ("raw", "tables", "figures"))
os.makedirs(FIG, exist_ok=True)

COL = {"WT": "#1f77b4", "KO": "#d62728"}
MARK = {"original": "o", "RESUB": "^"}
PRIMARY_PANELS = ["APM_core", "IFN_core", "Immune_core", "CD8_T", "Effector", "NK", "Leukocyte",
                  "Checkpoint", "MHC_II", "Claudin_1_3_4_7", "Paper_TJ_trio_Cldn1_Cldn7_Ocln",
                  "TJ_claudin_program"]
FOCAL = ["Tacstd2", "Epcam", "Cldn1", "Cldn3", "Cldn4", "Cldn7", "Ocln", "Cgn",
         "Cd8a", "Cd8b1", "Gzmb", "Prf1", "Nkg7", "Ifng", "Cxcl9", "Cxcl10",
         "B2m", "H2-K1", "H2-D1", "Tap1", "Psmb9", "Stat1", "Irf1", "Isg15"]


def strip(ax, values, sheet, title, note=None):
    for i, g in enumerate(["WT", "KO"]):
        sel = sheet.genotype.values == g
        xs = i + np.linspace(-0.12, 0.12, sel.sum())
        for x, v, w in zip(xs, values[sel], sheet.library_wave.values[sel]):
            ax.scatter(x, v, c=COL[g], marker=MARK[w], s=26, edgecolor="k", linewidth=0.4, zorder=3)
        ax.hlines(np.mean(values[sel]), i - 0.25, i + 0.25, color="k", lw=1.2, zorder=4)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["WT", "KO"]); ax.set_xlim(-0.6, 1.6)
    ax.set_title(title + ("\n" + note if note else ""), fontsize=7)
    ax.tick_params(labelsize=7)


def fig_gsea():
    d = pd.read_csv(os.path.join(TAB, "gsea_paper_fig4A_reproduction.csv"))
    d["label"] = d.paper_label + " [" + d.msigdb_set.str.replace("_", " ").str.slice(0, 34) + "]"
    y = np.arange(len(d))[::-1]
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    ax.barh(y + 0.2, d.paper_NES_WT_vs_KO, 0.38, color="#999999", label="paper (Fig 4A)")
    ax.barh(y - 0.2, d.repro_NES_desktop_emulation, 0.38, color="#2ca02c", label="reproduced (GSEA emulation)")
    for yy, (_, r) in zip(y, d.iterrows()):
        ax.text(2.25 if r.repro_NES_desktop_emulation > 0 else -2.25, yy,
                f"p={r.repro_nominal_p_desktop_emulation:.3f}  FDR={r.repro_FDR_within_collection:.2f}",
                va="center", ha="left" if r.repro_NES_desktop_emulation > 0 else "right", fontsize=6.5)
    ax.set_yticks(y); ax.set_yticklabels(d.label, fontsize=6.8)
    ax.axvline(0, color="k", lw=0.6)
    ax.set_xlim(-4.6, 4.6)
    ax.set_xlabel("NES (positive = enriched in Trop2-WT)")
    ax.set_title("GSE334497 paper Fig 4A vs reproduction\n5 WT vs 5 KO, exact 252-labeling null; "
                 "text = reproduced nominal p, FDR within MSigDB sub-collection", fontsize=8)
    ax.legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig1_gsea_paper_vs_reproduction.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_focal(norm, sheet, ann, de):
    sym2id = {}
    for gid, g in ann.gene_name.items():
        if isinstance(g, str) and g not in sym2id:
            sym2id[g] = gid
    lg = np.log2(norm + 1)
    fig, axes = plt.subplots(4, 6, figsize=(11, 8.2))
    for ax, g in zip(axes.flat, FOCAL):
        gid = sym2id[g]
        r = de.loc[gid]
        strip(ax, lg.loc[gid].values, sheet, g,
              f"log2FC {r.log2FC_KO_vs_WT:+.2f}; p={r.pvalue:.2g}\npadj={r.padj_BH:.2g}")
    fig.suptitle("GSE334497 focal genes, log2(DESeq2-normalized + 1); DESeq2 Wald ~genotype, n=5 vs 5\n"
                 "circle = original library wave, triangle = RESUB wave; bar = group mean", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(FIG, "fig2_focal_genes_per_sample.png"), dpi=150)
    plt.close(fig)


def fig_panels(scores, tests, sheet):
    t = tests.set_index("panel")
    fig, axes = plt.subplots(2, 6, figsize=(11, 5.2))
    for ax, p in zip(axes.flat, PRIMARY_PANELS):
        r = t.loc[p]
        strip(ax, scores[p].values, sheet, p.replace("Paper_TJ_trio_Cldn1_Cldn7_Ocln", "Cldn1/Cldn7/Ocln"),
              f"diff {r.diff_KO_minus_WT:+.2f}; perm p={r.perm_exact_p:.3f}\nwave-adj p={r.wave_adjusted_p:.3f}")
    fig.suptitle("Panel scores = mean log2(norm+1) over panel genes; diff = KO - WT; exact permutation over "
                 "252 labelings\ncircle = original library wave, triangle = RESUB wave", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(os.path.join(FIG, "fig3_panel_scores.png"), dpi=150)
    plt.close(fig)


def fig_epithelial(scores, sheet, adj):
    a = adj.set_index(["panel", "covariate"])
    pairs = [("Claudin_1_3_4_7", "DIAG_epithelial_differentiation"),
             ("Paper_TJ_trio_Cldn1_Cldn7_Ocln", "DIAG_epithelial_differentiation"),
             ("Paper_TJ_trio_Cldn1_Cldn7_Ocln", "DIAG_epithelial_content_Krt8_Krt18")]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.7))
    for ax, (yv, xv) in zip(axes, pairs):
        for g in ["WT", "KO"]:
            for w in ["original", "RESUB"]:
                m = (sheet.genotype.values == g) & (sheet.library_wave.values == w)
                ax.scatter(scores.loc[m, xv], scores.loc[m, yv], c=COL[g], marker=MARK[w], s=34,
                           edgecolor="k", linewidth=0.4, label=f"{g} {w}")
        r = a.loc[(yv, xv)]
        ax.set_xlabel(xv.replace("DIAG_", "") + " score", fontsize=8)
        ax.set_ylabel(yv + " score", fontsize=8)
        ax.set_title(f"r={r.r_score_vs_covariate:.2f}; genotype coef adjusted for x: "
                     f"{r.genotype_coef_KO_vs_WT:+.2f} (p={r.genotype_p:.2f})", fontsize=7.5)
        ax.tick_params(labelsize=7)
    axes[0].legend(fontsize=6.5)
    fig.suptitle("Post hoc diagnostic: claudin scores vs epithelial-state scores (OLS score ~ genotype + x)",
                 fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(os.path.join(FIG, "fig4_claudin_vs_epithelial_state.png"), dpi=150)
    plt.close(fig)


def fig_forest():
    f = pd.read_csv(os.path.join(TAB, "focal_genes_all_models.csv")).set_index("gene_name")
    genes = [g for g in FOCAL if g in f.index][::-1]
    models = [("M1_all10_genotype", "M1 all 10", "k", 0.0),
              ("M2_all10_wave_plus_genotype", "M2 + wave", "#ff7f0e", 0.18),
              ("M3_drop_RESUB170R_outlier", "M3 drop RESUB-170R", "#9467bd", -0.18),
              ("M4_drop_control170", "M4 drop control170", "#8c564b", -0.32),
              ("M5_drop_skin_positive_WT", "M5 drop skin+ WT (2 vs 5)", "#17becf", 0.32)]
    fig, ax = plt.subplots(figsize=(6.4, 8.4))
    y = np.arange(len(genes))
    for key, lab, c, off in models:
        lfc = f.loc[genes, f"{key}_log2FC"].values
        se = f.loc[genes, f"{key}_lfcSE"].values
        ax.errorbar(lfc, y + off, xerr=1.96 * se, fmt="o", ms=3, color=c, lw=0.8, label=lab)
    ax.axvline(0, color="grey", lw=0.6)
    ax.set_yticks(y); ax.set_yticklabels(genes, fontsize=7.5)
    ax.set_xlabel("DESeq2 log2FC KO vs WT (Wald, 95% CI)")
    ax.set_title("Focal genes across sensitivity models", fontsize=9)
    ax.legend(fontsize=7, loc="lower right")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig5_focal_gene_forest_sensitivity.png"), dpi=150)
    plt.close(fig)


def fig_contamination(sheet):
    per = pd.read_csv(os.path.join(TAB, "qc_tissue_contamination_per_sample.csv"), index_col=0).loc[sheet.index]
    order = list(per.sort_values(["genotype", "skin_epidermis_score"], ascending=[False, False]).index)
    fig, axes = plt.subplots(1, 4, figsize=(13, 3.6), gridspec_kw={"width_ratios": [1.6, 1, 1, 1]})
    ax = axes[0]
    x = np.arange(len(order))
    ax.bar(x - 0.2, per.loc[order, "skin_epidermis_score"], 0.4, color=[COL[g] for g in per.loc[order, "genotype"]],
           edgecolor="k", lw=0.4, label="epidermis score")
    ax.bar(x + 0.2, per.loc[order, "skeletal_muscle_score"], 0.4, color="none",
           edgecolor=[COL[g] for g in per.loc[order, "genotype"]], hatch="///", lw=0.8, label="muscle score")
    ax.set_xticks(x); ax.set_xticklabels(order, rotation=60, ha="right", fontsize=6.5)
    ax.set_ylabel("mean log2(norm+1)", fontsize=8); ax.legend(fontsize=6.5)
    ax.set_title("Epidermis (solid) and skeletal-muscle (hatched) carry-over", fontsize=8)
    for ax, g in zip(axes[1:], ["Cldn1", "Cldn7", "Tacstd2"]):
        for gt in ["WT", "KO"]:
            for w in ["original", "RESUB"]:
                m = (per.genotype == gt) & (per.library_wave == w)
                ax.scatter(per.loc[m, "skin_epidermis_score"], np.log2(per.loc[m, f"norm_{g}"] + 1), c=COL[gt],
                           marker=MARK[w], s=34, edgecolor="k", linewidth=0.4)
        ax.set_xlabel("epidermis score", fontsize=8); ax.set_ylabel(f"log2({g} norm + 1)", fontsize=8)
        ax.set_title(g, fontsize=8); ax.tick_params(labelsize=7)
    fig.suptitle("Post hoc: 3/5 WT and 0/5 KO tumours contain epidermis (Krt1/Flg2/Krt77/Dsg1a/Dsc3); "
                 "blue = WT, red = KO", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(os.path.join(FIG, "fig6_tissue_contamination.png"), dpi=150)
    plt.close(fig)


def main():
    norm = pd.read_csv(os.path.join(RAW, "GSE334497_normalized_counts.csv.gz"), index_col=0)
    sheet = pd.read_csv(os.path.join(TAB, "sample_sheet.csv")).set_index("library_name").loc[norm.columns]
    ann = pd.read_csv(os.path.join(RAW, "gene_annotation_ensembl102.tsv.gz"), sep="\t", index_col=0)
    de = pd.read_csv(os.path.join(TAB, "de_deseq2_M1_all10_genotype.csv.gz")).set_index("ensembl_gene_id")
    scores = pd.read_csv(os.path.join(TAB, "panel_scores_per_sample_log2.csv"), index_col=0).loc[norm.columns]
    tests = pd.read_csv(os.path.join(TAB, "panel_score_tests.csv"))
    adj = pd.read_csv(os.path.join(TAB, "tj_scores_adjusted_for_epithelial_state.csv"))
    fig_gsea()
    fig_focal(norm, sheet, ann, de)
    fig_panels(scores, tests, sheet)
    fig_epithelial(scores, sheet, adj)
    fig_forest()
    fig_contamination(sheet)
    print(sorted(os.listdir(FIG)))


if __name__ == "__main__":
    main()
