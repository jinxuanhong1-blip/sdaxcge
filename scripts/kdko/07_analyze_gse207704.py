#!/usr/bin/env python3
"""GSE207704 -- human CLDN4 CRISPR knockout in two breast cancer lines.

Design (from the GEO SOFT record):
    T47D WT        GSM6310640-41   T47D:CLDN4-/-   GSM6310642-43
    MCF7 WT        GSM6310644-45   MCF7:CLDN4-/-   GSM6310646-47
    416 bp of the CLDN4 CDS including the start codon removed with Cas9 + gRNA pair.

CRITICAL DATA LIMIT: the series has 8 samples (2 replicates per group), but the
only open processed file, GSE207704_CLDN4_RNAseq.txt.gz, publishes one FPKM
column per GROUP (4 columns), i.e. the replicates are already collapsed. There
are no per-GSM supplementary files and no counts. Therefore NO sample-level
p-value can be computed from public data for this series, and none is reported.

What can legitimately be done, and is done here:
  * per-cell-line log2 fold change (KO vs WT) from the group-mean FPKMs;
  * cross-cell-line reproducibility of those fold changes (T47D vs MCF7);
  * gene-set inference at GENE level: competitive Mann-Whitney of the set's fold
    changes against all other genes, and a self-contained Wilcoxon signed-rank
    test of the set's fold changes against zero. n for these tests is a number
    of GENES, not a number of biological replicates -- labelled as such
    everywhere in the output.
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

GSE = "GSE207704"
COLS = {"MCF7_KO": "MCF7_CLDN4KO_FPKM (fpkm)", "MCF7_WT": "MCF7_WT_FPKM (fpkm)",
        "T47D_KO": "T47D_CLDN4KO_FPKM (fpkm)", "T47D_WT": "T47D_WT_FPKM (fpkm)"}
PSEUDO = 0.5  # FPKM pseudocount before log2


def load() -> pd.DataFrame:
    df = pd.read_csv(os.path.join(L.DATA, GSE, f"{GSE}_CLDN4_RNAseq.txt.gz"),
                     sep="\t", low_memory=False)
    df = df[df["gene_short_name"].notna()]
    for c in COLS.values():
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=list(COLS.values()))
    df["_mean"] = df[list(COLS.values())].mean(axis=1)
    df = (df.sort_values("_mean", ascending=False)
            .drop_duplicates("gene_short_name")
            .set_index("gene_short_name"))
    fpkm = df[list(COLS.values())].rename(
        columns={v: k for k, v in COLS.items()})
    # require expression in at least one condition of at least one cell line
    return fpkm[fpkm.max(axis=1) >= 1.0]


def selfcontained_and_competitive(fc: pd.Series, sets: dict[str, list[str]],
                                  tag: str) -> pd.DataFrame:
    rows = []
    fc = fc.dropna().astype(float)
    for name, genes in sets.items():
        present = [g for g in genes if g in fc.index]
        rec = {"contrast": tag, "set_name": name,
               "n_genes_measured": len(present),
               "inference_unit": "genes (NOT biological replicates)"}
        if len(present) < 3:
            rows.append({**rec, "median_log2FC_in_set": np.nan,
                         "median_log2FC_background": np.nan,
                         "competitive_mannwhitney_p": np.nan,
                         "selfcontained_wilcoxon_p": np.nan,
                         "frac_genes_up": np.nan})
            continue
        s_in, s_out = fc.loc[present], fc.drop(index=present)
        u = stats.mannwhitneyu(s_in, s_out, alternative="two-sided")
        try:
            w = stats.wilcoxon(s_in, alternative="two-sided").pvalue
        except ValueError:
            w = np.nan
        rows.append({**rec,
                     "median_log2FC_in_set": float(s_in.median()),
                     "median_log2FC_background": float(s_out.median()),
                     "competitive_mannwhitney_p": float(u.pvalue),
                     "selfcontained_wilcoxon_p": float(w),
                     "frac_genes_up": float((s_in > 0).mean())})
    df = pd.DataFrame(rows)
    df["competitive_q_BH"] = L.bh_fdr(df["competitive_mannwhitney_p"].to_numpy())
    df["selfcontained_q_BH"] = L.bh_fdr(df["selfcontained_wilcoxon_p"].to_numpy())
    return df


def main() -> None:
    os.makedirs(L.RES, exist_ok=True)
    fpkm = load()
    sets = L.load_gene_sets("human")
    lg = np.log2(fpkm + PSEUDO)
    fc = pd.DataFrame({
        "log2FC_T47D_KO_vs_WT": lg["T47D_KO"] - lg["T47D_WT"],
        "log2FC_MCF7_KO_vs_WT": lg["MCF7_KO"] - lg["MCF7_WT"]})
    fc["log2FC_mean_bothlines"] = fc.mean(axis=1)
    fc["same_direction_both_lines"] = (
        np.sign(fc["log2FC_T47D_KO_vs_WT"]) == np.sign(fc["log2FC_MCF7_KO_vs_WT"]))
    out = fpkm.join(fc)
    out.insert(0, "note", "group-mean FPKM only; no per-replicate data in GEO")
    print(f"{GSE}: {len(out)} genes, 4 group-mean FPKM columns "
          "(no per-replicate values published)")

    # ---- KO validity -------------------------------------------------------
    val = []
    for line in ["T47D", "MCF7"]:
        val.append({"cell_line": line, "gene": "CLDN4",
                    "FPKM_WT": float(fpkm.loc["CLDN4", f"{line}_WT"]),
                    "FPKM_KO": float(fpkm.loc["CLDN4", f"{line}_KO"]),
                    "log2FC_KO_vs_WT": float(fc.loc["CLDN4",
                                                    f"log2FC_{line}_KO_vs_WT"]),
                    "pct_of_WT": float(100 * fpkm.loc["CLDN4", f"{line}_KO"] /
                                       fpkm.loc["CLDN4", f"{line}_WT"]),
                    "n_replicates_in_GEO": 2,
                    "p_value": "not computable: replicates collapsed in the "
                               "only open processed file"})
    L.write_tsv(pd.DataFrame(val), f"{GSE}_ko_validity.tsv", index=False)
    for v in val:
        print(f"  KO validity {v['cell_line']}: CLDN4 FPKM {v['FPKM_WT']:.2f} -> "
              f"{v['FPKM_KO']:.2f} ({v['pct_of_WT']:.1f}% of WT, "
              f"log2FC {v['log2FC_KO_vs_WT']:+.2f})")

    L.write_tsv(out.sort_values("log2FC_mean_bothlines"),
                f"{GSE}_gene_log2FC_CLDN4KO.tsv")

    # ---- cross-cell-line reproducibility -----------------------------------
    a, b = fc["log2FC_T47D_KO_vs_WT"], fc["log2FC_MCF7_KO_vs_WT"]
    r, pr = stats.pearsonr(a, b)
    rho, prho = stats.spearmanr(a, b)
    frac = float(fc["same_direction_both_lines"].mean())
    binom = stats.binomtest(int(fc["same_direction_both_lines"].sum()),
                            len(fc), 0.5, alternative="greater").pvalue
    L.write_tsv(pd.DataFrame([{
        "comparison": "T47D_vs_MCF7_CLDN4KO_log2FC", "pearson_r": r,
        "pearson_p": pr, "spearman_rho": rho, "spearman_p": prho,
        "frac_same_direction": frac, "binomial_p_same_direction": binom,
        "n_genes": len(fc)}]), f"{GSE}_cellline_concordance.tsv", index=False)
    print(f"  T47D vs MCF7 log2FC: r={r:.3f} (p={L.fmt_p(pr)}), "
          f"{100*frac:.1f}% same direction (binomial p={L.fmt_p(binom)}), "
          f"n={len(fc)} genes")

    # ---- gene-set inference ------------------------------------------------
    allsets = []
    for tag, col in [("T47D_CLDN4KO_vs_WT", "log2FC_T47D_KO_vs_WT"),
                     ("MCF7_CLDN4KO_vs_WT", "log2FC_MCF7_KO_vs_WT"),
                     ("mean_of_both_lines", "log2FC_mean_bothlines")]:
        allsets.append(selfcontained_and_competitive(fc[col], sets, tag))
    setdf = pd.concat(allsets, ignore_index=True)
    L.write_tsv(setdf, f"{GSE}_set_stats_CLDN4KO.tsv", index=False)

    focus = ["CLDN4", "TACSTD2", "EPCAM", "CLDN1", "CLDN3", "CLDN7", "OCLN",
             "TJP1", "TJP2", "DSG2", "CDH1", "VIM", "CD274", "B2M", "HLA-A",
             "HLA-B", "TAP1", "PSMB9", "NLRC5", "STAT1", "STAT2", "IRF7",
             "IRF9", "ISG15", "IFIT1", "IFIT3", "MX1", "OAS1", "OASL",
             "CXCL9", "CXCL10", "IDO1", "GBP1", "GBP2", "NFKBIA", "CXCL8", "IL6"]
    L.write_tsv(out.reindex([g for g in focus if g in out.index])[
        ["T47D_WT", "T47D_KO", "MCF7_WT", "MCF7_KO", "log2FC_T47D_KO_vs_WT",
         "log2FC_MCF7_KO_vs_WT", "log2FC_mean_bothlines",
         "same_direction_both_lines"]], f"{GSE}_focus_genes_CLDN4KO.tsv")

    make_figures(fpkm, fc, setdf, sets, r)
    print("done")


def make_figures(fpkm, fc, setdf, sets, r) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(17, 10))

    ax = axes[0, 0]
    genes = ["CLDN4", "TACSTD2", "EPCAM", "CLDN3", "CLDN7", "CLDN1"]
    genes = [g for g in genes if g in fpkm.index]
    x = np.arange(len(genes))
    w = 0.2
    for i, (col, c, lab) in enumerate([
            ("T47D_WT", "#4477aa", "T47D WT"), ("T47D_KO", "#cc3311", "T47D CLDN4-/-"),
            ("MCF7_WT", "#88ccee", "MCF7 WT"), ("MCF7_KO", "#ee7733", "MCF7 CLDN4-/-")]):
        ax.bar(x + (i - 1.5) * w, np.log2(fpkm.loc[genes, col] + PSEUDO),
               width=w, color=c, label=lab)
    ax.set_xticks(x)
    ax.set_xticklabels(genes, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("log2(FPKM + 0.5)")
    ax.set_title("GSE207704: CLDN4 KO validity and TACSTD2 response\n"
                 "(group means; GEO publishes no per-replicate values)")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[0, 1]
    ax.scatter(fc["log2FC_T47D_KO_vs_WT"], fc["log2FC_MCF7_KO_vs_WT"],
               s=4, alpha=.2, c="#999999")
    for g in ["CLDN4", "TACSTD2", "EPCAM", "CLDN3", "CLDN7", "ISG15", "STAT1"]:
        if g in fc.index:
            ax.scatter(fc.loc[g, "log2FC_T47D_KO_vs_WT"],
                       fc.loc[g, "log2FC_MCF7_KO_vs_WT"], s=70, c="#cc3311", zorder=5)
            ax.annotate(g, (fc.loc[g, "log2FC_T47D_KO_vs_WT"],
                            fc.loc[g, "log2FC_MCF7_KO_vs_WT"]), fontsize=9)
    ax.axhline(0, lw=.5, c="k")
    ax.axvline(0, lw=.5, c="k")
    ax.set_xlabel("log2FC T47D CLDN4-/- vs WT")
    ax.set_ylabel("log2FC MCF7 CLDN4-/- vs WT")
    ax.set_title(f"Cross-cell-line reproducibility\nPearson r={r:.3f}, "
                 f"n={len(fc)} genes")

    ax = axes[0, 2]
    order = ["TARGET_GENES", "CLAUDIN_FAMILY", "TIGHT_JUNCTION_CORE",
             "DESMOSOME_ADHERENS", "EPITHELIAL_IDENTITY", "EMT_CORE",
             "IFN_ALPHA_TYPE1", "IFN_GAMMA_TYPE2", "ANTIGEN_PRESENTATION_MHC1",
             "CHEMOKINE_T_RECRUITMENT", "ICI_CHECKPOINT_GENES",
             "TIS_18GENE_AYERS", "NFKB_INFLAMMATORY"]
    s = setdf[setdf["contrast"] == "mean_of_both_lines"].set_index(
        "set_name").reindex(order)
    y = np.arange(len(order))
    ax.barh(y, s["median_log2FC_in_set"].fillna(0),
            color=["#cc3311" if v > 0 else "#4477aa"
                   for v in s["median_log2FC_in_set"].fillna(0)])
    ax.set_yticks(y)
    ax.set_yticklabels([n.replace("_", " ").title() for n in order], fontsize=8)
    ax.invert_yaxis()
    ax.axvline(0, c="k", lw=.6)
    ax.set_xlabel("median log2FC of set genes (mean of both lines)")
    ax.set_title("Set-level shift; * = competitive q<0.05\n(gene-level test, "
                 "not replicate-level)")
    for i, (dd, q) in enumerate(zip(s["median_log2FC_in_set"],
                                    s["competitive_q_BH"])):
        if pd.notna(q) and q < 0.05:
            ax.text(dd + (.02 if dd > 0 else -.02), i, "*", va="center",
                    ha="left" if dd > 0 else "right", fontsize=13)

    ax = axes[1, 0]
    for i, name in enumerate(["CLAUDIN_FAMILY", "TIGHT_JUNCTION_CORE",
                              "IFN_ALPHA_TYPE1", "IFN_GAMMA_TYPE2",
                              "ANTIGEN_PRESENTATION_MHC1", "EPITHELIAL_IDENTITY"]):
        present = [g for g in sets[name] if g in fc.index]
        ax.scatter(np.full(len(present), i) + np.random.uniform(-.14, .14, len(present)),
                   fc.loc[present, "log2FC_mean_bothlines"], s=18, alpha=.7,
                   c="#cc3311")
        ax.hlines(fc.loc[present, "log2FC_mean_bothlines"].median(),
                  i - .28, i + .28, colors="k", lw=2)
    ax.set_xticks(range(6))
    ax.set_xticklabels(["CLAUDIN\nFAMILY", "TIGHT\nJUNCTION", "IFN-I", "IFN-II",
                        "MHC-I\nAPM", "EPITHELIAL\nIDENTITY"], fontsize=7)
    ax.axhline(0, lw=.6, c="k")
    ax.set_ylabel("log2FC (CLDN4 KO / WT), mean of both lines")
    ax.set_title("Per-gene distribution within sets\n(black bar = median)")

    ax = axes[1, 1]
    genes = [g for g in ["TACSTD2", "EPCAM", "CLDN1", "CLDN3", "CLDN7", "OCLN",
                         "TJP1", "TJP2", "CDH1", "VIM", "DSG2"] if g in fc.index]
    yy = np.arange(len(genes))
    ax.barh(yy - .2, fc.loc[genes, "log2FC_T47D_KO_vs_WT"], height=.4,
            color="#cc3311", label="T47D")
    ax.barh(yy + .2, fc.loc[genes, "log2FC_MCF7_KO_vs_WT"], height=.4,
            color="#ee7733", label="MCF7")
    ax.set_yticks(yy)
    ax.set_yticklabels(genes, fontsize=8)
    ax.invert_yaxis()
    ax.axvline(0, c="k", lw=.6)
    ax.set_xlabel("log2FC (CLDN4 KO / WT)")
    ax.set_title("Junction genes and TACSTD2, per cell line")
    ax.legend(frameon=False, fontsize=9)

    ax = axes[1, 2]
    genes = [g for g in ["ISG15", "IFIT1", "IFIT3", "MX1", "OAS1", "OASL",
                         "STAT1", "STAT2", "IRF7", "IRF9", "B2M", "HLA-A",
                         "HLA-B", "TAP1", "PSMB9", "CD274", "CXCL10", "GBP1"]
             if g in fc.index]
    yy = np.arange(len(genes))
    ax.barh(yy - .2, fc.loc[genes, "log2FC_T47D_KO_vs_WT"], height=.4,
            color="#cc3311", label="T47D")
    ax.barh(yy + .2, fc.loc[genes, "log2FC_MCF7_KO_vs_WT"], height=.4,
            color="#ee7733", label="MCF7")
    ax.set_yticks(yy)
    ax.set_yticklabels(genes, fontsize=8)
    ax.invert_yaxis()
    ax.axvline(0, c="k", lw=.6)
    ax.set_xlabel("log2FC (CLDN4 KO / WT)")
    ax.set_title("IFN and antigen-presentation genes")
    ax.legend(frameon=False, fontsize=9)

    fig.suptitle("GSE207704 — human CLDN4 CRISPR knockout, T47D and MCF7 breast "
                 "cancer lines (group-mean FPKM; no replicate-level p-values "
                 "obtainable from GEO)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    out = os.path.join(L.RES, f"fig_{GSE}_cldn4ko_T47D_MCF7.png")
    fig.savefig(out, dpi=160)
    print(f"  wrote {os.path.basename(out)}")
    plt.close(fig)


if __name__ == "__main__":
    np.random.seed(0)
    main()
