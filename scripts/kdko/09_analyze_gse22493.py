#!/usr/bin/env python3
"""GSE22493 -- human SKOV-3-IP-Luc ovarian cancer, CLDN4 siRNA-lentivirus silencing.

Design (from the GEO series matrix):
    GSM558700-702, three two-colour arrays (platform GPL10555, BWH Human Release
    3.0 oligos). Channel 1 (Cy3) = "CLDN4 overexpression (control)";
    channel 2 (Cy5) = "CLDN4 knockdown". The deposited values are normalised
    sample/control log ratios, one column per array, with dye-swap arrays used by
    the submitters to control hybridisation bias.

Two design facts that limit what this series can support, both handled below:
  1. The reference arm is a CLDN4-OVEREXPRESSING line, not the parental line, so
     the contrast is CLDN4-low vs CLDN4-high, a wider swing than a plain
     knockdown but not attributable to knockdown alone.
  2. TACSTD2 and EPCAM are absent from GPL10555, so this series cannot test the
     reciprocal CLDN4 -> TACSTD2 question at all.

Sign convention verified empirically: the CLDN4 probe is negative (-1.74, -0.71,
one NaN), consistent with log2(knockdown / control).
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

GSE = "GSE22493"
ARRAYS = ["GSM558700", "GSM558701", "GSM558702"]


def load() -> pd.DataFrame:
    expr, _ = L.read_series_matrix(
        os.path.join(L.DATA, GSE, f"{GSE}_series_matrix.txt.gz"))
    plat = L.read_platform_table(os.path.join(L.DATA, GSE, f"{GSE}_family.soft.gz"))
    sym = plat["ORF"].replace("", np.nan)
    expr.index = expr.index.astype(str)
    sym.index = sym.index.astype(str)
    mat = L.collapse_to_symbol(expr[ARRAYS], sym, how="mean")
    return mat.dropna(thresh=2)  # need >=2 of 3 arrays for a one-sample test


def main() -> None:
    os.makedirs(L.RES, exist_ok=True)
    mat = load()
    sets = L.load_gene_sets("human")
    print(f"{GSE}: {mat.shape[0]} genes with >=2 of 3 arrays "
          "(log2 knockdown/control ratios)")

    # ---- array-to-array agreement: this series' key quality metric ----------
    agree = []
    for i in range(3):
        for j in range(i + 1, 3):
            a, b = mat[ARRAYS[i]], mat[ARRAYS[j]]
            ok = a.notna() & b.notna()
            r, p = stats.pearsonr(a[ok], b[ok])
            rho, prho = stats.spearmanr(a[ok], b[ok])
            agree.append({"array_1": ARRAYS[i], "array_2": ARRAYS[j],
                          "pearson_r": r, "pearson_p": p, "spearman_rho": rho,
                          "spearman_p": prho, "n_genes": int(ok.sum())})
    ag = pd.DataFrame(agree)
    L.write_tsv(ag, f"{GSE}_array_agreement.tsv", index=False)
    print("  inter-array log-ratio agreement (Pearson r): "
          + ", ".join(f"{r:.3f}" for r in ag["pearson_r"]))

    # ---- per-gene one-sample test of the log ratio against zero ------------
    n_obs = mat.notna().sum(axis=1)
    mean_lr = mat.mean(axis=1, skipna=True)
    t = np.full(len(mat), np.nan)
    p = np.full(len(mat), np.nan)
    for i, g in enumerate(mat.index):
        v = mat.loc[g].dropna().to_numpy(dtype=float)
        if len(v) >= 3 and np.ptp(v) > 0:
            tt = stats.ttest_1samp(v, 0.0)
            t[i], p[i] = tt.statistic, tt.pvalue
    gs = pd.DataFrame({
        "contrast": "CLDN4_siRNA_KD_vs_CLDN4_overexpressing_control",
        "mean_log2ratio_KD_over_control": mean_lr,
        "n_arrays_with_value": n_obs,
        "t_stat_vs_zero": t, "p_onesample_t": p,
        "q_onesample_BH": L.bh_fdr(p),
        "n_biological_replicates": 3,
        "caveat": "reference arm is CLDN4-overexpressing, not parental",
    }, index=mat.index)
    for k, a in enumerate(ARRAYS):
        gs[f"log2ratio_{a}"] = mat[a]
    L.write_tsv(gs.sort_values("p_onesample_t"), f"{GSE}_gene_stats_CLDN4KD.tsv")

    val = {"gene": "CLDN4",
           "mean_log2ratio_KD_over_control": float(mean_lr.get("CLDN4", np.nan)),
           "n_arrays_with_value": int(n_obs.get("CLDN4", 0)),
           "per_array": ";".join(
               f"{a}={mat.loc['CLDN4', a]:.3f}" if pd.notna(mat.loc["CLDN4", a])
               else f"{a}=NA" for a in ARRAYS),
           "p_onesample_t": float(gs.loc["CLDN4", "p_onesample_t"])
           if "CLDN4" in gs.index else np.nan,
           "TACSTD2_on_platform": "TACSTD2" in mat.index,
           "EPCAM_on_platform": "EPCAM" in mat.index}
    L.write_tsv(pd.DataFrame([val]), f"{GSE}_kd_validity.tsv", index=False)
    print(f"  CLDN4 log2 ratio (KD/control): {val['per_array']}  "
          f"mean={val['mean_log2ratio_KD_over_control']:.3f}")
    print(f"  TACSTD2 on platform: {val['TACSTD2_on_platform']}  "
          f"(so the reciprocal read-out is NOT testable here)")

    # ---- set level ---------------------------------------------------------
    # Only gene-level tests are meaningful here: with three arrays whose log
    # ratios barely agree (see array_agreement), a sample-level set test on
    # n=3 ratios would be reporting noise.
    comp = L.competitive_from_ranking(mean_lr, sets)
    comp.insert(0, "contrast", "CLDN4_KD_vs_CLDN4_OE_control")
    # self-contained: is the set's mean log ratio different from 0?
    sc = []
    for name, genes in sets.items():
        present = [g for g in genes if g in mean_lr.index]
        if len(present) < 3:
            sc.append({"set_name": name, "selfcontained_wilcoxon_p": np.nan,
                       "median_log2ratio_in_set": np.nan})
            continue
        v = mean_lr.loc[present].dropna()
        sc.append({"set_name": name,
                   "selfcontained_wilcoxon_p": float(
                       stats.wilcoxon(v, alternative="two-sided").pvalue),
                   "median_log2ratio_in_set": float(v.median())})
    comp = comp.merge(pd.DataFrame(sc), on="set_name", how="left")
    comp["selfcontained_q_BH"] = L.bh_fdr(
        comp["selfcontained_wilcoxon_p"].to_numpy())
    comp["inference_unit"] = "genes; 3 arrays only, see array_agreement file"
    L.write_tsv(comp, f"{GSE}_set_stats_CLDN4KD.tsv", index=False)

    focus = ["CLDN4", "CLDN1", "CLDN3", "CLDN7", "OCLN", "TJP1", "TJP2", "CDH1",
             "VIM", "DSG2", "B2M", "HLA-A", "HLA-B", "TAP1", "PSMB9", "STAT1",
             "STAT2", "IRF7", "IFIT1", "MX1", "OAS1", "CXCL9", "CXCL10",
             "IDO1", "GBP1", "NFKBIA", "IL6"]
    L.write_tsv(gs.reindex([g for g in focus if g in gs.index]),
                f"{GSE}_focus_genes_CLDN4KD.tsv")

    make_figures(mat, gs, comp, sets, ag)
    print("done")


def make_figures(mat, gs, comp, sets, ag) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    ax = axes[0, 0]
    ax.scatter(mat[ARRAYS[0]], mat[ARRAYS[1]], s=5, alpha=.2, c="#999999")
    if "CLDN4" in mat.index:
        ax.scatter(mat.loc["CLDN4", ARRAYS[0]], mat.loc["CLDN4", ARRAYS[1]],
                   s=90, c="#cc3311", zorder=5)
        ax.annotate("CLDN4", (mat.loc["CLDN4", ARRAYS[0]],
                              mat.loc["CLDN4", ARRAYS[1]]), fontsize=10)
    ax.axhline(0, lw=.5, c="k")
    ax.axvline(0, lw=.5, c="k")
    ax.set_xlabel(f"log2 ratio {ARRAYS[0]}")
    ax.set_ylabel(f"log2 ratio {ARRAYS[1]}")
    r0 = ag.iloc[0]["pearson_r"]
    ax.set_title(f"GSE22493 replicate agreement\nPearson r={r0:.3f} "
                 f"(all pairs: "
                 f"{', '.join(f'{v:.2f}' for v in ag['pearson_r'])})")

    ax = axes[0, 1]
    genes = [g for g in ["CLDN4", "CLDN1", "CLDN3", "CLDN7", "OCLN", "TJP1",
                         "TJP2", "CDH1", "VIM", "DSG2"] if g in mat.index]
    for i, g in enumerate(genes):
        v = mat.loc[g, ARRAYS].astype(float)
        ax.scatter(np.full(3, i), v, c="#4477aa", s=40)
        ax.hlines(np.nanmean(v), i - .3, i + .3, colors="#cc3311", lw=2)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=45, ha="right", fontsize=8)
    ax.axhline(0, lw=.6, c="k")
    ax.set_ylabel("log2(CLDN4 KD / CLDN4-OE control)")
    ax.set_title("Junction genes, per array\n(red bar = mean of 3 arrays)")

    ax = axes[1, 0]
    order = ["CLAUDIN_FAMILY", "TIGHT_JUNCTION_CORE", "DESMOSOME_ADHERENS",
             "EPITHELIAL_IDENTITY", "EMT_CORE", "IFN_ALPHA_TYPE1",
             "IFN_GAMMA_TYPE2", "ANTIGEN_PRESENTATION_MHC1",
             "CHEMOKINE_T_RECRUITMENT", "NFKB_INFLAMMATORY"]
    s = comp.set_index("set_name").reindex(order)
    y = np.arange(len(order))
    ax.barh(y, s["median_stat_in_set"].fillna(0),
            color=["#cc3311" if v > 0 else "#4477aa"
                   for v in s["median_stat_in_set"].fillna(0)])
    ax.set_yticks(y)
    ax.set_yticklabels([n.replace("_", " ").title() for n in order], fontsize=8)
    ax.invert_yaxis()
    ax.axvline(0, c="k", lw=.6)
    ax.set_xlabel("median log2 ratio of set genes (KD / CLDN4-OE control)")
    ax.set_title("Set-level shift; * = competitive q<0.05")
    for i, (dd, q) in enumerate(zip(s["median_stat_in_set"], s["q_BH"])):
        if pd.notna(q) and q < 0.05:
            ax.text(dd + (.02 if dd > 0 else -.02), i, "*", va="center",
                    ha="left" if dd > 0 else "right", fontsize=13)

    ax = axes[1, 1]
    v = gs["mean_log2ratio_KD_over_control"].dropna()
    ax.hist(v, bins=80, color="#bbbbbb")
    if "CLDN4" in gs.index:
        ax.axvline(gs.loc["CLDN4", "mean_log2ratio_KD_over_control"],
                   c="#cc3311", lw=2, label="CLDN4")
    ax.axvline(0, c="k", lw=.8)
    ax.axvline(float(v.median()), c="#4477aa", lw=1.5, ls="--",
               label=f"median = {v.median():.2f}")
    ax.set_xlabel("mean log2 ratio (KD / CLDN4-OE control)")
    ax.set_ylabel("genes")
    ax.set_title(f"Global log-ratio distribution, n={len(v)} genes\n"
                 "off-centre median indicates residual normalisation bias")
    ax.legend(frameon=False, fontsize=9)

    fig.suptitle("GSE22493 — human SKOV-3 CLDN4 siRNA silencing vs a "
                 "CLDN4-overexpressing reference (3 two-colour arrays; "
                 "TACSTD2 not on the platform)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    out = os.path.join(L.RES, f"fig_{GSE}_cldn4_kd_SKOV3.png")
    fig.savefig(out, dpi=160)
    print(f"  wrote {os.path.basename(out)}")
    plt.close(fig)


if __name__ == "__main__":
    main()
