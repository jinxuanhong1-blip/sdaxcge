#!/usr/bin/env python3
"""GSE334497 -- mouse 4T1 Trop2 KO vs WT mammary tumours (BALB/c, immunocompetent).

Design (from the GEO SOFT record, verified in scripts/kdko/03_fetch_datasets.py):
    Trop2 KO : GSM9789678-82   library names KO162, KO164, KO165, KO172, RESUB-KO163R
    WT       : GSM9789683-87   library names RESUB-171R/170R/169R/168R, control170
    n = 5 vs 5, whole-tumour frozen sections, 3 weeks in vivo growth.

This is the single most informative public set for the mechanistic question,
because the host is immunocompetent, so immune-compartment transcripts are real.

Caveat handled explicitly below: library naming splits the samples into a
"RESUB" re-sequenced group (4/5 WT, 1/5 KO) and a non-RESUB group (1/5 WT,
4/5 KO). Prep group is therefore largely confounded with genotype and is
tested for, not ignored.
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

GSE = "GSE334497"
KO = ["KO162", "KO164", "KO165", "KO172", "RESUB-KO163R"]
WT = ["RESUB-171R", "RESUB-170R", "RESUB-169R", "RESUB-168R", "control170"]
GSM = {"KO162": "GSM9789678", "KO164": "GSM9789679", "KO165": "GSM9789680",
       "KO172": "GSM9789681", "RESUB-KO163R": "GSM9789682",
       "RESUB-171R": "GSM9789683", "RESUB-170R": "GSM9789684",
       "RESUB-169R": "GSM9789685", "RESUB-168R": "GSM9789686",
       "control170": "GSM9789687"}


def load() -> pd.DataFrame:
    raw = pd.read_csv(
        os.path.join(L.DATA, GSE, f"{GSE}_normalized_counts.csv.gz"), index_col=0)
    raw.index = raw.index.astype(str)
    missing = [c for c in KO + WT if c not in raw.columns]
    if missing:
        raise SystemExit(f"columns missing from matrix: {missing}")
    raw = raw[KO + WT]
    sym = L.map_ensembl_to_symbol(list(raw.index), "mouse",
                                  os.path.join(L.DATA, "ensmusg_symbol_cache.json"))
    raw = raw.loc[[i for i in raw.index if i in sym]]
    raw.index = [sym[i] for i in raw.index]
    # DESeq2-style normalised counts -> log2; collapse duplicate symbols by max mean
    lg = np.log2(raw + 1.0)
    lg = L.collapse_to_symbol(lg, pd.Series(lg.index, index=lg.index))
    # drop genes with essentially no signal in either group
    keep = (lg[KO].mean(axis=1) > 1) | (lg[WT].mean(axis=1) > 1)
    return lg[keep]


def main() -> None:
    os.makedirs(L.RES, exist_ok=True)
    mat = load()
    print(f"{GSE}: {mat.shape[0]} genes x {mat.shape[1]} samples "
          f"(KO n={len(KO)}, WT n={len(WT)})")

    gstats = L.per_gene_tests(mat, KO, WT, "KO", "WT")
    gstats.insert(0, "contrast", "Trop2_KO_vs_WT")
    L.write_tsv(gstats.sort_values("p_welch"), f"{GSE}_gene_stats_KOvsWT.tsv")

    sets = L.load_gene_sets("mouse")
    sstats = L.set_level_tests(mat, gstats, sets, KO, WT, "Trop2_KO", "WT")
    L.write_tsv(sstats, f"{GSE}_set_stats_KOvsWT.tsv", index=False)

    # ---- targeted read-outs ------------------------------------------------
    focus = ["Tacstd2", "Cldn4", "Cldn3", "Cldn7", "Cldn1", "Cldn6", "Epcam",
             "Cd274", "Cd8a", "Gzmb", "Prf1", "Ifng", "Cxcl9", "Cxcl10",
             "Stat1", "B2m", "Tap1", "H2-K1", "Psmb9", "Nlrc5", "Isg15",
             "Ifit1", "Ifit3", "Mx1", "Oasl2", "Irf7", "Gbp2", "Ptprc",
             "Ccl5", "Nkg7", "Klrd1", "Ido1", "Lag3", "Pdcd1", "Ctla4",
             "Ocln", "Tjp1", "Dsg2", "Cdh1", "Vim"]
    tf = gstats.reindex([g for g in focus if g in gstats.index])
    L.write_tsv(tf, f"{GSE}_focus_genes_KOvsWT.tsv")

    # ---- KO validity check -------------------------------------------------
    validity = {}
    if "Tacstd2" in mat.index:
        v = mat.loc["Tacstd2"]
        t, p = stats.ttest_ind(v[KO], v[WT], equal_var=False)
        validity = {"gene": "Tacstd2",
                    "mean_log2_KO": float(v[KO].mean()),
                    "mean_log2_WT": float(v[WT].mean()),
                    "log2FC_KO_vs_WT": float(v[KO].mean() - v[WT].mean()),
                    "pct_of_WT": float(100 * (2 ** v[KO].mean()) / (2 ** v[WT].mean())),
                    "welch_t": float(t), "welch_p": float(p),
                    "n_KO": len(KO), "n_WT": len(WT)}
        print("  KO validity:", {k: (round(x, 4) if isinstance(x, float) else x)
                                 for k, x in validity.items()})
    pd.DataFrame([validity]).to_csv(
        os.path.join(L.RES, f"{GSE}_ko_validity.tsv"), sep="\t", index=False)

    # ---- prep-group (RESUB) confound diagnostics ---------------------------
    resub = [s for s in mat.columns if s.startswith("RESUB")]
    nonresub = [s for s in mat.columns if not s.startswith("RESUB")]
    top = mat.loc[mat.var(axis=1).sort_values(ascending=False).index[:2000]]
    z = top.sub(top.mean(axis=1), axis=0).div(
        top.std(axis=1, ddof=1).replace(0, np.nan), axis=0).dropna()
    u, s, vt = np.linalg.svd(z.to_numpy(), full_matrices=False)
    pcs = pd.DataFrame(vt[:3].T, index=z.columns, columns=["PC1", "PC2", "PC3"])
    var_expl = (s ** 2 / (s ** 2).sum())[:3]
    diag = []
    for pc in ["PC1", "PC2", "PC3"]:
        tg, pg = stats.ttest_ind(pcs.loc[KO, pc], pcs.loc[WT, pc], equal_var=False)
        tb, pb = stats.ttest_ind(pcs.loc[resub, pc], pcs.loc[nonresub, pc],
                                 equal_var=False)
        diag.append({"component": pc,
                     "var_explained": float(var_expl[["PC1", "PC2", "PC3"].index(pc)]),
                     "genotype_t": float(tg), "genotype_p": float(pg),
                     "prepgroup_RESUB_t": float(tb), "prepgroup_RESUB_p": float(pb),
                     "n_KO": len(KO), "n_WT": len(WT),
                     "n_RESUB": len(resub), "n_nonRESUB": len(nonresub)})
    diag_df = pd.DataFrame(diag)
    L.write_tsv(diag_df, f"{GSE}_pca_confound_diagnostics.tsv", index=False)
    pcs.assign(genotype=["KO" if s in KO else "WT" for s in pcs.index],
               prep_group=["RESUB" if s.startswith("RESUB") else "original"
                           for s in pcs.index],
               gsm=[GSM[s] for s in pcs.index]).to_csv(
        os.path.join(L.RES, f"{GSE}_pca_coords.tsv"), sep="\t")

    # ---- figures -----------------------------------------------------------
    make_figures(mat, gstats, sstats, sets, pcs, var_expl)
    print("done")


def make_figures(mat, gstats, sstats, sets, pcs, var_expl) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(17, 10))

    ax = axes[0, 0]
    for name, cols, c in [("WT", WT, "#4477aa"), ("Trop2 KO", KO, "#cc3311")]:
        for g_i, g in enumerate(["Tacstd2", "Cldn4", "Cldn3", "Cldn7"]):
            if g not in mat.index:
                continue
            x = np.full(len(cols), g_i) + (0.18 if name == "Trop2 KO" else -0.18)
            ax.scatter(x + np.random.uniform(-.05, .05, len(cols)),
                       mat.loc[g, cols], color=c, s=34, alpha=.85,
                       label=name if g_i == 0 else None)
    ax.set_xticks(range(4))
    ax.set_xticklabels(["Tacstd2", "Cldn4", "Cldn3", "Cldn7"])
    ax.set_ylabel("log2(normalised count + 1)")
    ax.set_title("GSE334497 4T1: target and claudins\n(n=5 KO vs 5 WT)")
    ax.legend(frameon=False, fontsize=9)

    ax = axes[0, 1]
    q = gstats["q_welch_BH"].fillna(1.0)
    ax.scatter(gstats["log2FC"], -np.log10(gstats["p_welch"].clip(lower=1e-300)),
               s=5, c=np.where(q < 0.1, "#cc3311", "#bbbbbb"), alpha=.6)
    for g in ["Tacstd2", "Cldn4", "Cxcl9", "Gzmb", "Cd8a", "Ifng", "Cd274", "B2m"]:
        if g in gstats.index and not np.isnan(gstats.loc[g, "p_welch"]):
            ax.annotate(g, (gstats.loc[g, "log2FC"],
                            -np.log10(max(gstats.loc[g, "p_welch"], 1e-300))),
                        fontsize=8)
    ax.axvline(0, lw=.5, c="k")
    ax.set_xlabel("log2FC (Trop2 KO / WT)")
    ax.set_ylabel("-log10 Welch p")
    ax.set_title(f"Per-gene DE\n{int((q<0.1).sum())} genes at BH q<0.1 "
                 f"of {len(gstats)}")

    ax = axes[0, 2]
    order = ["TARGET_GENES", "CLAUDIN_FAMILY", "TIGHT_JUNCTION_CORE",
             "IFN_ALPHA_TYPE1", "IFN_GAMMA_TYPE2", "ANTIGEN_PRESENTATION_MHC1",
             "T_CELL_CYTOTOXICITY", "T_CELL_INFILTRATION_BROAD",
             "CHEMOKINE_T_RECRUITMENT", "ICI_CHECKPOINT_GENES",
             "TIS_18GENE_AYERS", "IFNG_6GENE_AYERS", "NFKB_INFLAMMATORY",
             "EMT_CORE", "EPITHELIAL_IDENTITY", "DESMOSOME_ADHERENS"]
    s = sstats.set_index("set_name").reindex(order)
    y = np.arange(len(s))
    cols = ["#cc3311" if v > 0 else "#4477aa" for v in s["sample_delta"].fillna(0)]
    ax.barh(y, s["sample_delta"].fillna(0), color=cols)
    ax.set_yticks(y)
    ax.set_yticklabels([n.replace("_", " ").title() for n in order], fontsize=8)
    ax.invert_yaxis()
    ax.axvline(0, c="k", lw=.6)
    ax.set_xlabel("Δ mean-z set score (KO − WT)")
    ax.set_title("Set score shift; * = sample-level p<0.05")
    for i, (d, p) in enumerate(zip(s["sample_delta"], s["sample_p"])):
        if pd.notna(p) and p < 0.05:
            ax.text(d + (0.03 if d > 0 else -0.03), i, "*", va="center",
                    ha="left" if d > 0 else "right", fontsize=13)

    ax = axes[1, 0]
    show = ["IFN_ALPHA_TYPE1", "IFN_GAMMA_TYPE2", "T_CELL_CYTOTOXICITY",
            "TIS_18GENE_AYERS", "CLAUDIN_FAMILY", "TIGHT_JUNCTION_CORE"]
    for i, name in enumerate(show):
        sc = L.set_score_matrix(mat, sets[name])
        if sc is None:
            continue
        ax.scatter(np.full(len(WT), i) - .16 + np.random.uniform(-.05, .05, len(WT)),
                   sc[WT], color="#4477aa", s=32,
                   label="WT" if i == 0 else None)
        ax.scatter(np.full(len(KO), i) + .16 + np.random.uniform(-.05, .05, len(KO)),
                   sc[KO], color="#cc3311", s=32,
                   label="Trop2 KO" if i == 0 else None)
    ax.set_xticks(range(len(show)))
    ax.set_xticklabels([n.replace("_", "\n") for n in show], fontsize=7)
    ax.axhline(0, lw=.5, c="k")
    ax.set_ylabel("per-sample mean-z score")
    ax.set_title("Per-sample set scores (real replicates)")
    ax.legend(frameon=False, fontsize=9)

    ax = axes[1, 1]
    for gt, c in [("WT", "#4477aa"), ("KO", "#cc3311")]:
        idx = WT if gt == "WT" else KO
        for s_ in idx:
            m = "o" if not s_.startswith("RESUB") else "^"
            ax.scatter(pcs.loc[s_, "PC1"], pcs.loc[s_, "PC2"], c=c, marker=m, s=80)
            ax.annotate(s_.replace("RESUB-", "R:"),
                        (pcs.loc[s_, "PC1"], pcs.loc[s_, "PC2"]), fontsize=6)
    ax.set_xlabel(f"PC1 ({100*var_expl[0]:.0f}%)")
    ax.set_ylabel(f"PC2 ({100*var_expl[1]:.0f}%)")
    ax.set_title("PCA, top-2000 variable genes\ncircle=original library, "
                 "triangle=RESUB (red=KO)")

    ax = axes[1, 2]
    imm = ["Cd8a", "Cd3e", "Gzmb", "Prf1", "Ifng", "Cxcl9", "Cxcl10", "Ccl5",
           "Nkg7", "B2m", "Tap1", "Psmb9", "Stat1", "Cd274", "Ptprc", "Isg15",
           "Ifit1", "Ifit3", "Mx1", "Irf7"]
    imm = [g for g in imm if g in gstats.index]
    fc = gstats.loc[imm, "log2FC"]
    ax.barh(np.arange(len(imm)), fc,
            color=["#cc3311" if v > 0 else "#4477aa" for v in fc])
    ax.set_yticks(np.arange(len(imm)))
    ax.set_yticklabels(imm, fontsize=8)
    ax.invert_yaxis()
    ax.axvline(0, c="k", lw=.6)
    ax.set_xlabel("log2FC (Trop2 KO / WT)")
    ax.set_title("Immune / IFN / APM genes")
    for i, g in enumerate(imm):
        p = gstats.loc[g, "p_welch"]
        if pd.notna(p) and p < 0.05:
            ax.text(fc[g] + (.05 if fc[g] > 0 else -.05), i, "*", va="center",
                    ha="left" if fc[g] > 0 else "right", fontsize=12)

    fig.suptitle("GSE334497 — 4T1 Trop2 knockout vs wild-type mammary tumours "
                 "(BALB/c, immunocompetent), n=5 vs 5", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    out = os.path.join(L.RES, f"fig_{GSE}_trop2ko_4T1.png")
    fig.savefig(out, dpi=160)
    print(f"  wrote {os.path.basename(out)}")
    plt.close(fig)


if __name__ == "__main__":
    np.random.seed(0)
    main()
