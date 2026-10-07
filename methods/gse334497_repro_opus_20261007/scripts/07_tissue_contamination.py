"""Post hoc: epidermis / skeletal-muscle carry-over in the GSE334497 tumour sections.

Krt1, Krtdap and Lypd5 head the KO-vs-WT DESeq2 list (all lower in KO). These are
suprabasal-epidermis genes that 4T1 cells do not make, so this script asks which
tumours contain skin (and skeletal muscle), and whether the tight-junction / claudin
and immune results survive removing the skin-positive tumours.

Outputs (tables/):
  qc_tissue_contamination_per_sample.csv   marker counts, scores and skin flag per tumour
  contamination_within_WT_spearman.csv     focal genes vs skin score within the 5 WT tumours
  paper_S4D_style_tests.csv                 Cldn1/Cldn7/Ocln etc. on normalized counts, with/without skin+ WT
  panel_score_tests_skin_free_WT.csv        panel scores, 2 skin-free WT vs 5 KO (exact, 21 labelings)
  gsea_paper_sets_leading_edge.csv          leading-edge genes of the paper's WT-enriched sets
  gsea_paper_sets_skin_free_WT.csv          paper / context sets: without skin-positive WT, and with
                                            contamination-tracking members pruned (all 10 tumours)
  contamination_tracking_genes_rho0.8.csv   genes pruned for the latter
"""
import importlib
import itertools
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "scripts"))
from gene_sets import (CONTEXT_SETS, PANELS, PAPER_FIG4A, SKELETAL_MUSCLE, SKIN_EPIDERMIS,  # noqa: E402
                       SKIN_FLAG_GENES, SKIN_FLAG_MIN_SUM)

gsea = importlib.import_module("04_gsea")
RAW, CACHE, TAB = (os.path.join(HERE, d) for d in ("raw", "cache", "tables"))
RHO_TRACKING = 0.8

FOCAL = ["Tacstd2", "Epcam", "Cldn1", "Cldn3", "Cldn4", "Cldn7", "Ocln", "Cgn", "Tjp2", "Marveld2",
         "Krt5", "Krt14", "Krt8", "Krt18", "Cdh1", "Esrp1", "Grhl2", "Rab25",
         "Ptprc", "Cd8a", "Cd8b1", "Gzmb", "Prf1", "Nkg7", "Ifng", "Cxcl9", "B2m", "Stat1"]


def exact_perm(x, is_ko):
    obs = x[is_ko].mean() - x[~is_ko].mean()
    diffs = []
    for ko in itertools.combinations(range(len(x)), int(is_ko.sum())):
        m = np.zeros(len(x), bool); m[list(ko)] = True
        diffs.append(x[m].mean() - x[~m].mean())
    diffs = np.array(diffs)
    return float(np.mean(np.abs(diffs) >= abs(obs) - 1e-12)), len(diffs)


def leading_edge(metric, genes, members):
    order = np.argsort(-metric)
    ranked, r = genes[order], metric[order]
    inset = np.isin(ranked, list(members))
    w = np.abs(r) * inset
    run = np.cumsum(w) / w.sum() - np.cumsum(~inset) / (~inset).sum()
    if abs(run.max()) >= abs(run.min()):
        k = int(run.argmax()); return list(ranked[:k + 1][inset[:k + 1]]), "WT"
    k = int(run.argmin()); return list(ranked[k:][inset[k:]]), "KO"


def main():
    norm = pd.read_csv(os.path.join(RAW, "GSE334497_normalized_counts.csv.gz"), index_col=0)
    sheet = pd.read_csv(os.path.join(TAB, "sample_sheet.csv")).set_index("library_name").loc[norm.columns]
    ann = pd.read_csv(os.path.join(RAW, "gene_annotation_ensembl102.tsv.gz"), sep="\t", index_col=0)
    chip = pd.read_csv(os.path.join(CACHE, "Mouse_Ensembl_Gene_ID_Human_Orthologs_MSigDB.v2024.1.Hs.chip"),
                       sep="\t", index_col=0)
    m1 = pd.read_csv(os.path.join(TAB, "de_deseq2_M1_all10_genotype.csv.gz")).set_index("ensembl_gene_id")
    m5 = pd.read_csv(os.path.join(TAB, "de_deseq2_M5_drop_skin_positive_WT.csv.gz")).set_index("ensembl_gene_id")
    scores = pd.read_csv(os.path.join(TAB, "panel_scores_per_sample_log2.csv"), index_col=0).loc[norm.columns]

    sym2id = {}
    for gid, g in ann.gene_name.items():
        if isinstance(g, str) and g not in sym2id:
            sym2id[g] = gid
    lg = np.log2(norm + 1)

    def score(genes):
        return lg.loc[[sym2id[g] for g in genes if g in sym2id]].mean()

    per = pd.DataFrame(index=norm.columns)
    per["genotype"] = sheet.genotype
    per["library_wave"] = sheet.library_wave
    per["skin_flag_sum_norm"] = norm.loc[[sym2id[g] for g in SKIN_FLAG_GENES]].sum()
    per["skin_positive"] = per.skin_flag_sum_norm >= SKIN_FLAG_MIN_SUM
    per["skin_epidermis_score"] = score(SKIN_EPIDERMIS)
    per["skeletal_muscle_score"] = score(SKELETAL_MUSCLE)
    for g in SKIN_FLAG_GENES + ["Krt10", "Lor", "Krtdap", "Krt5", "Krt14", "Acta1", "Ckm", "Myh1",
                                "Tacstd2", "Cldn1", "Cldn3", "Cldn4", "Cldn7", "Ocln"]:
        per[f"norm_{g}"] = norm.loc[sym2id[g]].round(1)
    per.rename_axis("library_name").to_csv(os.path.join(TAB, "qc_tissue_contamination_per_sample.csv"))
    ct = pd.crosstab(per.genotype, per.skin_positive)
    fisher_p = stats.fisher_exact(ct.reindex(columns=[True, False], fill_value=0).values)[1]

    wt = per.genotype == "WT"
    rows = []
    for g in FOCAL:
        x = lg.loc[sym2id[g], wt.values].values
        for sc in ["skin_epidermis_score", "skeletal_muscle_score"]:
            rho, p = stats.spearmanr(per.loc[wt, sc].values, x)
            rows.append({"gene": g, "score": sc, "n_WT": int(wt.sum()), "spearman_rho_within_WT": rho,
                         "spearman_p": p})
    sp = pd.DataFrame(rows)
    sp.to_csv(os.path.join(TAB, "contamination_within_WT_spearman.csv"), index=False)

    # Paper S4D style: DESeq2-normalized counts, unpaired two-group comparison.
    rows = []
    ko = per.genotype == "KO"
    for g in ["Tacstd2", "Cldn1", "Cldn3", "Cldn4", "Cldn7", "Ocln"]:
        v = norm.loc[sym2id[g]]
        for lab, wmask in [("all_WT", wt), ("skin_free_WT", wt & ~per.skin_positive)]:
            a, b = v[ko].values, v[wmask].values
            rows.append({"gene": g, "WT_group": lab, "n_WT": len(b), "n_KO": len(a),
                         "mean_norm_WT": b.mean(), "mean_norm_KO": a.mean(),
                         "log2_ratio_KO_over_WT_means": np.log2(a.mean() / b.mean()),
                         "welch_p_norm_counts": stats.ttest_ind(a, b, equal_var=False).pvalue,
                         "student_t_p_norm_counts": stats.ttest_ind(a, b).pvalue,
                         "mwu_exact_p_norm_counts": stats.mannwhitneyu(a, b, alternative="two-sided",
                                                                       method="exact").pvalue,
                         "welch_p_log2_norm": stats.ttest_ind(np.log2(a + 1), np.log2(b + 1),
                                                              equal_var=False).pvalue})
    pd.DataFrame(rows).to_csv(os.path.join(TAB, "paper_S4D_style_tests.csv"), index=False)

    keep = ~per.skin_positive.values
    is_ko = (per.genotype.values == "KO")[keep]
    rows = []
    cols = list(PANELS) + ["DIAG_epithelial_differentiation", "DIAG_epithelial_content_Krt8_Krt18",
                           "DIAG_basal_keratins_Krt5_Krt14"]
    for p in cols:
        x = scores[p].values[keep]
        k, w = x[is_ko], x[~is_ko]
        pp, nperm = exact_perm(x, is_ko)
        all_diff = scores[p][per.genotype == "KO"].mean() - scores[p][per.genotype == "WT"].mean()
        rows.append({"panel": p, "n_WT_skin_free": int((~is_ko).sum()), "n_KO": int(is_ko.sum()),
                     "diff_KO_minus_WT_all10": all_diff, "diff_KO_minus_skin_free_WT": k.mean() - w.mean(),
                     "welch_p": stats.ttest_ind(k, w, equal_var=False).pvalue,
                     "mwu_exact_p": stats.mannwhitneyu(k, w, alternative="two-sided", method="exact").pvalue,
                     "perm_exact_p": pp, "n_permutations": nperm,
                     "min_attainable_perm_p": 1 / nperm})
    pst = pd.DataFrame(rows)
    pst.to_csv(os.path.join(TAB, "panel_score_tests_skin_free_WT.csv"), index=False)

    hs = gsea.collapse_to_human(norm, chip)
    genes = np.asarray(hs.index)
    X = hs.values.astype(float)
    is_wt = (sheet.genotype == "WT").values
    metric = gsea.s2n(X, is_wt)
    names = {c for _, _, _, cands in PAPER_FIG4A for c in cands} | set(CONTEXT_SETS)
    target = {}
    for coll in gsea.COLLECTIONS:
        target.update({k: v for k, v in gsea.read_gmt(os.path.join(CACHE, f"{coll}.v2024.1.Hs.symbols.gmt")).items()
                       if k in names})
    skin_hs = set(chip.loc[[sym2id[g] for g in SKIN_EPIDERMIS if sym2id.get(g) in chip.index], "Gene Symbol"])
    muscle_hs = set(chip.loc[[sym2id[g] for g in SKELETAL_MUSCLE if sym2id.get(g) in chip.index], "Gene Symbol"])
    hs_to_mouse = chip["Gene Symbol"].dropna()
    le_rows = []
    for label, nes_rep, _, cands in PAPER_FIG4A:
        if nes_rep <= 0:
            continue
        for c in cands:
            le, side = leading_edge(metric, genes, target[c])
            for g in le:
                ids = [i for i in hs_to_mouse.index[hs_to_mouse == g] if i in m1.index]
                gid = max(ids, key=lambda i: m1.loc[i, "baseMean"]) if ids else None
                le_rows.append({
                    "paper_label": label, "msigdb_set": c, "leading_edge_side": side, "human_gene": g,
                    "mouse_gene": ann.loc[gid, "gene_name"] if gid else None,
                    "s2n_WT_vs_KO": float(metric[genes == g][0]),
                    "M1_log2FC_KO_vs_WT": m1.loc[gid, "log2FC_KO_vs_WT"] if gid else np.nan,
                    "M1_p": m1.loc[gid, "pvalue"] if gid else np.nan,
                    "M5_skin_free_log2FC_KO_vs_WT": m5.loc[gid, "log2FC_KO_vs_WT"] if gid in m5.index else np.nan,
                    "M5_skin_free_p": m5.loc[gid, "pvalue"] if gid in m5.index else np.nan,
                    "skin_epidermis_marker": g in skin_hs, "skeletal_muscle_marker": g in muscle_hs,
                    "norm_skin_positive_WT_mean": norm.loc[gid, per.index[per.skin_positive]].mean() if gid else np.nan,
                    "norm_skin_free_WT_mean": norm.loc[gid, per.index[wt & ~per.skin_positive]].mean() if gid else np.nan,
                    "norm_KO_mean": norm.loc[gid, per.index[per.genotype == "KO"]].mean() if gid else np.nan,
                })
    le_tab = pd.DataFrame(le_rows)
    le_tab.to_csv(os.path.join(TAB, "gsea_paper_sets_leading_edge.csv"), index=False)

    res_sf, _ = gsea.run_phenotype_gsea(X[:, keep], list(genes), is_wt[keep], target, label="skin_free_WT")
    full = pd.read_csv(os.path.join(TAB, "gsea_phenotype_perm_all_sets.csv.gz")).drop_duplicates("set").set_index("set")
    stat_m5 = -m5["wald_stat"].dropna()
    pr = pd.DataFrame({"hs": chip["Gene Symbol"].reindex(stat_m5.index), "stat": stat_m5}).dropna()
    pr = pr.loc[pr.groupby("hs")["stat"].apply(lambda s: s.abs().idxmax())]
    r_pre = gsea.run_preranked(pr["stat"].values, pr["hs"].values, target).set_index("set")
    out = []
    paper_nes = {c: nes for _, nes, _, cands in PAPER_FIG4A for c in cands}
    for c in sorted(target, key=lambda s: (s not in paper_nes, s)):
        a = res_sf[res_sf.set == c]
        out.append({"msigdb_set": c, "in_paper_fig4A": c in paper_nes, "paper_NES": paper_nes.get(c, np.nan),
                    "NES_all10": full.loc[c, "NES"] if c in full.index else np.nan,
                    "p_all10": full.loc[c, "nominal_p"] if c in full.index else np.nan,
                    "NES_skin_free_WT_2v5": a.NES.iloc[0] if len(a) else np.nan,
                    "p_skin_free_WT_2v5_21_labelings": a.nominal_p.iloc[0] if len(a) else np.nan,
                    "NES_preranked_M5_DESeq2": r_pre.loc[c, "NES"] if c in r_pre.index else np.nan,
                    "p_preranked_M5_DESeq2_anticonservative": r_pre.loc[c, "nominal_p"] if c in r_pre.index else np.nan})
    sf = pd.DataFrame(out)

    # Prune set members that track either contamination score across all 10 tumours, keep the
    # ranking unchanged, and re-run the exact 252-labeling GSEA on all samples.
    rk = np.apply_along_axis(stats.rankdata, 1, X)
    rk = (rk - rk.mean(1, keepdims=True)) / rk.std(1, keepdims=True).clip(1e-12)
    tracking = np.zeros(len(genes), bool)
    for sc in ["skin_epidermis_score", "skeletal_muscle_score"]:
        s = stats.rankdata(per.loc[hs.columns, sc].values)
        s = (s - s.mean()) / s.std()
        tracking |= (rk @ s) / len(s) >= RHO_TRACKING
    tracked = set(genes[tracking])
    pruned = {c: v - tracked for c, v in target.items()}
    res_pr, _ = gsea.run_phenotype_gsea(X, list(genes), is_wt, pruned, min_size=10, label="pruned")
    res_pr = res_pr.set_index("set")
    sf["n_members_in_data"] = [len(target[c] & set(genes)) for c in sf.msigdb_set]
    sf["n_contamination_tracking_removed"] = [len(target[c] & tracked) for c in sf.msigdb_set]
    sf["NES_all10_pruned"] = [res_pr.NES.get(c, np.nan) for c in sf.msigdb_set]
    sf["p_all10_pruned"] = [res_pr.nominal_p.get(c, np.nan) for c in sf.msigdb_set]
    sf.to_csv(os.path.join(TAB, "gsea_paper_sets_skin_free_WT.csv"), index=False)
    pd.DataFrame({"human_gene": sorted(tracked)}).to_csv(
        os.path.join(TAB, "contamination_tracking_genes_rho0.8.csv"), index=False)
    print(f"contamination-tracking genes (Spearman >= {RHO_TRACKING} with skin or muscle score, n=10): "
          f"{len(tracked)} of {len(genes)}")

    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
    print(per.round(2).to_string())
    print("skin-positive by genotype:\n", ct, "\nFisher exact p =", round(fisher_p, 4))
    print(sp[sp.score == "skin_epidermis_score"].round(3).to_string(index=False))
    print(pst.round(4).to_string(index=False))
    summ = le_tab.groupby("msigdb_set").agg(
        n_leading_edge=("human_gene", "size"), n_skin=("skin_epidermis_marker", "sum"),
        n_muscle=("skeletal_muscle_marker", "sum"),
        n_lost_when_skin_free=("M5_skin_free_log2FC_KO_vs_WT", lambda s: int((s > -0.3).sum())))
    print(summ.to_string())
    print(le_tab[["msigdb_set", "human_gene", "mouse_gene", "M1_log2FC_KO_vs_WT", "M5_skin_free_log2FC_KO_vs_WT",
                  "norm_skin_positive_WT_mean", "norm_skin_free_WT_mean", "norm_KO_mean"]].round(2).to_string(index=False))
    print(sf.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
