"""Sample-level gene-panel scores, exact tests and sensitivity analyses (KO vs WT).

Score = mean over panel genes of log2(DESeq2-normalized count + 1); the KO-WT
difference of this score equals the panel-average log2 fold change. A z-score
version (genes standardised across the 10 samples, then averaged) is also given.

Tests per panel (n = 5 KO vs 5 WT):
  * Welch two-sample t-test
  * exact two-sided Mann-Whitney U
  * exact permutation test of the mean difference over all C(10,5) = 252 labelings
  * OLS score ~ genotype + library_wave (putative RESUB submission wave)
  * leave-one-sample-out range of the difference
"""
import itertools
import os
import sys

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "scripts"))
from gene_sets import HUMAN_TO_MOUSE_NOTES, PANELS  # noqa: E402

RAW, CACHE, TAB = (os.path.join(HERE, d) for d in ("raw", "cache", "tables"))

# Defined after the primary results showed Epcam/keratins lower in KO; used only
# as diagnostics for tumour-epithelial content / state, not as primary end points.
DIAGNOSTIC_PANELS = {
    "DIAG_epithelial_content_Krt8_Krt18": ["Krt8", "Krt18"],
    "DIAG_epithelial_differentiation": ["Cdh1", "Esrp1", "Esrp2", "Grhl2", "Ovol1", "Rab25",
                                        "St14", "Krt7", "Krt19"],
    "DIAG_basal_keratins_Krt5_Krt14": ["Krt5", "Krt14"],
    "DIAG_mesenchymal": ["Vim", "Zeb1", "Snai2", "Twist1", "Cdh2"],
    "DIAG_proliferation": ["Mki67", "Top2a", "Ccnb1", "Cdk1", "Birc5", "Mcm2"],
}


def exact_perm_p(x, labels):
    obs = x[labels == "KO"].mean() - x[labels == "WT"].mean()
    idx = np.arange(len(x))
    diffs = []
    for ko in itertools.combinations(idx, int((labels == "KO").sum())):
        m = np.zeros(len(x), bool); m[list(ko)] = True
        diffs.append(x[m].mean() - x[~m].mean())
    diffs = np.array(diffs)
    return float(np.mean(np.abs(diffs) >= abs(obs) - 1e-12)), len(diffs)


def test_score(score, sheet):
    lab = sheet["genotype"].values
    ko, wt = score[lab == "KO"], score[lab == "WT"]
    diff = ko.mean() - wt.mean()
    welch = stats.ttest_ind(ko, wt, equal_var=False)
    mwu = stats.mannwhitneyu(ko, wt, alternative="two-sided", method="exact")
    perm_p, n_perm = exact_perm_p(score, lab)
    df = pd.DataFrame({"s": score, "genotype": pd.Categorical(lab, ["WT", "KO"]),
                       "wave": sheet["library_wave"].values})
    ols = smf.ols("s ~ genotype + wave", df).fit()
    loo = []
    for i in range(len(score)):
        keep = np.arange(len(score)) != i
        k, w = score[keep][lab[keep] == "KO"], score[keep][lab[keep] == "WT"]
        loo.append(k.mean() - w.mean())
    return {
        "n_KO": len(ko), "n_WT": len(wt),
        "mean_WT": wt.mean(), "mean_KO": ko.mean(), "diff_KO_minus_WT": diff,
        "welch_t": welch.statistic, "welch_p": welch.pvalue,
        "mwu_U": mwu.statistic, "mwu_exact_p": mwu.pvalue,
        "perm_exact_p": perm_p, "n_permutations": n_perm,
        "wave_adjusted_diff": ols.params["genotype[T.KO]"], "wave_adjusted_p": ols.pvalues["genotype[T.KO]"],
        "loo_diff_min": min(loo), "loo_diff_max": max(loo),
        "loo_sign_stable": bool(np.all(np.sign(loo) == np.sign(diff))),
    }


def main():
    norm = pd.read_csv(os.path.join(RAW, "GSE334497_normalized_counts.csv.gz"), index_col=0)
    sheet = pd.read_csv(os.path.join(TAB, "sample_sheet.csv")).set_index("library_name").loc[norm.columns]
    ann = pd.read_csv(os.path.join(RAW, "gene_annotation_ensembl102.tsv.gz"), sep="\t", index_col=0)
    full_ann = pd.read_csv(os.path.join(CACHE, "ensembl102_mouse_genes.tsv"), sep="\t", index_col=0)
    de = pd.read_csv(os.path.join(TAB, "de_deseq2_M1_all10_genotype.csv.gz")).set_index("ensembl_gene_id")

    lg = np.log2(norm + 1)
    sym = ann["gene_name"]
    sym_to_id = {}
    for gid, g in sym.items():
        if isinstance(g, str) and g not in sym_to_id:
            sym_to_id[g] = gid

    panels = dict(PANELS); panels.update(DIAGNOSTIC_PANELS)
    membership, scores_log, scores_z = [], {}, {}
    z = lg.sub(lg.mean(axis=1), axis=0).div(lg.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    for name, genes in panels.items():
        present = [g for g in genes if g in sym_to_id]
        for g in genes:
            gid = sym_to_id.get(g)
            row = {"panel": name, "gene": g, "in_matrix": gid is not None, "ensembl_gene_id": gid}
            if gid is not None and gid in de.index:
                r = de.loc[gid]
                row.update({"baseMean": r.baseMean, "deseq2_log2FC_KO_vs_WT": r.log2FC_KO_vs_WT,
                            "deseq2_p": r.pvalue, "deseq2_padj": r.padj_BH,
                            "mean_norm_WT": norm.loc[gid, sheet.genotype == "WT"].mean(),
                            "mean_norm_KO": norm.loc[gid, sheet.genotype == "KO"].mean()})
            membership.append(row)
        ids = [sym_to_id[g] for g in present]
        scores_log[name] = lg.loc[ids].mean()
        scores_z[name] = z.loc[ids].mean()
    memb = pd.DataFrame(membership)
    memb.to_csv(os.path.join(TAB, "geneset_membership.csv"), index=False)

    S = pd.DataFrame(scores_log)
    SZ = pd.DataFrame(scores_z)
    epi = S["DIAG_epithelial_content_Krt8_Krt18"]
    epi_diff = S["DIAG_epithelial_differentiation"]
    S["TJ_claudin_program_minus_epithelial_content"] = S["TJ_claudin_program"] - epi
    S["Claudin_1_3_4_7_minus_epithelial_content"] = S["Claudin_1_3_4_7"] - epi
    S["Claudin_1_3_4_7_minus_epithelial_differentiation"] = S["Claudin_1_3_4_7"] - epi_diff
    S["CD8_T_minus_Leukocyte"] = S["CD8_T"] - S["Leukocyte"]
    S["Effector_minus_Leukocyte"] = S["Effector"] - S["Leukocyte"]
    S["Immune_core_minus_Leukocyte"] = S["Immune_core"] - S["Leukocyte"]
    out = S.join(sheet[["gsm", "genotype", "library_wave"]])
    out.index.name = "library_name"
    out.to_csv(os.path.join(TAB, "panel_scores_per_sample_log2.csv"))
    SZ.join(sheet[["genotype"]]).rename_axis("library_name").to_csv(
        os.path.join(TAB, "panel_scores_per_sample_zmean.csv"))

    rows = []
    for name in S.columns:
        res = test_score(S[name].values, sheet)
        res["panel"] = name
        res["score_type"] = "mean_log2_norm"
        res["n_genes_used"] = int(memb.loc[(memb.panel == name) & memb.in_matrix].shape[0]) if name in panels else np.nan
        if name in SZ.columns:
            zres = test_score(SZ[name].values, sheet)
            res["zmean_diff_KO_minus_WT"] = zres["diff_KO_minus_WT"]
            res["zmean_welch_p"] = zres["welch_p"]
        rows.append(res)
    T = pd.DataFrame(rows)
    T["is_diagnostic_posthoc"] = T["panel"].str.startswith("DIAG_") | T["panel"].str.contains("minus_epithelial")
    cols = ["panel", "is_diagnostic_posthoc", "score_type", "n_genes_used", "n_WT", "n_KO", "mean_WT", "mean_KO",
            "diff_KO_minus_WT", "welch_t", "welch_p", "mwu_U", "mwu_exact_p", "perm_exact_p",
            "n_permutations", "wave_adjusted_diff", "wave_adjusted_p", "loo_diff_min", "loo_diff_max",
            "loo_sign_stable", "zmean_diff_KO_minus_WT", "zmean_welch_p"]
    T = T[cols]
    T.to_csv(os.path.join(TAB, "panel_score_tests.csv"), index=False)

    adj = []
    for name in ["TJ_claudin_program", "Claudin_1_3_4_7", "Paper_TJ_trio_Cldn1_Cldn7_Ocln"]:
        for cov in ["DIAG_epithelial_content_Krt8_Krt18", "DIAG_epithelial_differentiation"]:
            df = pd.DataFrame({"s": S[name], "c": S[cov],
                               "genotype": pd.Categorical(sheet.genotype, ["WT", "KO"])})
            fit = smf.ols("s ~ genotype + c", df).fit()
            adj.append({"panel": name, "covariate": cov,
                        "genotype_coef_KO_vs_WT": fit.params["genotype[T.KO]"],
                        "genotype_p": fit.pvalues["genotype[T.KO]"],
                        "covariate_coef": fit.params["c"], "covariate_p": fit.pvalues["c"],
                        "r_score_vs_covariate": np.corrcoef(S[name], S[cov])[0, 1]})
    pd.DataFrame(adj).to_csv(os.path.join(TAB, "tj_scores_adjusted_for_epithelial_state.csv"), index=False)

    cl = full_ann[full_ann.gene_name.str.match(r"^Cldn\d", na=False)].copy()
    cl["in_deposited_matrix"] = cl.index.isin(norm.index)
    for c in ["baseMean", "log2FC_KO_vs_WT", "lfcSE", "pvalue", "padj_BH"]:
        cl[c] = de[c].reindex(cl.index)
    cl["mean_norm_WT"] = norm.reindex(cl.index).loc[:, (sheet.genotype == "WT").values].mean(axis=1)
    cl["mean_norm_KO"] = norm.reindex(cl.index).loc[:, (sheet.genotype == "KO").values].mean(axis=1)
    cl["n_WT"], cl["n_KO"] = 5, 5
    cl["test"] = np.where(cl.in_deposited_matrix, "DESeq2 Wald, ~genotype, recovered counts", "not in GEO matrix")
    cl = cl.sort_values(["in_deposited_matrix", "baseMean"], ascending=[False, False])
    cl.rename_axis("ensembl_gene_id").to_csv(os.path.join(TAB, "claudin_family_deseq2.csv"))

    pd.Series(HUMAN_TO_MOUSE_NOTES).rename("mouse_substitute").rename_axis("human_symbol").to_csv(
        os.path.join(TAB, "human_to_mouse_panel_notes.csv"))

    pd.set_option("display.width", 250)
    show = ["panel", "n_genes_used", "diff_KO_minus_WT", "welch_p", "mwu_exact_p", "perm_exact_p",
            "wave_adjusted_diff", "wave_adjusted_p", "loo_diff_min", "loo_diff_max", "loo_sign_stable"]
    print(T[show].round(4).to_string(index=False))
    print(pd.DataFrame(adj).round(4).to_string(index=False))
    print(cl[["gene_name", "in_deposited_matrix", "baseMean", "log2FC_KO_vs_WT", "pvalue", "padj_BH"]].round(4).to_string())


if __name__ == "__main__":
    main()
