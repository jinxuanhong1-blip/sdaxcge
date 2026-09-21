#!/usr/bin/env python3
"""Compare scTenifoldKnk CLDN4 shifts to the observed CLDN4-high vs low contrast.

Pre-specified before looking at the knockout hits:
- scTenifoldKnk distance / Z / FC is a regulatory-manifold shift, not an expression logFC.
  The knocked-out gene is excluded from every enrichment test.
- Primary observed number is the cohort-adjusted Q4 vs Q1 family OLS (positive = higher in CLDN4-high).
- Primary KO test is a two-sided Wilcoxon of meta-Z (Stouffer across cohorts) for family vs other genes
  in the shared universe (gene in at least two cohorts).
- Signed agreement uses the denoised outgoing edge WT[CLDN4, target], not the manifold distance.
- A ribosomal control knockout on GSE131907 uses the same count matrix as the CLDN4 run.
  Preference order: RPL13A, RPLP0, RPS18 if that gene is in the analyzed universe and outside
  IFN/MHC/TJ. If none are, the alphabetically first background gene matching ^RPL or ^RPS.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from contrastlib import ROOT, ols_matrix

FAMILIES = {
    "IFN": ("in_IFN", -1),
    "MHC-I/APM": ("in_MHC", -1),
    "TJ": ("in_TJ", 1),
}
CONTROL_CANDIDATES = ["RPL13A", "RPLP0", "RPS18"]


def bh(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out.tolist()
    order = np.argsort(p[ok])
    ranked = p[ok][order]
    m = len(ranked)
    adj = np.empty(m)
    prev = 1.0
    for i in range(m - 1, -1, -1):
        prev = min(prev, ranked[i] * m / (i + 1))
        adj[i] = prev
    back = np.empty(m)
    back[order] = np.clip(adj, 0, 1)
    out[np.flatnonzero(ok)] = back
    return out.tolist()


def wilcox(a: np.ndarray, b: np.ndarray):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 3 or len(b) < 3:
        return dict(n_family=len(a), n_background=len(b), rank_biserial=np.nan, p=np.nan, median_family=np.nan, median_background=np.nan)
    res = stats.mannwhitneyu(a, b, alternative="two-sided")
    r = 2.0 * float(res.statistic) / (len(a) * len(b)) - 1.0
    return dict(
        n_family=int(len(a)),
        n_background=int(len(b)),
        rank_biserial=r,
        p=float(res.pvalue),
        median_family=float(np.median(a)),
        median_background=float(np.median(b)),
    )


def load_ko(path: Path, gko: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    dr = pd.read_csv(path / "diffregulation.tsv", sep="\t")
    edges = pd.read_csv(path / "outgoing_edges.tsv", sep="\t")
    meta = pd.read_csv(path / "run_meta.tsv", sep="\t").iloc[0]
    dr["gene"] = dr["gene"].astype(str)
    edges["gene"] = edges["gene"].astype(str)
    dr = dr.loc[dr["gene"] != gko].copy()
    return dr, edges, meta


def stouffer(z_by_cohort: dict[str, pd.Series]) -> pd.DataFrame:
    genes = sorted(set().union(*[set(s.index) for s in z_by_cohort.values()]))
    rows = []
    for g in genes:
        zs = []
        dists = []
        for ds, s in z_by_cohort.items():
            if g in s.index:
                zs.append(float(s.loc[g]))
        if not zs:
            continue
        k = len(zs)
        z = float(np.sum(zs) / np.sqrt(k))
        p = float(2 * stats.norm.sf(abs(z)))
        rows.append({"gene": g, "k_cohorts": k, "meta_Z": z, "meta_p": p})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", type=Path, default=Path("/tmp/c4_work"))
    ap.add_argument("--cohorts", nargs="+", default=["GSE123902", "GSE131907", "GSE205335", "GSE189357"])
    args = ap.parse_args()
    out = ROOT / "results" / "tables"
    figdir = ROOT / "results" / "figures"
    out.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    pack = np.load(args.work / "logcpm_de_units.npz", allow_pickle=True)
    genes = pack["genes"].astype(str)
    logcpm = pack["logcpm"]
    dataset = pack["dataset"].astype(str)
    quartile = pack["quartile"].astype(str)
    qmask = np.isin(quartile, ["Q1", "Q4"])
    logfc, pvals, df = ols_matrix(logcpm[:, qmask], dataset[qmask], quartile[qmask] == "Q4")
    observed = pd.DataFrame({"gene": genes, "obs_logFC_high_minus_low": logfc, "obs_p": pvals})

    universes = []
    zmaps = {}
    edge_maps = {}
    metas = []
    per_cohort_tests = []
    for ds in args.cohorts:
        uni = pd.read_csv(args.work / ds / "gene_universe.tsv", sep="\t")
        uni["dataset"] = ds
        universes.append(uni)
        dr, edges, meta = load_ko(args.work / ds / "knk_CLDN4", "CLDN4")
        meta = meta.copy()
        meta["dataset"] = ds
        metas.append(meta)
        dr.to_csv(out / f"knk_diffregulation_{ds}.tsv", sep="\t", index=False)
        edges.to_csv(out / f"knk_edges_{ds}.tsv", sep="\t", index=False)
        zmaps[ds] = dr.set_index("gene")["Z"]
        edge_maps[ds] = edges.set_index("gene")["weight"]
        merged = uni.merge(dr, on="gene", how="left").merge(observed, on="gene", how="left")
        merged = merged.loc[merged["gene"] != "CLDN4"]
        for fam, (col, _sign) in FAMILIES.items():
            w = wilcox(merged.loc[merged[col], "Z"], merged.loc[~merged["in_family"], "Z"])
            w.update({"dataset": ds, "family": fam, "test": "wilcoxon_Z_family_vs_nonfamily"})
            per_cohort_tests.append(w)

    meta_z = stouffer(zmaps)
    # mean edge across cohorts where the gene was in that network
    edge_rows = []
    all_edge_genes = sorted(set().union(*[set(s.index) for s in edge_maps.values()]))
    for g in all_edge_genes:
        if g == "CLDN4":
            continue
        ws = [float(edge_maps[ds].loc[g]) for ds in edge_maps if g in edge_maps[ds].index]
        edge_rows.append({"gene": g, "k_edge": len(ws), "mean_edge": float(np.mean(ws))})
    edge_df = pd.DataFrame(edge_rows)

    uni_all = pd.concat(universes, ignore_index=True)
    flags = (
        uni_all.groupby("gene")[["in_IFN", "in_MHC", "in_TJ", "in_family"]]
        .any()
        .reset_index()
    )
    tab = meta_z.merge(edge_df, on="gene", how="outer").merge(flags, on="gene", how="left").merge(observed, on="gene", how="left")
    for col in ["in_IFN", "in_MHC", "in_TJ", "in_family"]:
        tab[col] = tab[col].fillna(False).astype(bool)
    primary = tab.loc[(tab["k_cohorts"] >= 2) & (tab["gene"] != "CLDN4")].copy()
    primary.to_csv(out / "knk_meta_genes.tsv", sep="\t", index=False)
    top = primary.sort_values("meta_Z", ascending=False).head(40)
    top.to_csv(out / "knk_top_genes.tsv", sep="\t", index=False)

    fam_obs = pd.read_csv(out / "family_q4q1_ols.tsv", sep="\t")
    locked = pd.read_csv(ROOT / "data" / "locked_family_q4q1_ols.tsv", sep="\t")
    locked = locked.rename(columns={"logFC": "locked_logFC", "p": "locked_p", "n_genes": "locked_n_genes_present"})

    tests = []
    family_rows = []
    wilcox_ps = []
    family_order = []
    for fam, (col, expected_sign) in FAMILIES.items():
        sub_f = primary.loc[primary[col]]
        sub_b = primary.loc[~primary["in_family"]]
        w = wilcox(sub_f["meta_Z"], sub_b["meta_Z"])
        wilcox_ps.append(w["p"])
        family_order.append(fam)
        obs = fam_obs.loc[fam_obs["family"] == fam].iloc[0]
        lock = locked.loc[locked["family"] == fam].iloc[0]
        both = sub_f["mean_edge"].notna() & sub_f["obs_logFC_high_minus_low"].notna()
        ee = sub_f.loc[both, "mean_edge"].to_numpy()
        ff = sub_f.loc[both, "obs_logFC_high_minus_low"].to_numpy()
        nz = (ee != 0) & (ff != 0) & np.isfinite(ee) & np.isfinite(ff)
        n_nz = int(nz.sum())
        n_agree = int(np.sum(np.sign(ee[nz]) == np.sign(ff[nz]))) if n_nz else 0
        binom_p = float(stats.binomtest(n_agree, n_nz, 0.5).pvalue) if n_nz else np.nan
        family_rows.append(
            {
                "family": fam,
                "expected_sign_high_minus_low": expected_sign,
                "observed_logFC_high_minus_low": obs["logFC_high_minus_low"],
                "observed_p": obs["p"],
                "observed_n_q1": obs["n_q1"],
                "observed_n_q4": obs["n_q4"],
                "observed_n_genes_present": obs["n_genes_present"],
                "locked_logFC": lock["locked_logFC"],
                "locked_p": lock["locked_p"],
                "abs_diff_vs_locked_logFC": abs(obs["logFC_high_minus_low"] - lock["locked_logFC"]),
                "observed_matches_expected_sign": bool(np.sign(obs["logFC_high_minus_low"]) == expected_sign),
                "ko_n_family": w["n_family"],
                "ko_n_background": w["n_background"],
                "ko_median_Z_family": w["median_family"],
                "ko_median_Z_background": w["median_background"],
                "ko_rank_biserial": w["rank_biserial"],
                "ko_wilcoxon_p": w["p"],
                "edge_sign_n": n_nz,
                "edge_sign_n_agree": n_agree,
                "edge_sign_concordance": (n_agree / n_nz) if n_nz else np.nan,
                "edge_sign_binom_p": binom_p,
                "mean_edge_in_family": float(np.nanmean(ee)) if len(ee) else np.nan,
                "mean_gene_logFC_in_family": float(np.nanmean(ff)) if len(ff) else np.nan,
            }
        )
    qvals = bh(wilcox_ps)
    for row, q in zip(family_rows, qvals):
        row["ko_wilcoxon_q_bh3"] = q
    fam_tab = pd.DataFrame(family_rows)
    fam_tab.to_csv(out / "family_ko_vs_observed.tsv", sep="\t", index=False)

    # global tests on the primary gene set
    sp_abs = stats.spearmanr(primary["meta_Z"], np.abs(primary["obs_logFC_high_minus_low"]), nan_policy="omit")
    sp_edge = stats.spearmanr(primary["mean_edge"], primary["obs_logFC_high_minus_low"], nan_policy="omit")
    ee = primary["mean_edge"].to_numpy()
    ff = primary["obs_logFC_high_minus_low"].to_numpy()
    nz = np.isfinite(ee) & np.isfinite(ff) & (ee != 0) & (ff != 0)
    n_nz = int(nz.sum())
    n_agree = int(np.sum(np.sign(ee[nz]) == np.sign(ff[nz])))
    z = primary["meta_Z"].to_numpy()
    finite = np.isfinite(z) & np.isfinite(ff)
    zf = z[finite]
    ff_f = ff[finite]
    hi = zf >= np.quantile(zf, 0.9)
    up_high = ff_f >= np.quantile(ff_f, 0.9)
    up_low = ff_f <= np.quantile(ff_f, 0.1)
    big = np.abs(ff_f) >= np.quantile(np.abs(ff_f), 0.9)

    def fisher(a, b):
        tp = int(np.sum(a & b))
        fn = int(np.sum(a & ~b))
        fp = int(np.sum(~a & b))
        tn = int(np.sum(~a & ~b))
        odds, p = stats.fisher_exact([[tp, fn], [fp, tn]])
        return tp, int(a.sum()), int(b.sum()), float(odds), float(p)

    tests.append({"test": "spearman_metaZ_vs_abs_obs_logFC", "n": int(sp_abs.correlation.size if False else primary[["meta_Z", "obs_logFC_high_minus_low"]].dropna().shape[0]), "stat": float(sp_abs.statistic), "p": float(sp_abs.pvalue)})
    tests.append({"test": "spearman_mean_edge_vs_obs_logFC", "n": int(primary[["mean_edge", "obs_logFC_high_minus_low"]].dropna().shape[0]), "stat": float(sp_edge.statistic), "p": float(sp_edge.pvalue)})
    tests.append({"test": "edge_sign_concordance_all_genes", "n": n_nz, "stat": n_agree / n_nz if n_nz else np.nan, "p": float(stats.binomtest(n_agree, n_nz, 0.5).pvalue) if n_nz else np.nan, "n_agree": n_agree})
    for label, mask in [("top_decile_abs_logFC", big), ("top_decile_up_in_high", up_high), ("top_decile_up_in_low", up_low)]:
        tp, na, nb, odds, p = fisher(hi, mask)
        tests.append({"test": f"fisher_ko_top_decile_vs_{label}", "n": int(finite.sum()), "stat": odds, "p": p, "n_overlap": tp, "n_ko": na, "n_observed": nb})

    # Sensitivity: family genes at least as variable as the least-variable background gene.
    # Force-added low-variance family members would otherwise lower the family Z.
    sens_rows = []
    for ds in args.cohorts:
        uni = pd.read_csv(args.work / ds / "gene_universe.tsv", sep="\t")
        dr = pd.read_csv(out / f"knk_diffregulation_{ds}.tsv", sep="\t")
        merged = uni.merge(dr, on="gene", how="inner")
        merged = merged.loc[merged["gene"] != "CLDN4"]
        cutoff = float(merged.loc[merged["role"] == "background", "variance"].min())
        background_z = merged.loc[merged["role"] == "background", "Z"]
        for fam, (col, _sign) in FAMILIES.items():
            fam_z = merged.loc[merged[col] & (merged["variance"] >= cutoff), "Z"]
            w = wilcox(fam_z, background_z)
            w.update({
                "dataset": ds,
                "family": fam,
                "n_family_all": int(merged[col].sum()),
                "variance_cutoff": cutoff,
                "test": "wilcoxon_Z_family_variance_at_least_background_min",
            })
            sens_rows.append(w)
    pd.DataFrame(sens_rows).to_csv(out / "knk_family_wilcoxon_variance_matched.tsv", sep="\t", index=False)

    fdr_rows = []
    for ds in args.cohorts:
        dr = pd.read_csv(out / f"knk_diffregulation_{ds}.tsv", sep="\t")
        sig = dr.loc[dr["p.adj"] < 0.05].copy()
        sig["dataset"] = ds
        fdr_rows.append(sig)
    fdr = pd.concat(fdr_rows, ignore_index=True)
    fdr = fdr.merge(observed, on="gene", how="left")
    fdr = fdr.merge(flags, on="gene", how="left")
    fdr.to_csv(out / "knk_fdr05_genes.tsv", sep="\t", index=False)

    # control knockout, GSE131907 only, same analyzed universe as CLDN4
    uni131 = pd.read_csv(args.work / "GSE131907" / "gene_universe.tsv", sep="\t")
    blocked = set(uni131.loc[uni131["in_family"], "gene"]) | {"CLDN4"}
    present = set(uni131["gene"])
    ctrl_name = None
    ctrl_rule = None
    for cand in CONTROL_CANDIDATES:
        if cand in present and cand not in blocked:
            ctrl_name = cand
            ctrl_rule = "named_candidate"
            break
    if ctrl_name is None:
        ribo = uni131.loc[
            (uni131["role"] == "background")
            & uni131["gene"].str.match(r"^(RPL|RPS)")
            & ~uni131["gene"].isin(blocked)
        ].sort_values("gene")
        if ribo.empty:
            raise SystemExit("no RPL/RPS control gene in the GSE131907 background universe")
        ctrl_name = str(ribo.iloc[0]["gene"])
        ctrl_rule = "alphabetical_RPL_or_RPS_background"
    ctrl_dir = args.work / "GSE131907" / f"knk_{ctrl_name}"
    if ctrl_dir is not None and (ctrl_dir / "diffregulation.tsv").exists():
        dr_c, _, meta_c = load_ko(ctrl_dir, ctrl_name)
        meta_c = meta_c.copy()
        meta_c["dataset"] = "GSE131907"
        meta_c["gKO"] = ctrl_name
        metas.append(meta_c)
        merged = uni131.merge(dr_c, on="gene", how="inner")
        merged = merged.loc[~merged["gene"].isin([ctrl_name, "CLDN4"])]
        dr_cl = pd.read_csv(out / "knk_diffregulation_GSE131907.tsv", sep="\t")
        cl_merged = uni131.merge(dr_cl, on="gene", how="inner")
        cl_merged = cl_merged.loc[cl_merged["gene"] != "CLDN4"]
        for fam, (col, _sign) in FAMILIES.items():
            wc = wilcox(merged.loc[merged[col], "Z"], merged.loc[~merged["in_family"], "Z"])
            tests.append({"test": f"control_{ctrl_name}_GSE131907_{fam}_wilcoxon", "n": wc["n_family"], "stat": wc["rank_biserial"], "p": wc["p"], "median_family": wc["median_family"], "median_background": wc["median_background"]})
            wcl = wilcox(cl_merged.loc[cl_merged[col], "Z"], cl_merged.loc[~cl_merged["in_family"], "Z"])
            tests.append({"test": f"CLDN4_GSE131907_{fam}_wilcoxon", "n": wcl["n_family"], "stat": wcl["rank_biserial"], "p": wcl["p"], "median_family": wcl["median_family"], "median_background": wcl["median_background"]})
        both = dr_c.merge(dr_cl[["gene", "Z"]], on="gene", suffixes=("_control", "_CLDN4"))
        both = both.loc[~both["gene"].isin([ctrl_name, "CLDN4"])]
        rho, rp = stats.spearmanr(both["Z_control"], both["Z_CLDN4"])
        tests.append({"test": f"spearman_Z_{ctrl_name}_vs_CLDN4_GSE131907", "n": int(len(both)), "stat": float(rho), "p": float(rp)})
        dr_c.to_csv(out / f"knk_diffregulation_GSE131907_{ctrl_name}.tsv", sep="\t", index=False)
    else:
        tests.append({"test": "control_knockout", "n": 0, "stat": np.nan, "p": np.nan, "note": f"{ctrl_name} via {ctrl_rule}: output not found"})

    pd.DataFrame(per_cohort_tests).to_csv(out / "knk_family_wilcoxon_by_cohort.tsv", sep="\t", index=False)
    pd.DataFrame(tests).to_csv(out / "comparison_tests.tsv", sep="\t", index=False)
    pd.DataFrame(metas).to_csv(out / "knk_run_meta.tsv", sep="\t", index=False)

    # figure: observed logFC and KO shift
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.6))
    names = fam_tab["family"].tolist()
    x = np.arange(len(names))
    axes[0].bar(x, fam_tab["observed_logFC_high_minus_low"], color="#4C78A8", label="this run")
    axes[0].scatter(x, fam_tab["locked_logFC"], color="#F58518", zorder=3, label="locked OLS")
    axes[0].axhline(0, color="black", lw=0.6)
    axes[0].set_xticks(x, names)
    axes[0].set_ylabel("logFC CLDN4-high minus low")
    axes[0].set_title("Observed malignant contrast")
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].bar(x, fam_tab["ko_rank_biserial"], color="#54A24B")
    axes[1].axhline(0, color="black", lw=0.6)
    axes[1].set_xticks(x, names)
    axes[1].set_ylabel("rank-biserial of KO meta-Z")
    axes[1].set_title("scTenifoldKnk shift vs other genes")
    fig.tight_layout()
    fig.savefig(figdir / "family_ko_vs_observed.png", dpi=140)
    fig.savefig(figdir / "family_ko_vs_observed.pdf")
    plt.close(fig)

    summary = {
        "n_primary_genes": int(len(primary)),
        "control_gene": ctrl_name,
        "control_rule": ctrl_rule,
        "families": family_rows,
        "tests": tests,
        "run_meta": [m.to_dict() for m in metas],
    }
    (ROOT / "results" / "summary.json").write_text(json.dumps(summary, indent=2, default=float))
    print(fam_tab.to_string(index=False), flush=True)
    print(pd.DataFrame(tests).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
