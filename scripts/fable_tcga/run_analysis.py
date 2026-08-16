#!/usr/bin/env python3
"""TCGA LUAD / LUSC: TACSTD2 (TROP2) and CLDN4 vs immune contexture and survival.

Inputs (all open access, fetched by download_data.py):
  * Xena GDC hub STAR TPM matrices (log2(TPM+1), GENCODE v36)
  * TCGA-CDR curated survival endpoints + clinical covariates
  * MCP-counter marker-gene signatures (Becht et al. 2016)

Analyses (per cohort: LUAD, LUSC; primary tumors, one sample per patient):
  1. Tumor vs adjacent-normal expression (unpaired Mann-Whitney U + paired
     Wilcoxon on patients with matched pairs).
  2. Immune deconvolution: MCP-counter population scores (mean log2 TPM of
     marker genes), T-cell-inflamed GEP (Ayers et al., JCI 2017; mean of
     18 per-gene z-scores), cytolytic activity CYT (mean log2 of GZMA, PRF1),
     plus individual checkpoint / effector genes. Spearman correlations with
     TACSTD2 and CLDN4, BH-FDR corrected within each gene x cohort family.
  3. Quartile contrast: Q4 vs Q1 of each target gene, Mann-Whitney U on key
     immune scores.
  4. Survival: OS and PFI (TCGA-CDR). Kaplan-Meier median split + log-rank;
     Cox PH per-SD of log2 TPM, unadjusted and adjusted for age, sex and
     AJCC stage (III/IV vs I/II).

NOTE: TCGA patients were accrued before immune-checkpoint-inhibitor
approvals and were NOT ICI treated; immune associations here must not be
read as evidence of ICI benefit.

Usage:  python3 scripts/fable_tcga/run_analysis.py
Outputs: results/fable_tcga/  (CSV tables + PNG figures)
"""

import gzip
import json
import os
import warnings

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from statsmodels.stats.multitest import multipletests

DATA_DIR = os.environ.get("FABLE_TCGA_DATA", "/tmp/fable_tcga_data")
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(REPO, "results", "fable_tcga")
os.makedirs(OUT, exist_ok=True)

COHORTS = ["LUAD", "LUSC"]
TARGETS = ["TACSTD2", "CLDN4"]

# Ayers et al. JCI 2017 T-cell-inflamed 18-gene GEP
GEP18 = [
    "CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6",
    "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2",
    "PSMB10", "STAT1", "TIGIT",
]
CYT_GENES = ["GZMA", "PRF1"]
CHECKPOINTS = [
    "CD274", "PDCD1", "PDCD1LG2", "CTLA4", "LAG3", "TIGIT", "HAVCR2",
    "CD8A", "GZMB", "CXCL9",
]
KEY_SCORES = ["T cells", "CD8 T cells", "Cytotoxic lymphocytes", "GEP18", "CYT", "CD274"]
# tight-junction programme genes co-analysed with TACSTD2 (TROP2 physically
# engages claudins; Nakatsukasa et al. Am J Pathol 2010)
TJ_GENES = ["CLDN1", "CLDN4", "CLDN7", "F11R", "PARD3"]


def load_mcp_signatures():
    """MCP-counter marker genes; resolve symbols via GENCODE v36 Ensembl IDs
    where available (fixes a truncated 'C' entry = CAVIN2 in the upstream
    file and keeps symbols consistent with the expression matrix)."""
    probemap = pd.read_csv(
        os.path.join(DATA_DIR, "gencode.v36.probemap"), sep="\t", usecols=["id", "gene"]
    )
    ens2sym = {i.split(".")[0]: g for i, g in zip(probemap["id"], probemap["gene"])}
    sig = pd.read_csv(os.path.join(DATA_DIR, "mcpcounter_genes.txt"), sep="\t")
    sig.columns = [c.strip('"') for c in sig.columns]
    pops = {}
    for _, row in sig.iterrows():
        sym = ens2sym.get(str(row["ENSEMBL ID"]), str(row["HUGO symbols"]))
        pops.setdefault(row["Cell population"], set()).add(sym)
    return {pop: sorted(genes) for pop, genes in pops.items()}


def load_expression(cohort, symbols_needed):
    """Return log2(TPM+1) DataFrame (rows=gene symbols, cols=samples)."""
    probemap = pd.read_csv(
        os.path.join(DATA_DIR, "gencode.v36.probemap"), sep="\t", usecols=["id", "gene"]
    )
    id2sym = dict(zip(probemap["id"], probemap["gene"]))
    path = os.path.join(DATA_DIR, f"TCGA-{cohort}.star_tpm.tsv.gz")
    expr = pd.read_csv(path, sep="\t", index_col=0)
    expr.index = expr.index.map(lambda i: id2sym.get(i, i))
    expr = expr[expr.index.isin(symbols_needed)]
    # collapse duplicate symbols: keep the row with the highest mean expression
    expr = expr.loc[expr.mean(axis=1).groupby(expr.index).idxmax()]
    return expr


def split_samples(columns):
    """Return (tumor_samples one per patient, normal_samples)."""
    tumors, normals = {}, {}
    for s in sorted(columns):
        code = s.split("-")[3][:2]
        patient = "-".join(s.split("-")[:3])
        if code == "01" and patient not in tumors:
            tumors[patient] = s
        elif code == "11" and patient not in normals:
            normals[patient] = s
    return tumors, normals


def bh(pvals):
    return multipletests(pvals, method="fdr_bh")[1]


def partial_spearman(x, y, Z):
    """Rank-based partial correlation of x and y controlling for one or more
    covariates Z (1-D array or 2-D array with covariates in columns).
    Ranks of x and y are residualized on the ranks of Z (plus intercept);
    Pearson correlation of the residuals is returned with a t-test on
    n - 2 - k degrees of freedom."""
    Z = np.asarray(Z, dtype=float)
    if Z.ndim == 1:
        Z = Z[:, None]
    xr, yr = stats.rankdata(x), stats.rankdata(y)
    Zr = np.column_stack([stats.rankdata(Z[:, j]) for j in range(Z.shape[1])])
    design = np.column_stack([np.ones(len(xr)), Zr])
    rx = xr - design @ np.linalg.lstsq(design, xr, rcond=None)[0]
    ry = yr - design @ np.linalg.lstsq(design, yr, rcond=None)[0]
    r = np.corrcoef(rx, ry)[0, 1]
    n, k = len(xr), Zr.shape[1]
    t = r * np.sqrt((n - 2 - k) / (1 - r**2))
    p = 2 * stats.t.sf(abs(t), df=n - 2 - k)
    return r, p


def stage_group(s):
    if not isinstance(s, str) or "Stage" not in s:
        return np.nan
    s = s.replace("Stage ", "").rstrip("AB")
    return {"I": 0, "II": 0, "III": 1, "IV": 1}.get(s, np.nan)


def main():
    mcp = load_mcp_signatures()
    needed = set(TARGETS) | set(GEP18) | set(CYT_GENES) | set(CHECKPOINTS) | set(TJ_GENES)
    for genes in mcp.values():
        needed |= set(genes)

    purity_df = pd.read_csv(os.path.join(DATA_DIR, "tcga_absolute_purity.txt"), sep="\t")
    purity_df["patient"] = purity_df["array"].str.slice(0, 12)
    purity_df = purity_df[purity_df["array"].str.endswith("-01")]
    purity = purity_df.drop_duplicates("patient").set_index("patient")["purity"]

    cdr = pd.read_csv(os.path.join(DATA_DIR, "tcga_cdr_survival.tsv"), sep="\t")
    # one row per patient, preferring the primary-tumor (-01) sample row
    cdr["_is01"] = (~cdr["sample"].str.endswith("-01")).astype(int)
    cdr = (cdr.sort_values(["_PATIENT", "_is01", "sample"])
              .drop_duplicates("_PATIENT")
              .set_index("_PATIENT"))

    summary = {}
    tn_rows, corr_rows, quart_rows, lr_rows, cox_rows = [], [], [], [], []
    figs = {}

    for cohort in COHORTS:
        expr = load_expression(cohort, needed)
        tumors, normals = split_samples(expr.columns)
        tum = expr[list(tumors.values())]
        tum.columns = list(tumors.keys())  # patient IDs
        nor = expr[list(normals.values())]
        nor.columns = list(normals.keys())

        missing = sorted(needed - set(expr.index))
        summary[cohort] = {
            "n_tumor": tum.shape[1],
            "n_normal": nor.shape[1],
            "signature_genes_missing_from_matrix": missing,
        }

        # ---- 1. tumor vs normal -------------------------------------------
        for g in TARGETS:
            t, n = tum.loc[g], nor.loc[g]
            u, p_unp = stats.mannwhitneyu(t, n, alternative="two-sided")
            paired = sorted(set(t.index) & set(n.index))
            if len(paired) >= 10:
                w_stat, p_pair = stats.wilcoxon(t[paired], n[paired])
            else:
                p_pair = np.nan
            tn_rows.append({
                "cohort": cohort, "gene": g,
                "median_log2tpm_tumor": t.median(), "median_log2tpm_normal": n.median(),
                "delta_median": t.median() - n.median(),
                "n_tumor": len(t), "n_normal": len(n),
                "p_mannwhitney_unpaired": p_unp,
                "n_paired": len(paired), "p_wilcoxon_paired": p_pair,
            })

        # ---- 2. immune features -------------------------------------------
        feats = {}
        for pop, genes in mcp.items():
            present = [g for g in genes if g in tum.index]
            feats[pop] = tum.loc[present].mean(axis=0)
        gep = tum.loc[[g for g in GEP18 if g in tum.index]]
        feats["GEP18"] = ((gep.T - gep.mean(axis=1)) / gep.std(axis=1)).mean(axis=1)
        feats["CYT"] = tum.loc[[g for g in CYT_GENES if g in tum.index]].mean(axis=0)
        for g in CHECKPOINTS:
            if g in tum.index:
                feats[g] = tum.loc[g]
        feat_df = pd.DataFrame(feats)

        pur = purity.reindex(tum.columns).dropna()
        summary[cohort]["n_with_absolute_purity"] = len(pur)
        for g in TARGETS:
            rho_pur, p_pur = stats.spearmanr(tum.loc[g, pur.index], pur)
            summary[cohort][f"{g}_vs_purity_spearman_rho"] = round(float(rho_pur), 3)
            summary[cohort][f"{g}_vs_purity_p"] = float(p_pur)

        for g in TARGETS:
            x = tum.loc[g]
            block = []
            for fname, fvals in feat_df.items():
                rho, p = stats.spearmanr(x, fvals)
                pr, pp = partial_spearman(
                    x[pur.index].values, fvals[pur.index].values, pur.values)
                block.append({
                    "cohort": cohort, "gene": g, "immune_feature": fname,
                    "spearman_rho": rho, "p_value": p, "n": len(x),
                    "purity_adj_partial_rho": pr, "purity_adj_p": pp,
                    "n_purity_adj": len(pur),
                })
            block = pd.DataFrame(block)
            block["fdr_bh"] = bh(block["p_value"])
            block["purity_adj_fdr_bh"] = bh(block["purity_adj_p"])
            corr_rows.append(block)

            # quartile contrast Q4 vs Q1
            q1, q3 = x.quantile(0.25), x.quantile(0.75)
            lo, hi = x[x <= q1].index, x[x >= q3].index
            for f in KEY_SCORES:
                a, b = feat_df.loc[hi, f], feat_df.loc[lo, f]
                _, p = stats.mannwhitneyu(a, b, alternative="two-sided")
                quart_rows.append({
                    "cohort": cohort, "gene": g, "immune_feature": f,
                    "median_Q4_high": a.median(), "median_Q1_low": b.median(),
                    "n_Q4": len(a), "n_Q1": len(b), "p_mannwhitney": p,
                })

        # spearman between the two targets
        rho_tt, p_tt = stats.spearmanr(tum.loc["TACSTD2"], tum.loc["CLDN4"])
        summary[cohort]["TACSTD2_vs_CLDN4_spearman_rho"] = round(float(rho_tt), 3)
        summary[cohort]["TACSTD2_vs_CLDN4_p"] = float(p_tt)

        # ---- 3. survival ----------------------------------------------------
        cdr_c = cdr[cdr["cancer type abbreviation"] == cohort]
        clin = cdr_c.loc[cdr_c.index.intersection(tum.columns)].copy()
        clin["age"] = pd.to_numeric(clin["age_at_initial_pathologic_diagnosis"], errors="coerce")
        clin["male"] = (clin["gender"] == "MALE").astype(float)
        clin["stage34"] = clin["ajcc_pathologic_tumor_stage"].map(stage_group)
        summary[cohort]["n_with_cdr_clinical"] = len(clin)

        km_data = {}
        for endpoint in ["OS", "PFI"]:
            ev, tt = endpoint, f"{endpoint}.time"
            for g in TARGETS:
                df = clin[[ev, tt, "age", "male", "stage34"]].copy()
                df["expr"] = tum.loc[g, df.index].values
                df = df.dropna(subset=[ev, tt])
                df = df[df[tt] > 0]
                df["z"] = (df["expr"] - df["expr"].mean()) / df["expr"].std()
                med = df["expr"].median()
                df["high"] = (df["expr"] > med).astype(int)

                lr = logrank_test(
                    df.loc[df.high == 1, tt], df.loc[df.high == 0, tt],
                    df.loc[df.high == 1, ev], df.loc[df.high == 0, ev],
                )
                lr_rows.append({
                    "cohort": cohort, "gene": g, "endpoint": endpoint,
                    "n_high": int(df.high.sum()), "n_low": int((1 - df.high).sum()),
                    "events_high": int(df.loc[df.high == 1, ev].sum()),
                    "events_low": int(df.loc[df.high == 0, ev].sum()),
                    "logrank_p": lr.p_value,
                })
                km_data[(endpoint, g)] = df[[tt, ev, "high"]].rename(
                    columns={tt: "time", ev: "event"})

                for label, covs in [("unadjusted", ["z"]),
                                    ("adjusted", ["z", "age", "male", "stage34"])]:
                    sub = df[[tt, ev] + covs].dropna()
                    cph = CoxPHFitter()
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        cph.fit(sub, duration_col=tt, event_col=ev)
                    s = cph.summary.loc["z"]
                    cox_rows.append({
                        "cohort": cohort, "gene": g, "endpoint": endpoint,
                        "model": label, "n": len(sub),
                        "n_events": int(sub[ev].sum()),
                        "HR_per_SD": np.exp(s["coef"]),
                        "CI95_low": np.exp(s["coef lower 95%"]),
                        "CI95_high": np.exp(s["coef upper 95%"]),
                        "p_value": s["p"],
                    })

        figs[cohort] = {"tum": tum, "nor": nor, "feat": feat_df, "km": km_data,
                        "pur": pur}

    # ---- pooled LUAD+LUSC purity-partial analysis ---------------------------
    # Primary statistic: rank-based partial Spearman controlling ABSOLUTE
    # purity AND histology (LUAD/LUSC indicator). GEP18 z-scores were computed
    # within cohort, so pooling with a cohort covariate keeps them centred.
    pool_expr = pd.concat(
        [figs[c]["tum"].loc[sorted(needed & set(figs[c]["tum"].index))] for c in COHORTS],
        axis=1)
    pool_feat = pd.concat([figs[c]["feat"] for c in COHORTS], axis=0)
    pool_cohort = pd.Series(
        np.concatenate([[c] * figs[c]["tum"].shape[1] for c in COHORTS]),
        index=pool_expr.columns)
    pool_pur = pd.concat([figs[c]["pur"] for c in COHORTS])
    idx = pool_pur.index  # patients with purity available
    is_lusc = (pool_cohort[idx] == "LUSC").astype(float).values
    Z = np.column_stack([pool_pur.values, is_lusc])
    summary["POOLED"] = {
        "n_tumor_total": int(pool_expr.shape[1]),
        "n_with_purity_used_for_partial": int(len(idx)),
        "covariates": ["ABSOLUTE purity", "cohort (LUSC vs LUAD)"],
    }
    for g in TARGETS:
        x = pool_expr.loc[g]
        block = []
        for fname, fvals in pool_feat.items():
            rho, p = stats.spearmanr(x, fvals.loc[x.index])
            pr, pp = partial_spearman(x[idx].values, fvals[idx].values, Z)
            block.append({
                "cohort": "POOLED", "gene": g, "immune_feature": fname,
                "spearman_rho": rho, "p_value": p, "n": len(x),
                "purity_adj_partial_rho": pr, "purity_adj_p": pp,
                "n_purity_adj": len(idx),
            })
        block = pd.DataFrame(block)
        block["fdr_bh"] = bh(block["p_value"])
        block["purity_adj_fdr_bh"] = bh(block["purity_adj_p"])
        corr_rows.append(block)

    # ---- tight-junction genes vs TACSTD2 ------------------------------------
    tj_rows = []
    for scope in COHORTS + ["POOLED"]:
        if scope == "POOLED":
            mat, pidx, zmat = pool_expr, idx, Z
        else:
            mat = figs[scope]["tum"]
            pidx = figs[scope]["pur"].index
            zmat = figs[scope]["pur"].values
        x = mat.loc["TACSTD2"]
        block = []
        for tj in TJ_GENES:
            y = mat.loc[tj]
            rho, p = stats.spearmanr(x, y)
            pr, pp = partial_spearman(x[pidx].values, y[pidx].values, zmat)
            block.append({
                "cohort": scope, "gene": "TACSTD2", "tj_gene": tj,
                "spearman_rho": rho, "p_value": p, "n": len(x),
                "purity_adj_partial_rho": pr, "purity_adj_p": pp,
                "n_purity_adj": len(pidx),
            })
        block = pd.DataFrame(block)
        block["fdr_bh"] = bh(block["p_value"])
        block["purity_adj_fdr_bh"] = bh(block["purity_adj_p"])
        tj_rows.append(block)
    tj_df = pd.concat(tj_rows, ignore_index=True)
    tj_df.to_csv(os.path.join(OUT, "tj_correlations.csv"), index=False)

    # ---- tables -------------------------------------------------------------
    pd.DataFrame(tn_rows).to_csv(os.path.join(OUT, "tumor_vs_normal.csv"), index=False)
    corr_df = pd.concat(corr_rows, ignore_index=True)
    corr_df.to_csv(os.path.join(OUT, "immune_correlations.csv"), index=False)
    pd.DataFrame(quart_rows).to_csv(
        os.path.join(OUT, "quartile_immune_comparison.csv"), index=False)
    pd.DataFrame(lr_rows).to_csv(os.path.join(OUT, "survival_logrank.csv"), index=False)
    pd.DataFrame(cox_rows).to_csv(os.path.join(OUT, "survival_cox.csv"), index=False)
    with open(os.path.join(OUT, "cohort_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    # ---- figures --------------------------------------------------------------
    plt.rcParams.update({"font.size": 9, "figure.dpi": 150})

    # fig 1: tumor vs normal
    fig, axes = plt.subplots(1, 4, figsize=(11, 3.2), sharey=False)
    i = 0
    for cohort in COHORTS:
        for g in TARGETS:
            ax = axes[i]
            t = figs[cohort]["tum"].loc[g]
            n = figs[cohort]["nor"].loc[g]
            ax.boxplot([n, t], tick_labels=[f"N ({len(n)})", f"T ({len(t)})"],
                       showfliers=False, widths=0.55)
            for j, v in enumerate([n, t], start=1):
                ax.scatter(np.random.normal(j, 0.06, len(v)), v, s=3, alpha=0.25,
                           color="steelblue")
            p = [r for r in tn_rows if r["cohort"] == cohort and r["gene"] == g][0][
                "p_mannwhitney_unpaired"]
            ax.set_title(f"{cohort} {g}\nMWU p={p:.1e}")
            if i == 0:
                ax.set_ylabel("log2(TPM+1)")
            i += 1
    fig.suptitle("Tumor vs adjacent normal (Xena GDC STAR TPM)", y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig1_tumor_vs_normal.png"), bbox_inches="tight")
    plt.close(fig)

    # fig 2: correlation heatmaps (marginal + tumor-purity-adjusted)
    piv = corr_df.pivot_table(index="immune_feature", values="spearman_rho",
                              columns=["cohort", "gene"])
    order = piv.mean(axis=1).sort_values().index
    panels = [
        ("spearman_rho", "fdr_bh", "Marginal Spearman"),
        ("purity_adj_partial_rho", "purity_adj_fdr_bh", "Purity-adjusted (ABSOLUTE)"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    fig.subplots_adjust(wspace=0.55)
    for ax, (vcol, fcol, title) in zip(axes, panels):
        piv = corr_df.pivot_table(index="immune_feature", values=vcol,
                                  columns=["cohort", "gene"]).loc[order]
        fdr = corr_df.pivot_table(index="immune_feature", values=fcol,
                                  columns=["cohort", "gene"]).loc[order]
        im = ax.imshow(piv.values, cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")
        ax.set_xticks(range(piv.shape[1]))
        ax.set_xticklabels([f"{c}\n{g}" for c, g in piv.columns])
        ax.set_yticks(range(piv.shape[0]))
        ax.set_yticklabels(piv.index)
        for r in range(piv.shape[0]):
            for c in range(piv.shape[1]):
                star = "*" if fdr.values[r, c] < 0.05 else ""
                ax.text(c, r, f"{piv.values[r, c]:.2f}{star}", ha="center",
                        va="center", fontsize=7)
        ax.set_title(title)
    fig.colorbar(im, ax=axes, label="rho", shrink=0.5)
    fig.suptitle("Immune features vs TACSTD2 / CLDN4  (* BH-FDR < 0.05)", y=0.98)
    fig.savefig(os.path.join(OUT, "fig2_immune_correlation_heatmap.png"),
                bbox_inches="tight")
    plt.close(fig)

    # figs 3-4: KM curves
    for endpoint, fname in [("OS", "fig3_km_os.png"), ("PFI", "fig4_km_pfi.png")]:
        fig, axes = plt.subplots(2, 2, figsize=(9, 7.5))
        for r, cohort in enumerate(COHORTS):
            for c, g in enumerate(TARGETS):
                ax = axes[r][c]
                df = figs[cohort]["km"][(endpoint, g)]
                for grp, lab, col in [(1, "high (>median)", "firebrick"),
                                      (0, "low (<=median)", "steelblue")]:
                    sub = df[df.high == grp]
                    kmf = KaplanMeierFitter()
                    kmf.fit(sub["time"] / 365.25, sub["event"],
                            label=f"{lab} n={len(sub)}")
                    kmf.plot_survival_function(ax=ax, ci_show=False, color=col)
                p = [x for x in lr_rows if x["cohort"] == cohort and x["gene"] == g
                     and x["endpoint"] == endpoint][0]["logrank_p"]
                ax.set_title(f"{cohort} {g} — {endpoint}, log-rank p={p:.3f}")
                ax.set_xlabel("years")
                ax.set_ylabel(f"{endpoint} probability")
                ax.set_ylim(0, 1.02)
        fig.tight_layout()
        fig.savefig(os.path.join(OUT, fname), bbox_inches="tight")
        plt.close(fig)

    # fig 5: scatter vs cytotoxic lymphocytes
    fig, axes = plt.subplots(2, 2, figsize=(8, 7))
    for r, cohort in enumerate(COHORTS):
        for c, g in enumerate(TARGETS):
            ax = axes[r][c]
            x = figs[cohort]["tum"].loc[g]
            y = figs[cohort]["feat"]["Cytotoxic lymphocytes"]
            rho, p = stats.spearmanr(x, y)
            ax.scatter(x, y, s=5, alpha=0.35, color="darkslategray")
            ax.set_xlabel(f"{g} log2(TPM+1)")
            ax.set_ylabel("MCP-counter cytotoxic lymphocytes")
            ax.set_title(f"{cohort}: rho={rho:.2f}, p={p:.1e}")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig5_scatter_cytotoxic.png"), bbox_inches="tight")
    plt.close(fig)

    # fig 6: tight-junction genes vs TACSTD2 (marginal + purity-adjusted)
    fig, axes = plt.subplots(1, 2, figsize=(8, 4.2))
    fig.subplots_adjust(wspace=0.45)
    for ax, (vcol, fcol, title) in zip(axes, [
            ("spearman_rho", "fdr_bh", "Marginal Spearman"),
            ("purity_adj_partial_rho", "purity_adj_fdr_bh",
             "Purity(-and-cohort)-adjusted")]):
        piv = tj_df.pivot_table(index="tj_gene", values=vcol, columns="cohort")
        piv = piv[["LUAD", "LUSC", "POOLED"]].loc[TJ_GENES]
        fdr = tj_df.pivot_table(index="tj_gene", values=fcol, columns="cohort")
        fdr = fdr[["LUAD", "LUSC", "POOLED"]].loc[TJ_GENES]
        im = ax.imshow(piv.values, cmap="RdBu_r", vmin=-0.8, vmax=0.8, aspect="auto")
        ax.set_xticks(range(3))
        ax.set_xticklabels(piv.columns)
        ax.set_yticks(range(len(TJ_GENES)))
        ax.set_yticklabels(piv.index)
        for r in range(piv.shape[0]):
            for c in range(piv.shape[1]):
                star = "*" if fdr.values[r, c] < 0.05 else ""
                ax.text(c, r, f"{piv.values[r, c]:.2f}{star}", ha="center",
                        va="center", fontsize=8)
        ax.set_title(title, fontsize=9)
    fig.colorbar(im, ax=axes, label="rho", shrink=0.7)
    fig.suptitle("Tight-junction genes vs TACSTD2  (* BH-FDR < 0.05)", y=1.0)
    fig.savefig(os.path.join(OUT, "fig6_tj_correlations.png"), bbox_inches="tight")
    plt.close(fig)

    print(json.dumps(summary, indent=2))
    print("done ->", OUT)


if __name__ == "__main__":
    main()
