#!/usr/bin/env python3
"""TCGA-OV: TACSTD2 (TROP2) and CLDN4 vs immune contexture and survival.

Analog of the B1 slice (fable_tcga, TCGA LUAD/LUSC) applied to TCGA-OV
(ovarian serous cystadenocarcinoma).

Inputs (all open access, fetched by download_data.py):
  * Xena GDC hub STAR TPM matrix (log2(TPM+1), GENCODE v36)
  * TCGA-CDR curated survival endpoints + clinical covariates
  * MCP-counter marker-gene signatures (Becht et al. 2016)
  * ABSOLUTE tumor purity (PanCanAtlas open supplement)

Analyses (primary tumors only, one sample per patient):
  1. Tumor vs adjacent-normal expression: NOT POSSIBLE in TCGA-OV — the
     RNA-seq matrix contains no solid-tissue-normal (-11) samples (only
     422 primary -01 and 7 recurrent -02). Recorded as such; no synthetic
     comparison is produced.
  2. Immune deconvolution: MCP-counter population scores (mean log2 TPM of
     marker genes), T-cell-inflamed GEP (Ayers et al., JCI 2017; mean of
     18 per-gene z-scores), cytolytic activity CYT (mean log2 of GZMA, PRF1),
     plus individual checkpoint / effector genes. Spearman correlations with
     TACSTD2 and CLDN4, BH-FDR corrected within each gene family; rank-based
     partial correlations adjusting for ABSOLUTE tumor purity.
  3. Quartile contrast: Q4 vs Q1 of each target gene, Mann-Whitney U on key
     immune scores.
  4. Survival: OS and PFI (TCGA-CDR). Kaplan-Meier median split + log-rank;
     Cox PH per-SD of log2 TPM, unadjusted and adjusted for age and FIGO
     stage (III/IV vs I/II from clinical_stage; TCGA-OV has no AJCC staging
     and is all-female, so the sex covariate used in B1 is dropped).

NOTE: TCGA patients were accrued before immune-checkpoint-inhibitor
approvals and were NOT ICI treated; immune associations here must not be
read as evidence of ICI benefit.

Usage:  python3 scripts/w200/B1_OV/run_analysis.py
Outputs: results/w200/B1_OV/  (CSV tables + PNG figures)
"""

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

DATA_DIR = os.environ.get("W200_B1_OV_DATA", "/tmp/w200_b1_ov_data")
REPO = os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
)
OUT = os.path.join(REPO, "results", "w200", "B1_OV")
os.makedirs(OUT, exist_ok=True)

COHORT = "OV"
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


def load_expression(symbols_needed):
    """Return log2(TPM+1) DataFrame (rows=gene symbols, cols=samples)."""
    probemap = pd.read_csv(
        os.path.join(DATA_DIR, "gencode.v36.probemap"), sep="\t", usecols=["id", "gene"]
    )
    id2sym = dict(zip(probemap["id"], probemap["gene"]))
    path = os.path.join(DATA_DIR, f"TCGA-{COHORT}.star_tpm.tsv.gz")
    expr = pd.read_csv(path, sep="\t", index_col=0)
    expr.index = expr.index.map(lambda i: id2sym.get(i, i))
    expr = expr[expr.index.isin(symbols_needed)]
    # collapse duplicate symbols: keep the row with the highest mean expression
    expr = expr.loc[expr.mean(axis=1).groupby(expr.index).idxmax()]
    return expr


def split_samples(columns):
    """Return (tumor_samples one per patient, normal_samples, code_counts)."""
    tumors, normals, counts = {}, {}, {}
    for s in sorted(columns):
        code = s.split("-")[3][:2]
        counts[code] = counts.get(code, 0) + 1
        patient = "-".join(s.split("-")[:3])
        if code == "01" and patient not in tumors:
            tumors[patient] = s
        elif code == "11" and patient not in normals:
            normals[patient] = s
    return tumors, normals, counts


def bh(pvals):
    return multipletests(pvals, method="fdr_bh")[1]


def partial_spearman(x, y, z):
    """Spearman correlation of x and y controlling for z (rank-based partial
    correlation). Returns (rho, p) with p from a t-test on n-3 df."""
    xr, yr, zr = (stats.rankdata(v) for v in (x, y, z))
    rxy = np.corrcoef(xr, yr)[0, 1]
    rxz = np.corrcoef(xr, zr)[0, 1]
    ryz = np.corrcoef(yr, zr)[0, 1]
    r = (rxy - rxz * ryz) / np.sqrt((1 - rxz**2) * (1 - ryz**2))
    n = len(x)
    t = r * np.sqrt((n - 3) / (1 - r**2))
    p = 2 * stats.t.sf(abs(t), df=n - 3)
    return r, p


def stage_group(s):
    """FIGO clinical stage -> 1 for III/IV, 0 for I/II (OV substages end
    in A/B/C, e.g. 'Stage IIIC')."""
    if not isinstance(s, str) or "Stage" not in s:
        return np.nan
    s = s.replace("Stage ", "").rstrip("ABC")
    return {"I": 0, "II": 0, "III": 1, "IV": 1}.get(s, np.nan)


def main():
    mcp = load_mcp_signatures()
    needed = set(TARGETS) | set(GEP18) | set(CYT_GENES) | set(CHECKPOINTS)
    for genes in mcp.values():
        needed |= set(genes)

    purity_df = pd.read_csv(os.path.join(DATA_DIR, "tcga_absolute_purity.txt"), sep="\t")
    purity_df["patient"] = purity_df["array"].str.slice(0, 12)
    purity_df = purity_df[purity_df["array"].str.endswith("-01")]
    purity = purity_df.drop_duplicates("patient").set_index("patient")["purity"].dropna()

    cdr = pd.read_csv(os.path.join(DATA_DIR, "tcga_cdr_survival.tsv"), sep="\t")
    cdr["_is01"] = (~cdr["sample"].str.endswith("-01")).astype(int)
    cdr = (cdr.sort_values(["_PATIENT", "_is01", "sample"])
              .drop_duplicates("_PATIENT")
              .set_index("_PATIENT"))

    expr = load_expression(needed)
    tumors, normals, code_counts = split_samples(expr.columns)
    tum = expr[list(tumors.values())]
    tum.columns = list(tumors.keys())  # patient IDs

    missing = sorted(needed - set(expr.index))
    summary = {
        "cohort": COHORT,
        "n_tumor_primary_01": tum.shape[1],
        "sample_type_code_counts_in_matrix": code_counts,
        "n_adjacent_normal_11": len(normals),
        "tumor_vs_normal": (
            "NOT PERFORMED: TCGA-OV STAR-TPM matrix contains no solid-tissue-"
            "normal (-11) samples, so the B1 tumor-vs-normal comparison is "
            "impossible for this cohort."
        ),
        "signature_genes_missing_from_matrix": missing,
    }

    corr_rows, quart_rows, lr_rows, cox_rows = [], [], [], []

    # ---- immune features -----------------------------------------------------
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
    summary["n_with_absolute_purity"] = len(pur)
    for g in TARGETS:
        rho_pur, p_pur = stats.spearmanr(tum.loc[g, pur.index], pur)
        summary[f"{g}_vs_purity_spearman_rho"] = round(float(rho_pur), 3)
        summary[f"{g}_vs_purity_p"] = float(p_pur)

    for g in TARGETS:
        x = tum.loc[g]
        block = []
        for fname, fvals in feat_df.items():
            rho, p = stats.spearmanr(x, fvals)
            pr, pp = partial_spearman(
                x[pur.index].values, fvals[pur.index].values, pur.values)
            block.append({
                "cohort": COHORT, "gene": g, "immune_feature": fname,
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
                "cohort": COHORT, "gene": g, "immune_feature": f,
                "median_Q4_high": a.median(), "median_Q1_low": b.median(),
                "n_Q4": len(a), "n_Q1": len(b), "p_mannwhitney": p,
            })

    rho_tt, p_tt = stats.spearmanr(tum.loc["TACSTD2"], tum.loc["CLDN4"])
    summary["TACSTD2_vs_CLDN4_spearman_rho"] = round(float(rho_tt), 3)
    summary["TACSTD2_vs_CLDN4_p"] = float(p_tt)

    # ---- survival ------------------------------------------------------------
    cdr_c = cdr[cdr["cancer type abbreviation"] == COHORT]
    clin = cdr_c.loc[cdr_c.index.intersection(tum.columns)].copy()
    clin["age"] = pd.to_numeric(clin["age_at_initial_pathologic_diagnosis"], errors="coerce")
    # TCGA-OV: ajcc_pathologic_tumor_stage is entirely missing; FIGO stage
    # lives in clinical_stage. All patients are female, so no sex covariate.
    clin["stage34"] = clin["clinical_stage"].map(stage_group)
    summary["n_with_cdr_clinical"] = len(clin)
    summary["n_stage_III_IV"] = int((clin["stage34"] == 1).sum())
    summary["n_stage_I_II"] = int((clin["stage34"] == 0).sum())

    km_data = {}
    for endpoint in ["OS", "PFI"]:
        ev, tt = endpoint, f"{endpoint}.time"
        for g in TARGETS:
            df = clin[[ev, tt, "age", "stage34"]].copy()
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
                "cohort": COHORT, "gene": g, "endpoint": endpoint,
                "n_high": int(df.high.sum()), "n_low": int((1 - df.high).sum()),
                "events_high": int(df.loc[df.high == 1, ev].sum()),
                "events_low": int(df.loc[df.high == 0, ev].sum()),
                "logrank_p": lr.p_value,
            })
            km_data[(endpoint, g)] = df[[tt, ev, "high"]].rename(
                columns={tt: "time", ev: "event"})

            for label, covs in [("unadjusted", ["z"]),
                                ("adjusted", ["z", "age", "stage34"])]:
                sub = df[[tt, ev] + covs].dropna()
                cph = CoxPHFitter()
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    cph.fit(sub, duration_col=tt, event_col=ev)
                s = cph.summary.loc["z"]
                cox_rows.append({
                    "cohort": COHORT, "gene": g, "endpoint": endpoint,
                    "model": label, "n": len(sub),
                    "n_events": int(sub[ev].sum()),
                    "HR_per_SD": np.exp(s["coef"]),
                    "CI95_low": np.exp(s["coef lower 95%"]),
                    "CI95_high": np.exp(s["coef upper 95%"]),
                    "p_value": s["p"],
                })

    # ---- tables ---------------------------------------------------------------
    corr_df = pd.concat(corr_rows, ignore_index=True)
    corr_df.to_csv(os.path.join(OUT, "immune_correlations.csv"), index=False)
    pd.DataFrame(quart_rows).to_csv(
        os.path.join(OUT, "quartile_immune_comparison.csv"), index=False)
    pd.DataFrame(lr_rows).to_csv(os.path.join(OUT, "survival_logrank.csv"), index=False)
    pd.DataFrame(cox_rows).to_csv(os.path.join(OUT, "survival_cox.csv"), index=False)
    with open(os.path.join(OUT, "cohort_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    # ---- figures ----------------------------------------------------------------
    plt.rcParams.update({"font.size": 9, "figure.dpi": 150})

    # fig 1: target expression distributions (tumor-vs-normal impossible: no
    # -11 samples in TCGA-OV; shown honestly as tumor-only distributions)
    fig, axes = plt.subplots(1, 2, figsize=(7, 3.2))
    for ax, g in zip(axes, TARGETS):
        v = tum.loc[g]
        ax.hist(v, bins=40, color="steelblue", alpha=0.8)
        ax.axvline(v.median(), color="firebrick", lw=1,
                   label=f"median={v.median():.2f}")
        ax.set_xlabel(f"{g} log2(TPM+1)")
        ax.set_ylabel("primary tumors")
        ax.set_title(f"OV {g} (n={len(v)})")
        ax.legend(fontsize=7)
    fig.suptitle("TCGA-OV primary tumors — no adjacent-normal RNA-seq exists",
                 y=1.04)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig1_target_expression_distributions.png"),
                bbox_inches="tight")
    plt.close(fig)

    # fig 2: correlation heatmaps (marginal + tumor-purity-adjusted)
    piv0 = corr_df.pivot_table(index="immune_feature", values="spearman_rho",
                               columns="gene")
    order = piv0.mean(axis=1).sort_values().index
    panels = [
        ("spearman_rho", "fdr_bh", "Marginal Spearman"),
        ("purity_adj_partial_rho", "purity_adj_fdr_bh", "Purity-adjusted (ABSOLUTE)"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(8.5, 8))
    fig.subplots_adjust(wspace=0.9)
    for ax, (vcol, fcol, title) in zip(axes, panels):
        piv = corr_df.pivot_table(index="immune_feature", values=vcol,
                                  columns="gene").loc[order]
        fdr = corr_df.pivot_table(index="immune_feature", values=fcol,
                                  columns="gene").loc[order]
        im = ax.imshow(piv.values, cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")
        ax.set_xticks(range(piv.shape[1]))
        ax.set_xticklabels(piv.columns)
        ax.set_yticks(range(piv.shape[0]))
        ax.set_yticklabels(piv.index)
        for r in range(piv.shape[0]):
            for c in range(piv.shape[1]):
                star = "*" if fdr.values[r, c] < 0.05 else ""
                ax.text(c, r, f"{piv.values[r, c]:.2f}{star}", ha="center",
                        va="center", fontsize=7)
        ax.set_title(title)
    fig.colorbar(im, ax=axes, label="rho", shrink=0.5)
    fig.suptitle("TCGA-OV: immune features vs TACSTD2 / CLDN4  (* BH-FDR < 0.05)",
                 y=0.98)
    fig.savefig(os.path.join(OUT, "fig2_immune_correlation_heatmap.png"),
                bbox_inches="tight")
    plt.close(fig)

    # figs 3-4: KM curves
    for endpoint, fname in [("OS", "fig3_km_os.png"), ("PFI", "fig4_km_pfi.png")]:
        fig, axes = plt.subplots(1, 2, figsize=(9, 3.9))
        for c, g in enumerate(TARGETS):
            ax = axes[c]
            df = km_data[(endpoint, g)]
            for grp, lab, col in [(1, "high (>median)", "firebrick"),
                                  (0, "low (<=median)", "steelblue")]:
                sub = df[df.high == grp]
                kmf = KaplanMeierFitter()
                kmf.fit(sub["time"] / 365.25, sub["event"],
                        label=f"{lab} n={len(sub)}")
                kmf.plot_survival_function(ax=ax, ci_show=False, color=col)
            p = [x for x in lr_rows if x["gene"] == g
                 and x["endpoint"] == endpoint][0]["logrank_p"]
            ax.set_title(f"OV {g} — {endpoint}, log-rank p={p:.3f}")
            ax.set_xlabel("years")
            ax.set_ylabel(f"{endpoint} probability")
            ax.set_ylim(0, 1.02)
        fig.tight_layout()
        fig.savefig(os.path.join(OUT, fname), bbox_inches="tight")
        plt.close(fig)

    # fig 5: scatter vs cytotoxic lymphocytes
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.7))
    for c, g in enumerate(TARGETS):
        ax = axes[c]
        x = tum.loc[g]
        y = feat_df["Cytotoxic lymphocytes"]
        rho, p = stats.spearmanr(x, y)
        ax.scatter(x, y, s=5, alpha=0.35, color="darkslategray")
        ax.set_xlabel(f"{g} log2(TPM+1)")
        ax.set_ylabel("MCP-counter cytotoxic lymphocytes")
        ax.set_title(f"OV: rho={rho:.2f}, p={p:.1e}")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig5_scatter_cytotoxic.png"), bbox_inches="tight")
    plt.close(fig)

    print(json.dumps(summary, indent=2))
    print("done ->", OUT)


if __name__ == "__main__":
    main()
