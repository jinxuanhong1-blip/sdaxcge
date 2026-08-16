#!/usr/bin/env python3
"""CPTAC LUAD/LSCC TACSTD2/CLDN4 protein vs immune contexture (treatment-naive)."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from scipy import stats
from statsmodels.stats.multitest import multipletests

from config import (
    CLDN4_ENSG,
    CPTAC_COHORTS,
    CPTAC_PROTEIN_SUFFIX,
    DATA,
    FIGDIR,
    IMMUNE_SCORE_COLS,
    PROTEIN_IMMUNE_MARKERS,
    RESULTS,
    TACSTD2_ENSG,
)

sns.set_theme(style="ticks", context="talk")
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["figure.dpi"] = 140


def strip_ensg(x: str) -> str:
    return str(x).split(".")[0]


def load_protein(cohort: str) -> pd.DataFrame:
    p = DATA / "cptac" / f"{cohort}{CPTAC_PROTEIN_SUFFIX}"
    df = pd.read_csv(p, sep="\t")
    df["ensg"] = df["idx"].map(strip_ensg)
    df = df.drop_duplicates("ensg").set_index("ensg")
    return df.drop(columns=["idx"])


def load_table(cohort: str, kind: str, index_col: str) -> pd.DataFrame:
    df = pd.read_csv(DATA / "cptac" / f"{cohort}_{kind}.txt", sep="\t")
    if kind == "phenotype":
        df = df.rename(columns={"idx": "case_id"})
    df["case_id"] = df["case_id"].astype(str)
    return df.set_index("case_id")


def extract_targets(prot: pd.DataFrame) -> pd.DataFrame:
    rows = {}
    for name, ensg in {"TACSTD2": TACSTD2_ENSG, "CLDN4": CLDN4_ENSG}.items():
        if ensg in prot.index:
            rows[name] = pd.to_numeric(prot.loc[ensg], errors="coerce")
        else:
            rows[name] = pd.Series(np.nan, index=prot.columns)
    for name, ensg in PROTEIN_IMMUNE_MARKERS.items():
        col = f"prot_{name}"
        if ensg in prot.index:
            rows[col] = pd.to_numeric(prot.loc[ensg], errors="coerce")
        else:
            rows[col] = pd.Series(np.nan, index=prot.columns)
    out = pd.DataFrame(rows)
    out.index.name = "case_id"
    return out


def missingness_table(joined: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cohort, sub in joined.groupby("cohort"):
        n = len(sub)
        for gene in ("TACSTD2", "CLDN4"):
            n_obs = int(sub[gene].notna().sum())
            n_na = int(sub[gene].isna().sum())
            rows.append(
                {
                    "cohort": cohort,
                    "gene": gene,
                    "n_tumor": n,
                    "n_observed": n_obs,
                    "n_NA": n_na,
                    "pct_NA": 100.0 * n_na / n if n else np.nan,
                    "min": float(sub[gene].min(skipna=True)) if n_obs else np.nan,
                    "median": float(sub[gene].median(skipna=True)) if n_obs else np.nan,
                    "max": float(sub[gene].max(skipna=True)) if n_obs else np.nan,
                }
            )
    return pd.DataFrame(rows)


def spearman_block(df: pd.DataFrame, xcols: list[str], ycols: list[str], cohort: str) -> pd.DataFrame:
    rows = []
    for x in xcols:
        for y in ycols:
            if x not in df.columns or y not in df.columns:
                continue
            pair = df[[x, y]].apply(pd.to_numeric, errors="coerce").dropna()
            n = len(pair)
            if n < 8:
                rho = p = np.nan
            else:
                rho, p = stats.spearmanr(pair[x], pair[y])
            rows.append(
                {
                    "cohort": cohort,
                    "x": x,
                    "y": y,
                    "n": n,
                    "spearman_rho": rho,
                    "p": p,
                }
            )
    out = pd.DataFrame(rows)
    if out["p"].notna().any():
        mask = out["p"].notna()
        out.loc[mask, "fdr_bh"] = multipletests(out.loc[mask, "p"], method="fdr_bh")[1]
    else:
        out["fdr_bh"] = np.nan
    return out


def mw_high_low(df: pd.DataFrame, gene: str, y: str) -> dict:
    pair = df[[gene, y]].apply(pd.to_numeric, errors="coerce").dropna()
    if len(pair) < 10:
        return {"n": len(pair), "U": np.nan, "p": np.nan, "median_high": np.nan, "median_low": np.nan}
    med = pair[gene].median()
    hi = pair.loc[pair[gene] > med, y]
    lo = pair.loc[pair[gene] <= med, y]
    if len(hi) < 3 or len(lo) < 3:
        return {"n": len(pair), "U": np.nan, "p": np.nan, "median_high": np.nan, "median_low": np.nan}
    U, p = stats.mannwhitneyu(hi, lo, alternative="two-sided")
    return {
        "n": len(pair),
        "n_high": int(len(hi)),
        "n_low": int(len(lo)),
        "U": float(U),
        "p": float(p),
        "median_high": float(hi.median()),
        "median_low": float(lo.median()),
    }


def clnd4_missing_vs_immune(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cohort, sub in df.groupby("cohort"):
        miss = sub["CLDN4"].isna()
        for y in ["ESTIMATE_ImmuneScore", "xCell_immune_score", "CIBERSORT_T_cell_CD8+", "Age"]:
            if y not in sub.columns:
                continue
            a = pd.to_numeric(sub.loc[miss, y], errors="coerce").dropna()
            b = pd.to_numeric(sub.loc[~miss, y], errors="coerce").dropna()
            if len(a) < 5 or len(b) < 5:
                p = np.nan
            else:
                p = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
            rows.append(
                {
                    "cohort": cohort,
                    "variable": y,
                    "n_CLDN4_NA": int(len(a)),
                    "n_CLDN4_obs": int(len(b)),
                    "median_NA": float(a.median()) if len(a) else np.nan,
                    "median_obs": float(b.median()) if len(b) else np.nan,
                    "mannwhitney_p": p,
                }
            )
    return pd.DataFrame(rows)


def survival_split(df: pd.DataFrame, gene: str, time_col: str, event_col: str) -> dict:
    pair = df[[gene, time_col, event_col]].apply(pd.to_numeric, errors="coerce").dropna()
    pair = pair[pair[time_col] > 0]
    n_event = int(pair[event_col].sum())
    if len(pair) < 20 or n_event < 5:
        return {
            "n": int(len(pair)),
            "n_event": n_event,
            "logrank_p": np.nan,
            "cox_hr": np.nan,
            "cox_p": np.nan,
            "note": "too_few_events_or_samples",
        }
    med = pair[gene].median()
    hi = pair[gene] > med
    lr = logrank_test(
        pair.loc[hi, time_col],
        pair.loc[~hi, time_col],
        event_observed_A=pair.loc[hi, event_col],
        event_observed_B=pair.loc[~hi, event_col],
    )
    cph_df = pair[[time_col, event_col, gene]].copy()
    cph_df[gene] = stats.zscore(cph_df[gene])
    cph = CoxPHFitter()
    try:
        cph.fit(cph_df, duration_col=time_col, event_col=event_col)
        hr = float(np.exp(cph.params_[gene]))
        cp = float(cph.summary.loc[gene, "p"])
    except Exception as e:
        hr, cp = np.nan, np.nan
        note = f"cox_failed:{e}"
    else:
        note = "ok"
    return {
        "n": int(len(pair)),
        "n_event": n_event,
        "n_high": int(hi.sum()),
        "n_low": int((~hi).sum()),
        "median_split": float(med),
        "logrank_p": float(lr.p_value),
        "cox_hr_per_SD": hr,
        "cox_p": cp,
        "note": note,
    }


def plot_missingness(miss: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    plot_df = miss.copy()
    plot_df["label"] = plot_df["cohort"] + " " + plot_df["gene"]
    colors = ["#2a6f97" if g == "TACSTD2" else "#c1121f" for g in plot_df["gene"]]
    ax.bar(plot_df["label"], plot_df["pct_NA"], color=colors)
    for i, r in plot_df.iterrows():
        ax.text(
            list(plot_df["label"]).index(r["label"]),
            r["pct_NA"] + 0.8,
            f"{int(r['n_NA'])}/{int(r['n_tumor'])}",
            ha="center",
            fontsize=11,
        )
    ax.set_ylabel("% samples with protein NA")
    ax.set_title("CPTAC tumor protein missingness")
    ax.set_ylim(0, max(35, plot_df["pct_NA"].max() + 8))
    sns.despine(ax=ax)
    fig.tight_layout()
    fig.savefig(FIGDIR / "cptac_missingness.png")
    fig.savefig(FIGDIR / "cptac_missingness.pdf")
    plt.close(fig)


def plot_scatters(joined: pd.DataFrame) -> None:
    pairs = [
        ("TACSTD2", "ESTIMATE_ImmuneScore"),
        ("CLDN4", "ESTIMATE_ImmuneScore"),
        ("TACSTD2", "xCell_immune_score"),
        ("CLDN4", "CIBERSORT_T_cell_CD8+"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 9.2))
    for ax, (x, y) in zip(axes.ravel(), pairs):
        for cohort, sub in joined.groupby("cohort"):
            pair = sub[[x, y]].apply(pd.to_numeric, errors="coerce").dropna()
            ax.scatter(pair[x], pair[y], s=22, alpha=0.7, label=f"{cohort} n={len(pair)}")
            if len(pair) >= 8:
                rho, p = stats.spearmanr(pair[x], pair[y])
                ax.text(
                    0.02,
                    0.98 if cohort == "LUAD" else 0.88,
                    f"{cohort} ρ={rho:.2f} p={p:.2e}",
                    transform=ax.transAxes,
                    va="top",
                    fontsize=10,
                )
        ax.set_xlabel(f"{x} protein (log2)")
        ax.set_ylabel(y)
        ax.legend(frameon=False, fontsize=9)
    fig.suptitle("CPTAC treatment-naive tumors (no ICI labels)", fontsize=14)
    fig.tight_layout()
    fig.savefig(FIGDIR / "cptac_scatter_immune.png")
    fig.savefig(FIGDIR / "cptac_scatter_immune.pdf")
    plt.close(fig)


def plot_corr_heatmap(corr: pd.DataFrame) -> None:
    genes = ["TACSTD2", "CLDN4"]
    for cohort in CPTAC_COHORTS:
        sub = corr[(corr["cohort"] == cohort) & (corr["x"].isin(genes))]
        if sub.empty:
            continue
        mat = sub.pivot(index="y", columns="x", values="spearman_rho")
        mat.index = [s.replace("CIBERSORT_", "CB_").replace("ESTIMATE_", "EST_").replace("xCell_", "xC_") for s in mat.index]
        fig, ax = plt.subplots(figsize=(6.4, 8.8))
        sns.heatmap(
            mat,
            ax=ax,
            cmap="RdBu_r",
            center=0,
            vmin=-0.6,
            vmax=0.6,
            annot=True,
            fmt=".2f",
            annot_kws={"size": 8},
            cbar_kws={"label": "Spearman ρ"},
        )
        ax.set_title(f"{cohort} protein vs freeze immune scores")
        fig.savefig(FIGDIR / f"cptac_{cohort}_corr_heatmap.png", bbox_inches="tight")
        fig.savefig(FIGDIR / f"cptac_{cohort}_corr_heatmap.pdf", bbox_inches="tight")
        plt.close(fig)


def plot_km(joined: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 8.6))
    axes = axes.ravel()
    i = 0
    for cohort in CPTAC_COHORTS:
        sub = joined[joined["cohort"] == cohort]
        for gene in ("TACSTD2", "CLDN4"):
            ax = axes[i]
            i += 1
            pair = sub[[gene, "OS_days", "OS_event"]].apply(pd.to_numeric, errors="coerce").dropna()
            pair = pair[pair["OS_days"] > 0]
            if len(pair) < 10:
                ax.set_title(f"{cohort} {gene} OS: insufficient data")
                continue
            med = pair[gene].median()
            km = KaplanMeierFitter()
            for lab, mask in (("high", pair[gene] > med), ("low", pair[gene] <= med)):
                if mask.sum() < 3:
                    continue
                km.fit(
                    pair.loc[mask, "OS_days"],
                    pair.loc[mask, "OS_event"],
                    label=f"{lab} n={int(mask.sum())} ev={int(pair.loc[mask,'OS_event'].sum())}",
                )
                km.plot(ax=ax, ci_show=False)
            if (pair[gene] > med).sum() >= 3 and (pair[gene] <= med).sum() >= 3:
                lr = logrank_test(
                    pair.loc[pair[gene] > med, "OS_days"],
                    pair.loc[pair[gene] <= med, "OS_days"],
                    event_observed_A=pair.loc[pair[gene] > med, "OS_event"],
                    event_observed_B=pair.loc[pair[gene] <= med, "OS_event"],
                )
                ax.set_title(f"{cohort} {gene} OS median-split logrank p={lr.p_value:.3g}")
            else:
                ax.set_title(f"{cohort} {gene} OS")
            ax.set_xlabel("OS days")
            ax.set_ylabel("Survival")
    fig.suptitle("CPTAC OS is not ICI outcome (treatment-naive resection cohort)", fontsize=13)
    fig.tight_layout()
    fig.savefig(FIGDIR / "cptac_os_km.png")
    fig.savefig(FIGDIR / "cptac_os_km.pdf")
    plt.close(fig)


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGDIR.mkdir(parents=True, exist_ok=True)
    frames = []
    miss_rows = []
    for cohort in CPTAC_COHORTS:
        prot = load_protein(cohort)
        targets = extract_targets(prot)
        ph = load_table(cohort, "phenotype", "case_id")
        sv = load_table(cohort, "survival", "case_id")
        meta = load_table(cohort, "meta", "case_id")
        # drop phenotype header-like data_type row if present
        for tbl in (ph, sv, meta):
            tbl.drop(index=[i for i in tbl.index if str(i).lower() == "data_type"], inplace=True, errors="ignore")
        j = targets.join(ph, how="left").join(sv, how="left").join(meta, how="left")
        j["cohort"] = cohort
        frames.append(j.reset_index())
        miss_rows.append(
            {
                "cohort": cohort,
                "n_protein_samples": int(prot.shape[1]),
                "n_genes": int(prot.shape[0]),
                "TACSTD2_in_matrix": TACSTD2_ENSG in prot.index,
                "CLDN4_in_matrix": CLDN4_ENSG in prot.index,
            }
        )
    joined = pd.concat(frames, ignore_index=True)
    # coerce numeric phenotype/survival
    for c in IMMUNE_SCORE_COLS + ["OS_days", "OS_event", "PFS_days", "PFS_event", "Age"]:
        if c in joined.columns:
            joined[c] = pd.to_numeric(joined[c], errors="coerce")

    miss = missingness_table(joined)
    miss.to_csv(RESULTS / "cptac_missingness.tsv", sep="\t", index=False)

    ycols = [c for c in IMMUNE_SCORE_COLS if c in joined.columns]
    prot_mark = [c for c in joined.columns if c.startswith("prot_")]
    corr_parts = []
    for cohort, sub in joined.groupby("cohort"):
        corr_parts.append(spearman_block(sub, ["TACSTD2", "CLDN4"], ycols, cohort))
        corr_parts.append(spearman_block(sub, ["TACSTD2", "CLDN4"], prot_mark, cohort))
    # combined, adjusting for cohort via residualization is overkill; report pooled Spearman too
    corr_parts.append(spearman_block(joined, ["TACSTD2", "CLDN4"], ycols, "LUAD+LSCC"))
    corr = pd.concat(corr_parts, ignore_index=True)
    corr.to_csv(RESULTS / "cptac_gene_immune_correlations.tsv", sep="\t", index=False)

    mw_rows = []
    for cohort, sub in joined.groupby("cohort"):
        for gene in ("TACSTD2", "CLDN4"):
            for y in ["ESTIMATE_ImmuneScore", "xCell_immune_score", "CIBERSORT_T_cell_CD8+"]:
                rec = mw_high_low(sub, gene, y)
                rec.update({"cohort": cohort, "gene": gene, "immune_score": y})
                mw_rows.append(rec)
    pd.DataFrame(mw_rows).to_csv(RESULTS / "cptac_highlow_mannwhitney.tsv", sep="\t", index=False)

    clnd4_miss = clnd4_missing_vs_immune(joined)
    clnd4_miss.to_csv(RESULTS / "cptac_cldn4_missingness_vs_immune.tsv", sep="\t", index=False)

    surv_rows = []
    for cohort, sub in joined.groupby("cohort"):
        for gene in ("TACSTD2", "CLDN4"):
            for time_col, event_col, endpoint in (
                ("OS_days", "OS_event", "OS"),
                ("PFS_days", "PFS_event", "PFS"),
            ):
                rec = survival_split(sub, gene, time_col, event_col)
                rec.update({"cohort": cohort, "gene": gene, "endpoint": endpoint})
                surv_rows.append(rec)
    pd.DataFrame(surv_rows).to_csv(RESULTS / "cptac_survival.tsv", sep="\t", index=False)

    # sample-level table used for stats (no extra identifiers invented)
    keep = (
        ["case_id", "cohort", "TACSTD2", "CLDN4"]
        + [c for c in ycols]
        + [c for c in ["OS_days", "OS_event", "PFS_days", "PFS_event", "Age", "Sex", "Stage"] if c in joined.columns]
        + prot_mark
    )
    joined[keep].to_csv(RESULTS / "cptac_joined_tumor_table.tsv", sep="\t", index=False)

    plot_missingness(miss)
    plot_scatters(joined)
    plot_corr_heatmap(corr[corr["y"].isin(ycols)])
    plot_km(joined)

    summary = {
        "ici_labels": "none — CPTAC LUAD/LSCC discovery cohorts are treatment-naive surgical resections",
        "protein_matrix": "gene abundance log2 reference-intensity normalized Tumor, freeze v1.2",
        "immune_scores": "same-freeze phenotype table (ESTIMATE/xCell/CIBERSORT; RNA-derived)",
        "n_LUAD": int((joined.cohort == "LUAD").sum()),
        "n_LSCC": int((joined.cohort == "LSCC").sum()),
        "missingness": miss.to_dict(orient="records"),
    }
    (RESULTS / "cptac_run_summary.json").write_text(json.dumps(summary, indent=2))
    print(miss.to_string(index=False))
    print("wrote", RESULTS)


if __name__ == "__main__":
    main()
