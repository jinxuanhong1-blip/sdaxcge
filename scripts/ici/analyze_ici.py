#!/usr/bin/env python3
"""Reproducible target-gene analysis in open lung-cancer ICI cohorts."""

from __future__ import annotations

import gzip
from pathlib import Path
import re

from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu, spearmanr


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "results" / "ici" / "data"
OUT = ROOT / "results" / "ici"
GENES = ("TACSTD2", "CLDN4")
ENSG = {
    "TACSTD2": "ENSG00000184292",
    "CLDN4": "ENSG00000189143",
    "CD3D": "ENSG00000167286",
    "CD3E": "ENSG00000198851",
    "CD8A": "ENSG00000153563",
}


def parse_soft_samples(path: Path) -> pd.DataFrame:
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    with gzip.open(path, "rt", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current:
                    records.append(current)
                current = {"gsm": line.split(" = ", 1)[1]}
            elif current is not None and line.startswith("!Sample_title = "):
                current["title"] = line.split(" = ", 1)[1]
            elif current is not None and line.startswith("!Sample_characteristics_ch1 = "):
                value = line.split(" = ", 1)[1]
                if ": " in value:
                    key, item = value.split(": ", 1)
                    current[key.strip().lower()] = item.strip()
    if current:
        records.append(current)
    return pd.DataFrame(records)


def bh_adjust(values: pd.Series) -> pd.Series:
    p = values.astype(float).to_numpy()
    order = np.argsort(p)
    adjusted = np.empty(len(p))
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    adjusted[order] = np.minimum.accumulate(ranked[::-1])[::-1]
    return pd.Series(np.minimum(adjusted, 1.0), index=values.index)


def median_or(expression: pd.Series, positive: pd.Series) -> tuple[float, float, float, float]:
    high = expression >= expression.median()
    table = np.array(
        [
            [(high & positive).sum(), (high & ~positive).sum()],
            [(~high & positive).sum(), (~high & ~positive).sum()],
        ],
        dtype=float,
    )
    _, fisher_p = fisher_exact(table.astype(int))
    corrected = table + (0.5 if (table == 0).any() else 0)
    odds_ratio = corrected[0, 0] * corrected[1, 1] / (
        corrected[0, 1] * corrected[1, 0]
    )
    se = np.sqrt(np.sum(1 / corrected))
    low, high_ci = np.exp(np.log(odds_ratio) + np.array([-1, 1]) * 1.96 * se)
    return odds_ratio, low, high_ci, fisher_p


def response_test(
    cohort: str,
    endpoint: str,
    expression: pd.DataFrame,
    positive: pd.Series,
) -> list[dict[str, float | str | int]]:
    positive = positive.astype(bool).reindex(expression.columns)
    output = []
    for gene in GENES:
        values = expression.loc[gene].astype(float)
        case = values[positive]
        control = values[~positive]
        u_stat, p_value = mannwhitneyu(case, control, alternative="two-sided")
        odds, ci_low, ci_high, fisher_p = median_or(values, positive)
        output.append(
            {
                "cohort": cohort,
                "endpoint": endpoint,
                "gene": gene,
                "n_positive": len(case),
                "n_negative": len(control),
                "median_positive": case.median(),
                "median_negative": control.median(),
                "median_difference": case.median() - control.median(),
                "rank_biserial": 2 * u_stat / (len(case) * len(control)) - 1,
                "mann_whitney_p": p_value,
                "high_vs_low_response_or": odds,
                "or_ci95_low": ci_low,
                "or_ci95_high": ci_high,
                "fisher_p": fisher_p,
            }
        )
    return output


def load_response_cohorts() -> tuple[list[dict], dict[str, tuple[pd.DataFrame, pd.Series, str]]]:
    cohorts: dict[str, tuple[pd.DataFrame, pd.Series, str]] = {}
    rows: list[dict] = []

    counts = pd.read_csv(DATA / "GSE126044_counts.txt.gz", sep="\t", index_col=0)
    expression = np.log2(counts.div(counts.sum(axis=0), axis=1) * 1_000_000 + 0.5)
    meta = parse_soft_samples(DATA / "GSE126044_family.soft.gz")
    meta["sample"] = meta["title"].str.replace("RNA-seq_", "", regex=False)
    response = meta.set_index("sample")["patient response"].eq("responder")
    cohorts["GSE126044 anti-PD-1 response"] = (expression, response, "response")
    rows += response_test("GSE126044", "GEO responder label", expression, response)

    expression = pd.read_csv(
        DATA / "GSE166449_Raw_gene_TPM_matrix.txt.gz", sep="\t", index_col=0
    )
    meta = parse_soft_samples(DATA / "GSE166449_family.soft.gz")
    # GEO lists samples in the same order as the 22 matrix columns.
    response = pd.Series(
        meta["title"].str.contains("_Responder", case=True).to_numpy(),
        index=expression.columns,
    )
    cohorts["GSE166449 immunotherapy response"] = (expression, response, "response")
    rows += response_test("GSE166449", "GEO responder label", expression, response)

    expression = pd.read_csv(
        DATA / "GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz",
        sep="\t",
        index_col=0,
    )
    meta = pd.read_excel(DATA / "GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx")
    meta = meta[meta["Sample"].astype(str).isin(expression.columns)].set_index("Sample")
    mpr = meta["Pathologic Response"].astype(str).str.startswith("MPR")
    rows += response_test("GSE207422", "MPR vs NMPR", expression, mpr)
    orr = meta["RECIST"].isin(["CR", "PR"])
    rows += response_test("GSE207422", "RECIST ORR (CR/PR) vs SD", expression, orr)
    cohorts["GSE207422 neoadjuvant MPR"] = (expression, mpr, "MPR")
    return rows, cohorts


def load_gse135222() -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(
        DATA / "GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz", sep="\t", index_col=0
    )
    raw.index = raw.index.str.split(".").str[0]
    wanted = {symbol: ENSG[symbol] for symbol in (*GENES, "CD3D", "CD3E", "CD8A")}
    expression = pd.DataFrame(
        {symbol: np.log2(raw.loc[ensembl].astype(float) + 1) for symbol, ensembl in wanted.items()}
    ).T
    meta = parse_soft_samples(DATA / "GSE135222_family.soft.gz")
    meta["sample"] = meta["title"].str.replace(" ", "", regex=False)
    meta = meta.set_index("sample").reindex(expression.columns)
    meta["pfs_event"] = pd.to_numeric(meta["progression-free survival (pfs)"])
    meta["pfs_days"] = pd.to_numeric(meta["pfs.time"])
    return expression, meta


def survival_analysis(expression: pd.DataFrame, meta: pd.DataFrame) -> list[dict]:
    rows = []
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    for ax, gene in zip(axes, GENES):
        values = expression.loc[gene]
        z = (values - values.mean()) / values.std(ddof=1)
        model_data = pd.DataFrame(
            {"time": meta["pfs_days"], "event": meta["pfs_event"], "expression_z": z}
        ).dropna()
        model = CoxPHFitter().fit(
            model_data, duration_col="time", event_col="event", formula="expression_z"
        )
        summary = model.summary.loc["expression_z"]
        high = values >= values.median()
        km_high = KaplanMeierFitter().fit(
            meta.loc[high, "pfs_days"],
            meta.loc[high, "pfs_event"],
            label=f"High (n={high.sum()})",
        )
        km_low = KaplanMeierFitter().fit(
            meta.loc[~high, "pfs_days"],
            meta.loc[~high, "pfs_event"],
            label=f"Low (n={(~high).sum()})",
        )
        km_high.plot_survival_function(ax=ax, ci_show=False)
        km_low.plot_survival_function(ax=ax, ci_show=False)
        logrank = logrank_test(
            meta.loc[high, "pfs_days"],
            meta.loc[~high, "pfs_days"],
            event_observed_A=meta.loc[high, "pfs_event"],
            event_observed_B=meta.loc[~high, "pfs_event"],
        )
        ax.set(title=f"GSE135222: {gene}", xlabel="PFS (days)", ylabel="Probability")
        ax.text(0.98, 0.95, f"log-rank p={logrank.p_value:.3g}", ha="right", va="top", transform=ax.transAxes)
        rows.append(
            {
                "cohort": "GSE135222",
                "endpoint": "PFS",
                "gene": gene,
                "n": len(model_data),
                "events": int(model_data["event"].sum()),
                "cox_hr_per_sd": summary["exp(coef)"],
                "cox_ci95_low": summary["exp(coef) lower 95%"],
                "cox_ci95_high": summary["exp(coef) upper 95%"],
                "cox_p": summary["p"],
                "median_split_logrank_p": logrank.p_value,
            }
        )
    fig.tight_layout()
    fig.savefig(OUT / "gse135222_pfs_km.png", dpi=180)
    plt.close(fig)
    return rows


def response_figure(
    cohorts: dict[str, tuple[pd.DataFrame, pd.Series, str]], response_stats: pd.DataFrame
) -> None:
    fig, axes = plt.subplots(len(cohorts), 2, figsize=(9, 10))
    rng = np.random.default_rng(20260816)
    for row_index, (label, (expression, positive, endpoint)) in enumerate(cohorts.items()):
        positive = positive.reindex(expression.columns).astype(bool)
        for column_index, gene in enumerate(GENES):
            ax = axes[row_index, column_index]
            groups = [expression.loc[gene, ~positive], expression.loc[gene, positive]]
            ax.boxplot(groups, tick_labels=[f"No (n={len(groups[0])})", f"Yes (n={len(groups[1])})"])
            for x, values in enumerate(groups, 1):
                ax.scatter(
                    x + rng.uniform(-0.07, 0.07, len(values)),
                    values,
                    color="#276FBF",
                    s=18,
                    alpha=0.75,
                )
            cohort = label.split()[0]
            matched = response_stats[
                (response_stats["cohort"] == cohort)
                & (response_stats["gene"] == gene)
                & (
                    (response_stats["endpoint"].str.startswith(endpoint))
                    if endpoint == "MPR"
                    else True
                )
            ].iloc[0]
            ax.set_title(f"{label}\n{gene}; MW p={matched['mann_whitney_p']:.3g}")
            ax.set_ylabel("log2 expression")
            ax.set_xlabel(endpoint)
    fig.tight_layout()
    fig.savefig(OUT / "response_expression.png", dpi=180)
    plt.close(fig)


def tcell_proxy(
    response_cohorts: dict[str, tuple[pd.DataFrame, pd.Series, str]],
    gse135_expression: pd.DataFrame,
) -> pd.DataFrame:
    expressions = {name: value[0] for name, value in response_cohorts.items()}
    expressions["GSE135222 anti-PD-(L)1"] = gse135_expression
    rows = []
    fig, axes = plt.subplots(2, 2, figsize=(8, 7))
    for ax, (label, expression) in zip(axes.flat, expressions.items()):
        immune_genes = [gene for gene in ("CD3D", "CD3E", "CD8A") if gene in expression.index]
        standardized = expression.loc[immune_genes].apply(
            lambda row: (row - row.mean()) / row.std(ddof=1), axis=1
        )
        proxy = standardized.mean(axis=0)
        target = expression.loc["TACSTD2"].astype(float)
        rho, p_value = spearmanr(target, proxy)
        rows.append(
            {
                "cohort": label.split()[0],
                "n": len(target),
                "proxy_genes": ",".join(immune_genes),
                "spearman_rho": rho,
                "spearman_p": p_value,
            }
        )
        ax.scatter(target, proxy, color="#9C2C77", alpha=0.8)
        ax.set(
            title=f"{label}\nrho={rho:.2f}, p={p_value:.3g}",
            xlabel="TACSTD2 expression",
            ylabel="T-cell proxy (mean z-score)",
        )
    fig.tight_layout()
    fig.savefig(OUT / "tacstd2_tcell_proxy.png", dpi=180)
    plt.close(fig)
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    response_rows, response_cohorts = load_response_cohorts()
    response_stats = pd.DataFrame(response_rows)
    response_stats["mann_whitney_q_bh"] = bh_adjust(response_stats["mann_whitney_p"])
    response_stats["fisher_q_bh"] = bh_adjust(response_stats["fisher_p"])
    response_stats.to_csv(OUT / "response_statistics.tsv", sep="\t", index=False)
    response_figure(response_cohorts, response_stats)

    gse135_expression, gse135_meta = load_gse135222()
    survival_stats = pd.DataFrame(survival_analysis(gse135_expression, gse135_meta))
    survival_stats["cox_q_bh"] = bh_adjust(survival_stats["cox_p"])
    survival_stats.to_csv(OUT / "survival_statistics.tsv", sep="\t", index=False)

    proxy_stats = tcell_proxy(response_cohorts, gse135_expression)
    proxy_stats["spearman_q_bh"] = bh_adjust(proxy_stats["spearman_p"])
    proxy_stats.to_csv(OUT / "tcell_proxy_statistics.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
