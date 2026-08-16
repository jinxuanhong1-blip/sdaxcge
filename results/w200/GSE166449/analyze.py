#!/usr/bin/env python3
"""Test TACSTD2 and CLDN4 expression against pembrolizumab response in GSE166449."""

from __future__ import annotations

import csv
import gzip
import itertools
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
from scipy import stats


HERE = Path(__file__).resolve().parent
DATA_DIR = HERE / "data"
MATRIX_PATH = DATA_DIR / "GSE166449_Raw_gene_TPM_matrix.txt.gz"
METADATA_PATH = DATA_DIR / "GSE166449_series_matrix.txt.gz"
GENES = ("TACSTD2", "CLDN4")
RNG_SEED = 166449
N_BOOTSTRAP = 20_000


def read_metadata_row(key: str) -> list[str]:
    with gzip.open(METADATA_PATH, "rt") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if row and row[0] == key:
                return row[1:]
    raise ValueError(f"Missing metadata row: {key}")


def load_samples() -> pd.DataFrame:
    accessions = read_metadata_row("!Sample_geo_accession")
    titles = read_metadata_row("!Sample_title")
    sample_ids = read_metadata_row("!Sample_description")
    if not (len(accessions) == len(titles) == len(sample_ids) == 22):
        raise ValueError("Expected 22 aligned samples in GEO metadata")

    response = []
    for title in titles:
        if title.startswith("Immunotherapy_Responder"):
            response.append("Responder")
        elif title.startswith("Immunotherapy_nonResponder"):
            response.append("Non-responder")
        else:
            raise ValueError(f"Unrecognized response label in title: {title}")

    samples = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "geo_accession": accessions,
            "geo_title": titles,
            "response": response,
        }
    )
    samples["response_binary"] = (samples["response"] == "Responder").astype(int)
    if samples["response"].value_counts().to_dict() != {
        "Non-responder": 15,
        "Responder": 7,
    }:
        raise ValueError("Expected 7 responders and 15 non-responders")
    return samples


def load_expression(samples: pd.DataFrame) -> pd.DataFrame:
    matrix = pd.read_csv(MATRIX_PATH, sep="\t", index_col="Gene")
    if matrix.index.duplicated().any():
        raise ValueError("Gene symbols are not unique")
    if matrix.columns.tolist() != samples["sample_id"].tolist():
        raise ValueError("Expression columns do not match GEO sample metadata")
    missing = sorted(set(GENES) - set(matrix.index))
    if missing:
        raise ValueError(f"Missing target genes: {missing}")

    # The deposited values are log2(TPM + 1): low positive values map exactly
    # back to TPM increments of 0.01, despite the source filename saying TPM.
    target = matrix.loc[list(GENES)].T.reset_index(names="sample_id")
    return samples.merge(target, on="sample_id", validate="one_to_one")


def bh_adjust(p_values: list[float]) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def bootstrap_auc(
    responder: np.ndarray, non_responder: np.ndarray, rng: np.random.Generator
) -> tuple[float, float]:
    aucs = np.empty(N_BOOTSTRAP)
    for i in range(N_BOOTSTRAP):
        r = rng.choice(responder, size=len(responder), replace=True)
        nr = rng.choice(non_responder, size=len(non_responder), replace=True)
        aucs[i] = stats.mannwhitneyu(r, nr, method="asymptotic").statistic / (
            len(r) * len(nr)
        )
    return tuple(np.quantile(aucs, [0.025, 0.975]))


def hedges_g(responder: np.ndarray, non_responder: np.ndarray) -> float:
    n_r, n_nr = len(responder), len(non_responder)
    pooled_sd = math.sqrt(
        ((n_r - 1) * responder.var(ddof=1) + (n_nr - 1) * non_responder.var(ddof=1))
        / (n_r + n_nr - 2)
    )
    d = (responder.mean() - non_responder.mean()) / pooled_sd
    correction = 1 - 3 / (4 * (n_r + n_nr) - 9)
    return correction * d


def logistic_or(values: np.ndarray, outcome: np.ndarray) -> tuple[float, float, float, float]:
    standardized = (values - values.mean()) / values.std(ddof=1)
    fit = sm.GLM(
        outcome, sm.add_constant(standardized), family=sm.families.Binomial()
    ).fit()
    beta = fit.params[1]
    ci_low, ci_high = fit.conf_int()[1]
    return math.exp(beta), math.exp(ci_low), math.exp(ci_high), fit.pvalues[1]


def summarize_feature(
    data: pd.DataFrame, feature: str, family: str, rng: np.random.Generator
) -> dict[str, float | str | int]:
    responder = data.loc[data["response_binary"] == 1, feature].to_numpy(float)
    non_responder = data.loc[data["response_binary"] == 0, feature].to_numpy(float)
    welch = stats.ttest_ind(responder, non_responder, equal_var=False)
    mean_diff = responder.mean() - non_responder.mean()
    se_diff = math.sqrt(
        responder.var(ddof=1) / len(responder)
        + non_responder.var(ddof=1) / len(non_responder)
    )
    ci = stats.t.interval(0.95, df=welch.df, loc=mean_diff, scale=se_diff)
    mw = stats.mannwhitneyu(
        responder, non_responder, alternative="two-sided", method="exact"
    )
    auc = mw.statistic / (len(responder) * len(non_responder))
    auc_ci = bootstrap_auc(responder, non_responder, rng)
    odds_ratio, or_low, or_high, logistic_p = logistic_or(
        data[feature].to_numpy(float), data["response_binary"].to_numpy(int)
    )

    return {
        "feature": feature,
        "analysis_family": family,
        "n_responder": len(responder),
        "n_non_responder": len(non_responder),
        "mean_responder": responder.mean(),
        "sd_responder": responder.std(ddof=1),
        "median_responder": np.median(responder),
        "mean_non_responder": non_responder.mean(),
        "sd_non_responder": non_responder.std(ddof=1),
        "median_non_responder": np.median(non_responder),
        "mean_difference_R_minus_NR": mean_diff,
        "mean_difference_ci95_low": ci[0],
        "mean_difference_ci95_high": ci[1],
        "geometric_mean_ratio_TPM_plus_1": 2**mean_diff,
        "hedges_g": hedges_g(responder, non_responder),
        "welch_t_p_two_sided": welch.pvalue,
        "mann_whitney_u": mw.statistic,
        "mann_whitney_p_exact_two_sided": mw.pvalue,
        "auc_higher_value_predicts_response": auc,
        "auc_bootstrap_ci95_low": auc_ci[0],
        "auc_bootstrap_ci95_high": auc_ci[1],
        "odds_ratio_response_per_1sd": odds_ratio,
        "odds_ratio_ci95_low": or_low,
        "odds_ratio_ci95_high": or_high,
        "logistic_wald_p_two_sided": logistic_p,
    }


def exact_composite_permutation_p(data: pd.DataFrame) -> float:
    values = data["two_gene_z_mean"].to_numpy()
    observed = values[data["response_binary"].to_numpy() == 1].mean() - values[
        data["response_binary"].to_numpy() == 0
    ].mean()
    n = len(values)
    n_r = int(data["response_binary"].sum())
    extreme = 0
    total = 0
    all_sum = values.sum()
    for indices in itertools.combinations(range(n), n_r):
        responder_sum = values[list(indices)].sum()
        difference = responder_sum / n_r - (all_sum - responder_sum) / (n - n_r)
        extreme += abs(difference) >= abs(observed) - 1e-12
        total += 1
    return extreme / total


def make_plots(data: pd.DataFrame) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    palette = {"Non-responder": "#C44E52", "Responder": "#4C72B0"}
    order = ["Non-responder", "Responder"]

    long = data.melt(
        id_vars=["sample_id", "response"],
        value_vars=[*GENES, "two_gene_z_mean"],
        var_name="feature",
        value_name="value",
    )
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    labels = {
        "TACSTD2": "TACSTD2",
        "CLDN4": "CLDN4",
        "two_gene_z_mean": "Equal-weight two-gene score",
    }
    for ax, feature in zip(axes, labels):
        subset = long[long["feature"] == feature]
        sns.boxplot(
            data=subset,
            x="response",
            y="value",
            order=order,
            hue="response",
            palette=palette,
            legend=False,
            width=0.55,
            showfliers=False,
            ax=ax,
        )
        sns.stripplot(
            data=subset,
            x="response",
            y="value",
            order=order,
            color="black",
            jitter=0.12,
            size=6,
            ax=ax,
        )
        ax.set_title(labels[feature])
        ax.set_xlabel("")
        ax.set_ylabel("log2(TPM + 1)" if feature in GENES else "Mean within-gene z-score")
        ax.tick_params(axis="x", rotation=20)
    fig.tight_layout()
    fig.savefig(HERE / "expression_by_response.png", dpi=200, bbox_inches="tight")
    fig.savefig(HERE / "expression_by_response.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 6))
    sns.scatterplot(
        data=data,
        x="TACSTD2",
        y="CLDN4",
        hue="response",
        hue_order=order,
        palette=palette,
        s=90,
        ax=ax,
    )
    ax.set_xlabel("TACSTD2, log2(TPM + 1)")
    ax.set_ylabel("CLDN4, log2(TPM + 1)")
    ax.legend(title="")
    fig.tight_layout()
    fig.savefig(HERE / "tacstd2_cldn4_scatter.png", dpi=200, bbox_inches="tight")
    fig.savefig(HERE / "tacstd2_cldn4_scatter.pdf", bbox_inches="tight")
    plt.close(fig)


def fmt(value: float, digits: int = 3) -> str:
    return f"{value:.{digits}f}"


def write_report(
    stats_table: pd.DataFrame,
    spearman_rho: float,
    spearman_p: float,
    composite_permutation_p: float,
) -> None:
    indexed = stats_table.set_index("feature")
    tac = indexed.loc["TACSTD2"]
    cldn = indexed.loc["CLDN4"]
    composite = indexed.loc["two_gene_z_mean"]

    report = f"""# GSE166449: TACSTD2 and CLDN4 vs pembrolizumab response

## Bottom line

This small cohort provides **no persuasive evidence** that pretreatment `TACSTD2`,
`CLDN4`, or their equal-weight two-gene score is associated with pembrolizumab
response. Point estimates are weak and uncertain. If anything, `TACSTD2` is
numerically higher in responders, opposite to a simple “high TROP2 means
resistance” hypothesis, but the confidence interval includes effects in both
directions.

## Cohort and endpoint

- 22 advanced lung adenocarcinoma patients: 7 responders and 15 non-responders.
- All received pembrolizumab; all biopsies were pretreatment.
- The source paper assessed response using RECIST 1.1 and defined CR/PR as
  responder and SD/PD as non-responder.
- GEO titles provide only the binary responder label, not patient-level RECIST
  categories, follow-up, survival, treatment line, PD-L1, histology details, or
  clinical covariates.
- The deposited matrix values are `log2(TPM + 1)` (although the source filename
  says TPM). Tests use these deposited log-scale values.

## Results

| Feature | Mean R | Mean NR | Difference R−NR (95% CI) | Exact Wilcoxon/Mann–Whitney p | BH q (2 genes) | AUC, high predicts R (bootstrap 95% CI) |
|---|---:|---:|---:|---:|---:|---:|
| TACSTD2 | {fmt(tac.mean_responder)} | {fmt(tac.mean_non_responder)} | {fmt(tac.mean_difference_R_minus_NR)} ({fmt(tac.mean_difference_ci95_low)}, {fmt(tac.mean_difference_ci95_high)}) | {fmt(tac.mann_whitney_p_exact_two_sided)} | {fmt(tac.mann_whitney_bh_q_two_genes)} | {fmt(tac.auc_higher_value_predicts_response)} ({fmt(tac.auc_bootstrap_ci95_low)}, {fmt(tac.auc_bootstrap_ci95_high)}) |
| CLDN4 | {fmt(cldn.mean_responder)} | {fmt(cldn.mean_non_responder)} | {fmt(cldn.mean_difference_R_minus_NR)} ({fmt(cldn.mean_difference_ci95_low)}, {fmt(cldn.mean_difference_ci95_high)}) | {fmt(cldn.mann_whitney_p_exact_two_sided)} | {fmt(cldn.mann_whitney_bh_q_two_genes)} | {fmt(cldn.auc_higher_value_predicts_response)} ({fmt(cldn.auc_bootstrap_ci95_low)}, {fmt(cldn.auc_bootstrap_ci95_high)}) |

The equal-weight score (mean of within-cohort z-scored `TACSTD2` and `CLDN4`)
also does not separate groups: AUC {fmt(composite.auc_higher_value_predicts_response)}
(bootstrap 95% CI {fmt(composite.auc_bootstrap_ci95_low)}–{fmt(composite.auc_bootstrap_ci95_high)});
exact label-permutation p={fmt(composite_permutation_p)}. This score is
exploratory and was not externally validated.

`TACSTD2` and `CLDN4` are moderately correlated in this cohort (Spearman
rho={fmt(spearman_rho)}, p={fmt(spearman_p)}), so treating them as independent
signals would overstate the information available.

## Interpretation

- `TACSTD2`: exact two-sided p={fmt(tac.mann_whitney_p_exact_two_sided)};
  AUC={fmt(tac.auc_higher_value_predicts_response)}. This is compatible with
  chance and is not evidence of predictive utility.
- `CLDN4`: exact two-sided p={fmt(cldn.mann_whitney_p_exact_two_sided)};
  AUC={fmt(cldn.auc_higher_value_predicts_response)}, essentially random.
- Neither gene survives even a minimal two-gene multiplicity correction.
- These are treatment-outcome associations in a single-arm cohort. They cannot
  establish a treatment-specific predictive biomarker because there is no
  untreated/control arm; prognostic effects cannot be separated from
  pembrolizumab interaction effects.
- With only 7 responders, estimates and bootstrap intervals are wide. Cutpoint
  searching, multivariable fitting, or reporting an in-sample optimized model
  would be especially prone to overfitting, so none was used.
- Bulk RNA-seq mixes tumor expression with purity and cell-composition effects.
  No available covariates permit adjustment for those factors or FFPE versus
  fresh tissue.

## Files

- `sample_level_expression.csv`: response labels and both genes per sample.
- `association_statistics.csv`: descriptive statistics, effect sizes, tests,
  AUCs, and univariable odds ratios.
- `expression_by_response.*` and `tacstd2_cldn4_scatter.*`: plots.
- `analysis_manifest.txt`: input checksums and software versions.
- `analyze.py`: complete reproducible analysis.

## Sources

1. GEO [GSE166449](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166449).
2. Lee et al. *Cell* (2021), PMID 33857424,
   DOI [10.1016/j.cell.2021.03.030](https://doi.org/10.1016/j.cell.2021.03.030).
"""
    (HERE / "REPORT.md").write_text(report)


def write_manifest() -> None:
    import hashlib
    import platform

    import matplotlib
    import scipy
    import seaborn
    import statsmodels

    lines = [
        f"python={platform.python_version()}",
        f"pandas={pd.__version__}",
        f"numpy={np.__version__}",
        f"scipy={scipy.__version__}",
        f"matplotlib={matplotlib.__version__}",
        f"seaborn={seaborn.__version__}",
        f"statsmodels={statsmodels.__version__}",
        f"random_seed={RNG_SEED}",
        f"bootstrap_replicates={N_BOOTSTRAP}",
    ]
    for path in (MATRIX_PATH, METADATA_PATH):
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        lines.append(f"sha256  {digest}  data/{path.name}")
    (HERE / "analysis_manifest.txt").write_text("\n".join(lines) + "\n")


def main() -> None:
    samples = load_samples()
    data = load_expression(samples)
    for gene in GENES:
        data[f"{gene}_z"] = (data[gene] - data[gene].mean()) / data[gene].std(ddof=1)
    data["two_gene_z_mean"] = data[[f"{gene}_z" for gene in GENES]].mean(axis=1)

    rng = np.random.default_rng(RNG_SEED)
    rows = [summarize_feature(data, gene, "primary_gene", rng) for gene in GENES]
    rows.append(summarize_feature(data, "two_gene_z_mean", "exploratory_composite", rng))
    stats_table = pd.DataFrame(rows)
    primary_mask = stats_table["analysis_family"] == "primary_gene"
    stats_table.loc[primary_mask, "mann_whitney_bh_q_two_genes"] = bh_adjust(
        stats_table.loc[primary_mask, "mann_whitney_p_exact_two_sided"].tolist()
    )

    spearman = stats.spearmanr(data["TACSTD2"], data["CLDN4"])
    composite_permutation_p = exact_composite_permutation_p(data)

    output_columns = [
        "sample_id",
        "geo_accession",
        "geo_title",
        "response",
        "response_binary",
        "TACSTD2",
        "CLDN4",
        "TACSTD2_z",
        "CLDN4_z",
        "two_gene_z_mean",
    ]
    data[output_columns].to_csv(HERE / "sample_level_expression.csv", index=False)
    stats_table.to_csv(HERE / "association_statistics.csv", index=False)
    pd.DataFrame(
        [
            {
                "comparison": "TACSTD2_vs_CLDN4",
                "spearman_rho": spearman.statistic,
                "spearman_p_two_sided": spearman.pvalue,
            }
        ]
    ).to_csv(HERE / "gene_correlation.csv", index=False)

    make_plots(data)
    write_report(
        stats_table,
        spearman.statistic,
        spearman.pvalue,
        composite_permutation_p,
    )
    write_manifest()


if __name__ == "__main__":
    main()
