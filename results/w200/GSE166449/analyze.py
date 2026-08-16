#!/usr/bin/env python3
"""Test TACSTD2 and CLDN4 against pembrolizumab response in GSE166449."""

from __future__ import annotations

import csv
import gzip
import hashlib
import itertools
import json
import math
import platform
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
from scipy import stats


HERE = Path(__file__).resolve().parent
DATA_DIR = HERE / "data"
TABLE_DIR = HERE / "tables"
FIGURE_DIR = HERE / "figures"
MATRIX_PATH = DATA_DIR / "GSE166449_Raw_gene_TPM_matrix.txt.gz"
METADATA_PATH = DATA_DIR / "GSE166449_series_matrix.txt.gz"
SOFT_PATH = DATA_DIR / "GSE166449_family.soft.gz"
GENES = ("TACSTD2", "CLDN4")
TCELL_GENES = ("CD3D", "CD3E", "CD8A")
RNG_SEED = 166449
N_BOOTSTRAP = 20_000


def read_metadata_row(key: str) -> list[str]:
    with gzip.open(METADATA_PATH, "rt") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if row and row[0] == key:
                return row[1:]
    raise ValueError(f"Missing metadata row: {key}")


def parse_soft_titles() -> list[str]:
    titles: list[str] = []
    with gzip.open(SOFT_PATH, "rt", errors="replace") as handle:
        for raw in handle:
            if raw.startswith("!Sample_title = "):
                titles.append(raw.split(" = ", 1)[1].strip())
    if len(titles) != 22:
        raise ValueError(f"Expected 22 SOFT sample titles, found {len(titles)}")
    return titles


def load_samples() -> pd.DataFrame:
    accessions = read_metadata_row("!Sample_geo_accession")
    titles = read_metadata_row("!Sample_title")
    sample_ids = read_metadata_row("!Sample_description")
    if not (len(accessions) == len(titles) == len(sample_ids) == 22):
        raise ValueError("Expected 22 aligned samples in GEO metadata")
    if titles != parse_soft_titles():
        raise ValueError("Series matrix titles do not match family.soft titles")

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


def load_expression(samples: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    matrix = pd.read_csv(MATRIX_PATH, sep="\t", index_col="Gene")
    if matrix.index.duplicated().any():
        raise ValueError("Gene symbols are not unique")
    if matrix.columns.tolist() != samples["sample_id"].tolist():
        raise ValueError("Expression columns do not match GEO sample metadata")
    wanted = list(GENES + TCELL_GENES)
    missing = sorted(set(wanted) - set(matrix.index))
    if missing:
        raise ValueError(f"Missing genes: {missing}")

    # Deposited values are log2(TPM + 1): low positives map to TPM steps of 0.01.
    target = matrix.loc[wanted].T.reset_index(names="sample_id")
    data = samples.merge(target, on="sample_id", validate="one_to_one")
    return data, matrix


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


def median_split_or(values: np.ndarray, outcome: np.ndarray) -> dict[str, float | int]:
    high = values >= np.median(values)
    table = np.array(
        [
            [int(((high) & (outcome == 1)).sum()), int(((high) & (outcome == 0)).sum())],
            [int(((~high) & (outcome == 1)).sum()), int(((~high) & (outcome == 0)).sum())],
        ]
    )
    _, fisher_p = stats.fisher_exact(table)
    corrected = table + (0.5 if (table == 0).any() else 0)
    odds_ratio = (corrected[0, 0] * corrected[1, 1]) / (
        corrected[0, 1] * corrected[1, 0]
    )
    se = math.sqrt(np.sum(1 / corrected))
    low, high_ci = np.exp(np.log(odds_ratio) + np.array([-1.0, 1.0]) * 1.96 * se)
    return {
        "median_split_high_n": int(high.sum()),
        "median_split_low_n": int((~high).sum()),
        "high_responder_n": int(table[0, 0]),
        "high_non_responder_n": int(table[0, 1]),
        "low_responder_n": int(table[1, 0]),
        "low_non_responder_n": int(table[1, 1]),
        "median_split_or_high_vs_low": float(odds_ratio),
        "median_split_or_ci95_low": float(low),
        "median_split_or_ci95_high": float(high_ci),
        "median_split_fisher_p_two_sided": float(fisher_p),
    }


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
    median_or = median_split_or(
        data[feature].to_numpy(float), data["response_binary"].to_numpy(int)
    )
    return {
        "feature": feature,
        "analysis_family": family,
        "n_responder": len(responder),
        "n_non_responder": len(non_responder),
        "mean_responder": responder.mean(),
        "sd_responder": responder.std(ddof=1),
        "median_responder": float(np.median(responder)),
        "mean_non_responder": non_responder.mean(),
        "sd_non_responder": non_responder.std(ddof=1),
        "median_non_responder": float(np.median(non_responder)),
        "mean_difference_R_minus_NR": mean_diff,
        "mean_difference_ci95_low": ci[0],
        "mean_difference_ci95_high": ci[1],
        "median_difference_R_minus_NR": float(np.median(responder) - np.median(non_responder)),
        "geometric_mean_ratio_TPM_plus_1": (
            2**mean_diff if family == "primary_gene" else np.nan
        ),
        "hedges_g": hedges_g(responder, non_responder),
        "rank_biserial_R_gt_NR": 2 * auc - 1,
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
        **median_or,
    }


def exact_mean_difference_permutation_p(values: np.ndarray, outcome: np.ndarray) -> float:
    observed = values[outcome == 1].mean() - values[outcome == 0].mean()
    n = len(values)
    n_r = int(outcome.sum())
    extreme = 0
    total = 0
    all_sum = values.sum()
    for indices in itertools.combinations(range(n), n_r):
        responder_sum = values[list(indices)].sum()
        difference = responder_sum / n_r - (all_sum - responder_sum) / (n - n_r)
        extreme += abs(difference) >= abs(observed) - 1e-12
        total += 1
    return extreme / total


def add_derived_scores(data: pd.DataFrame) -> pd.DataFrame:
    for gene in GENES + TCELL_GENES:
        data[f"{gene}_z"] = (data[gene] - data[gene].mean()) / data[gene].std(ddof=1)
    data["two_gene_z_mean"] = data[[f"{gene}_z" for gene in GENES]].mean(axis=1)
    data["tcell_proxy"] = data[[f"{gene}_z" for gene in TCELL_GENES]].mean(axis=1)
    return data


def correlation_table(data: pd.DataFrame) -> pd.DataFrame:
    pairs = [
        ("TACSTD2", "CLDN4", "epithelial_coexpression"),
        ("TACSTD2", "tcell_proxy", "immune_cold_claim"),
        ("CLDN4", "tcell_proxy", "immune_cold_claim"),
        ("two_gene_z_mean", "tcell_proxy", "exploratory"),
    ]
    rows = []
    for left, right, family in pairs:
        rho, p_value = stats.spearmanr(data[left], data[right])
        rows.append(
            {
                "left": left,
                "right": right,
                "analysis_family": family,
                "n": len(data),
                "spearman_rho": rho,
                "spearman_p_two_sided": p_value,
            }
        )
    return pd.DataFrame(rows)


def make_plots(data: pd.DataFrame) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    palette = {"Non-responder": "#C44E52", "Responder": "#4C72B0"}
    order = ["Non-responder", "Responder"]
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    long = data.melt(
        id_vars=["sample_id", "response"],
        value_vars=[*GENES, "two_gene_z_mean", "tcell_proxy"],
        var_name="feature",
        value_name="value",
    )
    fig, axes = plt.subplots(1, 4, figsize=(18, 5))
    labels = {
        "TACSTD2": ("TACSTD2", "log2(TPM + 1)"),
        "CLDN4": ("CLDN4", "log2(TPM + 1)"),
        "two_gene_z_mean": ("Equal-weight two-gene score", "Mean within-gene z-score"),
        "tcell_proxy": ("T-cell proxy", "Mean z of CD3D/CD3E/CD8A"),
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
        ax.set_title(labels[feature][0])
        ax.set_xlabel("")
        ax.set_ylabel(labels[feature][1])
        ax.tick_params(axis="x", rotation=20)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "expression_by_response.png", dpi=200, bbox_inches="tight")
    fig.savefig(FIGURE_DIR / "expression_by_response.pdf", bbox_inches="tight")
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
    fig.savefig(FIGURE_DIR / "tacstd2_cldn4_scatter.png", dpi=200, bbox_inches="tight")
    fig.savefig(FIGURE_DIR / "tacstd2_cldn4_scatter.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    for ax, gene in zip(axes, GENES):
        sns.scatterplot(
            data=data,
            x=gene,
            y="tcell_proxy",
            hue="response",
            hue_order=order,
            palette=palette,
            s=90,
            ax=ax,
        )
        ax.set_xlabel(f"{gene}, log2(TPM + 1)")
        ax.set_ylabel("T-cell proxy (mean z-score)")
        ax.legend(title="")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "tcell_proxy_scatter.png", dpi=200, bbox_inches="tight")
    fig.savefig(FIGURE_DIR / "tcell_proxy_scatter.pdf", bbox_inches="tight")
    plt.close(fig)


def fmt(value: float, digits: int = 3) -> str:
    return f"{value:.{digits}f}"


def write_writeup(
    stats_table: pd.DataFrame,
    correlations: pd.DataFrame,
    composite_permutation_p: float,
    tcell_permutation_p: float,
) -> None:
    indexed = stats_table.set_index("feature")
    tac = indexed.loc["TACSTD2"]
    cldn = indexed.loc["CLDN4"]
    composite = indexed.loc["two_gene_z_mean"]
    tcell = indexed.loc["tcell_proxy"]
    corr = correlations.set_index(["left", "right"])

    writeup = f"""# GSE166449: TACSTD2 and CLDN4 vs pembrolizumab response

# GSE166449：TACSTD2 与 CLDN4 对 pembrolizumab 应答

**Verdict / 结论:** this 22-patient single-arm LUAD pembrolizumab cohort provides
**no persuasive evidence** that pretreatment `TACSTD2`, `CLDN4`, their
equal-weight two-gene score, or a bulk T-cell proxy is associated with
response. `TACSTD2` is numerically **higher** in responders, opposite a simple
high-TROP2 resistance story, but the interval includes both directions.

**结论：** 这个 22 例单臂肺腺癌 pembrolizumab 队列**不能提供有说服力的证据**
表明治疗前 `TACSTD2`、`CLDN4`、等权双基因评分或 bulk T 细胞代理指标与应答
相关。`TACSTD2` 在应答者中数值上**更高**，与“高 TROP2 即耐药”的简单假说
相反，但置信区间同时包含两个方向。

---

## 1. Cohort / 队列

| Item | Value |
|------|--------|
| GEO | [GSE166449](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166449) |
| Paper | Lee et al., *Cell* 2021, PMID [33857424](https://pubmed.ncbi.nlm.nih.gov/33857424/), DOI [10.1016/j.cell.2021.03.030](https://doi.org/10.1016/j.cell.2021.03.030) |
| Disease | advanced lung adenocarcinoma |
| Treatment | pembrolizumab (paper); GEO titles say only “immunotherapy” |
| Biopsy | pretreatment |
| Labels | 7 `Immunotherapy_Responder` vs 15 `Immunotherapy_nonResponder` |
| Response rule in the paper | RECIST 1.1; CR/PR = responder, SD/PD = non-responder |
| Public covariates | none (no PD-L1, line, survival, or patient-level RECIST category) |

GEO itself does not name the drug. The pembrolizumab assignment comes from the
Cell 2021 STAR Methods description of the Samsung Medical Center cohort.

GEO 本身未写明药物。pembrolizumab 来自 Cell 2021 STAR Methods 对三星医疗中心
队列的描述。

## 2. Methods (pre-specified) / 方法（预先指定）

- **Scale.** Use the deposited matrix as `log2(TPM + 1)`. Do not re-log.
- **Primary test.** Exact two-sided Mann–Whitney U for `TACSTD2` and `CLDN4`.
- **Effect size.** Mean difference with Welch 95% CI; Hedges’ g; AUC =
  P(R > NR); rank-biserial = 2·AUC − 1; univariable OR per 1 SD.
- **Secondary, descriptive.** Median-split Fisher exact OR (high vs low).
  This loses information and is not used to claim a cutoff.
- **Exploratory score.** Mean of within-cohort z-scored `TACSTD2` and `CLDN4`.
  Exact label-permutation p for the mean difference.
- **Context, not a hunt.** Bulk T-cell proxy = mean z of `CD3D`, `CD3E`,
  `CD8A`. Spearman vs each gene; same response tests as a power/context check.
- **Multiplicity.** BH q across the two primary gene tests only.
- **Not done.** Cutpoint search, multivariable fitting, gene-set fishing, or
  pooling with other GEO series.

未因结果更接近预期假说而改检验、改分组或改基因。

## 3. Results / 结果

### 3.1 Primary genes vs response

| Feature | Median R | Median NR | Mean R−NR (95% CI) | Exact MW p | BH q | AUC (bootstrap 95% CI) | High/low OR (95% CI) | Fisher p |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| TACSTD2 | {fmt(tac.median_responder)} | {fmt(tac.median_non_responder)} | {fmt(tac.mean_difference_R_minus_NR)} ({fmt(tac.mean_difference_ci95_low)}, {fmt(tac.mean_difference_ci95_high)}) | {fmt(tac.mann_whitney_p_exact_two_sided)} | {fmt(tac.mann_whitney_bh_q_two_genes)} | {fmt(tac.auc_higher_value_predicts_response)} ({fmt(tac.auc_bootstrap_ci95_low)}, {fmt(tac.auc_bootstrap_ci95_high)}) | {fmt(tac.median_split_or_high_vs_low)} ({fmt(tac.median_split_or_ci95_low)}, {fmt(tac.median_split_or_ci95_high)}) | {fmt(tac.median_split_fisher_p_two_sided)} |
| CLDN4 | {fmt(cldn.median_responder)} | {fmt(cldn.median_non_responder)} | {fmt(cldn.mean_difference_R_minus_NR)} ({fmt(cldn.mean_difference_ci95_low)}, {fmt(cldn.mean_difference_ci95_high)}) | {fmt(cldn.mann_whitney_p_exact_two_sided)} | {fmt(cldn.mann_whitney_bh_q_two_genes)} | {fmt(cldn.auc_higher_value_predicts_response)} ({fmt(cldn.auc_bootstrap_ci95_low)}, {fmt(cldn.auc_bootstrap_ci95_high)}) | {fmt(cldn.median_split_or_high_vs_low)} ({fmt(cldn.median_split_or_ci95_low)}, {fmt(cldn.median_split_or_ci95_high)}) | {fmt(cldn.median_split_fisher_p_two_sided)} |

`TACSTD2` median difference R−NR = {fmt(tac.median_difference_R_minus_NR)}.
`CLDN4` median difference R−NR = {fmt(cldn.median_difference_R_minus_NR)}.

### 3.2 Exploratory two-gene score

Equal-weight z-mean: AUC {fmt(composite.auc_higher_value_predicts_response)}
(bootstrap 95% CI {fmt(composite.auc_bootstrap_ci95_low)}–{fmt(composite.auc_bootstrap_ci95_high)});
exact label-permutation p={fmt(composite_permutation_p)}. Not validated.

### 3.3 T-cell proxy and co-expression

The T-cell proxy does not separate response groups: AUC
{fmt(tcell.auc_higher_value_predicts_response)}
(bootstrap 95% CI {fmt(tcell.auc_bootstrap_ci95_low)}–{fmt(tcell.auc_bootstrap_ci95_high)});
exact MW p={fmt(tcell.mann_whitney_p_exact_two_sided)};
exact permutation p={fmt(tcell_permutation_p)}.
This is a context check: the cohort is too small to recover even a conventional
immune-infiltration marker.

| Pair | Spearman rho | p |
|---|---:|---:|
| TACSTD2 vs CLDN4 | {fmt(corr.loc[("TACSTD2", "CLDN4")].spearman_rho)} | {fmt(corr.loc[("TACSTD2", "CLDN4")].spearman_p_two_sided)} |
| TACSTD2 vs T-cell proxy | {fmt(corr.loc[("TACSTD2", "tcell_proxy")].spearman_rho)} | {fmt(corr.loc[("TACSTD2", "tcell_proxy")].spearman_p_two_sided)} |
| CLDN4 vs T-cell proxy | {fmt(corr.loc[("CLDN4", "tcell_proxy")].spearman_rho)} | {fmt(corr.loc[("CLDN4", "tcell_proxy")].spearman_p_two_sided)} |

`TACSTD2` and `CLDN4` are moderately correlated, so they are not independent
signals. Neither gene shows an immune-cold correlation here (`rho ≈ 0`).

## 4. Honest reading / 诚实解读

1. **No association.** Neither primary gene, the two-gene score, nor the T-cell
   proxy reaches even an unadjusted 0.05 threshold.
2. **Direction.** `TACSTD2` is higher in responders. That is the opposite of
   “high TROP2 predicts pembrolizumab resistance” and also opposite the
   GSE126044 / GSE207422 direction reported elsewhere in this repository.
   `CLDN4` is essentially random and does not support a CLDN4-high → worse
   response claim in this cohort.
3. **Not a predictive biomarker test.** There is no control arm, so prognostic
   and treatment-interaction effects cannot be separated.
4. **Power.** With 7 vs 15, only large effects are detectable. A null p-value
   is not proof of no effect. The wide AUC intervals (crossing 0.5) are the
   honest result.
5. **Bulk RNA-seq.** Expression mixes tumor, stroma, and purity. No covariates
   are available for adjustment. The T-cell proxy is not a cell count.

## 5. Files / 文件

- `tables/sample_level_expression.csv`
- `tables/association_statistics.csv`
- `tables/gene_correlation.csv`
- `figures/expression_by_response.*`
- `figures/tacstd2_cldn4_scatter.*`
- `figures/tcell_proxy_scatter.*`
- `summary.json`
- `file_manifest.json`
- `analyze.py`, `download.py`

Reproduce:

```bash
python3 results/w200/GSE166449/download.py
python3 results/w200/GSE166449/analyze.py
# or: python3 scripts/w200/GSE166449.py
```
"""
    (HERE / "WRITEUP.md").write_text(writeup)
    (HERE / "REPORT.md").write_text(
        "# See WRITEUP.md\n\n"
        "The bilingual honest report is `WRITEUP.md`. "
        "This file is kept as a pointer for the first-pass path.\n"
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_manifest() -> None:
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
    for path in (MATRIX_PATH, METADATA_PATH, SOFT_PATH):
        lines.append(f"sha256  {sha256(path)}  data/{path.name}")
    (HERE / "analysis_manifest.txt").write_text("\n".join(lines) + "\n")


def write_summary(
    stats_table: pd.DataFrame,
    correlations: pd.DataFrame,
    composite_permutation_p: float,
    tcell_permutation_p: float,
) -> None:
    indexed = stats_table.set_index("feature")
    corr = correlations.set_index(["left", "right"])

    def gene_block(name: str) -> dict:
        row = indexed.loc[name]
        return {
            "median_responder": row["median_responder"],
            "median_non_responder": row["median_non_responder"],
            "median_difference_R_minus_NR": row["median_difference_R_minus_NR"],
            "mean_difference_R_minus_NR": row["mean_difference_R_minus_NR"],
            "mean_difference_ci95": [
                row["mean_difference_ci95_low"],
                row["mean_difference_ci95_high"],
            ],
            "hedges_g": row["hedges_g"],
            "mann_whitney_p_exact_two_sided": row["mann_whitney_p_exact_two_sided"],
            "mann_whitney_bh_q_two_genes": row.get("mann_whitney_bh_q_two_genes"),
            "auc_higher_expression_predicts_response": row[
                "auc_higher_value_predicts_response"
            ],
            "auc_bootstrap_ci95": [
                row["auc_bootstrap_ci95_low"],
                row["auc_bootstrap_ci95_high"],
            ],
            "median_split_or_high_vs_low": row["median_split_or_high_vs_low"],
            "median_split_or_ci95": [
                row["median_split_or_ci95_low"],
                row["median_split_or_ci95_high"],
            ],
            "median_split_fisher_p_two_sided": row["median_split_fisher_p_two_sided"],
        }

    summary = {
        "task": "GSE166449_TACSTD2_CLDN4_vs_pembrolizumab_response",
        "cohort": {
            "disease": "advanced lung adenocarcinoma",
            "treatment_from_paper": "pembrolizumab",
            "treatment_from_geo": "immunotherapy",
            "biopsy_timing": "pretreatment",
            "response_definition": "RECIST 1.1; CR/PR responder, SD/PD non-responder",
            "n_total": 22,
            "n_responder": 7,
            "n_non_responder": 15,
        },
        "expression_scale": "deposited log2(TPM + 1); no additional transform",
        "primary_gene_results": {gene: gene_block(gene) for gene in GENES},
        "exploratory_equal_weight_two_gene_score": {
            "auc_higher_score_predicts_response": indexed.loc[
                "two_gene_z_mean", "auc_higher_value_predicts_response"
            ],
            "auc_bootstrap_ci95": [
                indexed.loc["two_gene_z_mean", "auc_bootstrap_ci95_low"],
                indexed.loc["two_gene_z_mean", "auc_bootstrap_ci95_high"],
            ],
            "exact_label_permutation_p_two_sided": composite_permutation_p,
        },
        "tcell_proxy": {
            "genes": list(TCELL_GENES),
            "auc_higher_score_predicts_response": indexed.loc[
                "tcell_proxy", "auc_higher_value_predicts_response"
            ],
            "mann_whitney_p_exact_two_sided": indexed.loc[
                "tcell_proxy", "mann_whitney_p_exact_two_sided"
            ],
            "exact_label_permutation_p_two_sided": tcell_permutation_p,
            "spearman_vs_TACSTD2": {
                "rho": corr.loc[("TACSTD2", "tcell_proxy")].spearman_rho,
                "p_two_sided": corr.loc[("TACSTD2", "tcell_proxy")].spearman_p_two_sided,
            },
            "spearman_vs_CLDN4": {
                "rho": corr.loc[("CLDN4", "tcell_proxy")].spearman_rho,
                "p_two_sided": corr.loc[("CLDN4", "tcell_proxy")].spearman_p_two_sided,
            },
        },
        "gene_correlation": {
            "TACSTD2_vs_CLDN4_spearman_rho": corr.loc[
                ("TACSTD2", "CLDN4")
            ].spearman_rho,
            "p_two_sided_unadjusted": corr.loc[
                ("TACSTD2", "CLDN4")
            ].spearman_p_two_sided,
            "status": "descriptive_exploratory",
        },
        "honest_verdict": {
            "supports_TACSTD2_response_association": False,
            "supports_CLDN4_response_association": False,
            "supports_two_gene_response_association": False,
            "supports_TACSTD2_immune_cold_in_this_cohort": False,
            "supports_CLDN4_high_worse_response_in_this_cohort": False,
            "supports_treatment_specific_predictive_biomarker": False,
            "statement": (
                "No persuasive association was detected for TACSTD2, CLDN4, the "
                "equal-weight score, or the CD3D/CD3E/CD8A T-cell proxy. TACSTD2 "
                "was numerically higher in responders, opposite a simple "
                "high-TROP2 resistance hypothesis, with wide uncertainty. "
                "TACSTD2 vs T-cell proxy rho is near zero. A single-arm cohort "
                "cannot estimate a biomarker-by-treatment interaction."
            ),
            "limitations": [
                "7 responders and 15 non-responders",
                "single-arm observational treatment-outcome association",
                "GEO names immunotherapy; pembrolizumab is from the paper",
                "binary labels only; no patient-level RECIST category",
                "no public clinical covariates for adjustment",
                "bulk RNA-seq is sensitive to tumor purity and cell composition",
                "T-cell proxy is not a measured cell count",
                "exploratory two-gene score has no external validation",
                "median-split OR is descriptive and loses information",
            ],
        },
    }
    (HERE / "summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n"
    )


def cleanup_legacy_root_outputs() -> None:
    for name in (
        "sample_level_expression.csv",
        "association_statistics.csv",
        "gene_correlation.csv",
        "expression_by_response.png",
        "expression_by_response.pdf",
        "tacstd2_cldn4_scatter.png",
        "tacstd2_cldn4_scatter.pdf",
    ):
        path = HERE / name
        if path.exists():
            path.unlink()


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    samples = load_samples()
    data, _matrix = load_expression(samples)
    data = add_derived_scores(data)

    rng = np.random.default_rng(RNG_SEED)
    rows = [summarize_feature(data, gene, "primary_gene", rng) for gene in GENES]
    rows.append(summarize_feature(data, "two_gene_z_mean", "exploratory_composite", rng))
    rows.append(summarize_feature(data, "tcell_proxy", "context_tcell_proxy", rng))
    stats_table = pd.DataFrame(rows)
    primary_mask = stats_table["analysis_family"] == "primary_gene"
    stats_table.loc[primary_mask, "mann_whitney_bh_q_two_genes"] = bh_adjust(
        stats_table.loc[primary_mask, "mann_whitney_p_exact_two_sided"].tolist()
    )

    correlations = correlation_table(data)
    composite_permutation_p = exact_mean_difference_permutation_p(
        data["two_gene_z_mean"].to_numpy(), data["response_binary"].to_numpy()
    )
    tcell_permutation_p = exact_mean_difference_permutation_p(
        data["tcell_proxy"].to_numpy(), data["response_binary"].to_numpy()
    )

    output_columns = [
        "sample_id",
        "geo_accession",
        "geo_title",
        "response",
        "response_binary",
        "TACSTD2",
        "CLDN4",
        "CD3D",
        "CD3E",
        "CD8A",
        "TACSTD2_z",
        "CLDN4_z",
        "two_gene_z_mean",
        "tcell_proxy",
    ]
    data[output_columns].to_csv(TABLE_DIR / "sample_level_expression.csv", index=False)
    stats_table.to_csv(TABLE_DIR / "association_statistics.csv", index=False)
    correlations.to_csv(TABLE_DIR / "gene_correlation.csv", index=False)

    make_plots(data)
    write_writeup(
        stats_table,
        correlations,
        composite_permutation_p,
        tcell_permutation_p,
    )
    write_summary(
        stats_table,
        correlations,
        composite_permutation_p,
        tcell_permutation_p,
    )
    write_manifest()
    cleanup_legacy_root_outputs()


if __name__ == "__main__":
    main()
