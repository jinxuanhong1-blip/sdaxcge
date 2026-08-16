#!/usr/bin/env python3
"""Analyze TACSTD2/CLDN4 against outcome in open SCLC ICI RNA cohorts."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import urllib.parse
import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import mannwhitneyu, spearmanr


RECORD_ID = "21614171"  # SCAPeSCLC v1.3.5, pinned rather than "latest"
BASE_URL = f"https://zenodo.org/api/records/{RECORD_ID}/files"
FILES = {
    "expression": "D4. Patient-Level Gene Expression (Log2 Normalized).csv",
    "survival": "D6. Patient Survival Data, Time to Event Intervals and Censoring.csv",
}
EXPECTED_MD5 = {
    # Checked against the Zenodo record API at analysis time.
    "expression": "08cbc8f402da6aa308864d562e1f84b9",
    "survival": "f4403a55f1c1a62361ad8f85878e0052",
}
TARGETS = ("TACSTD2", "CLDN4")


def download(name: str, cache_dir: Path) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    filename = FILES[name]
    destination = cache_dir / filename
    if not destination.exists():
        url = f"{BASE_URL}/{urllib.parse.quote(filename, safe='')}/content"
        with urllib.request.urlopen(url) as response:
            destination.write_bytes(response.read())
    digest = hashlib.md5(destination.read_bytes()).hexdigest()  # nosec: integrity only
    if digest != EXPECTED_MD5[name]:
        raise RuntimeError(
            f"Checksum mismatch for {filename}: expected {EXPECTED_MD5[name]}, got {digest}"
        )
    return destination


def clean_columns(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    frame.columns = [column.lstrip("\ufeff") for column in frame.columns]
    return frame


def rank_biserial(x_positive: np.ndarray, x_negative: np.ndarray) -> float:
    """Rank-biserial correlation, positive when the positive group is higher."""
    u = mannwhitneyu(x_positive, x_negative, alternative="two-sided").statistic
    return 2 * u / (len(x_positive) * len(x_negative)) - 1


def compare_groups(
    data: pd.DataFrame,
    endpoint: str,
    positive_label: str,
    stratum: str,
    gene: str = "TACSTD2",
) -> dict[str, object]:
    selected = data.loc[data[endpoint].notna(), [endpoint, gene]]
    positive = selected.loc[selected[endpoint].eq(True), gene].to_numpy(float)
    negative = selected.loc[selected[endpoint].eq(False), gene].to_numpy(float)
    if not len(positive) or not len(negative):
        return {
            "analysis": endpoint,
            "stratum": stratum,
            "positive_label": positive_label,
            "n_positive": len(positive),
            "n_negative": len(negative),
            "median_positive": np.nan,
            "median_negative": np.nan,
            "median_difference": np.nan,
            "rank_biserial": np.nan,
            "mann_whitney_p": np.nan,
        }
    test = mannwhitneyu(positive, negative, alternative="two-sided")
    return {
        "analysis": endpoint,
        "stratum": stratum,
        "positive_label": positive_label,
        "n_positive": len(positive),
        "n_negative": len(negative),
        "median_positive": np.median(positive),
        "median_negative": np.median(negative),
        "median_difference": np.median(positive) - np.median(negative),
        "rank_biserial": rank_biserial(positive, negative),
        "mann_whitney_p": test.pvalue,
    }


def adjusted_logistic(data: pd.DataFrame, endpoint: str) -> dict[str, object]:
    selected = data.loc[data[endpoint].notna(), [endpoint, "TACSTD2", "Trial"]].copy()
    selected["TACSTD2_z"] = (
        selected["TACSTD2"] - selected["TACSTD2"].mean()
    ) / selected["TACSTD2"].std(ddof=1)
    selected["IMfirst"] = (selected["Trial"] == "IMfirst").astype(int)
    design = sm.add_constant(selected[["TACSTD2_z", "IMfirst"]])
    fit = sm.GLM(
        selected[endpoint].astype(int), design, family=sm.families.Binomial()
    ).fit()
    beta = fit.params["TACSTD2_z"]
    ci_low, ci_high = fit.conf_int().loc["TACSTD2_z"]
    return {
        "analysis": endpoint,
        "stratum": "Pooled, cohort-adjusted logistic",
        "positive_label": (
            "CR/PR" if endpoint == "objective_response" else "TTP >=12 months"
        ),
        "n_positive": int(selected[endpoint].sum()),
        "n_negative": int(selected[endpoint].eq(False).sum()),
        "median_positive": np.nan,
        "median_negative": np.nan,
        "median_difference": np.nan,
        "rank_biserial": np.nan,
        "mann_whitney_p": np.nan,
        "odds_ratio_per_sd": math.exp(beta),
        "or_95ci_low": math.exp(ci_low),
        "or_95ci_high": math.exp(ci_high),
        "logistic_p": fit.pvalues["TACSTD2_z"],
    }


def holm_adjust(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index)
    valid = values.dropna().sort_values()
    count = len(valid)
    running = 0.0
    for rank, (index, value) in enumerate(valid.items()):
        adjusted = min(1.0, (count - rank) * value)
        running = max(running, adjusted)
        result.loc[index] = running
    return result


def format_number(value: object, digits: int = 3) -> str:
    if value is None or pd.isna(value):
        return "NA"
    return f"{float(value):.{digits}g}"


def make_response_plot(data: pd.DataFrame, destination: Path) -> None:
    plot_data = data.loc[data["objective_response"].notna()].copy()
    order = [
        ("CANTABRICO", False),
        ("CANTABRICO", True),
        ("IMfirst", False),
        ("IMfirst", True),
    ]
    labels = ["CANTABRICO\nSD/PD", "CANTABRICO\nCR/PR", "IMfirst\nSD/PD", "IMfirst\nCR/PR"]
    colors = ["#9e9e9e", "#2f6f9f", "#9e9e9e", "#2f6f9f"]
    groups = [
        plot_data.loc[
            (plot_data["Trial"] == trial)
            & (plot_data["objective_response"] == response),
            "TACSTD2",
        ].to_numpy()
        for trial, response in order
    ]
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    boxes = ax.boxplot(groups, patch_artist=True, widths=0.58, showfliers=False)
    for patch, color in zip(boxes["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.65)
    rng = np.random.default_rng(20260816)
    for position, values in enumerate(groups, start=1):
        jitter = rng.uniform(-0.14, 0.14, len(values))
        ax.scatter(
            np.full(len(values), position) + jitter,
            values,
            s=24,
            c="#222222",
            alpha=0.75,
            linewidths=0,
        )
    ax.set_xticks(range(1, 5), labels)
    ax.set_ylabel("Patient-level TACSTD2 (log2 normalized)")
    ax.set_title("Pretreatment TACSTD2 by objective response")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(destination, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", type=Path, default=Path("results/w200/SCLC_open")
    )
    parser.add_argument(
        "--cache", type=Path, default=Path(".cache/sclc_open/SCAPeSCLC_v1.3.5")
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    expression_path = download("expression", args.cache)
    survival_path = download("survival", args.cache)
    expression = clean_columns(pd.read_csv(expression_path))
    survival = clean_columns(pd.read_csv(survival_path))

    # The expression file already contains the clinical fields, but independently
    # merge the pinned survival table and assert one record per patient.
    expression_genes = [gene for gene in TARGETS if gene in expression.columns]
    survival_fields = [
        "Patient ID",
        "Trial",
        "Immunotherapy",
        "RECIST",
        "TTP",
        "TTP Censored",
        "PFS",
        "PFS Censored",
        "OS",
        "OS Censored",
    ]
    # "Trial" is in the expression/clinical table, while outcome rows carry drug.
    survival_for_merge = survival[
        [field for field in survival_fields if field in survival.columns]
    ].copy()
    expression["Trial"] = np.where(
        expression["Patient ID"].str.startswith("CAN"), "CANTABRICO", "IMfirst"
    )
    data = expression[
        ["Patient ID", "Trial", "Immunotherapy", *expression_genes]
    ].merge(
        survival_for_merge.drop(columns=["Trial", "Immunotherapy"], errors="ignore"),
        on="Patient ID",
        validate="one_to_one",
    )
    data["objective_response"] = data["RECIST"].map(
        {"CR": True, "PR": True, "SD": False, "PD": False}
    )
    ttp = pd.to_numeric(data["TTP"], errors="coerce")
    ttp_censored = data["TTP Censored"].str.lower().eq("yes")
    data["long_term_benefit"] = np.where(
        ttp >= 12,
        True,
        np.where((ttp < 12) & ~ttp_censored, False, np.nan),
    )
    data["long_term_benefit"] = data["long_term_benefit"].map(
        {1.0: True, 0.0: False}
    )

    availability_rows = []
    for gene in TARGETS:
        availability_rows.append(
            {
                "gene": gene,
                "available": gene in expression.columns,
                "action": (
                    "analyzed"
                    if gene in expression.columns
                    else "skipped: absent from the ~1,800-gene GeoMx CTA/custom panel"
                ),
            }
        )
    pd.DataFrame(availability_rows).to_csv(
        args.output / "target_availability.tsv", sep="\t", index=False
    )

    if "TACSTD2" not in data:
        raise RuntimeError("TACSTD2 was unexpectedly absent from the pinned release")

    results: list[dict[str, object]] = []
    endpoint_labels = {
        "objective_response": "CR/PR",
        "long_term_benefit": "TTP >=12 months",
    }
    for endpoint, label in endpoint_labels.items():
        results.append(compare_groups(data, endpoint, label, "Pooled"))
        for trial in ("CANTABRICO", "IMfirst"):
            results.append(
                compare_groups(
                    data.loc[data["Trial"] == trial], endpoint, label, trial
                )
            )
        results.append(adjusted_logistic(data, endpoint))
    tests = pd.DataFrame(results)
    tests["holm_p_across_reported_tests"] = holm_adjust(
        tests["mann_whitney_p"].combine_first(tests["logistic_p"])
    )
    tests.to_csv(args.output / "association_tests.tsv", sep="\t", index=False)

    export_fields = [
        "Patient ID",
        "Trial",
        "Immunotherapy",
        "RECIST",
        "objective_response",
        "TTP",
        "TTP Censored",
        "long_term_benefit",
        "TACSTD2",
    ]
    data[export_fields].to_csv(
        args.output / "analysis_input_patient_level.tsv", sep="\t", index=False
    )

    spearman = spearmanr(
        data["TACSTD2"], pd.to_numeric(data["TTP"]), nan_policy="omit"
    )
    summary = {
        "record": f"https://doi.org/10.5281/zenodo.{RECORD_ID}",
        "record_version": "1.3.5",
        "n_patients": int(len(data)),
        "n_objective_response_evaluable": int(data["objective_response"].notna().sum()),
        "n_long_term_benefit_evaluable": int(data["long_term_benefit"].notna().sum()),
        "targets_requested": list(TARGETS),
        "targets_analyzed": expression_genes,
        "targets_skipped": [gene for gene in TARGETS if gene not in expression.columns],
        "tacstd2_ttp_spearman_rho": spearman.statistic,
        "tacstd2_ttp_spearman_p": spearman.pvalue,
    }
    (args.output / "summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n"
    )
    make_response_plot(data, args.output / "TACSTD2_by_RECIST_response.png")

    primary = tests.loc[
        (tests["analysis"] == "objective_response")
        & (tests["stratum"] == "Pooled")
    ].iloc[0]
    adjusted = tests.loc[
        (tests["analysis"] == "objective_response")
        & (tests["stratum"] == "Pooled, cohort-adjusted logistic")
    ].iloc[0]
    long_term = tests.loc[
        (tests["analysis"] == "long_term_benefit")
        & (tests["stratum"] == "Pooled")
    ].iloc[0]
    report = f"""# Open SCLC ICI RNA: TACSTD2/CLDN4 versus outcome

## Result

In the two analyzable open cohorts, patient-level pretreatment **TACSTD2 was not
associated with RECIST objective response**. Among {int(primary.n_positive)}
CR/PR and {int(primary.n_negative)} SD/PD patients, median log2-normalized
expression was {format_number(primary.median_positive)} versus
{format_number(primary.median_negative)} (difference
{format_number(primary.median_difference)}; rank-biserial
{format_number(primary.rank_biserial)}; two-sided Mann–Whitney
P={format_number(primary.mann_whitney_p)}). The cohort-adjusted odds ratio per
1-SD higher TACSTD2 was {format_number(adjusted.odds_ratio_per_sd)} (95% CI
{format_number(adjusted.or_95ci_low)}–{format_number(adjusted.or_95ci_high)};
P={format_number(adjusted.logistic_p)}).

**CLDN4 was not measured by the deposited GeoMx CTA/custom panel and was
skipped.** No proxy gene was substituted.

The secondary long-term-benefit comparison was also null: median TACSTD2 was
{format_number(long_term.median_positive)} in {int(long_term.n_positive)}
patients with TTP ≥12 months and {format_number(long_term.median_negative)} in
{int(long_term.n_negative)} evaluable patients without it (Mann–Whitney
P={format_number(long_term.mann_whitney_p)}). As a continuous check, Spearman
ρ for TACSTD2 versus TTP was {format_number(spearman.statistic)}
(P={format_number(spearman.pvalue)}).

## Included open cohorts

- **GSE261345 / CANTABRICO:** 26 ES-SCLC patients, pretreatment spatial RNA,
  first-line durvalumab plus platinum–etoposide.
- **GSE261348 / IMfirst:** 32 ES-SCLC patients, pretreatment spatial RNA,
  first-line atezolizumab plus platinum–etoposide.
- Expression and outcomes were taken from the pinned open harmonization
  [SCAPeSCLC v1.3.5](https://doi.org/10.5281/zenodo.{RECORD_ID}), derived from
  those GEO deposits. Analysis is at the patient level, avoiding ROI
  pseudoreplication.

## Interpretation limits

- Both studies are single-arm chemoimmunotherapy cohorts. This analysis is an
  **outcome association**, not evidence that TACSTD2 predicts benefit specific
  to PD-L1 blockade rather than chemotherapy or prognosis.
- Objective response (CR/PR versus SD/PD) was the primary endpoint. The
  prespecified study-style secondary endpoint was time to progression at least
  12 months; patients censored before 12 months were excluded from that binary
  endpoint.
- Tests are exploratory, two-sided, and unadjusted for clinical covariates
  beyond cohort in logistic models. `association_tests.tsv` also reports a
  conservative Holm adjustment across every inferential row shown.
- The assay is a targeted ~1,800-gene panel, not whole-transcriptome RNA-seq.
  Patient values are harmonized patient-level summaries of spatial ROIs.

## Controlled/request-only cohorts skipped

See `cohort_audit.tsv`. IMpower133 was deliberately excluded because its
expression and linked clinical data are EGA-controlled. Other cohorts were
skipped when sample-level expression linked to response was request-only,
controlled, or not deposited; plots and aggregate tables were not digitized.

## Reproduce

```bash
python3 scripts/analyze_sclc_open.py --output results/w200/SCLC_open
```

The script pins the open record version and verifies source-file MD5 checksums.
Search/audit cutoff: 2026-08-16.
"""
    (args.output / "README.md").write_text(report)


if __name__ == "__main__":
    main()
