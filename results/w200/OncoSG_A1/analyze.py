#!/usr/bin/env python3
"""Reproduce TACSTD2–immune associations in OncoSG LUAD."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import tempfile
import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests


STUDY_ID = "luad_oncosg_2020"
DATAHUB_REVISION = "165bd77077b03038f9c2ee104959eb474770b2a9"
BASE_URL = (
    f"https://media.githubusercontent.com/media/cBioPortal/datahub/{DATAHUB_REVISION}/"
    f"public/{STUDY_ID}"
)
EXPRESSION_URL = (
    f"{BASE_URL}/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt"
)
CLINICAL_URL = f"{BASE_URL}/data_clinical_sample.txt"
IMMUNE_COLUMNS = {
    "IMSIG_B_CELLS": "B cells",
    "IMSIG_INTERFERON": "Interferon",
    "IMSIG_MACROPHAGES": "Macrophages",
    "IMSIG_MONOCYTES": "Monocytes",
    "IMSIG_NEUTROPHILS": "Neutrophils",
    "IMSIG_NK_CELLS": "NK cells",
    "IMSIG_PLASMA_CELLS": "Plasma cells",
    "IMSIG_T_CELLS": "T cells",
}
CONTROL_COLUMNS = {
    "IMSIG_PROLIFERATION": "Proliferation",
    "IMSIG_TRANSLATION": "Translation",
}


def download(url: str, destination: Path) -> None:
    """Download a public source file atomically."""
    partial = destination.with_suffix(destination.suffix + ".part")
    with urllib.request.urlopen(url, timeout=120) as response, partial.open("wb") as out:
        while chunk := response.read(1024 * 1024):
            out.write(chunk)
    partial.replace(destination)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def load_tacstd2(expression_path: Path) -> pd.Series:
    """Read only the TACSTD2 row from the public gene-by-sample matrix."""
    for chunk in pd.read_csv(expression_path, sep="\t", chunksize=2_000):
        row = chunk.loc[chunk["Hugo_Symbol"].eq("TACSTD2")]
        if len(row) == 1:
            values = pd.to_numeric(
                row.drop(columns=["Hugo_Symbol", "Entrez_Gene_Id"]).iloc[0],
                errors="coerce",
            )
            values.name = "TACSTD2_Z"
            values.index.name = "SAMPLE_ID"
            return values
        if len(row) > 1:
            raise ValueError("TACSTD2 occurs more than once in the expression matrix")
    raise ValueError("TACSTD2 is absent from the expression matrix")


def load_clinical(clinical_path: Path) -> pd.DataFrame:
    clinical = pd.read_csv(
        clinical_path,
        sep="\t",
        comment="#",
        na_values=["NA", "NaN", ""],
        keep_default_na=True,
    ).set_index("SAMPLE_ID")
    columns = ["PURITY", *IMMUNE_COLUMNS, *CONTROL_COLUMNS]
    return clinical[columns].apply(pd.to_numeric, errors="coerce")


def residualize_rank(values: pd.Series, covariate: pd.Series) -> np.ndarray:
    ranked_values = stats.rankdata(values, method="average")
    ranked_covariate = stats.rankdata(covariate, method="average")
    design = np.column_stack([np.ones(len(values)), ranked_covariate])
    fitted = design @ np.linalg.lstsq(design, ranked_values, rcond=None)[0]
    return ranked_values - fitted


def partial_spearman(frame: pd.DataFrame, outcome: str) -> tuple[float, float, float, float]:
    """Partial Spearman rho and Fisher 95% CI, controlling for purity."""
    x_resid = residualize_rank(frame["TACSTD2_Z"], frame["PURITY"])
    y_resid = residualize_rank(frame[outcome], frame["PURITY"])
    rho = float(stats.pearsonr(x_resid, y_resid).statistic)
    df = len(frame) - 3
    if abs(rho) == 1:
        p_value = 0.0
    else:
        statistic = rho * math.sqrt(df / (1 - rho**2))
        p_value = float(2 * stats.t.sf(abs(statistic), df))
    fisher_se = 1 / math.sqrt(len(frame) - 4)
    fisher_z = np.arctanh(np.clip(rho, -0.999999, 0.999999))
    ci_low, ci_high = np.tanh(
        [fisher_z - 1.96 * fisher_se, fisher_z + 1.96 * fisher_se]
    )
    return rho, p_value, float(ci_low), float(ci_high)


def calculate_associations(merged: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, float | int | str | bool]] = []
    all_columns = {**IMMUNE_COLUMNS, **CONTROL_COLUMNS}
    for column, label in all_columns.items():
        complete = merged[["TACSTD2_Z", "PURITY", column]].dropna()
        raw = stats.spearmanr(complete["TACSTD2_Z"], complete[column])
        partial_rho, partial_p, ci_low, ci_high = partial_spearman(complete, column)
        rows.append(
            {
                "signature": label,
                "source_column": column,
                "is_immune": column in IMMUNE_COLUMNS,
                "n": len(complete),
                "raw_spearman_rho": raw.statistic,
                "raw_p_value": raw.pvalue,
                "partial_spearman_rho": partial_rho,
                "partial_ci95_low": ci_low,
                "partial_ci95_high": ci_high,
                "partial_p_value": partial_p,
            }
        )
    results = pd.DataFrame(rows)
    for prefix in ("raw", "partial"):
        results[f"{prefix}_q_value"] = np.nan
        immune = results["is_immune"]
        results.loc[immune, f"{prefix}_q_value"] = multipletests(
            results.loc[immune, f"{prefix}_p_value"], method="fdr_bh"
        )[1]
    return results.sort_values(
        ["is_immune", "partial_spearman_rho"], ascending=[False, False]
    )


def plot_forest(results: pd.DataFrame, output: Path) -> None:
    plot_data = results.loc[results["is_immune"]].sort_values(
        "partial_spearman_rho"
    )
    y = np.arange(len(plot_data))
    rho = plot_data["partial_spearman_rho"].to_numpy()
    errors = np.vstack(
        [
            rho - plot_data["partial_ci95_low"].to_numpy(),
            plot_data["partial_ci95_high"].to_numpy() - rho,
        ]
    )
    significant = plot_data["partial_q_value"].lt(0.05).to_numpy()
    colors = np.where(significant, "#0072B2", "#888888")
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for index in range(len(plot_data)):
        ax.errorbar(
            rho[index],
            y[index],
            xerr=errors[:, index, None],
            fmt="o",
            color=colors[index],
            capsize=3,
        )
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_yticks(y, plot_data["signature"])
    ax.set_xlabel("Partial Spearman rho (adjusted for purity)")
    ax.set_title("OncoSG LUAD: TACSTD2 vs IMSIG immune scores")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_residuals(merged: pd.DataFrame, results: pd.DataFrame, output: Path) -> None:
    ordered = results.loc[results["is_immune"]].sort_values(
        "partial_spearman_rho", ascending=False
    )
    fig, axes = plt.subplots(2, 4, figsize=(13, 6.8))
    for ax, row in zip(axes.flat, ordered.itertuples(index=False), strict=True):
        complete = merged[["TACSTD2_Z", "PURITY", row.source_column]].dropna()
        x_resid = residualize_rank(complete["TACSTD2_Z"], complete["PURITY"])
        y_resid = residualize_rank(complete[row.source_column], complete["PURITY"])
        sns.regplot(
            x=x_resid,
            y=y_resid,
            ax=ax,
            scatter_kws={"s": 18, "alpha": 0.55, "edgecolor": "none"},
            line_kws={"color": "#D55E00", "linewidth": 1.5},
            ci=None,
        )
        ax.set_title(
            f"{row.signature}\n"
            rf"$\rho_p$={row.partial_spearman_rho:.2f}, "
            f"q={row.partial_q_value:.3g}"
        )
        ax.set_xlabel("TACSTD2 rank residual")
        ax.set_ylabel("IMSIG rank residual")
    fig.suptitle("Purity-adjusted rank residuals", y=1.01, fontsize=14)
    fig.tight_layout()
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def write_report(
    output_dir: Path,
    merged: pd.DataFrame,
    results: pd.DataFrame,
    purity_result: stats.SignificanceResult,
    provenance: dict,
) -> None:
    immune = results.loc[results["is_immune"]]
    significant = immune.loc[immune["partial_q_value"].lt(0.05)]
    strongest = immune.iloc[immune["partial_spearman_rho"].abs().argmax()]
    direction = "positive" if strongest["partial_spearman_rho"] > 0 else "negative"
    significant_text = (
        ", ".join(significant["signature"]) if len(significant) else "none"
    )
    report = f"""# OncoSG LUAD: TACSTD2 versus immune scores after purity adjustment

## Result

The public cBioPortal/Datahub matrix was available, so the cohort was analyzed rather
than skipped. There were {len(merged)} expression samples; {int(merged["PURITY"].notna().sum())}
had purity and {int(merged[list(IMMUNE_COLUMNS)].notna().all(axis=1).sum())} had all eight
IMSIG immune scores. Complete-case counts for each test are reported in
`associations.tsv`.

The strongest purity-adjusted association was {strongest["signature"]}
(partial Spearman rho = {strongest["partial_spearman_rho"]:.3f},
95% CI {strongest["partial_ci95_low"]:.3f} to {strongest["partial_ci95_high"]:.3f},
BH q = {strongest["partial_q_value"]:.3g}; n = {int(strongest["n"])}), a {direction}
association. Immune associations passing 5% BH FDR: {significant_text}.

TACSTD2 expression itself correlated with purity (Spearman rho =
{purity_result.statistic:.3f}, p = {purity_result.pvalue:.3g}), supporting the requested
purity adjustment.

## Method

- Source: cBioPortal Datahub study `{STUDY_ID}` (OncoSG, Nat Genet 2020).
- Exposure: TACSTD2 (Entrez 4070) z-score from log RNA-seq V2 RSEM, standardized
  relative to all samples.
- Outcomes: the eight sample-level IMSIG immune scores supplied with the study.
  Proliferation and translation were analyzed as controls but excluded from immune-test
  multiple-testing correction.
- Primary statistic: partial Spearman correlation. TACSTD2, each score, and purity were
  rank transformed; each ranked variable was residualized on ranked purity; the
  residual correlation was tested with df = n - 3.
- Multiplicity: Benjamini-Hochberg correction across the eight immune outcomes.
- Confidence intervals: Fisher-z approximation for the partial correlation.
- Missing values: per-outcome complete cases; no imputation.

This is an association analysis of bulk tumors. Purity adjustment reduces one source of
mixture confounding but does not establish cell-intrinsic regulation or causality.

## Files

- `associations.tsv`: raw and purity-adjusted results.
- `analysis_samples.tsv`: joined public values used in the tests.
- `partial_spearman_forest.png`: primary estimates and 95% intervals.
- `purity_adjusted_residuals.png`: rank-residual diagnostic panels.
- `provenance.json`: source URLs, checksums, and counts.
- `analyze.py`: complete reproduction script.

Run `python3 analyze.py` from this directory. Downloads are cached outside the
repository by default; use `--cache-dir PATH` to select another cache.
"""
    (output_dir / "README.md").write_text(report, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=Path(tempfile.gettempdir()) / "oncosg_a1_cache",
    )
    args = parser.parse_args()
    output_dir = Path(__file__).resolve().parent
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    sources = {
        "expression": (EXPRESSION_URL, args.cache_dir / "expression.tsv"),
        "clinical": (CLINICAL_URL, args.cache_dir / "clinical_sample.tsv"),
    }
    for _, (url, path) in sources.items():
        if not path.exists():
            download(url, path)

    tacstd2 = load_tacstd2(sources["expression"][1])
    clinical = load_clinical(sources["clinical"][1])
    merged = clinical.join(tacstd2, how="inner")
    if merged.empty:
        raise ValueError("No sample IDs overlap between expression and clinical matrices")
    results = calculate_associations(merged)
    purity_complete = merged[["TACSTD2_Z", "PURITY"]].dropna()
    purity_result = stats.spearmanr(
        purity_complete["TACSTD2_Z"], purity_complete["PURITY"]
    )

    results.to_csv(output_dir / "associations.tsv", sep="\t", index=False)
    merged.reset_index().to_csv(
        output_dir / "analysis_samples.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
        float_format="%.6g",
    )
    plot_forest(results, output_dir / "partial_spearman_forest.png")
    plot_residuals(merged, results, output_dir / "purity_adjusted_residuals.png")
    provenance = {
        "study_id": STUDY_ID,
        "study_url": f"https://www.cbioportal.org/study/summary?id={STUDY_ID}",
        "datahub_revision": DATAHUB_REVISION,
        "sources": {
            name: {"url": url, "sha256": sha256(path), "bytes": path.stat().st_size}
            for name, (url, path) in sources.items()
        },
        "expression_samples": int(len(tacstd2)),
        "clinical_samples": int(len(clinical)),
        "overlap_samples": int(len(merged)),
        "tacstd2_entrez_id": 4070,
        "expression_profile": "rna_seq_v2_mrna_median_all_sample_Zscores",
    }
    (output_dir / "provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )
    write_report(output_dir, merged, results, purity_result, provenance)


if __name__ == "__main__":
    main()
