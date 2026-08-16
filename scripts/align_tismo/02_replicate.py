#!/usr/bin/env python3
"""Replicate the claim: Tacstd2 is higher in ICB-treated than baseline TISMO cohorts.

Claim under test (user-supplied):
    64 tumor models, baseline vs ICB, paired; 49/64 up; mean 1.05 -> 1.3;
    paired p = 5.8e-5.

The script reproduces that headline analysis exactly, then runs the sensitivity
analyses the headline design does not support (clustering of cohorts within cell
lines and studies, lung-only subset, robustness to influential cohorts).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tismo_client as tc

GENE = "Tacstd2"
BOOT_SEED = 20240816
N_BOOT = 20000


# --------------------------------------------------------------------------- #
# cohort table construction
# --------------------------------------------------------------------------- #
def load_cell_line_cancer_types() -> dict[str, str]:
    rows = json.loads((tc.DATA_DIR / "cellLineMeta.json").read_text())
    return {r["cellLine"]: (r.get("cancerType") or "").strip() for r in rows}


def build_cohort_table(sample_df: pd.DataFrame, cancer_types: dict[str, str]) -> pd.DataFrame:
    """Collapse per-sample values to one baseline/ICB pair per TISMO cohort.

    TISMO labels each cohort `<cellline>_<study>_<condition>_<regimen>(n=N)` and
    flags the study's own control arm with Baseline=1. Every treated (Baseline=0)
    arm in the in vivo Gene module carries an ICB agent, so the contrast is
    within-study ICB vs that study's matched control.
    """
    df = sample_df.copy()
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    df = df.dropna(subset=["value"])
    if "cohort" not in df.columns:
        df["cohort"] = df["cell_line"].map(tc.cohort_key)

    records = []
    for cohort, sub in df.groupby("cohort", sort=True):
        base = sub[sub["Baseline"] == 1]["value"]
        icb = sub[sub["Baseline"] == 0]["value"]
        if base.empty or icb.empty:
            continue
        cell_line = cohort.split("_")[0]
        records.append(
            {
                "cohort": cohort,
                "cell_line": cell_line,
                "cancer_type": cancer_types.get(cell_line, "unknown"),
                "study": sub["GSE_ID"].iloc[0],
                "regimen": ";".join(sorted(sub[sub["Baseline"] == 0]["Mouse_treatment"].unique())),
                "n_baseline": int(base.size),
                "n_icb": int(icb.size),
                "mean_baseline": float(base.mean()),
                "mean_icb": float(icb.mean()),
                "delta": float(icb.mean() - base.mean()),
                # TISMO's own precomputed within-cohort DESeq2 statistic.
                "tismo_within_cohort_p": _single_p(sub[sub["Baseline"] == 0]["pvalue"]),
            }
        )
    out = pd.DataFrame.from_records(records)
    return out.sort_values("cohort").reset_index(drop=True)


def _single_p(series: pd.Series) -> float:
    vals = pd.to_numeric(series, errors="coerce").dropna().unique()
    return float(np.min(vals)) if len(vals) else float("nan")


# --------------------------------------------------------------------------- #
# statistics
# --------------------------------------------------------------------------- #
def paired_summary(tab: pd.DataFrame, label: str) -> dict:
    """Every paired statistic a reader might mean by 'paired p', reported together."""
    base = tab["mean_baseline"].to_numpy(float)
    icb = tab["mean_icb"].to_numpy(float)
    d = icb - base
    n = len(d)
    up = int((d > 0).sum())
    down = int((d < 0).sum())
    tie = int((d == 0).sum())

    res: dict = {
        "subset": label,
        "n_cohorts": n,
        "n_cell_lines": int(tab["cell_line"].nunique()),
        "n_studies": int(tab["study"].nunique()),
        "n_cancer_types": int(tab["cancer_type"].nunique()),
        "n_up": up,
        "n_down": down,
        "n_tie": tie,
        "frac_up": up / n if n else float("nan"),
        "mean_baseline": float(base.mean()) if n else float("nan"),
        "mean_icb": float(icb.mean()) if n else float("nan"),
        "mean_delta": float(d.mean()) if n else float("nan"),
        "median_baseline": float(np.median(base)) if n else float("nan"),
        "median_icb": float(np.median(icb)) if n else float("nan"),
        "median_delta": float(np.median(d)) if n else float("nan"),
    }

    if n >= 2 and np.std(d) > 0:
        t_stat, t_p = stats.ttest_rel(icb, base)
        res["paired_t_stat"] = float(t_stat)
        res["paired_t_p"] = float(t_p)
        res["cohens_dz"] = float(d.mean() / d.std(ddof=1))
    else:
        res["paired_t_stat"] = res["paired_t_p"] = res["cohens_dz"] = float("nan")

    nz = d[d != 0]
    if len(nz) >= 1:
        w_stat, w_p = stats.wilcoxon(icb, base, zero_method="wilcox", alternative="two-sided")
        res["wilcoxon_stat"] = float(w_stat)
        res["wilcoxon_p"] = float(w_p)
        m = len(nz)
        res["rank_biserial"] = float(2 * w_stat / (m * (m + 1)) - 1)
        res["sign_test_p"] = float(stats.binomtest(up, up + down, 0.5).pvalue)
    else:
        res["wilcoxon_stat"] = res["wilcoxon_p"] = float("nan")
        res["rank_biserial"] = res["sign_test_p"] = float("nan")
    return res


def cluster_bootstrap(tab: pd.DataFrame, cluster_col: str, n_boot: int = N_BOOT) -> dict:
    """Resample whole cell lines/studies, because cohorts within one are not independent."""
    rng = np.random.default_rng(BOOT_SEED)
    groups = [g["delta"].to_numpy(float) for _, g in tab.groupby(cluster_col, sort=True)]
    k = len(groups)
    if k < 2:
        return {"cluster": cluster_col, "n_clusters": k, "ci_low": float("nan"),
                "ci_high": float("nan"), "p_two_sided": float("nan")}

    observed = float(np.concatenate(groups).mean())
    means = np.empty(n_boot)
    for b in range(n_boot):
        pick = rng.integers(0, k, size=k)
        means[b] = np.concatenate([groups[i] for i in pick]).mean()

    lo, hi = np.percentile(means, [2.5, 97.5])
    # Percentile-bootstrap two-sided p: how often the resampled mean crosses zero.
    centred = means - observed
    p = float((np.abs(centred) >= abs(observed)).mean())
    return {
        "cluster": cluster_col,
        "n_clusters": k,
        "observed_mean_delta": observed,
        "ci_low": float(lo),
        "ci_high": float(hi),
        "p_two_sided": max(p, 1.0 / n_boot),
    }


def collapse(tab: pd.DataFrame, key: str) -> pd.DataFrame:
    """One pair per cell line (or study): average its cohorts first."""
    def one_or_mixed(s: pd.Series) -> str:
        return s.iloc[0] if s.nunique() == 1 else "mixed"

    aggs = {
        "mean_baseline": ("mean_baseline", "mean"),
        "mean_icb": ("mean_icb", "mean"),
        "cell_line": ("cell_line", one_or_mixed),
        "study": ("study", one_or_mixed),
        "cancer_type": ("cancer_type", one_or_mixed),
        "n_cohorts": ("cohort", "size"),
    }
    aggs.pop(key, None)
    g = tab.groupby(key, sort=True).agg(**aggs).reset_index()
    g["cohort"] = g[key]
    g["delta"] = g["mean_icb"] - g["mean_baseline"]
    return g


def leave_one_cell_line_out(tab: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for line in sorted(tab["cell_line"].unique()):
        kept = tab[tab["cell_line"] != line]
        s = paired_summary(kept, f"drop_{line}")
        rows.append({
            "dropped_cell_line": line,
            "n_cohorts_removed": int((tab["cell_line"] == line).sum()),
            "n_cohorts_left": s["n_cohorts"],
            "n_up": s["n_up"],
            "mean_delta": s["mean_delta"],
            "paired_t_p": s["paired_t_p"],
            "wilcoxon_p": s["wilcoxon_p"],
        })
    return pd.DataFrame(rows).sort_values("wilcoxon_p", ascending=False)


# --------------------------------------------------------------------------- #
def main() -> int:
    tc.TABLES_DIR.mkdir(parents=True, exist_ok=True)
    cancer_types = load_cell_line_cancer_types()

    raw = tc.load_expression_csv(tc.DATA_DIR / f"vivo_expression_{GENE}.csv", gene=GENE)
    cohorts = build_cohort_table(raw, cancer_types)
    cohorts.to_csv(tc.TABLES_DIR / f"cohort_level_{GENE}.csv", index=False)

    # How many ICB cohorts does TISMO expose in total, across all genes fetched?
    all_cohorts: set[str] = set()
    for path in sorted(tc.DATA_DIR.glob("reference_genes/*.csv")) + \
            sorted(tc.ensure_null_genes().glob("*.csv")):
        try:
            other = tc.load_expression_csv(path)
        except Exception:
            continue
        all_cohorts |= set(other["cohort"].unique())
    missing_for_gene = sorted(all_cohorts - set(cohorts["cohort"])) if all_cohorts else []

    # ---- provenance / honest denominators -------------------------------- #
    raw_num = raw.copy()
    raw_num["value"] = pd.to_numeric(raw_num["value"], errors="coerce")
    provenance = {
        "gene": GENE,
        "source": "TISMO Gene module, in vivo (tismo.cistrome.org -> tismo.pku-genomics.org)",
        "value_units": "quantile-normalised, ComBat-corrected log-scale TPM as served by TISMO",
        "n_rows_returned": int(len(raw)),
        "n_rows_dropped_malformed": int(raw.attrs.get("dropped_malformed_rows", 0)),
        "n_rows_dropped_wrong_gene": int(raw.attrs.get("dropped_wrong_gene_rows", 0)),
        "n_unique_samples": int(raw["Samples"].nunique()),
        "n_baseline_samples": int(raw.loc[raw.Baseline == 1, "Samples"].nunique()),
        "n_icb_samples": int(raw.loc[raw.Baseline == 0, "Samples"].nunique()),
        "n_samples_in_multiple_cohorts": int(
            (raw.groupby("Samples")["cohort"].nunique() > 1).sum()
        ),
        "n_cohorts": int(cohorts.shape[0]),
        "n_icb_cohorts_tismo_exposes_overall": len(all_cohorts) or None,
        "cohorts_absent_for_this_gene": missing_for_gene,
        "n_cell_lines": int(cohorts["cell_line"].nunique()),
        "n_studies": int(cohorts["study"].nunique()),
        "n_cancer_types": int(cohorts["cancer_type"].nunique()),
        "cell_lines": sorted(cohorts["cell_line"].unique().tolist()),
        "studies": sorted(cohorts["study"].unique().tolist()),
        "frac_samples_exactly_zero": float((raw_num["value"] == 0).mean()),
        "frac_cohorts_with_baseline_below_0.5": float((cohorts["mean_baseline"] < 0.5).mean()),
    }
    (tc.TABLES_DIR / "provenance.json").write_text(json.dumps(provenance, indent=2))

    # ---- primary replication + subsets ----------------------------------- #
    summaries = [paired_summary(cohorts, "ALL cohorts (headline design)")]

    lung = cohorts[cohorts["cancer_type"].str.contains("Lung", case=False, na=False)]
    if len(lung):
        summaries.append(paired_summary(lung, "Lung cohorts only"))
    non_lung = cohorts[~cohorts["cancer_type"].str.contains("Lung", case=False, na=False)]
    summaries.append(paired_summary(non_lung, "Non-lung cohorts"))

    for ct, sub in cohorts.groupby("cancer_type", sort=True):
        summaries.append(paired_summary(sub, f"cancer type: {ct}"))

    by_line = collapse(cohorts, "cell_line")
    summaries.append(paired_summary(by_line, "collapsed to one pair per cell line"))
    by_study = collapse(cohorts, "study")
    summaries.append(paired_summary(by_study, "collapsed to one pair per study"))

    summary_df = pd.DataFrame(summaries)
    summary_df.to_csv(tc.TABLES_DIR / f"paired_summaries_{GENE}.csv", index=False)

    by_line.to_csv(tc.TABLES_DIR / f"cell_line_level_{GENE}.csv", index=False)
    by_study.to_csv(tc.TABLES_DIR / f"study_level_{GENE}.csv", index=False)

    # ---- clustering-aware inference -------------------------------------- #
    boot = pd.DataFrame([
        cluster_bootstrap(cohorts, "cell_line"),
        cluster_bootstrap(cohorts, "study"),
    ])
    boot.to_csv(tc.TABLES_DIR / f"cluster_bootstrap_{GENE}.csv", index=False)

    loo = leave_one_cell_line_out(cohorts)
    loo.to_csv(tc.TABLES_DIR / f"leave_one_cell_line_out_{GENE}.csv", index=False)

    # ---- TISMO's own within-cohort statistic ----------------------------- #
    has_p = cohorts.dropna(subset=["tismo_within_cohort_p"])
    tismo_view = {
        "n_cohorts_with_tismo_p": int(len(has_p)),
        "n_cohorts_without_tismo_p": int(len(cohorts) - len(has_p)),
        "n_tismo_p_below_0.05": int((has_p["tismo_within_cohort_p"] < 0.05).sum()),
        "n_tismo_p_below_0.05_and_up": int(
            ((has_p["tismo_within_cohort_p"] < 0.05) & (has_p["delta"] > 0)).sum()
        ),
        "n_tismo_p_below_0.05_and_down": int(
            ((has_p["tismo_within_cohort_p"] < 0.05) & (has_p["delta"] < 0)).sum()
        ),
    }
    (tc.TABLES_DIR / "tismo_within_cohort_stats.json").write_text(json.dumps(tismo_view, indent=2))

    # ---- console report --------------------------------------------------- #
    pd.set_option("display.width", 200)
    print("=" * 78)
    print(f"TISMO {GENE}: baseline vs ICB, paired at the cohort level")
    print("=" * 78)
    print(json.dumps(provenance, indent=2))
    print()
    cols = ["subset", "n_cohorts", "n_cell_lines", "n_studies", "n_up", "n_down",
            "mean_baseline", "mean_icb", "mean_delta", "paired_t_p", "wilcoxon_p", "sign_test_p"]
    print(summary_df[cols].to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print()
    print("Cluster bootstrap on the mean delta:")
    print(boot.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print()
    print("Leave-one-cell-line-out (sorted by weakest resulting p):")
    print(loo.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print()
    print("TISMO's own within-cohort DESeq2 statistic:")
    print(json.dumps(tismo_view, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
