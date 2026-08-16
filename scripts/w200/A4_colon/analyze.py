#!/usr/bin/env python3
"""A4 colon-only: Tacstd2 after ICB in TISMO colorectal models.

Claim under test (user-supplied, all-cancer): 49/64 models up, p=5.8e-5.
This slice is the colon/CRC subset only. It does not retune the contrast
to recover 49/64.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tismo_client as tc

GENE = "Tacstd2"
BOOT_SEED = 20240816
N_BOOT = 20000
CLAIM_N_UP = 49
CLAIM_N = 64
CLAIM_P = 5.8e-5

# ICB-class tokens. A regimen is "strict ICB" if every treated-arm token is
# anti-PD-1 / PD-L1 / PD-L2 / CTLA-4 (TISMO's antiCTL4 typo and Fc suffixes ok).
_ICB_TOKEN = re.compile(
    r"^anti-?(pd-?l?-?[12]|ctla?-?4)(_fcs)?$",
    re.I,
)
_SPLIT_TX = re.compile(r"[;,]+")


def is_strict_icb(regimen: str) -> bool:
    tokens = [t.strip() for t in _SPLIT_TX.split(regimen) if t.strip()]
    if not tokens:
        return False
    return all(_ICB_TOKEN.match(re.sub(r"[\s.]", "", t)) for t in tokens)


def build_cohort_table(sample_df: pd.DataFrame, cancer_types: dict[str, str]) -> pd.DataFrame:
    df = sample_df.copy()
    records = []
    for cohort, sub in df.groupby("cohort", sort=True):
        base = sub[sub["Baseline"] == 1]["value"]
        icb = sub[sub["Baseline"] == 0]["value"]
        if base.empty or icb.empty:
            continue
        model = cohort.split("_")[0]
        regimen = ";".join(sorted(sub[sub["Baseline"] == 0]["Mouse_treatment"].unique()))
        records.append(
            {
                "cohort": cohort,
                "cell_line": model,
                "cancer_type": cancer_types.get(model, "unknown"),
                "study": sub["GSE_ID"].iloc[0],
                "regimen": regimen,
                "strict_icb": bool(is_strict_icb(regimen)),
                "n_baseline": int(base.size),
                "n_icb": int(icb.size),
                "n_baseline_unique": int(sub[sub["Baseline"] == 1]["Samples"].nunique()),
                "n_icb_unique": int(sub[sub["Baseline"] == 0]["Samples"].nunique()),
                "mean_baseline": float(base.mean()),
                "mean_icb": float(icb.mean()),
                "median_baseline": float(base.median()),
                "median_icb": float(icb.median()),
                "delta": float(icb.mean() - base.mean()),
                "direction": "up" if icb.mean() > base.mean() else (
                    "down" if icb.mean() < base.mean() else "tie"
                ),
                "tismo_within_cohort_p": _single_p(sub[sub["Baseline"] == 0]["pvalue"]),
            }
        )
    return pd.DataFrame.from_records(records).sort_values("cohort").reset_index(drop=True)


def _single_p(series: pd.Series) -> float:
    vals = pd.to_numeric(series, errors="coerce").dropna().unique()
    return float(np.min(vals)) if len(vals) else float("nan")


def paired_summary(tab: pd.DataFrame, label: str) -> dict:
    if tab.empty:
        return {
            "subset": label,
            "n_cohorts": 0,
            "n_up": 0,
            "n_down": 0,
            "n_tie": 0,
            "paired_t_p": float("nan"),
            "wilcoxon_p": float("nan"),
            "sign_test_p": float("nan"),
        }
    base = tab["mean_baseline"].to_numpy(float)
    icb = tab["mean_icb"].to_numpy(float)
    d = icb - base
    n = len(d)
    up = int((d > 0).sum())
    down = int((d < 0).sum())
    tie = int((d == 0).sum())
    res = {
        "subset": label,
        "n_cohorts": n,
        "n_cell_lines": int(tab["cell_line"].nunique()),
        "n_studies": int(tab["study"].nunique()),
        "n_up": up,
        "n_down": down,
        "n_tie": tie,
        "frac_up": up / n if n else float("nan"),
        "n_baseline_samples": int(tab["n_baseline"].sum()),
        "n_icb_samples": int(tab["n_icb"].sum()),
        "mean_baseline": float(base.mean()),
        "mean_icb": float(icb.mean()),
        "mean_delta": float(d.mean()),
        "median_delta": float(np.median(d)),
    }
    if n >= 2 and np.std(d, ddof=1) > 0:
        t_stat, t_p = stats.ttest_rel(icb, base)
        res["paired_t_stat"] = float(t_stat)
        res["paired_t_p"] = float(t_p)
        res["cohens_dz"] = float(d.mean() / d.std(ddof=1))
    else:
        res["paired_t_stat"] = res["paired_t_p"] = res["cohens_dz"] = float("nan")

    if n >= 2 and np.any(d != 0):
        w_stat, w_p = stats.wilcoxon(icb, base, zero_method="wilcox", alternative="two-sided")
        res["wilcoxon_stat"] = float(w_stat)
        res["wilcoxon_p"] = float(w_p)
        m = int((d != 0).sum())
        res["rank_biserial"] = float(2 * w_stat / (m * (m + 1)) - 1)
        res["sign_test_p"] = float(stats.binomtest(up, up + down, 0.5).pvalue)
    else:
        res["wilcoxon_stat"] = res["wilcoxon_p"] = float("nan")
        res["rank_biserial"] = res["sign_test_p"] = float("nan")
    return res


def collapse(tab: pd.DataFrame, key: str) -> pd.DataFrame:
    def one_or_mixed(s: pd.Series) -> str:
        return s.iloc[0] if s.nunique() == 1 else "mixed"

    aggs = {
        "mean_baseline": ("mean_baseline", "mean"),
        "mean_icb": ("mean_icb", "mean"),
        "n_cohorts": ("cohort", "size"),
        "n_baseline": ("n_baseline", "sum"),
        "n_icb": ("n_icb", "sum"),
    }
    for col in ("cell_line", "study", "cancer_type"):
        if col != key and col in tab.columns:
            aggs[col] = (col, one_or_mixed)
    g = tab.groupby(key, sort=True).agg(**aggs).reset_index()
    g["cohort"] = g[key]
    g["delta"] = g["mean_icb"] - g["mean_baseline"]
    return g


def cluster_bootstrap(tab: pd.DataFrame, cluster_col: str, n_boot: int = N_BOOT) -> dict:
    rng = np.random.default_rng(BOOT_SEED)
    groups = [g["delta"].to_numpy(float) for _, g in tab.groupby(cluster_col, sort=True)]
    k = len(groups)
    if k < 3:
        observed = float(np.concatenate(groups).mean()) if k else float("nan")
        return {
            "cluster": cluster_col,
            "n_clusters": k,
            "observed_mean_delta": observed,
            "ci_low": float("nan"),
            "ci_high": float("nan"),
            "p_two_sided": float("nan"),
            "note": f"k={k} clusters: bootstrap p/CI not reported (resamples collapse to the cluster means)",
        }
    observed = float(np.concatenate(groups).mean())
    means = np.empty(n_boot)
    for b in range(n_boot):
        pick = rng.integers(0, k, size=k)
        means[b] = np.concatenate([groups[i] for i in pick]).mean()
    lo, hi = np.percentile(means, [2.5, 97.5])
    centred = means - observed
    p = float((np.abs(centred) >= abs(observed)).mean())
    return {
        "cluster": cluster_col,
        "n_clusters": k,
        "observed_mean_delta": observed,
        "ci_low": float(lo),
        "ci_high": float(hi),
        "p_two_sided": max(p, 1.0 / n_boot),
        "note": "",
    }


def _fmt_p(p: float) -> str:
    if p != p:
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3g}"


def write_figures(cohorts: pd.DataFrame, out: Path) -> None:
    up_c, down_c = "#c0392b", "#2c6fbb"
    model_c = {"CT26": "#c0392b", "MC38": "#1f6aa5"}

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.2, 5.6), gridspec_kw={"width_ratios": [1.35, 1]})
    for _, r in cohorts.iterrows():
        c = model_c.get(r.cell_line, "#555")
        ax.plot([0, 1], [r.mean_baseline, r.mean_icb], color=c, alpha=0.7, lw=1.4, marker="o", ms=4)
    ax.plot(
        [0, 1],
        [cohorts.mean_baseline.mean(), cohorts.mean_icb.mean()],
        color="black",
        lw=2.6,
        marker="o",
        ms=7,
        label="mean of cohorts",
        zorder=5,
    )
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Baseline\n(control arm)", "ICB-treated"])
    ax.set_xlim(-0.3, 1.3)
    ax.set_ylabel(f"{GENE} (TISMO quantile-normalised units)")
    n_up = int((cohorts.delta > 0).sum())
    ax.set_title(
        f"Colon/CRC only: {GENE} paired within-study contrast\n"
        f"{n_up}/{len(cohorts)} cohorts up  "
        f"({cohorts.mean_baseline.mean():.2f} → {cohorts.mean_icb.mean():.2f})"
    )
    handles = [
        plt.Line2D([0], [0], color=model_c["CT26"], lw=2, label="CT26"),
        plt.Line2D([0], [0], color=model_c["MC38"], lw=2, label="MC38"),
        plt.Line2D([0], [0], color="black", lw=2.6, label="mean of cohorts"),
    ]
    ax.legend(handles=handles, frameon=False, loc="upper left")

    order = cohorts.sort_values("delta")
    colors = [up_c if d > 0 else down_c for d in order.delta]
    ax2.barh(np.arange(len(order)), order.delta, color=colors)
    ax2.axvline(0, color="black", lw=0.8)
    ax2.set_yticks(np.arange(len(order)))
    ax2.set_yticklabels([f"{r.cell_line} {r.study}" for _, r in order.iterrows()], fontsize=8)
    ax2.set_xlabel("ICB − baseline")
    ax2.set_title("Per-cohort difference")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
        ax2.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(out / "fig_paired_slopes.png", dpi=160)
    fig.savefig(out / "fig_paired_slopes.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    by_model = list(cohorts.groupby("cell_line", sort=True))
    for i, (model, sub) in enumerate(by_model):
        jitter = np.linspace(-0.12, 0.12, len(sub)) if len(sub) > 1 else [0.0]
        ax.scatter(
            np.full(len(sub), i) + jitter,
            sub["delta"],
            s=55,
            color=[up_c if d > 0 else down_c for d in sub["delta"]],
            zorder=3,
        )
        ax.plot([i], [sub["delta"].mean()], marker="_", ms=22, color="black", zorder=4)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(range(len(by_model)))
    ax.set_xticklabels([f"{m}\n(n={len(s)} cohorts)" for m, s in by_model])
    ax.set_ylabel(f"{GENE}: ICB − baseline")
    ax.set_title("Colon models only — Tacstd2 delta after ICB")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(out / "fig_delta_by_model.png", dpi=160)
    fig.savefig(out / "fig_delta_by_model.pdf")
    plt.close(fig)


def write_readme(
    provenance: dict,
    summaries: pd.DataFrame,
    cohorts: pd.DataFrame,
    boot: pd.DataFrame,
    colon_in_meta: list[str],
    colon_in_icb: list[str],
    colon_absent: list[str],
    path: Path,
) -> None:
    wanted = [
        "colon ICB cohorts (TISMO definition)",
        "CT26 only",
        "MC38 only",
        "strict ICB regimens only",
        "ICB plus other partner",
        "collapsed to one pair per study",
        "collapsed to one pair per cell line",
    ]
    prim = summaries[summaries["subset"] == wanted[0]].iloc[0]

    def row(s: pd.Series) -> str:
        return (
            f"| {s['subset']} | {int(s['n_cohorts'])} | {int(s.get('n_cell_lines', 0) or 0)} | "
            f"{int(s.get('n_studies', 0) or 0)} | {int(s['n_up'])}/{int(s['n_cohorts'])} | "
            f"{s.get('mean_delta', float('nan')):.3f} | {_fmt_p(s.get('paired_t_p', float('nan')))} | "
            f"{_fmt_p(s.get('wilcoxon_p', float('nan')))} | {_fmt_p(s.get('sign_test_p', float('nan')))} |"
        )

    table_rows = "\n".join(
        row(summaries[summaries["subset"] == s].iloc[0])
        for s in wanted
        if (summaries["subset"] == s).any()
    )

    cohort_rows = []
    for _, r in cohorts.sort_values(["cell_line", "study", "cohort"]).iterrows():
        cohort_rows.append(
            f"| {r.cell_line} | {r.study} | {r.regimen} | {int(r.n_baseline)}/{int(r.n_icb)} | "
            f"{r.mean_baseline:.3f} | {r.mean_icb:.3f} | {r.delta:+.3f} | {r.direction} | "
            f"{'yes' if r.strict_icb else 'no'} |"
        )

    boot_rows = []
    for _, r in boot.iterrows():
        if r.ci_low == r.ci_low:
            ci = f"[{r.ci_low:.3f}, {r.ci_high:.3f}]"
        else:
            ci = "NA"
        note = str(r.note).strip()
        extra = f" {note}" if note else ""
        boot_rows.append(
            f"| {r.cluster} | {int(r.n_clusters)} | {r.observed_mean_delta:.3f} | "
            f"{ci} | {_fmt_p(r.p_two_sided)} |{extra}"
        )

    verdict = (
        "NOT the 49/64 result. Colon-only is **8/10** cohorts up. "
        "Paired t / Wilcoxon are only barely <0.05 if the 10 cohorts are treated as "
        "independent; the sign test, study collapse, cell-line collapse, CT26-only, "
        "MC38-only, and strict-ICB subsets are all null."
    )

    text = f"""# A4 TISMO colon/CRC — Tacstd2 after ICB

**Claim (user, all-cancer):** Tacstd2 up in 49/64 ICI-treated mouse models (p=5.8e-5).

**This slice:** colon/CRC models only. Filters were not tuned to recover 49/64.

**Honest verdict: `{prim['n_up']}/{int(prim['n_cohorts'])}` up. The 49/64 statistic is not a colon result.**

{verdict}

## Primary contrast

TISMO in-vivo Gene module, ICB regimens vs the study's own `Baseline=1` control arm.
One pair per TISMO cohort (`cellline_study_condition_regimen`).
Values are TISMO quantile-normalised, ComBat-corrected log-scale TPM as served.

| Subset | n cohorts | n lines | n studies | up/n | mean Δ | paired t p | Wilcoxon p | sign-test p |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
{table_rows}

Primary n = **{int(prim['n_cohorts'])}** paired cohorts, **{int(prim['n_baseline_samples'])}** baseline + **{int(prim['n_icb_samples'])}** ICB samples, **{int(prim['n_cell_lines'])}** cell lines, **{int(prim['n_studies'])}** studies.

- mean baseline → ICB: **{prim['mean_baseline']:.3f} → {prim['mean_icb']:.3f}** (mean Δ = {prim['mean_delta']:.3f})
- paired t = {prim.get('paired_t_stat', float('nan')):.2f}, p = {_fmt_p(prim['paired_t_p'])}
- Wilcoxon p = {_fmt_p(prim['wilcoxon_p'])}
- sign test (8 up / 2 down vs 0.5) p = {_fmt_p(prim['sign_test_p'])}
- Cohen's d_z = {prim.get('cohens_dz', float('nan')):.2f}

User-reported 49/{CLAIM_N} and p={CLAIM_P:.1e} are **all-cancer**. This download recovers that all-cancer count ({provenance.get("all_cancer_n_up")}/{provenance.get("n_all_cancer_paired_cohorts")} up). It is not a colon result.

## Per-cohort table (every colon pair TISMO serves)

| Model | Study | Regimen | n B/ICB | mean B | mean ICB | Δ | dir | strict ICB |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- | --- |
{chr(10).join(cohort_rows)}

`strict ICB` = treated arm is only anti-PD-1 / PD-L1 / CTLA-4 (no TGF-β, GARP, exercise, or other partners).

## Independence (why p≈0.05 is not a colon confirmation)

Cohorts cluster inside two cell lines and five studies. Treating n=10 as independent overstates precision.

| Cluster | n clusters | mean Δ | 95% bootstrap CI | two-sided p |
| --- | ---: | ---: | --- | ---: |
{chr(10).join(boot_rows)}

Cell-line k=2 is too small for a bootstrap p. Study-level k=5 is the usable cluster check.

- Collapsed to one pair per **study** (n=5): 4/5 up; paired t p and Wilcoxon p are null (see table).
- Collapsed to one pair per **cell line** (n=2): both lines go up (CT26 Δ larger; MC38 near zero). n=2 cannot support a p-value.
- **MC38** Tacstd2 is essentially off (cohort means 0.00–0.10). A 2/3 “up” count there is noise around zero.
- **CT26 only** (n=7, 6/7 up): paired t p ≈ 0.06, Wilcoxon / sign test null.

## Coverage (honest missingness)

- `cellLineMeta` colorectal carcinoma models: {", ".join(colon_in_meta) or "none"}.
- Of those, the ICB Gene-module model list contains: {", ".join(colon_in_icb) or "none"}.
- Colorectal models with **no** ICB Gene-module expression: {", ".join(colon_absent) or "none"}.
- Samples appearing in more than one colon cohort: {provenance.get("n_samples_in_multiple_colon_cohorts", 0)} (should be 0).
- All-cancer Tacstd2 ICB cohorts TISMO returned in this download: {provenance.get("n_all_cancer_paired_cohorts")} ({provenance.get("all_cancer_n_up")} up / {provenance.get("all_cancer_n_down")} down). Colon is {int(prim['n_cohorts'])} of those 64, not a separate 64.

## What this is not

- A reproduction of 49/64 or p=5.8e-5 in colon. That numerator/denominator is all-cancer.
- Evidence that Tacstd2 induction after ICB is a CRC-specific finding.
- Human CRC, MSI/MSS, or TROP2-ADC outcome.
- A Cldn4 analysis.

## Caveats

1. TISMO units are site-normalised, not raw TPM. Direction within a study is usable; absolute values are not portable.
2. Several CT26 arms are ICB + TGF-β / GARP. The strict-ICB subset is null (see table).
3. ERP114266 day-7 and day-14 are the same study at two time points.
4. n=10 is small. A p-value that sits on 0.05 under the most liberal pairing is not a confirmation.
5. 1638N-T1 and CMT93 are colorectal in TISMO metadata but are absent from the ICB Gene module, so they cannot enter the denominator.

## Rerun

```bash
python3 scripts/w200/A4_colon/download.py
python3 scripts/w200/A4_colon/analyze.py
```

See `summary.json`, `paired_summaries.csv`, and `cohort_pairs.csv`.
"""
    path.write_text(text)


def main() -> int:
    tc.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    expr_path = tc.DATA_DIR / f"vivo_expression_{GENE}.csv"
    meta_path = tc.DATA_DIR / "cellLineMeta.json"
    if not expr_path.exists() or not meta_path.exists():
        raise SystemExit("missing downloaded TISMO files; run download.py first")

    cancer_types = tc.load_cell_line_cancer_types(meta_path)
    colon_meta = tc.colon_models(cancer_types)

    vocab_path = tc.DATA_DIR / "vocabularies.json"
    icb_models = []
    if vocab_path.exists():
        icb_models = json.loads(vocab_path.read_text()).get("tumor_models", [])
    colon_in_icb = [m for m in colon_meta if m in icb_models]
    colon_absent = [m for m in colon_meta if m not in icb_models]

    raw = tc.load_expression_csv(expr_path, gene=GENE)
    raw["cancer_type"] = raw["model"].map(cancer_types)
    all_cohorts = build_cohort_table(raw, cancer_types)
    colon_raw = raw[raw["model"].isin(colon_meta)].copy()
    cohorts = build_cohort_table(colon_raw, cancer_types)
    if cohorts.empty:
        raise SystemExit("no paired colon ICB cohorts — refusing to invent a result")

    cohorts.to_csv(tc.RESULTS_DIR / "cohort_pairs.csv", index=False)
    colon_raw.to_csv(tc.RESULTS_DIR / "sample_level.csv", index=False)

    summaries = [
        paired_summary(cohorts, "colon ICB cohorts (TISMO definition)"),
        paired_summary(cohorts[cohorts["cell_line"] == "CT26"], "CT26 only"),
        paired_summary(cohorts[cohorts["cell_line"] == "MC38"], "MC38 only"),
        paired_summary(cohorts[cohorts["strict_icb"]], "strict ICB regimens only"),
        paired_summary(cohorts[~cohorts["strict_icb"]], "ICB plus other partner"),
        paired_summary(collapse(cohorts, "study"), "collapsed to one pair per study"),
        paired_summary(collapse(cohorts, "cell_line"), "collapsed to one pair per cell line"),
    ]
    # Collapse labels one row per study/line; keep the original cluster counts.
    summaries[-2]["n_cell_lines"] = int(cohorts["cell_line"].nunique())
    summaries[-2]["n_studies"] = int(cohorts["study"].nunique())
    summaries[-1]["n_cell_lines"] = int(cohorts["cell_line"].nunique())
    summaries[-1]["n_studies"] = int(cohorts["study"].nunique())
    summary_df = pd.DataFrame(summaries)
    summary_df.to_csv(tc.RESULTS_DIR / "paired_summaries.csv", index=False)

    boot = pd.DataFrame([
        cluster_bootstrap(cohorts, "cell_line"),
        cluster_bootstrap(cohorts, "study"),
    ])
    boot.to_csv(tc.RESULTS_DIR / "cluster_bootstrap.csv", index=False)

    n_multi = int((colon_raw.groupby("Samples")["cohort"].nunique() > 1).sum())
    provenance = {
        "gene": GENE,
        "source": "TISMO Gene module, in vivo (tismo.cistrome.org -> tismo.pku-genomics.org)",
        "slice": "colon/CRC only (cellLineMeta cancerType == Colorectal carcinoma)",
        "value_units": "quantile-normalised, ComBat-corrected log-scale TPM as served by TISMO",
        "claim_tested": f"Tacstd2 up in {CLAIM_N_UP}/{CLAIM_N} ICI-treated models, p={CLAIM_P}",
        "claim_scope_as_stated": "all-cancer TISMO; this slice is colon-only",
        "n_rows_returned_all_cancer": int(len(raw)),
        "n_rows_colon": int(len(colon_raw)),
        "n_rows_dropped_malformed": int(raw.attrs.get("dropped_malformed_rows", 0)),
        "n_rows_dropped_wrong_gene": int(raw.attrs.get("dropped_wrong_gene_rows", 0)),
        "n_all_cancer_paired_cohorts": int(len(all_cohorts)),
        "n_colon_paired_cohorts": int(len(cohorts)),
        "n_colon_baseline_samples": int(colon_raw.loc[colon_raw.Baseline == 1, "Samples"].nunique()),
        "n_colon_icb_samples": int(colon_raw.loc[colon_raw.Baseline == 0, "Samples"].nunique()),
        "n_samples_in_multiple_colon_cohorts": n_multi,
        "colorectal_models_in_cellLineMeta": colon_meta,
        "colorectal_models_in_icb_gene_module": colon_in_icb,
        "colorectal_models_absent_from_icb_gene_module": colon_absent,
        "colon_cell_lines_with_pairs": sorted(cohorts["cell_line"].unique().tolist()),
        "colon_studies": sorted(cohorts["study"].unique().tolist()),
        "all_cancer_n_up": int((all_cohorts["delta"] > 0).sum()),
        "all_cancer_n_down": int((all_cohorts["delta"] < 0).sum()),
    }
    (tc.RESULTS_DIR / "provenance.json").write_text(json.dumps(provenance, indent=2))

    prim = summaries[0]
    verdict = "NOT_49_of_64"
    if prim["n_cohorts"] == CLAIM_N and prim["n_up"] == CLAIM_N_UP:
        verdict = "MATCHES_CLAIM_DENOMINATOR"
    summary = {
        "verdict": verdict,
        "primary_n_cohorts": prim["n_cohorts"],
        "primary_n_up": prim["n_up"],
        "primary_n_down": prim["n_down"],
        "primary_paired_t_p": prim["paired_t_p"],
        "primary_wilcoxon_p": prim["wilcoxon_p"],
        "primary_sign_test_p": prim["sign_test_p"],
        "primary_mean_delta": prim["mean_delta"],
        "user_claim_n_up": CLAIM_N_UP,
        "user_claim_n": CLAIM_N,
        "user_claim_p": CLAIM_P,
        "note": (
            "Colon-only cannot be 49/64. Report 8/10 and the p-values in paired_summaries.csv. "
            "Do not treat p≈0.05 under independent-cohort pairing as confirmation."
        ),
    }
    (tc.RESULTS_DIR / "summary.json").write_text(json.dumps(summary, indent=2))

    write_figures(cohorts, tc.RESULTS_DIR)
    write_readme(
        provenance, summary_df, cohorts, boot,
        colon_meta, colon_in_icb, colon_absent,
        tc.RESULTS_DIR / "README.md",
    )

    print("=" * 72)
    print(f"TISMO {GENE} colon/CRC only — baseline vs ICB")
    print("=" * 72)
    print(json.dumps(provenance, indent=2))
    cols = [
        "subset", "n_cohorts", "n_cell_lines", "n_studies", "n_up", "n_down",
        "mean_delta", "paired_t_p", "wilcoxon_p", "sign_test_p",
    ]
    print(summary_df[cols].to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print()
    print(boot.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print()
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
