#!/usr/bin/env python3
"""Per-cohort CLDN4-high vs ICI response odds ratios and random-effects meta.

Primary question (claim B5): is CLDN4-high associated with ICI non-response at
OR ≈ 0.42 across 11 named GBM/NSCLC/melanoma/RCC/urothelial cohorts?

This script does not invent missing cohorts. GBM is absent from the public
ORCESTRA/PredictIO processed objects and is not silently filled in.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm

ROOT = Path(__file__).resolve().parents[1]
DER = ROOT / "data" / "derived"
OUT = ROOT / "results" / "claim_B5"
OUT.mkdir(parents=True, exist_ok=True)

# Haldane–Anscombe continuity correction applied when any 2x2 cell is 0.
CORR = 0.5


def recist_or(recist: pd.Series) -> pd.Series:
    """CR/PR = R, PD = NR, everything else (SD, mixed, missing) = NA."""
    r = recist.astype(str).str.strip().str.upper()
    out = pd.Series(np.nan, index=recist.index, dtype=object)
    out[r.isin(["CR", "PR"])] = "R"
    out[r.isin(["PD", "PD*", "PD "])] = "NR"
    return out


def twobytwo(high: np.ndarray, resp: np.ndarray) -> dict:
    """high True/False, resp 'R'/'NR'. Returns a,b,c,d and OR on log scale."""
    a = int(np.sum(high & (resp == "R")))  # high responders
    b = int(np.sum(high & (resp == "NR")))  # high nonresponders
    c = int(np.sum(~high & (resp == "R")))  # low responders
    d = int(np.sum(~high & (resp == "NR")))  # low nonresponders
    cells = np.array([a, b, c, d], dtype=float)
    used_corr = bool((cells == 0).any())
    if used_corr:
        cells = cells + CORR
    aa, bb, cc, dd = cells
    lor = np.log((aa * dd) / (bb * cc))
    se = np.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    z = lor / se
    p = 2 * stats.norm.sf(abs(z))
    return {
        "n_high_R": a,
        "n_high_NR": b,
        "n_low_R": c,
        "n_low_NR": d,
        "n_high": a + b,
        "n_low": c + d,
        "n_R": a + c,
        "n_NR": b + d,
        "n": a + b + c + d,
        "continuity_correction": used_corr,
        "logOR": float(lor),
        "SE": float(se),
        "OR": float(np.exp(lor)),
        "OR_lo": float(np.exp(lor - 1.96 * se)),
        "OR_hi": float(np.exp(lor + 1.96 * se)),
        "p": float(p),
    }


def median_split(x: np.ndarray) -> np.ndarray:
    return x >= np.median(x)


def tertile_high_low(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    q1, q2 = np.quantile(x, [1 / 3, 2 / 3])
    high = x >= q2
    low = x <= q1
    keep = high | low
    return high[keep], keep


def quartile_high_low(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    q1, q3 = np.quantile(x, [0.25, 0.75])
    high = x >= q3
    low = x <= q1
    keep = high | low
    return high[keep], keep


def re_meta(rows: list[dict]) -> dict:
    """DerSimonian–Laird random-effects on logOR; also reports IV fixed effect."""
    if not rows:
        return {"k": 0}
    yi = np.array([r["logOR"] for r in rows], dtype=float)
    vi = np.array([r["SE"] ** 2 for r in rows], dtype=float)
    w = 1.0 / vi
    ybar = np.sum(w * yi) / np.sum(w)
    Q = float(np.sum(w * (yi - ybar) ** 2))
    k = len(rows)
    df = k - 1
    c = float(np.sum(w) - np.sum(w**2) / np.sum(w)) if k > 1 else np.nan
    tau2 = max(0.0, (Q - df) / c) if k > 1 and c > 0 else 0.0
    I2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 and k > 1 else 0.0
    p_het = float(stats.chi2.sf(Q, df)) if k > 1 else np.nan
    wstar = 1.0 / (vi + tau2)
    yre = float(np.sum(wstar * yi) / np.sum(wstar))
    se_re = float(np.sqrt(1.0 / np.sum(wstar)))
    z = yre / se_re
    p = 2 * stats.norm.sf(abs(z))
    yfe = float(ybar)
    se_fe = float(np.sqrt(1.0 / np.sum(w)))
    zfe = yfe / se_fe
    pfe = 2 * stats.norm.sf(abs(zfe))
    return {
        "k": k,
        "n_total": int(sum(r["n"] for r in rows)),
        "method": "DerSimonian-Laird random effects on logOR",
        "OR": float(np.exp(yre)),
        "OR_lo": float(np.exp(yre - 1.96 * se_re)),
        "OR_hi": float(np.exp(yre + 1.96 * se_re)),
        "logOR": yre,
        "SE": se_re,
        "p": float(p),
        "FE_OR": float(np.exp(yfe)),
        "FE_OR_lo": float(np.exp(yfe - 1.96 * se_fe)),
        "FE_OR_hi": float(np.exp(yfe + 1.96 * se_fe)),
        "FE_p": float(pfe),
        "tau2": tau2,
        "Q": Q,
        "I2": I2,
        "p_heterogeneity": p_het,
    }


def mannwhitney(x_r: np.ndarray, x_nr: np.ndarray) -> dict:
    if len(x_r) < 2 or len(x_nr) < 2:
        return {"U_p": np.nan, "cliffs_delta": np.nan, "median_R": np.nan, "median_NR": np.nan}
    u, p = stats.mannwhitneyu(x_r, x_nr, alternative="two-sided")
    # Cliff's delta: (n_R>n_NR - n_R<n_NR) / (nR*nNR)
    gt = np.sum(x_r[:, None] > x_nr[None, :])
    lt = np.sum(x_r[:, None] < x_nr[None, :])
    delta = (gt - lt) / (len(x_r) * len(x_nr))
    return {
        "U_p": float(p),
        "cliffs_delta": float(delta),
        "median_R": float(np.median(x_r)),
        "median_NR": float(np.median(x_nr)),
    }


def logistic_or_per_sd(x: np.ndarray, y: np.ndarray) -> dict:
    """y is 1=R, 0=NR. OR per +1 SD CLDN4."""
    if y.sum() < 2 or (1 - y).sum() < 2 or np.nanstd(x) == 0:
        return {"logit_OR_per_SD": np.nan, "logit_OR_lo": np.nan, "logit_OR_hi": np.nan, "logit_p": np.nan}
    z = (x - np.mean(x)) / np.std(x, ddof=1)
    X = sm.add_constant(z)
    try:
        fit = sm.Logit(y, X).fit(disp=False)
        b = float(fit.params[1])
        se = float(fit.bse[1])
        return {
            "logit_OR_per_SD": float(np.exp(b)),
            "logit_OR_lo": float(np.exp(b - 1.96 * se)),
            "logit_OR_hi": float(np.exp(b + 1.96 * se)),
            "logit_p": float(fit.pvalues[1]),
        }
    except Exception:
        return {"logit_OR_per_SD": np.nan, "logit_OR_lo": np.nan, "logit_OR_hi": np.nan, "logit_p": np.nan}


def analyze_cohort(df: pd.DataFrame, response_col: str, split: str) -> dict | None:
    sub = df[df[response_col].isin(["R", "NR"]) & np.isfinite(df["CLDN4"])].copy()
    if len(sub) < 8:
        return None
    x = sub["CLDN4"].to_numpy()
    y = sub[response_col].to_numpy()
    if len(np.unique(y)) < 2:
        return None
    iqr = float(np.percentile(x, 75) - np.percentile(x, 25))
    xmin = float(np.min(x))
    # floor: values within 0.05 of the study minimum
    frac_floor = float(np.mean(np.isclose(x, xmin, atol=0.05))) if np.isfinite(xmin) else np.nan
    if split == "median":
        high = median_split(x)
        keep = np.ones(len(x), dtype=bool)
    elif split == "tertile":
        high, keep = tertile_high_low(x)
        y = y[keep]
        x = x[keep]
    elif split == "quartile":
        high, keep = quartile_high_low(x)
        y = y[keep]
        x = x[keep]
    else:
        raise ValueError(split)
    if high.sum() < 2 or (~high).sum() < 2 or len(np.unique(y)) < 2:
        return None
    orstats = twobytwo(high, y)
    mw = mannwhitney(x[y == "R"], x[y == "NR"])
    logit = logistic_or_per_sd(sub["CLDN4"].to_numpy(), (sub[response_col] == "R").astype(int).to_numpy())
    out = {
        "study": df["study"].iloc[0],
        "indication": df["indication"].iloc[0],
        "io_class": df["io_class"].iloc[0],
        "in_scope": bool(df["in_scope"].iloc[0]),
        "split": split,
        "response_def": response_col,
        "cldn4_iqr": iqr,
        "cldn4_frac_floor": frac_floor,
        "low_dynamic_range": iqr < 0.5,
        **orstats,
        **mw,
        **logit,
    }
    return out


def forest(rows: list[dict], summary: dict, title: str, path: Path) -> None:
    rows = sorted(rows, key=lambda r: (r["indication"], r["study"]))
    labels = [f"{r['study']} ({r['indication']}, n={r['n']})" for r in rows]
    ors = [r["OR"] for r in rows]
    lo = [r["OR_lo"] for r in rows]
    hi = [r["OR_hi"] for r in rows]
    n = len(rows) + (1 if summary.get("k") else 0)
    fig_h = max(3.5, 0.42 * n + 1.6)
    fig, ax = plt.subplots(figsize=(8.2, fig_h))
    ys = list(range(len(rows), 0, -1))
    ax.axvline(1.0, color="0.5", lw=1, ls="--")
    for y, o, l, h in zip(ys, ors, lo, hi):
        ax.plot([l, h], [y, y], color="0.15", lw=1.2)
        ax.plot(o, y, "s", color="0.15", ms=5)
    if summary.get("k"):
        ax.plot([summary["OR_lo"], summary["OR_hi"]], [0, 0], color="#8B1E3F", lw=1.6)
        ax.plot(summary["OR"], 0, "D", color="#8B1E3F", ms=7)
        labels = labels + [
            f"RE meta k={summary['k']}  OR={summary['OR']:.2f} "
            f"({summary['OR_lo']:.2f}–{summary['OR_hi']:.2f})"
        ]
        ys = ys + [0]
    ax.set_yticks(ys)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Odds ratio (CLDN4-high vs low; response = R)")
    ax.set_xscale("log")
    xmin = min(lo + ([summary["OR_lo"]] if summary.get("k") else []))
    xmax = max(hi + ([summary["OR_hi"]] if summary.get("k") else []))
    ax.set_xlim(max(0.02, xmin / 1.4), min(20, xmax * 1.4))
    ax.set_title(title, fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    df = pd.read_csv(DER / "cldn4_clinical.tsv", sep="\t")
    df["recist_bin"] = recist_or(df["recist"])

    # Primary: PredictIO binary response (R vs NR).
    # Sensitivity: RECIST CR/PR vs PD (SD dropped).
    all_rows = []
    for study, g in df.groupby("study", sort=False):
        for split in ("median", "tertile", "quartile"):
            for rcol in ("response", "recist_bin"):
                row = analyze_cohort(g, rcol, split)
                if row:
                    all_rows.append(row)
    res = pd.DataFrame(all_rows)
    res.to_csv(OUT / "per_cohort_or.tsv", sep="\t", index=False)

    def pick(frame: pd.DataFrame, **kw) -> list[dict]:
        q = frame
        for k, v in kw.items():
            if k == "in_scope":
                q = q[q["in_scope"] == v]
            elif k == "io_class":
                q = q[q["io_class"] == v]
            elif k == "split":
                q = q[q["split"] == v]
            elif k == "response_def":
                q = q[q["response_def"] == v]
            elif k == "drop_low_var":
                if v:
                    q = q[~q["low_dynamic_range"]]
            elif k == "min_n":
                q = q[q["n"] >= v]
        return q.to_dict("records")

    analyses = {
        "primary_pd1_inscope_median_predictio": {
            "title": "Primary: in-scope PD-1/PD-L1, median split, PredictIO R/NR",
            "rows": pick(
                res,
                in_scope=True,
                io_class="PD-1/PD-L1",
                split="median",
                response_def="response",
                min_n=10,
            ),
        },
        "pd1_inscope_drop_lowvar": {
            "title": "PD-1/PD-L1 in-scope, drop IQR<0.5 (CLDN4 nearly constant)",
            "rows": pick(
                res,
                in_scope=True,
                io_class="PD-1/PD-L1",
                split="median",
                response_def="response",
                drop_low_var=True,
                min_n=10,
            ),
        },
        "pd1_inscope_recist_crpr_vs_pd": {
            "title": "PD-1/PD-L1 in-scope, median, RECIST CR/PR vs PD (SD dropped)",
            "rows": pick(
                res,
                in_scope=True,
                io_class="PD-1/PD-L1",
                split="median",
                response_def="recist_bin",
                min_n=10,
            ),
        },
        "pd1_inscope_tertile": {
            "title": "PD-1/PD-L1 in-scope, tertile high vs low",
            "rows": pick(
                res,
                in_scope=True,
                io_class="PD-1/PD-L1",
                split="tertile",
                response_def="response",
                min_n=8,
            ),
        },
        "pd1_plus_ctla4_inscope": {
            "title": "In-scope PD-1/PD-L1 + CTLA-4, median, PredictIO R/NR",
            "rows": pick(
                res,
                in_scope=True,
                split="median",
                response_def="response",
                min_n=10,
            ),
        },
        "all_studies_with_cldn4": {
            "title": "All studies with CLDN4 (includes gastric/pancreas), median",
            "rows": pick(res, split="median", response_def="response", min_n=10),
        },
        "epithelial_only": {
            "title": "Epithelial in-scope (NSCLC + RCC + urothelial; drop melanoma)",
            "rows": [
                r
                for r in pick(
                    res,
                    in_scope=True,
                    io_class="PD-1/PD-L1",
                    split="median",
                    response_def="response",
                    min_n=10,
                )
                if r["indication"] in ("NSCLC", "RCC", "urothelial")
            ],
        },
    }

    summaries = {}
    for key, spec in analyses.items():
        rows = spec["rows"]
        smry = re_meta(rows)
        smry["title"] = spec["title"]
        smry["studies"] = [r["study"] for r in rows]
        summaries[key] = smry
        if rows:
            forest(rows, smry, spec["title"], OUT / f"forest_{key}.png")

        # leave-one-out for primary
        if key == "primary_pd1_inscope_median_predictio" and len(rows) >= 3:
            loo = []
            for i, dropped in enumerate(rows):
                rest = rows[:i] + rows[i + 1 :]
                s = re_meta(rest)
                s["dropped"] = dropped["study"]
                loo.append(s)
            pd.DataFrame(loo).to_csv(OUT / "leave_one_out_primary.tsv", sep="\t", index=False)

    # indication-stratified primary
    prim = analyses["primary_pd1_inscope_median_predictio"]["rows"]
    by_ind = {}
    for ind in ("NSCLC", "melanoma", "RCC", "urothelial", "GBM"):
        sub = [r for r in prim if r["indication"] == ind]
        by_ind[ind] = re_meta(sub)
        by_ind[ind]["studies"] = [r["study"] for r in sub]
    summaries["by_indication"] = by_ind

    # claim comparison
    claim = {
        "claimed_k": 11,
        "claimed_indications": ["GBM", "NSCLC", "melanoma", "RCC", "urothelial"],
        "claimed_OR": 0.42,
        "named_public_pd1_cohorts_with_cldn4": summaries["primary_pd1_inscope_median_predictio"].get("studies", []),
        "reconstructed_k": summaries["primary_pd1_inscope_median_predictio"].get("k", 0),
        "reconstructed_OR": summaries["primary_pd1_inscope_median_predictio"].get("OR"),
        "reconstructed_OR_CI": [
            summaries["primary_pd1_inscope_median_predictio"].get("OR_lo"),
            summaries["primary_pd1_inscope_median_predictio"].get("OR_hi"),
        ],
        "gbm_public_processed_paired": False,
        "gbm_note": (
            "No GBM cohort is present in the public ORCESTRA/PredictIO ICB objects. "
            "Zhao et al. Nat Med 2019 (PMID 30742119) deposited raw reads in SRA PRJNA482620 "
            "and stated processed data are available on request — not used. "
            "GSE121810 series matrix is metadata-only (3.2 KB). "
            "GSE154795 is scRNA-seq of immune cells, not bulk tumour CLDN4. "
            "GSE264695 is neoadjuvant pembrolizumab bulk RNA-seq without a public "
            "paired RECIST R/NR table in the GEO supplement we could use here."
        ),
        "published_source_of_OR_0_42": None,
        "verdict": None,
    }
    rec_or = claim["reconstructed_OR"]
    rec_lo = claim["reconstructed_OR_CI"][0]
    rec_hi = claim["reconstructed_OR_CI"][1]
    named = claim["named_public_pd1_cohorts_with_cldn4"]
    claim["named_set_includes_gbm"] = False
    claim["named_set_is_the_claimed_11"] = False
    if rec_or is None:
        claim["verdict"] = "NOT RECONSTRUCTED: no eligible public cohorts."
    else:
        contains = rec_lo <= 0.42 <= rec_hi
        close = abs(np.log(rec_or) - np.log(0.42)) < np.log(1.25)
        claim["claim_OR_inside_reconstructed_CI"] = bool(contains)
        claim["reconstructed_OR_within_25pct_of_claim"] = bool(close)
        claim["verdict"] = (
            f"NOT REPRODUCED. The claim specifies 11 cohorts spanning GBM/NSCLC/"
            f"melanoma/RCC/urothelial with CLDN4-high OR=0.42. No published source "
            f"naming those 11 was found, and GBM is not reconstructable from public "
            f"processed expression+response files. Named public PD-1/PD-L1 cohorts "
            f"with usable CLDN4 + binary response: k={claim['reconstructed_k']} "
            f"({', '.join(named)}). Pooled DerSimonian–Laird OR = {rec_or:.2f} "
            f"(95% CI {rec_lo:.2f}–{rec_hi:.2f}, p={summaries['primary_pd1_inscope_median_predictio'].get('p'):.3g}). "
            f"The claimed OR 0.42 "
            f"{'falls inside' if contains else 'is outside'} this CI. "
            f"Point estimate is {'in the claimed direction (OR<1)' if rec_or < 1 else 'in the opposite direction (OR>1)'}."
        )

    (OUT / "meta_summaries.json").write_text(json.dumps(summaries, indent=2, default=str))
    (OUT / "claim_verdict.json").write_text(json.dumps(claim, indent=2, default=str))

    # compact primary table
    prim_df = pd.DataFrame(prim)
    if len(prim_df):
        cols = [
            "study",
            "indication",
            "n",
            "n_R",
            "n_NR",
            "n_high_R",
            "n_high_NR",
            "n_low_R",
            "n_low_NR",
            "OR",
            "OR_lo",
            "OR_hi",
            "p",
            "continuity_correction",
            "cldn4_iqr",
            "cldn4_frac_floor",
            "low_dynamic_range",
            "U_p",
            "cliffs_delta",
            "logit_OR_per_SD",
            "logit_p",
        ]
        prim_df[cols].to_csv(OUT / "primary_cohort_table.tsv", sep="\t", index=False)

    print(json.dumps({k: {kk: summaries[k].get(kk) for kk in ("k", "OR", "OR_lo", "OR_hi", "p", "I2", "studies")}
                      if isinstance(summaries[k], dict) and "k" in summaries[k] else summaries[k]
                      for k in summaries}, indent=2, default=str))
    print("VERDICT:", claim["verdict"])


if __name__ == "__main__":
    main()
