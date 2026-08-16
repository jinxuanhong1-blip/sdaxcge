#!/usr/bin/env python3
"""TACSTD2 / CLDN4 vs ICI pathologic response and recurrence in leftover bulk tumour RNA-seq.

Primary leftover cohorts (whole transcriptome, genes present, outcome recovered from the
open-access trial papers because GEO itself only stores the randomisation arm):

  * GSE253564 — pre-treatment FFPE RNA-seq, neoadjuvant durvalumab +/- SBRT (n=32)
  * GSE248378 — post-treatment / surgical RNA-seq from the same trial (n=29)

Expression is the submitter FPKM matrix, analysed as log2(FPKM + 1). Tests are two-sided
Mann–Whitney U (binary MPR / recurrence), Spearman (continuous % viable-cancer reduction),
and log-rank + univariate Cox (PFS). Cliff's delta is the effect-size companion to the U test.
No values are imputed; a test is simply omitted when a label is missing.
"""
from __future__ import annotations

import csv
import gzip
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "opus_geo_leftover"
DATA = Path("/tmp/geo_dl/data")
GENES = ("TACSTD2", "CLDN4")
CONTROLS = ("CD8A", "EPCAM", "MKI67")


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """Cliff's δ: P(a > b) − P(a < b). Positive means group a is stochastically larger."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if a.size == 0 or b.size == 0:
        return float("nan")
    gt = 0
    lt = 0
    for x in a:
        gt += int(np.sum(x > b))
        lt += int(np.sum(x < b))
    return (gt - lt) / (a.size * b.size)


def read_fpkm(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", compression="gzip")
    gene_col = df.columns[0]
    df = df.set_index(gene_col)
    if "Entrez.ID" in df.columns:
        df = df.drop(columns=["Entrez.ID"])
    df.index = df.index.astype(str).str.strip()
    return df.apply(pd.to_numeric, errors="coerce")


def log2p1(s: pd.Series) -> pd.Series:
    return np.log2(s.astype(float) + 1.0)


def mwu_row(gene: str, values: pd.Series, labels: pd.Series, pos_name: str, neg_name: str) -> dict:
    aligned = pd.concat([values.rename("x"), labels.rename("y")], axis=1).dropna()
    pos = aligned.loc[aligned["y"] == 1, "x"].to_numpy()
    neg = aligned.loc[aligned["y"] == 0, "x"].to_numpy()
    if pos.size < 3 or neg.size < 3:
        return {
            "gene": gene,
            "test": f"Mann-Whitney U: {pos_name} vs {neg_name}",
            "n_pos": int(pos.size),
            "n_neg": int(neg.size),
            "median_pos": float(np.median(pos)) if pos.size else None,
            "median_neg": float(np.median(neg)) if neg.size else None,
            "U": None,
            "p": None,
            "cliffs_delta": None,
            "note": "too few labelled samples",
        }
    u, p = stats.mannwhitneyu(pos, neg, alternative="two-sided")
    return {
        "gene": gene,
        "test": f"Mann-Whitney U: {pos_name} vs {neg_name}",
        "n_pos": int(pos.size),
        "n_neg": int(neg.size),
        "median_pos": float(np.median(pos)),
        "median_neg": float(np.median(neg)),
        "iqr_pos": float(np.subtract(*np.percentile(pos, [75, 25]))),
        "iqr_neg": float(np.subtract(*np.percentile(neg, [75, 25]))),
        "U": float(u),
        "p": float(p),
        "cliffs_delta": float(cliffs_delta(pos, neg)),
        "note": "",
    }


def spearman_row(gene: str, values: pd.Series, other: pd.Series, label: str) -> dict:
    aligned = pd.concat([values.rename("x"), other.rename("y")], axis=1).dropna()
    if len(aligned) < 5:
        return {
            "gene": gene,
            "test": f"Spearman vs {label}",
            "n": int(len(aligned)),
            "rho": None,
            "p": None,
            "note": "too few labelled samples",
        }
    rho, p = stats.spearmanr(aligned["x"], aligned["y"])
    return {
        "gene": gene,
        "test": f"Spearman vs {label}",
        "n": int(len(aligned)),
        "rho": float(rho),
        "p": float(p),
        "note": "",
    }


def survival_row(gene: str, values: pd.Series, time: pd.Series, event: pd.Series) -> list[dict]:
    aligned = pd.concat(
        [values.rename("x"), time.rename("t"), event.rename("e")], axis=1
    ).dropna()
    aligned["t"] = pd.to_numeric(aligned["t"], errors="coerce")
    aligned["e"] = pd.to_numeric(aligned["e"], errors="coerce")
    aligned = aligned.dropna()
    out = []
    if len(aligned) < 8 or aligned["e"].sum() < 3:
        return [
            {
                "gene": gene,
                "test": "log-rank median split (PFS)",
                "n": int(len(aligned)),
                "n_events": int(aligned["e"].sum()) if len(aligned) else 0,
                "p": None,
                "note": "too few events",
            }
        ]
    med = float(aligned["x"].median())
    high = aligned["x"] >= med
    lr = logrank_test(
        aligned.loc[high, "t"],
        aligned.loc[~high, "t"],
        event_observed_A=aligned.loc[high, "e"],
        event_observed_B=aligned.loc[~high, "e"],
    )
    out.append(
        {
            "gene": gene,
            "test": "log-rank median split (PFS)",
            "n": int(len(aligned)),
            "n_high": int(high.sum()),
            "n_low": int((~high).sum()),
            "n_events": int(aligned["e"].sum()),
            "median_cut": med,
            "p": float(lr.p_value),
            "test_statistic": float(lr.test_statistic),
            "note": "high = expression >= cohort median",
        }
    )
    cox_df = aligned.rename(columns={"x": "expr", "t": "T", "e": "E"})
    cph = CoxPHFitter()
    try:
        cph.fit(cox_df, duration_col="T", event_col="E")
        s = cph.summary.loc["expr"]
        out.append(
            {
                "gene": gene,
                "test": "Cox PH, expression as continuous covariate (PFS)",
                "n": int(len(aligned)),
                "n_events": int(aligned["e"].sum()),
                "hr": float(s["exp(coef)"]),
                "hr_ci_low": float(s["exp(coef) lower 95%"]),
                "hr_ci_high": float(s["exp(coef) upper 95%"]),
                "p": float(s["p"]),
                "note": "HR per +1 log2(FPKM+1)",
            }
        )
    except Exception as exc:  # noqa: BLE001
        out.append(
            {
                "gene": gene,
                "test": "Cox PH, expression as continuous covariate (PFS)",
                "n": int(len(aligned)),
                "p": None,
                "note": f"Cox fit failed: {exc}",
            }
        )
    return out


def bh_fdr(pvals: list[float | None]) -> list[float | None]:
    pairs = [(i, p) for i, p in enumerate(pvals) if p is not None and not math.isnan(p)]
    out: list[float | None] = [None] * len(pvals)
    if not pairs:
        return out
    pairs.sort(key=lambda t: t[1])
    m = len(pairs)
    prev = 1.0
    ranked = []
    for rank, (i, p) in enumerate(pairs, start=1):
        ranked.append((i, min(1.0, p * m / rank)))
    for i, q in reversed(ranked):
        prev = min(prev, q)
        out[i] = prev
    # recompute standard BH (monotone from largest p)
    q = [None] * len(pvals)
    prev = 1.0
    for rank, (i, p) in zip(range(m, 0, -1), reversed(pairs)):
        prev = min(prev, p * m / rank)
        q[i] = prev
    return q


def analyse_series(acc: str, matrix: Path, clin: pd.DataFrame) -> tuple[pd.DataFrame, list[dict]]:
    expr = read_fpkm(matrix)
    sub = clin[clin["series"] == acc].copy()
    # Map expression columns onto the clinical sample_title.
    colmap = {c: c for c in expr.columns}
    missing = [t for t in sub["sample_title"] if t not in colmap]
    if missing:
        # GSE248378 titles match columns exactly in the downloaded file; fail loudly if not.
        raise RuntimeError(f"{acc}: expression columns missing for {missing}")
    sub = sub.set_index("sample_title")
    rows = []
    per_sample = sub.reset_index()[["series", "gsm", "sample_title", "patient", "arm_canonical",
                                    "arm_conflict", "mpr", "path_response_pct",
                                    "path_response_2group", "recurrence_event", "pfs_months",
                                    "histology", "egfr_status"]].copy()
    for gene in GENES + CONTROLS:
        if gene not in expr.index:
            rows.append({"series": acc, "gene": gene, "test": "presence", "note": "gene absent from matrix"})
            continue
        values = log2p1(expr.loc[gene, per_sample["sample_title"]])
        values.index = per_sample["sample_title"].to_numpy()
        per_sample[f"{gene}_log2fpkm1"] = values.to_numpy()
        mpr = pd.to_numeric(sub["mpr"], errors="coerce")
        mpr.index = sub.index
        rec = pd.to_numeric(sub["recurrence_event"], errors="coerce")
        rec.index = sub.index
        pct = pd.to_numeric(sub["path_response_pct"], errors="coerce")
        pct.index = sub.index
        # Pathology Response in Table S1 is signed reduction (more negative = deeper response).
        # Flip the sign so Spearman >0 means higher expression with deeper response.
        depth = -pct
        pfs = pd.to_numeric(sub["pfs_months"], errors="coerce")
        pfs.index = sub.index

        for rec_row in (
            mwu_row(gene, values, mpr, "MPR", "no-MPR"),
            mwu_row(gene, values, rec, "recurrence", "no-recurrence"),
            spearman_row(gene, values, depth, "pathologic response depth (−Pathology Response %)"),
        ):
            rec_row["series"] = acc
            rec_row["timepoint"] = "pre" if acc == "GSE253564" else "post"
            rows.append(rec_row)
        for rec_row in survival_row(gene, values, pfs, rec):
            rec_row["series"] = acc
            rec_row["timepoint"] = "pre" if acc == "GSE253564" else "post"
            rows.append(rec_row)

        # Sensitivity: dual-therapy arm only (Arm2), excluding the one GEO arm conflict.
        arm2 = sub["arm_canonical"].eq("arm2") & sub["arm_conflict"].fillna(0).astype(int).eq(0)
        if gene in GENES:
            r = mwu_row(gene, values[arm2.index][arm2], mpr[arm2], "MPR", "no-MPR")
            r.update({"series": acc, "timepoint": "pre" if acc == "GSE253564" else "post",
                      "test": "Mann-Whitney U: MPR vs no-MPR (Arm2 only)"})
            rows.append(r)

        # Sensitivity: is a TACSTD2/CLDN4–MPR gap independent of proliferation / epithelium?
        if gene in GENES and "MKI67" in expr.index and acc == "GSE253564":
            for cov_name in ("MKI67", "EPCAM"):
                if cov_name not in expr.index:
                    continue
                cov = log2p1(expr.loc[cov_name, values.index])
                aligned = pd.concat([values.rename("y"), cov.rename("x"), mpr.rename("m")], axis=1).dropna()
                if len(aligned) < 10:
                    continue
                slope, intercept, *_ = stats.linregress(aligned["x"], aligned["y"])
                resid = aligned["y"] - (intercept + slope * aligned["x"])
                r = mwu_row(gene, resid, aligned["m"], "MPR", "no-MPR")
                r.update(
                    {
                        "series": acc,
                        "timepoint": "pre",
                        "test": f"Mann-Whitney U on residual after {cov_name} OLS: MPR vs no-MPR",
                    }
                )
                rows.append(r)

    # Gene–gene Spearman as a tumour-epithelial sanity check.
    if all(g in expr.index for g in GENES):
        a = log2p1(expr.loc["TACSTD2", per_sample["sample_title"]])
        b = log2p1(expr.loc["CLDN4", per_sample["sample_title"]])
        rho, p = stats.spearmanr(a, b)
        rows.append(
            {
                "series": acc,
                "timepoint": "pre" if acc == "GSE253564" else "post",
                "gene": "TACSTD2 vs CLDN4",
                "test": "Spearman gene–gene (sanity)",
                "n": int(len(a)),
                "rho": float(rho),
                "p": float(p),
                "note": "expected positive in bulk tumour",
            }
        )
    return per_sample, rows


def paired_delta(pre: pd.DataFrame, post: pd.DataFrame) -> list[dict]:
    both = sorted(set(pre["patient"]) & set(post["patient"]))
    rows = []
    for gene in GENES:
        col = f"{gene}_log2fpkm1"
        if col not in pre.columns or col not in post.columns:
            continue
        p = pre.set_index("patient")
        q = post.set_index("patient")
        delta = q.loc[both, col] - p.loc[both, col]
        mpr = pd.to_numeric(p.loc[both, "mpr"], errors="coerce")
        rec = mwu_row(gene, delta, mpr, "MPR", "no-MPR")
        rec.update({"series": "paired GSE253564→GSE248378", "timepoint": "post−pre",
                    "test": "Mann-Whitney U on Δlog2(FPKM+1): MPR vs no-MPR"})
        rows.append(rec)
        depth = -pd.to_numeric(p.loc[both, "path_response_pct"], errors="coerce")
        rec2 = spearman_row(gene, delta, depth, "pathologic response depth")
        rec2.update({"series": "paired GSE253564→GSE248378", "timepoint": "post−pre"})
        rows.append(rec2)
    return rows


def main() -> int:
    clin = pd.read_csv(OUT / "clinical_annotation.csv")
    pre_expr, pre_rows = analyse_series(
        "GSE253564",
        DATA / "GSE253564" / "GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz",
        clin,
    )
    post_expr, post_rows = analyse_series(
        "GSE248378",
        DATA / "GSE248378" / "GSE248378_Durva_Post_FPKMs.txt.gz",
        clin,
    )
    rows = pre_rows + post_rows + paired_delta(pre_expr, post_expr)

    # BH-FDR on the primary leftover tests only (two genes × two series × MPR U-test).
    primary = [
        i
        for i, r in enumerate(rows)
        if r.get("gene") in GENES
        and r.get("test") == "Mann-Whitney U: MPR vs no-MPR"
        and r.get("series") in {"GSE253564", "GSE248378"}
    ]
    q = bh_fdr([rows[i].get("p") for i in primary])
    for i, qi in zip(primary, q):
        rows[i]["q_bh_primary_mpr"] = qi

    stats_df = pd.DataFrame(rows)
    stats_df.to_csv(OUT / "tacstd2_cldn4_stats.csv", index=False)
    pre_expr.to_csv(OUT / "GSE253564_sample_level.csv", index=False)
    post_expr.to_csv(OUT / "GSE248378_sample_level.csv", index=False)

    # Compact human-readable summary of the primary tests.
    keep = stats_df[stats_df["gene"].isin(GENES) | stats_df["gene"].eq("TACSTD2 vs CLDN4")]
    print(keep[["series", "timepoint", "gene", "test", "n_pos", "n_neg", "n", "median_pos",
                "median_neg", "cliffs_delta", "rho", "hr", "p", "q_bh_primary_mpr"]].to_string(index=False))
    print(f"\n[analyze] wrote {len(stats_df)} tests")
    _write_figures(pre_expr, post_expr)
    return 0


def _write_figures(pre: pd.DataFrame, post: pd.DataFrame) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.8), constrained_layout=True)
    for ax, gene in zip(axes, GENES):
        col = f"{gene}_log2fpkm1"
        a = pre.loc[pre["mpr"] == 1, col].dropna()
        b = pre.loc[pre["mpr"] == 0, col].dropna()
        ax.boxplot([b, a], tick_labels=["no-MPR", "MPR"], widths=0.55)
        rng = np.random.default_rng(0)
        ax.scatter(1 + rng.uniform(-0.08, 0.08, size=len(b)), b, s=14, color="#444444", alpha=0.7, zorder=3)
        ax.scatter(2 + rng.uniform(-0.08, 0.08, size=len(a)), a, s=14, color="#b33", alpha=0.7, zorder=3)
        ax.set_title(f"GSE253564 pre-treatment {gene}")
        ax.set_ylabel("log2(FPKM+1)")
    fig.savefig(OUT / "GSE253564_MPR_boxplots.png", dpi=140)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.8), constrained_layout=True)
    for ax, gene in zip(axes, GENES):
        col = f"{gene}_log2fpkm1"
        a = post.loc[post["recurrence_event"] == 1, col].dropna()
        b = post.loc[post["recurrence_event"] == 0, col].dropna()
        ax.boxplot([b, a], tick_labels=["no recurrence", "recurrence"], widths=0.55)
        rng = np.random.default_rng(1)
        ax.scatter(1 + rng.uniform(-0.08, 0.08, size=len(b)), b, s=14, color="#444444", alpha=0.7, zorder=3)
        ax.scatter(2 + rng.uniform(-0.08, 0.08, size=len(a)), a, s=14, color="#b33", alpha=0.7, zorder=3)
        ax.set_title(f"GSE248378 post-treatment {gene}")
        ax.set_ylabel("log2(FPKM+1)")
    fig.savefig(OUT / "GSE248378_recurrence_boxplots.png", dpi=140)
    plt.close(fig)
    print("[analyze] wrote boxplots")


if __name__ == "__main__":
    raise SystemExit(main())
