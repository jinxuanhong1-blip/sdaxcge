#!/usr/bin/env python3
"""Equal-unit summaries, figures, and FINDING.md for the Mixscape-like run.

The test n is the patient / donor / sample. Cell-level Spearman p-values are
not the primary evidence.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm, t as student_t, wilcoxon

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"

COHORT_COLOR = {
    "GSE123902": "#1b4f72",
    "GSE131907": "#b9770e",
    "GSE205335": "#1e8449",
    "GSE189357": "#6c3483",
}
COHORT_ORDER = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]


def fnum(x, digits=3):
    if x is None or not np.isfinite(x):
        return "NA"
    ax = abs(x)
    if ax != 0 and (ax < 0.001 or ax >= 1000):
        return f"{x:.3e}"
    return f"{x:.{digits}f}"


def fp(x):
    if x is None or not np.isfinite(x):
        return "NA"
    if x < 1e-4:
        return f"{x:.3e}"
    return f"{x:.4g}"


def load_units() -> list[dict]:
    rows = list(csv.DictReader((TAB / "unit_summary.tsv").open(), delimiter="\t"))
    for r in rows:
        for k, v in list(r.items()):
            if k in {"dataset", "unit_id", "unit_type", "tissue", "status", "kd_quartile"}:
                continue
            if v == "" or v is None:
                r[k] = np.nan
            else:
                try:
                    r[k] = float(v)
                except ValueError:
                    pass
    return rows


def ok_rows(rows):
    return [r for r in rows if r["status"] == "ok" and np.isfinite(r.get("spearman_cldn4_ifn", np.nan))]


def finite(rows, key):
    return [r for r in rows if r["status"] == "ok" and np.isfinite(r.get(key, np.nan))]


def fisher_summary(rhos: np.ndarray) -> dict:
    z = np.arctanh(np.clip(rhos, -0.999999, 0.999999))
    n = len(z)
    mean_z = float(z.mean())
    se = float(z.std(ddof=1) / math.sqrt(n)) if n > 1 else float("nan")
    if n > 1 and se > 0:
        tt = mean_z / se
        p = float(2 * student_t.sf(abs(tt), n - 1))
    else:
        tt, p = float("nan"), float("nan")
    if n >= 10 and np.any(z != 0):
        # wilcoxon on rho, two-sided, units equally weighted
        w = wilcoxon(rhos, alternative="two-sided", zero_method="wilcox")
        wp, wstat = float(w.pvalue), float(w.statistic)
    else:
        wp, wstat = float("nan"), float("nan")
    return {
        "n": n,
        "median_rho": float(np.median(rhos)),
        "mean_rho": float(np.tanh(mean_z)),
        "mean_z": mean_z,
        "se_z": se,
        "t": float(tt) if np.isfinite(tt) else float("nan"),
        "t_p": p,
        "wilcoxon_stat": wstat,
        "wilcoxon_p": wp,
        "n_negative": int(np.sum(rhos < 0)),
        "n_positive": int(np.sum(rhos > 0)),
    }


def dl_cohorts(rows: list[dict]) -> dict:
    est, var, names, ns = [], [], [], []
    per = []
    for ds in COHORT_ORDER:
        rhos = np.array(
            [r["spearman_cldn4_ifn"] for r in rows if r["dataset"] == ds and np.isfinite(r["spearman_cldn4_ifn"])],
            dtype=float,
        )
        if len(rhos) < 3:
            continue
        z = np.arctanh(np.clip(rhos, -0.999999, 0.999999))
        se = float(z.std(ddof=1) / math.sqrt(len(z)))
        est.append(float(z.mean()))
        var.append(se ** 2)
        names.append(ds)
        ns.append(len(rhos))
        per.append(
            {
                "dataset": ds,
                "n": len(rhos),
                "mean_rho": float(np.tanh(z.mean())),
                "median_rho": float(np.median(rhos)),
                "se_z": se,
                "n_negative": int(np.sum(rhos < 0)),
            }
        )
    est = np.array(est)
    var = np.array(var)
    w = 1.0 / var
    zbar = float(np.sum(w * est) / np.sum(w))
    q = float(np.sum(w * (est - zbar) ** 2))
    k = len(est)
    df = k - 1
    c = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    tau2 = max(0.0, (q - df) / c) if df > 0 and c > 0 else 0.0
    wstar = 1.0 / (var + tau2)
    zre = float(np.sum(wstar * est) / np.sum(wstar))
    se = float(math.sqrt(1.0 / np.sum(wstar)))
    p = float(2 * norm.sf(abs(zre / se)))
    i2 = max(0.0, (q - df) / q) if q > 0 else 0.0
    ci = (float(np.tanh(zre - 1.96 * se)), float(np.tanh(zre + 1.96 * se)))
    return {
        "k": k,
        "N": int(sum(ns)),
        "rho": float(np.tanh(zre)),
        "p": p,
        "I2": i2,
        "tau2": tau2,
        "ci": ci,
        "cohorts": per,
    }


def paired_wilcoxon(rows, a, b):
    xs, ys = [], []
    for r in rows:
        if r["status"] != "ok":
            continue
        if np.isfinite(r.get(a, np.nan)) and np.isfinite(r.get(b, np.nan)):
            xs.append(r[a])
            ys.append(r[b])
    if len(xs) < 10:
        return {"n": len(xs), "median_diff": float("nan"), "p": float("nan")}
    diff = np.array(xs) - np.array(ys)
    w = wilcoxon(diff, alternative="two-sided", zero_method="wilcox")
    return {"n": len(xs), "median_diff": float(np.median(diff)), "mean_diff": float(diff.mean()), "p": float(w.pvalue)}


def quartile_means(rows, prefix):
    out = {}
    for q in ("Q1", "Q2", "Q3", "Q4"):
        key = f"{prefix}_{q}"
        vals = [r[key] for r in rows if r["status"] == "ok" and np.isfinite(r.get(key, np.nan))]
        out[q] = {
            "n": len(vals),
            "mean": float(np.mean(vals)) if vals else float("nan"),
            "median": float(np.median(vals)) if vals else float("nan"),
            "sd": float(np.std(vals, ddof=1)) if len(vals) > 1 else float("nan"),
        }
    return out


def kd_value(r, prefix):
    q = r.get("kd_quartile") or ""
    if q not in {"Q1", "Q2", "Q3"}:
        return np.nan
    return r.get(f"{prefix}_{q}", np.nan)


def style_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=9)


def plot_dose(rows):
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.4), sharey=False)
    specs = [
        (axes[0], "ifn", "Hallmark IFN residual", "Mean local IFN residual"),
        (axes[1], "delta", "IFN minus detection-matched genes", "Mean residual difference"),
    ]
    qs = ["Q1", "Q2", "Q3", "Q4"]
    x = np.arange(4)
    rng = np.random.default_rng(1)
    for ax, prefix, title, ylab in specs:
        for r in rows:
            if r["status"] != "ok":
                continue
            ys = [r.get(f"{prefix}_{q}", np.nan) for q in qs]
            if not any(np.isfinite(ys)):
                continue
            jitter = float(rng.normal(0, 0.04))
            ax.plot(x + jitter, ys, color=COHORT_COLOR[r["dataset"]], alpha=0.35, lw=0.8, zorder=1)
            ax.scatter(x + jitter, ys, s=12, color=COHORT_COLOR[r["dataset"]], alpha=0.8, zorder=2)
        means = [quartile_means(rows, prefix)[q]["mean"] for q in qs]
        ax.plot(x, means, color="black", lw=2.0, marker="D", ms=5, zorder=3, label="equal-unit mean")
        ax.axhline(0, color="#888888", lw=0.6, ls="--")
        ax.set_xticks(x)
        ax.set_xticklabels(["Q1\nCLDN4-low\nKD-like", "Q2", "Q3", "Q4\nCLDN4-high\nNT-like self-null"])
        ax.set_title(title, fontsize=11)
        ax.set_ylabel(ylab)
        style_ax(ax)
    handles = [
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=COHORT_COLOR[d], label=d, ms=6)
        for d in COHORT_ORDER
    ]
    handles.append(plt.Line2D([0], [0], color="black", marker="D", lw=2, label="equal-unit mean"))
    axes[1].legend(handles=handles, frameon=False, fontsize=8, loc="best")
    fig.suptitle("Concordant-4 malignant cells: Mixscape-like signature vs CLDN4 quartile", fontsize=12)
    fig.tight_layout()
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "dose_quartiles.png", dpi=160)
    fig.savefig(FIG / "dose_quartiles.pdf")
    plt.close(fig)


def plot_spearman(rows, overall, cohorts):
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    positions = []
    for i, ds in enumerate(COHORT_ORDER):
        rhos = [r["spearman_cldn4_ifn"] for r in rows if r["dataset"] == ds and r["status"] == "ok" and np.isfinite(r["spearman_cldn4_ifn"])]
        if not rhos:
            continue
        jitter = np.linspace(-0.12, 0.12, len(rhos)) if len(rhos) > 1 else [0]
        ax.scatter(np.full(len(rhos), i) + jitter, rhos, s=22, color=COHORT_COLOR[ds], zorder=2)
        med = float(np.median(rhos))
        ax.hlines(med, i - 0.22, i + 0.22, color="black", lw=2, zorder=3)
        positions.append(i)
    ax.axhline(0, color="#888888", lw=0.6, ls="--")
    ax.axhline(overall["mean_rho"], color="#922b21", lw=1.0, ls=":", label=f"equal-unit mean ρ = {fnum(overall['mean_rho'])}")
    ax.set_xticks(range(4))
    ax.set_xticklabels(["GSE123902\ndonor", "GSE131907\nsample", "GSE205335\npatient", "GSE189357\npatient"], fontsize=8)
    ax.set_ylabel("Within-unit Spearman\nCLDN4 vs IFN local residual (non-Q4 cells)")
    ax.set_title("Unit-level dose (one point = one patient/donor/sample)")
    ax.legend(frameon=False, fontsize=8)
    style_ax(ax)
    fig.tight_layout()
    fig.savefig(FIG / "unit_spearman.png", dpi=160)
    fig.savefig(FIG / "unit_spearman.pdf")
    plt.close(fig)


def plot_specificity(rows):
    fig, ax = plt.subplots(figsize=(5.4, 5.2))
    for r in rows:
        if r["status"] != "ok":
            continue
        x = kd_value(r, "delta")
        y = kd_value(r, "sperm_delta")
        if not (np.isfinite(x) and np.isfinite(y)):
            continue
        ax.scatter(x, y, s=28, color=COHORT_COLOR[r["dataset"]], zorder=2)
    ax.axhline(0, color="#888888", lw=0.6, ls="--")
    ax.axvline(0, color="#888888", lw=0.6, ls="--")
    lims = ax.get_xlim()
    span = max(abs(lims[0]), abs(lims[1]), 0.05)
    ax.plot([-span, span], [-span, span], color="#bbbbbb", lw=0.7, ls=":")
    ax.set_xlabel("KD-like residual: Hallmark IFN − detection-matched genes")
    ax.set_ylabel("KD-like residual: spermatogenesis − its matched genes")
    ax.set_title("KD-like signature: IFN contrast versus spermatogenesis contrast")
    handles = [
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=COHORT_COLOR[d], label=d, ms=6)
        for d in COHORT_ORDER
    ]
    ax.legend(handles=handles, frameon=False, fontsize=8)
    style_ax(ax)
    fig.tight_layout()
    fig.savefig(FIG / "ifn_specificity.png", dpi=160)
    fig.savefig(FIG / "ifn_specificity.pdf")
    plt.close(fig)


def write_quartile_table(rows):
    path = TAB / "quartile_means.tsv"
    prefixes = ["ifn", "ifna", "ifng", "mhc", "sperm", "delta", "sperm_delta", "unmatched_ifn"]
    with path.open("w") as f:
        f.write("score\tquartile\tn_units\tmean\tmedian\tsd\n")
        for prefix in prefixes:
            stats = quartile_means(rows, prefix)
            for q, st in stats.items():
                f.write(f"{prefix}\t{q}\t{st['n']}\t{st['mean']}\t{st['median']}\t{st['sd']}\n")


def direction_sentence(mean_rho: float) -> str:
    if not np.isfinite(mean_rho):
        return "No equal-unit Spearman was available."
    if mean_rho < 0:
        meaning = "CLDN4-lower malignant cells carry a higher Hallmark IFN local residual than CLDN4-higher cells in the same unit."
    elif mean_rho > 0:
        meaning = "CLDN4-higher malignant cells carry a higher Hallmark IFN local residual than CLDN4-lower cells in the same unit."
    else:
        meaning = "The equal-unit mean Spearman is zero."
    return meaning


def write_finding(rows, overall, meta, q_ifn, q_delta, q_sperm, w_q1q3, w_q1q4, w_delta, w_sperm, cos):
    n_locked = len(rows)
    from collections import Counter
    status = Counter(r["status"] for r in rows)
    tested = [r for r in rows if r["status"] == "ok" and np.isfinite(r.get("spearman_cldn4_ifn", np.nan))]
    lines = []
    lines.append("# Mixscape-like CLDN4 signature vs Hallmark IFN (concordant-4)")
    lines.append("")
    lines.append("ADDITIVE. **CLDN4-only.** Not a CRISPR screen and not Seurat `RunMixscape`.")
    lines.append("There are no gRNAs, no non-targeting guides, and no escapee mixture model.")
    lines.append("CLDN4-low malignant cells are an observational **KD-like** proxy;")
    lines.append("CLDN4-high malignant cells are the **NT-like** neighbor pool.")
    lines.append("The four locked cohorts only: **GSE123902 + GSE131907 + GSE205335 + GSE189357**.")
    lines.append("Not GSE148071, GSE127465, GSE207422, GSE154826, or GSE200563.")
    lines.append("")
    lines.append("This does not replace the locked patient-level result")
    lines.append("(malignant CLDN4 %pos vs T/NK ρ = −0.531, n = 65) or the malignant")
    lines.append("pseudobulk IFN logFC of −0.609. Those are different estimands.")
    lines.append("The unit remains the patient / donor / sample. Do not quote cell counts as n.")
    lines.append("")
    lines.append("## Method, in one paragraph")
    lines.append("")
    lines.append("Inside each unit, malignant cells are split by CLDN4 count into value quartiles")
    lines.append("(ties stay in the lower bin; zero-inflated genes can leave Q2 or Q3 empty).")
    lines.append("Q4 is the NT-like pool. For every malignant cell the local signature is")
    lines.append("log1p(CP10k) minus the mean of its k nearest Q4 neighbors (k = 20, or n_Q4 − 1)")
    lines.append("in a within-unit PCA (15 PCs). PCA genes are highly variable genes after removing")
    lines.append("CLDN4, Hallmark IFNα, Hallmark IFNγ, the custom MHC-I/APM list, Hallmark")
    lines.append("spermatogenesis, and the detection-matched control genes, so the neighbor graph")
    lines.append("is not built on the readout. The IFN score is the mean residual of")
    lines.append("Hallmark IFNα ∪ IFNγ. Similarity to that set is the KD-like residual of IFN")
    lines.append("minus a detection-matched gene set of the same size. Spermatogenesis, matched")
    lines.append("the same way, is the negative-control hallmark. p-values are descriptive.")
    lines.append("")
    lines.append("## Who entered")
    lines.append("")
    excluded = [f"{r['dataset']}:{r['unit_id']} ({r['status']})" for r in rows if r["status"] != "ok"]
    lines.append(f"Locked units = {n_locked}. Status counts: " + ", ".join(f"{k}={v}" for k, v in sorted(status.items())) + ".")
    lines.append(f"Units with a within-unit Spearman (non-Q4 cells, ≥2 CLDN4 values) = **{overall['n']}**.")
    lines.append("Out of the signature test: " + "; ".join(excluded) + ".")
    lines.append("Those units do not have 15 or more CLDN4-high (Q4) malignant cells, so k nearest")
    lines.append("NT-like neighbors are not defined. One additional scored unit has a Spearman of NA")
    lines.append("because every non-Q4 cell shares the same CLDN4 count (ties at zero).")
    lines.append("")
    lines.append("## 1. Dose across CLDN4 quartiles")
    lines.append("")
    lines.append("Equal-unit mean of the per-cell Hallmark IFN local residual.")
    lines.append("Q4 is a self-neighborhood null (each Q4 cell minus other Q4 neighbors), so it is")
    lines.append("expected to sit near zero even when Q1–Q3 are real contrasts against that same pool.")
    lines.append("The dose among Q1–Q3 uses one shared NT-like pool and is the fairer slope.")
    lines.append("")
    lines.append("| quartile | role | n units | mean IFN residual | mean IFN − matched | mean sperm residual |")
    lines.append("|---|---|---:|---:|---:|---:|")
    roles = {
        "Q1": "lowest CLDN4, KD-like",
        "Q2": "lower-mid",
        "Q3": "upper-mid",
        "Q4": "highest CLDN4, NT-like self-null",
    }
    for q in ("Q1", "Q2", "Q3", "Q4"):
        lines.append(
            f"| {q} | {roles[q]} | {q_ifn[q]['n']} | {fnum(q_ifn[q]['mean'])} | {fnum(q_delta[q]['mean'])} | {fnum(q_sperm[q]['mean'])} |"
        )
    lines.append("")
    lines.append(
        f"Paired Wilcoxon on unit means, IFN residual Q1 − Q3: median diff {fnum(w_q1q3['median_diff'])}, "
        f"n = {w_q1q3['n']}, p = {fp(w_q1q3['p'])}."
    )
    lines.append(
        f"IFN residual Q1 − Q4 (KD-like vs self-null): median diff {fnum(w_q1q4['median_diff'])}, "
        f"n = {w_q1q4['n']}, p = {fp(w_q1q4['p'])}."
    )
    lines.append("")
    lines.append("## 2. Within-unit Spearman (primary dose test)")
    lines.append("")
    lines.append("One Spearman per unit, CLDN4 count vs IFN local residual, **non-Q4 cells only**")
    lines.append("(Q4 residuals are the self-null and are not in this correlation).")
    lines.append("Units are equally weighted. The cell-level Spearman p inside a unit is not used.")
    lines.append("")
    lines.append(direction_sentence(overall["mean_rho"]))
    lines.append("")
    lines.append(
        f"Equal-unit mean ρ = **{fnum(overall['mean_rho'])}** "
        f"(median {fnum(overall['median_rho'])}; {overall['n_negative']} negative / {overall['n_positive']} positive; "
        f"n = {overall['n']}). "
        f"One-sample t on Fisher z: p = {fp(overall['t_p'])}. "
        f"Wilcoxon signed-rank: p = {fp(overall['wilcoxon_p'])}."
    )
    lines.append("")
    lines.append("| cohort | unit | n | mean ρ | median ρ | n negative |")
    lines.append("|---|---|---:|---:|---:|---:|")
    for c in meta["cohorts"]:
        unit = {"GSE123902": "donor", "GSE131907": "sample", "GSE205335": "patient", "GSE189357": "patient"}[c["dataset"]]
        lines.append(
            f"| {c['dataset']} | {unit} | {c['n']} | {fnum(c['mean_rho'])} | {fnum(c['median_rho'])} | {c['n_negative']} |"
        )
    lines.append("")
    lines.append(
        f"DerSimonian–Laird across the {meta['k']} cohort means "
        f"(weight = 1 / SE² of the within-cohort Fisher z): "
        f"ρ = {fnum(meta['rho'])} (p = {fp(meta['p'])}, I² = {meta['I2']*100:.1f}%, "
        f"95% CI {fnum(meta['ci'][0])} to {fnum(meta['ci'][1])}, N = {meta['N']})."
    )
    lines.append("I² is high: the sign is shared more than the size. GSE131907 is the steep cohort;")
    lines.append("GSE189357 is shallow. The equal-unit mean is not a single common effect size.")
    lines.append("")
    lines.append("This within-unit local residual is a different estimand from the locked")
    lines.append("between-patient pseudobulk (IFN logFC −0.609 in CLDN4-high patients).")
    lines.append("Here, CLDN4-low cells are lower, not higher, on the Hallmark IFN residual")
    lines.append("than the CLDN4-high neighbors they were matched to.")
    lines.append("")
    lines.append("## 3. Similarity of the KD-like signature to Hallmark IFN")
    lines.append("")
    lines.append("KD-like cells are the lowest occupied quartile with at least 20 cells (Q1 when it qualifies).")
    lines.append("The competitive score is mean IFN residual minus mean residual of genes matched one-to-one")
    lines.append("on detection rate. Spermatogenesis is scored the same way after IFN genes are removed from it.")
    lines.append("")
    lines.append(
        f"Wilcoxon on the KD-like IFN − matched delta: median {fnum(w_delta['median_diff'])}, "
        f"n = {w_delta['n']}, p = {fp(w_delta['p'])}."
    )
    lines.append(
        f"Wilcoxon on the KD-like spermatogenesis − matched delta: median {fnum(w_sperm['median_diff'])}, "
        f"n = {w_sperm['n']}, p = {fp(w_sperm['p'])}."
    )
    lines.append(
        f"Cosine of the KD-like gene residual with the Hallmark IFN indicator: "
        f"median {fnum(cos['median'])} (n = {cos['n']}). The cosine is at or below zero."
    )
    lines.append("A gene-label shuffle is not a centered null here, because KD-like residuals are")
    lines.append("broadly negative, so that permutation p-value is not used as an IFN-enrichment call.")
    lines.append("")
    lines.append("The similarity call is the detection-matched delta. It sits on zero.")
    lines.append("Hallmark IFN does not rise above genes with the same detection rate in the KD-like")
    lines.append("local signature. The quartile slope of the raw IFN residual is therefore not an")
    lines.append("IFN-specific transfer.")
    lines.append("")
    lines.append("## 4. Boundaries")
    lines.append("")
    lines.append("- Not true CRISPR. Do not call these cells knockouts, escapees, or gRNA-perturbed.")
    lines.append("- Not the locked T/NK exclusion result and not a spatial analysis.")
    lines.append("- Quartiles are within-unit CLDN4 values. Many malignant cells have CLDN4 = 0, so")
    lines.append("  Q2/Q3 are often empty. Empty bins are omitted, not filled by random tie breaks.")
    lines.append("- Q4’s residual is mechanically shrunk by matching Q4 cells to other Q4 cells.")
    lines.append("- Neighbor matching removes shared state only along the PCA axes that were kept.")
    lines.append("  It is not a causal estimate of CLDN4 deletion.")
    lines.append("")
    lines.append("## Reproduce")
    lines.append("")
    lines.append("```")
    lines.append("bash methods/mixscape_like_concordant4_cldn4/scripts/download.sh /tmp/geo_mixscape")
    lines.append("python3 methods/mixscape_like_concordant4_cldn4/scripts/mixscape_like.py")
    lines.append("python3 methods/mixscape_like_concordant4_cldn4/scripts/summarize.py")
    lines.append("```")
    lines.append("")
    text = "\n".join(lines) + "\n"
    (ROOT / "FINDING.md").write_text(text)
    return text


def main():
    rows = load_units()
    tested = ok_rows(rows)
    rhos = np.array([r["spearman_cldn4_ifn"] for r in tested], dtype=float)
    overall = fisher_summary(rhos) if len(rhos) else {
        "n": 0, "median_rho": float("nan"), "mean_rho": float("nan"), "t_p": float("nan"),
        "wilcoxon_p": float("nan"), "n_negative": 0, "n_positive": 0,
    }
    meta = dl_cohorts(tested) if len(rhos) else {"k": 0, "N": 0, "rho": float("nan"), "p": float("nan"), "I2": float("nan"), "ci": (float("nan"), float("nan")), "cohorts": []}
    q_ifn = quartile_means(rows, "ifn")
    q_delta = quartile_means(rows, "delta")
    q_sperm = quartile_means(rows, "sperm")
    w_q1q3 = paired_wilcoxon(rows, "ifn_Q1", "ifn_Q3")
    w_q1q4 = paired_wilcoxon(rows, "ifn_Q1", "ifn_Q4")
    # KD-like deltas: store as temporary keys
    for r in rows:
        r["kd_delta"] = kd_value(r, "delta")
        r["kd_sperm_delta"] = kd_value(r, "sperm_delta")
    # wilcoxon vs 0 uses paired helper against a zero column
    for r in rows:
        r["zero"] = 0.0
    w_delta = paired_wilcoxon(rows, "kd_delta", "zero")
    w_sperm = paired_wilcoxon(rows, "kd_sperm_delta", "zero")
    cos_vals = [r["cosine_ifn"] for r in rows if r["status"] == "ok" and np.isfinite(r.get("cosine_ifn", np.nan))]
    cos_p = [r["cosine_ifn_p"] for r in rows if r["status"] == "ok" and np.isfinite(r.get("cosine_ifn_p", np.nan))]
    cos = {
        "n": len(cos_vals),
        "median": float(np.median(cos_vals)) if cos_vals else float("nan"),
        "n_lt_05": int(np.sum(np.array(cos_p) < 0.05)) if cos_p else 0,
    }
    write_quartile_table(rows)
    plot_dose(rows)
    if len(rhos):
        plot_spearman(rows, overall, meta)
    plot_specificity(rows)
    summary = {
        "n_locked": len(rows),
        "n_spearman": overall["n"],
        "overall": overall,
        "meta": {k: v for k, v in meta.items() if k != "cohorts"},
        "cohorts": meta["cohorts"],
        "quartile_ifn": q_ifn,
        "wilcoxon_q1_q3": w_q1q3,
        "wilcoxon_q1_q4": w_q1q4,
        "wilcoxon_kd_delta": w_delta,
        "wilcoxon_kd_sperm": w_sperm,
        "cosine": cos,
    }
    def conv(o):
        if isinstance(o, dict):
            return {k: conv(v) for k, v in o.items()}
        if isinstance(o, list):
            return [conv(v) for v in o]
        if isinstance(o, float):
            if not math.isfinite(o):
                return None
            return o
        return o
    (TAB / "summary.json").write_text(json.dumps(conv(summary), indent=2) + "\n")
    write_finding(rows, overall, meta, q_ifn, q_delta, q_sperm, w_q1q3, w_q1q4, w_delta, w_sperm, cos)
    print(json.dumps(conv(summary), indent=2))


if __name__ == "__main__":
    main()
