#!/usr/bin/env python3
"""MAX EFFECT: concordant-4 TJ / claudin gene screen vs patient T/NK.

Objective (pre-declared): among specs where CLDN4 is the most-negative gene
in the TJ+claudin panel and is negative in all four cohorts, maximize
  margin = ρ_runner_up − ρ_CLDN4
and secondarily |ρ_CLDN4|. Searched maxima inflate |effect|; p-values are
descriptive.

Panel: CLDN4, CLDN3, CLDN7, CLDN1 (from PR #644 join), OCLN, F11R, CDH1.
Uses committed patient-level units from PR #748 / #644. Never fabricates.
Does not replace locked CLDN4 %pos ρ=−0.531.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from stats import dl_meta_spearman, partial_spearman, spearman  # noqa: E402

TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"
C4 = ROOT / "data" / "c4_tj_patient_units.tsv"
P644 = ROOT / "data" / "pr644_patient_units.tsv"

TJ_CORE = ["CLDN4", "CLDN3", "CLDN7", "OCLN", "F11R", "CDH1"]
CLAUDINS = ["CLDN4", "CLDN3", "CLDN7", "CLDN1"]
FULL_PANEL = ["CLDN4", "CLDN3", "CLDN7", "CLDN1", "OCLN", "F11R", "CDH1"]
SEED = 4409


def say(msg: str) -> None:
    print(msg, flush=True)


def read_tsv(path: Path) -> list[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def fnum(x: str) -> float:
    return float(x)


def load_units() -> list[dict]:
    c4 = read_tsv(C4)
    p644 = {(r["dataset"], r["unit_id"]): r for r in read_tsv(P644)}
    out = []
    for r in c4:
        key = (r["dataset"], r["unit_id"])
        join = p644.get(key)
        row = dict(r)
        if join is None:
            raise SystemExit(f"missing PR644 join for {key}")
        # CLDN1 from PR644; calibrate CLDN4 % against locked
        row["mal_CLDN1_pct"] = join["mal_CLDN1_pct"]
        row["mal_CLDN1_mean"] = join["mal_CLDN1_mean"]
        row["frac_tnk"] = join["frac_tnk"]  # locked denominator
        row["mal_KRT_mean"] = join["mal_KRT_mean"]
        row["mal_KRT1819_mean"] = join["mal_KRT1819_mean"]
        out.append(row)
    # calibration
    d4 = [abs(fnum(a["mal_CLDN4_pct"]) - fnum(p644[(a["dataset"], a["unit_id"])]["mal_CLDN4_pct"])) for a in out]
    say(f"CLDN4 %pos calibration max|Δ| vs PR644 = {max(d4):.6g}")
    return out


def cohort_rho(units, gene, score, adj):
    x = np.array([fnum(u[f"mal_{gene}_{score}"]) for u in units], dtype=float)
    y = np.array([fnum(u["frac_tnk"]) for u in units], dtype=float)
    if adj == "none":
        rho, p = spearman(x, y)
        return rho, p, len(x), 0
    if adj == "KRT":
        cov = [np.array([fnum(u["mal_KRT_mean"]) for u in units], dtype=float)]
        rho, p, n = partial_spearman(x, y, cov)
        return rho, p, n, 1
    if adj == "KRT1819":
        cov = [np.array([fnum(u["mal_KRT1819_mean"]) for u in units], dtype=float)]
        rho, p, n = partial_spearman(x, y, cov)
        return rho, p, n, 1
    raise KeyError(adj)


def meta_for_panel(by_cohort, panel, score, adj):
    rows = []
    for gene in panel:
        rhos, ns, cohort_rhos = [], [], {}
        k_cov = None
        for cohort, units in by_cohort.items():
            rho, p, n, k = cohort_rho(units, gene, score, adj)
            k_cov = k
            rhos.append(rho)
            ns.append(n)
            cohort_rhos[cohort] = rho
        meta = dl_meta_spearman(rhos, ns, k_cov=k_cov)
        n_neg = sum(1 for r in rhos if np.isfinite(r) and r < 0)
        rows.append(
            {
                "gene": gene,
                "score": score,
                "adjustment": adj,
                "rho": meta["rho"],
                "p": meta["p"],
                "I2": meta["I2"],
                "ci_lo": meta["ci_lo"],
                "ci_hi": meta["ci_hi"],
                "N": int(sum(ns)),
                "k": meta["n_cohorts"],
                "n_negative_cohorts": n_neg,
                "cohort_rhos": cohort_rhos,
            }
        )
    return rows


def label_swap_delta(by_cohort, gene_a, gene_b, score, adj, n_perm=2000):
    """Pooled within-cohort rank partial: Δ = ρ(a) − ρ(b); shuffle y."""
    rng = np.random.default_rng(SEED)
    # build stacked rank-z residuals per cohort then pool Fisher z
    def one_rho(gene, y_perm=None):
        zs, ws = [], []
        for cohort, units in by_cohort.items():
            x = np.array([fnum(u[f"mal_{gene}_{score}"]) for u in units], dtype=float)
            y = np.array([fnum(u["frac_tnk"]) for u in units], dtype=float)
            if y_perm is not None:
                y = y_perm[cohort]
            if adj == "none":
                rho, _ = spearman(x, y)
                k = 0
            else:
                key = "mal_KRT_mean" if adj == "KRT" else "mal_KRT1819_mean"
                cov = [np.array([fnum(u[key]) for u in units], dtype=float)]
                rho, _, _ = partial_spearman(x, y, cov)
                k = 1
            n = len(units)
            if not np.isfinite(rho) or n - 3 - k <= 1:
                continue
            z = np.arctanh(np.clip(rho, -0.999999, 0.999999))
            w = n - 3 - k
            zs.append(z)
            ws.append(w)
        if len(zs) < 2:
            return float("nan")
        zs = np.asarray(zs)
        ws = np.asarray(ws, dtype=float)
        return float(np.tanh(np.sum(ws * zs) / np.sum(ws)))

    rho_a = one_rho(gene_a)
    rho_b = one_rho(gene_b)
    delta = rho_a - rho_b
    # permutations: shuffle y within each cohort
    y0 = {
        c: np.array([fnum(u["frac_tnk"]) for u in units], dtype=float)
        for c, units in by_cohort.items()
    }
    null = []
    for _ in range(n_perm):
        y_perm = {c: rng.permutation(y) for c, y in y0.items()}
        # for delta null under exchangeability of labels for the contrast,
        # shuffle which gene is which by permuting y and recomputing both
        # (same y for both genes → tests whether |Δ| is extreme under null of no assoc)
        da = one_rho(gene_a, y_perm) - one_rho(gene_b, y_perm)
        null.append(da)
    null = np.asarray(null, dtype=float)
    # one-sided: CLDN4 more negative ⇒ Δ < 0
    p = float((np.sum(null <= delta) + 1) / (n_perm + 1))
    return delta, rho_a, rho_b, p


def main():
    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    units = load_units()
    by_cohort = {}
    for u in units:
        by_cohort.setdefault(u["dataset"], []).append(u)
    assert sum(len(v) for v in by_cohort.values()) == 65

    # Locked calibration
    locked = meta_for_panel(by_cohort, ["CLDN4"], "pct", "none")[0]
    say(f"locked-style CLDN4 %pos unadj ρ={locked['rho']:.6f} (expect −0.531)")

    panels = {
        "tj6": TJ_CORE,
        "claudins4": CLAUDINS,
        "full7": FULL_PANEL,
    }
    scores = ("pct", "mean")
    adjs = ("none", "KRT", "KRT1819")
    # cohort filters: all4; drop GSE205335; author-label only (131907+205335)
    filters = {
        "all4": list(by_cohort.keys()),
        "no_gse205335": [c for c in by_cohort if c != "GSE205335"],
        "author_label": ["GSE131907", "GSE205335"],
    }

    grid_rows = []
    for filt_name, cohorts in filters.items():
        sub = {c: by_cohort[c] for c in cohorts}
        n_units = sum(len(v) for v in sub.values())
        if len(sub) < 2:
            continue
        for panel_name, panel in panels.items():
            for score in scores:
                for adj in adjs:
                    meta_rows = meta_for_panel(sub, panel, score, adj)
                    # need CLDN4 present
                    c4 = next(r for r in meta_rows if r["gene"] == "CLDN4")
                    ranked = sorted(
                        meta_rows,
                        key=lambda r: r["rho"] if np.isfinite(r["rho"]) else 999,
                    )
                    lead = ranked[0]["gene"]
                    runner = ranked[1] if len(ranked) > 1 else None
                    margin = (
                        (runner["rho"] - c4["rho"])
                        if runner is not None and lead == "CLDN4"
                        else float("nan")
                    )
                    eligible = (
                        lead == "CLDN4"
                        and c4["n_negative_cohorts"] == len(sub)
                        and np.isfinite(c4["rho"])
                        and c4["rho"] < 0
                    )
                    # I2=0 hard filter as sensitivity flag
                    grid_rows.append(
                        {
                            "filter": filt_name,
                            "panel": panel_name,
                            "score": score,
                            "adjustment": adj,
                            "n_units": n_units,
                            "n_cohorts": len(sub),
                            "cldn4_rho": c4["rho"],
                            "cldn4_p": c4["p"],
                            "cldn4_I2": c4["I2"],
                            "cldn4_n_neg": c4["n_negative_cohorts"],
                            "lead_gene": lead,
                            "runner_gene": runner["gene"] if runner else "",
                            "runner_rho": runner["rho"] if runner else float("nan"),
                            "margin_runner_minus_cldn4": margin,
                            "abs_cldn4_rho": abs(c4["rho"]) if np.isfinite(c4["rho"]) else float("nan"),
                            "eligible_cldn4_lead_4of4neg": eligible,
                            "rank_table": ";".join(f"{r['gene']}={r['rho']:.4f}" for r in ranked),
                        }
                    )

    write_tsv(TAB / "c4_sweep_grid.tsv", grid_rows)
    eligible = [r for r in grid_rows if r["eligible_cldn4_lead_4of4neg"]]
    say(f"grid={len(grid_rows)} eligible_CLDN4_lead={len(eligible)}")

    # Primary objective: max margin among eligible; tie-break |ρ|
    def sort_key(r):
        return (
            r["margin_runner_minus_cldn4"] if np.isfinite(r["margin_runner_minus_cldn4"]) else -999,
            r["abs_cldn4_rho"] if np.isfinite(r["abs_cldn4_rho"]) else -999,
        )

    winners = sorted(eligible, key=sort_key, reverse=True)
    # also I2=0 subset
    i2zero = [r for r in winners if r["cldn4_I2"] == 0.0 or r["cldn4_I2"] < 1e-12]
    # locked-like baseline: all4, tj6/full, pct, none
    baseline = [
        r
        for r in grid_rows
        if r["filter"] == "all4" and r["panel"] == "full7" and r["score"] == "pct" and r["adjustment"] == "none"
    ][0]

    primary = winners[0] if winners else None
    primary_i2 = i2zero[0] if i2zero else None

    # Head-to-head perm on primary and on locked-like KRT partial (PR748 default)
    h2h = []
    for label, spec in (
        ("primary_max_margin", primary),
        ("i2zero_max_margin", primary_i2),
        (
            "pr748_krt_pct_full7",
            next(
                r
                for r in grid_rows
                if r["filter"] == "all4"
                and r["panel"] == "full7"
                and r["score"] == "pct"
                and r["adjustment"] == "KRT"
            ),
        ),
    ):
        if spec is None:
            continue
        panel = panels[spec["panel"]]
        sub = {c: by_cohort[c] for c in filters[spec["filter"]]}
        meta_rows = meta_for_panel(sub, panel, spec["score"], spec["adjustment"])
        others = [g for g in panel if g != "CLDN4"]
        for other in others:
            delta, ra, rb, p = label_swap_delta(
                sub, "CLDN4", other, spec["score"], spec["adjustment"], n_perm=2000
            )
            h2h.append(
                {
                    "spec_label": label,
                    "filter": spec["filter"],
                    "panel": spec["panel"],
                    "score": spec["score"],
                    "adjustment": spec["adjustment"],
                    "other": other,
                    "rho_cldn4": ra,
                    "rho_other": rb,
                    "delta_cldn4_minus_other": delta,
                    "perm_p_delta_le": p,
                    "separable_p05": bool(delta < 0 and p < 0.05),
                }
            )
    write_tsv(TAB / "c4_head_to_head.tsv", h2h)

    # Rank table for primary
    rank_out = []
    if primary:
        sub = {c: by_cohort[c] for c in filters[primary["filter"]]}
        meta_rows = meta_for_panel(sub, panels[primary["panel"]], primary["score"], primary["adjustment"])
        ranked = sorted(meta_rows, key=lambda r: r["rho"] if np.isfinite(r["rho"]) else 999)
        for i, r in enumerate(ranked, 1):
            rank_out.append(
                {
                    "spec": "primary_max_margin",
                    "rank": i,
                    "gene": r["gene"],
                    "rho": r["rho"],
                    "p": r["p"],
                    "I2": r["I2"],
                    "ci_lo": r["ci_lo"],
                    "ci_hi": r["ci_hi"],
                    "n_negative_cohorts": r["n_negative_cohorts"],
                    "is_cldn4": r["gene"] == "CLDN4",
                }
            )
    write_tsv(TAB / "c4_primary_rank.tsv", rank_out)

    summary = {
        "n_grid": len(grid_rows),
        "n_eligible_cldn4_lead": len(eligible),
        "n_eligible_i2zero": len(i2zero),
        "baseline_unadj_pct_full7": baseline,
        "primary_max_margin": primary,
        "primary_i2zero_max_margin": primary_i2,
        "locked_cldn4_pct_unadj": {
            "rho": locked["rho"],
            "p": locked["p"],
            "I2": locked["I2"],
            "N": locked["N"],
        },
        "note": "Part2 specificity screen. Searched maximum; do not force into Part1.",
    }
    with (TAB / "c4_sweep_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2, default=lambda o: float(o) if isinstance(o, (np.floating,)) else o)

    # Figure: primary ranks
    if rank_out:
        fig, ax = plt.subplots(figsize=(7.2, 4.0))
        genes = [r["gene"] for r in rank_out][::-1]
        rhos = [r["rho"] for r in rank_out][::-1]
        colors = ["#b45309" if g == "CLDN4" else "#374151" for g in genes]
        ax.barh(range(len(genes)), rhos, color=colors)
        ax.set_yticks(range(len(genes)))
        ax.set_yticklabels(genes)
        ax.axvline(0, color="#9ca3af", lw=1)
        ax.set_xlabel("DL Spearman ρ vs patient T/NK")
        ax.set_title(
            f"C4 max-margin: {primary['score']}/{primary['adjustment']}/"
            f"{primary['filter']}  margin={primary['margin_runner_minus_cldn4']:.3f}"
        )
        fig.tight_layout()
        fig.savefig(FIG / "c4_primary_ranks.png", dpi=150)
        fig.savefig(FIG / "c4_primary_ranks.pdf")
        plt.close(fig)

    # Forest of eligible margins
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    top = winners[:15]
    labels = [
        f"{r['filter']}|{r['panel']}|{r['score']}|{r['adjustment']}" for r in top
    ][::-1]
    margins = [r["margin_runner_minus_cldn4"] for r in top][::-1]
    ax.barh(range(len(labels)), margins, color="#b45309")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("margin = ρ_runner − ρ_CLDN4 (eligible CLDN4 leads)")
    ax.set_title("Concordant-4 TJ screen: top margins")
    fig.tight_layout()
    fig.savefig(FIG / "c4_top_margins.png", dpi=150)
    fig.savefig(FIG / "c4_top_margins.pdf")
    plt.close(fig)

    say("=== PRIMARY ===")
    if primary:
        say(
            f"{primary['filter']} {primary['panel']} {primary['score']} {primary['adjustment']} "
            f"CLDN4 ρ={primary['cldn4_rho']:.3f} runner={primary['runner_gene']} "
            f"ρ={primary['runner_rho']:.3f} margin={primary['margin_runner_minus_cldn4']:.3f}"
        )
    say("done concordant-4 sweep")


if __name__ == "__main__":
    main()
