#!/usr/bin/env python3
"""Multi-sample NicheNet: CLDN4-high vs low malignant → same-patient T/NK.

Cohorts: GSE131907 + GSE205335 (PR #320 winning public merge) + GSE207422.
Patient is the unit. CLDN4-only (no dual-high). Documented NicheNet-v2 prior
+ MultiNicheNet-style ranks. R `nichenetr` / `multinichenetr` are not used.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gene_sets import (  # noqa: E402
    CYTOTOXICITY,
    DETECT_FRAC,
    EXHAUSTION,
    EXTRA_TNK,
    IFN,
    MIN_MALIG_PER_STATE,
    MIN_MALIG_PROGRAM,
    MIN_TNK,
    THIN_Q_N,
    THIN_TAIL,
)
from lib_io import CACHE, DATA, ROOT, log, log1p_cp10k  # noqa: E402
from lib_nichenet import (  # noqa: E402
    ligand_activity,
    ligand_activity_continuous,
    prioritization_table,
    receptors_for,
)
from lib_stats import (  # noqa: E402
    fdr_bh,
    fisher_z_pool,
    rank_biserial_q4q1,
    spearman_safe,
    wilcoxon_paired,
)

OUT = ROOT / "results"
FIG = ROOT / "figures"
COHORTS = ("GSE131907", "GSE205335", "GSE207422")
META_COLS = {
    "barcode",
    "patient",
    "sample",
    "cohort",
    "is_tumor",
    "is_malignant",
    "is_tnk",
    "celltype",
    "note",
    "ncount",
}
WINNING = ("GSE131907", "GSE205335")


def gene_cols(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c not in META_COLS]


def add_logs(df: pd.DataFrame) -> pd.DataFrame:
    genes = gene_cols(df)
    ncount = df["ncount"].to_numpy()
    logs = {f"log_{g}": log1p_cp10k(df[g].to_numpy(), ncount) for g in genes}
    return pd.concat([df, pd.DataFrame(logs, index=df.index)], axis=1)


def mean_score(df: pd.DataFrame, genes: list[str]) -> pd.Series:
    cols = [f"log_{g}" for g in genes if f"log_{g}" in df.columns]
    if not cols:
        return pd.Series(np.nan, index=df.index)
    return df[cols].mean(axis=1)


def load_prior() -> tuple[pd.DataFrame, pd.DataFrame]:
    lr = pd.read_csv(DATA / "lr_network.tsv", sep="\t")
    lt = pd.read_parquet(CACHE / "prior_ligand_target.parquet")
    if lt.shape[0] < lt.shape[1]:
        log(f"warning: prior shape {lt.shape} (expected genes x ligands)")
    return lr, lt


def process_cohort(name: str, lr: pd.DataFrame) -> dict:
    path = CACHE / f"{name}_panel.parquet"
    raw = pd.read_parquet(path)
    n_all = pd.read_csv(CACHE / f"{name}_patient_n.tsv", sep="\t")
    n_map = n_all.set_index("patient")["n_cells"].to_dict()
    raw = add_logs(raw)
    mal = raw[raw["is_malignant"]].copy()
    tnk = raw[raw["is_tnk"]].copy()
    if "log_CLDN4" not in mal.columns:
        raise SystemExit(f"{name}: CLDN4 missing")
    thr = float(mal["log_CLDN4"].median())
    mal["cldn4_high"] = mal["log_CLDN4"] >= thr
    log(
        f"[{name}] cells={len(raw)} mal={len(mal)} tnk={len(tnk)} "
        f"patients={raw['patient'].nunique()} CLDN4_med={thr:.3f} "
        f"high={int(mal['cldn4_high'].sum())} low={int((~mal['cldn4_high']).sum())}"
    )
    ligands = sorted(set(lr["from"].astype(str)) & {c.replace("log_", "") for c in mal.columns if c.startswith("log_")})
    receptors = sorted(set(lr["to"].astype(str)) & {c.replace("log_", "") for c in tnk.columns if c.startswith("log_")})

    rows = []
    lig_rows = []
    tnk_gene_rows = []
    tnk_genes = sorted(
        {g for g in (CYTOTOXICITY + IFN + EXHAUSTION + EXTRA_TNK) if f"log_{g}" in tnk.columns}
    )
    for pid, g in raw.groupby("patient"):
        mg = mal[mal["patient"] == pid]
        tg = tnk[tnk["patient"] == pid]
        hi = mg[mg["cldn4_high"]]
        lo = mg[~mg["cldn4_high"]]
        rec = {
            "cohort": name,
            "patient": pid,
            "note": str(g["note"].iloc[0]) if len(g) else "",
            "n_cells": int(len(g)),
            "n_malignant": int(len(mg)),
            "n_cldn4_high": int(len(hi)),
            "n_cldn4_low": int(len(lo)),
            "n_tnk": int(len(tg)),
            "n_cells_tumor": int(n_map.get(pid, len(g))),
            "tnk_frac": (
                float(len(tg) / n_map[pid]) if n_map.get(pid, 0) else np.nan
            ),
        }
        if len(mg):
            rec["mal_CLDN4_mean"] = float(mg["log_CLDN4"].mean())
            rec["mal_CLDN4_pct"] = float((mg["CLDN4"] > 0).mean()) if "CLDN4" in mg.columns else np.nan
        if len(tg):
            rec["tnk_cyto"] = float(mean_score(tg, CYTOTOXICITY).mean())
            rec["tnk_ifn"] = float(mean_score(tg, IFN).mean())
            rec["tnk_exh"] = float(mean_score(tg, EXHAUSTION).mean())
        rec["paired"] = (
            rec["n_cldn4_high"] >= MIN_MALIG_PER_STATE
            and rec["n_cldn4_low"] >= MIN_MALIG_PER_STATE
            and rec["n_tnk"] >= MIN_TNK
        )
        rec["program"] = rec["n_malignant"] >= MIN_MALIG_PROGRAM and rec["n_tnk"] >= MIN_TNK
        rows.append(rec)
        if rec["paired"]:
            for L in ligands:
                col = f"log_{L}"
                if col not in mg.columns:
                    continue
                lig_rows.append(
                    {
                        "cohort": name,
                        "patient": pid,
                        "ligand": L,
                        "mean_high": float(hi[col].mean()) if len(hi) else np.nan,
                        "mean_low": float(lo[col].mean()) if len(lo) else np.nan,
                        "frac_high": float((hi[L] > 0).mean()) if len(hi) and L in hi.columns else np.nan,
                        "frac_low": float((lo[L] > 0).mean()) if len(lo) and L in lo.columns else np.nan,
                    }
                )
        if rec["program"] and len(tg):
            grow = {"cohort": name, "patient": pid}
            for gene in tnk_genes:
                grow[gene] = float(tg[f"log_{gene}"].mean())
            tnk_gene_rows.append(grow)

    patients = pd.DataFrame(rows).sort_values(["cohort", "patient"])
    lig_pat = pd.DataFrame(lig_rows)
    tnk_pat = pd.DataFrame(tnk_gene_rows)
    return {
        "name": name,
        "raw": raw,
        "mal": mal,
        "tnk": tnk,
        "thr": thr,
        "ligands": ligands,
        "receptors": receptors,
        "patients": patients,
        "lig_pat": lig_pat,
        "tnk_pat": tnk_pat,
        "tnk_genes": tnk_genes,
    }


def ligand_de_table(lig_pat: pd.DataFrame) -> pd.DataFrame:
    if lig_pat.empty:
        return pd.DataFrame()
    rows = []
    for (cohort, lig), g in lig_pat.groupby(["cohort", "ligand"]):
        delta = g["mean_high"] - g["mean_low"]
        w = wilcoxon_paired(g["mean_high"].to_numpy(), g["mean_low"].to_numpy())
        rows.append(
            {
                "cohort": cohort,
                "ligand": lig,
                "n_patients": int(len(g)),
                "median_delta": float(delta.median()),
                "mean_delta": float(delta.mean()),
                "paired_p": w["p"],
                "frac_patients_expressed_high": float((g["frac_high"] >= DETECT_FRAC).mean()),
                "median_frac_high": float(g["frac_high"].median()),
                "median_frac_low": float(g["frac_low"].median()),
            }
        )
    out = pd.DataFrame(rows)
    out["paired_fdr"] = out.groupby("cohort")["paired_p"].transform(fdr_bh)
    return out


def pooled_ligand_de(lig_pat: pd.DataFrame) -> pd.DataFrame:
    if lig_pat.empty:
        return pd.DataFrame()
    rows = []
    for lig, g in lig_pat.groupby("ligand"):
        delta = g["mean_high"] - g["mean_low"]
        w = wilcoxon_paired(g["mean_high"].to_numpy(), g["mean_low"].to_numpy())
        rows.append(
            {
                "ligand": lig,
                "n_patients": int(g["patient"].nunique()),
                "n_cohorts": int(g["cohort"].nunique()),
                "median_delta": float(delta.median()),
                "mean_delta": float(delta.mean()),
                "paired_p": w["p"],
                "frac_patients_expressed_high": float((g["frac_high"] >= DETECT_FRAC).mean()),
            }
        )
    out = pd.DataFrame(rows)
    out["paired_fdr"] = fdr_bh(out["paired_p"])
    return out


def potential_ligands(de: pd.DataFrame, lr: pd.DataFrame, receptors: set[str], min_frac: float = 0.25) -> list[str]:
    keep = de[de["frac_patients_expressed_high"] >= min_frac]["ligand"]
    out = []
    for L in keep:
        recs = set(receptors_for(lr, L)) & receptors
        if recs:
            out.append(L)
    return sorted(set(out))


def empirical_tnk_de(tnk_pat: pd.DataFrame, patients: pd.DataFrame) -> pd.Series:
    """Patient-level T/NK gene contrast: CLDN4-high vs low patients (median split)."""
    prog = patients[patients["program"]].copy()
    if prog.empty or tnk_pat.empty:
        return pd.Series(dtype=float)
    med = prog["mal_CLDN4_pct"].median()
    prog["hi"] = prog["mal_CLDN4_pct"] >= med
    merged = tnk_pat.merge(prog[["cohort", "patient", "hi"]], on=["cohort", "patient"])
    genes = [c for c in tnk_pat.columns if c not in {"cohort", "patient"}]
    scores = {}
    for g in genes:
        a = merged.loc[~merged["hi"], g].dropna()
        b = merged.loc[merged["hi"], g].dropna()
        if len(a) < 3 or len(b) < 3:
            continue
        # signed rank-biserial of high minus low (negative = lower in CLDN4-high patients)
        from lib_stats import mannwhitney

        mw = mannwhitney(a.to_numpy(), b.to_numpy())
        scores[g] = mw["rank_biserial"]
    return pd.Series(scores, dtype=float)


def program_tests(patients: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cohort, g in patients.groupby("cohort"):
        prog = g[g["program"]]
        for metric in ["tnk_frac", "tnk_cyto", "tnk_ifn", "tnk_exh"]:
            sp = spearman_safe(prog["mal_CLDN4_pct"], prog[metric])
            q = rank_biserial_q4q1(prog[metric], prog["mal_CLDN4_pct"])
            thin = (len(prog) < THIN_Q_N) or (q["n_Q1"] < THIN_TAIL) or (q["n_Q4"] < THIN_TAIL)
            rows.append(
                {
                    "cohort": cohort,
                    "metric": metric,
                    "n_program": int(len(prog)),
                    "spearman_rho": sp["rho"],
                    "spearman_p": sp["p"],
                    "q4q1_r": q["rank_biserial"],
                    "q4q1_p": q["p"],
                    "n_Q1": q["n_Q1"],
                    "n_Q4": q["n_Q4"],
                    "median_Q1": q["median_Q1"],
                    "median_Q4": q["median_Q4"],
                    "thin_q4q1": thin,
                }
            )
    return pd.DataFrame(rows)


def pool_program(tests: pd.DataFrame, metric: str, cohorts: tuple[str, ...]) -> dict:
    sub = tests[(tests["metric"] == metric) & (tests["cohort"].isin(cohorts))]
    rows = [{"n": int(r.n_program), "rho": r.spearman_rho} for r in sub.itertuples()]
    out = fisher_z_pool(rows)
    out["metric"] = metric
    out["cohorts"] = "+".join(cohorts)
    return out


def receptor_means(cohorts: list[dict], receptors: list[str]) -> pd.Series:
    acc = {r: [] for r in receptors}
    for C in cohorts:
        tnk = C["tnk"]
        for r in receptors:
            col = f"log_{r}"
            if col in tnk.columns:
                acc[r].append(float(tnk[col].mean()))
    return pd.Series({k: float(np.mean(v)) if v else np.nan for k, v in acc.items()})


def fmt(x, nd=3):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    if isinstance(x, float):
        if abs(x) >= 0.001 or x == 0:
            return f"{x:.{nd}g}" if nd == 3 else f"{x:.{nd}f}"
        return f"{x:.2e}"
    return str(x)


def pfmt(p) -> str:
    if p is None or (isinstance(p, float) and not np.isfinite(p)):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3g}"


def write_figures(patients: pd.DataFrame, tests: pd.DataFrame, prio: pd.DataFrame, de_all: pd.DataFrame, pools: list[dict]) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    colors = {"GSE131907": "#1f77b4", "GSE205335": "#ff7f0e", "GSE207422": "#2ca02c"}

    # 1. inclusion n
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    labels = []
    vals = []
    for c in COHORTS:
        g = patients[patients["cohort"] == c]
        labels += [f"{c}\nprogram", f"{c}\npaired"]
        vals += [int(g["program"].sum()), int(g["paired"].sum())]
    ax.bar(range(len(vals)), vals, color=["#444"] * 2 + ["#888"] * 2 + ["#bbb"] * 2)
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("patients")
    ax.set_title("Honest n: program (≥20 mal + ≥20 T/NK) vs paired high/low")
    fig.tight_layout()
    fig.savefig(FIG / "n_patients_inclusion.png", dpi=160)
    fig.savefig(FIG / "n_patients_inclusion.pdf")
    plt.close(fig)

    # 2. CLDN4 %pos vs T/NK fraction
    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    for c, g in patients[patients["program"]].groupby("cohort"):
        ax.scatter(g["mal_CLDN4_pct"], g["tnk_frac"], s=28, alpha=0.85, label=c, color=colors[c])
    ax.set_xlabel("malignant CLDN4 %pos")
    ax.set_ylabel("same-patient T/NK fraction")
    ax.set_title("CLDN4-only vs T/NK fraction (patient unit)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "cldn4_vs_tnk_frac.png", dpi=160)
    fig.savefig(FIG / "cldn4_vs_tnk_frac.pdf")
    plt.close(fig)

    # 3. programs
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.4), sharex=True)
    for ax, metric, title in zip(
        axes, ["tnk_cyto", "tnk_ifn", "tnk_exh"], ["cytotoxicity", "IFN", "exhaustion"]
    ):
        for c, g in patients[patients["program"]].groupby("cohort"):
            ax.scatter(g["mal_CLDN4_pct"], g[metric], s=22, alpha=0.85, label=c, color=colors[c])
        ax.set_title(title)
        ax.set_xlabel("mal CLDN4 %pos")
    axes[0].set_ylabel("T/NK mean log1p(CP10k)")
    axes[2].legend(fontsize=7)
    fig.suptitle("Same-patient T/NK programs vs malignant CLDN4", y=1.02)
    fig.tight_layout()
    fig.savefig(FIG / "cldn4_vs_tnk_programs.png", dpi=160)
    fig.savefig(FIG / "cldn4_vs_tnk_programs.pdf")
    plt.close(fig)

    # 4. Spearman forest T/NK frac
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    sub = tests[tests["metric"] == "tnk_frac"].copy()
    ylabels, rhos, lo, hi, ys = [], [], [], [], []
    y = 0
    for r in sub.itertuples():
        ylabels.append(f"{r.cohort} n={r.n_program}")
        rhos.append(r.spearman_rho)
        # crude 95% via Fisher z
        n = max(int(r.n_program), 4)
        z = np.arctanh(np.clip(r.spearman_rho, -0.999, 0.999))
        se = 1 / np.sqrt(n - 3)
        lo.append(np.tanh(z - 1.96 * se))
        hi.append(np.tanh(z + 1.96 * se))
        ys.append(y)
        y += 1
    for pool in pools:
        if pool["metric"] != "tnk_frac":
            continue
        ylabels.append(f"{pool['cohorts']} N={pool['N']}")
        rhos.append(pool["rho"])
        z = np.arctanh(np.clip(pool["rho"], -0.999, 0.999))
        se = 1 / np.sqrt(max(pool["N"] - 3 * pool["k"], 4))
        lo.append(np.tanh(z - 1.96 * se))
        hi.append(np.tanh(z + 1.96 * se))
        ys.append(y)
        y += 1
    ax.axvline(0, color="0.6", lw=0.8)
    ax.errorbar(rhos, ys, xerr=[np.array(rhos) - np.array(lo), np.array(hi) - np.array(rhos)], fmt="o", color="#222")
    ax.set_yticks(ys)
    ax.set_yticklabels(ylabels, fontsize=8)
    ax.set_xlabel("Spearman ρ (CLDN4 %pos vs T/NK fraction)")
    ax.set_title("Patient-level Spearman (descriptive p)")
    fig.tight_layout()
    fig.savefig(FIG / "forest_spearman_tnk_frac.png", dpi=160)
    fig.savefig(FIG / "forest_spearman_tnk_frac.pdf")
    plt.close(fig)

    # 5. Q4 vs Q1 boxes for winning+207422 T/NK frac
    fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.4), sharey=True)
    for ax, c in zip(axes, COHORTS):
        g = patients[(patients["cohort"] == c) & (patients["program"])]
        from lib_stats import quartile_tails

        low, high = quartile_tails(g.set_index("patient")["mal_CLDN4_pct"])
        q1 = g[g["patient"].isin(low)]["tnk_frac"]
        q4 = g[g["patient"].isin(high)]["tnk_frac"]
        ax.boxplot([q1.to_numpy(), q4.to_numpy()], tick_labels=["Q1", "Q4"], widths=0.55)
        ax.scatter([1] * len(q1), q1, s=16, color=colors[c], alpha=0.8)
        ax.scatter([2] * len(q4), q4, s=16, color=colors[c], alpha=0.8)
        ax.set_title(f"{c}\n{len(q1)} vs {len(q4)}")
    axes[0].set_ylabel("T/NK fraction")
    fig.suptitle("Q4 vs Q1 malignant CLDN4 %pos (tails only)", y=1.02)
    fig.tight_layout()
    fig.savefig(FIG / "q4q1_tnk_frac_box.png", dpi=160)
    fig.savefig(FIG / "q4q1_tnk_frac_box.pdf")
    plt.close(fig)

    # 6. ligand activity / prioritization
    if not prio.empty:
        top = prio.head(20)
        fig, ax = plt.subplots(figsize=(7.2, 5.2))
        ax.barh(top["ligand"][::-1], top["prioritization_score"][::-1], color="#4c72b0")
        ax.set_xlabel("MultiNicheNet-style prioritization (equal-weight scaled)")
        ax.set_title("Top ligands: CLDN4-high vs low malignant → T/NK")
        fig.tight_layout()
        fig.savefig(FIG / "ligand_prioritization_top.png", dpi=160)
        fig.savefig(FIG / "ligand_prioritization_top.pdf")
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(6.4, 5.0))
        ax.barh(top["ligand"][::-1], top["pearson"][::-1], color="#dd8452")
        ax.set_xlabel("NicheNet-v2 ligand activity (Pearson vs empirical T/NK DE)")
        ax.set_title("Ligand activity on receiver DE (patient unit)")
        fig.tight_layout()
        fig.savefig(FIG / "ligand_activity_empirical.png", dpi=160)
        fig.savefig(FIG / "ligand_activity_empirical.pdf")
        plt.close(fig)

    # 7. paired delta forest top
    if not de_all.empty:
        topd = de_all.reindex(de_all["median_delta"].abs().sort_values(ascending=False).index).head(15)
        fig, ax = plt.subplots(figsize=(6.8, 4.8))
        y = np.arange(len(topd))
        ax.axvline(0, color="0.6", lw=0.8)
        ax.scatter(topd["median_delta"], y, c=["#c44" if v < 0 else "#4c4" for v in topd["median_delta"]])
        ax.set_yticks(y)
        ax.set_yticklabels([f"{r.ligand} n={int(r.n_patients)}" for r in topd.itertuples()], fontsize=8)
        ax.set_xlabel("median paired Δ (CLDN4-high − low malignant)")
        ax.set_title("Sender ligand DE (patient-paired)")
        fig.tight_layout()
        fig.savefig(FIG / "ligand_delta_top.png", dpi=160)
        fig.savefig(FIG / "ligand_delta_top.pdf")
        plt.close(fig)

    # 8. extra: per-cohort paired n vs dropped
    fig, ax = plt.subplots(figsize=(5.6, 3.4))
    x = np.arange(len(COHORTS))
    prog_n = [int(patients[patients["cohort"] == c]["program"].sum()) for c in COHORTS]
    pair_n = [int(patients[patients["cohort"] == c]["paired"].sum()) for c in COHORTS]
    all_n = [int((patients["cohort"] == c).sum()) for c in COHORTS]
    ax.bar(x - 0.25, all_n, 0.25, label="any tumor extract", color="#bbb")
    ax.bar(x, prog_n, 0.25, label="program n", color="#888")
    ax.bar(x + 0.25, pair_n, 0.25, label="paired high/low", color="#444")
    ax.set_xticks(x)
    ax.set_xticklabels(COHORTS)
    ax.set_ylabel("patients")
    ax.legend(fontsize=8)
    ax.set_title("Extra: who is dropped from each n")
    fig.tight_layout()
    fig.savefig(FIG / "extra_n_dropped.png", dpi=160)
    fig.savefig(FIG / "extra_n_dropped.pdf")
    plt.close(fig)


def write_finding(
    patients: pd.DataFrame,
    tests: pd.DataFrame,
    pools: list[dict],
    prio: pd.DataFrame,
    de_all: pd.DataFrame,
    ntab: pd.DataFrame,
    method_note: str,
) -> None:
    prog = patients[patients["program"]]
    paired = patients[patients["paired"]]
    win = [p for p in pools if p["cohorts"] == "GSE131907+GSE205335"]
    three = [p for p in pools if p["cohorts"] == "GSE131907+GSE205335+GSE207422"]

    def pool_row(items, metric):
        hits = [p for p in items if p["metric"] == metric]
        return hits[0] if hits else {"k": 0, "N": 0, "rho": np.nan, "p": np.nan, "I2": np.nan}

    w_frac = pool_row(win, "tnk_frac")
    t_frac = pool_row(three, "tnk_frac")
    w_cyto = pool_row(win, "tnk_cyto")
    t_cyto = pool_row(three, "tnk_cyto")
    w_ifn = pool_row(win, "tnk_ifn")
    t_ifn = pool_row(three, "tnk_ifn")

    top = prio.head(15) if not prio.empty else pd.DataFrame()
    lines = [
        "# FINDING — Multi-sample NicheNet, CLDN4-high vs low malignant → same-patient T/NK",
        "",
        "ADDITIVE **CLDN4-only**. No dual-high TACSTD2×CLDN4. Patient is the unit.",
        "Sender = author-malignant cells split by CLDN4 `log1p(CP10k)` at the **cohort-wide",
        "malignant median**. Receiver = **same-patient** T/NK. Run on the PR #320 winning",
        "public merge **GSE131907+GSE205335** plus **GSE207422** as a third cohort.",
        "",
        f"**Method:** {method_note}",
        "",
        "## Verdict",
        "",
    ]
    # honest verdict from numbers
    rho3 = t_frac.get("rho", np.nan)
    p3 = t_frac.get("p", np.nan)
    if np.isfinite(rho3) and rho3 < 0 and np.isfinite(p3) and p3 < 0.05:
        verdict = (
            f"Malignant CLDN4 %pos is negatively associated with same-patient T/NK fraction "
            f"in the three-cohort patient pool (ρ={rho3:.3f}, p={pfmt(p3)}, N={t_frac['N']}). "
            "Ligand ranks below are prior + paired sender DE, not proof of causation."
        )
    elif np.isfinite(rho3) and rho3 < 0:
        verdict = (
            f"Direction is CLDN4-high → fewer T/NK (three-cohort ρ={rho3:.3f}, p={pfmt(p3)}, "
            f"N={t_frac['N']}) but the signed program tests are not uniformly significant. "
            "Ligand activity is a documented prior score, not a causal induction test."
        )
    else:
        verdict = (
            f"The three-cohort patient pool does **not** support a clean CLDN4-high → T/NK-low "
            f"program (ρ={fmt(rho3)}, p={pfmt(p3)}, N={t_frac.get('N', 0)}). "
            "Ligand ranks are still reported as the reusable pipeline output."
        )
    lines += [verdict, ""]
    lines += [
        "## Honest n",
        "",
        "| Item | n | Note |",
        "| --- | ---: | --- |",
        f"| Program patients (all 3) | **{int(prog['patient'].nunique())}** | ≥20 malignant + ≥20 T/NK; unit of Spearman / Q4Q1 |",
        f"| Paired patients (all 3) | **{int(paired['patient'].nunique())}** | ≥10 CLDN4-high + ≥10 low + ≥20 T/NK; unit of ligand Δ |",
        f"| GSE131907 program / paired | {int(prog[prog.cohort=='GSE131907'].shape[0])} / {int(paired[paired.cohort=='GSE131907'].shape[0])} | tumor-origin samples pooled per `patient_id` |",
        f"| GSE205335 program / paired | {int(prog[prog.cohort=='GSE205335'].shape[0])} / {int(paired[paired.cohort=='GSE205335'].shape[0])} | tumor samples pooled; normals dropped |",
        f"| GSE207422 program / paired | {int(prog[prog.cohort=='GSE207422'].shape[0])} / {int(paired[paired.cohort=='GSE207422'].shape[0])} | post-tx DRMref; pCR counted as MPR in the note column |",
        f"| Winning merge program N | **{int(prog[prog.cohort.isin(WINNING)].shape[0])}** | GSE131907+GSE205335 (PR #320 given) |",
        f"| Three-cohort program N | **{int(prog.shape[0])}** | +GSE207422 |",
        "",
        "Cells are counts, not n. GSE131907 is patient-pooled (not the sample-level T/NK extract in PR #320).",
        "Q4 vs Q1 uses quartile **tails only**; thin tails (n<8 or a tail <3) are flagged and not treated as a pool.",
        "",
        "## Patient-level T/NK vs malignant CLDN4 %pos",
        "",
        "| cohort | n | ρ T/NK frac (p) | Q4 vs Q1 r (p; n_Q1/n_Q4) | ρ cyto (p) | ρ IFN (p) |",
        "| --- | ---: | --- | --- | --- | --- |",
    ]
    for c in COHORTS:
        tf = tests[(tests.cohort == c) & (tests.metric == "tnk_frac")].iloc[0]
        cy = tests[(tests.cohort == c) & (tests.metric == "tnk_cyto")].iloc[0]
        inf = tests[(tests.cohort == c) & (tests.metric == "tnk_ifn")].iloc[0]
        thin = " thin" if tf.thin_q4q1 else ""
        lines.append(
            f"| {c} | {int(tf.n_program)} | {tf.spearman_rho:.3f} ({pfmt(tf.spearman_p)}) | "
            f"{tf.q4q1_r:.3f} ({pfmt(tf.q4q1_p)}; {int(tf.n_Q1)}/{int(tf.n_Q4)}){thin} | "
            f"{cy.spearman_rho:.3f} ({pfmt(cy.spearman_p)}) | {inf.spearman_rho:.3f} ({pfmt(inf.spearman_p)}) |"
        )
    lines += [
        f"| **GSE131907+GSE205335** | {w_frac['N']} | **{fmt(w_frac['rho'])} ({pfmt(w_frac['p'])}, I²={fmt(w_frac['I2'])}%)** | — | "
        f"{fmt(w_cyto['rho'])} ({pfmt(w_cyto['p'])}) | {fmt(w_ifn['rho'])} ({pfmt(w_ifn['p'])}) |",
        f"| **+GSE207422 (3 cohorts)** | {t_frac['N']} | **{fmt(t_frac['rho'])} ({pfmt(t_frac['p'])}, I²={fmt(t_frac['I2'])}%)** | — | "
        f"{fmt(t_cyto['rho'])} ({pfmt(t_cyto['p'])}) | {fmt(t_ifn['rho'])} ({pfmt(t_ifn['p'])}) |",
        "",
        "p-values are descriptive. Fisher-z pool uses patient n per cohort (k=2 or 3).",
        "",
        "## Ligand-activity table (primary)",
        "",
        "Full table: `results/ligand_activity.tsv`. Rank = MultiNicheNet-style equal-weight",
        "mean of min-max scaled (1) paired sender Δ, (2) NicheNet-v2 activity vs empirical",
        "T/NK DE (CLDN4-high vs low **patients**), (3) receptor expression in T/NK,",
        "(4) fraction of paired patients with ligand detected in ≥10% of CLDN4-high senders.",
        "",
    ]
    if top.empty:
        lines.append("_(no potential ligands passed the expression/receptor gate)_")
    else:
        lines += [
            "| rank | ligand | n patients | n cohorts | median Δ | paired p | activity r | prio |",
            "| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
        for r in top.itertuples():
            lines.append(
                f"| {int(r.rank)} | {r.ligand} | {int(r.n_patients)} | {int(r.n_cohorts)} | "
                f"{r.median_delta:.3f} | {pfmt(r.paired_p)} | {r.pearson:.3f} | {r.prioritization_score:.3f} |"
            )
    lines += [
        "",
        "Activity Pearson is prior-structure (unsigned). A high IFN/MHC ligand score does",
        "**not** mean CLDN4-high cells induce that program — test direction on the patient rows above.",
        "",
        "## What this is not",
        "",
        "- Not R `nichenetr` or `multinichenetr` (packages absent; Python prior scoring is used).",
        "- Not dual-high TACSTD2×CLDN4. TACSTD2 is never a gate.",
        "- Not CellChat / LIANA. Those live in other folders.",
        "- Not a full-transcriptome NicheNet run. Background = panel ∩ prior ∩ T/NK-expressed genes.",
        "- GSE131907 has no ICI/MPR labels. GSE205335 RECIST is not substituted for MPR.",
        "- Cells are not replicates.",
        "",
        "## Extra figures",
        "",
        "- `figures/n_patients_inclusion.png`",
        "- `figures/cldn4_vs_tnk_frac.png`",
        "- `figures/cldn4_vs_tnk_programs.png`",
        "- `figures/forest_spearman_tnk_frac.png`",
        "- `figures/q4q1_tnk_frac_box.png`",
        "- `figures/ligand_prioritization_top.png`",
        "- `figures/ligand_activity_empirical.png`",
        "- `figures/ligand_delta_top.png`",
        "- `figures/extra_n_dropped.png`",
        "",
        "## Reproduce",
        "",
        "```bash",
        "cd methods/multinichenet_cldn4",
        "python3 scripts/00_download.py",
        "python3 scripts/01_convert_prior.py",
        "python3 scripts/02_extract_panels.py",
        "python3 scripts/03_analyze.py",
        "```",
        "",
    ]
    text = "\n".join(lines) + "\n"
    (ROOT / "FINDING.md").write_text(text)
    (OUT / "FINDING.md").write_text(text)
    log(f"wrote {ROOT / 'FINDING.md'}")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    lr, lt = load_prior()
    method_note = (
        "Documented NicheNet-v2 ligand–target prior (`ligand_target_matrix_nsga2r_final`, "
        "Zenodo 10.5281/zenodo.7074291) scored in Python (`rdata`); MultiNicheNet-style "
        "equal-weight prioritization (Bonte/Browaeys vignette). "
        "R packages `nichenetr` and `multinichenetr` were **not** installed."
    )
    log(method_note)

    packed = [process_cohort(c, lr) for c in COHORTS]
    patients = pd.concat([C["patients"] for C in packed], ignore_index=True)
    lig_pat = pd.concat([C["lig_pat"] for C in packed], ignore_index=True)
    tnk_pat = pd.concat([C["tnk_pat"] for C in packed], ignore_index=True)
    patients.to_csv(OUT / "patient_metrics.tsv", sep="\t", index=False)
    lig_pat.to_csv(OUT / "ligand_patient_means.tsv", sep="\t", index=False)

    tests = program_tests(patients)
    tests.to_csv(OUT / "program_tests.tsv", sep="\t", index=False)
    pools = []
    for metric in ["tnk_frac", "tnk_cyto", "tnk_ifn", "tnk_exh"]:
        pools.append(pool_program(tests, metric, WINNING))
        pools.append(pool_program(tests, metric, COHORTS))
    pd.DataFrame(pools).to_csv(OUT / "program_pools.tsv", sep="\t", index=False)

    de_coh = ligand_de_table(lig_pat)
    de_all = pooled_ligand_de(lig_pat)
    de_coh.to_csv(OUT / "ligand_de_by_cohort.tsv", sep="\t", index=False)
    de_all.to_csv(OUT / "ligand_de_pooled.tsv", sep="\t", index=False)

    rec_set = set()
    lig_set = set()
    for C in packed:
        rec_set.update(C["receptors"])
        lig_set.update(C["ligands"])
    pots = potential_ligands(de_all, lr, rec_set, min_frac=0.25)
    log(f"potential ligands (pooled, ≥25% patients expressed in high): {len(pots)}")

    # background: T/NK-expressed panel genes present in prior
    tnk_all = pd.concat([C["tnk"] for C in packed], ignore_index=True)
    background = []
    for g in gene_cols(tnk_all):
        if g in lt.index and float((tnk_all[g] > 0).mean()) >= DETECT_FRAC:
            background.append(g)
    log(f"T/NK background genes in prior: {len(background)}")

    act_ifn = ligand_activity(lt, set(IFN), background, pots)
    act_cyto = ligand_activity(lt, set(CYTOTOXICITY), background, pots)
    act_exh = ligand_activity(lt, set(EXHAUSTION), background, pots)
    emp = empirical_tnk_de(tnk_pat, patients)
    emp.to_csv(OUT / "tnk_empirical_de_rankbiserial.tsv", sep="\t")
    act_emp = ligand_activity_continuous(lt, emp, pots)
    act_ifn.to_csv(OUT / "ligand_activity_ifn.tsv", sep="\t", index=False)
    act_cyto.to_csv(OUT / "ligand_activity_cyto.tsv", sep="\t", index=False)
    act_exh.to_csv(OUT / "ligand_activity_exh.tsv", sep="\t", index=False)
    act_emp.to_csv(OUT / "ligand_activity_empirical.tsv", sep="\t", index=False)

    rec_score = receptor_means(packed, sorted(rec_set))
    # receptor score per ligand = max receptor mean among its T/NK receptors
    rec_for_lig = {}
    rec_names = {}
    for L in pots:
        recs = [r for r in receptors_for(lr, L) if r in rec_score.index]
        rec_for_lig[L] = float(np.nanmax(rec_score.reindex(recs))) if recs else np.nan
        rec_names[L] = ",".join(recs[:6])
    rec_series = pd.Series(rec_for_lig)

    if act_emp.empty:
        activity = act_ifn.rename(columns={"pearson": "pearson", "pearson_p": "pearson_p"})
        log("empirical activity empty; falling back to IFN gene-set activity")
    else:
        activity = act_emp
    prio = prioritization_table(de_all[de_all["ligand"].isin(pots)].copy(), activity, rec_series)
    if not prio.empty:
        prio["receptors"] = prio["ligand"].map(rec_names)
        if not act_ifn.empty:
            prio = prio.merge(
                act_ifn[["ligand", "pearson", "auroc"]].rename(
                    columns={"pearson": "activity_ifn", "auroc": "auroc_ifn"}
                ),
                on="ligand",
                how="left",
            )
        if not act_cyto.empty:
            prio = prio.merge(
                act_cyto[["ligand", "pearson"]].rename(columns={"pearson": "activity_cyto"}),
                on="ligand",
                how="left",
            )
        if not act_exh.empty:
            prio = prio.merge(
                act_exh[["ligand", "pearson"]].rename(columns={"pearson": "activity_exh"}),
                on="ligand",
                how="left",
            )
    prio.to_csv(OUT / "ligand_activity.tsv", sep="\t", index=False)
    if not prio.empty:
        prio.head(40).to_csv(OUT / "top_ligands.tsv", sep="\t", index=False)

    ntab = patients[
        [
            "cohort",
            "patient",
            "n_malignant",
            "n_cldn4_high",
            "n_cldn4_low",
            "n_tnk",
            "tnk_frac",
            "mal_CLDN4_mean",
            "mal_CLDN4_pct",
            "paired",
            "program",
            "note",
        ]
    ]
    ntab.to_csv(OUT / "n_table.tsv", sep="\t", index=False)

    write_figures(patients, tests, prio, de_all, pools)
    write_finding(patients, tests, pools, prio, de_all, ntab, method_note)

    summary = {
        "method": method_note,
        "n_program": int(patients["program"].sum()),
        "n_paired": int(patients["paired"].sum()),
        "n_potential_ligands": int(len(pots)),
        "n_background": int(len(background)),
        "pools": pools,
        "cldn4_only": True,
        "dual_high": False,
        "r_nichenetr": False,
        "r_multinichenetr": False,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    log("ok")
    return 0


if __name__ == "__main__":
    sys.exit(main())
