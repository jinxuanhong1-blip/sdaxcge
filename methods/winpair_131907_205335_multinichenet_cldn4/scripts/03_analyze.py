#!/usr/bin/env python3
"""Multi-sample NicheNet (MultiNicheNet-style) on the winning pair.

Sender = CLDN4-high vs low malignant (CLDN4 only; no dual-high).
Receiver = same-patient T/NK. Patient is the unit.
Cohorts = GSE131907 + GSE205335 only. GSE207422 is not used (PR #334 NS).

nichenetr / multinichenetr are not installed. Ligand activity is the
documented NicheNet-v2 prior score (Pearson / AUROC / AUPR of prior
target weights vs gene-set membership; Browaeys et al. Nat Methods 2020
and MultiNicheNet Nat Commun 2023 prioritization logic).
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
from scipy import stats
from sklearn.metrics import average_precision_score, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
CACHE = Path("/tmp/winpair_multinichenet_cldn4")
OUT = ROOT / "results"
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "scripts"))
from gene_sets import CYTOTOXICITY, EXHAUSTION, IFN  # noqa: E402

DETECT_FRAC = 0.10
TOP_N = 15
EMPIRICAL_P = 0.15
PRIMARY_SETS = ("a_priori_ifn", "a_priori_cytotoxicity")
DATASETS = ("GSE131907", "GSE205335")


def log(msg: str) -> None:
    print(msg, flush=True)


def fmt_p(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    return f"{value:.2e}" if value < 0.001 else f"{value:.3g}"


def fmt_num(value: float, digits: int = 3) -> str:
    if not np.isfinite(value):
        return "NA"
    return f"{value:+.{digits}f}" if value != 0 else f"{value:.{digits}f}"


def spearman_safe(x: pd.Series, y: pd.Series) -> dict:
    m = x.notna() & y.notna()
    n = int(m.sum())
    if n < 4:
        return {"n": n, "rho": np.nan, "p": np.nan}
    rho, p = stats.spearmanr(x[m], y[m])
    return {"n": n, "rho": float(rho), "p": float(p)}


def wilcoxon_paired(a: pd.Series, b: pd.Series) -> dict:
    x = a.to_numpy(dtype=float)
    y = b.to_numpy(dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    n = int(len(x))
    if n < 3 or np.allclose(x, y):
        return {
            "n": n,
            "mean_high": float(np.mean(x)) if n else np.nan,
            "mean_low": float(np.mean(y)) if n else np.nan,
            "median_delta": float(np.median(x - y)) if n else np.nan,
            "mean_delta": float(np.mean(x - y)) if n else np.nan,
            "W": np.nan,
            "p": np.nan,
        }
    try:
        res = stats.wilcoxon(x, y, zero_method="wilcox", alternative="two-sided")
        w, p = float(res.statistic), float(res.pvalue)
    except ValueError:
        w, p = np.nan, np.nan
    return {
        "n": n,
        "mean_high": float(np.mean(x)),
        "mean_low": float(np.mean(y)),
        "median_delta": float(np.median(x - y)),
        "mean_delta": float(np.mean(x - y)),
        "W": w,
        "p": p,
    }


def mwu(a: pd.Series, b: pd.Series) -> dict:
    xa = a.dropna().to_numpy(dtype=float)
    xb = b.dropna().to_numpy(dtype=float)
    na, nb = len(xa), len(xb)
    if na < 2 or nb < 2:
        return {"n_a": na, "n_b": nb, "mean_a": np.nan, "mean_b": np.nan, "U": np.nan, "p": np.nan}
    res = stats.mannwhitneyu(xb, xa, alternative="two-sided", method="auto")
    return {
        "n_a": na,
        "n_b": nb,
        "mean_a": float(np.mean(xa)),
        "mean_b": float(np.mean(xb)),
        "median_a": float(np.median(xa)),
        "median_b": float(np.median(xb)),
        "U": float(res.statistic),
        "p": float(res.pvalue),
    }


def bh_fdr(pvals: pd.Series) -> pd.Series:
    p = pvals.to_numpy(dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return pd.Series(out, index=pvals.index)
    pv = p[ok]
    n = len(pv)
    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    tmp = np.empty(n, dtype=float)
    tmp[order] = q
    out[ok] = tmp
    return pd.Series(out, index=pvals.index)


def minmax(s: pd.Series) -> pd.Series:
    v = s.astype(float)
    lo, hi = v.min(), v.max()
    if not np.isfinite(lo) or not np.isfinite(hi) or np.isclose(hi, lo):
        return pd.Series(0.5, index=s.index)
    return (v - lo) / (hi - lo)


def ligand_activity(
    lt: pd.DataFrame, geneset: set[str], background: list[str], ligands: list[str]
) -> pd.DataFrame:
    genes = [g for g in background if g in lt.index]
    if len(genes) < 20 or not ligands:
        return pd.DataFrame()
    y = np.array([1 if g in geneset else 0 for g in genes], dtype=int)
    if y.sum() < 2 or y.sum() > len(y) - 2:
        return pd.DataFrame()
    use = [L for L in ligands if L in lt.columns]
    X = lt.loc[genes, use].to_numpy(dtype=float)
    rows = []
    for j, L in enumerate(use):
        s = X[:, j]
        if not np.isfinite(s).all() or np.allclose(s, s[0]):
            continue
        pear = stats.pearsonr(s, y)
        try:
            auroc = float(roc_auc_score(y, s))
        except ValueError:
            auroc = np.nan
        try:
            aupr = float(average_precision_score(y, s))
        except ValueError:
            aupr = np.nan
        rows.append(
            {
                "ligand": L,
                "pearson": float(pear.statistic),
                "pearson_p": float(pear.pvalue),
                "auroc": auroc,
                "aupr": aupr,
                "n_background": len(genes),
                "n_geneset": int(y.sum()),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values("pearson", ascending=False).reset_index(drop=True)


def program_score(stats: pd.DataFrame, genes: list[str], compartment: str) -> pd.Series:
    sub = stats[(stats["compartment"] == compartment) & (stats["gene"].isin(genes))]
    if sub.empty:
        return pd.Series(dtype=float)
    return sub.groupby(["dataset", "patient"])["mean_log1p"].mean()


def md_table(df: pd.DataFrame, cols: list[str], fmt: dict[str, str] | None = None) -> list[str]:
    if df is None or df.empty:
        return ["_(empty)_", ""]
    fmt = fmt or {}
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for rec in df.itertuples(index=False):
        row = []
        for c in cols:
            v = getattr(rec, c)
            if c in fmt and pd.notna(v):
                row.append(fmt[c].format(v))
            elif isinstance(v, float) and pd.notna(v):
                row.append(f"{v:.3g}")
            else:
                row.append(str(v))
        lines.append("| " + " | ".join(row) + " |")
    lines.append("")
    return lines


def write_finding(summary: dict, ntab: pd.DataFrame, tests: pd.DataFrame, table: pd.DataFrame) -> None:
    n_elig = int(summary["n"]["eligible_patients_combined"])
    n_131 = int(summary["n"]["eligible_GSE131907"])
    n_205 = int(summary["n"]["eligible_GSE205335"])
    rho_ifn = summary["patient_tests"].get("spearman_mal_CLDN4__tnk_ifn", {})
    rho_cy = summary["patient_tests"].get("spearman_mal_CLDN4__tnk_cyto", {})
    rho_frac = summary["patient_tests"].get("spearman_mal_CLDN4__tnk_frac", {})

    def tp(key: str) -> str:
        d = summary["patient_tests"].get(key, {})
        if "rho" in d:
            return f"ρ={fmt_num(d.get('rho', np.nan))} p={fmt_p(d.get('p', np.nan))} (n={d.get('n', 'NA')})"
        return f"p={fmt_p(d.get('p', np.nan))}"

    signed_ok = (
        np.isfinite(rho_ifn.get("p", np.nan))
        and rho_ifn.get("p", 1) < 0.05
        and rho_ifn.get("rho", 0) < 0
    ) or (
        np.isfinite(rho_cy.get("p", np.nan))
        and rho_cy.get("p", 1) < 0.05
        and rho_cy.get("rho", 0) < 0
    )
    if signed_ok:
        verdict = (
            f"On the winning pair (eligible n={n_elig}: {n_131} GSE131907 + {n_205} GSE205335), "
            f"patient malignant CLDN4 vs T/NK IFN is {tp('spearman_mal_CLDN4__tnk_ifn')}; "
            f"vs cytotoxicity {tp('spearman_mal_CLDN4__tnk_cyto')}. "
            "Ligand ranks below are prior + multi-sample DE, not a causal induction claim."
        )
    else:
        verdict = (
            f"The signed T/NK-program claim is **not supported** on the winning pair "
            f"(eligible n=**{n_elig}**: {n_131} GSE131907 + {n_205} GSE205335). "
            f"Malignant CLDN4 vs T/NK IFN {tp('spearman_mal_CLDN4__tnk_ifn')}; "
            f"vs cytotoxicity {tp('spearman_mal_CLDN4__tnk_cyto')}; "
            f"vs T/NK fraction {tp('spearman_mal_CLDN4__tnk_frac')}. "
            "The ligand-activity table is a NicheNet-v2 prior ranking plus patient-paired "
            "CLDN4-high vs low ligand DE. Prior recovery of an IFN list is **not** evidence "
            "that CLDN4-high cells induce or repress IFN in these patients."
        )

    prim = table[table["geneset"].isin(PRIMARY_SETS)].copy() if not table.empty else table
    ifn = (
        prim[prim["geneset"] == "a_priori_ifn"]
        .sort_values("prioritization", ascending=False)
        .head(12)
        if not prim.empty
        else prim
    )
    cyto = (
        prim[prim["geneset"] == "a_priori_cytotoxicity"]
        .sort_values("prioritization", ascending=False)
        .head(12)
        if not prim.empty
        else prim
    )
    act_ifn = (
        prim[prim["geneset"] == "a_priori_ifn"].sort_values("pearson", ascending=False).head(8)
        if not prim.empty
        else prim
    )

    lines = [
        "# FINDING — MultiNicheNet / multi-sample NicheNet on winning pair GSE131907+GSE205335",
        "",
        f"**Verdict:** {verdict}",
        "",
        "ADDITIVE. **CLDN4 only.** No dual-high. **GSE207422 was not used** (NicheNet #334 NS). "
        "Sender = CLDN4-high vs low malignant. Receiver = **same-patient T/NK**. "
        "**Patient is the unit.** `nichenetr` / `multinichenetr` were not installed; "
        "ligand activity is the published NicheNet-v2 ligand–target prior (Zenodo 7074291).",
        "",
        "---",
        "",
        "## English",
        "",
        "### Honest n",
        "",
        "| Item | n | Note |",
        "| --- | ---: | --- |",
    ]
    for rec in ntab.itertuples(index=False):
        lines.append(f"| {rec.item} | {rec.n} | {rec.note} |")
    lines += [
        "",
        "GSE131907 is treatment-naive (Kim et al., *Nat Commun* 2020). It has **no ICI / MPR**. "
        "GSE205335 is palliative ICI biopsy/effusion (Hu / Ahn / Lee). RECIST is recorded and "
        "**not** used as MPR. Tumor samples are pooled per patient. "
        "A patient enters the paired ligand DE if it has ≥10 CLDN4-high, ≥10 CLDN4-low malignant "
        "cells and ≥20 T/NK.",
        "",
        "### E1 / E2 — patient T/NK programs vs malignant CLDN4",
        "",
        "Unit = eligible patient. Cells are not n.",
        "",
        "| Test | n | Result |",
        "| --- | ---: | --- |",
    ]
    for rec in tests.itertuples(index=False):
        if rec.test.startswith("spearman") or rec.test.startswith("paired") or rec.test.startswith("Q4"):
            if "rho" in tests.columns and pd.notna(getattr(rec, "rho", np.nan)):
                lines.append(
                    f"| {rec.test} | {int(rec.n) if pd.notna(rec.n) else 'NA'} | "
                    f"ρ={fmt_num(rec.rho)} p={fmt_p(rec.p)} |"
                )
            elif "median_delta" in tests.columns and pd.notna(getattr(rec, "median_delta", np.nan)):
                lines.append(
                    f"| {rec.test} | {int(rec.n) if pd.notna(rec.n) else 'NA'} | "
                    f"Δ={fmt_num(rec.median_delta)} p={fmt_p(rec.p)} |"
                )
            elif "p" in tests.columns:
                extra = ""
                if "mean_a" in tests.columns and pd.notna(getattr(rec, "mean_a", np.nan)):
                    extra = f" mean {fmt_num(rec.mean_a)} vs {fmt_num(rec.mean_b)};"
                lines.append(
                    f"| {rec.test} | {getattr(rec, 'n', getattr(rec, 'n_a', 'NA'))} | "
                    f"{extra} p={fmt_p(rec.p)} |"
                )
    lines += [
        "",
        "### E3 — ligand activity (unsigned NicheNet-v2 prior)",
        "",
        "Primary sender = CLDN4-high malignant (CLDN4 only); receiver = same-patient T/NK. "
        "Rank by Pearson of prior target scores vs gene-set membership on T/NK-expressed "
        "background genes (panel ∩ prior; **not** a full-transcriptome `nichenetr` run).",
        "",
        "Full table: [`results/ligand_activity_table.tsv`](results/ligand_activity_table.tsv) "
        "(also `ligand_activity_all.tsv`).",
        "",
        "**A priori IFN. Top 8 by Pearson (combined potential ligands):**",
        "",
    ]
    lines += md_table(
        act_ifn,
        ["rank_activity", "ligand", "pearson", "auroc", "pearson_p", "n_patients_de"],
        {
            "pearson": "{:.3f}",
            "auroc": "{:.3f}",
            "pearson_p": "{:.2e}",
            "rank_activity": "{:.0f}",
            "n_patients_de": "{:.0f}",
        },
    )
    lines += [
        "**MultiNicheNet-style prioritization (IFN set).** Combines scaled ligand activity, "
        "patient-paired CLDN4-high vs low ligand logFC, and T/NK receptor coverage. Top 12:",
        "",
    ]
    lines += md_table(
        ifn,
        [
            "rank_priority",
            "ligand",
            "pearson",
            "auroc",
            "mean_delta_high_minus_low",
            "paired_p",
            "frac_patients_receptor",
            "prioritization",
        ],
        {
            "rank_priority": "{:.0f}",
            "pearson": "{:.3f}",
            "auroc": "{:.3f}",
            "mean_delta_high_minus_low": "{:+.3f}",
            "paired_p": "{:.3g}",
            "frac_patients_receptor": "{:.2f}",
            "prioritization": "{:.3f}",
        },
    )
    lines += [
        "**A priori cytotoxicity. Top 12 by prioritization:**",
        "",
    ]
    lines += md_table(
        cyto,
        [
            "rank_priority",
            "ligand",
            "pearson",
            "auroc",
            "mean_delta_high_minus_low",
            "paired_p",
            "prioritization",
        ],
        {
            "rank_priority": "{:.0f}",
            "pearson": "{:.3f}",
            "auroc": "{:.3f}",
            "mean_delta_high_minus_low": "{:+.3f}",
            "paired_p": "{:.3g}",
            "prioritization": "{:.3f}",
        },
    )
    lines += [
        "A high Pearson on the IFN list means the ligand sits next to IFN genes in the "
        "published prior (ISG / MHC ligands). It does **not** mean IFN is up or down in "
        "CLDN4-high patients. Direction is E1/E2 only.",
        "",
        "### E4 — multi-sample (patient-paired) ligand DE",
        "",
        "For each potential ligand, mean `log1p(CP10k)` in CLDN4-high vs CLDN4-low malignant "
        "cells is compared with a **paired Wilcoxon across patients** (combined eligible set, "
        "and each dataset alone). That is the MultiNicheNet sender-DE term. "
        "Receivers are the same patient's T/NK, so receiver gene-set membership does not "
        "split by sender state within a patient; receiver association is the between-patient "
        "Spearman of T/NK gene means vs malignant CLDN4 (empirical sets).",
        "",
        "### What this is not",
        "",
        "- Not GSE207422 and not a re-run of PR #334.",
        "- Not TACSTD2∩CLDN4 dual-high senders.",
        "- Not `nichenetr` / `multinichenetr` R (packages absent; documented v2 prior used).",
        "- Not CellChat / LIANA communication probability.",
        "- Not inferCNV/CopyKAT recomputed malignant IDs (author labels).",
        "- Not ICI / MPR on GSE131907 (none labeled). RECIST on GSE205335 is not MPR.",
        "- Cells are not the sample size.",
        "",
        "### Files",
        "",
        "`results/ligand_activity_table.tsv`, `ligand_activity_all.tsv`, `top_ligands_primary.tsv`, "
        "`ligand_de_paired.tsv`, `patient_level_tests.tsv`, `n_table.tsv`, `sample_metrics.tsv`, "
        "`fig1`–`fig7`.",
        "",
        "---",
        "",
        "## 中文",
        "",
        f"**结论：** {verdict}",
        "",
        f"只做胜出对 **GSE131907 + GSE205335**。**不用 GSE207422**（#334 NS）。"
        f"发送端只用 CLDN4（恶性细胞全局中位数拆高/低），不用 TACSTD2 双高。"
        f"接收端是**同一患者** T/NK。推断单位是患者："
        f"**合格 n={n_elig}（131907: {n_131}；205335: {n_205}）**。",
        "配体活性是公开 NicheNet-v2 先验（Pearson / AUROC），加上患者配对的高 vs 低配体 DE "
        "（MultiNicheNet 式排序）。先验回收 IFN 基因集不是本队列里 CLDN4 诱导/抑制 IFN 的证据。"
        "GSE131907 无 ICI/MPR；GSE205335 的 RECIST 不代替 MPR。",
        "",
    ]
    text = "\n".join(lines)
    (ROOT / "FINDING.md").write_text(text)
    (OUT / "FINDING.md").write_text(text)


def main() -> int:
    pat = pd.read_csv(CACHE / "patient_inventory.tsv", sep="\t")
    stats_long = pd.read_parquet(CACHE / "patient_gene_stats.parquet")
    extract_meta = json.loads((CACHE / "extract_meta.json").read_text())
    lr = pd.read_csv(DATA / "lr_network.tsv", sep="\t")
    lt = pd.read_parquet(CACHE / "prior_ligand_target.parquet")
    if lt.shape[0] < lt.shape[1]:
        log(f"warning: lt shape {lt.shape} (expected genes x ligands)")

    elig = pat[pat["eligible_paired"]].copy()
    log(f"eligible patients: {len(elig)} / {len(pat)}")
    print(elig.groupby("dataset").size().to_string(), flush=True)

    # Patient program scores
    ifn_tnk = program_score(stats_long, IFN, "tnk").rename("tnk_ifn")
    cyto_tnk = program_score(stats_long, CYTOTOXICITY, "tnk").rename("tnk_cyto")
    exh_tnk = program_score(stats_long, EXHAUSTION, "tnk").rename("tnk_exh")
    sample = elig.set_index(["dataset", "patient"]).join([ifn_tnk, cyto_tnk, exh_tnk]).reset_index()
    sample["tnk_ifn_minus_cyto"] = sample["tnk_ifn"] - sample["tnk_cyto"]
    sample.to_csv(OUT / "sample_metrics.tsv", sep="\t", index=False)

    tests = []
    for ds, sub in list(sample.groupby("dataset")) + [("combined", sample)]:
        for x, y in [
            ("mal_CLDN4_mean", "tnk_ifn"),
            ("mal_CLDN4_mean", "tnk_cyto"),
            ("mal_CLDN4_mean", "tnk_ifn_minus_cyto"),
            ("mal_CLDN4_mean", "tnk_exh"),
            ("mal_CLDN4_mean", "tnk_frac"),
            ("mal_CLDN4_pct_pos", "tnk_ifn"),
            ("mal_CLDN4_pct_pos", "tnk_cyto"),
            ("mal_CLDN4_pct_pos", "tnk_frac"),
        ]:
            if x not in sub.columns or y not in sub.columns:
                continue
            sp = spearman_safe(sub[x], sub[y])
            tests.append({"test": f"spearman_{x}__{y}", "dataset": ds, **sp})
        # Q4 vs Q1 on T/NK fraction and programs
        if len(sub) >= 8 and "mal_CLDN4_mean" in sub.columns:
            ranks = sub["mal_CLDN4_mean"].rank(method="average")
            try:
                qs = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
            except ValueError:
                qs = None
            if qs is not None and qs.nunique() == 4:
                for metric in ["tnk_ifn", "tnk_cyto", "tnk_frac"]:
                    w = mwu(sub.loc[qs == "Q1", metric], sub.loc[qs == "Q4", metric])
                    tests.append({"test": f"Q4_vs_Q1_{metric}", "dataset": ds, **w, "n": w["n_a"] + w["n_b"]})

    # Potential ligands: expressed in >=10% of CLDN4-high cells in the eligible set
    # Use patient-weighted mean of frac_pos (unweighted across patients, then require mean>=0.10
    # in at least one dataset OR combined mean>=0.10). Also require a T/NK receptor.
    ligands_all = sorted(set(lr["from"].astype(str)) & set(lt.columns.astype(str)))
    receptors_all = sorted(set(lr["to"].astype(str)))

    def expressed(compartment: str, genes: list[str], dataset: str | None = None) -> list[str]:
        sub = stats_long[stats_long["compartment"] == compartment]
        sub = sub[sub["gene"].isin(genes)]
        keys = elig[["dataset", "patient"]]
        sub = sub.merge(keys, on=["dataset", "patient"], how="inner")
        if dataset:
            sub = sub[sub["dataset"] == dataset]
        if sub.empty:
            return []
        # cell-weighted detection: sum(frac*n)/sum(n)
        g = sub.groupby("gene").apply(
            lambda d: float(np.average(d["frac_pos"], weights=d["n_cells"])) if d["n_cells"].sum() else 0.0,
            include_groups=False,
        )
        return sorted(g[g >= DETECT_FRAC].index.astype(str))

    expr_lig_high = expressed("mal_high", ligands_all)
    expr_lig_low = expressed("mal_low", ligands_all)
    expr_rec_tnk = expressed("tnk", receptors_all)
    potential = sorted(
        set(lr.loc[lr["from"].isin(expr_lig_high) & lr["to"].isin(expr_rec_tnk), "from"])
    )
    pot_by_ds = {}
    for ds in DATASETS:
        h = expressed("mal_high", ligands_all, ds)
        r = expressed("tnk", receptors_all, ds)
        pot_by_ds[ds] = sorted(set(lr.loc[lr["from"].isin(h) & lr["to"].isin(r), "from"]))
    log(
        f"expressed ligands high={len(expr_lig_high)} low={len(expr_lig_low)} "
        f"T/NK receptors={len(expr_rec_tnk)} potential={len(potential)}"
    )

    # Background: T/NK-expressed genes in prior
    tnk_expr = expressed("tnk", sorted(set(lt.index.astype(str))))
    background = sorted(set(tnk_expr) | set(CYTOTOXICITY) | set(IFN) | set(EXHAUSTION))
    background = [g for g in background if g in lt.index]
    log(f"background genes {len(background)}")

    # Empirical T/NK genes vs malignant CLDN4 (patient Spearman, combined)
    tnk_wide = stats_long[stats_long["compartment"] == "tnk"].merge(
        elig[["dataset", "patient", "mal_CLDN4_mean"]], on=["dataset", "patient"], how="inner"
    )
    emp_rows = []
    for gene, sub in tnk_wide.groupby("gene"):
        sp = spearman_safe(sub["mal_CLDN4_mean"], sub["mean_log1p"])
        emp_rows.append({"gene": gene, **sp})
    emp = pd.DataFrame(emp_rows).sort_values("p")
    emp["in_ifn"] = emp["gene"].isin(IFN)
    emp["in_cyto"] = emp["gene"].isin(CYTOTOXICITY)
    emp["in_exh"] = emp["gene"].isin(EXHAUSTION)
    emp.to_csv(OUT / "tnk_gene_vs_mal_CLDN4.tsv", sep="\t", index=False)
    emp_up = set(emp.loc[(emp["p"] < EMPIRICAL_P) & (emp["rho"] > 0), "gene"])
    emp_down = set(emp.loc[(emp["p"] < EMPIRICAL_P) & (emp["rho"] < 0), "gene"])

    sets = {
        "a_priori_ifn": set(IFN) & set(lt.index),
        "a_priori_cytotoxicity": set(CYTOTOXICITY) & set(lt.index),
        "a_priori_exhaustion_sensitivity": set(EXHAUSTION) & set(lt.index),
        "empirical_TNK_up_with_CLDN4": emp_up & set(lt.index),
        "empirical_TNK_down_with_CLDN4": emp_down & set(lt.index),
    }

    activities = []
    settings = [("combined", potential)] + [(f"potential_{ds}", pot_by_ds[ds]) for ds in DATASETS]
    for set_name, gs in sets.items():
        for setting, ligs in settings:
            act = ligand_activity(lt, gs, background, ligs)
            if act.empty:
                continue
            act.insert(0, "setting", setting)
            act.insert(1, "geneset", set_name)
            activities.append(act)
    act_all = pd.concat(activities, ignore_index=True) if activities else pd.DataFrame()
    if not act_all.empty:
        act_all.to_csv(OUT / "ligand_activity_all.tsv", sep="\t", index=False)

    # Paired ligand DE (high vs low) per patient
    hi = stats_long[stats_long["compartment"] == "mal_high"][["dataset", "patient", "gene", "mean_log1p", "frac_pos", "n_cells"]]
    lo = stats_long[stats_long["compartment"] == "mal_low"][["dataset", "patient", "gene", "mean_log1p", "frac_pos", "n_cells"]]
    paired = hi.merge(lo, on=["dataset", "patient", "gene"], suffixes=("_high", "_low"))
    paired = paired.merge(elig[["dataset", "patient"]], on=["dataset", "patient"], how="inner")
    paired["delta"] = paired["mean_log1p_high"] - paired["mean_log1p_low"]
    focus_ligands = sorted(set(potential) | set().union(*pot_by_ds.values()))
    de_rows = []
    for scope, dsub in [("combined", paired)] + [(ds, paired[paired["dataset"] == ds]) for ds in DATASETS]:
        for L in focus_ligands:
            sub = dsub[dsub["gene"] == L]
            if sub.empty:
                continue
            w = wilcoxon_paired(sub["mean_log1p_high"], sub["mean_log1p_low"])
            de_rows.append({"dataset": scope, "ligand": L, **w})
    de = pd.DataFrame(de_rows)
    if not de.empty:
        de["padj"] = de.groupby("dataset")["p"].transform(bh_fdr)
        de = de.sort_values(["dataset", "p"])
        de.to_csv(OUT / "ligand_de_paired.tsv", sep="\t", index=False)

    # Receptor coverage: fraction of eligible patients with any receptor expressed in T/NK (>=10%)
    rec_map = lr.groupby("from")["to"].apply(lambda s: sorted(set(s.astype(str)))).to_dict()
    tnk_rec = stats_long[
        (stats_long["compartment"] == "tnk") & (stats_long["gene"].isin(receptors_all))
    ].merge(elig[["dataset", "patient"]], on=["dataset", "patient"], how="inner")
    rec_cov = {}
    n_pat = int(elig[["dataset", "patient"]].drop_duplicates().shape[0])
    for L in focus_ligands:
        recs = rec_map.get(L, [])
        if not recs or tnk_rec.empty:
            rec_cov[L] = 0.0
            continue
        sub = tnk_rec[tnk_rec["gene"].isin(recs)]
        if sub.empty:
            rec_cov[L] = 0.0
            continue
        hit = sub.groupby(["dataset", "patient"])["frac_pos"].max() >= DETECT_FRAC
        rec_cov[L] = float(hit.mean()) if len(hit) else 0.0

    # Build primary ligand-activity table (combined, a priori sets)
    table_rows = []
    if not act_all.empty:
        prim_act = act_all[act_all["setting"] == "combined"]
        de_comb = de[de["dataset"] == "combined"] if not de.empty else pd.DataFrame()
        for rec in prim_act.itertuples(index=False):
            L = rec.ligand
            drow = de_comb[de_comb["ligand"] == L]
            table_rows.append(
                {
                    "ligand": L,
                    "geneset": rec.geneset,
                    "pearson": rec.pearson,
                    "pearson_p": rec.pearson_p,
                    "auroc": rec.auroc,
                    "aupr": rec.aupr,
                    "n_background": rec.n_background,
                    "n_geneset": rec.n_geneset,
                    "mean_delta_high_minus_low": float(drow["mean_delta"].iloc[0]) if len(drow) else np.nan,
                    "median_delta_high_minus_low": float(drow["median_delta"].iloc[0]) if len(drow) else np.nan,
                    "paired_p": float(drow["p"].iloc[0]) if len(drow) else np.nan,
                    "paired_padj": float(drow["padj"].iloc[0]) if len(drow) else np.nan,
                    "n_patients_de": int(drow["n"].iloc[0]) if len(drow) else 0,
                    "frac_patients_receptor": rec_cov.get(L, np.nan),
                    "potential_GSE131907": L in pot_by_ds["GSE131907"],
                    "potential_GSE205335": L in pot_by_ds["GSE205335"],
                }
            )
    table = pd.DataFrame(table_rows)
    if not table.empty:
        chunks = []
        for gs, sub in table.groupby("geneset"):
            sub = sub.copy()
            sub["score_activity"] = minmax(sub["pearson"])
            sub["score_de"] = minmax(sub["mean_delta_high_minus_low"])
            sub["score_receptor"] = minmax(sub["frac_patients_receptor"])
            # MultiNicheNet-like: average of scaled activity, signed DE, receptor coverage
            sub["prioritization"] = sub[["score_activity", "score_de", "score_receptor"]].mean(axis=1)
            sub["rank_activity"] = sub["pearson"].rank(ascending=False, method="min")
            sub["rank_priority"] = sub["prioritization"].rank(ascending=False, method="min")
            chunks.append(sub)
        table = pd.concat(chunks, ignore_index=True)
        table = table.sort_values(["geneset", "rank_priority"])
        table.to_csv(OUT / "ligand_activity_table.tsv", sep="\t", index=False)
        primary = table[table["geneset"].isin(PRIMARY_SETS)].copy()
        primary.to_csv(OUT / "top_ligands_primary.tsv", sep="\t", index=False)
        top = []
        for gs, sub in primary.groupby("geneset"):
            chunk = sub.sort_values("pearson", ascending=False).head(TOP_N).copy()
            top.append(chunk)
        if top:
            pd.concat(top, ignore_index=True).to_csv(OUT / "top_ligands.tsv", sep="\t", index=False)

    tests_df = pd.DataFrame(tests)
    tests_df.to_csv(OUT / "patient_level_tests.tsv", sep="\t", index=False)

    n_tbl = pd.DataFrame(
        [
            {
                "item": "datasets",
                "n": 2,
                "note": "GSE131907 + GSE205335 only; GSE207422 excluded",
            },
            {
                "item": "GSE131907_cells",
                "n": extract_meta["GSE131907"]["n_cells"],
                "note": "Kim et al. author barcodes",
            },
            {
                "item": "GSE131907_malignant",
                "n": extract_meta["GSE131907"]["n_malignant"],
                "note": "tumor-origin ∩ {Malignant cells, tS1, tS2, tS3}",
            },
            {
                "item": "GSE131907_cldn4_high_low",
                "n": extract_meta["GSE131907"]["n_cldn4_high"],
                "note": f"high; low={extract_meta['GSE131907']['n_cldn4_low']}; median={extract_meta['GSE131907']['cldn4_threshold']:.3f}",
            },
            {
                "item": "GSE131907_TNK",
                "n": extract_meta["GSE131907"]["n_tnk"],
                "note": "tumor-origin T lymphocytes + NK",
            },
            {
                "item": "GSE205335_cells",
                "n": extract_meta["GSE205335"]["n_cells"],
                "note": "Hu/Ahn/Lee processed UMI",
            },
            {
                "item": "GSE205335_malignant",
                "n": extract_meta["GSE205335"]["n_malignant"],
                "note": "tumor sample ∩ lineage.sub == Malignant cells",
            },
            {
                "item": "GSE205335_cldn4_high_low",
                "n": extract_meta["GSE205335"]["n_cldn4_high"],
                "note": f"high; low={extract_meta['GSE205335']['n_cldn4_low']}; median={extract_meta['GSE205335']['cldn4_threshold']:.3f}",
            },
            {
                "item": "GSE205335_TNK",
                "n": extract_meta["GSE205335"]["n_tnk"],
                "note": "tumor sample ∩ lineage.total == T/NK cells",
            },
            {
                "item": "eligible_patients_GSE131907",
                "n": int((elig["dataset"] == "GSE131907").sum()),
                "note": "≥10 high, ≥10 low, ≥20 T/NK; patient pooled",
            },
            {
                "item": "eligible_patients_GSE205335",
                "n": int((elig["dataset"] == "GSE205335").sum()),
                "note": "same gate; RECIST not used as MPR",
            },
            {
                "item": "eligible_patients_combined",
                "n": int(len(elig)),
                "note": "unit of inference",
            },
            {
                "item": "potential_ligands_combined",
                "n": int(len(potential)),
                "note": "CLDN4-high ligand + T/NK receptor + v2 prior",
            },
            {
                "item": "potential_ligands_GSE131907",
                "n": int(len(pot_by_ds["GSE131907"])),
                "note": "dataset-specific 10% rule",
            },
            {
                "item": "potential_ligands_GSE205335",
                "n": int(len(pot_by_ds["GSE205335"])),
                "note": "dataset-specific 10% rule",
            },
            {
                "item": "background_genes",
                "n": int(len(background)),
                "note": "T/NK-expressed ∩ prior targets (panel-restricted)",
            },
            {
                "item": "a_priori_ifn_in_prior",
                "n": int(len(sets["a_priori_ifn"])),
                "note": "",
            },
            {
                "item": "a_priori_cyto_in_prior",
                "n": int(len(sets["a_priori_cytotoxicity"])),
                "note": "",
            },
            {
                "item": "empirical_TNK_up_with_CLDN4",
                "n": int(len(emp_up)),
                "note": f"patient Spearman p<{EMPIRICAL_P}, ρ>0",
            },
            {
                "item": "empirical_TNK_down_with_CLDN4",
                "n": int(len(emp_down)),
                "note": f"patient Spearman p<{EMPIRICAL_P}, ρ<0",
            },
        ]
    )
    n_tbl.to_csv(OUT / "n_table.tsv", sep="\t", index=False)

    def grab(test: str, dataset: str = "combined") -> dict:
        hit = tests_df[(tests_df["test"] == test) & (tests_df["dataset"] == dataset)]
        if hit.empty:
            return {}
        return {k: (None if pd.isna(v) else v) for k, v in hit.iloc[0].to_dict().items() if k not in {"test", "dataset"}}

    summary = {
        "datasets": ["GSE131907", "GSE205335"],
        "excluded": {"GSE207422": "PR #334 NS; not re-run; not added"},
        "sender_definition": "CLDN4-only malignant high (log1p CP10k >= dataset malignant median)",
        "receiver_definition": "same-patient T/NK (author labels)",
        "unit": "patient",
        "method": "documented NicheNet-v2 ligand-target prior + MultiNicheNet-like prioritization",
        "packages_missing": ["nichenetr", "multinichenetr"],
        "n": {
            "eligible_patients_combined": int(len(elig)),
            "eligible_GSE131907": int((elig["dataset"] == "GSE131907").sum()),
            "eligible_GSE205335": int((elig["dataset"] == "GSE205335").sum()),
            "potential_ligands": int(len(potential)),
            "background_genes": int(len(background)),
        },
        "cutoffs": {
            "GSE131907_median_CLDN4": extract_meta["GSE131907"]["cldn4_threshold"],
            "GSE205335_median_CLDN4": extract_meta["GSE205335"]["cldn4_threshold"],
            "detect_frac": DETECT_FRAC,
        },
        "patient_tests": {
            "spearman_mal_CLDN4__tnk_ifn": grab("spearman_mal_CLDN4_mean__tnk_ifn"),
            "spearman_mal_CLDN4__tnk_cyto": grab("spearman_mal_CLDN4_mean__tnk_cyto"),
            "spearman_mal_CLDN4__tnk_frac": grab("spearman_mal_CLDN4_mean__tnk_frac"),
            "spearman_mal_CLDN4_pct__tnk_frac": grab("spearman_mal_CLDN4_pct_pos__tnk_frac"),
        },
        "empirical_sets": {
            "up": sorted(emp_up),
            "down": sorted(emp_down),
            "p_cutoff": EMPIRICAL_P,
        },
        "caveats": [
            "GSE207422 not used",
            "No dual-high / TACSTD2 gate",
            "Prior is unsigned; direction is patient-level only",
            "Background is panel-restricted",
            "Author malignant labels, not CNV re-call",
            "GSE131907 has no ICI/MPR; GSE205335 RECIST is not MPR",
        ],
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    # -------- figures (extra) --------
    colors_ds = {"GSE131907": "#1b9e77", "GSE205335": "#d95f02"}

    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.6))
    for ax, col, title in [
        (axes[0], "n_malignant", "Malignant cells / patient"),
        (axes[1], "n_tnk", "T/NK cells / patient"),
        (axes[2], "mal_CLDN4_mean", "Malignant CLDN4 mean"),
    ]:
        for i, ds in enumerate(DATASETS):
            sub = elig[elig["dataset"] == ds]
            rng = np.random.default_rng(1)
            ax.scatter(
                rng.uniform(-0.12, 0.12, len(sub)) + i,
                sub[col],
                c=colors_ds[ds],
                s=28,
                alpha=0.85,
                label=f"{ds} n={len(sub)}",
            )
        ax.set_xticks([0, 1], ["GSE131907", "GSE205335"])
        ax.set_title(title, fontsize=9)
        ax.legend(fontsize=6, frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "fig1_n_cohort.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2))
    for ax, y, title in [
        (axes[0], "tnk_ifn", "T/NK IFN vs malignant CLDN4"),
        (axes[1], "tnk_cyto", "T/NK cytotoxicity vs malignant CLDN4"),
    ]:
        for ds, sub in sample.groupby("dataset"):
            ax.scatter(sub["mal_CLDN4_mean"], sub[y], c=colors_ds[ds], s=36, label=f"{ds} n={len(sub)}")
        ax.set_xlabel("Malignant CLDN4 mean log1p(CP10k)")
        ax.set_ylabel(f"{y} mean log1p(CP10k)")
        ax.set_title(title, fontsize=9)
        ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "fig2_cldn4_vs_tnk_programs.png", dpi=160)
    plt.close(fig)

    if not table.empty:
        fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.4))
        for ax, gs, title in [
            (axes[0], "a_priori_ifn", "Prior activity vs IFN set"),
            (axes[1], "a_priori_cytotoxicity", "Prior activity vs cytotoxicity set"),
        ]:
            sub = table[table["geneset"] == gs].sort_values("pearson", ascending=False).head(12)
            if sub.empty:
                ax.axis("off")
                continue
            ax.barh(sub["ligand"][::-1], sub["pearson"][::-1], color="#4C72B0")
            ax.set_xlabel("Pearson (prior vs gene-set membership)")
            ax.set_title(title, fontsize=9)
        fig.tight_layout()
        fig.savefig(OUT / "fig3_ligand_activity_top.png", dpi=160)
        plt.close(fig)

        fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.4))
        for ax, gs, title in [
            (axes[0], "a_priori_ifn", "MultiNicheNet-like priority (IFN)"),
            (axes[1], "a_priori_cytotoxicity", "MultiNicheNet-like priority (cytotoxicity)"),
        ]:
            sub = table[table["geneset"] == gs].sort_values("prioritization", ascending=False).head(12)
            if sub.empty:
                ax.axis("off")
                continue
            ax.barh(sub["ligand"][::-1], sub["prioritization"][::-1], color="#8c564b")
            ax.set_xlabel("Prioritization (scaled activity + DE + receptor)")
            ax.set_title(title, fontsize=9)
        fig.tight_layout()
        fig.savefig(OUT / "fig4_prioritization.png", dpi=160)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(6.2, 5.2))
        sub = table[table["geneset"] == "a_priori_ifn"]
        sc = ax.scatter(
            sub["mean_delta_high_minus_low"],
            sub["pearson"],
            c=sub["frac_patients_receptor"],
            s=28,
            cmap="viridis",
        )
        fig.colorbar(sc, ax=ax, label="T/NK receptor coverage")
        show = sub.reindex(sub["prioritization"].sort_values(ascending=False).index).head(10)
        for r in show.itertuples():
            ax.annotate(r.ligand, (r.mean_delta_high_minus_low, r.pearson), fontsize=7)
        ax.axvline(0, color="0.7", lw=0.8)
        ax.set_xlabel("Paired ligand Δ (CLDN4-high − low)")
        ax.set_ylabel("IFN prior Pearson")
        ax.set_title("Extra: activity vs sender DE (IFN set)", fontsize=9)
        fig.tight_layout()
        fig.savefig(OUT / "fig5_activity_vs_de.png", dpi=160)
        plt.close(fig)

    if not de.empty:
        a = de[de["dataset"] == "GSE131907"][["ligand", "mean_delta"]].rename(columns={"mean_delta": "d_131907"})
        b = de[de["dataset"] == "GSE205335"][["ligand", "mean_delta"]].rename(columns={"mean_delta": "d_205335"})
        both = a.merge(b, on="ligand", how="inner")
        if not both.empty:
            fig, ax = plt.subplots(figsize=(5.4, 5.2))
            ax.scatter(both["d_131907"], both["d_205335"], s=16, c="#555555")
            lim = [
                min(both["d_131907"].min(), both["d_205335"].min()) - 0.02,
                max(both["d_131907"].max(), both["d_205335"].max()) + 0.02,
            ]
            ax.plot(lim, lim, ls="--", c="0.7", lw=1)
            both = both.assign(delta=both["d_205335"] - both["d_131907"])
            for r in both.reindex(both["d_131907"].abs().sort_values(ascending=False).index).head(8).itertuples():
                ax.annotate(r.ligand, (r.d_131907, r.d_205335), fontsize=7)
            ax.set_xlabel("Ligand Δ high−low (GSE131907)")
            ax.set_ylabel("Ligand Δ high−low (GSE205335)")
            ax.set_title("Extra: dataset concordance of sender DE", fontsize=9)
            fig.tight_layout()
            fig.savefig(OUT / "fig6_dataset_concordance.png", dpi=160)
            plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    ax.scatter(sample["mal_CLDN4_mean"], sample["tnk_frac"], c=sample["dataset"].map(colors_ds), s=36)
    for ds, sub in sample.groupby("dataset"):
        ax.scatter([], [], c=colors_ds[ds], s=36, label=f"{ds} n={len(sub)}")
    ax.set_xlabel("Malignant CLDN4 mean log1p(CP10k)")
    ax.set_ylabel("T/NK fraction")
    ax.set_title("Extra: CLDN4 vs same-patient T/NK fraction", fontsize=9)
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "fig7_cldn4_vs_tnk_frac.png", dpi=160)
    plt.close(fig)

    write_finding(summary, n_tbl, tests_df, table if not table.empty else pd.DataFrame())
    log(f"wrote results to {OUT}")
    print(json.dumps(summary["n"], indent=2), flush=True)
    print(json.dumps(summary["patient_tests"], indent=2), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
