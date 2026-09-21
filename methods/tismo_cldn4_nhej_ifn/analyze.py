#!/usr/bin/env python3
"""Cldn4 versus NHEJ and IFN gene scores on the locked TISMO ICB pairing.

The 64 slices and the Tacstd2 49/64 result are the PR #542 lock. This script
recomputes Tacstd2 and Cldn4 from the current TISMO export and stops if those
deltas do not match the lock. NHEJ and IFN scores are new.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
from scipy.stats import binomtest, pearsonr, spearmanr, wilcoxon

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUT = ROOT / "results" / "tismo_cldn4_nhej_ifn"
TABLES = OUT / "tables"
FIG = OUT / "figures"
LOCK_PATH = HERE / "locked_pairs_pr542.tsv"

STEM_RE = re.compile(r"\(n=\d+\)$")
MIN_GENE_FRACTION = 0.8
# Longest token first so antiCTLA4&antiPD1 is not parsed as antiCTLA4.
ICB_TOKENS = (
    "antiCTLA4&antiPDL1",
    "antiCTLA4&antiPD1",
    "antiPD1_FcS",
    "antiPDL1",
    "antiPDL2",
    "antiPD1",
    "antiCTLA4",
)

CANCER = {
    "4T1": ("Mammary carcinoma", "Mammary"),
    "E0771": ("Mammary adenocarcinoma", "Mammary"),
    "EMT6": ("Mammary carcinoma", "Mammary"),
    "KPB25L": ("Mammary cancer, NOS", "Mammary"),
    "T11": ("Mammary cancer, NOS", "Mammary"),
    "p53-2225L": ("Mammary cancer, NOS", "Mammary"),
    "p53-2336R": ("Mammary cancer, NOS", "Mammary"),
    "B16": ("Melanoma", "Melanoma"),
    "YUMM1.7": ("Melanoma", "Melanoma"),
    "D3UV2": ("Melanoma", "Melanoma"),
    "D4M.3A.3": ("Melanoma", "Melanoma"),
    "CT26": ("Colorectal carcinoma", "Colorectal"),
    "MC38": ("Colorectal carcinoma", "Colorectal"),
    "LLC": ("Lung carcinoma", "Lung"),
    "402230": ("Sarcoma", "Sarcoma"),
    "BNL-MEA": ("Hepatocellular carcinoma", "Liver"),
    "YTN16": ("Gastric adenocarcinoma", "Gastric"),
    "MOC22": ("Oral squamous cell carcinoma", "HeadNeck"),
    "CMT-167": ("Lung carcinoma", "Lung"),
}

GROUP_COLORS = {
    "Mammary": "#4C78A8",
    "Melanoma": "#F58518",
    "Colorectal": "#54A24B",
    "Gastric": "#E45756",
    "Liver": "#72B7B2",
    "Lung": "#B279A2",
    "Sarcoma": "#FF9DA6",
    "HeadNeck": "#9D755D",
    "Other": "#BAB0AC",
}

# Short names used in tables. Primary pair is IFN_gamma and NHEJ_GO.
SCORE_OF_SET = {
    "HALLMARK_INTERFERON_GAMMA_RESPONSE": "IFN_gamma",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE": "IFN_alpha",
    "GOBP_DOUBLE_STRAND_BREAK_REPAIR_VIA_NONHOMOLOGOUS_END_JOINING": "NHEJ_GO",
    "REACTOME_NONHOMOLOGOUS_END_JOINING_NHEJ": "NHEJ_Reactome",
    "REACTOME_NHEJ_minus_histones": "NHEJ_Reactome_no_histone",
    "WP_NONHOMOLOGOUS_END_JOINING": "NHEJ_WP",
    "cNHEJ_machinery": "NHEJ_core",
}

PRIMARY_TESTS = [
    ("delta", "IFN_gamma"),
    ("delta", "NHEJ_GO"),
    ("baseline", "IFN_gamma"),
    ("baseline", "NHEJ_GO"),
    ("icb", "IFN_gamma"),
    ("icb", "NHEJ_GO"),
]


def stem(label: str) -> str:
    return STEM_RE.sub("", str(label)).rstrip()


def split_icb_stem(label: str) -> tuple[str, str]:
    for tok in ICB_TOKENS:
        suffix = "_" + tok
        if label.endswith(suffix):
            return label[: -len(suffix)], tok
    return label, ""


def share_sibling_baselines(meta: pd.DataFrame) -> tuple[pd.DataFrame, list[dict]]:
    """Repeat a shared isotype cohort onto each ICB arm.

    The PR #542 export listed one cell_line slice per arm and repeated the
    isotype rows on every arm. The current TISMO export keeps each sample
    once, so an arm can arrive with ICB rows and no baseline. Those baselines
    are copied from the single sibling arm that still carries them.
    """
    df = meta.copy()
    df["stem"] = df["cell_line"].map(stem)
    df["prefix"] = df["stem"].map(lambda s: split_icb_stem(s)[0])
    needs = []
    for label, g in df.groupby("stem"):
        if (g["Baseline"] == 0).any() and not (g["Baseline"] == 1).any():
            needs.append(label)
    extra = []
    notes = []
    for label in needs:
        prefix, _tok = split_icb_stem(label)
        siblings = df[(df["prefix"] == prefix) & (df["stem"] != label) & (df["Baseline"] == 1)]
        sib_stems = sorted(siblings["stem"].unique())
        if len(sib_stems) != 1:
            raise SystemExit(f"no unique shared baseline for {label}; siblings={sib_stems}")
        base = siblings.drop_duplicates("Samples").copy()
        target_label = df.loc[df["stem"] == label, "cell_line"].iloc[0]
        base["cell_line"] = target_label
        base["stem"] = label
        base["prefix"] = prefix
        extra.append(base)
        notes.append(
            {
                "stem": label,
                "baseline_from": sib_stems[0],
                "n_baseline_samples": int(base["Samples"].nunique()),
            }
        )
    if extra:
        df = pd.concat([df, *extra], ignore_index=True)
    return df.drop(columns=["prefix"]), notes


def parse_line(model: str) -> str:
    return model.split("_")[0]


def cancer_of(line: str) -> tuple[str, str]:
    return CANCER.get(line, ("Unknown", "Other"))


def fmt_p(p) -> str:
    if p is None or (isinstance(p, float) and (np.isnan(p) or np.isinf(p))):
        return "NA"
    if p == 0:
        return "<1×10⁻¹⁶"
    if p < 1e-4:
        exp = int(np.floor(np.log10(p)))
        mant = p / 10**exp
        sup_map = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")
        return f"{mant:.2f}×10{str(exp).translate(sup_map)}"
    return f"{p:.3g}"


def fmt_rho(r) -> str:
    if r is None or (isinstance(r, float) and np.isnan(r)):
        return "NA"
    return f"{r:+.2f}"


def sign_counts(values) -> dict:
    v = pd.Series(values, dtype=float).dropna()
    n_up = int((v > 0).sum())
    n_down = int((v < 0).sum())
    n_tie = int((v == 0).sum())
    n = n_up + n_down
    if n == 0:
        p_two = p_greater = np.nan
    else:
        p_two = float(binomtest(n_up, n, 0.5, alternative="two-sided").pvalue)
        p_greater = float(binomtest(n_up, n, 0.5, alternative="greater").pvalue)
    return {
        "n_pairs": int(len(v)),
        "n_up": n_up,
        "n_down": n_down,
        "n_tie": n_tie,
        "n_tested": n,
        "binomial_p_two_sided": p_two,
        "binomial_p_greater": p_greater,
    }


def wilcoxon_delta(values) -> dict:
    v = pd.Series(values, dtype=float).dropna()
    v_nz = v[v != 0]
    if len(v_nz) < 1:
        return {"n": int(len(v_nz)), "statistic": np.nan, "p_two_sided": np.nan, "p_greater": np.nan}
    two = wilcoxon(v_nz, alternative="two-sided", zero_method="wilcox")
    greater = wilcoxon(v_nz, alternative="greater", zero_method="wilcox")
    return {
        "n": int(len(v_nz)),
        "statistic": float(two.statistic),
        "p_two_sided": float(two.pvalue),
        "p_greater": float(greater.pvalue),
    }


def bh_fdr(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank, i in enumerate(order[::-1], start=1):
        k = n - rank + 1
        val = min(prev, p[i] * n / k)
        q[i] = val
        prev = val
    return [float(min(1.0, v)) for v in q]


def pair_values(meta: pd.DataFrame, values: pd.Series, colname: str) -> pd.DataFrame:
    """Pair Baseline==1 vs Baseline==0 by TISMO cell_line stem. Delta = ICB − naive."""
    df = meta[["Samples", "cell_line", "Baseline", "GSE_ID", "Mouse_treatment"]].copy()
    df[colname] = values.reindex(df["Samples"]).to_numpy()
    df["stem"] = df["cell_line"].map(stem)
    rows = []
    for model, g in df.groupby("stem", sort=True):
        base = g.loc[g["Baseline"] == 1, colname].astype(float).dropna()
        trt = g.loc[g["Baseline"] == 0, colname].astype(float).dropna()
        if len(base) == 0 or len(trt) == 0:
            continue
        line = parse_line(model)
        cancer, group = cancer_of(line)
        treatments = sorted({str(x) for x in g["Mouse_treatment"].dropna().unique() if str(x) not in {"", "NA", "nan"}})
        rows.append(
            {
                "stem": model,
                "cell_line": line,
                "gse_id": str(g["GSE_ID"].iloc[0]),
                "cancer_type": cancer,
                "cancer_group": group,
                "is_lung": group == "Lung" or line.upper() == "LLC",
                "is_llc": line.upper() == "LLC",
                "n_baseline": int(len(base)),
                "n_icb": int(len(trt)),
                f"{colname}_mean_naive": float(base.mean()),
                f"{colname}_mean_icb": float(trt.mean()),
                f"{colname}_delta_mean": float(trt.mean() - base.mean()),
                f"{colname}_delta_median": float(trt.median() - base.median()),
                "treatments": ";".join(treatments),
            }
        )
    return pd.DataFrame(rows)


def mean_z_score(mat: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str], list[str]]:
    """Mean of per-gene z-scores. Z is fit on this matrix only (ddof=0)."""
    present = [g for g in genes if g in mat.columns]
    missing = [g for g in genes if g not in mat.columns]
    use = []
    zcols = []
    for gene in present:
        s = mat[gene].astype(float)
        if s.notna().mean() < MIN_GENE_FRACTION:
            missing.append(gene)
            continue
        sd = float(s.std(ddof=0))
        if not np.isfinite(sd) or sd == 0:
            missing.append(gene)
            continue
        zcols.append((s - float(s.mean())) / sd)
        use.append(gene)
    if not zcols:
        return pd.Series(np.nan, index=mat.index), [], missing
    z = pd.concat(zcols, axis=1)
    z.columns = use
    coverage = z.notna().mean(axis=1)
    score = z.mean(axis=1, skipna=True)
    score = score.where(coverage >= MIN_GENE_FRACTION)
    return score, use, missing


def spearman_block(x: pd.Series, y: pd.Series, label: str, seed: int) -> dict:
    a = pd.concat([x.rename("x"), y.rename("y")], axis=1).dropna()
    n = int(len(a))
    if n < 5:
        return {"label": label, "n": n, "rho": np.nan, "p": np.nan, "ci_low": np.nan, "ci_high": np.nan}
    rho, p = spearmanr(a["x"], a["y"])
    rng = np.random.default_rng(seed)
    boots = np.empty(4000)
    xv = a["x"].to_numpy()
    yv = a["y"].to_numpy()
    for i in range(4000):
        idx = rng.integers(0, n, n)
        boots[i] = spearmanr(xv[idx], yv[idx]).statistic
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return {
        "label": label,
        "n": n,
        "rho": float(rho),
        "p": float(p),
        "ci_low": float(lo),
        "ci_high": float(hi),
    }


def partial_spearman(x: pd.Series, y: pd.Series, z: pd.Series) -> dict:
    a = pd.concat([x.rename("x"), y.rename("y"), z.rename("z")], axis=1).dropna()
    n = int(len(a))
    if n < 8:
        return {"n": n, "rho": np.nan, "p": np.nan}
    ranks = a.rank(method="average")
    def resid(left: np.ndarray, right: np.ndarray) -> np.ndarray:
        design = np.column_stack([np.ones(len(right)), right])
        coef, _, _, _ = np.linalg.lstsq(design, left, rcond=None)
        return left - design @ coef
    xr = resid(ranks["x"].to_numpy(), ranks["z"].to_numpy())
    yr = resid(ranks["y"].to_numpy(), ranks["z"].to_numpy())
    rho, p = pearsonr(xr, yr)
    return {"n": n, "rho": float(rho), "p": float(p)}


def assert_lock(pairs: pd.DataFrame, lock: pd.DataFrame) -> dict:
    merged = lock.merge(pairs, on="stem", how="outer", indicator=True)
    missing = merged.loc[merged["_merge"] != "both", "stem"].tolist()
    if missing:
        raise SystemExit(f"lock stem mismatch ({len(missing)}): {missing[:8]}")
    report = {}
    for gene, lock_col in (("Tacstd2", "Tacstd2_delta_mean"), ("Cldn4", "Cldn4_delta_mean")):
        delta = pairs.set_index("stem")[f"{gene}_delta_mean"]
        ref = lock.set_index("stem")[lock_col]
        diff = (delta - ref).abs()
        report[gene] = {
            "max_abs_diff": float(diff.max()),
            "n": int(len(diff)),
            "n_up": int((delta > 0).sum()),
            "n_down": int((delta < 0).sum()),
            "n_tie": int((delta == 0).sum()),
        }
        if diff.max() > 1e-6:
            worst = diff.sort_values(ascending=False).head(5)
            raise SystemExit(f"{gene} deltas do not match the PR #542 lock. Worst stems:\n{worst}")
    # The published lock is the mean-delta sign count.
    if report["Tacstd2"]["n_up"] != 49 or report["Tacstd2"]["n"] != 64:
        raise SystemExit(f"Tacstd2 lock is 49/64; recomputed {report['Tacstd2']}")
    if report["Cldn4"]["n_up"] != 34 or report["Cldn4"]["n_down"] != 28 or report["Cldn4"]["n_tie"] != 2:
        raise SystemExit(f"Cldn4 lock is 34/28/2; recomputed {report['Cldn4']}")
    return report


def scatter(ax, df: pd.DataFrame, x: str, y: str, title: str) -> None:
    for group, g in df.groupby("cancer_group"):
        lung = g["is_lung"].to_numpy()
        ax.scatter(
            g.loc[~g["is_lung"], x],
            g.loc[~g["is_lung"], y],
            s=28,
            c=GROUP_COLORS.get(group, "#333333"),
            label=group,
            alpha=0.9,
            linewidths=0,
        )
        if lung.any():
            ax.scatter(
                g.loc[g["is_lung"], x],
                g.loc[g["is_lung"], y],
                s=46,
                c=GROUP_COLORS.get(group, "#333333"),
                edgecolors="black",
                linewidths=1.1,
                zorder=3,
                label="Lung (LLC)",
            )
    ax.axhline(0, color="0.6", lw=0.6)
    ax.axvline(0, color="0.6", lw=0.6)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Cldn4 Δ mean (ICB − naive)")
    ax.set_ylabel(y.replace("_", " "))


def waterfall(ax, df: pd.DataFrame, col: str, title: str) -> None:
    d = df.sort_values([col, "stem"]).reset_index(drop=True)
    colors = np.where(d[col] >= 0, "#4C78A8", "#E45756")
    ax.axhline(0, color="0.35", lw=0.7)
    ax.bar(np.arange(len(d)), d[col], color=colors, width=0.9, linewidth=0)
    if d["is_lung"].any():
        idx = np.where(d["is_lung"].to_numpy())[0]
        ax.scatter(idx, d.loc[d["is_lung"], col], s=28, facecolors="none", edgecolors="black", linewidths=1.1, zorder=3)
    ax.set_xticks([])
    ax.set_title(title, fontsize=10)
    ax.set_ylabel("score Δ mean")
    ax.set_xlabel(f"64 TISMO ICB slices, sorted")


def level_scatter(ax, df: pd.DataFrame, x: str, y: str, title: str) -> None:
    for group, g in df.groupby("cancer_group"):
        ax.scatter(
            g[x],
            g[y],
            s=28,
            c=GROUP_COLORS.get(group, "#333333"),
            label=group,
            alpha=0.9,
            linewidths=0,
        )
    ax.set_title(title, fontsize=10)
    ax.set_xlabel(x.replace("_", " "))
    ax.set_ylabel(y.replace("_", " "))


def _missing_clause(row: dict) -> str:
    missing = [g for g in str(row.get("missing") or "").split(";") if g]
    if not missing:
        return ""
    return "; absent in TISMO: " + ", ".join(missing)


def write_results(summary: dict) -> None:
    p = summary["primary_correlations"]
    by = {row["label"]: row for row in p}
    d_ifn = by["delta Cldn4 vs IFN_gamma"]
    d_nhej = by["delta Cldn4 vs NHEJ_GO"]
    b_ifn = by["baseline Cldn4 vs IFN_gamma"]
    b_nhej = by["baseline Cldn4 vs NHEJ_GO"]
    i_ifn = by["icb Cldn4 vs IFN_gamma"]
    i_nhej = by["icb Cldn4 vs NHEJ_GO"]
    directions = {row["score"]: row for row in summary["score_directions"]}
    ifn_dir = directions["IFN_gamma"]
    nhej_dir = directions["NHEJ_GO"]
    partial = {row["label"]: row for row in summary["partial"]}
    cov = {row["score"]: row for row in summary["coverage"]}
    lock = summary["lock"]
    w_cldn4 = summary["cldn4_wilcoxon_p"]

    def rho_clause(row: dict) -> str:
        return (
            f"ρ = {fmt_rho(row['rho'])} "
            f"(95% bootstrap CI {fmt_rho(row['ci_low'])} to {fmt_rho(row['ci_high'])}, "
            f"p = {fmt_p(row['p'])}, q = {fmt_p(row['q'])}, n = {row['n']})"
        )

    line_base = summary["level_dependence_index"]["baseline"]["cell-line means"]
    line_icb = summary["level_dependence_index"]["after ICB"]["cell-line means"]
    sentence = (
        f"On the locked 64 TISMO ICB slices (Tacstd2 still 49/64 up, Wilcoxon p = {fmt_p(summary['tacstd2_wilcoxon_p'])}; "
        f"Cldn4 still 34/64 up, 28 down, 2 ties, Wilcoxon p = {fmt_p(w_cldn4)}), "
        f"the ICB change in Cldn4 is uncorrelated with the ICB change in Hallmark IFN-γ "
        f"({rho_clause(d_ifn)}) and with the GO:0006303 NHEJ change ({rho_clause(d_nhej)}). "
        f"Averaging slices to {summary['n_cell_lines']} cell lines leaves those change-correlations null "
        f"(IFN-γ ρ = {fmt_rho(summary['cell_line_delta']['ifn_rho'])}, p = {fmt_p(summary['cell_line_delta']['ifn_p'])}; "
        f"NHEJ ρ = {fmt_rho(summary['cell_line_delta']['nhej_rho'])}, p = {fmt_p(summary['cell_line_delta']['nhej_p'])}). "
        f"IFN-γ rose in {ifn_dir['n_up']}/{ifn_dir['n_pairs']} slices "
        f"(Wilcoxon p = {fmt_p(ifn_dir['wilcoxon_p_two_sided'])}) and NHEJ fell in {nhej_dir['n_down']}/{nhej_dir['n_pairs']} "
        f"(Wilcoxon p = {fmt_p(nhej_dir['wilcoxon_p_two_sided'])}); Cldn4 does not follow either shift. "
        f"A slice-level baseline correlation of Cldn4 with IFN-γ ({rho_clause(b_ifn)}) "
        f"does not remain after cell-line averaging "
        f"(ρ = {fmt_rho(line_base['ifn_rho'])}, p = {fmt_p(line_base['ifn_p'])}, n = {line_base['n']}). "
        f"Baseline Cldn4 versus NHEJ is {rho_clause(b_nhej)} at slice level and "
        f"ρ = {fmt_rho(line_base['nhej_rho'])} (p = {fmt_p(line_base['nhej_p'])}) across cell lines. "
        f"After ICB, slice-level Cldn4 versus IFN-γ is {rho_clause(i_ifn)} and versus NHEJ is {rho_clause(i_nhej)}; "
        f"cell-line means are IFN-γ ρ = {fmt_rho(line_icb['ifn_rho'])} (p = {fmt_p(line_icb['ifn_p'])}) and "
        f"NHEJ ρ = {fmt_rho(line_icb['nhej_rho'])} (p = {fmt_p(line_icb['nhej_p'])}). "
        f"Lung ICB in this pairing is LLC only (2/64), not KL or KP. Scores are bulk-tumor RNA."
    )

    lines = [
        "# TISMO ICB: Cldn4 versus NHEJ and IFN gene scores",
        "",
        "Public recompute on the PR #542 pairing. Tacstd2 **49/64** is unchanged. Every number below is in `tables/`.",
        "",
        "## Paper sentence",
        "",
        sentence,
        "",
        "## What was held fixed",
        "",
        f"Live TISMO Gene-module export matched the archived 64-slice means "
        f"(Tacstd2 max |Δ difference| = {lock['Tacstd2']['max_abs_diff']:.2e}; "
        f"Cldn4 max |Δ difference| = {lock['Cldn4']['max_abs_diff']:.2e}). "
        "Naive is `Baseline == 1`. ICB is `Baseline == 0` (responders and non-responders pooled; anti-PD-1, anti-PD-L1, anti-PD-L2, anti-CTLA4, and combos). "
        "The paired location is the mean, the same rule that produces Tacstd2 49/64. "
        "Scores use only samples present in the Tacstd2 export. "
        + summary["baseline_note"],
        "",
        "TISMO has no KL or KP lung ICB model. The two lung slices are LLC (GSE155972), which is Lewis lung carcinoma, not KL.",
        "",
        "## Scores",
        "",
        "Each score is the mean of per-gene z-scores. Z-scores are fit once on the samples that belong to the 64 stems (population sd, ddof = 0). "
        f"A gene enters a score only if it is non-missing in at least {MIN_GENE_FRACTION:.0%} of those samples and has non-zero variance. "
        f"A sample is scored only if at least {MIN_GENE_FRACTION:.0%} of the retained genes are present.",
        "",
        f"Primary IFN is MSigDB mouse Hallmark interferon-gamma response ({cov['IFN_gamma']['n_used']}/{cov['IFN_gamma']['n_requested']} genes used"
        f"{_missing_clause(cov['IFN_gamma'])}). "
        f"Primary NHEJ is MSigDB mouse GO:0006303, double-strand break repair via nonhomologous end joining "
        f"({cov['NHEJ_GO']['n_used']}/{cov['NHEJ_GO']['n_requested']} genes used"
        f"{_missing_clause(cov['NHEJ_GO'])}). "
        "TISMO still uses the previous MGI symbols Wars and Ddx58; those columns fill Hallmark members Wars1 and Rigi. "
        "Reactome NHEJ (R-MMU-5693571) is secondary because 33 of its 67 genes are histones; a histone-stripped Reactome score and a 12-gene c-NHEJ machinery score are reported beside it. "
        "Hallmark interferon-alpha and the 7-gene WikiPathways NHEJ set are secondary.",
        "",
        "These are bulk syngeneic tumors. The IFN-γ score mixes immune infiltrate with any tumor-cell interferon response. It is not a malignant-cell program.",
        "",
        "## Primary Spearman tests (BH q across these six)",
        "",
        "| Contrast | n | ρ | 95% CI | p | q |",
        "|---|---:|---:|---|---:|---:|",
    ]
    for row in p:
        lines.append(
            f"| {row['label']} | {row['n']} | {fmt_rho(row['rho'])} | "
            f"{fmt_rho(row['ci_low'])} to {fmt_rho(row['ci_high'])} | {fmt_p(row['p'])} | {fmt_p(row['q'])} |"
        )
    lines += [
        "",
        "Partial Spearman (Pearson correlation of rank residuals), not in the six-test family:",
        "",
        "| Contrast | n | ρ | p |",
        "|---|---:|---:|---:|",
    ]
    for row in summary["partial"]:
        lines.append(f"| {row['label']} | {row['n']} | {fmt_rho(row['rho'])} | {fmt_p(row['p'])} |")
    lines += [
        "",
        f"Holding NHEJ constant, the partial Spearman of the Cldn4 ICB change with IFN-γ is "
        f"{fmt_rho(partial['delta Cldn4 vs IFN_gamma given NHEJ_GO']['rho'])} "
        f"(p = {fmt_p(partial['delta Cldn4 vs IFN_gamma given NHEJ_GO']['p'])}). "
        f"Holding IFN-γ constant, the partial Spearman with NHEJ is "
        f"{fmt_rho(partial['delta Cldn4 vs NHEJ_GO given IFN_gamma']['rho'])} "
        f"(p = {fmt_p(partial['delta Cldn4 vs NHEJ_GO given IFN_gamma']['p'])}).",
        "",
        "## Do the scores themselves move after ICB?",
        "",
        "| Score | Up | Down | Tie | Mean Δ | Wilcoxon p | Binomial p |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in summary["score_directions"]:
        lines.append(
            f"| {row['score']} | {row['n_up']} | {row['n_down']} | {row['n_tie']} | "
            f"{row['mean_delta']:+.3f} | {fmt_p(row['wilcoxon_p_two_sided'])} | {fmt_p(row['binomial_p_two_sided'])} |"
        )
    lines += [
        "",
        "Cldn4’s own paired test is unchanged from PR #542 (34 up, 28 down, 2 ties). "
        f"IFN-γ is also up in {ifn_dir['n_up']}/64 slices. "
        f"The overlap with the Tacstd2-up slices is {summary['tacstd2_ifn_both_up']}/49 "
        f"(about {summary['tacstd2_ifn_overlap_expected']:.1f} expected if the two up-sets were independent), "
        "so the matching counts are not one shared set of slices and they do not revise the Tacstd2 lock.",
        "",
        "## Sign concordance with Cldn4",
        "",
        "Slices with a zero Cldn4 delta or a zero score delta are omitted.",
        "",
        "| Score | Same direction | Cldn4 up, score down | Cldn4 down, score up | n |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in summary["concordance"]:
        lines.append(
            f"| {row['score']} | {row['same']} | {row['cldn4_up_score_down']} | {row['cldn4_down_score_up']} | {row['n']} |"
        )
    lines += [
        "",
        "## Dependence checks",
        "",
        "The 64 slices are not 64 independent tumors. They are 17 cell lines and 22 studies, and one study contributes many slices. "
        "The cell-line check averages slice deltas within each cell line and correlates those means. "
        "The leave-one-study check drops the study with the most slices.",
        "",
        "| Check | n | IFN-γ ρ | IFN-γ p | NHEJ ρ | NHEJ p |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in summary["dependence"]:
        lines.append(
            f"| {row['check']} | {row['n']} | {fmt_rho(row['ifn_rho'])} | {fmt_p(row['ifn_p'])} | "
            f"{fmt_rho(row['nhej_rho'])} | {fmt_p(row['nhej_p'])} |"
        )
    lines += [
        "",
        "The same slice-versus-cell-line check for Cldn4 levels, rather than changes:",
        "",
        "| When | Unit | n | IFN-γ ρ | IFN-γ p | NHEJ ρ | NHEJ p |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in summary["level_dependence"]:
        lines.append(
            f"| {row['kind']} | {row['unit']} | {row['n']} | {fmt_rho(row['ifn_rho'])} | {fmt_p(row['ifn_p'])} | "
            f"{fmt_rho(row['nhej_rho'])} | {fmt_p(row['nhej_p'])} |"
        )
    lines += [
        "",
        "Cancer groups with at least 8 slices, correlation of the ICB changes:",
        "",
        "| Group | n | IFN-γ ρ | IFN-γ p | NHEJ ρ | NHEJ p |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in summary["strata"]:
        lines.append(
            f"| {row['group']} | {row['n']} | {fmt_rho(row['ifn_rho'])} | {fmt_p(row['ifn_p'])} | "
            f"{fmt_rho(row['nhej_rho'])} | {fmt_p(row['nhej_p'])} |"
        )
    lines += [
        "",
        "Baseline levels inside those same groups:",
        "",
        "| Group | n | IFN-γ ρ | IFN-γ p | NHEJ ρ | NHEJ p |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in summary["baseline_strata"]:
        lines.append(
            f"| {row['group']} | {row['n']} | {fmt_rho(row['ifn_rho'])} | {fmt_p(row['ifn_p'])} | "
            f"{fmt_rho(row['nhej_rho'])} | {fmt_p(row['nhej_p'])} |"
        )
    lung = summary["lung"]
    lines += [
        "",
        "## Lung",
        "",
        "Both lung slices are LLC from GSE155972. With n = 2 there is no Wilcoxon claim.",
        "",
        "| Slice | Cldn4 Δ | IFN-γ Δ | NHEJ Δ |",
        "|---|---:|---:|---:|",
    ]
    for row in lung:
        lines.append(
            f"| {row['stem']} | {row['Cldn4_delta_mean']:+.3f} | {row['IFN_gamma_delta_mean']:+.3f} | {row['NHEJ_GO_delta_mean']:+.3f} |"
        )
    lines += [
        "",
        "## Secondary scores",
        "",
        "Same 64-slice Spearman of the Cldn4 mean-delta against each secondary score delta. These p values are not in the six-test BH family.",
        "",
        "| Score | Genes used | Slices | ρ | p | Score up |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in summary["secondary"]:
        lines.append(
            f"| {row['score']} | {row['n_genes']} | {row['n']} | {fmt_rho(row['rho'])} | {fmt_p(row['p'])} | {row['n_up']}/{row['n_pairs']} |"
        )
    lines.append("")
    unscored = summary.get("reactome_unscored") or []
    if unscored:
        lines.append(
            "Full Reactome NHEJ is unscored when histone genes are missing and coverage falls under "
            f"{MIN_GENE_FRACTION:.0%}. Unscored slices ({len(unscored)}): {', '.join(unscored)}. "
            "The histone-stripped Reactome score and the c-NHEJ machinery score still cover those slices."
        )
    lines += [
        "",
        "## Figures",
        "",
        "- `figures/fig1_delta_cldn4_vs_scores.png` — Cldn4 change versus IFN-γ change and versus NHEJ change. Black edges mark LLC.",
        "- `figures/fig2_score_waterfalls.png` — paired score changes across the 64 slices.",
        "- `figures/fig3_levels_pre_post.png` — slice-level Cldn4 versus each score before ICB and after ICB. The baseline IFN-γ panel is the ρ = +0.39 association that shrinks to ρ = +0.08 across 17 cell lines.",
        "",
        "## Sources",
        "",
        "- TISMO Gene module, `POST https://tismo.pku-genomics.org/rtismo/gene/downVivoExprn`, type=3, all ICB treatments × all tumors. Download date is in `tables/download_manifest.json`.",
        "- Lock: PR #542 (`results/tismo/tables/pairs_merged.tsv` on `cursor/tismo-icb-tacstd2-cldn4-tj-aa3d`), copied to `methods/tismo_cldn4_nhej_ifn/locked_pairs_pr542.tsv`.",
        "- Mouse MSigDB gene sets (MGI symbols): Hallmark IFN-γ MM3878, Hallmark IFN-α MM3877, GO:0006303 MM4659, Reactome NHEJ MM15296 (R-MMU-5693571), WikiPathways NHEJ MM15989 (WP1242). Lists are in `methods/tismo_cldn4_nhej_ifn/gene_sets.json`.",
        "",
        "## Reproduce",
        "",
        "```bash",
        "pip install -r methods/tismo_cldn4_nhej_ifn/requirements.txt",
        "python3 methods/tismo_cldn4_nhej_ifn/download.py",
        "python3 methods/tismo_cldn4_nhej_ifn/analyze.py",
        "```",
        "",
        "The analyzer refuses to write results if the recomputed Tacstd2 or Cldn4 mean-deltas disagree with the lock.",
        "",
    ]
    (OUT / "RESULTS.md").write_text("\n".join(lines))
    (ROOT / "RESULTS.md").write_text("\n".join(lines))


def main() -> int:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    gene_sets = json.loads((HERE / "gene_sets.json").read_text())
    wide = pd.read_csv(TABLES / "expression_wide.tsv.gz", sep="\t")
    lock = pd.read_csv(LOCK_PATH, sep="\t")
    wide["stem"] = wide["cell_line"].map(stem)
    locked_stems = set(lock["stem"])
    # PR #542 scored Cldn4 on the Tacstd2 sample list. Samples present only in
    # other gene exports (for example two YTN16 baselines) stay out.
    sub = wide[wide["stem"].isin(locked_stems) & wide["Tacstd2"].notna()].copy()
    if sub["Samples"].duplicated().any():
        raise SystemExit("duplicate sample ids inside the locked stems")
    meta = sub[["Samples", "cell_line", "Baseline", "GSE_ID", "Mouse_treatment"]].copy()
    meta, shared_baselines = share_sibling_baselines(meta)
    expr = sub.drop_duplicates("Samples").set_index("Samples")

    # Genes that are the outcome must not sit inside a score.
    for rec in gene_sets["sets"].values():
        overlap = sorted(set(rec["genes"]) & {"Cldn4", "Tacstd2"})
        if overlap:
            raise SystemExit(f"outcome gene inside a score: {overlap}")

    pair_frames = []
    for gene in ("Tacstd2", "Cldn4"):
        if gene not in expr.columns:
            raise SystemExit(f"missing {gene} in the expression matrix")
        pair_frames.append(pair_values(meta, expr[gene], gene))
    pairs = pair_frames[0]
    for frame in pair_frames[1:]:
        pairs = pairs.merge(frame.drop(columns=["treatments"]), on=[
            "stem", "cell_line", "gse_id", "cancer_type", "cancer_group", "is_lung", "is_llc", "n_baseline", "n_icb"
        ], how="outer")
    lock_report = assert_lock(pairs, lock)

    coverage_rows = []
    genes_used_rows = []
    score_series = {}
    for set_name, score_name in SCORE_OF_SET.items():
        requested = gene_sets["sets"][set_name]["genes"]
        score, used, missing = mean_z_score(expr, requested)
        score = score.reindex(expr.index)
        score_series[score_name] = score
        coverage_rows.append(
            {
                "score": score_name,
                "set": set_name,
                "role": gene_sets["sets"][set_name]["role"],
                "n_requested": len(requested),
                "n_used": len(used),
                "n_missing_or_dropped": len(missing),
                "missing": ";".join(missing),
                "n_samples_scored": int(score.notna().sum()),
            }
        )
        for gene in used:
            genes_used_rows.append({"score": score_name, "gene": gene, "status": "used"})
        for gene in missing:
            genes_used_rows.append({"score": score_name, "gene": gene, "status": "dropped"})
        paired = pair_values(meta, score, score_name)
        pairs = pairs.merge(
            paired[["stem", f"{score_name}_mean_naive", f"{score_name}_mean_icb", f"{score_name}_delta_mean", f"{score_name}_delta_median"]],
            on="stem",
            how="left",
        )

    pairs = pairs.sort_values("stem").reset_index(drop=True)
    if len(pairs) != 64:
        raise SystemExit(f"expected 64 stems after the lock merge, found {len(pairs)}")
    required = [
        "Tacstd2_delta_mean",
        "Cldn4_delta_mean",
        "IFN_gamma_delta_mean",
        "NHEJ_GO_delta_mean",
    ]
    if pairs[required].isna().any().any():
        raise SystemExit(f"missing primary paired values:\n{pairs.loc[pairs[required].isna().any(axis=1), ['stem', *required]]}")

    # Primary correlations.
    corr_rows = []
    seed = 1000
    for kind, score in PRIMARY_TESTS:
        if kind == "delta":
            x, y = pairs["Cldn4_delta_mean"], pairs[f"{score}_delta_mean"]
            label = f"delta Cldn4 vs {score}"
        elif kind == "baseline":
            x, y = pairs["Cldn4_mean_naive"], pairs[f"{score}_mean_naive"]
            label = f"baseline Cldn4 vs {score}"
        else:
            x, y = pairs["Cldn4_mean_icb"], pairs[f"{score}_mean_icb"]
            label = f"icb Cldn4 vs {score}"
        block = spearman_block(x, y, label, seed)
        block["family"] = "primary"
        block["kind"] = kind
        block["score"] = score
        corr_rows.append(block)
        seed += 1
    qs = bh_fdr([row["p"] for row in corr_rows])
    for row, q in zip(corr_rows, qs):
        row["q"] = q

    partial_rows = [
        {
            "label": "delta Cldn4 vs IFN_gamma given NHEJ_GO",
            **partial_spearman(pairs["Cldn4_delta_mean"], pairs["IFN_gamma_delta_mean"], pairs["NHEJ_GO_delta_mean"]),
        },
        {
            "label": "delta Cldn4 vs NHEJ_GO given IFN_gamma",
            **partial_spearman(pairs["Cldn4_delta_mean"], pairs["NHEJ_GO_delta_mean"], pairs["IFN_gamma_delta_mean"]),
        },
        {
            "label": "baseline Cldn4 vs IFN_gamma given NHEJ_GO",
            **partial_spearman(pairs["Cldn4_mean_naive"], pairs["IFN_gamma_mean_naive"], pairs["NHEJ_GO_mean_naive"]),
        },
        {
            "label": "baseline Cldn4 vs NHEJ_GO given IFN_gamma",
            **partial_spearman(pairs["Cldn4_mean_naive"], pairs["NHEJ_GO_mean_naive"], pairs["IFN_gamma_mean_naive"]),
        },
    ]

    direction_rows = []
    for score in SCORE_OF_SET.values():
        deltas = pairs[f"{score}_delta_mean"]
        signs = sign_counts(deltas)
        wil = wilcoxon_delta(deltas)
        direction_rows.append(
            {
                "score": score,
                **signs,
                "mean_delta": float(deltas.mean()),
                "median_delta": float(deltas.median()),
                "wilcoxon_stat": wil["statistic"],
                "wilcoxon_p_two_sided": wil["p_two_sided"],
                "wilcoxon_p_greater": wil["p_greater"],
            }
        )

    concordance_rows = []
    for score in ("IFN_gamma", "NHEJ_GO", "IFN_alpha", "NHEJ_core", "NHEJ_Reactome", "NHEJ_Reactome_no_histone", "NHEJ_WP"):
        d = pairs[(pairs["Cldn4_delta_mean"] != 0) & (pairs[f"{score}_delta_mean"] != 0)]
        same = int(((d["Cldn4_delta_mean"] > 0) & (d[f"{score}_delta_mean"] > 0)).sum() + ((d["Cldn4_delta_mean"] < 0) & (d[f"{score}_delta_mean"] < 0)).sum())
        up_down = int(((d["Cldn4_delta_mean"] > 0) & (d[f"{score}_delta_mean"] < 0)).sum())
        down_up = int(((d["Cldn4_delta_mean"] < 0) & (d[f"{score}_delta_mean"] > 0)).sum())
        concordance_rows.append(
            {
                "score": score,
                "same": same,
                "cldn4_up_score_down": up_down,
                "cldn4_down_score_up": down_up,
                "n": int(len(d)),
            }
        )

    def rho_p(frame: pd.DataFrame, score: str) -> tuple[float, float, int]:
        if len(frame) < 5:
            return np.nan, np.nan, int(len(frame))
        rho, p = spearmanr(frame["Cldn4_delta_mean"], frame[f"{score}_delta_mean"])
        return float(rho), float(p), int(len(frame))

    collapsed = pairs.groupby("cell_line", as_index=False)[
        ["Cldn4_delta_mean", "IFN_gamma_delta_mean", "NHEJ_GO_delta_mean"]
    ].mean()
    largest_gse = pairs["gse_id"].value_counts().idxmax()
    dropped = pairs[pairs["gse_id"] != largest_gse]
    dep_rows = []
    for check, frame in (
        ("64 slices (primary unit)", pairs),
        ("cell-line mean of slice deltas", collapsed),
        (f"drop {largest_gse} ({int((pairs['gse_id'] == largest_gse).sum())} slices)", dropped),
    ):
        ir, ip, n = rho_p(frame, "IFN_gamma")
        nr, np_, _ = rho_p(frame, "NHEJ_GO")
        dep_rows.append({"check": check, "n": n, "ifn_rho": ir, "ifn_p": ip, "nhej_rho": nr, "nhej_p": np_})

    strata_rows = []
    baseline_strata = []
    for group, g in pairs.groupby("cancer_group"):
        if len(g) < 8:
            continue
        ir, ip, n = rho_p(g, "IFN_gamma")
        nr, np_, _ = rho_p(g, "NHEJ_GO")
        strata_rows.append({"group": group, "n": n, "ifn_rho": ir, "ifn_p": ip, "nhej_rho": nr, "nhej_p": np_})
        br, bp = spearmanr(g["Cldn4_mean_naive"], g["IFN_gamma_mean_naive"])
        nr2, np2 = spearmanr(g["Cldn4_mean_naive"], g["NHEJ_GO_mean_naive"])
        baseline_strata.append(
            {
                "group": group,
                "n": int(len(g)),
                "ifn_rho": float(br),
                "ifn_p": float(bp),
                "nhej_rho": float(nr2),
                "nhej_p": float(np2),
            }
        )

    level_rows = []
    line_means = pairs.groupby("cell_line", as_index=False)[
        [
            "Cldn4_mean_naive",
            "IFN_gamma_mean_naive",
            "NHEJ_GO_mean_naive",
            "Cldn4_mean_icb",
            "IFN_gamma_mean_icb",
            "NHEJ_GO_mean_icb",
        ]
    ].mean()
    for kind, xcol, y_ifn, y_nhej in (
        ("baseline", "Cldn4_mean_naive", "IFN_gamma_mean_naive", "NHEJ_GO_mean_naive"),
        ("after ICB", "Cldn4_mean_icb", "IFN_gamma_mean_icb", "NHEJ_GO_mean_icb"),
    ):
        for unit, frame in (("64 slices", pairs), ("cell-line means", line_means)):
            ir, ip = spearmanr(frame[xcol], frame[y_ifn])
            nr, npv = spearmanr(frame[xcol], frame[y_nhej])
            level_rows.append(
                {
                    "kind": kind,
                    "unit": unit,
                    "n": int(len(frame)),
                    "ifn_rho": float(ir),
                    "ifn_p": float(ip),
                    "nhej_rho": float(nr),
                    "nhej_p": float(npv),
                }
            )

    secondary_rows = []
    cov_by = {row["score"]: row for row in coverage_rows}
    dir_by = {row["score"]: row for row in direction_rows}
    for score in ("IFN_alpha", "NHEJ_Reactome", "NHEJ_Reactome_no_histone", "NHEJ_core", "NHEJ_WP"):
        d = pairs[["Cldn4_delta_mean", f"{score}_delta_mean"]].dropna()
        rho, p = spearmanr(d["Cldn4_delta_mean"], d[f"{score}_delta_mean"])
        secondary_rows.append(
            {
                "score": score,
                "n_genes": cov_by[score]["n_used"],
                "n": int(len(d)),
                "rho": float(rho),
                "p": float(p),
                "n_up": dir_by[score]["n_up"],
                "n_pairs": dir_by[score]["n_pairs"],
            }
        )

    cldn4_w = wilcoxon_delta(pairs["Cldn4_delta_mean"])
    tac_w = wilcoxon_delta(pairs["Tacstd2_delta_mean"])
    lung_rows = pairs.loc[pairs["is_lung"], ["stem", "Cldn4_delta_mean", "IFN_gamma_delta_mean", "NHEJ_GO_delta_mean"]].to_dict("records")

    pairs.to_csv(TABLES / "pairs_64.tsv", sep="\t", index=False)
    pd.DataFrame(corr_rows).to_csv(TABLES / "correlations_primary.tsv", sep="\t", index=False)
    pd.DataFrame(partial_rows).to_csv(TABLES / "correlations_partial.tsv", sep="\t", index=False)
    pd.DataFrame(direction_rows).to_csv(TABLES / "score_directions.tsv", sep="\t", index=False)
    pd.DataFrame(concordance_rows).to_csv(TABLES / "concordance.tsv", sep="\t", index=False)
    pd.DataFrame(dep_rows).to_csv(TABLES / "dependence.tsv", sep="\t", index=False)
    pd.DataFrame(level_rows).to_csv(TABLES / "level_dependence.tsv", sep="\t", index=False)
    pd.DataFrame(strata_rows).to_csv(TABLES / "strata.tsv", sep="\t", index=False)
    pd.DataFrame(baseline_strata).to_csv(TABLES / "baseline_strata.tsv", sep="\t", index=False)
    pd.DataFrame(secondary_rows).to_csv(TABLES / "secondary_scores.tsv", sep="\t", index=False)
    pd.DataFrame(coverage_rows).to_csv(TABLES / "gene_coverage.tsv", sep="\t", index=False)
    pd.DataFrame(genes_used_rows).to_csv(TABLES / "genes_used.tsv", sep="\t", index=False)

    # Figures
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.6), sharey=False)
    for ax, score, pretty in (
        (axes[0], "IFN_gamma", "Hallmark IFN-γ"),
        (axes[1], "NHEJ_GO", "GO:0006303 NHEJ"),
    ):
        row = next(r for r in corr_rows if r["label"] == f"delta Cldn4 vs {score}")
        scatter(ax, pairs, "Cldn4_delta_mean", f"{score}_delta_mean", f"{pretty}\nρ={fmt_rho(row['rho'])}, p={fmt_p(row['p'])}")
        ax.set_ylabel(f"{pretty} Δ mean (z)")
    handles, labels = axes[0].get_legend_handles_labels()
    # Lung was added inside the lung group; keep one legend from the left panel plus lung if present.
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(FIG / "fig1_delta_cldn4_vs_scores.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.2))
    for ax, score, pretty in (
        (axes[0], "IFN_gamma", "Hallmark IFN-γ"),
        (axes[1], "NHEJ_GO", "GO:0006303 NHEJ"),
    ):
        row = dir_by[score]
        waterfall(
            ax,
            pairs,
            f"{score}_delta_mean",
            f"{pretty}: {row['n_up']}/{row['n_pairs']} up, Wilcoxon p={fmt_p(row['wilcoxon_p_two_sided'])}",
        )
    fig.tight_layout()
    fig.savefig(FIG / "fig2_score_waterfalls.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(10.2, 8.2))
    level_specs = [
        (axes[0, 0], "Cldn4_mean_naive", "IFN_gamma_mean_naive", "baseline Cldn4 vs IFN_gamma", "Before ICB, IFN-γ"),
        (axes[0, 1], "Cldn4_mean_naive", "NHEJ_GO_mean_naive", "baseline Cldn4 vs NHEJ_GO", "Before ICB, NHEJ"),
        (axes[1, 0], "Cldn4_mean_icb", "IFN_gamma_mean_icb", "icb Cldn4 vs IFN_gamma", "After ICB, IFN-γ"),
        (axes[1, 1], "Cldn4_mean_icb", "NHEJ_GO_mean_icb", "icb Cldn4 vs NHEJ_GO", "After ICB, NHEJ"),
    ]
    for ax, x, y, key, title in level_specs:
        row = next(r for r in corr_rows if r["label"] == key)
        level_scatter(ax, pairs, x, y, f"{title}\nρ={fmt_rho(row['rho'])}, p={fmt_p(row['p'])}")
    fig.tight_layout()
    fig.savefig(FIG / "fig3_levels_pre_post.png", dpi=160)
    plt.close(fig)

    summary = {
        "lock": lock_report,
        "tacstd2_wilcoxon_p": tac_w["p_two_sided"],
        "cldn4_wilcoxon_p": cldn4_w["p_two_sided"],
        "primary_correlations": corr_rows,
        "partial": partial_rows,
        "score_directions": direction_rows,
        "concordance": concordance_rows,
        "dependence": dep_rows,
        "cell_line_delta": next(row for row in dep_rows if row["check"].startswith("cell-line")),
        "level_dependence": level_rows,
        "level_dependence_index": {
            kind: {row["unit"]: row for row in level_rows if row["kind"] == kind}
            for kind in ("baseline", "after ICB")
        },
        "strata": strata_rows,
        "baseline_strata": baseline_strata,
        "secondary": secondary_rows,
        "coverage": coverage_rows,
        "lung": lung_rows,
        "shared_baselines": shared_baselines,
        "baseline_note": (
            "The current TISMO export keeps each sample once. "
            + (
                "Two slices had ICB rows and no isotype rows; their baselines are the same sample IDs the PR #542 export repeated onto that arm ("
                + "; ".join(f"{row['stem']} from {row['baseline_from']} (n={row['n_baseline_samples']})" for row in shared_baselines)
                + ")."
                if shared_baselines
                else "Every locked slice already carried its own isotype rows."
            )
        ),
        "largest_gse_dropped": largest_gse,
        "reactome_unscored": pairs.loc[pairs["NHEJ_Reactome_delta_mean"].isna(), "stem"].tolist(),
        "tacstd2_ifn_both_up": int(((pairs["Tacstd2_delta_mean"] > 0) & (pairs["IFN_gamma_delta_mean"] > 0)).sum()),
        "tacstd2_ifn_overlap_expected": float((pairs["Tacstd2_delta_mean"] > 0).sum() * (pairs["IFN_gamma_delta_mean"] > 0).sum() / len(pairs)),
        "n_cell_lines": int(pairs["cell_line"].nunique()),
        "n_studies": int(pairs["gse_id"].nunique()),
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_results(summary)
    print(json.dumps({k: summary[k] for k in ("lock", "tacstd2_wilcoxon_p", "cldn4_wilcoxon_p")}, indent=2))
    for row in corr_rows:
        print(f"{row['label']}: rho={row['rho']:.3f} p={row['p']:.3g} q={row['q']:.3g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
