#!/usr/bin/env python3
"""TISMO Tacstd2/Cldn4 vs ICI response and immune infiltrate.

Maximize odds ratio, Cliff's delta, and |Spearman ρ| on two honest units:

- mouse: one unique SRA run (SRX)
- study: one GEO series (GSE), summarized by the median of its baseline mice

Arm (TISMO cell-line × series × treatment group) is the unit that actually
carries the responder label. It is reported separately and is not called a
study. Within-line residual correlations are a different estimand and are
reported separately.

Nothing is imputed. Median splits that collapse (zero-inflated scores whose
median is the floor) are ineligible for OR and Cliff's delta. Expression is
the TISMO `value` column, untransformed.
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "tismo_max_effect"
OUT = ROOT / "results" / "tismo_max_effect"
TABLES = OUT / "tables"
FIGS = OUT / "figures"

META_BASE = "https://tismo.pku-genomics.org/tismo"
R_BASE = "https://tismo.pku-genomics.org/rtismo"
ICB_TREATMENTS = [
    "antiCTLA4",
    "antiCTLA4&antiPD1",
    "antiCTLA4&antiPDL1",
    "antiPD1",
    "antiPDL1",
    "antiPDL2",
]

GENES = ("Tacstd2", "Cldn4")
RNG = np.random.default_rng(0)
N_BOOT = 2000

# Non-immune rows in the TISMO infiltrate export.
STROMAL_TOKENS = (
    "Fibroblast",
    "Endothelial",
    "Cancer associated",
    "Vessels",
    "Lymphatics",
    "progenitor",
    "Hematopoietic stem",
    "Common lymphoid",
    "Common myeloid",
)

CANCER = {
    "4T1": "Mammary",
    "E0771": "Mammary",
    "EMT6": "Mammary",
    "KPB25L": "Mammary",
    "T11": "Mammary",
    "p53-2225L": "Mammary",
    "p53-2336R": "Mammary",
    "B16": "Melanoma",
    "YUMM1.7": "Melanoma",
    "D3UV2": "Melanoma",
    "D4M.3A.3": "Melanoma",
    "CT26": "Colorectal",
    "MC38": "Colorectal",
    "LLC": "Lung",
    "402230": "Sarcoma",
    "BNL-MEA": "Liver",
    "YTN16": "Gastric",
    "MOC22": "HeadNeck",
}

CD8_FEATURES = (
    "T CD8_TIMER",
    "T CD8_CIBERSORT_abs",
    "T CD8_EPIC",
    "T CD8_quanTIseq",
    "T CD8_xCell",
    "CD8 T_mMCPcounter",
)

MIN_N = {
    "mouse": 30,
    "study": 10,
    "arm": 12,
    "line": 8,
    "mouse_line_residual": 30,
}


def stem_label(value: str) -> str:
    return re.sub(r"\(n=\d+\)$", "", str(value)).rstrip()


def is_immune(name: str) -> bool:
    return not any(token in name for token in STROMAL_TOKENS)


def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    """Cliff's δ. Positive means values in x tend to exceed values in y."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[~np.isnan(x)]
    y = y[~np.isnan(y)]
    if len(x) == 0 or len(y) == 0:
        return float("nan")
    u = stats.mannwhitneyu(x, y, alternative="two-sided").statistic
    return float(2 * u / (len(x) * len(y)) - 1)


def spearman(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = ~(np.isnan(x) | np.isnan(y))
    x, y = x[mask], y[mask]
    if len(x) < 3 or np.unique(x).size < 2 or np.unique(y).size < 2:
        return float("nan"), float("nan"), int(len(x))
    rho, p = stats.spearmanr(x, y)
    return float(rho), float(p), int(len(x))


def median_split_or(x, y):
    """OR for (y >= median) given (x >= median). Ineligible when a margin collapses."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = ~(np.isnan(x) | np.isnan(y))
    x, y = x[mask], y[mask]
    n = int(len(x))
    empty = {
        "n": n,
        "a": 0,
        "b": 0,
        "c": 0,
        "d": 0,
        "n_high_x": 0,
        "n_low_x": 0,
        "n_high_y": 0,
        "n_low_y": 0,
        "OR": float("nan"),
        "OR_low": float("nan"),
        "OR_high": float("nan"),
        "fisher_p": float("nan"),
        "eligible": False,
        "reason": "too few",
    }
    if n < 8:
        return empty
    mx = float(np.median(x))
    my = float(np.median(y))
    hx = x >= mx
    hy = y >= my
    a = int((hx & hy).sum())
    b = int((hx & ~hy).sum())
    c = int((~hx & hy).sum())
    d = int((~hx & ~hy).sum())
    n_hx, n_lx = int(hx.sum()), int((~hx).sum())
    n_hy, n_ly = int(hy.sum()), int((~hy).sum())
    out = {
        "n": n,
        "a": a,
        "b": b,
        "c": c,
        "d": d,
        "median_x": mx,
        "median_y": my,
        "n_high_x": n_hx,
        "n_low_x": n_lx,
        "n_high_y": n_hy,
        "n_low_y": n_ly,
    }
    # Both sides of each split must be real groups. A median of 0 on a
    # zero-inflated score puts every observation in the high group.
    margins_ok = min(n_hx, n_lx, n_hy, n_ly) >= 3 and min(n_hx, n_lx) >= 0.2 * n and min(n_hy, n_ly) >= 0.2 * n
    if not margins_ok or min(a, b, c, d) < 1:
        out.update(OR=float("nan"), OR_low=float("nan"), OR_high=float("nan"), fisher_p=float("nan"), eligible=False, reason="split collapsed or empty cell")
        return out
    odd = (a / b) / (c / d)
    log_or = math.log(odd)
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    _, fisher_p = stats.fisher_exact([[a, b], [c, d]])
    out.update(
        OR=float(odd),
        OR_low=float(math.exp(log_or - 1.96 * se)),
        OR_high=float(math.exp(log_or + 1.96 * se)),
        fisher_p=float(fisher_p),
        eligible=True,
        reason="ok",
    )
    return out


def bootstrap_spearman(x, y, n_boot=N_BOOT):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = ~(np.isnan(x) | np.isnan(y))
    x, y = x[mask], y[mask]
    n = len(x)
    if n < 8:
        return float("nan"), float("nan")
    rhos = []
    for _ in range(n_boot):
        idx = RNG.integers(0, n, n)
        if np.unique(x[idx]).size < 2 or np.unique(y[idx]).size < 2:
            continue
        rhos.append(stats.spearmanr(x[idx], y[idx]).statistic)
    if len(rhos) < 100:
        return float("nan"), float("nan")
    lo, hi = np.percentile(rhos, [2.5, 97.5])
    return float(lo), float(hi)


def bootstrap_cliff(high, low, n_boot=N_BOOT):
    high = np.asarray(high, dtype=float)
    low = np.asarray(low, dtype=float)
    high = high[~np.isnan(high)]
    low = low[~np.isnan(low)]
    if len(high) < 5 or len(low) < 5:
        return float("nan"), float("nan")
    deltas = []
    for _ in range(n_boot):
        h = high[RNG.integers(0, len(high), len(high))]
        l = low[RNG.integers(0, len(low), len(low))]
        deltas.append(cliffs_delta(h, l))
    lo, hi = np.percentile(deltas, [2.5, 97.5])
    return float(lo), float(hi)


def leave_one_out_spearman(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = ~(np.isnan(x) | np.isnan(y))
    x, y = x[mask], y[mask]
    rhos = []
    for i in range(len(x)):
        rho, _, n = spearman(np.delete(x, i), np.delete(y, i))
        if n >= 8 and rho == rho:
            rhos.append(rho)
    if not rhos:
        return float("nan"), float("nan")
    return float(min(rhos)), float(max(rhos))


def prepare_gene(path: Path, gene: str, immune: pd.DataFrame) -> pd.DataFrame:
    df = pd.read_csv(path)
    if df["value"].dtype == object:
        df["value"] = pd.to_numeric(df["value"], errors="coerce")
    # Duplicate SRX rows are the same value copied onto two arms that share controls.
    value_nunique = df.groupby("Samples")["value"].nunique(dropna=True)
    if int(value_nunique.max()) > 1:
        raise RuntimeError(f"{gene}: a sample has conflicting expression values")
    df = df.drop_duplicates("Samples").dropna(subset=["value"]).copy()
    df["stem"] = df["cell_line"].map(stem_label)
    df["line"] = df["stem"].str.split("_").str[0]
    df["cancer"] = df["line"].map(CANCER).fillna("Other")
    df["gene"] = gene
    keep = ["Samples", "value", "stem", "line", "cancer", "Responder", "Baseline", "GSE_ID", "Mouse_treatment", "gene"]
    df = df[keep].merge(immune, left_on="Samples", right_index=True, how="left")
    return df


def load_immune(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path)
    raw["value"] = pd.to_numeric(raw["value"], errors="coerce")
    # Identical duplicate sample × feature rows are collapsed. Conflicting values abort.
    nunique = raw.groupby(["Samples", "geneID"])["value"].nunique(dropna=True)
    if int((nunique > 1).sum()) > 0:
        raise RuntimeError("infiltrate table has conflicting values for a sample × feature")
    raw = raw.groupby(["Samples", "geneID"], as_index=False)["value"].mean()
    wide = raw.pivot(index="Samples", columns="geneID", values="value")
    cols = [c for c in wide.columns if is_immune(c)]
    return wide[cols]


def aggregate(df: pd.DataFrame, unit: str, features: list[str]) -> pd.DataFrame:
    base = df[df["Baseline"] == 1].copy()
    if unit == "mouse":
        return base
    if unit == "study":
        keys = ["GSE_ID"]
    elif unit == "line":
        keys = ["line"]
    elif unit == "arm":
        keys = ["stem"]
    else:
        raise ValueError(unit)
    rows = []
    for key, g in base.groupby(keys):
        rec = {
            "unit_id": key if isinstance(key, str) else key[0],
            "value": float(g["value"].median()),
            "n_mice": int(len(g)),
            "line": g["line"].iloc[0] if g["line"].nunique() == 1 else "mixed",
            "cancer": g["cancer"].iloc[0] if g["cancer"].nunique() == 1 else "mixed",
            "GSE_ID": g["GSE_ID"].iloc[0] if g["GSE_ID"].nunique() == 1 else "mixed",
        }
        for feat in features:
            rec[feat] = float(g[feat].median()) if g[feat].notna().any() else float("nan")
        rows.append(rec)
    return pd.DataFrame(rows)


def score_feature(x, y, gene, feature, unit, stratum, min_n) -> dict | None:
    rho, p, n = spearman(x, y)
    if n < min_n or rho != rho:
        return None
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = ~(np.isnan(x) | np.isnan(y))
    x, y = x[mask], y[mask]
    n_unique_x = int(np.unique(np.round(x, 10)).size)
    n_unique_y = int(np.unique(np.round(y, 10)).size)
    # A Spearman of 1.0 on a handful of distinct values is a tie pattern
    # (several studies at 0, 0), not a maximized effect.
    rho_eligible = n_unique_x >= 5 and n_unique_y >= 5
    mx = float(np.median(x))
    high = y[x >= mx]
    low = y[x < mx]
    min_arm = 5 if unit != "mouse" else 15
    cliff_ok = len(high) >= min_arm and len(low) >= min_arm
    delta = cliffs_delta(high, low) if cliff_ok else float("nan")
    odd = median_split_or(x, y)
    rec = {
        "gene": gene,
        "feature": feature,
        "unit": unit,
        "stratum": stratum,
        "cohort": "baseline",
        "n": n,
        "rho": rho,
        "abs_rho": abs(rho) if rho_eligible else float("nan"),
        "rho_raw": rho,
        "rho_eligible": rho_eligible,
        "n_unique_x": n_unique_x,
        "n_unique_y": n_unique_y,
        "spearman_p": p,
        "cliff": delta,
        "abs_cliff": abs(delta) if delta == delta else float("nan"),
        "cliff_eligible": bool(cliff_ok),
        "n_high_gene": int(len(high)),
        "n_low_gene": int(len(low)),
        "OR": odd["OR"],
        "OR_low": odd["OR_low"],
        "OR_high": odd["OR_high"],
        "fisher_p": odd["fisher_p"],
        "OR_eligible": bool(odd["eligible"]),
        "a_high_high": odd["a"],
        "b_high_low": odd["b"],
        "c_low_high": odd["c"],
        "d_low_low": odd["d"],
        "OR_reason": odd["reason"],
    }
    if odd["eligible"]:
        rec["OR_extreme"] = max(odd["OR"], 1 / odd["OR"])
    else:
        rec["OR_extreme"] = float("nan")
    return rec


def infiltrate_grid(genes: dict[str, pd.DataFrame], features: list[str]) -> pd.DataFrame:
    rows = []
    for gene, df in genes.items():
        tables = {
            "mouse": aggregate(df, "mouse", features),
            "study": aggregate(df, "study", features),
            "line": aggregate(df, "line", features),
        }
        for unit, tab in tables.items():
            strata = [("pan-cancer", tab)]
            if unit in ("mouse", "study") and "cancer" in tab.columns:
                for cancer, sub in tab.groupby("cancer"):
                    if cancer in ("mixed", "Other"):
                        continue
                    strata.append((str(cancer), sub))
            for stratum, sub in strata:
                min_n = MIN_N[unit]
                if stratum != "pan-cancer":
                    min_n = max(6, MIN_N[unit] // 2) if unit == "study" else MIN_N[unit]
                for feat in features:
                    rec = score_feature(sub["value"], sub[feat], gene, feat, unit, stratum, min_n)
                    if rec is not None:
                        rows.append(rec)
        # Within-line residual: mouse association after subtracting the line median.
        base = df[df["Baseline"] == 1].copy()
        counts = base["line"].value_counts()
        base = base[base["line"].isin(counts[counts >= 3].index)]
        for feat in features:
            x = base["value"] - base.groupby("line")["value"].transform("median")
            y = base[feat] - base.groupby("line")[feat].transform("median")
            rec = score_feature(x, y, gene, feat, "mouse_line_residual", "pan-cancer", MIN_N["mouse_line_residual"])
            if rec is not None:
                rec["n_lines"] = int(base.loc[base["value"].notna() & base[feat].notna(), "line"].nunique())
                rows.append(rec)
    return pd.DataFrame(rows)


def arm_table(df: pd.DataFrame) -> pd.DataFrame:
    treated = df[df["Responder"].isin(["Responders", "Non-responders"])].copy()
    rows = []
    for stem, g in treated.groupby("stem"):
        n_r = int((g["Responder"] == "Responders").sum())
        n_nr = int((g["Responder"] == "Non-responders").sum())
        if n_r > 0 and n_nr > 0:
            label = "mixed"
            med = float("nan")
        elif n_r > 0:
            label = "R"
            med = float(g["value"].median())
        else:
            label = "NR"
            med = float(g["value"].median())
        rows.append(
            {
                "stem": stem,
                "label": label,
                "median": med,
                "n_mice": int(len(g)),
                "n_R": n_r,
                "n_NR": n_nr,
                "line": g["line"].iloc[0],
                "cancer": g["cancer"].iloc[0],
                "GSE_ID": g["GSE_ID"].iloc[0],
                "treatment": g["Mouse_treatment"].iloc[0],
            }
        )
    return pd.DataFrame(rows)


def response_grid(genes: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, dict]:
    rows = []
    mixed_detail = {}
    for gene, df in genes.items():
        arms = arm_table(df)
        pure = arms[arms["label"].isin(["R", "NR"])].copy()
        strata = [("pan-cancer", pure)]
        for cancer, sub in pure.groupby("cancer"):
            strata.append((str(cancer), sub))
        for stratum, sub in strata:
            r = sub.loc[sub["label"] == "R", "median"]
            nr = sub.loc[sub["label"] == "NR", "median"]
            if len(sub) < 6 or len(r) < 4 or len(nr) < 4:
                eligible_groups = False
            else:
                eligible_groups = True
            rho, p, n = spearman(sub["median"], (sub["label"] == "R").astype(float))
            delta = cliffs_delta(r.to_numpy(), nr.to_numpy()) if len(r) and len(nr) else float("nan")
            odd = median_split_or(sub["median"], (sub["label"] == "R").astype(float))
            # Response is already binary. A median of 1 (more responders than not)
            # collapses the outcome split. Require a real gene split and both outcomes.
            if odd["n_low_x"] < 3 or odd["n_high_x"] < 3 or min(odd["n_high_y"], odd["n_low_y"]) < 1:
                odd["eligible"] = False
            rows.append(
                {
                    "gene": gene,
                    "unit": "arm",
                    "stratum": stratum,
                    "n_arms": int(len(sub)),
                    "n_R": int(len(r)),
                    "n_NR": int(len(nr)),
                    "rho": rho,
                    "abs_rho": abs(rho) if rho == rho else float("nan"),
                    "spearman_p": p,
                    "cliff_R_vs_NR": delta if eligible_groups else float("nan"),
                    "abs_cliff": abs(delta) if eligible_groups and delta == delta else float("nan"),
                    "cliff_eligible": bool(eligible_groups),
                    "OR_response_given_high_gene": odd["OR"] if odd["eligible"] else float("nan"),
                    "OR_low": odd["OR_low"] if odd["eligible"] else float("nan"),
                    "OR_high": odd["OR_high"] if odd["eligible"] else float("nan"),
                    "fisher_p": odd["fisher_p"] if odd["eligible"] else float("nan"),
                    "OR_eligible": bool(odd["eligible"]),
                    "a_high_R": odd["a"],
                    "b_high_NR": odd["b"],
                    "c_low_R": odd["c"],
                    "d_low_NR": odd["d"],
                    "median_R": float(r.median()) if len(r) else float("nan"),
                    "median_NR": float(nr.median()) if len(nr) else float("nan"),
                }
            )
        # Study-level: median gene vs responder fraction among treated mice.
        treated = df[df["Responder"].isin(["Responders", "Non-responders"])]
        study_rows = []
        for gse, g in treated.groupby("GSE_ID"):
            study_rows.append(
                {
                    "GSE_ID": gse,
                    "median": float(g["value"].median()),
                    "response_rate": float((g["Responder"] == "Responders").mean()),
                    "n_mice": int(len(g)),
                    "n_lines": int(g["line"].nunique()),
                }
            )
        st = pd.DataFrame(study_rows)
        rho, p, n = spearman(st["median"], st["response_rate"])
        odd = median_split_or(st["median"], st["response_rate"])
        rows.append(
            {
                "gene": gene,
                "unit": "study",
                "stratum": "pan-cancer",
                "n_arms": int(n),
                "n_R": int((st["response_rate"] >= st["response_rate"].median()).sum()) if len(st) else 0,
                "n_NR": int((st["response_rate"] < st["response_rate"].median()).sum()) if len(st) else 0,
                "rho": rho,
                "abs_rho": abs(rho) if rho == rho else float("nan"),
                "spearman_p": p,
                "cliff_R_vs_NR": float("nan"),
                "abs_cliff": float("nan"),
                "cliff_eligible": False,
                "OR_response_given_high_gene": odd["OR"] if odd["eligible"] else float("nan"),
                "OR_low": odd["OR_low"] if odd["eligible"] else float("nan"),
                "OR_high": odd["OR_high"] if odd["eligible"] else float("nan"),
                "fisher_p": odd["fisher_p"] if odd["eligible"] else float("nan"),
                "OR_eligible": bool(odd["eligible"]),
                "a_high_R": odd["a"],
                "b_high_NR": odd["b"],
                "c_low_R": odd["c"],
                "d_low_NR": odd["d"],
                "median_R": float("nan"),
                "median_NR": float("nan"),
            }
        )
        # The only arms with mouse-level responder labels.
        mixed = arms[arms["label"] == "mixed"]
        detail = []
        for _, arm in mixed.iterrows():
            g = treated[treated["stem"] == arm["stem"]]
            r = g.loc[g["Responder"] == "Responders", "value"]
            nr = g.loc[g["Responder"] == "Non-responders", "value"]
            detail.append(
                {
                    "stem": arm["stem"],
                    "GSE_ID": arm["GSE_ID"],
                    "line": arm["line"],
                    "n_R": int(len(r)),
                    "n_NR": int(len(nr)),
                    "median_R": float(r.median()),
                    "median_NR": float(nr.median()),
                    "cliff_R_vs_NR": cliffs_delta(r.to_numpy(), nr.to_numpy()),
                }
            )
        mixed_detail[gene] = detail
    return pd.DataFrame(rows), mixed_detail


def pick_max(df: pd.DataFrame, column: str, mask: pd.Series) -> dict | None:
    sub = df.loc[mask & df[column].notna()]
    if sub.empty:
        return None
    return sub.loc[sub[column].idxmax()].to_dict()


def fmt_p(p) -> str:
    if p is None or (isinstance(p, float) and (math.isnan(p) or math.isinf(p))):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3g}"


def fmt(x, digits=3) -> str:
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return "NA"
    return f"{x:.{digits}f}"


def annotate(ax, xs, ys, labels):
    for x, y, lab in zip(xs, ys, labels):
        ax.annotate(lab, (x, y), textcoords="offset points", xytext=(3, 3), fontsize=7, color="#333333")


def plot_study_scatter(tab: pd.DataFrame, feature: str, gene: str, path: Path, rho: float, n: int):
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    colors = tab["cancer"].map(
        {
            "Mammary": "#b45309",
            "Melanoma": "#1d4ed8",
            "Colorectal": "#15803d",
            "Lung": "#b91c1c",
            "Sarcoma": "#7c3aed",
            "Liver": "#0f766e",
            "Gastric": "#be185d",
            "HeadNeck": "#44403c",
            "mixed": "#737373",
        }
    ).fillna("#737373")
    ax.scatter(tab["value"], tab[feature], s=28 + 4 * tab["n_mice"], c=colors, zorder=3)
    annotate(ax, tab["value"], tab[feature], tab["unit_id"].str.replace("GSE", ""))
    ax.set_xlabel(f"Baseline {gene} (TISMO value, study median)")
    ax.set_ylabel(feature.replace("_", " "))
    ax.set_title(f"{gene} vs {feature}\nstudy unit, Spearman ρ = {rho:.3f}, n = {n} studies")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_quadrants(tab: pd.DataFrame, feature: str, gene: str, path: Path, cells: dict):
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    mx = float(tab["value"].median())
    my = float(tab[feature].median())
    ax.axvline(mx, color="#a3a3a3", lw=1)
    ax.axhline(my, color="#a3a3a3", lw=1)
    ax.scatter(tab["value"], tab[feature], s=36, c="#0f172a", zorder=3)
    annotate(ax, tab["value"], tab[feature], tab["unit_id"].str.replace("GSE", ""))
    ax.set_xlabel(f"Baseline {gene} (study median)")
    ax.set_ylabel(feature.replace("_", " "))
    ax.set_title(
        f"{gene} high/low vs {feature} high/low\n"
        f"OR = {cells['OR']:.0f}  table HH/HL/LH/LL = "
        f"{cells['a']}/{cells['b']}/{cells['c']}/{cells['d']}  n = {len(tab)} studies"
    )
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_top_bars(specs: pd.DataFrame, path: Path):
    sub = specs[(specs["unit"] == "study") & (specs["stratum"] == "pan-cancer")].nlargest(12, "abs_rho")
    sub = sub.iloc[::-1]
    fig, ax = plt.subplots(figsize=(8, 5.2))
    colors = ["#b91c1c" if g == "Cldn4" else "#1d4ed8" for g in sub["gene"]]
    ax.barh(np.arange(len(sub)), sub["rho"], color=colors)
    labels = [f"{g} · {f}" for g, f in zip(sub["gene"], sub["feature"])]
    ax.set_yticks(np.arange(len(sub)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.axvline(0, color="#111", lw=0.8)
    ax.set_xlabel("Spearman ρ at the study unit (baseline medians)")
    ax.set_title("Largest pan-cancer study-level |ρ|")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def write_report(summary: dict, path: Path):
    h = summary["headline"]
    rho = h["max_abs_rho_study"]
    cliff = h["max_abs_cliff_study"]
    odd = h["max_OR_study"]
    rho_m = h["max_abs_rho_mouse"]
    odd_m = h["max_OR_mouse"]
    cliff_m = h["max_abs_cliff_mouse"]
    resp = h["response_pancancer"]
    resp_max = h["response_max_abs_cliff"]
    resid = h["max_abs_rho_within_line"]
    lines = []
    lines.append("# TISMO Tacstd2/Cldn4 vs ICI response and immune infiltrate")
    lines.append("")
    lines.append("## Paper sentence")
    lines.append("")
    sens = summary["sensitivity_n_mice_ge_8"]
    lines.append(
        f"Across TISMO baseline tumors, the largest pan-cancer study-level association "
        f"with at least five distinct values on each axis is "
        f"**{rho['gene']} vs {rho['feature']}**, Spearman ρ = **{rho['rho']:.3f}** "
        f"(95% bootstrap CI {fmt(rho['rho_low'])} to {fmt(rho['rho_high'])}; "
        f"n = **{rho['n']} studies**; leave-one-study-out ρ {fmt(rho['loo_min'])} to {fmt(rho['loo_max'])}). "
        f"The largest study-level odds ratio is **{odd['gene']} vs {odd['feature']}**, "
        f"OR = **{odd['OR']:.0f}** (Woolf 95% CI {fmt(odd['OR_low'], 1)} to {fmt(odd['OR_high'], 0)}; "
        f"Fisher p = {fmt_p(odd['fisher_p'])}; "
        f"high/high, high/low, low/high, low/low = "
        f"{odd['a_high_high']}/{odd['b_high_low']}/{odd['c_low_high']}/{odd['d_low_low']}; "
        f"n = {odd['n']} studies). "
        f"Cliff's δ for high vs low {cliff['gene']} on {cliff['feature']} is **{cliff['cliff']:.3f}** "
        f"(n_high = {cliff['n_high_gene']}, n_low = {cliff['n_low_gene']}). "
        f"All three point the same way: higher bulk Tacstd2/Cldn4 tracks higher granulocyte infiltrate. "
        f"Restricting to studies with at least 8 baseline mice leaves ρ = {sens['rho_winner']['rho']:.3f} "
        f"(n = {sens['rho_winner']['n']}) and the Tacstd2 granulocyte OR at {sens['or_winner']['OR']:.0f} "
        f"({sens['or_winner']['table']}, n = {sens['or_winner']['n']})."
    )
    lines.append("")
    lines.append(
        f"小鼠单位上的最大 |ρ| 是 {rho_m['gene']} vs {rho_m['feature']}，ρ = {rho_m['rho']:.3f}（n = {rho_m['n']} 只小鼠）。"
        f"这个小鼠相关嵌在细胞系差异里。扣掉细胞系中位数之后，最大 |ρ| 是 {resid['abs_rho']:.3f}"
        f"（{resid['gene']} vs {resid['feature']}，ρ = {resid['rho']:.3f}，n = {resid['n']} 只小鼠，{int(resid['n_lines'])} 个细胞系）。"
        f"ICI 应答在全癌种 arm 单位上接近空：Tacstd2 OR = {fmt(resp['Tacstd2']['OR'])}，"
        f"Cldn4 OR = {fmt(resp['Cldn4']['OR'])}。"
    )
    lines.append("")
    lines.append("## Units")
    lines.append("")
    lines.append("Two units were eligible for the maximum.")
    lines.append("")
    lines.append("- **Mouse.** One unique SRX run. Baseline mice that TISMO copies onto two arms sharing an isotype control are kept once.")
    lines.append("- **Study.** One source series (GEO `GSE` or ArrayExpress `ERP`/`RU` accession). The gene and the infiltrate score are the median of that series' baseline mice. A 3-mouse series counts the same as a 91-mouse series. The sensitivity below repeats the winners after dropping series with fewer than 8 baseline mice.")
    lines.append("")
    lines.append("Responder labels are an arm property. In the Tacstd2 export, 62 of 64 treated arms are uniformly responders or uniformly non-responders. Only CT26 anti-PD-1 (GSE139475, 5 vs 4 mice) and YUMM1.7 BRAF-inhibitor + anti-PD-L2 (GSE103725, 2 vs 2 mice) carry mouse-level labels. A mouse-level response test would repeat the arm label. It is not used as a maximized effect.")
    lines.append("")
    lines.append("Cell-line medians and within-line residuals are reported so the study-level correlation is not mistaken for a within-mouse effect. TISMO has no KL or KP lung line. LLC is the only lung carcinoma in this ICB table.")
    lines.append("")
    lines.append("## Immune infiltrate, study unit")
    lines.append("")
    lines.append(
        f"Search space: 2 genes × {summary['n_immune_features']} immune scores at the study unit. "
        f"{summary['n_specs_study_pancancer']} rows had enough non-missing studies to score, and "
        f"{summary['n_rho_eligible_study']} of those had at least five distinct values on each axis. "
        f"Stromal, endothelial, fibroblast, and progenitor scores were excluded before ranking."
    )
    lines.append("")
    lines.append("| Metric | Winner | Estimate | n | Direction |")
    lines.append("|---|---|---|---|---|")
    lines.append(f"| max \\|ρ\\| | {rho['gene']} vs {rho['feature']} | ρ = {rho['rho']:.3f} (p = {fmt_p(rho['spearman_p'])}) | {rho['n']} studies | higher gene, higher infiltrate |")
    lines.append(f"| max \\|Cliff δ\\| | {cliff['gene']} vs {cliff['feature']} | δ = {cliff['cliff']:.3f} (CI {fmt(cliff['cliff_low'])} to {fmt(cliff['cliff_high'])}) | {cliff['n_high_gene']} vs {cliff['n_low_gene']} studies | high-gene studies have higher infiltrate |")
    lines.append(f"| max OR | {odd['gene']} vs {odd['feature']} | OR = {odd['OR']:.0f} ({odd['a_high_high']}/{odd['b_high_low']}/{odd['c_low_high']}/{odd['d_low_low']}) | {odd['n']} studies | high gene co-occurs with high infiltrate |")
    lines.append("")
    lines.append("Odds ratios use one pre-specified cut: the median of each variable inside that spec. A split is ineligible when either side has fewer than 3 observations, when either side is under 20% of n, or when any cell of the 2×2 is 0. Zero-inflated scores such as mMCPcounter neutrophils often have a median of 0, so their OR is left undefined even when Spearman ρ is large. No outcome-guided cutpoint was scanned.")
    lines.append("")
    lines.append("## Immune infiltrate, mouse unit")
    lines.append("")
    lines.append("| Metric | Winner | Estimate | n |")
    lines.append("|---|---|---|---|---|")
    lines.append(f"| max \\|ρ\\| | {rho_m['gene']} vs {rho_m['feature']} | ρ = {rho_m['rho']:.3f} (p = {fmt_p(rho_m['spearman_p'])}) | {rho_m['n']} mice |")
    lines.append(f"| max \\|Cliff δ\\| | {cliff_m['gene']} vs {cliff_m['feature']} | δ = {cliff_m['cliff']:.3f} | {cliff_m['n_high_gene']} vs {cliff_m['n_low_gene']} mice |")
    lines.append(f"| max OR | {odd_m['gene']} vs {odd_m['feature']} | OR = {odd_m['OR']:.2f} ({odd_m['a_high_high']}/{odd_m['b_high_low']}/{odd_m['c_low_high']}/{odd_m['d_low_low']}; Fisher p = {fmt_p(odd_m['fisher_p'])}) | {odd_m['n']} mice |")
    lines.append("")
    lines.append(
        f"Mouse-level p-values treat SRX runs as independent. They are descriptive. "
        f"After subtracting each cell line's median, the largest |ρ| is {resid['abs_rho']:.3f} "
        f"({resid['gene']} vs {resid['feature']}, ρ = {resid['rho']:.3f}, n = {resid['n']} mice across "
        f"{int(resid['n_lines'])} lines). The large infiltrate association is a between-study fact."
    )
    lines.append("")
    lines.append("## CD8 reference")
    lines.append("")
    cd8_best = max(summary["cd8_reference"], key=lambda row: abs(row["study_rho"]) if row["study_rho"] == row["study_rho"] else -1)
    lines.append(
        f"CD8 scores stayed in the search. At the study unit the strongest is "
        f"{cd8_best['gene']} vs {cd8_best['feature']}, ρ = {cd8_best['study_rho']:+.3f} "
        f"(n = {cd8_best['study_n']}). It is positive and smaller than the granulocyte ρ. "
        f"This table does not support CD8 exclusion."
    )
    lines.append("")
    lines.append("| Gene | Score | Study ρ | n studies | Mouse ρ | n mice |")
    lines.append("|---|---|---|---|---|---|")
    for row in summary["cd8_reference"]:
        lines.append(
            f"| {row['gene']} | {row['feature']} | {fmt(row['study_rho'])} | {row['study_n']} | {fmt(row['mouse_rho'])} | {row['mouse_n']} |"
        )
    lines.append("")
    lines.append("## ICI response")
    lines.append("")
    lines.append("Pure arms only. An arm is one TISMO comparison group. Mixed arms are held out of this table and listed below.")
    lines.append("")
    lines.append("| Gene | Unit | Stratum | n | ρ (gene vs response) | Cliff δ (responders vs non-responders) | OR for response given high gene | 2×2 R/NR among high, R/NR among low |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for row in summary["response_rows"]:
        lines.append(
            f"| {row['gene']} | {row['unit']} | {row['stratum']} | {row['n_arms']} | {fmt(row['rho'])} | "
            f"{fmt(row['cliff_R_vs_NR'])} | {fmt(row['OR_response_given_high_gene'])} | "
            f"{row['a_high_R']}/{row['b_high_NR']}/{row['c_low_R']}/{row['d_low_NR']} |"
        )
    lines.append("")
    if resp_max is not None:
        lines.append(
            f"The largest eligible response |Cliff's δ| is {resp_max['gene']} in {resp_max['stratum']} arms: "
            f"δ = {resp_max['cliff_R_vs_NR']:.3f} (responders minus non-responders), "
            f"{resp_max['n_R']} responder arms vs {resp_max['n_NR']} non-responder arms, "
            f"Spearman ρ = {fmt(resp_max['rho'])}. "
            f"Positive δ means responder arms have the higher gene median. "
            f"That is the opposite of a high-gene resistance marker. "
            f"Pan-cancer arm ORs sit near 1."
        )
        lines.append("")
    lines.append("Mouse-level labels inside the only two mixed arms:")
    lines.append("")
    lines.append("| Gene | Arm | Mice R vs NR | Median R | Median NR | Cliff δ |")
    lines.append("|---|---|---|---|---|---|")
    for row in summary["mixed_arms"]:
        lines.append(
            f"| {row['gene']} | {row['stem']} | {row['n_R']} vs {row['n_NR']} | {row['median_R']:.3f} | {row['median_NR']:.3f} | {row['cliff_R_vs_NR']:.3f} |"
        )
    lines.append("")
    lines.append("## What was not done")
    lines.append("")
    lines.append("- The locked Tacstd2 49/64 ICB-versus-control count was not recomputed and is not this result.")
    lines.append("- LLC is not labeled KL. No KP lung line is in this table.")
    lines.append("- Cutpoints were not tuned to maximize OR.")
    lines.append("- Mouse-level response tests that repeat an arm label were not treated as independent mice.")
    lines.append("")
    lines.append("## Methods")
    lines.append("")
    lines.append("Expression: TISMO `POST /rtismo/gene/downVivoExprn` for Tacstd2 and Cldn4, all six ICB treatment classes and all vivo models. Infiltrate: `POST /rtismo/gene/downICBTreated` with `type=3`. That export returns every deconvolution score (TIMER, CIBERSORT absolute, EPIC, quanTIseq, xCell, mMCPcounter) rather than the one score named in the request. Scores were joined on SRX. Spearman ρ is the rank correlation. A ρ is eligible for the maximum only when both variables have at least five distinct values, so a tie pattern such as several melanoma series at zero does not rank as ρ = 1. Cliff's δ is `2U/(n1 n2) − 1` from the Mann–Whitney U of the high-gene group versus the low-gene group. OR is the median-split odds ratio with a Woolf interval and a two-sided Fisher exact test. Bootstrap intervals are percentile intervals from 2000 resamples (seed 0).")
    lines.append("")
    lines.append("Code: `scripts/tismo_max_effect.py`. Tables: `results/tismo_max_effect/tables/`.")
    lines.append("")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _post(url: str, fields: dict) -> bytes:
    import time
    from urllib.error import HTTPError, URLError
    from urllib.parse import urlencode
    from urllib.request import Request, urlopen

    body = urlencode(fields).encode()
    last: Exception | None = None
    for attempt in range(5):
        try:
            req = Request(
                url,
                data=body,
                headers={
                    "Content-Type": "application/x-www-form-urlencoded",
                    "User-Agent": "tismo-max-effect/1.0 (public TISMO API)",
                    "Origin": "https://tismo.pku-genomics.org",
                    "Referer": "https://tismo.pku-genomics.org/",
                },
            )
            with urlopen(req, timeout=180) as resp:
                return resp.read()
        except (HTTPError, URLError, TimeoutError, OSError) as err:
            last = err
            time.sleep(2**attempt)
    raise RuntimeError(f"POST {url} failed: {last}")


def _post_json(url: str, fields: dict | None = None) -> dict:
    payload = json.loads(_post(url, fields or {}).decode())
    if payload.get("status") not in (None, 600200):
        raise RuntimeError(f"{url} status={payload.get('status')}: {payload.get('msg')}")
    return payload


def download_inputs() -> None:
    """Fetch the public TISMO tables. Existing files are reused."""
    vivo = DATA / "vivo"
    vivo.mkdir(parents=True, exist_ok=True)
    needed = [vivo / f"{gene}.csv" for gene in GENES] + [DATA / "infiltrate_icb.csv"]
    if all(path.exists() and path.stat().st_size > 1000 for path in needed):
        return
    models = [
        row["name"]
        for row in _post_json(f"{META_BASE}/gene/getVivoCohort", {})["data"]
        if row["name"] != "All"
    ]
    for gene in GENES:
        dest = vivo / f"{gene}.csv"
        if dest.exists() and dest.stat().st_size > 1000:
            continue
        raw = _post(
            f"{R_BASE}/gene/downVivoExprn",
            {
                "filename": f"{gene}.csv",
                "type": "3",
                "gene": gene,
                "icbList": json.dumps(ICB_TREATMENTS),
                "tumorList": json.dumps(models),
            },
        )
        text = raw.decode("utf-8", errors="replace")
        if not text.startswith("Samples"):
            raise RuntimeError(f"unexpected expression payload for {gene}: {text[:160]}")
        dest.write_text(text, encoding="utf-8")
    infil = DATA / "infiltrate_icb.csv"
    if not infil.exists() or infil.stat().st_size < 1000:
        raw = _post(
            f"{R_BASE}/gene/downICBTreated",
            {
                "filename": "infiltrationICB.csv",
                "immune": "TIMER",
                "type": "3",
                "immuneInfiltrate": "T cell CD8+",
                "icbList": json.dumps(ICB_TREATMENTS),
                "tumorList": json.dumps(models),
            },
        )
        text = raw.decode("utf-8", errors="replace")
        if not text.startswith("Samples"):
            raise RuntimeError(f"unexpected infiltrate payload: {text[:160]}")
        infil.write_bytes(raw)


def main():
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    download_inputs()
    immune = load_immune(DATA / "infiltrate_icb.csv")
    features = list(immune.columns)
    genes = {gene: prepare_gene(DATA / "vivo" / f"{gene}.csv", gene, immune) for gene in GENES}

    specs = infiltrate_grid(genes, features)
    specs.to_csv(TABLES / "infiltrate_specs.tsv", sep="\t", index=False)

    primary = specs["stratum"].eq("pan-cancer")
    study = specs["unit"].eq("study") & primary
    mouse = specs["unit"].eq("mouse") & primary
    residual = specs["unit"].eq("mouse_line_residual") & primary

    max_rho_study = pick_max(specs, "abs_rho", study)
    max_cliff_study = pick_max(specs, "abs_cliff", study & specs["cliff_eligible"])
    max_or_study = pick_max(specs, "OR_extreme", study & specs["OR_eligible"])
    max_rho_mouse = pick_max(specs, "abs_rho", mouse)
    max_cliff_mouse = pick_max(specs, "abs_cliff", mouse & specs["cliff_eligible"])
    max_or_mouse = pick_max(specs, "OR_extreme", mouse & specs["OR_eligible"])
    max_rho_resid = pick_max(specs, "abs_rho", residual)

    # Attach uncertainty for the three study-level winners and confirm the 2×2 from the raw medians.
    def study_frame(gene: str, feature: str) -> pd.DataFrame:
        tab = aggregate(genes[gene], "study", [feature])
        return tab.dropna(subset=["value", feature]).reset_index(drop=True)

    rho_tab = study_frame(max_rho_study["gene"], max_rho_study["feature"])
    rho_lo, rho_hi = bootstrap_spearman(rho_tab["value"], rho_tab[max_rho_study["feature"]])
    loo_min, loo_max = leave_one_out_spearman(rho_tab["value"], rho_tab[max_rho_study["feature"]])
    max_rho_study["rho_low"] = rho_lo
    max_rho_study["rho_high"] = rho_hi
    max_rho_study["loo_min"] = loo_min
    max_rho_study["loo_max"] = loo_max
    rho_tab.to_csv(TABLES / "study_max_rho.tsv", sep="\t", index=False)

    cliff_tab = study_frame(max_cliff_study["gene"], max_cliff_study["feature"])
    mx = float(cliff_tab["value"].median())
    high = cliff_tab.loc[cliff_tab["value"] >= mx, max_cliff_study["feature"]]
    low = cliff_tab.loc[cliff_tab["value"] < mx, max_cliff_study["feature"]]
    c_lo, c_hi = bootstrap_cliff(high.to_numpy(), low.to_numpy())
    max_cliff_study["cliff_low"] = c_lo
    max_cliff_study["cliff_high"] = c_hi
    cliff_tab.to_csv(TABLES / "study_max_cliff.tsv", sep="\t", index=False)

    or_tab = study_frame(max_or_study["gene"], max_or_study["feature"])
    or_check = median_split_or(or_tab["value"], or_tab[max_or_study["feature"]])
    if not or_check["eligible"] or abs(or_check["OR"] - max_or_study["OR"]) > 1e-6:
        raise RuntimeError("OR recomputation does not match the ranked spec")
    or_tab.to_csv(TABLES / "study_max_or.tsv", sep="\t", index=False)

    def sized_sensitivity(tab: pd.DataFrame, feature: str, min_mice: int = 8) -> dict:
        sub = tab[tab["n_mice"] >= min_mice]
        rho_s, p_s, n_s = spearman(sub["value"], sub[feature])
        odd_s = median_split_or(sub["value"], sub[feature])
        return {
            "rho": rho_s,
            "spearman_p": p_s,
            "n": n_s,
            "OR": odd_s["OR"] if odd_s["eligible"] else None,
            "table": f"{odd_s['a']}/{odd_s['b']}/{odd_s['c']}/{odd_s['d']}",
            "fisher_p": odd_s["fisher_p"] if odd_s["eligible"] else None,
            "eligible": bool(odd_s["eligible"]),
        }

    sens_rho = sized_sensitivity(rho_tab, max_rho_study["feature"])
    sens_or = sized_sensitivity(or_tab, max_or_study["feature"])

    plot_study_scatter(
        rho_tab,
        max_rho_study["feature"],
        max_rho_study["gene"],
        FIGS / "study_max_rho.png",
        max_rho_study["rho"],
        max_rho_study["n"],
    )
    plot_quadrants(
        or_tab,
        max_or_study["feature"],
        max_or_study["gene"],
        FIGS / "study_max_or.png",
        {"OR": or_check["OR"], "a": or_check["a"], "b": or_check["b"], "c": or_check["c"], "d": or_check["d"]},
    )
    plot_top_bars(specs, FIGS / "study_top_rho.png")

    cd8_rows = []
    for gene in GENES:
        for feat in CD8_FEATURES:
            if feat not in features:
                continue
            mouse_hit = specs[(specs.gene == gene) & (specs.feature == feat) & (specs.unit == "mouse") & (specs.stratum == "pan-cancer")]
            study_hit = specs[(specs.gene == gene) & (specs.feature == feat) & (specs.unit == "study") & (specs.stratum == "pan-cancer")]
            cd8_rows.append(
                {
                    "gene": gene,
                    "feature": feat,
                    "study_rho": float(study_hit.iloc[0]["rho"]) if len(study_hit) else float("nan"),
                    "study_n": int(study_hit.iloc[0]["n"]) if len(study_hit) else 0,
                    "mouse_rho": float(mouse_hit.iloc[0]["rho"]) if len(mouse_hit) else float("nan"),
                    "mouse_n": int(mouse_hit.iloc[0]["n"]) if len(mouse_hit) else 0,
                }
            )
    pd.DataFrame(cd8_rows).to_csv(TABLES / "cd8_reference.tsv", sep="\t", index=False)

    resp, mixed = response_grid(genes)
    resp.to_csv(TABLES / "response_specs.tsv", sep="\t", index=False)
    # Keep eligible pan-cancer rows and any stratum that clears the arm floor.
    show = resp[(resp["stratum"] == "pan-cancer") | (resp["cliff_eligible"]) | (resp["OR_eligible"])].copy()
    resp_max = pick_max(resp, "abs_cliff", resp["cliff_eligible"] & resp["unit"].eq("arm"))
    pancancer = {}
    for gene in GENES:
        hit = resp[(resp.gene == gene) & (resp.unit == "arm") & (resp.stratum == "pan-cancer")].iloc[0]
        pancancer[gene] = {"OR": float(hit["OR_response_given_high_gene"]), "cliff": float(hit["cliff_R_vs_NR"]), "rho": float(hit["rho"])}

    mixed_rows = []
    for gene, detail in mixed.items():
        for row in detail:
            mixed_rows.append({"gene": gene, **row})
    pd.DataFrame(mixed_rows).to_csv(TABLES / "mixed_arm_mice.tsv", sep="\t", index=False)

    # Line-level companion for the winning feature, so the nesting is visible.
    line_notes = []
    for label, spec in (("rho", max_rho_study), ("cliff", max_cliff_study), ("OR", max_or_study)):
        hit = specs[
            (specs.gene == spec["gene"])
            & (specs.feature == spec["feature"])
            & (specs.unit == "line")
            & (specs.stratum == "pan-cancer")
        ]
        if len(hit):
            line_notes.append(
                {
                    "paired_with": label,
                    "gene": spec["gene"],
                    "feature": spec["feature"],
                    "n_lines": int(hit.iloc[0]["n"]),
                    "rho": float(hit.iloc[0]["rho"]),
                    "cliff": float(hit.iloc[0]["cliff"]) if hit.iloc[0]["cliff"] == hit.iloc[0]["cliff"] else None,
                    "OR": float(hit.iloc[0]["OR"]) if hit.iloc[0]["OR_eligible"] else None,
                }
            )

    def clean(rec: dict) -> dict:
        out = {}
        for k, v in rec.items():
            if isinstance(v, (np.floating, float)):
                out[k] = None if math.isnan(float(v)) or math.isinf(float(v)) else float(v)
            elif isinstance(v, (np.integer,)):
                out[k] = int(v)
            elif isinstance(v, (np.bool_,)):
                out[k] = bool(v)
            else:
                out[k] = v
        return out

    headline = {
        "max_abs_rho_study": clean(max_rho_study),
        "max_abs_cliff_study": clean(max_cliff_study),
        "max_OR_study": clean(max_or_study),
        "max_abs_rho_mouse": clean(max_rho_mouse),
        "max_abs_cliff_mouse": clean(max_cliff_mouse),
        "max_OR_mouse": clean(max_or_mouse),
        "max_abs_rho_within_line": clean(max_rho_resid),
        "response_pancancer": pancancer,
        "response_max_abs_cliff": clean(resp_max) if resp_max else None,
    }
    summary = {
        "n_immune_features": len(features),
        "n_specs": int(len(specs)),
        "n_specs_study_pancancer": int(study.sum()),
        "n_rho_eligible_study": int((study & specs["rho_eligible"]).sum()),
        "n_specs_mouse_pancancer": int(mouse.sum()),
        "excluded_nonimmune_rule": list(STROMAL_TOKENS),
        "genes": {
            gene: {
                "n_unique_samples": int(df["Samples"].nunique()),
                "n_baseline_mice": int((df["Baseline"] == 1).sum()),
                "n_treated_mice": int(df["Responder"].isin(["Responders", "Non-responders"]).sum()),
                "n_studies_baseline": int(df.loc[df["Baseline"] == 1, "GSE_ID"].nunique()),
                "n_lines_baseline": int(df.loc[df["Baseline"] == 1, "line"].nunique()),
                "lung_lines": sorted(df.loc[df["cancer"] == "Lung", "line"].unique()),
            }
            for gene, df in genes.items()
        },
        "headline": headline,
        "sensitivity_n_mice_ge_8": {
            "rho_winner": {
                "gene": max_rho_study["gene"],
                "feature": max_rho_study["feature"],
                "rho": sens_rho["rho"],
                "n": sens_rho["n"],
                "spearman_p": sens_rho["spearman_p"],
            },
            "or_winner": {
                "gene": max_or_study["gene"],
                "feature": max_or_study["feature"],
                "OR": sens_or["OR"],
                "n": sens_or["n"],
                "table": sens_or["table"],
                "fisher_p": sens_or["fisher_p"],
            },
        },
        "line_level_same_features": line_notes,
        "cd8_reference": cd8_rows,
        "response_rows": json.loads(show.to_json(orient="records")),
        "mixed_arms": mixed_rows,
    }
    # json roundtrip converts numpy leftovers inside response_rows
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_report(summary, OUT / "REPORT.md")

    # Invariants that keep the headline honest.
    assert max_rho_study["unit"] == "study" and max_rho_study["stratum"] == "pan-cancer"
    assert max_or_study["OR_eligible"] is True or max_or_study["OR_eligible"] == True
    assert max_or_study["a_high_high"] >= 1 and max_or_study["d_low_low"] >= 1
    assert set(genes["Tacstd2"]["cancer"].unique()) <= set(CANCER.values()) | {"Other"}
    assert "KL" not in set(genes["Tacstd2"]["line"]) and "KP" not in set(genes["Tacstd2"]["line"])
    assert int((genes["Tacstd2"]["line"] == "LLC").sum()) > 0
    print(json.dumps({k: headline[k] for k in ("max_abs_rho_study", "max_abs_cliff_study", "max_OR_study", "max_abs_rho_mouse", "max_OR_mouse", "max_abs_rho_within_line")}, indent=2))
    print("response pancancer", pancancer)
    print("response max cliff", None if resp_max is None else {k: resp_max[k] for k in ("gene", "stratum", "cliff_R_vs_NR", "n_R", "n_NR", "rho")})
    print("wrote", OUT)


if __name__ == "__main__":
    main()
