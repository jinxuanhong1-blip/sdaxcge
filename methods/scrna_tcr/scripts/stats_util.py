"""Patient-level tests. Cells are never the unit of an MPR / residual contrast."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """Cliff's δ: P(a>b) - P(a<b). Positive => a tends larger than b."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if a.size == 0 or b.size == 0:
        return float("nan")
    # (2U)/(n1 n2) - 1 with U = #pairs a>b + 0.5 ties
    u = stats.mannwhitneyu(a, b, alternative="two-sided").statistic
    return float(2.0 * u / (a.size * b.size) - 1.0)


def mw_row(a: pd.Series, b: pd.Series, metric: str, contrast: str) -> dict:
    a = pd.to_numeric(a, errors="coerce").dropna()
    b = pd.to_numeric(b, errors="coerce").dropna()
    n_a, n_b = int(a.size), int(b.size)
    if n_a < 3 or n_b < 3:
        return {
            "metric": metric,
            "contrast": contrast,
            "n_a": n_a,
            "n_b": n_b,
            "median_a": float(a.median()) if n_a else float("nan"),
            "median_b": float(b.median()) if n_b else float("nan"),
            "mean_a": float(a.mean()) if n_a else float("nan"),
            "mean_b": float(b.mean()) if n_b else float("nan"),
            "mannwhitney_U": float("nan"),
            "mannwhitney_p_twosided": float("nan"),
            "cliffs_delta_a_minus_b": float("nan"),
            "note": "n<3 in a group; test not run",
        }
    res = stats.mannwhitneyu(a, b, alternative="two-sided")
    return {
        "metric": metric,
        "contrast": contrast,
        "n_a": n_a,
        "n_b": n_b,
        "median_a": float(a.median()),
        "median_b": float(b.median()),
        "mean_a": float(a.mean()),
        "mean_b": float(b.mean()),
        "q25_a": float(a.quantile(0.25)),
        "q75_a": float(a.quantile(0.75)),
        "q25_b": float(b.quantile(0.25)),
        "q75_b": float(b.quantile(0.75)),
        "mannwhitney_U": float(res.statistic),
        "mannwhitney_p_twosided": float(res.pvalue),
        "cliffs_delta_a_minus_b": cliffs_delta(a.to_numpy(), b.to_numpy()),
        "note": "",
    }


def spearman_row(x: pd.Series, y: pd.Series, metric_x: str, metric_y: str) -> dict:
    d = pd.DataFrame({"x": pd.to_numeric(x, errors="coerce"), "y": pd.to_numeric(y, errors="coerce")})
    d = d.dropna()
    n = int(len(d))
    if n < 5:
        return {
            "metric_x": metric_x,
            "metric_y": metric_y,
            "n": n,
            "spearman_rho": float("nan"),
            "spearman_p": float("nan"),
            "note": "n<5; correlation not run",
        }
    rho, p = stats.spearmanr(d["x"], d["y"])
    return {
        "metric_x": metric_x,
        "metric_y": metric_y,
        "n": n,
        "spearman_rho": float(rho),
        "spearman_p": float(p),
        "note": "",
    }


def paired_mw(a: pd.Series, b: pd.Series, metric: str) -> dict:
    d = pd.DataFrame({"a": pd.to_numeric(a, errors="coerce"), "b": pd.to_numeric(b, errors="coerce")})
    d = d.dropna()
    n = int(len(d))
    if n < 5:
        return {
            "metric": metric,
            "n_pairs": n,
            "median_a": float("nan"),
            "median_b": float("nan"),
            "median_a_minus_b": float("nan"),
            "wilcoxon_p_twosided": float("nan"),
            "note": "n<5 pairs; test not run",
        }
    try:
        res = stats.wilcoxon(d["a"], d["b"], alternative="two-sided", zero_method="wilcox")
        p = float(res.pvalue)
    except ValueError as e:
        p = float("nan")
        note = f"wilcoxon failed: {e}"
    else:
        note = ""
    return {
        "metric": metric,
        "n_pairs": n,
        "median_a": float(d["a"].median()),
        "median_b": float(d["b"].median()),
        "median_a_minus_b": float((d["a"] - d["b"]).median()),
        "wilcoxon_p_twosided": p,
        "note": note,
    }


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def write_json(obj, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2) + "\n")
