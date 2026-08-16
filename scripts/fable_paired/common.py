"""Shared statistics and plotting helpers for the fable_paired analysis."""
from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Sequence

import numpy as np
from scipy import stats

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# Consistent, colour-blind friendly palette.
PALETTE = {
    "responder": "#2c7fb8",
    "non_responder": "#d95f0e",
    "pre": "#8c96c6",
    "on": "#88419d",
    "post": "#810f7c",
    "TACSTD2": "#1b7837",
    "CLDN4": "#762a83",
}


@dataclass
class GroupComparison:
    dataset: str
    gene: str
    comparison: str
    group_hi: str
    group_lo: str
    n_hi: int
    n_lo: int
    median_hi: float
    median_lo: float
    mean_hi: float
    mean_lo: float
    log2fc_hi_vs_lo: float
    auc: float          # P(hi > lo); 0.5 = no separation
    test: str
    statistic: float
    p_value: float

    def as_row(self) -> dict:
        return asdict(self)


def _auc(hi: np.ndarray, lo: np.ndarray) -> float:
    """Probability that a random 'hi' value exceeds a random 'lo' value (AUC)."""
    if len(hi) == 0 or len(lo) == 0:
        return float("nan")
    try:
        u, _ = stats.mannwhitneyu(hi, lo, alternative="two-sided")
        return float(u) / (len(hi) * len(lo))
    except ValueError:
        return float("nan")


def compare_groups(
    dataset: str,
    gene: str,
    comparison: str,
    values_hi: Sequence[float],
    values_lo: Sequence[float],
    label_hi: str,
    label_lo: str,
    log_scale_input: bool = True,
) -> GroupComparison:
    """Two-group comparison on (log-scale) expression via Mann-Whitney U.

    ``log2fc_hi_vs_lo`` is the difference of means when inputs are already on a
    log2 scale, otherwise log2(mean_hi+eps) - log2(mean_lo+eps).
    """
    hi = np.asarray(values_hi, dtype=float)
    lo = np.asarray(values_lo, dtype=float)
    hi = hi[~np.isnan(hi)]
    lo = lo[~np.isnan(lo)]

    if log_scale_input:
        log2fc = float(np.mean(hi) - np.mean(lo)) if len(hi) and len(lo) else float("nan")
    else:
        eps = 1e-9
        log2fc = float(np.log2(np.mean(hi) + eps) - np.log2(np.mean(lo) + eps))

    if len(hi) >= 1 and len(lo) >= 1 and (len(hi) + len(lo)) >= 3:
        try:
            u, p = stats.mannwhitneyu(hi, lo, alternative="two-sided")
        except ValueError:
            u, p = float("nan"), float("nan")
    else:
        u, p = float("nan"), float("nan")

    return GroupComparison(
        dataset=dataset,
        gene=gene,
        comparison=comparison,
        group_hi=label_hi,
        group_lo=label_lo,
        n_hi=int(len(hi)),
        n_lo=int(len(lo)),
        median_hi=float(np.median(hi)) if len(hi) else float("nan"),
        median_lo=float(np.median(lo)) if len(lo) else float("nan"),
        mean_hi=float(np.mean(hi)) if len(hi) else float("nan"),
        mean_lo=float(np.mean(lo)) if len(lo) else float("nan"),
        log2fc_hi_vs_lo=log2fc,
        auc=_auc(hi, lo),
        test="Mann-Whitney U",
        statistic=float(u),
        p_value=float(p),
    )


def paired_test(values_pre: Sequence[float], values_on: Sequence[float]):
    """Wilcoxon signed-rank on paired samples. Returns (stat, p, median_delta)."""
    pre = np.asarray(values_pre, dtype=float)
    on = np.asarray(values_on, dtype=float)
    delta = on - pre
    delta = delta[~np.isnan(delta)]
    if len(delta) < 3 or np.all(delta == 0):
        return float("nan"), float("nan"), float(np.median(delta)) if len(delta) else float("nan")
    try:
        stat, p = stats.wilcoxon(delta)
    except ValueError:
        stat, p = float("nan"), float("nan")
    return float(stat), float(p), float(np.median(delta))


def strip_box(ax, groups: dict[str, np.ndarray], colors: dict[str, str], ylabel: str, title: str):
    """Draw a box + jittered strip plot for the provided named groups."""
    labels = list(groups.keys())
    data = [np.asarray(groups[k], dtype=float) for k in labels]
    bp = ax.boxplot(data, positions=range(len(labels)), widths=0.55,
                    showfliers=False, patch_artist=True)
    for patch, lab in zip(bp["boxes"], labels):
        patch.set_facecolor(colors.get(lab, "#cccccc"))
        patch.set_alpha(0.35)
    for med in bp["medians"]:
        med.set_color("black")
    rng = np.random.default_rng(0)
    for i, (lab, arr) in enumerate(zip(labels, data)):
        x = np.full(len(arr), i) + rng.uniform(-0.14, 0.14, len(arr))
        ax.scatter(x, arr, s=26, color=colors.get(lab, "#333333"),
                   edgecolor="white", linewidth=0.5, zorder=3)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=0)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)


def pstars(p: float) -> str:
    if p != p:  # NaN
        return "n/a"
    if p < 1e-3:
        return "***"
    if p < 1e-2:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def savefig(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig ] {path}")
