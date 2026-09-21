#!/usr/bin/env python3
"""Signed IFN/APM concordance for the true public CLDN4-loss sets.

Inputs are the locked Claim C4 gene effects (PR #88), not a re-fit.
Positive log2FC is KD/KO versus control and matches the thesis that
CLDN4 loss opens IFN / MHC-I / APM.

Included contrasts
  GSE207704  T47D and MCF-7 CRISPR KO     cancer, breast
  GSE22493   SKOV-3 siRNA, probes mapped  cancer, ovary
  GSE50927   naive lung KO vs WT only     non-cancer, lung

GSE22493 is dropped if either panel is under 50% observed. VILI arms
of GSE50927 are not part of this match wave.
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
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon

ROOT = Path(__file__).resolve().parents[2]
LOCK = ROOT / "methods" / "cldn4_kd_match_wave" / "locked"
OUT_TAB = ROOT / "methods" / "cldn4_kd_match_wave" / "tables"
OUT_FIG = ROOT / "methods" / "cldn4_kd_match_wave" / "figures"
FINDING = ROOT / "methods" / "cldn4_kd_match_wave" / "FINDING.md"

B_RESAMPLE = 10000
RNG_SEED = 50927
MAP_MIN_COVERAGE = 0.50
Z95 = 1.959963984540054

CONTRASTS = [
    {
        "dataset": "GSE50927",
        "contrast": "baseline KO vs WT",
        "organ": "lung",
        "species": "mouse",
        "perturbation": "germline Cldn4 KO",
        "class": "non-cancer",
        "short": "GSE50927   lung baseline",
        "note": (
            "Author edgeR naive KO vs WT (GSM1232581 vs GSM1232580). "
            "GEO lists one library per arm. VILI arms are excluded. Not a tumor."
        ),
    },
    {
        "dataset": "GSE207704",
        "contrast": "T47D KO vs WT",
        "organ": "breast",
        "species": "human",
        "perturbation": "CRISPR CLDN4 KO",
        "class": "cancer",
        "short": "GSE207704  T47D breast",
        "note": (
            "Deposited condition-mean FPKM. Many MHC-I/APM genes are absent "
            "from the file. No replicate-level SE."
        ),
    },
    {
        "dataset": "GSE207704",
        "contrast": "MCF-7 KO vs WT",
        "organ": "breast",
        "species": "human",
        "perturbation": "CRISPR CLDN4 KO",
        "class": "cancer",
        "short": "GSE207704  MCF-7 breast",
        "note": (
            "Deposited condition-mean FPKM. Panel coverage is incomplete. "
            "No replicate-level SE."
        ),
    },
    {
        "dataset": "GSE22493",
        "contrast": "SKOV-3 KD",
        "organ": "ovary",
        "species": "human",
        "perturbation": "lentiviral CLDN4 siRNA",
        "class": "cancer",
        "short": "GSE22493   SKOV-3 ovary",
        "note": (
            "Included because GPL10555 probes map to symbols "
            "(median-collapsed; G1P2 read as ISG15). "
            "Control is CLDN4 overexpression, not scramble. Three two-colour arrays."
        ),
    },
]

MODULES = ["IFN_ISG", "MHC_APM"]
MODULE_TITLE = {"IFN_ISG": "IFN / ISG", "MHC_APM": "MHC-I / APM"}
CANCER_COLOR = "#B85C38"
NONCANCER_COLOR = "#2C6E8A"
RANGE_COLOR = "#5C534C"
THESIS = "up_after_CLDN4_loss"


def read_tsv(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fields})


def fnum(value: float) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    return format(float(value), ".10g")


def pnum(text: str) -> str:
    if text in ("", "NA", "nan", None):
        return ""
    return text


def call_vs_thesis(direction: str) -> str:
    if direction == "up":
        return "concordant"
    if direction == "down":
        return "discordant"
    if direction in ("weak_up", "weak_down", "null"):
        return "unsigned"
    raise RuntimeError(f"unexpected locked call: {direction}")


def resample_median(values: np.ndarray, rng: np.random.Generator) -> tuple[float, float, float, float]:
    observed = np.asarray(values, dtype=float)
    observed = observed[np.isfinite(observed)]
    n = observed.size
    if n == 0:
        raise RuntimeError("empty panel")
    point = float(np.median(observed))
    draws = rng.integers(0, n, size=(B_RESAMPLE, n))
    medians = np.median(observed[draws], axis=1)
    lo, hi = np.quantile(medians, [0.025, 0.975])
    sd = float(np.std(medians, ddof=1))
    return point, float(lo), float(hi), sd


def der_simonian_laird(effects: list[float], ses: list[float]) -> dict:
    y = np.asarray(effects, dtype=float)
    se = np.asarray(ses, dtype=float)
    if np.any(~np.isfinite(se)) or np.any(se <= 0):
        raise RuntimeError("SE must be positive and finite")
    variance = se**2
    weight = 1.0 / variance
    fixed = float(np.sum(weight * y) / np.sum(weight))
    Q = float(np.sum(weight * (y - fixed) ** 2))
    k = int(y.size)
    df = k - 1
    C = float(np.sum(weight) - np.sum(weight**2) / np.sum(weight))
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0
    re_weight = 1.0 / (variance + tau2)
    mu = float(np.sum(re_weight * y) / np.sum(re_weight))
    se_mu = float(np.sqrt(1.0 / np.sum(re_weight)))
    I2 = max(0.0, (Q - df) / Q) if Q > 0 else 0.0
    return {
        "mu": mu,
        "se": se_mu,
        "lo": mu - Z95 * se_mu,
        "hi": mu + Z95 * se_mu,
        "tau2": tau2,
        "Q": Q,
        "I2": I2,
        "k": k,
    }


def load_inputs() -> tuple[list[dict], dict]:
    effects = read_tsv(LOCK / "gene_effects_ifn_apm.tsv")
    summary_rows = read_tsv(LOCK / "geneset_summary_ifn_apm.tsv")
    summary = {(row["dataset"], row["contrast"], row["gene_set"]): row for row in summary_rows}
    return effects, summary


def panel_values(effects: list[dict], dataset: str, contrast: str, module: str) -> list[tuple[str, float]]:
    found = []
    for row in effects:
        if row["dataset"] != dataset or row["contrast"] != contrast or row["gene_set"] != module:
            continue
        raw = row["log2FC"]
        if raw in ("", "NA", "nan"):
            continue
        found.append((row["human_symbol"], float(raw)))
    return found


def build_rows(effects: list[dict], summary: dict) -> tuple[list[dict], list[dict]]:
    rng = np.random.default_rng(RNG_SEED)
    rows = []
    excluded = []
    for meta in CONTRASTS:
        for module in MODULES:
            key = (meta["dataset"], meta["contrast"], module)
            locked = summary[key]
            measured = panel_values(effects, meta["dataset"], meta["contrast"], module)
            values = np.array([value for _, value in measured], dtype=float)
            n_in_set = int(locked["n_in_set"])
            n_observed = int(values.size)
            coverage = n_observed / n_in_set
            if meta["dataset"] == "GSE22493" and coverage < MAP_MIN_COVERAGE:
                excluded.append({**meta, "module": module, "coverage": coverage})
                continue
            point, lo, hi, sd = resample_median(values, rng)
            n_up = int(np.sum(values > 0))
            n_down = int(np.sum(values < 0))
            mean_fc = float(np.mean(values))
            locked_median = float(locked["median_log2FC"])
            locked_mean = float(locked["mean_log2FC"])
            locked_n = int(locked["n_observed"])
            locked_up = int(locked["n_up"])
            if n_observed != locked_n or n_up != locked_up:
                raise RuntimeError(f"count drift {key}: {n_observed}/{n_up} vs {locked_n}/{locked_up}")
            if not math.isclose(point, locked_median, rel_tol=0, abs_tol=1e-9):
                raise RuntimeError(f"median drift {key}: {point} vs {locked_median}")
            if not math.isclose(mean_fc, locked_mean, rel_tol=0, abs_tol=1e-9):
                raise RuntimeError(f"mean drift {key}: {mean_fc} vs {locked_mean}")
            direction = locked["direction_call"]
            mapping_gate = "pass" if meta["dataset"] != "GSE22493" else "pass_mapped"
            if meta["dataset"] == "GSE22493":
                mapping_gate = "pass_mapped"
            rows.append(
                {
                    "accession": meta["dataset"],
                    "contrast": meta["contrast"],
                    "class": meta["class"],
                    "organ": meta["organ"],
                    "species": meta["species"],
                    "perturbation": meta["perturbation"],
                    "module": module,
                    "cldn4_log2FC": float(locked["CLDN4_log2FC"]),
                    "n_in_set": n_in_set,
                    "n_observed": n_observed,
                    "coverage": coverage,
                    "n_up": n_up,
                    "n_down": n_down,
                    "median_log2FC": point,
                    "mean_log2FC": mean_fc,
                    "resample_p2_5": lo,
                    "resample_p97_5": hi,
                    "resample_sd": sd,
                    "signed_balance": (n_up - n_down) / n_observed,
                    "median_sign": 1 if point > 0 else (-1 if point < 0 else 0),
                    "thesis": THESIS,
                    "thesis_sign_match": "yes" if point > 0 else "no",
                    "locked_direction_call": direction,
                    "call_vs_thesis": call_vs_thesis(direction),
                    "wilcoxon_p_greater": pnum(locked["wilcoxon_p_greater"]),
                    "wilcoxon_p_less": pnum(locked["wilcoxon_p_less"]),
                    "sign_p_greater": pnum(locked["sign_p_greater"]),
                    "mapping_gate": mapping_gate,
                    "short": meta["short"],
                    "note": meta["note"],
                }
            )
    if excluded:
        raise RuntimeError(f"GSE22493 failed the mapping gate: {excluded}")
    return rows, excluded


def stratum_rows(rows: list[dict]) -> list[dict]:
    out = []
    for module in MODULES:
        for klass in ("non-cancer", "cancer"):
            block = [row for row in rows if row["module"] == module and row["class"] == klass]
            effects = [row["median_log2FC"] for row in block]
            if klass == "cancer":
                dl = der_simonian_laird(effects, [row["resample_sd"] for row in block])
                method = "median_of_medians_plus_descriptive_DL"
                re_mu, re_lo, re_hi = dl["mu"], dl["lo"], dl["hi"]
                I2, tau2, Q, k = dl["I2"], dl["tau2"], dl["Q"], dl["k"]
            else:
                only = block[0]
                method = "single_contrast"
                re_mu = re_lo = re_hi = ""
                I2 = tau2 = Q = ""
                k = 1
                dl = None
            out.append(
                {
                    "class": klass,
                    "module": module,
                    "n_contrasts": len(block),
                    "n_median_positive": sum(row["median_sign"] > 0 for row in block),
                    "n_median_negative": sum(row["median_sign"] < 0 for row in block),
                    "n_call_concordant": sum(row["call_vs_thesis"] == "concordant" for row in block),
                    "n_call_discordant": sum(row["call_vs_thesis"] == "discordant" for row in block),
                    "n_call_unsigned": sum(row["call_vs_thesis"] == "unsigned" for row in block),
                    "median_of_medians": float(np.median(effects)),
                    "min_median": min(effects),
                    "max_median": max(effects),
                    "dl_mu": "" if dl is None else dl["mu"],
                    "dl_wald_lo": "" if dl is None else dl["lo"],
                    "dl_wald_hi": "" if dl is None else dl["hi"],
                    "dl_I2": "" if dl is None else dl["I2"],
                    "dl_tau2": "" if dl is None else dl["tau2"],
                    "dl_Q": "" if dl is None else dl["Q"],
                    "method": method,
                    "dl_note": (
                        ""
                        if dl is None
                        else (
                            "Descriptive DerSimonian-Laird using gene-resample SDs as weights. "
                            "Not a replicate-level confidence interval and not a confirmatory test."
                        )
                    ),
                    "_re_mu": re_mu,
                    "_I2": I2,
                    "_tau2": tau2,
                    "_Q": Q,
                    "_k": k,
                    "_re_lo": re_lo,
                    "_re_hi": re_hi,
                }
            )
    return out


def gene_tables(effects: list[dict]) -> tuple[list[dict], list[dict], dict]:
    by_key: dict[tuple, dict] = {}
    symbols = {module: set() for module in MODULES}
    for row in effects:
        if row["gene_set"] not in MODULES:
            continue
        raw = row["log2FC"]
        value = None if raw in ("", "NA", "nan") else float(raw)
        key = (row["gene_set"], row["human_symbol"])
        slot = by_key.setdefault(
            key,
            {
                "module": row["gene_set"],
                "human_symbol": row["human_symbol"],
                "lung": None,
                "T47D": None,
                "MCF7": None,
                "SKOV3": None,
            },
        )
        if row["dataset"] == "GSE50927":
            slot["lung"] = value
        elif row["contrast"] == "T47D KO vs WT":
            slot["T47D"] = value
        elif row["contrast"] == "MCF-7 KO vs WT":
            slot["MCF7"] = value
        elif row["contrast"] == "SKOV-3 KD":
            slot["SKOV3"] = value
        symbols[row["gene_set"]].add(row["human_symbol"])

    signed = []
    cross = []
    for module in MODULES:
        names = sorted(symbols[module])
        for name in names:
            slot = by_key[(module, name)]
            cancer_vals = [slot[key] for key in ("T47D", "MCF7", "SKOV3") if slot[key] is not None]
            lung = slot["lung"]
            cancer_median = float(np.median(cancer_vals)) if cancer_vals else None
            rec = {
                "module": module,
                "human_symbol": name,
                "log2FC_lung_baseline": "" if lung is None else lung,
                "log2FC_T47D": "" if slot["T47D"] is None else slot["T47D"],
                "log2FC_MCF7": "" if slot["MCF7"] is None else slot["MCF7"],
                "log2FC_SKOV3": "" if slot["SKOV3"] is None else slot["SKOV3"],
                "n_cancer_observed": len(cancer_vals),
                "cancer_median_log2FC": "" if cancer_median is None else cancer_median,
                "lung_matches_thesis": "" if lung is None else ("yes" if lung > 0 else "no"),
                "cancer_matches_thesis": (
                    "" if cancer_median is None else ("yes" if cancer_median > 0 else "no")
                ),
            }
            signed.append(rec)
            if lung is None or cancer_median is None or lung == 0 or cancer_median == 0:
                continue
            lung_sign = 1 if lung > 0 else -1
            cancer_sign = 1 if cancer_median > 0 else -1
            cross.append(
                {
                    **rec,
                    "lung_sign": lung_sign,
                    "cancer_sign": cancer_sign,
                    "same_sign": "yes" if lung_sign == cancer_sign else "no",
                    "lung_up_cancer_down": "yes" if lung_sign > 0 and cancer_sign < 0 else "no",
                }
            )

    summary = []
    summary_obj = {}
    for module in MODULES:
        block = [row for row in cross if row["module"] == module]
        n = len(block)
        n_same = sum(row["same_sign"] == "yes" for row in block)
        n_split = sum(row["lung_up_cancer_down"] == "yes" for row in block)
        n_lung_up = sum(row["lung_sign"] > 0 for row in block)
        n_cancer_up = sum(row["cancer_sign"] > 0 for row in block)
        n_both_up = sum(row["lung_sign"] > 0 and row["cancer_sign"] > 0 for row in block)
        rec = {
            "module": module,
            "n_genes_in_lung_and_cancer": n,
            "n_same_sign": n_same,
            "n_opposite_sign": n - n_same,
            "frac_same_sign": n_same / n,
            "n_lung_up": n_lung_up,
            "n_cancer_median_up": n_cancer_up,
            "n_both_up": n_both_up,
            "n_lung_up_cancer_down": n_split,
            "frac_lung_up_cancer_down": n_split / n,
        }
        summary.append(rec)
        summary_obj[module] = rec
    return signed, cross, summary_obj


def lookup_gene(effects: list[dict], dataset: str, contrast: str, module: str, symbol: str) -> float | None:
    for row in effects:
        if (
            row["dataset"] == dataset
            and row["contrast"] == contrast
            and row["gene_set"] == module
            and row["human_symbol"] == symbol
            and row["log2FC"] not in ("", "NA", "nan")
        ):
            return float(row["log2FC"])
    return None


def forest_table(rows: list[dict], strata: list[dict]) -> list[dict]:
    table = []
    for row in rows:
        table.append(
            {
                "row_type": "contrast",
                "accession": row["accession"],
                "contrast": row["contrast"],
                "class": row["class"],
                "organ": row["organ"],
                "module": row["module"],
                "label": f"{'Non-cancer' if row['class'] == 'non-cancer' else 'Cancer':<11}{row['short']}",
                "estimate": row["median_log2FC"],
                "lo": row["resample_p2_5"],
                "hi": row["resample_p97_5"],
                "interval": "gene_resample_2.5_97.5",
                "n_up": row["n_up"],
                "n_observed": row["n_observed"],
                "call_vs_thesis": row["call_vs_thesis"],
                "thesis_sign_match": row["thesis_sign_match"],
            }
        )
    for stratum in strata:
        if stratum["class"] != "cancer":
            continue
        table.append(
            {
                "row_type": "cancer_median",
                "accession": "cancer_median",
                "contrast": "median of 3 cancer medians",
                "class": "cancer",
                "organ": "breast+ovary",
                "module": stratum["module"],
                "label": "Cancer      median of medians",
                "estimate": stratum["median_of_medians"],
                "lo": stratum["min_median"],
                "hi": stratum["max_median"],
                "interval": "min_max_of_contrast_medians",
                "n_up": "",
                "n_observed": stratum["n_contrasts"],
                "call_vs_thesis": (
                    f"{stratum['n_call_concordant']} concordant / "
                    f"{stratum['n_call_discordant']} discordant / "
                    f"{stratum['n_call_unsigned']} unsigned"
                ),
                "thesis_sign_match": "no" if stratum["median_of_medians"] < 0 else "yes",
            }
        )
    return table


def draw_forest(plot_rows: list[dict]) -> None:
    OUT_FIG.mkdir(parents=True, exist_ok=True)
    key_order = [
        ("contrast", "baseline KO vs WT"),
        ("contrast", "T47D KO vs WT"),
        ("contrast", "MCF-7 KO vs WT"),
        ("contrast", "SKOV-3 KD"),
        ("cancer_median", "median of 3 cancer medians"),
    ]
    y_pos = {key: float(len(key_order) - 1 - i) for i, key in enumerate(key_order)}

    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.15), sharey=True)
    all_lo = min(float(row["lo"]) for row in plot_rows)
    all_hi = max(float(row["hi"]) for row in plot_rows)
    pad = 0.08 * (all_hi - all_lo)
    xlim = (all_lo - pad, all_hi + pad)

    for ax, module in zip(axes, MODULES):
        block = [row for row in plot_rows if row["module"] == module]
        ax.axvspan(0, xlim[1], color="#E7F3EC", zorder=0)
        ax.axvspan(xlim[0], 0, color="#F8EFEA", zorder=0)
        ax.axvline(0, color="#2B2B2B", lw=0.9, zorder=1)
        for row in block:
            y = y_pos[(row["row_type"], row["contrast"])]
            est, lo, hi = float(row["estimate"]), float(row["lo"]), float(row["hi"])
            if row["row_type"] == "cancer_median":
                height = 0.16
                verts = [(lo, y), (est, y + height), (hi, y), (est, y - height)]
                ax.add_patch(
                    Polygon(verts, closed=True, facecolor=RANGE_COLOR, edgecolor=RANGE_COLOR, zorder=3)
                )
                ax.annotate(
                    f"{est:+.2f}   range of medians",
                    xy=(1.01, y),
                    xycoords=ax.get_yaxis_transform(),
                    va="center",
                    ha="left",
                    fontsize=8,
                    fontfamily="DejaVu Sans",
                    color="#333333",
                )
                continue
            color = NONCANCER_COLOR if row["class"] == "non-cancer" else CANCER_COLOR
            ax.plot([lo, hi], [y, y], color=color, lw=1.6, solid_capstyle="round", zorder=2)
            ax.plot([lo, lo], [y - 0.08, y + 0.08], color=color, lw=1.2, zorder=2)
            ax.plot([hi, hi], [y - 0.08, y + 0.08], color=color, lw=1.2, zorder=2)
            ax.scatter([est], [y], s=42, color=color, zorder=4, linewidths=0)
            ax.annotate(
                f"{est:+.2f}   {row['n_up']}/{row['n_observed']}   {row['call_vs_thesis']}",
                xy=(1.01, y),
                xycoords=ax.get_yaxis_transform(),
                va="center",
                ha="left",
                fontsize=8,
                fontfamily="DejaVu Sans",
                color="#333333",
            )
        ax.set_xlim(*xlim)
        ax.set_ylim(-0.55, len(key_order) - 0.45)
        ax.set_title(MODULE_TITLE[module], loc="left", fontsize=12, pad=8)
        ax.set_xlabel("Median log2FC after CLDN4 loss", fontsize=9)
        ax.tick_params(axis="x", labelsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.text(0.01, 1.02, "closes", transform=ax.transAxes, ha="left", va="bottom", fontsize=8, color="#8A4B32")
        ax.text(0.99, 1.02, "opens (thesis)", transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color="#1E6B45")

    labels = [
        "Non-cancer   GSE50927  lung",
        "Cancer       GSE207704 T47D breast",
        "Cancer       GSE207704 MCF-7 breast",
        "Cancer       GSE22493  SKOV-3 ovary",
        "Cancer       median of the 3 medians",
    ]
    axes[0].set_yticks([y_pos[key] for key in key_order])
    axes[0].set_yticklabels(labels, fontsize=8.5, fontfamily="DejaVu Sans")
    for tick, key in zip(axes[0].get_yticklabels(), key_order):
        tick.set_color(NONCANCER_COLOR if key[1].startswith("baseline") else CANCER_COLOR)
    axes[0].axhline(y_pos[key_order[0]] - 0.5, color="#D0D0D0", lw=0.6, zorder=0)
    axes[1].axhline(y_pos[key_order[0]] - 0.5, color="#D0D0D0", lw=0.6, zorder=0)

    legend_handles = [
        Line2D([0], [0], marker="o", color=NONCANCER_COLOR, lw=1.6, label="Non-cancer"),
        Line2D([0], [0], marker="o", color=CANCER_COLOR, lw=1.6, label="Cancer"),
        Line2D([0], [0], marker="D", color=RANGE_COLOR, lw=0, markersize=6, label="Cancer median (min–max)"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.46, -0.02),
        fontsize=8.5,
    )
    fig.suptitle(
        "True CLDN4-loss sets: signed IFN / APM versus the open-after-loss thesis",
        fontsize=12.5,
        x=0.01,
        ha="left",
    )
    fig.text(
        0.01,
        -0.06,
        "Thin bars resample genes inside the locked panel median (2.5–97.5%). "
        "They are not replicate confidence intervals. "
        "The diamond spans the min and max of the three cancer medians. "
        "GSE22493 is in because the array probes map. GSE50927 is the uninjured baseline only.",
        fontsize=7.5,
        color="#444444",
        ha="left",
    )
    fig.subplots_adjust(left=0.24, right=0.78, top=0.84, bottom=0.16, wspace=0.55)
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(OUT_FIG / f"forest_ifn_apm.{suffix}", dpi=160)
    plt.close(fig)


def fmt(value: float, digits: int = 2) -> str:
    return f"{value:+.{digits}f}"


def write_finding(rows: list[dict], strata: list[dict], cross_summary: dict, h2k1: float) -> None:
    by = {(row["accession"], row["contrast"], row["module"]): row for row in rows}
    lung_ifn = by[("GSE50927", "baseline KO vs WT", "IFN_ISG")]
    lung_apm = by[("GSE50927", "baseline KO vs WT", "MHC_APM")]
    cancer_ifn = next(s for s in strata if s["class"] == "cancer" and s["module"] == "IFN_ISG")
    cancer_apm = next(s for s in strata if s["class"] == "cancer" and s["module"] == "MHC_APM")

    def cell(accession: str, contrast: str, module: str) -> str:
        row = by[(accession, contrast, module)]
        return (
            f"{fmt(row['median_log2FC'])} "
            f"({row['n_up']}/{row['n_observed']}; {row['call_vs_thesis']})"
        )

    lines = [
        "# CLDN4 KD match wave: signed IFN/APM concordance",
        "",
        "Thesis sign, taken from the public C4 test: **CLDN4 loss opens IFN and MHC-I/APM**. "
        "The signed effect is the median log2 fold change of the prespecified panel "
        "(KD/KO versus matched control). Positive matches the thesis.",
        "",
        "This wave combines only the true CLDN4-loss transcriptomes. "
        "TACSTD2/TROP2 knockdowns are not included. "
        "GSE22493 stays in because GPL10555 probes map to gene symbols "
        f"(IFN/ISG {by[('GSE22493','SKOV-3 KD','IFN_ISG')]['n_observed']}/"
        f"{by[('GSE22493','SKOV-3 KD','IFN_ISG')]['n_in_set']}, "
        f"APM {by[('GSE22493','SKOV-3 KD','MHC_APM')]['n_observed']}/"
        f"{by[('GSE22493','SKOV-3 KD','MHC_APM')]['n_in_set']}; both above the 50% gate). "
        "GSE50927 contributes the uninjured baseline contrast only. VILI is not in the forest.",
        "",
        "Numbers are the locked Claim C4 effects (PR #88), recomputed here only to attach "
        "gene-resample intervals and the cancer-versus-non-cancer labels. "
        "Medians, means, and up-counts are asserted against that table.",
        "",
        "## Concordance",
        "",
        "| Class | Organ | Contrast | CLDN4 log2FC | IFN/ISG | MHC-I/APM |",
        "|---|---|---|---:|---|---|",
        f"| Non-cancer | Lung | GSE50927 naive KO vs WT | {fmt(lung_ifn['cldn4_log2FC'])} | {cell('GSE50927','baseline KO vs WT','IFN_ISG')} | {cell('GSE50927','baseline KO vs WT','MHC_APM')} |",
        f"| Cancer | Breast | GSE207704 T47D KO | {fmt(by[('GSE207704','T47D KO vs WT','IFN_ISG')]['cldn4_log2FC'])} | {cell('GSE207704','T47D KO vs WT','IFN_ISG')} | {cell('GSE207704','T47D KO vs WT','MHC_APM')} |",
        f"| Cancer | Breast | GSE207704 MCF-7 KO | {fmt(by[('GSE207704','MCF-7 KO vs WT','IFN_ISG')]['cldn4_log2FC'])} | {cell('GSE207704','MCF-7 KO vs WT','IFN_ISG')} | {cell('GSE207704','MCF-7 KO vs WT','MHC_APM')} |",
        f"| Cancer | Ovary | GSE22493 SKOV-3 KD | {fmt(by[('GSE22493','SKOV-3 KD','IFN_ISG')]['cldn4_log2FC'])} | {cell('GSE22493','SKOV-3 KD','IFN_ISG')} | {cell('GSE22493','SKOV-3 KD','MHC_APM')} |",
        "",
        "A locked call is **concordant** only when Claim C4 called the panel `up` "
        "(median above 0 and a one-sided Wilcoxon or sign test P ≤ 0.05). "
        "**Discordant** means the locked call was `down`. "
        "**Unsigned** means the median did not clear that bar (`weak_down` here).",
        "",
        f"Non-cancer lung is concordant on both panels: IFN/ISG median {fmt(lung_ifn['median_log2FC'])} "
        f"({lung_ifn['n_up']}/{lung_ifn['n_observed']} up) and MHC-I/APM median {fmt(lung_apm['median_log2FC'])} "
        f"({lung_apm['n_up']}/{lung_apm['n_observed']} up). "
        f"Classical H2-K1, the mouse stand-in for HLA-A, is {fmt(h2k1)} in that same baseline table, "
        "so the APM concordance is the panel median, not MHC-I heavy-chain induction.",
        "",
        "Every cancer median is negative, in breast and in ovary, on both panels. "
        f"IFN/ISG locked calls: {cancer_ifn['n_call_concordant']} concordant, "
        f"{cancer_ifn['n_call_discordant']} discordant, {cancer_ifn['n_call_unsigned']} unsigned. "
        f"MHC-I/APM locked calls: {cancer_apm['n_call_concordant']} concordant, "
        f"{cancer_apm['n_call_discordant']} discordant, {cancer_apm['n_call_unsigned']} unsigned. "
        "No cancer contrast matches the thesis call. "
        f"The median of the three cancer medians is {fmt(cancer_ifn['median_of_medians'])} for IFN/ISG "
        f"(range {fmt(cancer_ifn['min_median'])} to {fmt(cancer_ifn['max_median'])}) and "
        f"{fmt(cancer_apm['median_of_medians'])} for MHC-I/APM "
        f"(range {fmt(cancer_apm['min_median'])} to {fmt(cancer_apm['max_median'])}).",
        "",
        "T47D clears a down call on both panels. "
        f"Its IFN/ISG gene-resample interval still touches 0 "
        f"({fmt(by[('GSE207704','T47D KO vs WT','IFN_ISG')]['resample_p2_5'])} to "
        f"{fmt(by[('GSE207704','T47D KO vs WT','IFN_ISG')]['resample_p97_5'])}). "
        "That call is the locked Wilcoxon on the gene list (7 up, 14 down), "
        "which answers a different question than the median's resample interval. "
        "SKOV-3 clears a down call for APM and stays unsigned for IFN/ISG "
        "(median essentially zero, slightly negative). "
        "MCF-7 is unsigned on both, with a negative median and a thin gene list.",
        "",
        "## Cross-organ gene signs",
        "",
        "Genes counted below are observed in the lung baseline and in at least one cancer contrast. "
        "The cancer sign is the median of the cancer log2 fold changes available for that gene.",
        "",
        "| Module | Genes in both | Same sign | Lung up, cancer down | Both up |",
        "|---|---:|---:|---:|---:|",
    ]
    for module in MODULES:
        rec = cross_summary[module]
        lines.append(
            f"| {MODULE_TITLE[module]} | {rec['n_genes_in_lung_and_cancer']} | "
            f"{rec['n_same_sign']} | {rec['n_lung_up_cancer_down']} | {rec['n_both_up']} |"
        )
    lines += [
        "",
        "Same-sign agreement across organ class is the minority. "
        "The common pattern among shared genes is lung up and cancer down.",
        "",
        "## What the forest is",
        "",
        "Thin bars are 2.5–97.5% gene-resample intervals of the locked panel median "
        f"({B_RESAMPLE} draws, seed {RNG_SEED}). "
        "They describe whether the genes in that panel move together. "
        "They are not confidence intervals from biological replicates. "
        "GSE207704 deposits pooled FPKM. GSE50927 baseline is one naive KO library versus one naive WT library. "
        "GSE22493 has three arrays, but the interval still resamples genes, so it matches the other rows.",
        "",
        "The diamond is the median of the three cancer medians. Its width is the min-to-max of those medians. "
        "Cancer and non-cancer stay on separate rows.",
        "",
        "A DerSimonian–Laird summary that weights each cancer median by its gene-resample SD is in "
        "`tables/stratum_summary.tsv`. "
        f"IFN/ISG: {fmt(cancer_ifn['dl_mu'])} "
        f"(Wald {fmt(cancer_ifn['dl_wald_lo'])} to {fmt(cancer_ifn['dl_wald_hi'])}). "
        f"MHC-I/APM: {fmt(cancer_apm['dl_mu'])} "
        f"(Wald {fmt(cancer_apm['dl_wald_lo'])} to {fmt(cancer_apm['dl_wald_hi'])}). "
        "I² is 0 on both because those SDs are wide relative to the gap between medians. "
        "The weights are gene-resample SDs, so that summary stays in the table and off the figure.",
        "",
        "## Files",
        "",
        "- `tables/concordance_ifn_apm.tsv` — one row per contrast and panel",
        "- `tables/gene_signed.tsv` — gene-level log2FC by contrast",
        "- `tables/cross_organ_genes.tsv` — genes observed on both sides of the class split",
        "- `tables/cross_organ_summary.tsv`",
        "- `tables/stratum_summary.tsv`",
        "- `tables/forest_effects.tsv`",
        "- `figures/forest_ifn_apm.png`",
        "",
        "Reproduce: `python3 scripts/cldn4_kd_match_wave/build_concordance.py`",
        "",
    ]
    FINDING.write_text("\n".join(lines))


def assert_story(rows: list[dict], cross_summary: dict) -> None:
    if len(rows) != 8:
        raise RuntimeError(f"expected 8 rows, got {len(rows)}")
    if any("VILI" in row["contrast"] for row in rows):
        raise RuntimeError("VILI leaked into the match wave")
    classes = {row["class"] for row in rows}
    if classes != {"cancer", "non-cancer"}:
        raise RuntimeError(classes)
    lung = [row for row in rows if row["accession"] == "GSE50927"]
    if len(lung) != 2 or any(row["call_vs_thesis"] != "concordant" or row["median_sign"] != 1 for row in lung):
        raise RuntimeError("lung baseline should be concordant on both panels")
    cancer = [row for row in rows if row["class"] == "cancer"]
    if any(row["median_sign"] >= 0 for row in cancer):
        raise RuntimeError("a cancer median is not negative")
    if any(row["call_vs_thesis"] == "concordant" for row in cancer):
        raise RuntimeError("a cancer contrast was called concordant")
    if not any(row["accession"] == "GSE22493" for row in rows):
        raise RuntimeError("mapped GSE22493 missing")
    for row in rows:
        if row["cldn4_log2FC"] >= 0:
            raise RuntimeError("CLDN4 did not fall in a contrast")
    for module, rec in cross_summary.items():
        if rec["n_genes_in_lung_and_cancer"] < 5:
            raise RuntimeError(f"cross-organ overlap too small for {module}")
    png = OUT_FIG / "forest_ifn_apm.png"
    if not png.exists() or png.stat().st_size < 20000:
        raise RuntimeError("forest png missing or too small")


def dump_numeric(rows: list[dict]) -> list[dict]:
    numeric = [
        "cldn4_log2FC",
        "coverage",
        "median_log2FC",
        "mean_log2FC",
        "resample_p2_5",
        "resample_p97_5",
        "resample_sd",
        "signed_balance",
    ]
    out = []
    for row in rows:
        item = dict(row)
        for key in numeric:
            item[key] = fnum(item[key])
        for key in ("n_in_set", "n_observed", "n_up", "n_down", "median_sign"):
            item[key] = str(item[key])
        out.append(item)
    return out


def main() -> None:
    effects, summary = load_inputs()
    rows, _excluded = build_rows(effects, summary)
    strata = stratum_rows(rows)
    signed, cross, cross_summary = gene_tables(effects)
    h2k1 = lookup_gene(effects, "GSE50927", "baseline KO vs WT", "MHC_APM", "HLA-A")
    if h2k1 is None:
        raise RuntimeError("H2-K1 / HLA-A missing from lung APM panel")
    plot_rows = forest_table(rows, strata)
    draw_forest(plot_rows)
    write_finding(rows, strata, cross_summary, h2k1)

    conc_fields = [
        "accession",
        "contrast",
        "class",
        "organ",
        "species",
        "perturbation",
        "module",
        "cldn4_log2FC",
        "n_in_set",
        "n_observed",
        "coverage",
        "n_up",
        "n_down",
        "median_log2FC",
        "mean_log2FC",
        "resample_p2_5",
        "resample_p97_5",
        "resample_sd",
        "signed_balance",
        "median_sign",
        "thesis",
        "thesis_sign_match",
        "locked_direction_call",
        "call_vs_thesis",
        "wilcoxon_p_greater",
        "wilcoxon_p_less",
        "sign_p_greater",
        "mapping_gate",
        "note",
    ]
    write_tsv(OUT_TAB / "concordance_ifn_apm.tsv", dump_numeric(rows), conc_fields)

    stratum_fields = [
        "class",
        "module",
        "n_contrasts",
        "n_median_positive",
        "n_median_negative",
        "n_call_concordant",
        "n_call_discordant",
        "n_call_unsigned",
        "median_of_medians",
        "min_median",
        "max_median",
        "dl_mu",
        "dl_wald_lo",
        "dl_wald_hi",
        "dl_I2",
        "dl_tau2",
        "dl_Q",
        "method",
        "dl_note",
    ]
    stratum_out = []
    for row in strata:
        item = {key: row[key] for key in stratum_fields}
        for key in (
            "median_of_medians",
            "min_median",
            "max_median",
            "dl_mu",
            "dl_wald_lo",
            "dl_wald_hi",
            "dl_I2",
            "dl_tau2",
            "dl_Q",
        ):
            if item[key] != "":
                item[key] = fnum(item[key])
        stratum_out.append(item)
    write_tsv(OUT_TAB / "stratum_summary.tsv", stratum_out, stratum_fields)

    gene_fields = [
        "module",
        "human_symbol",
        "log2FC_lung_baseline",
        "log2FC_T47D",
        "log2FC_MCF7",
        "log2FC_SKOV3",
        "n_cancer_observed",
        "cancer_median_log2FC",
        "lung_matches_thesis",
        "cancer_matches_thesis",
    ]
    gene_out = []
    for row in signed:
        item = dict(row)
        for key in (
            "log2FC_lung_baseline",
            "log2FC_T47D",
            "log2FC_MCF7",
            "log2FC_SKOV3",
            "cancer_median_log2FC",
        ):
            if item[key] != "":
                item[key] = fnum(item[key])
        item["n_cancer_observed"] = str(item["n_cancer_observed"])
        gene_out.append(item)
    write_tsv(OUT_TAB / "gene_signed.tsv", gene_out, gene_fields)

    cross_fields = gene_fields + ["lung_sign", "cancer_sign", "same_sign", "lung_up_cancer_down"]
    cross_out = []
    for row in cross:
        item = {key: row[key] for key in cross_fields}
        for key in (
            "log2FC_lung_baseline",
            "log2FC_T47D",
            "log2FC_MCF7",
            "log2FC_SKOV3",
            "cancer_median_log2FC",
        ):
            if item[key] != "":
                item[key] = fnum(item[key])
        item["n_cancer_observed"] = str(item["n_cancer_observed"])
        cross_out.append(item)
    write_tsv(OUT_TAB / "cross_organ_genes.tsv", cross_out, cross_fields)

    cross_sum_fields = [
        "module",
        "n_genes_in_lung_and_cancer",
        "n_same_sign",
        "n_opposite_sign",
        "frac_same_sign",
        "n_lung_up",
        "n_cancer_median_up",
        "n_both_up",
        "n_lung_up_cancer_down",
        "frac_lung_up_cancer_down",
    ]
    cross_sum_out = []
    for module in MODULES:
        item = dict(cross_summary[module])
        item["frac_same_sign"] = fnum(item["frac_same_sign"])
        item["frac_lung_up_cancer_down"] = fnum(item["frac_lung_up_cancer_down"])
        for key in (
            "n_genes_in_lung_and_cancer",
            "n_same_sign",
            "n_opposite_sign",
            "n_lung_up",
            "n_cancer_median_up",
            "n_both_up",
            "n_lung_up_cancer_down",
        ):
            item[key] = str(item[key])
        cross_sum_out.append(item)
    write_tsv(OUT_TAB / "cross_organ_summary.tsv", cross_sum_out, cross_sum_fields)

    forest_fields = [
        "row_type",
        "accession",
        "contrast",
        "class",
        "organ",
        "module",
        "label",
        "estimate",
        "lo",
        "hi",
        "interval",
        "n_up",
        "n_observed",
        "call_vs_thesis",
        "thesis_sign_match",
    ]
    forest_out = []
    for row in plot_rows:
        item = dict(row)
        for key in ("estimate", "lo", "hi"):
            item[key] = fnum(item[key])
        item["n_up"] = str(item["n_up"])
        item["n_observed"] = str(item["n_observed"])
        forest_out.append(item)
    write_tsv(OUT_TAB / "forest_effects.tsv", forest_out, forest_fields)

    key = {
        "thesis": THESIS,
        "n_contrasts": 4,
        "gse22493_mapped": True,
        "gse50927_contrast": "baseline KO vs WT",
        "lung_ifn_median": lung_median(rows, "IFN_ISG"),
        "lung_apm_median": lung_median(rows, "MHC_APM"),
        "h2k1_log2FC": h2k1,
        "cancer_ifn_calls": next(
            s for s in strata if s["class"] == "cancer" and s["module"] == "IFN_ISG"
        ),
        "cancer_apm_calls": next(
            s for s in strata if s["class"] == "cancer" and s["module"] == "MHC_APM"
        ),
        "cross_organ": cross_summary,
    }
    # JSON-safe strata slices
    for name in ("cancer_ifn_calls", "cancer_apm_calls"):
        key[name] = {
            k: (None if isinstance(v, float) and math.isnan(v) else v)
            for k, v in key[name].items()
            if not k.startswith("_")
        }
    (OUT_TAB / "key_stats.json").write_text(json.dumps(key, indent=2) + "\n")
    assert_story(rows, cross_summary)
    print(json.dumps({"lung_h2k1": h2k1, "cross": cross_summary, "n_rows": len(rows)}, indent=2, default=str))


def lung_median(rows: list[dict], module: str) -> float:
    for row in rows:
        if row["accession"] == "GSE50927" and row["module"] == module:
            return row["median_log2FC"]
    raise RuntimeError(module)


if __name__ == "__main__":
    main()
