#!/usr/bin/env python3
"""A10: is GRHL1 a public correlate of TACSTD2 (TROP2) and CLDN4 in lung?

Honest test of the user claim that GRHL1 is a shared regulator/correlate of BOTH
TACSTD2 and CLDN4 in lung. We compute Spearman and Pearson correlations of GRHL1
with each gene across TCGA lung cohorts, and report GRHL2/GRHL3 and the two
claudins as context. We do not pick the cohort or method that flatters the claim;
every cohort/method is reported side by side.

Association is a co-expression correlation only. It is consistent with, but does
NOT prove, direct transcriptional regulation.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# Pre-specified honest thresholds for a co-expression "association" (Spearman rho,
# two-sided p). These are fixed before looking at the numbers.
STRONG = 0.50
MODERATE = 0.30
WEAK = 0.10
ALPHA = 0.05


def corr_block(x: np.ndarray, y: np.ndarray) -> dict:
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = int(len(x))
    if n < 8 or np.std(x) == 0 or np.std(y) == 0:
        return {
            "n": n,
            "spearman_rho": None,
            "spearman_p": None,
            "pearson_r": None,
            "pearson_p": None,
        }
    rho, p_s = stats.spearmanr(x, y)
    r, p_p = stats.pearsonr(x, y)
    return {
        "n": n,
        "spearman_rho": float(rho),
        "spearman_p": float(p_s),
        "pearson_r": float(r),
        "pearson_p": float(p_p),
    }


def strength_label(rho: float | None, p: float | None) -> str:
    if rho is None or p is None or not np.isfinite(rho):
        return "NOT_COMPUTED"
    if p > ALPHA:
        return "NOT_SIGNIFICANT"
    a = abs(rho)
    direction = "pos" if rho > 0 else "neg"
    if a >= STRONG:
        mag = "STRONG"
    elif a >= MODERATE:
        mag = "MODERATE"
    elif a >= WEAK:
        mag = "WEAK"
    else:
        mag = "NEGLIGIBLE"
    return f"{mag}_{direction}"


def supports_claim(rho: float | None, p: float | None) -> bool:
    """A single pair 'supports' the claim if it is a positive, significant
    association of at least MODERATE magnitude."""
    if rho is None or p is None or not np.isfinite(rho):
        return False
    return (rho >= MODERATE) and (p <= ALPHA)


def pair(df: pd.DataFrame, g1: str, g2: str) -> dict:
    b = corr_block(df[g1].to_numpy(float), df[g2].to_numpy(float))
    return {
        "gene_a": g1,
        "gene_b": g2,
        **b,
        "spearman_strength": strength_label(b["spearman_rho"], b["spearman_p"]),
    }


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--in-dir", default="results/w200/A10_GRHL1")
    p.add_argument("--out-dir", default="results/w200/A10_GRHL1")
    args = p.parse_args()
    indir = Path(args.in_dir)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(indir / "expression_tcga_lung.csv")
    for c in df.columns:
        if c not in ("sample", "cohort", "sample_type"):
            df[c] = pd.to_numeric(df[c], errors="coerce")

    # Cohort slices. PRIMARY = NSCLC primary tumors (LUAD+LUSC), the lung-tumor
    # setting the claim is about. Per-cohort tumor slices avoid LUAD/LUSC mixing
    # artifacts; normal lung is reported as context.
    tumor = df[df["sample_type"] == "primary_tumor"]
    slices = {
        "NSCLC_primary_tumor": tumor,
        "LUAD_primary_tumor": tumor[tumor["cohort"] == "LUAD"],
        "LUSC_primary_tumor": tumor[tumor["cohort"] == "LUSC"],
        "LUAD_normal": df[(df["cohort"] == "LUAD") & (df["sample_type"] == "normal")],
        "LUSC_normal": df[(df["cohort"] == "LUSC") & (df["sample_type"] == "normal")],
        "all_lung_normal": df[df["sample_type"] == "normal"],
    }
    PRIMARY = "NSCLC_primary_tumor"

    # Pairs we report. First two are the claim; the rest are honest context.
    claim_pairs = [("GRHL1", "TACSTD2"), ("GRHL1", "CLDN4")]
    context_pairs = [
        ("GRHL2", "TACSTD2"),
        ("GRHL2", "CLDN4"),
        ("GRHL3", "TACSTD2"),
        ("GRHL3", "CLDN4"),
        ("TACSTD2", "CLDN4"),
        ("GRHL1", "GRHL2"),
        ("GRHL1", "CLDN3"),
        ("GRHL1", "CLDN7"),
    ]

    rows: list[dict] = []
    per_slice: dict[str, dict] = {}
    for name, sub in slices.items():
        block = {"cohort": name, "n_samples": int(len(sub)), "pairs": []}
        for g1, g2 in claim_pairs + context_pairs:
            if g1 not in sub.columns or g2 not in sub.columns:
                continue
            rec = pair(sub, g1, g2)
            rec["cohort"] = name
            rec["is_claim_pair"] = (g1, g2) in claim_pairs
            rows.append(rec)
            block["pairs"].append(rec)
        per_slice[name] = block

    flat = pd.DataFrame(rows)[
        [
            "cohort",
            "gene_a",
            "gene_b",
            "is_claim_pair",
            "n",
            "spearman_rho",
            "spearman_p",
            "pearson_r",
            "pearson_p",
            "spearman_strength",
        ]
    ]
    flat.to_csv(out / "correlations.csv", index=False)

    # Honest verdict on the two claim pairs in the primary cohort.
    def get(cohort: str, g1: str, g2: str) -> dict:
        for r in per_slice[cohort]["pairs"]:
            if r["gene_a"] == g1 and r["gene_b"] == g2:
                return r
        return {}

    g1_trop2 = get(PRIMARY, "GRHL1", "TACSTD2")
    g1_cldn4 = get(PRIMARY, "GRHL1", "CLDN4")

    trop2_ok = supports_claim(g1_trop2.get("spearman_rho"), g1_trop2.get("spearman_p"))
    cldn4_ok = supports_claim(g1_cldn4.get("spearman_rho"), g1_cldn4.get("spearman_p"))

    if trop2_ok and cldn4_ok:
        verdict_label = "SUPPORTED_FOR_BOTH"
    elif trop2_ok or cldn4_ok:
        verdict_label = "PARTIAL_ONLY_ONE_TARGET"
    else:
        verdict_label = "NOT_SUPPORTED"

    luad_trop2 = get("LUAD_primary_tumor", "GRHL1", "TACSTD2")
    luad_cldn4 = get("LUAD_primary_tumor", "GRHL1", "CLDN4")
    lusc_trop2 = get("LUSC_primary_tumor", "GRHL1", "TACSTD2")
    lusc_cldn4 = get("LUSC_primary_tumor", "GRHL1", "CLDN4")

    statement = _statement(
        g1_trop2,
        g1_cldn4,
        trop2_ok,
        cldn4_ok,
        luad_trop2,
        luad_cldn4,
        lusc_trop2,
        lusc_cldn4,
        PRIMARY,
    )

    summary = {
        "task": "A10_GRHL1",
        "question": (
            "Is GRHL1 a public co-expression correlate of BOTH TACSTD2 (TROP2) "
            "and CLDN4 in lung?"
        ),
        "data": "TCGA LUAD+LUSC, UCSC Xena HiSeqV2 log2(norm_count+1)",
        "primary_cohort": PRIMARY,
        "thresholds": {
            "moderate_rho": MODERATE,
            "strong_rho": STRONG,
            "alpha": ALPHA,
            "support_rule": "positive Spearman rho >= 0.30 and p <= 0.05",
        },
        "claim_pairs_primary": {
            "GRHL1_vs_TACSTD2": g1_trop2,
            "GRHL1_vs_CLDN4": g1_cldn4,
        },
        "claim_pairs_by_histology": {
            "LUAD_primary_tumor": {
                "GRHL1_vs_TACSTD2": luad_trop2,
                "GRHL1_vs_CLDN4": luad_cldn4,
                "both_supported": supports_claim(
                    luad_trop2.get("spearman_rho"), luad_trop2.get("spearman_p")
                )
                and supports_claim(
                    luad_cldn4.get("spearman_rho"), luad_cldn4.get("spearman_p")
                ),
            },
            "LUSC_primary_tumor": {
                "GRHL1_vs_TACSTD2": lusc_trop2,
                "GRHL1_vs_CLDN4": lusc_cldn4,
                "both_supported": supports_claim(
                    lusc_trop2.get("spearman_rho"), lusc_trop2.get("spearman_p")
                )
                and supports_claim(
                    lusc_cldn4.get("spearman_rho"), lusc_cldn4.get("spearman_p")
                ),
            },
        },
        "honest_verdict": {
            "label": verdict_label,
            "applied_to": PRIMARY,
            "GRHL1_TACSTD2_supported": trop2_ok,
            "GRHL1_CLDN4_supported": cldn4_ok,
            "statement": statement,
        },
        "per_cohort": per_slice,
        "caveats": [
            "Co-expression correlation only; not evidence of direct transcriptional regulation.",
            "Primary verdict is mixed NSCLC. LUAD and LUSC are pre-specified sensitivities, not a search for a nicer number.",
            "Pooling LUAD+LUSC can cancel a within-histology GRHL1–CLDN4 correlation (different means / ranks across histologies).",
            "GRHL2/GRHL3 rows are context. They are not used to rescue or replace the GRHL1 claim.",
            "No cohort or correlation method was chosen to favor the claim.",
        ],
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    _figures(slices, PRIMARY, out, g1_trop2, g1_cldn4, per_slice)

    print(json.dumps(summary["honest_verdict"], indent=2))
    print()
    print(flat.to_string(index=False))
    return 0


def _statement(
    trop2,
    cldn4,
    trop2_ok,
    cldn4_ok,
    luad_trop2,
    luad_cldn4,
    lusc_trop2,
    lusc_cldn4,
    primary,
) -> str:
    def fmt(rec):
        if not rec or rec.get("spearman_rho") is None:
            return "not computed"
        return (
            f"rho={rec['spearman_rho']:.3f} (p={rec['spearman_p']:.2e}, "
            f"n={rec['n']}, {rec['spearman_strength']})"
        )

    parts = [
        f"In {primary}, GRHL1 vs TACSTD2: {fmt(trop2)}; "
        f"GRHL1 vs CLDN4: {fmt(cldn4)}."
    ]
    if trop2_ok and cldn4_ok:
        parts.append(
            "Both associations are positive, significant, and at least moderate, "
            "so the claim that GRHL1 tracks both TROP2 and CLDN4 in mixed NSCLC is supported."
        )
    elif trop2_ok and not cldn4_ok:
        parts.append(
            "On the pre-specified mixed-NSCLC slice, GRHL1 tracks TACSTD2 (TROP2) "
            "but not CLDN4 (near-zero, not significant). The shared-target claim is "
            "only half true here."
        )
    elif cldn4_ok and not trop2_ok:
        parts.append(
            "GRHL1 tracks CLDN4 but not TACSTD2 at the moderate+significant bar; "
            "the claim is only partially supported on mixed NSCLC."
        )
    else:
        parts.append(
            "Neither association clears the moderate+significant bar on mixed NSCLC, "
            "so this public slice does not support GRHL1 as a shared correlate of "
            "TROP2 and CLDN4."
        )
    parts.append(
        f"Histology split (pre-specified sensitivity, not a rescue): "
        f"LUAD GRHL1–TACSTD2 {fmt(luad_trop2)}; LUAD GRHL1–CLDN4 {fmt(luad_cldn4)}; "
        f"LUSC GRHL1–TACSTD2 {fmt(lusc_trop2)}; LUSC GRHL1–CLDN4 {fmt(lusc_cldn4)}. "
        "Pooling LUAD+LUSC can cancel a within-LUAD GRHL1–CLDN4 correlation. "
        "We still score the claim on mixed NSCLC because that was the primary cohort."
    )
    return " ".join(parts)


def _figures(slices, primary, out, trop2, cldn4, per_slice) -> None:
    sub = slices[primary]
    colors = {"LUAD": "#1f77b4", "LUSC": "#d62728"}
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    for ax, gy, rec in (
        (axes[0], "TACSTD2", trop2),
        (axes[1], "CLDN4", cldn4),
    ):
        for hist, g in sub.groupby("cohort"):
            ax.scatter(
                g["GRHL1"],
                g[gy],
                s=14,
                alpha=0.55,
                c=colors.get(hist, "#7f7f7f"),
                edgecolors="none",
                label=f"{hist} (n={len(g)})",
            )
        rho = rec.get("spearman_rho")
        title = f"GRHL1 vs {gy}\n"
        title += (
            f"mixed NSCLC Spearman rho = {rho:.3f} (n={rec.get('n')})"
            if rho is not None
            else "n/a"
        )
        ax.set_title(title)
        ax.set_xlabel("GRHL1  log2(norm_count+1)")
        ax.set_ylabel(f"{gy}  log2(norm_count+1)")
        ax.legend(loc="best", fontsize=8, frameon=False)
        ax.grid(True, alpha=0.25)
    fig.suptitle("TCGA NSCLC primary tumors: GRHL1 vs TROP2 / CLDN4")
    fig.tight_layout()
    fig.savefig(out / "fig_grhl1_vs_targets_nsclc.png", dpi=160)
    plt.close(fig)

    # Bar of the two claim pairs across the three tumor slices.
    cohorts = ["NSCLC_primary_tumor", "LUAD_primary_tumor", "LUSC_primary_tumor"]
    labels = ["NSCLC mixed", "LUAD", "LUSC"]
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    x = np.arange(len(cohorts))
    width = 0.36
    trop2_vals = []
    cldn4_vals = []
    for c in cohorts:
        t = next(r for r in per_slice[c]["pairs"] if r["gene_a"] == "GRHL1" and r["gene_b"] == "TACSTD2")
        d = next(r for r in per_slice[c]["pairs"] if r["gene_a"] == "GRHL1" and r["gene_b"] == "CLDN4")
        trop2_vals.append(t["spearman_rho"])
        cldn4_vals.append(d["spearman_rho"])
    ax.bar(x - width / 2, trop2_vals, width, color="#1f77b4", label="GRHL1–TACSTD2")
    ax.bar(x + width / 2, cldn4_vals, width, color="#d62728", label="GRHL1–CLDN4")
    ax.axhline(0.30, color="0.35", ls="--", lw=1, label="support bar (ρ=0.30)")
    ax.axhline(0, color="0.5", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Spearman ρ")
    ax.set_ylim(-0.15, 0.65)
    ax.set_title("GRHL1 vs TACSTD2 / CLDN4 by TCGA lung slice")
    ax.legend(loc="upper right", fontsize=8, frameon=False)
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out / "fig_grhl1_rho_by_histology.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
