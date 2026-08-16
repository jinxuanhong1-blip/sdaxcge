#!/usr/bin/env python3
"""Within-histology Fisher-z meta and forest plots. One gene × histology at a time."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lib import fisher_z, fisher_z_se, meta_fisher_z, tanh_rho


ROOT = Path(__file__).resolve().parents[1]


def _eligible(effects: pd.DataFrame, gene: str, hist: str, min_n: int, primary_only: bool) -> pd.DataFrame:
    sub = effects[(effects["gene"] == gene) & (effects["histology"] == hist)].copy()
    # Do not double-count overlapping GSE131907 site subsets or IIT+RWC combined.
    sub = sub[~sub["cohort"].astype(str).str.contains("combined|tumor_sites")]
    if primary_only:
        sub = sub[sub["primary_meta"] == 1]
    else:
        sub = sub[sub["n"] >= min_n]
        sub = sub[sub["rho"].notna()]
    return sub


def run_meta(effects: pd.DataFrame, min_n: int, primary_only: bool) -> pd.DataFrame:
    rows = []
    for gene in ("TACSTD2", "CLDN4"):
        for hist in ("LUAD", "LUSC"):
            sub = _eligible(effects, gene, hist, min_n, primary_only)
            meta = meta_fisher_z(sub.to_dict("records"), min_n=min_n)
            rows.append(
                {
                    "gene": gene,
                    "histology": hist,
                    "min_n": min_n,
                    "k": meta["k"],
                    "n_total": meta["n_total"],
                    "rho_re": meta["rho_re"],
                    "ci_re_low": meta["ci_re_low"],
                    "ci_re_high": meta["ci_re_high"],
                    "p_re": meta["p_re"],
                    "rho_fe": meta["rho_fe"],
                    "ci_fe_low": meta["ci_fe_low"],
                    "ci_fe_high": meta["ci_fe_high"],
                    "p_fe": meta["p_fe"],
                    "tau2": meta["tau2"],
                    "I2": meta["I2"],
                    "Q": meta["Q"],
                    "Q_p": meta["Q_p"],
                    "cohorts": ";".join(sub["cohort"].tolist()),
                }
            )
    return pd.DataFrame(rows)


def forest(effects: pd.DataFrame, meta: pd.DataFrame, gene: str, out: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 6.4), sharex=True)
    for ax, hist, color in zip(axes, ("LUAD", "LUSC"), ("#2c7fb8", "#d95f0e")):
        sub = _eligible(effects, gene, hist, min_n=4, primary_only=False)
        m = meta[(meta["gene"] == gene) & (meta["histology"] == hist)].iloc[0]
        labels = []
        rhos = []
        los = []
        his = []
        pooled = []
        for _, r in sub.iterrows():
            tag = "" if int(r["primary_meta"]) == 1 else "  [n<6, not pooled]"
            labels.append(f"{r['cohort']}  n={int(r['n'])}{tag}")
            rhos.append(r["rho"])
            los.append(r["ci_low"])
            his.append(r["ci_high"])
            pooled.append(int(r["primary_meta"]) == 1)
        if m["k"] >= 1 and np.isfinite(m["rho_re"]):
            labels.append(f"RE meta  k={int(m['k'])} n={int(m['n_total'])}")
            rhos.append(m["rho_re"])
            los.append(m["ci_re_low"])
            his.append(m["ci_re_high"])
        y = np.arange(len(labels))[::-1]
        for i, (lab, rho, lo, hi) in enumerate(zip(labels, rhos, los, his)):
            is_meta = lab.startswith("RE meta")
            is_pooled = (not is_meta) and (i < len(pooled) and pooled[i])
            ax.plot([lo, hi], [y[i], y[i]], color="black" if is_meta else color, lw=1.6, alpha=1.0 if (is_meta or is_pooled) else 0.55)
            ax.plot(
                rho,
                y[i],
                "D" if is_meta else ("o" if is_pooled else "o"),
                color="black" if is_meta else color,
                mfc=("black" if is_meta else (color if is_pooled else "white")),
                ms=7 if is_meta else 6,
            )
        ax.axvline(0, color="#666666", lw=0.8, ls="--")
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_title(f"{hist}  {gene} vs T/NK", fontsize=11)
        ax.set_xlabel("Spearman ρ (Fisher-z 95% CI)")
        ax.set_xlim(-1.05, 1.05)
        if m["k"] >= 1 and np.isfinite(m["rho_re"]):
            ax.text(
                0.02,
                0.02,
                f"RE ρ={m['rho_re']:+.2f}  p={m['p_re']:.3g}  I²={m['I2']:.0f}%",
                transform=ax.transAxes,
                fontsize=8,
                va="bottom",
            )
        else:
            ax.text(0.02, 0.02, "no primary cohort with n≥6", transform=ax.transAxes, fontsize=8)
    fig.suptitle(
        f"Public lung tumor scRNA · malignant/epithelial {gene} vs T/NK fraction\n"
        "split by author histology · sample/patient unit · not pooled across histology",
        fontsize=11,
    )
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=160)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--effects", default=str(ROOT / "results" / "cohort_effects.tsv"))
    ap.add_argument("--out-dir", default=str(ROOT / "results"))
    args = ap.parse_args()
    out_dir = Path(args.out_dir)
    effects = pd.read_csv(args.effects, sep="\t")
    meta6 = run_meta(effects, min_n=6, primary_only=True)
    meta4 = run_meta(effects, min_n=4, primary_only=False)
    meta6.to_csv(out_dir / "meta_by_histology.tsv", sep="\t", index=False)
    meta4.to_csv(out_dir / "meta_by_histology_n4_sensitivity.tsv", sep="\t", index=False)
    fig_dir = out_dir / "figures"
    forest(effects, meta6, "TACSTD2", fig_dir / "forest_TACSTD2_by_histology.png")
    forest(effects, meta6, "CLDN4", fig_dir / "forest_CLDN4_by_histology.png")
    (out_dir / "meta_summary.json").write_text(
        json.dumps(
            {
                "primary_min_n": 6,
                "meta": meta6.to_dict("records"),
            },
            indent=2,
        )
        + "\n"
    )
    print(meta6.to_string(index=False))


if __name__ == "__main__":
    main()
