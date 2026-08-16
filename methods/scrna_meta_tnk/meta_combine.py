#!/usr/bin/env python3
"""Random-effects / Stouffer / Fisher combination + forest plot."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lib_stats import fisher_combine, fisher_z, fisher_z_var, random_effects_dl, stouffer

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"


def _ci_rho(rho: float, n: int) -> tuple[float, float]:
    z = fisher_z(rho)
    se = np.sqrt(fisher_z_var(n))
    return float(np.tanh(z - 1.96 * se)), float(np.tanh(z + 1.96 * se))


def pool_gene(df: pd.DataFrame, gene: str) -> dict:
    sub = df[(df["gene"] == gene) & (df["primary"] == True) & (df["source"] != "skip")].copy()
    sub = sub[sub["n"] >= 4]
    sub = sub[np.isfinite(sub["rho"])]
    rhos = sub["rho"].tolist()
    ps = sub["p"].tolist()
    ns = sub["n"].astype(int).tolist()
    re = random_effects_dl(rhos, ns)
    st = stouffer(rhos, ps, ns)
    fi = fisher_combine(ps)
    return {
        "gene": gene,
        "n_cohorts": int(len(sub)),
        "n_patients_total": int(sub["n"].sum()),
        "cohorts": sub["cohort"].tolist(),
        "random_effects": re,
        "stouffer": st,
        "fisher": fi,
        "rows": sub.to_dict(orient="records"),
    }


def forest(df: pd.DataFrame, gene: str, pooled: dict, path: Path) -> None:
    sub = df[(df["gene"] == gene) & (df["primary"] == True) & (df["source"] != "skip")].copy()
    sub = sub[sub["n"] >= 4].sort_values("n", ascending=False)
    re = pooled["random_effects"]
    labels = []
    rhos, los, his, ns, ps = [], [], [], [], []
    for _, r in sub.iterrows():
        lo, hi = _ci_rho(float(r["rho"]), int(r["n"]))
        labels.append(f"{r['cohort']}  n={int(r['n'])}")
        rhos.append(float(r["rho"]))
        los.append(lo)
        his.append(hi)
        ns.append(int(r["n"]))
        ps.append(float(r["p"]))
    labels.append(f"RE pooled  N={re['n_patients_total']}")
    rhos.append(re["pooled_rho"])
    los.append(re["ci95_rho"][0])
    his.append(re["ci95_rho"][1])
    ns.append(re["n_patients_total"])
    ps.append(re["p"])

    fig_h = 1.0 + 0.42 * len(labels)
    fig, ax = plt.subplots(figsize=(8.6, fig_h))
    y = np.arange(len(labels))
    for i, (rho, lo, hi) in enumerate(zip(rhos, los, his)):
        color = "#1f4e79" if i == len(labels) - 1 else "#4a4a4a"
        marker = "D" if i == len(labels) - 1 else "o"
        ax.plot([lo, hi], [y[i], y[i]], color=color, lw=1.6, solid_capstyle="round")
        ax.plot(rho, y[i], marker=marker, color=color, ms=7 if i == len(labels) - 1 else 6)
        ax.text(1.02, y[i], f"ρ={rho:+.2f}  p={ps[i]:.3g}", va="center", ha="left", fontsize=8, family="monospace")
    ax.axvline(0, color="#888888", lw=0.8, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel(f"Spearman ρ  ({gene} in malignant/epithelial cells vs T/NK fraction)")
    ax.set_xlim(-1.05, 1.05)
    ax.set_title(
        f"{gene} vs T/NK  ·  RE ρ={re['pooled_rho']:+.3f} "
        f"[{re['ci95_rho'][0]:+.3f}, {re['ci95_rho'][1]:+.3f}]  "
        f"p={re['p']:.3g}  I²={re['I2']:.0f}%  N={re['n_patients_total']}",
        fontsize=10,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    df = pd.read_csv(OUT / "cohort_effects.tsv", sep="\t")
    # pandas may read True as bool already
    if df["primary"].dtype != bool:
        df["primary"] = df["primary"].astype(str).str.lower().isin(["true", "1", "yes"])
    summary = {}
    for gene in ("TACSTD2", "CLDN4"):
        pooled = pool_gene(df, gene)
        summary[gene] = pooled
        forest(df, gene, pooled, OUT / f"forest_{gene}.png")
        print(
            gene,
            "k=",
            pooled["n_cohorts"],
            "N=",
            pooled["n_patients_total"],
            "RE ρ=",
            f"{pooled['random_effects']['pooled_rho']:.3f}",
            "p=",
            f"{pooled['random_effects']['p']:.4g}",
            "I2=",
            f"{pooled['random_effects']['I2']:.1f}",
            "Stouffer z=",
            f"{pooled['stouffer']['z']:.3f}",
            "p=",
            f"{pooled['stouffer']['p']:.4g}",
        )
    (OUT / "meta_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    # compact table
    rows = []
    for gene, pooled in summary.items():
        re = pooled["random_effects"]
        rows.append(
            {
                "gene": gene,
                "n_cohorts": pooled["n_cohorts"],
                "n_patients_total": pooled["n_patients_total"],
                "pooled_rho": re["pooled_rho"],
                "ci95_lo": re["ci95_rho"][0],
                "ci95_hi": re["ci95_rho"][1],
                "p_RE": re["p"],
                "I2": re["I2"],
                "tau2": re["tau2"],
                "stouffer_z": pooled["stouffer"]["z"],
                "stouffer_p": pooled["stouffer"]["p"],
                "fisher_p": pooled["fisher"]["p"],
            }
        )
    pd.DataFrame(rows).to_csv(OUT / "meta_pooled.tsv", sep="\t", index=False)
    print("wrote", OUT / "meta_summary.json")


if __name__ == "__main__":
    main()
