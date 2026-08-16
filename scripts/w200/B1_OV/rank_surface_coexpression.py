#!/usr/bin/env python3
"""TCGA-OV: honest TACSTD2–CLDN4 rank among the cell-surface proteome.

Continuation of the B1_OV slice. Mirrors the B1_BRCA analog
(scripts/coexpression.py on cursor/tacstd2-cldn4-coexpression-brca-c503):
rank every surfaceome gene by Spearman co-expression with TACSTD2 and
report where CLDN4 actually lands. No pre-selected neighbourhood.

Definitions
-----------
* Surface-gene universe: the 2,886-protein "in silico human surfaceome"
  (Bausch-Fluck et al., PNAS 2018, table S3), intersected with genes
  measured in the TCGA-OV STAR-TPM matrix, with the anchor (TACSTD2)
  itself removed.
* Samples: TCGA-OV primary tumours (barcode sample-type ``01``), one
  sample per patient — same n as the rest of B1_OV.
* Expression: Xena GDC hub STAR TPM, log2(TPM+1), GENCODE v36 symbols.
* Primary metric: Spearman rho. Pearson is reported as a sensitivity.

Outputs land in results/w200/B1_OV/ (the w200 window is the top-200
companion file; the full ranking always covers the entire universe).

Usage:  python3 scripts/w200/B1_OV/rank_surface_coexpression.py
"""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA_DIR = os.environ.get("W200_B1_OV_DATA", "/tmp/w200_b1_ov_data")
REPO = os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
)
OUT = os.path.join(REPO, "results", "w200", "B1_OV")
os.makedirs(OUT, exist_ok=True)

ANCHOR = "TACSTD2"
FOCUS = "CLDN4"
WINDOW = 200


def load_expression() -> pd.DataFrame:
    """log2(TPM+1) DataFrame, rows = HGNC symbols, cols = samples."""
    probemap = pd.read_csv(
        os.path.join(DATA_DIR, "gencode.v36.probemap"), sep="\t", usecols=["id", "gene"]
    )
    id2sym = dict(zip(probemap["id"], probemap["gene"]))
    path = os.path.join(DATA_DIR, "TCGA-OV.star_tpm.tsv.gz")
    expr = pd.read_csv(path, sep="\t", index_col=0)
    expr.index = expr.index.map(lambda i: id2sym.get(i, i))
    # Collapse duplicate symbols (PAR / multi-mapped GENCODE IDs) by keeping
    # the highest-mean row. loc[symbol] cannot do this: it re-selects every
    # row sharing that symbol.
    expr = (
        expr.assign(_mean=expr.mean(axis=1))
        .sort_values("_mean", ascending=False)
        .drop(columns="_mean")
    )
    expr = expr[~expr.index.duplicated(keep="first")]
    return expr


def load_surfaceome() -> pd.DataFrame:
    path = os.path.join(DATA_DIR, "table_S3_surfaceome.xlsx")
    surf = pd.read_excel(
        path, sheet_name="in silico surfaceome only", engine="openpyxl", header=1
    )
    surf = surf.rename(columns={"UniProt gene": "gene"})
    surf = surf[surf["gene"].notna()].copy()
    surf["gene"] = surf["gene"].astype(str)
    return surf


def primary_tumors_one_per_patient(columns) -> list[str]:
    tumors = {}
    for s in sorted(columns):
        parts = s.split("-")
        if len(parts) < 4:
            continue
        code = parts[3][:2]
        patient = "-".join(parts[:3])
        if code == "01" and patient not in tumors:
            tumors[patient] = s
    return list(tumors.values())


def _corr_matrix_vs_vector(mat: np.ndarray, vec: np.ndarray) -> np.ndarray:
    vc = vec - vec.mean()
    mc = mat - mat.mean(axis=1, keepdims=True)
    num = mc @ vc
    den = np.sqrt((mc**2).sum(axis=1) * (vc**2).sum())
    with np.errstate(invalid="ignore", divide="ignore"):
        return num / den


def _pvalues(r: np.ndarray, n: int) -> np.ndarray:
    r = np.clip(r, -0.999999999, 0.999999999)
    with np.errstate(invalid="ignore", divide="ignore"):
        t = r * np.sqrt((n - 2) / (1 - r**2))
    return 2 * stats.t.sf(np.abs(t), df=n - 2)


def _bh_fdr(p: np.ndarray) -> np.ndarray:
    q = np.full_like(p, np.nan, dtype=float)
    ok = ~np.isnan(p)
    pv = p[ok]
    m = pv.size
    order = np.argsort(pv)
    ranked = pv[order]
    adj = ranked * m / (np.arange(m) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.clip(adj, 0, 1)
    q[ok] = out
    return q


def compute_ranking(expr: pd.DataFrame, universe: list[str], samples: list[str]) -> pd.DataFrame:
    sub = expr.loc[universe, samples]
    anchor_vec = expr.loc[ANCHOR, samples].to_numpy(dtype=float)
    mat = sub.to_numpy(dtype=float)
    n = len(samples)

    pear_r = _corr_matrix_vs_vector(mat, anchor_vec)
    anchor_rank = stats.rankdata(anchor_vec)
    mat_rank = np.apply_along_axis(stats.rankdata, 1, mat)
    spear_r = _corr_matrix_vs_vector(mat_rank, anchor_rank)

    df = pd.DataFrame(
        {
            "gene": universe,
            "spearman_r": spear_r,
            "spearman_p": _pvalues(spear_r, n),
            "pearson_r": pear_r,
            "pearson_p": _pvalues(pear_r, n),
            "mean_log2tpm": mat.mean(axis=1),
            "pct_expressed": (mat > 0).mean(axis=1) * 100.0,
        }
    )
    df["spearman_q"] = _bh_fdr(df["spearman_p"].to_numpy())
    df["pearson_q"] = _bh_fdr(df["pearson_p"].to_numpy())
    df = df.sort_values("spearman_r", ascending=False, kind="mergesort").reset_index(drop=True)
    df.insert(1, "spearman_rank", np.arange(1, len(df) + 1))
    df["pearson_rank"] = df["pearson_r"].rank(ascending=False, method="min").astype(int)
    return df


def main() -> int:
    expr = load_expression()
    if ANCHOR not in expr.index:
        raise SystemExit(f"anchor {ANCHOR} not in expression matrix")
    if FOCUS not in expr.index:
        raise SystemExit(f"focus {FOCUS} not in expression matrix")

    surf = load_surfaceome()
    surf_genes = sorted({g for g in surf["gene"] if g not in ("nan", "None", "")})
    n_surfaceome_table = len(surf_genes)

    samples = primary_tumors_one_per_patient(expr.columns)
    n = len(samples)

    present = [g for g in surf_genes if g in expr.index]
    missing_from_matrix = sorted(set(surf_genes) - set(present))
    candidates = [g for g in present if g != ANCHOR]
    stds = expr.loc[candidates, samples].std(axis=1)
    zero_var = sorted(stds.index[stds == 0].tolist())
    universe = [g for g in candidates if g not in set(zero_var)]
    n_universe = len(universe)
    if zero_var:
        print(f"[warn] dropped {len(zero_var)} zero-variance surface genes (rho undefined)")
    if FOCUS not in universe:
        raise SystemExit(
            f"{FOCUS} is not in the surfaceome ∩ matrix universe "
            f"(present={FOCUS in present}, in_table={FOCUS in surf_genes})"
        )

    df = compute_ranking(expr, universe, samples)
    df.to_csv(os.path.join(OUT, "coexpression_TACSTD2_surfaceome.csv"), index=False)
    df.head(WINDOW).to_csv(os.path.join(OUT, "top200.csv"), index=False)

    row = df[df["gene"] == FOCUS].iloc[0]
    s_rank = int(row["spearman_rank"])
    p_rank = int(row["pearson_rank"])
    report = {
        "gene": FOCUS,
        "in_universe": True,
        "is_top_spearman": s_rank == 1,
        "is_top_pearson": p_rank == 1,
        "n_samples": n,
        "n_surface_universe": n_universe,
        "spearman_rank": s_rank,
        "spearman_rho": float(row["spearman_r"]),
        "spearman_p": float(row["spearman_p"]),
        "spearman_q": float(row["spearman_q"]),
        "spearman_percentile": round(100 * (1 - (s_rank - 1) / n_universe), 3),
        "pearson_rank": p_rank,
        "pearson_r": float(row["pearson_r"]),
        "pearson_q": float(row["pearson_q"]),
    }

    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "cohort": "TCGA-OV",
        "expression_dataset": "Xena GDC hub TCGA-OV.star_tpm (log2(TPM+1), GENCODE v36)",
        "surfaceome_source": "Bausch-Fluck et al. 2018 in-silico surfaceome (table S3)",
        "anchor_gene": ANCHOR,
        "focus_gene": FOCUS,
        "sample_rule": "primary tumor (-01), one sample per patient",
        "n_samples": n,
        "n_surfaceome_table_unique_genes": n_surfaceome_table,
        "n_surfaceome_in_matrix": len(present),
        "n_surfaceome_missing_from_matrix": len(missing_from_matrix),
        "n_zero_variance_dropped": len(zero_var),
        "zero_variance_genes": zero_var,
        "surface_gene_universe": n_universe,
        "correlation_primary": "spearman",
        "window": WINDOW,
        "focus_result": report,
        "top10_spearman": df.head(10)[
            ["gene", "spearman_rank", "spearman_r", "spearman_q"]
        ].to_dict("records"),
        "surfaceome_genes_missing_from_matrix_head": missing_from_matrix[:40],
    }
    with open(os.path.join(OUT, "surface_rank_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    # ---- figures -----------------------------------------------------------
    plt.rcParams.update({"font.size": 9, "figure.dpi": 150})

    x = expr.loc[ANCHOR, samples].to_numpy(dtype=float)
    y = expr.loc[FOCUS, samples].to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(5.2, 5.0))
    ax.scatter(x, y, s=8, alpha=0.4, edgecolors="none", color="darkslategray")
    ax.set_xlabel(f"{ANCHOR} log2(TPM+1)")
    ax.set_ylabel(f"{FOCUS} log2(TPM+1)")
    ax.set_title(
        f"{ANCHOR} vs {FOCUS} (TCGA-OV, n={n})\n"
        f"Spearman ρ={row['spearman_r']:.3f}, rank #{s_rank} of {n_universe}"
    )
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig6_scatter_tacstd2_cldn4.png"), bbox_inches="tight")
    plt.close(fig)

    topn = df.head(20).copy()
    if FOCUS not in set(topn["gene"]):
        topn = pd.concat([topn, df[df["gene"] == FOCUS]], ignore_index=True)
    topn = topn.iloc[::-1]
    colors = ["#d62728" if g == FOCUS else "#4c78a8" for g in topn["gene"]]
    labels = [
        f"{g}  (#{int(r)})" if g == FOCUS else g
        for g, r in zip(topn["gene"], topn["spearman_rank"])
    ]
    fig, ax = plt.subplots(figsize=(6.2, 7.4))
    ax.barh(labels, topn["spearman_r"], color=colors)
    ax.set_xlabel(f"Spearman ρ with {ANCHOR}")
    ax.set_title(
        f"Top 20 surface-gene partners of {ANCHOR} in TCGA-OV\n"
        f"{FOCUS} (red) is rank #{s_rank} of {n_universe}, not in the top 20"
    )
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig7_top_surface_partners.png"), bbox_inches="tight")
    plt.close(fig)

    print(json.dumps(report, indent=2))
    print(f"top1={df.iloc[0]['gene']} rho={df.iloc[0]['spearman_r']:.4f}")
    print("done ->", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
