"""Honest co-expression ranking of a hub gene against the cell-surface proteome.

This is the "B1 analog" for TCGA-CESC: it reproduces, for the cervical-cancer
cohort, the panel-B1 style analysis of asking how a candidate surface marker
(default CLDN4) ranks among *all* surface genes when ordered by co-expression
with an anchor gene (default TACSTD2 / TROP2).

The point of the exercise is an *honest* ranking: we score every surface gene
in the surfaceome against the anchor and report where the focus gene actually
lands -- we do not pre-select a favourable neighbourhood or method. Both
Spearman (primary) and Pearson correlations are reported, and the focus gene's
rank is reported against the full surface-gene universe.

Definitions
-----------
* Surface-gene universe: the 2,886-protein "in silico human surfaceome"
  (Bausch-Fluck et al., PNAS 2018), intersected with genes measured in the
  expression matrix, with the anchor gene itself removed.
* Samples: TCGA-CESC primary tumours (barcode sample-type code ``01``) by
  default; configurable via --sample-types.
* Co-expression: correlation of the anchor gene's expression vector with each
  surface gene's expression vector across the selected samples.

"w200" convention
------------------
The output directory prefix ``w200`` denotes the reporting *window*: the number
of top-ranked partners written to the ``top<window>.csv`` companion file
(default 200). It does not alter the full ranking, which always covers the
entire surface-gene universe.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data"
DEFAULT_EXPR = DATA_DIR / "TCGA-CESC.HiSeqV2.gz"
DEFAULT_SURF = DATA_DIR / "table_S3_surfaceome.xlsx"


def load_expression(path: Path) -> pd.DataFrame:
    """Load the Xena matrix: rows = gene symbols, cols = TCGA sample barcodes."""
    expr = pd.read_csv(path, sep="\t", index_col=0)
    expr.index = expr.index.astype(str)
    # Collapse duplicate symbols (rare) by mean so lookups are unambiguous.
    if expr.index.duplicated().any():
        expr = expr.groupby(level=0).mean()
    return expr


def load_surfaceome(path: Path) -> pd.DataFrame:
    """Return the in-silico surfaceome table (header is on the 2nd row)."""
    surf = pd.read_excel(
        path, sheet_name="in silico surfaceome only", engine="openpyxl", header=1
    )
    surf = surf.rename(columns={"UniProt gene": "gene"})
    surf = surf[surf["gene"].notna()].copy()
    surf["gene"] = surf["gene"].astype(str)
    return surf


def select_samples(expr: pd.DataFrame, sample_types: list[str]) -> list[str]:
    """Keep columns whose TCGA sample-type code is in ``sample_types``."""
    keep = []
    for col in expr.columns:
        parts = col.split("-")
        code = parts[3][:2] if len(parts) >= 4 else ""
        if code in sample_types:
            keep.append(col)
    return keep


def _corr_matrix_vs_vector(mat: np.ndarray, vec: np.ndarray) -> np.ndarray:
    """Pearson correlation of each row of ``mat`` against ``vec`` (vectorised)."""
    vc = vec - vec.mean()
    mc = mat - mat.mean(axis=1, keepdims=True)
    num = mc @ vc
    den = np.sqrt((mc**2).sum(axis=1) * (vc**2).sum())
    with np.errstate(invalid="ignore", divide="ignore"):
        return num / den


def _pvalues(r: np.ndarray, n: int) -> np.ndarray:
    """Two-sided p-values for correlation coefficients (t approximation)."""
    r = np.clip(r, -0.999999999, 0.999999999)
    with np.errstate(invalid="ignore", divide="ignore"):
        t = r * np.sqrt((n - 2) / (1 - r**2))
    return 2 * stats.t.sf(np.abs(t), df=n - 2)


def _bh_fdr(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR, NaN-safe."""
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


def compute_ranking(
    expr: pd.DataFrame,
    anchor: str,
    universe: list[str],
    samples: list[str],
) -> pd.DataFrame:
    sub = expr.loc[universe, samples]
    anchor_vec = expr.loc[anchor, samples].to_numpy(dtype=float)
    mat = sub.to_numpy(dtype=float)
    n = len(samples)

    # Pearson (linear).
    pear_r = _corr_matrix_vs_vector(mat, anchor_vec)
    # Spearman (rank-based) = Pearson on rank-transformed data.
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
            "mean_log2_expr": mat.mean(axis=1),
            "pct_expressed": (mat > 0).mean(axis=1) * 100.0,
        }
    )
    df["spearman_q"] = _bh_fdr(df["spearman_p"].to_numpy())
    df["pearson_q"] = _bh_fdr(df["pearson_p"].to_numpy())

    df = df.sort_values("spearman_r", ascending=False, kind="mergesort").reset_index(drop=True)
    df.insert(1, "spearman_rank", np.arange(1, len(df) + 1))
    df["pearson_rank"] = df["pearson_r"].rank(ascending=False, method="min").astype(int)
    return df


def focus_report(df: pd.DataFrame, focus: str, n_universe: int) -> dict:
    row = df[df["gene"] == focus]
    if row.empty:
        return {"gene": focus, "in_universe": False}
    r = row.iloc[0]
    s_rank = int(r["spearman_rank"])
    p_rank = int(r["pearson_rank"])
    return {
        "gene": focus,
        "in_universe": True,
        "is_top_spearman": s_rank == 1,
        "is_top_pearson": p_rank == 1,
        "spearman_rank": s_rank,
        "spearman_r": float(r["spearman_r"]),
        "spearman_q": float(r["spearman_q"]),
        "spearman_percentile": round(100 * (1 - (s_rank - 1) / n_universe), 3),
        "pearson_rank": p_rank,
        "pearson_r": float(r["pearson_r"]),
        "pearson_q": float(r["pearson_q"]),
        "pearson_percentile": round(100 * (1 - (p_rank - 1) / n_universe), 3),
    }


def make_plots(
    expr: pd.DataFrame,
    df: pd.DataFrame,
    anchor: str,
    focus: str,
    samples: list[str],
    window: int,
    outdir: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # 1) Anchor vs focus scatter.
    if focus in expr.index:
        x = expr.loc[anchor, samples].to_numpy(dtype=float)
        y = expr.loc[focus, samples].to_numpy(dtype=float)
        row = df[df["gene"] == focus].iloc[0]
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.scatter(x, y, s=8, alpha=0.4, edgecolors="none")
        ax.set_xlabel(f"{anchor} log2(norm_count+1)")
        ax.set_ylabel(f"{focus} log2(norm_count+1)")
        ax.set_title(
            f"{anchor} vs {focus} (TCGA-CESC, n={len(samples)})\n"
            f"Spearman rho={row['spearman_r']:.3f}, rank #{int(row['spearman_rank'])}"
        )
        fig.tight_layout()
        fig.savefig(outdir / f"scatter_{anchor}_vs_{focus}.png", dpi=150)
        plt.close(fig)

    # 2) Top-N bar chart, focus highlighted.
    topn = df.head(min(20, window)).iloc[::-1]
    colors = ["#d62728" if g == focus else "#4c78a8" for g in topn["gene"]]
    fig, ax = plt.subplots(figsize=(6, 7))
    ax.barh(topn["gene"], topn["spearman_r"], color=colors)
    ax.set_xlabel(f"Spearman rho with {anchor}")
    ax.set_title(f"Top surface-gene co-expression with {anchor}\nTCGA-CESC (focus: {focus})")
    fig.tight_layout()
    fig.savefig(outdir / f"top_partners_{anchor}.png", dpi=150)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--anchor", default="TACSTD2", help="anchor / hub gene (default: TACSTD2)")
    ap.add_argument("--focus", default="CLDN4", help="gene whose rank we report (default: CLDN4)")
    ap.add_argument("--expression", type=Path, default=DEFAULT_EXPR)
    ap.add_argument("--surfaceome", type=Path, default=DEFAULT_SURF)
    ap.add_argument(
        "--sample-types",
        default="01",
        help="comma-separated TCGA sample-type codes to keep (default 01=primary tumour)",
    )
    ap.add_argument("--window", type=int, default=200, help="size of top-N companion file (default 200)")
    ap.add_argument("--outdir", type=Path, default=ROOT / "results" / "w200" / "B1_CESC")
    args = ap.parse_args(argv)

    sample_types = [s.strip() for s in args.sample_types.split(",") if s.strip()]
    args.outdir.mkdir(parents=True, exist_ok=True)

    print(f"[load] expression: {args.expression}")
    expr = load_expression(args.expression)
    print(f"       genes x samples = {expr.shape[0]} x {expr.shape[1]}")

    if args.anchor not in expr.index:
        raise SystemExit(f"anchor gene {args.anchor!r} not found in expression matrix")

    print(f"[load] surfaceome: {args.surfaceome}")
    surf = load_surfaceome(args.surfaceome)
    surf_genes = sorted({g for g in surf["gene"] if g not in ("nan", "None", "")})

    samples = select_samples(expr, sample_types)
    if not samples:
        raise SystemExit(f"no samples matched sample-type codes {sample_types}")
    print(f"[data] sample-type {sample_types}: {len(samples)} samples")

    universe = [g for g in surf_genes if g in expr.index and g != args.anchor]
    n_universe = len(universe)
    print(f"[data] surface-gene universe (present, anchor removed): {n_universe}")

    df = compute_ranking(expr, args.anchor, universe, samples)

    # Write full ranking + top-N window.
    full_csv = args.outdir / f"coexpression_{args.anchor}_surfaceome.csv"
    df.to_csv(full_csv, index=False)
    top_csv = args.outdir / f"top{args.window}.csv"
    df.head(args.window).to_csv(top_csv, index=False)

    report = focus_report(df, args.focus, n_universe)

    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "cohort": "TCGA-CESC",
        "expression_dataset": "UCSC Xena TCGA.CESC.sampleMap/HiSeqV2 (log2 norm_count+1)",
        "surfaceome_source": "Bausch-Fluck et al. 2018 in-silico surfaceome (table S3)",
        "anchor_gene": args.anchor,
        "focus_gene": args.focus,
        "sample_type_codes": sample_types,
        "n_samples": len(samples),
        "surface_gene_universe": n_universe,
        "correlation_primary": "spearman",
        "window": args.window,
        "focus_result": report,
        "top10_spearman": df.head(10)[["gene", "spearman_rank", "spearman_r", "spearman_q"]].to_dict("records"),
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2))

    # Human-readable answer.
    q = report
    if q.get("in_universe"):
        verdict = (
            f"NO -- {args.focus} is NOT the single top surface-gene co-expression "
            f"partner of {args.anchor}."
            if not q["is_top_spearman"]
            else f"YES -- {args.focus} is the top surface-gene co-expression partner of {args.anchor}."
        )
        top_gene = df.iloc[0]["gene"]
        lines = [
            f"# {args.focus} vs {args.anchor} co-expression in TCGA-CESC (honest ranking)",
            "",
            f"**Question:** Is `{args.focus}` the top co-expression partner of "
            f"`{args.anchor}` among surface genes?",
            "",
            f"**Answer:** {verdict}",
            "",
            f"- Cohort: TCGA-CESC, {len(samples)} primary-tumour samples "
            f"(sample-type {sample_types}).",
            f"- Surface-gene universe: {n_universe} surfaceome genes present in the "
            f"matrix (anchor removed).",
            f"- Primary metric: Spearman correlation.",
            "",
            "| metric | value |",
            "| --- | --- |",
            f"| {args.focus} Spearman rank | **#{q['spearman_rank']} of {n_universe}** "
            f"({q['spearman_percentile']}th percentile) |",
            f"| {args.focus} Spearman rho | {q['spearman_r']:.3f} (FDR q={q['spearman_q']:.2e}) |",
            f"| {args.focus} Pearson rank | #{q['pearson_rank']} of {n_universe} |",
            f"| {args.focus} Pearson rho | {q['pearson_r']:.3f} |",
            f"| Actual #1 (Spearman) | {top_gene} (rho={df.iloc[0]['spearman_r']:.3f}) |",
            "",
            "## Top 15 surface-gene partners of "
            f"{args.anchor} (Spearman)",
            "",
            "| rank | gene | spearman_rho | pearson_rho | FDR q |",
            "| --- | --- | --- | --- | --- |",
        ]
        for _, r in df.head(15).iterrows():
            star = "  <-- focus" if r["gene"] == args.focus else ""
            lines.append(
                f"| {int(r['spearman_rank'])} | {r['gene']}{star} | "
                f"{r['spearman_r']:.3f} | {r['pearson_r']:.3f} | {r['spearman_q']:.2e} |"
            )
    else:
        lines = [f"{args.focus} is not present in the surface-gene universe."]

    (args.outdir / "summary.md").write_text("\n".join(lines) + "\n")

    make_plots(expr, df, args.anchor, args.focus, samples, args.window, args.outdir)

    print("\n".join(lines))
    print(f"\n[done] outputs written to {args.outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
