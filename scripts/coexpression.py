"""Honest co-expression ranking of a hub gene against the cell-surface proteome.

This is the "B1 analog" for TCGA-PRAD: it reproduces, for the prostate
adenocarcinoma cohort, the panel-B1 style analysis of asking how a candidate
surface marker (default CLDN4) ranks among *all* surface genes when ordered
by co-expression with an anchor gene (default TACSTD2 / TROP2).

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
* Samples: TCGA-PRAD primary tumours (barcode sample-type code ``01``) by
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
DEFAULT_EXPR = DATA_DIR / "TCGA-PRAD.HiSeqV2.gz"
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

    # Zero-variance genes yield undefined correlations; drop them rather than
    # invent a rank. This is reported in the summary as n_undefined.
    finite = np.isfinite(df["spearman_r"].to_numpy()) & np.isfinite(df["pearson_r"].to_numpy())
    df = df.loc[finite].copy()
    df = df.sort_values("spearman_r", ascending=False, kind="mergesort").reset_index(drop=True)
    df.insert(1, "spearman_rank", np.arange(1, len(df) + 1))
    df["pearson_rank"] = df["pearson_r"].rank(ascending=False, method="min").astype(int)
    return df


def bootstrap_focus(
    expr: pd.DataFrame,
    anchor: str,
    universe: list[str],
    samples: list[str],
    focus: str,
    n_boot: int = 1000,
    seed: int = 20260816,
) -> dict:
    """Case-resample Spearman rho and how often ``focus`` is rank #1.

    Ranking the full surfaceome each resample is the honest stability check:
    a 0.001 lead can evaporate. Returns percentile CI and win-rate.
    """
    rng = np.random.default_rng(seed)
    genes = [focus] + [g for g in universe if g != focus]
    mat = expr.loc[genes, samples].to_numpy(dtype=float)
    anchor_vec = expr.loc[anchor, samples].to_numpy(dtype=float)
    n = len(samples)
    g = mat.shape[0]
    rhos = np.empty(n_boot, dtype=float)
    wins = 0
    # Track the runner-up when focus loses, for an honest near-tie note.
    usurper = {}
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        a = stats.rankdata(anchor_vec[idx])
        # Rank each gene on the resampled samples (G x N).
        ranked = np.apply_along_axis(stats.rankdata, 1, mat[:, idx])
        r = _corr_matrix_vs_vector(ranked, a)
        rhos[i] = r[0]
        k = int(np.nanargmax(r))
        if k == 0:
            wins += 1
        else:
            name = genes[k]
            usurper[name] = usurper.get(name, 0) + 1
    lo, med, hi = np.percentile(rhos, [2.5, 50, 97.5])
    top_usurpers = sorted(usurper.items(), key=lambda kv: -kv[1])[:5]
    return {
        "n_boot": n_boot,
        "seed": seed,
        "spearman_r_ci95": [float(lo), float(hi)],
        "spearman_r_median": float(med),
        "focus_rank1_rate": wins / n_boot,
        "top_usurpers": [{"gene": g, "n_wins": c} for g, c in top_usurpers],
    }


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
            f"{anchor} vs {focus} (TCGA-PRAD, n={len(samples)})\n"
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
    ax.set_title(f"Top surface-gene co-expression with {anchor}\nTCGA-PRAD (focus: {focus})")
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
    ap.add_argument("--n-boot", type=int, default=1000, help="bootstrap resamples for rank stability (0 to skip)")
    ap.add_argument("--outdir", type=Path, default=ROOT / "results" / "w200" / "B1_PRAD")
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
    n_ranked = len(df)
    n_undefined = n_universe - n_ranked
    if n_undefined:
        print(f"[data] dropped {n_undefined} zero-variance / undefined-correlation genes; ranked {n_ranked}")

    # Write full ranking + top-N window.
    full_csv = args.outdir / f"coexpression_{args.anchor}_surfaceome.csv"
    df.to_csv(full_csv, index=False)
    top_csv = args.outdir / f"top{args.window}.csv"
    df.head(args.window).to_csv(top_csv, index=False)

    report = focus_report(df, args.focus, n_ranked)
    boot = None
    if args.n_boot > 0 and report.get("in_universe"):
        print(f"[boot] {args.n_boot} case-resamples of the full surfaceome ranking")
        boot = bootstrap_focus(
            expr, args.anchor, list(df["gene"]), samples, args.focus, n_boot=args.n_boot
        )
        print(
            f"       {args.focus} Spearman rho 95% CI "
            f"[{boot['spearman_r_ci95'][0]:.3f}, {boot['spearman_r_ci95'][1]:.3f}]; "
            f"rank-1 rate {boot['focus_rank1_rate']:.3f}"
        )

    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "cohort": "TCGA-PRAD",
        "expression_dataset": "UCSC Xena TCGA.PRAD.sampleMap/HiSeqV2 (log2 norm_count+1)",
        "surfaceome_source": "Bausch-Fluck et al. 2018 in-silico surfaceome (table S3)",
        "anchor_gene": args.anchor,
        "focus_gene": args.focus,
        "sample_type_codes": sample_types,
        "n_samples": len(samples),
        "surface_gene_universe_present": n_universe,
        "surface_gene_universe": n_ranked,
        "n_undefined_correlation": n_undefined,
        "correlation_primary": "spearman",
        "window": args.window,
        "focus_result": report,
        "bootstrap": boot,
        "top10_spearman": df.head(10)[["gene", "spearman_rank", "spearman_r", "spearman_q"]].to_dict("records"),
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2))

    # Human-readable answer.
    q = report
    if q.get("in_universe"):
        if not q["is_top_spearman"]:
            verdict = (
                f"NO -- {args.focus} is NOT the single top surface-gene co-expression "
                f"partner of {args.anchor}."
            )
        elif boot is not None and boot["focus_rank1_rate"] < 0.5:
            verdict = (
                f"POINT ESTIMATE YES, STABILITY NO -- {args.focus} is Spearman #1 "
                f"in this sample, but the lead is not stable "
                f"(rank-1 in only {100 * boot['focus_rank1_rate']:.0f}% of bootstraps; "
                f"Pearson rank #{q['pearson_rank']})."
            )
        else:
            verdict = (
                f"YES -- {args.focus} is the top surface-gene co-expression partner of {args.anchor}."
            )
        top_gene = df.iloc[0]["gene"]
        lines = [
            f"# {args.focus} vs {args.anchor} co-expression in TCGA-PRAD (honest ranking)",
            "",
            f"**Question:** Is `{args.focus}` the top co-expression partner of "
            f"`{args.anchor}` among surface genes?",
            "",
            f"**Answer:** {verdict}",
            "",
            f"- Cohort: TCGA-PRAD, {len(samples)} primary-tumour samples "
            f"(sample-type {sample_types}).",
            f"- Surface-gene universe: {n_ranked} surfaceome genes with defined "
            f"correlation (anchor removed"
            + (f"; {n_undefined} zero-variance genes dropped" if n_undefined else "")
            + ").",
            f"- Primary metric: Spearman correlation.",
            "",
            "| metric | value |",
            "| --- | --- |",
            f"| {args.focus} Spearman rank | **#{q['spearman_rank']} of {n_ranked}** "
            f"({q['spearman_percentile']}th percentile) |",
            f"| {args.focus} Spearman rho | {q['spearman_r']:.3f} (FDR q={q['spearman_q']:.2e}) |",
            f"| {args.focus} Pearson rank | #{q['pearson_rank']} of {n_ranked} |",
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
        runner = df.iloc[1] if len(df) > 1 else None
        lines += [
            "",
            "## Honest caveats",
            "",
            "- Bulk-tumour mRNA co-expression (mixes tumour/stroma/immune), not protein "
            "or single-cell. Xena HiSeqV2 is log2(norm_count+1), not TPM.",
            "- Primary metric is Spearman, matching the BRCA B1 analog. Pearson can "
            "disagree: here Pearson rank is "
            f"#{q['pearson_rank']} (#{1} is {df.sort_values('pearson_r', ascending=False).iloc[0]['gene']}).",
        ]
        if runner is not None:
            delta = float(df.iloc[0]["spearman_r"] - runner["spearman_r"])
            lines.append(
                f"- The Spearman lead over #{2} {runner['gene']} is only "
                f"{delta:.4f}. Treat '#1' as a near-tie, not a unique winner."
            )
        if boot is not None:
            usurper_txt = ", ".join(
                f"{u['gene']} ({u['n_wins']}/{boot['n_boot']})" for u in boot["top_usurpers"][:3]
            ) or "none"
            lines += [
                f"- Bootstrap (n={boot['n_boot']}, seed {boot['seed']}): "
                f"{args.focus} Spearman rho 95% CI "
                f"[{boot['spearman_r_ci95'][0]:.3f}, {boot['spearman_r_ci95'][1]:.3f}]; "
                f"it is rank #1 in **{100 * boot['focus_rank1_rate']:.1f}%** of resamples. "
                f"Most frequent usurpers: {usurper_txt}.",
            ]
        lines += [
            "- Surface universe is the 2018 in-silico surfaceome (legacy HGNC symbols, "
            "e.g. PVRL4 = NECTIN4, PVRL2 = NECTIN2). Four zero-variance surface genes "
            "were dropped rather than ranked.",
            "- One primary-tumour aliquot per patient in this matrix (n=497 code-01; "
            "52 adjacent-normal and 1 metastatic aliquot exist and were not used).",
        ]
    else:
        lines = [f"{args.focus} is not present in the surface-gene universe."]

    (args.outdir / "summary.md").write_text("\n".join(lines) + "\n")

    make_plots(expr, df, args.anchor, args.focus, samples, args.window, args.outdir)

    print("\n".join(lines))
    print(f"\n[done] outputs written to {args.outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
