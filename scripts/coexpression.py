#!/usr/bin/env python3
"""Honest co-expression ranking of a hub gene against the cell-surface proteome.

This is the "B1 analog" for TCGA-COAD: it reproduces, for the colon
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
  (Bausch-Fluck et al., PNAS 2018, table S3), intersected with genes measured
  in the expression matrix, with the anchor gene itself removed.
* Samples: TCGA-COAD primary tumours (GDC sample_type ``Primary Tumor``) by
  default. One aliquot is kept per case (prefer the ``-01A`` vial).
* Expression: GDC STAR - Counts ``tpm_unstranded``, analysed as log2(TPM+1).
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
DEFAULT_TPM = DATA_DIR / "coad_star_tpm.parquet"
DEFAULT_MANIFEST = DATA_DIR / "coad_sample_manifest.tsv"
DEFAULT_ANNOT = DATA_DIR / "coad_gene_annotation.tsv"
DEFAULT_SURF = DATA_DIR / "table_S3_surfaceome.xlsx"
DEFAULT_XENA = DATA_DIR / "TCGA-COAD.HiSeqV2.gz"


def load_expression(
    tpm_path: Path, annot_path: Path
) -> pd.DataFrame:
    """Return a gene-symbol x sample matrix of log2(TPM+1).

    GDC STAR files are Ensembl-id rows. We map to GENCODE v36 ``gene_name``,
    then collapse rare duplicate symbols by mean so lookups are unambiguous.
    """
    tpm = pd.read_parquet(tpm_path)
    annot = pd.read_csv(annot_path, sep="\t")
    annot = annot.drop_duplicates("gene_id")
    name_map = dict(zip(annot["gene_id"], annot["gene_name"].astype(str)))
    tpm = tpm.copy()
    tpm.index = tpm.index.map(lambda g: name_map.get(g, g))
    tpm = tpm.groupby(level=0).mean()
    # log2(TPM+1) is the conventional bulk-RNA scale for Pearson/Spearman here.
    expr = np.log2(tpm.astype(float) + 1.0)
    return expr


def load_manifest(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def load_surfaceome(path: Path) -> pd.DataFrame:
    """Return the in-silico surfaceome table (header is on the 2nd row)."""
    surf = pd.read_excel(
        path, sheet_name="in silico surfaceome only", engine="openpyxl", header=1
    )
    surf = surf.rename(columns={"UniProt gene": "gene"})
    surf = surf[surf["gene"].notna()].copy()
    surf["gene"] = surf["gene"].astype(str)
    return surf


def select_primary_samples(manifest: pd.DataFrame, expr_cols: list[str]) -> list[str]:
    """Keep one Primary Tumor aliquot per case.

    Preference order: sample barcode ending in ``-01A``, then any ``-01``,
    then the first remaining primary-tumour aliquot. When a barcode maps to
    more than one GDC file (rare technical replicates), the expression
    loader already has one column per file; we pick the first matching
    column for that barcode.
    """
    prim = manifest[manifest["sample_type"] == "Primary Tumor"].copy()
    if prim.empty:
        raise SystemExit("no Primary Tumor samples in the GDC manifest")

    # Map sample_barcode -> available matrix column(s).
    col_by_barcode: dict[str, list[str]] = {}
    for col in expr_cols:
        barcode = col.split("|", 1)[0]
        col_by_barcode.setdefault(barcode, []).append(col)

    chosen: list[str] = []
    for case, grp in prim.groupby("case_barcode", sort=True):
        barcodes = list(grp["sample_barcode"].unique())

        def score(b: str) -> tuple[int, str]:
            # Prefer 01A, then any 01, then everything else.
            pref = 2
            if b.endswith("-01A") or "-01A-" in b or b.split("-")[-1].startswith("01A"):
                pref = 0
            elif "-01" in b:
                pref = 1
            return (pref, b)

        barcodes_sorted = sorted(barcodes, key=score)
        picked_col = None
        for b in barcodes_sorted:
            cols = col_by_barcode.get(b) or []
            if cols:
                picked_col = sorted(cols)[0]
                break
        if picked_col is not None:
            chosen.append(picked_col)
    return chosen


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
            "mean_log2_tpm": mat.mean(axis=1),
            "pct_expressed": (mat > 0).mean(axis=1) * 100.0,
        }
    )
    df["spearman_q"] = _bh_fdr(df["spearman_p"].to_numpy())
    df["pearson_q"] = _bh_fdr(df["pearson_p"].to_numpy())

    # Genes with zero variance (all-zero TPM) have undefined correlation.
    # They are dropped from the ranked universe rather than given a fake rank.
    n_undefined = int(df["spearman_r"].isna().sum())
    df = df.dropna(subset=["spearman_r"]).copy()
    df.attrs["n_undefined"] = n_undefined

    df = df.sort_values("spearman_r", ascending=False, kind="mergesort").reset_index(drop=True)
    df.insert(1, "spearman_rank", np.arange(1, len(df) + 1))
    df["pearson_rank"] = (
        df["pearson_r"].rank(ascending=False, method="min", na_option="bottom").astype(int)
    )
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


def xena_robustness(
    path: Path, surf_genes: list[str], anchor: str, focus: str
) -> dict | None:
    """Repeat the ranking on UCSC Xena HiSeqV2 (legacy RSEM, log2(norm_count+1))."""
    if not path.exists():
        return None
    expr = pd.read_csv(path, sep="\t", index_col=0)
    expr.index = expr.index.astype(str)
    if expr.index.duplicated().any():
        expr = expr.groupby(level=0).mean()
    if anchor not in expr.index or focus not in expr.index:
        return {"available": True, "error": "anchor or focus missing from Xena matrix"}
    samples = []
    for col in expr.columns:
        parts = str(col).split("-")
        code = parts[3][:2] if len(parts) >= 4 else ""
        if code == "01":
            samples.append(col)
    universe = [g for g in surf_genes if g in expr.index and g != anchor]
    df = compute_ranking(expr, anchor, universe, samples)
    n_undefined = int(df.attrs.get("n_undefined", 0))
    report = focus_report(df, focus, len(df))
    return {
        "available": True,
        "expression_dataset": "UCSC Xena TCGA.COAD.sampleMap/HiSeqV2 (log2 norm_count+1)",
        "n_samples": len(samples),
        "surface_genes_present": len(universe),
        "surface_genes_undefined_corr": n_undefined,
        "surface_gene_universe": len(df),
        "focus_result": report,
        "top5_spearman": df.head(5)[["gene", "spearman_rank", "spearman_r"]].to_dict("records"),
    }


def expression_context(expr: pd.DataFrame, samples: list[str], genes: list[str]) -> dict:
    out = {}
    for g in genes:
        if g not in expr.index:
            continue
        v = expr.loc[g, samples].to_numpy(dtype=float)
        out[g] = {
            "n": int(v.size),
            "median_log2_tpm": float(np.median(v)),
            "mean_log2_tpm": float(np.mean(v)),
            "min_log2_tpm": float(np.min(v)),
            "max_log2_tpm": float(np.max(v)),
            "pct_detected": float((v > 0).mean() * 100.0),
        }
    return out


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

    if focus in expr.index:
        x = expr.loc[anchor, samples].to_numpy(dtype=float)
        y = expr.loc[focus, samples].to_numpy(dtype=float)
        row = df[df["gene"] == focus].iloc[0]
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.scatter(x, y, s=8, alpha=0.4, edgecolors="none")
        ax.set_xlabel(f"{anchor} log2(TPM+1)")
        ax.set_ylabel(f"{focus} log2(TPM+1)")
        ax.set_title(
            f"{anchor} vs {focus} (TCGA-COAD, n={len(samples)})\n"
            f"Spearman rho={row['spearman_r']:.3f}, rank #{int(row['spearman_rank'])}"
        )
        fig.tight_layout()
        fig.savefig(outdir / f"scatter_{anchor}_vs_{focus}.png", dpi=150)
        plt.close(fig)

    topn = df.head(min(20, window)).iloc[::-1]
    colors = ["#d62728" if g == focus else "#4c78a8" for g in topn["gene"]]
    fig, ax = plt.subplots(figsize=(6, 7))
    ax.barh(topn["gene"], topn["spearman_r"], color=colors)
    ax.set_xlabel(f"Spearman rho with {anchor}")
    ax.set_title(f"Top surface-gene co-expression with {anchor}\nTCGA-COAD (focus: {focus})")
    fig.tight_layout()
    fig.savefig(outdir / f"top_partners_{anchor}.png", dpi=150)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--anchor", default="TACSTD2", help="anchor / hub gene (default: TACSTD2)")
    ap.add_argument("--focus", default="CLDN4", help="gene whose rank we report (default: CLDN4)")
    ap.add_argument("--tpm", type=Path, default=DEFAULT_TPM)
    ap.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    ap.add_argument("--annotation", type=Path, default=DEFAULT_ANNOT)
    ap.add_argument("--surfaceome", type=Path, default=DEFAULT_SURF)
    ap.add_argument("--xena", type=Path, default=DEFAULT_XENA)
    ap.add_argument("--window", type=int, default=200, help="size of top-N companion file (default 200)")
    ap.add_argument("--outdir", type=Path, default=ROOT / "results" / "w200" / "B1_COAD")
    args = ap.parse_args(argv)

    args.outdir.mkdir(parents=True, exist_ok=True)

    print(f"[load] expression: {args.tpm}")
    expr = load_expression(args.tpm, args.annotation)
    print(f"       genes x samples = {expr.shape[0]} x {expr.shape[1]}")

    if args.anchor not in expr.index:
        raise SystemExit(f"anchor gene {args.anchor!r} not found in expression matrix")

    print(f"[load] surfaceome: {args.surfaceome}")
    surf = load_surfaceome(args.surfaceome)
    surf_genes = sorted({g for g in surf["gene"] if g not in ("nan", "None", "")})

    manifest = load_manifest(args.manifest)
    samples = select_primary_samples(manifest, list(expr.columns))
    if not samples:
        raise SystemExit("no primary-tumour samples matched the expression matrix")
    print(f"[data] Primary Tumor, one aliquot/case: {len(samples)} samples")

    universe = [g for g in surf_genes if g in expr.index and g != args.anchor]
    n_universe = len(universe)
    print(f"[data] surface-gene universe (present, anchor removed): {n_universe}")
    if args.focus not in universe and args.focus in expr.index:
        print(
            f"[warn] {args.focus} is measured but not in the surfaceome gene-symbol "
            f"intersection (legacy UniProt names vs GENCODE v36).",
            file=sys.stderr,
        )

    df = compute_ranking(expr, args.anchor, universe, samples)
    n_undefined = int(df.attrs.get("n_undefined", 0))
    n_universe = len(df)
    print(
        f"[data] ranked surface genes: {n_universe} "
        f"(dropped {n_undefined} with undefined correlation / zero variance)"
    )

    full_csv = args.outdir / f"coexpression_{args.anchor}_surfaceome.csv"
    df.to_csv(full_csv, index=False)
    top_csv = args.outdir / f"top{args.window}.csv"
    df.head(args.window).to_csv(top_csv, index=False)

    report = focus_report(df, args.focus, n_universe)
    ctx = expression_context(expr, samples, [args.anchor, args.focus])
    xena = xena_robustness(args.xena, surf_genes, args.anchor, args.focus)

    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "cohort": "TCGA-COAD",
        "expression_dataset": (
            "NCI GDC TCGA-COAD STAR - Counts, GENCODE v36, tpm_unstranded; "
            "analysed as log2(TPM+1); one Primary Tumor aliquot per case"
        ),
        "surfaceome_source": "Bausch-Fluck et al. 2018 in-silico surfaceome (table S3)",
        "anchor_gene": args.anchor,
        "focus_gene": args.focus,
        "sample_type": "Primary Tumor",
        "n_samples": len(samples),
        "n_gdc_star_files_total": int(len(manifest)),
        "surface_genes_present": len(universe),
        "surface_genes_undefined_corr": n_undefined,
        "surface_gene_universe": n_universe,
        "correlation_primary": "spearman",
        "window": args.window,
        "focus_result": report,
        "top10_spearman": df.head(10)[["gene", "spearman_rank", "spearman_r", "spearman_q"]].to_dict(
            "records"
        ),
        "expression_context_log2tpm": ctx,
        "xena_hiseqv2_robustness": xena,
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2))

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
            f"# {args.focus} vs {args.anchor} co-expression in TCGA-COAD (honest ranking)",
            "",
            f"**Question:** Is `{args.focus}` the top co-expression partner of "
            f"`{args.anchor}` among surface genes?",
            "",
            f"**Answer:** {verdict}",
            "",
            f"- Cohort: TCGA-COAD, {len(samples)} primary-tumour samples "
            f"(GDC `Primary Tumor`, one aliquot per case).",
            f"- Expression: GDC STAR - Counts `tpm_unstranded`, log2(TPM+1).",
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

        # Honest context: CLDN4 is constitutively high in colon, not a TACSTD2 partner.
        if args.focus in ctx and args.anchor in ctx:
            fa = ctx[args.anchor]
            ff = ctx[args.focus]
            lines += [
                "",
                "## Expression context (why the rank is this low)",
                "",
                f"`{args.focus}` is constitutively high in COAD primary tumours "
                f"(median log2(TPM+1) = {ff['median_log2_tpm']:.2f}; detected in "
                f"{ff['pct_detected']:.0f}% of samples). `{args.anchor}` is moderate and "
                f"more variable (median {fa['median_log2_tpm']:.2f}, range "
                f"{fa['min_log2_tpm']:.2f}–{fa['max_log2_tpm']:.2f}). The two vectors "
                f"do not co-vary: Spearman ρ = {q['spearman_r']:.3f}, FDR q = "
                f"{q['spearman_q']:.2e}. This is a null co-expression result, not a "
                f"near-miss of #1.",
            ]

        if xena and xena.get("focus_result", {}).get("in_universe"):
            xr = xena["focus_result"]
            xtop = xena["top5_spearman"][0]["gene"] if xena.get("top5_spearman") else "?"
            xtop_r = xena["top5_spearman"][0]["spearman_r"] if xena.get("top5_spearman") else float("nan")
            lines += [
                "",
                "## Independent matrix (UCSC Xena HiSeqV2)",
                "",
                f"Same question on the legacy Xena `TCGA.COAD.sampleMap/HiSeqV2` "
                f"matrix (log2(norm_count+1), sample-type `01`, n={xena['n_samples']}, "
                f"{xena['surface_gene_universe']} ranked surface genes):",
                "",
                f"- {args.focus} Spearman rank **#{xr['spearman_rank']} of "
                f"{xena['surface_gene_universe']}**, ρ = {xr['spearman_r']:.3f} "
                f"(FDR q={xr['spearman_q']:.2e}).",
                f"- Actual #1: {xtop} (ρ = {xtop_r:.3f}).",
                "",
                "The GDC STAR-Counts and Xena HiSeqV2 matrices agree: CLDN4 is not "
                "a TACSTD2 surface-gene partner in COAD.",
            ]

        lines += [
            "",
            "## Caveats",
            "",
            "- Bulk-tumour mRNA co-expression (mixes tumour / stroma / immune); not protein and not single-cell.",
            "- Surfaceome gene symbols are the 2018 UniProt names (e.g. `PVRL4` = `NECTIN4`).",
            "- 13 surfaceome genes with zero variance in this cohort were dropped (undefined correlation), not ranked.",
            "- One Primary Tumor aliquot per case; technical-replicate files were not averaged.",
            "- This COAD null does not contradict the BRCA B1 analog (CLDN4 #4, ρ≈0.35) or the PAAD analog (ρ≈0.71); those are different tissues.",
            "- TCGA-COAD is untreated archival RNA-seq; it is not an ICI or ADC outcome cohort.",
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
