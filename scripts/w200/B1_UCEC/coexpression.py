"""Honest co-expression ranking of TACSTD2 against the cell-surface proteome.

B1 analog for TCGA-UCEC: score every in-silico surfaceome gene against the
anchor (default TACSTD2 / TROP2) and report where the focus gene (default
CLDN4) actually lands. No favourable neighbourhood or single metric is
pre-selected. Spearman is primary; Pearson is reported alongside.

Definitions
-----------
* Surface-gene universe: Bausch-Fluck et al., PNAS 2018, table S3
  ("in silico surfaceome only", 2,886 proteins). A surfaceome row is kept
  if its UniProt gene symbol is in the expression matrix, or — if that
  2018 symbol is absent — if its Ensembl gene id maps to a current GENCODE
  v36 symbol that is present. The anchor itself is removed. This is the
  BRCA analog's symbol intersection, plus a documented Ensembl fallback so
  renamed genes (e.g. PVRL4 → NECTIN4) are not silently dropped.
* Samples: TCGA-UCEC primary tumours (barcode sample-type code ``01``).
* Expression: Xena GDC STAR FPKM-UQ, log2(fpkm-uq + 1), Ensembl rows
  collapsed to symbols via the GENCODE v36 probeMap.

"w200" convention
-----------------
``w200`` is the reporting window (size of ``top200.csv``). It does not
change the full-universe ranking.
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = ROOT / "data"
DEFAULT_EXPR = DATA_DIR / "TCGA-UCEC.star_fpkm-uq.tsv.gz"
DEFAULT_MAP = DATA_DIR / "gencode.v36.annotation.gtf.gene.probemap"
DEFAULT_SURF = DATA_DIR / "table_S3_surfaceome.xlsx"


def load_expression(path: Path, probemap: Path) -> pd.DataFrame:
    """Load Xena STAR matrix and collapse Ensembl rows to HGNC symbols."""
    expr = pd.read_csv(path, sep="\t", index_col=0)
    expr.index = expr.index.astype(str)
    pm = pd.read_csv(probemap, sep="\t")
    id_to_symbol = dict(zip(pm["id"].astype(str), pm["gene"].astype(str)))
    symbols = expr.index.map(lambda i: id_to_symbol.get(i, ""))
    expr = expr.loc[symbols.astype(bool)].copy()
    expr.index = symbols[symbols.astype(bool)]
    expr.index.name = "gene"
    if expr.index.duplicated().any():
        expr = expr.groupby(level=0).mean()
    return expr


def load_surfaceome(path: Path) -> pd.DataFrame:
    """Return the in-silico surfaceome table (header is on the 2nd row)."""
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Unknown extension")
        surf = pd.read_excel(
            path, sheet_name="in silico surfaceome only", engine="openpyxl", header=1
        )
    surf = surf.rename(
        columns={"UniProt gene": "gene", "Ensembl gene": "ensembl_gene"}
    )
    surf = surf[surf["gene"].notna()].copy()
    surf["gene"] = surf["gene"].astype(str)
    return surf


def surface_universe(
    surf: pd.DataFrame, expr_index: pd.Index, ensembl_to_symbol: dict[str, str]
) -> tuple[list[str], dict]:
    """Unique surface symbols present in the matrix (anchor removed later).

    Preference: UniProt gene symbol if measured; else GENCODE symbol of the
    table-S3 Ensembl id. Returns (sorted unique symbols, mapping stats).
    """
    present: set[str] = set()
    n_symbol = 0
    n_ensembl_fallback = 0
    fallback_examples: list[dict] = []
    for gene, ens in zip(surf["gene"], surf["ensembl_gene"]):
        if gene in expr_index:
            present.add(gene)
            n_symbol += 1
            continue
        ens_id = "" if pd.isna(ens) else str(ens).split(".")[0]
        mapped = ensembl_to_symbol.get(ens_id, "")
        if mapped and mapped in expr_index:
            present.add(mapped)
            n_ensembl_fallback += 1
            if len(fallback_examples) < 15:
                fallback_examples.append({"surfaceome_symbol": gene, "mapped_symbol": mapped})
    stats_out = {
        "surfaceome_rows": int(len(surf)),
        "unique_uniprot_symbols": int(surf["gene"].nunique()),
        "matched_by_uniprot_symbol": n_symbol,
        "matched_by_ensembl_fallback": n_ensembl_fallback,
        "unique_symbols_in_matrix": len(present),
        "ensembl_fallback_examples": fallback_examples,
    }
    return sorted(present), stats_out


def select_samples(expr: pd.DataFrame, sample_types: list[str]) -> list[str]:
    """Keep columns whose TCGA sample-type code is in ``sample_types``."""
    keep = []
    for col in expr.columns:
        parts = str(col).split("-")
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
            "mean_log2_expr": mat.mean(axis=1),
            "pct_expressed": (mat > 0).mean(axis=1) * 100.0,
        }
    )
    df["spearman_q"] = _bh_fdr(df["spearman_p"].to_numpy())
    df["pearson_q"] = _bh_fdr(df["pearson_p"].to_numpy())

    # Zero-variance genes yield undefined correlations; they are not ranked.
    n_undefined = int(df["spearman_r"].isna().sum())
    df = df.dropna(subset=["spearman_r"]).copy()
    df.attrs["n_undefined"] = n_undefined

    df = df.sort_values("spearman_r", ascending=False, kind="mergesort").reset_index(drop=True)
    df.insert(1, "spearman_rank", np.arange(1, len(df) + 1))
    df["pearson_rank"] = (
        df["pearson_r"].rank(ascending=False, method="min", na_option="keep")
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
        "spearman_p": float(r["spearman_p"]),
        "spearman_q": float(r["spearman_q"]),
        "spearman_percentile": round(100 * (1 - (s_rank - 1) / n_universe), 3),
        "pearson_rank": p_rank,
        "pearson_r": float(r["pearson_r"]),
        "pearson_p": float(r["pearson_p"]),
        "pearson_q": float(r["pearson_q"]),
        "pearson_percentile": round(100 * (1 - (p_rank - 1) / n_universe), 3),
        "mean_log2_expr": float(r["mean_log2_expr"]),
        "pct_expressed": float(r["pct_expressed"]),
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

    if focus in expr.index:
        x = expr.loc[anchor, samples].to_numpy(dtype=float)
        y = expr.loc[focus, samples].to_numpy(dtype=float)
        row = df[df["gene"] == focus].iloc[0]
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.scatter(x, y, s=8, alpha=0.4, edgecolors="none")
        ax.set_xlabel(f"{anchor} log2(fpkm-uq+1)")
        ax.set_ylabel(f"{focus} log2(fpkm-uq+1)")
        ax.set_title(
            f"{anchor} vs {focus} (TCGA-UCEC, n={len(samples)})\n"
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
    ax.set_title(f"Top surface-gene co-expression with {anchor}\nTCGA-UCEC (focus: {focus})")
    fig.tight_layout()
    fig.savefig(outdir / f"top_partners_{anchor}.png", dpi=150)
    plt.close(fig)


def _table_rows(df: pd.DataFrame, focus: str, n: int) -> list[str]:
    show = df.head(n)
    if focus in df["gene"].values and focus not in set(show["gene"]):
        show = pd.concat([show, df[df["gene"] == focus]], ignore_index=True)
    lines = []
    for _, r in show.iterrows():
        star = "  <-- focus" if r["gene"] == focus else ""
        lines.append(
            f"| {int(r['spearman_rank'])} | {r['gene']}{star} | "
            f"{r['spearman_r']:.3f} | {r['pearson_r']:.3f} | {r['spearman_q']:.2e} |"
        )
    return lines


def write_summary_md(
    df: pd.DataFrame,
    report: dict,
    n_samples: int,
    n_universe: int,
    sample_types: list[str],
    n_patients: int,
    mapping_stats: dict,
    outdir: Path,
    anchor: str,
    focus: str,
    n_undefined: int = 0,
    sensitivities: dict | None = None,
) -> list[str]:
    q = report
    if q.get("in_universe"):
        verdict = (
            f"NO -- {focus} is NOT the single top surface-gene co-expression "
            f"partner of {anchor}."
            if not q["is_top_spearman"]
            else f"YES -- {focus} is the top surface-gene co-expression partner of {anchor}."
        )
        top_gene = df.iloc[0]["gene"]
        examples = ", ".join(
            f"{e['surfaceome_symbol']}→{e['mapped_symbol']}"
            for e in mapping_stats.get("ensembl_fallback_examples", [])[:8]
        )
        lines = [
            f"# {focus} vs {anchor} co-expression in TCGA-UCEC (honest ranking)",
            "",
            f"**Question:** Is `{focus}` the top co-expression partner of "
            f"`{anchor}` among surface genes?",
            "",
            f"**Answer:** {verdict}",
            "",
            f"- Cohort: TCGA-UCEC, {n_samples} primary-tumour samples "
            f"(sample-type {sample_types}; {n_patients} unique patients).",
            f"- Surface-gene universe: {n_universe} surfaceome genes with defined "
            f"Spearman (anchor removed; {n_undefined} zero-variance genes dropped).",
            f"- Primary metric: Spearman correlation.",
            f"- Expression: Xena GDC STAR FPKM-UQ, log2(fpkm-uq+1).",
            "",
            "| metric | value |",
            "| --- | --- |",
            f"| {focus} Spearman rank | **#{q['spearman_rank']} of {n_universe}** "
            f"({q['spearman_percentile']}th percentile) |",
            f"| {focus} Spearman rho | {q['spearman_r']:.3f} "
            f"(p={q['spearman_p']:.2e}, FDR q={q['spearman_q']:.2e}) |",
            f"| {focus} Pearson rank | #{int(q['pearson_rank'])} of {n_universe} |",
            f"| {focus} Pearson rho | {q['pearson_r']:.3f} |",
            f"| Actual #1 (Spearman) | {top_gene} (rho={df.iloc[0]['spearman_r']:.3f}) |",
            "",
            "## Top 25 surface-gene partners of "
            f"{anchor} (Spearman)",
            "",
            "| rank | gene | spearman_rho | pearson_rho | FDR q |",
            "| --- | --- | --- | --- | --- |",
        ]
        lines += _table_rows(df, focus, 25)
        if sensitivities:
            lines += ["", "## Sensitivities (same Spearman, not cherry-picked)", ""]
            for name, rec in sensitivities.items():
                lines.append(
                    f"- **{name}:** {focus} rank #{rec['spearman_rank']} of "
                    f"{rec['n_universe']} (rho={rec['spearman_r']:.3f}); "
                    f"#1 is {rec['top_gene']} (rho={rec['top_r']:.3f})."
                )
        lines += [
            "",
            "## Surfaceome mapping",
            "",
            f"- Table S3 rows: {mapping_stats['surfaceome_rows']}; unique UniProt symbols: "
            f"{mapping_stats['unique_uniprot_symbols']}.",
            f"- Matched by 2018 UniProt symbol: {mapping_stats['matched_by_uniprot_symbol']}.",
            f"- Recovered by Ensembl id → GENCODE v36 symbol: "
            f"{mapping_stats['matched_by_ensembl_fallback']} (e.g. {examples}).",
            "",
            "## Caveats",
            "",
            "- Bulk-tumour mRNA co-expression (mixes tumour/stroma/immune), not protein or single-cell.",
            "- UCEC mixes endometrioid and serous histologies; this ranking is not histology-stratified.",
            "- Surfaceome symbols are 2018 UniProt names; renamed genes (PVRL4→NECTIN4) enter via Ensembl ids.",
            "- Rank is among surfaceome genes only, not the full transcriptome.",
            "- Pearson and Spearman disagree on exact rank (CLDN4 is higher by Pearson); Spearman is primary.",
        ]
    else:
        lines = [f"{focus} is not present in the surface-gene universe."]
    (outdir / "summary.md").write_text("\n".join(lines) + "\n")
    return lines


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--anchor", default="TACSTD2")
    ap.add_argument("--focus", default="CLDN4")
    ap.add_argument("--expression", type=Path, default=DEFAULT_EXPR)
    ap.add_argument("--probemap", type=Path, default=DEFAULT_MAP)
    ap.add_argument("--surfaceome", type=Path, default=DEFAULT_SURF)
    ap.add_argument(
        "--sample-types",
        default="01",
        help="comma-separated TCGA sample-type codes (default 01=primary tumour)",
    )
    ap.add_argument("--window", type=int, default=200)
    ap.add_argument("--outdir", type=Path, default=ROOT / "results" / "w200" / "B1_UCEC")
    args = ap.parse_args(argv)

    sample_types = [s.strip() for s in args.sample_types.split(",") if s.strip()]
    args.outdir.mkdir(parents=True, exist_ok=True)

    print(f"[load] expression: {args.expression}")
    expr = load_expression(args.expression, args.probemap)
    print(f"       genes x samples = {expr.shape[0]} x {expr.shape[1]}")

    if args.anchor not in expr.index:
        raise SystemExit(f"anchor gene {args.anchor!r} not found in expression matrix")

    print(f"[load] surfaceome: {args.surfaceome}")
    surf = load_surfaceome(args.surfaceome)

    pm = pd.read_csv(args.probemap, sep="\t")
    ensembl_to_symbol = {
        str(i).split(".")[0]: str(g) for i, g in zip(pm["id"], pm["gene"])
    }

    samples = select_samples(expr, sample_types)
    if not samples:
        raise SystemExit(f"no samples matched sample-type codes {sample_types}")
    n_patients = len({"-".join(str(c).split("-")[:3]) for c in samples})
    print(f"[data] sample-type {sample_types}: {len(samples)} samples / {n_patients} patients")

    present, mapping_stats = surface_universe(surf, expr.index, ensembl_to_symbol)
    universe = [g for g in present if g != args.anchor]
    n_universe = len(universe)
    print(f"[data] surface-gene universe (present, anchor removed): {n_universe}")
    print(f"[data] mapping: {mapping_stats}")

    df = compute_ranking(expr, args.anchor, universe, samples)
    n_undefined = int(df.attrs.get("n_undefined", 0))
    n_universe = len(df)
    print(
        f"[data] ranked surface genes: {n_universe} "
        f"(dropped {n_undefined} with undefined Spearman / zero variance)"
    )

    # Independent spot-check of the focus pair (do not invent; recompute).
    if args.focus in expr.index:
        a = expr.loc[args.anchor, samples].to_numpy(dtype=float)
        b = expr.loc[args.focus, samples].to_numpy(dtype=float)
        spr = stats.spearmanr(a, b)
        per = stats.pearsonr(a, b)
        print(
            f"[check] scipy.stats {args.anchor} vs {args.focus}: "
            f"spearman={spr.statistic:.6f} p={spr.pvalue:.3e}; "
            f"pearson={per.statistic:.6f} p={per.pvalue:.3e}"
        )

    full_csv = args.outdir / f"coexpression_{args.anchor}_surfaceome.csv"
    df.to_csv(full_csv, index=False)
    top_csv = args.outdir / f"top{args.window}.csv"
    df.head(args.window).to_csv(top_csv, index=False)

    report = focus_report(df, args.focus, n_universe)

    # Sensitivities: BRCA-style symbol-only universe, and one sample per patient.
    symbol_only = sorted(
        {g for g in surf["gene"] if g in expr.index and g != args.anchor}
    )
    df_sym = compute_ranking(expr, args.anchor, symbol_only, samples)
    df_sym = df_sym.dropna(subset=["spearman_r"]) if "spearman_r" in df_sym else df_sym
    n_sym = len(df_sym)
    focus_sym = focus_report(df_sym, args.focus, n_sym)

    one_per: dict[str, str] = {}
    for col in samples:
        patient = "-".join(str(col).split("-")[:3])
        one_per.setdefault(patient, col)
    samples_1pp = list(one_per.values())
    df_1pp = compute_ranking(expr, args.anchor, list(df["gene"]), samples_1pp)
    n_1pp = len(df_1pp)
    focus_1pp = focus_report(df_1pp, args.focus, n_1pp)

    sensitivities = {
        "symbol_only_no_ensembl_fallback": {
            "n_universe": n_sym,
            "spearman_rank": focus_sym.get("spearman_rank"),
            "spearman_r": focus_sym.get("spearman_r"),
            "top_gene": str(df_sym.iloc[0]["gene"]) if n_sym else None,
            "top_r": float(df_sym.iloc[0]["spearman_r"]) if n_sym else None,
        },
        "one_sample_per_patient": {
            "n_samples": len(samples_1pp),
            "n_universe": n_1pp,
            "spearman_rank": focus_1pp.get("spearman_rank"),
            "spearman_r": focus_1pp.get("spearman_r"),
            "top_gene": str(df_1pp.iloc[0]["gene"]) if n_1pp else None,
            "top_r": float(df_1pp.iloc[0]["spearman_r"]) if n_1pp else None,
        },
    }

    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "cohort": "TCGA-UCEC",
        "expression_dataset": (
            "UCSC Xena GDC hub TCGA-UCEC.star_fpkm-uq (log2(fpkm-uq+1), GENCODE v36)"
        ),
        "surfaceome_source": "Bausch-Fluck et al. 2018 in-silico surfaceome (table S3)",
        "anchor_gene": args.anchor,
        "focus_gene": args.focus,
        "sample_type_codes": sample_types,
        "n_samples": len(samples),
        "n_patients": n_patients,
        "surface_gene_universe": n_universe,
        "n_undefined_correlation": n_undefined,
        "surfaceome_mapping": mapping_stats,
        "correlation_primary": "spearman",
        "window": args.window,
        "focus_result": report,
        "sensitivities": sensitivities,
        "top10_spearman": df.head(10)[
            ["gene", "spearman_rank", "spearman_r", "spearman_q"]
        ].to_dict("records"),
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    lines = write_summary_md(
        df,
        report,
        len(samples),
        n_universe,
        sample_types,
        n_patients,
        mapping_stats,
        args.outdir,
        args.anchor,
        args.focus,
        n_undefined=n_undefined,
        sensitivities=sensitivities,
    )
    make_plots(expr, df, args.anchor, args.focus, samples, args.window, args.outdir)

    print("\n".join(lines))
    print(f"\n[done] outputs written to {args.outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
