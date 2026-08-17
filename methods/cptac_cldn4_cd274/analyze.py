#!/usr/bin/env python3
"""CLDN4 protein vs CD274 and IFN proteins in CPTAC LUAD / LSCC TMT.

Public freeze v1.2 tumor protein only. Pairwise-complete Spearman ρ,
two-sided p, 2,000-resample bootstrap 95% CI (seed 20260817). Honest n.
IFN proteins are a locked list; absent or n<8 rows are reported as
not quantified, not as ρ=0.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from genes import (  # noqa: E402
    ALL_ENDPOINTS,
    CLASS,
    ENSEMBL,
    IFN_SCORE_MEMBERS,
    MHC1_SCORE_MEMBERS,
    PRIMARY,
)

SEED = 20260817
N_BOOT = 2000
MIN_N = 8
MIN_SCORE_GENES = 4
MIN_MHC1_GENES = 3  # B2M is absent in this freeze; HLA-A/B/C still form a 3-gene score


def find_row(index: pd.Index, prefixes: list[str]) -> str | None:
    hits: list[str] = []
    for prefix in prefixes:
        for i in index:
            s = str(i)
            if s == prefix or s.startswith(prefix + ".") or s.startswith(prefix + "|"):
                hits.append(s)
    hits = list(dict.fromkeys(hits))
    if not hits:
        return None
    if len(hits) > 1:
        # Prefer exact prefix+version over a different gene that shares a prefix.
        versioned = [h for h in hits if any(h == p or h.startswith(p + ".") for p in prefixes)]
        if len(versioned) == 1:
            return versioned[0]
        raise ValueError(f"multiple rows for {prefixes}: {hits}")
    return hits[0]


def load_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    return df.apply(pd.to_numeric, errors="coerce")


def extract_gene(mat: pd.DataFrame, symbol: str) -> tuple[pd.Series, str | None]:
    row = find_row(mat.index, ENSEMBL[symbol])
    if row is None:
        return pd.Series(np.nan, index=mat.columns, name=symbol), None
    s = pd.to_numeric(mat.loc[row], errors="coerce")
    s.name = symbol
    return s, row


def mean_z(series_list: list[pd.Series], min_genes: int) -> pd.Series:
    if not series_list:
        return pd.Series(dtype=float)
    rows = []
    for s in series_list:
        sd = s.std(ddof=0)
        if sd == 0 or not np.isfinite(sd):
            continue
        rows.append((s - s.mean()) / sd)
    if not rows:
        return pd.Series(dtype=float)
    z = pd.concat(rows, axis=1)
    n = z.notna().sum(axis=1)
    score = z.mean(axis=1, skipna=True)
    score[n < min_genes] = np.nan
    score.name = "score"
    return score


def spearman_boot(x: pd.Series, y: pd.Series, seed: int = SEED, n_boot: int = N_BOOT) -> dict:
    d = pd.concat([x, y], axis=1).dropna()
    n = int(len(d))
    rec = {
        "n": n,
        "rho": np.nan,
        "p": np.nan,
        "ci_lo": np.nan,
        "ci_hi": np.nan,
    }
    if n < MIN_N or d.iloc[:, 0].nunique() < 2 or d.iloc[:, 1].nunique() < 2:
        return rec
    rho, p = stats.spearmanr(d.iloc[:, 0], d.iloc[:, 1])
    rec["rho"] = float(rho)
    rec["p"] = float(p)
    rng = np.random.default_rng(seed)
    xv = d.iloc[:, 0].to_numpy()
    yv = d.iloc[:, 1].to_numpy()
    boots = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        r, _ = stats.spearmanr(xv[idx], yv[idx])
        if np.isfinite(r):
            boots.append(float(r))
    if boots:
        rec["ci_lo"] = float(np.percentile(boots, 2.5))
        rec["ci_hi"] = float(np.percentile(boots, 97.5))
    return rec


def cohort_protein(data: Path, cohort: str) -> Path:
    return (
        data
        / cohort
        / f"{cohort}_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"
    )


def analyze_cohort(mat: pd.DataFrame, cohort: str) -> tuple[pd.DataFrame, list[dict], list[dict], dict]:
    cldn4, cldn4_row = extract_gene(mat, "CLDN4")
    n_tumors = int(mat.shape[1])
    n_cldn4 = int(cldn4.notna().sum())

    coverage = []
    series: dict[str, pd.Series] = {"CLDN4": cldn4}
    rows_used: dict[str, str | None] = {"CLDN4": cldn4_row}

    for symbol in ALL_ENDPOINTS:
        s, row = extract_gene(mat, symbol)
        series[symbol] = s
        rows_used[symbol] = row
        n_gene = int(s.notna().sum())
        pair = pd.concat([cldn4, s], axis=1).dropna()
        coverage.append(
            {
                "cohort": cohort,
                "symbol": symbol,
                "class": CLASS[symbol],
                "ensembl_query": ";".join(ENSEMBL[symbol]),
                "row": row or "",
                "present": bool(row is not None),
                "n_quantified": n_gene,
                "n_pairwise_with_cldn4": int(len(pair)),
                "usable_n_ge_8": bool(row is not None and len(pair) >= MIN_N),
            }
        )

    corrs = []
    for symbol in ALL_ENDPOINTS:
        cov = next(c for c in coverage if c["symbol"] == symbol)
        rec = {
            "cohort": cohort,
            "predictor": "CLDN4_protein",
            "endpoint": f"{symbol}_protein",
            "symbol": symbol,
            "class": CLASS[symbol],
            "row": cov["row"],
            "present": cov["present"],
            "usable": cov["usable_n_ge_8"],
        }
        if cov["usable_n_ge_8"]:
            rec.update(spearman_boot(cldn4, series[symbol]))
        else:
            rec.update({"n": int(cov["n_pairwise_with_cldn4"]), "rho": np.nan, "p": np.nan, "ci_lo": np.nan, "ci_hi": np.nan})
        corrs.append(rec)

    # Extra scores from present members only.
    for score_name, members, min_genes in (
        ("IFN_core_protein", IFN_SCORE_MEMBERS, MIN_SCORE_GENES),
        ("MHC1_protein", MHC1_SCORE_MEMBERS, MIN_MHC1_GENES),
    ):
        used = []
        vecs = []
        for g in members:
            cov = next(c for c in coverage if c["symbol"] == g)
            if cov["usable_n_ge_8"]:
                used.append(g)
                vecs.append(series[g])
        score = mean_z(vecs, min_genes)
        pair_n = int(pd.concat([cldn4, score], axis=1).dropna().shape[0]) if len(score) else 0
        rec = {
            "cohort": cohort,
            "predictor": "CLDN4_protein",
            "endpoint": score_name,
            "symbol": score_name,
            "class": "score",
            "row": ",".join(used),
            "present": bool(len(used) >= min_genes),
            "usable": bool(len(used) >= min_genes and pair_n >= MIN_N),
            "members_used": ",".join(used),
            "n_members_used": len(used),
            "n_members_requested": len(members),
        }
        if rec["usable"]:
            rec.update(spearman_boot(cldn4, score))
            series[score_name] = score
        else:
            rec.update({"n": pair_n, "rho": np.nan, "p": np.nan, "ci_lo": np.nan, "ci_hi": np.nan})
        corrs.append(rec)

    core = pd.DataFrame({"CLDN4_protein": cldn4})
    core["cohort"] = cohort
    for symbol in PRIMARY + ["STAT1", "ISG15", "MX1", "HLA-A"]:
        if symbol in series:
            core[f"{symbol}_protein"] = series[symbol]
    for score_name in ("IFN_core_protein", "MHC1_protein"):
        if score_name in series:
            core[score_name] = series[score_name]
    core.index.name = "case_id"

    info = {
        "cohort": cohort,
        "n_protein_columns": n_tumors,
        "n_protein_rows": int(mat.shape[0]),
        "cldn4_row": cldn4_row,
        "n_cldn4_protein": n_cldn4,
        "n_cldn4_protein_na": n_tumors - n_cldn4,
        "cd274_row": rows_used["CD274"],
        "n_cd274_protein": int(series["CD274"].notna().sum()),
        "n_cldn4_cd274_pairwise": int(pd.concat([cldn4, series["CD274"]], axis=1).dropna().shape[0]),
    }
    return core, coverage, corrs, info


def style(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=8)


def fmt_p(p: float) -> str:
    if not np.isfinite(p):
        return "NA"
    if p < 1e-4:
        return f"{p:.1e}"
    return f"{p:.3g}"


def plot_scatter(cores: dict[str, pd.DataFrame], corrs: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.8))
    for ax, cohort in zip(axes, ("LUAD", "LSCC")):
        d = cores[cohort]
        x = d["CLDN4_protein"]
        y = d["CD274_protein"]
        pair = pd.concat([x, y], axis=1).dropna()
        ax.scatter(pair.iloc[:, 0], pair.iloc[:, 1], s=18, c="#2c5f8a", alpha=0.75, edgecolors="none")
        row = corrs[(corrs.cohort == cohort) & (corrs.symbol == "CD274")].iloc[0]
        ax.set_xlabel("CLDN4 protein (log2 TMT)", fontsize=8)
        ax.set_ylabel("CD274 protein (log2 TMT)", fontsize=8)
        rho = row["rho"]
        n = int(row["n"])
        p = row["p"]
        title = f"{cohort}  n={n}"
        if np.isfinite(rho):
            title += f"  ρ={rho:+.3f}  p={fmt_p(p)}"
        else:
            title += "  CD274 not usable"
        ax.set_title(title, fontsize=9)
        style(ax)
    fig.suptitle(
        "Extra scatter · CLDN4 vs CD274 protein · CPTAC public TMT freeze v1.2",
        fontsize=11,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_ifn_extra(cores: dict[str, pd.DataFrame], corrs: pd.DataFrame, path: Path) -> None:
    """Extra IFN-axis scatters for proteins that are actually present."""
    want = []
    for symbol in ("STAT1", "ISG15", "MX1", "HLA-A"):
        ok = False
        for cohort in ("LUAD", "LSCC"):
            sub = corrs[(corrs.cohort == cohort) & (corrs.symbol == symbol)]
            if len(sub) and bool(sub.iloc[0]["usable"]):
                ok = True
        if ok:
            want.append(symbol)
    if not want:
        return
    fig, axes = plt.subplots(2, len(want), figsize=(3.1 * len(want), 6.0), squeeze=False)
    for i, cohort in enumerate(("LUAD", "LSCC")):
        d = cores[cohort]
        for j, symbol in enumerate(want):
            ax = axes[i, j]
            col = f"{symbol}_protein"
            if col not in d.columns:
                ax.set_axis_off()
                continue
            pair = pd.concat([d["CLDN4_protein"], d[col]], axis=1).dropna()
            ax.scatter(pair.iloc[:, 0], pair.iloc[:, 1], s=14, c="#2c5f8a", alpha=0.7, edgecolors="none")
            row = corrs[(corrs.cohort == cohort) & (corrs.symbol == symbol)]
            rho = float(row.iloc[0]["rho"]) if len(row) else np.nan
            p = float(row.iloc[0]["p"]) if len(row) else np.nan
            n = int(row.iloc[0]["n"]) if len(row) else len(pair)
            ax.set_xlabel("CLDN4 protein", fontsize=7)
            ax.set_ylabel(f"{symbol} protein", fontsize=7)
            lab = f"{cohort} {symbol}\nn={n}"
            if np.isfinite(rho):
                lab += f"  ρ={rho:+.2f} p={fmt_p(p)}"
            ax.set_title(lab, fontsize=8)
            style(ax)
    fig.suptitle("Extra scatter · CLDN4 vs present IFN-axis proteins", fontsize=11)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/cptac_cldn4_cd274")
    ap.add_argument("--outdir", default="methods/cptac_cldn4_cd274")
    args = ap.parse_args()
    data = Path(args.data)
    out = Path(args.outdir)
    tables = out / "tables"
    figs = out / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figs.mkdir(parents=True, exist_ok=True)

    cores = {}
    cov_all = []
    corr_all = []
    infos = []
    for cohort in ("LUAD", "LSCC"):
        path = cohort_protein(data, cohort)
        if not path.exists():
            print(f"missing {path}", file=sys.stderr)
            return 1
        mat = load_matrix(path)
        core, coverage, corrs, info = analyze_cohort(mat, cohort)
        cores[cohort] = core
        cov_all.extend(coverage)
        corr_all.extend(corrs)
        infos.append(info)

    sample = pd.concat(cores.values(), axis=0)
    sample.to_csv(tables / "sample_scores.tsv", sep="\t")
    cov = pd.DataFrame(cov_all)
    corr = pd.DataFrame(corr_all)
    cov.to_csv(tables / "coverage.tsv", sep="\t", index=False)
    corr.to_csv(tables / "correlations.tsv", sep="\t", index=False)

    n_rows = []
    for info in infos:
        n_rows.append(
            {
                "cohort": info["cohort"],
                "n_tumors": info["n_protein_columns"],
                "n_cldn4_protein": info["n_cldn4_protein"],
                "n_cldn4_protein_na": info["n_cldn4_protein_na"],
                "n_cd274_protein": info["n_cd274_protein"],
                "n_cldn4_vs_cd274": info["n_cldn4_cd274_pairwise"],
                "cldn4_row": info["cldn4_row"],
                "cd274_row": info["cd274_row"],
            }
        )
    ntab = pd.DataFrame(n_rows)
    ntab.to_csv(tables / "n_table.tsv", sep="\t", index=False)

    plot_scatter(cores, corr, figs / "fig_cldn4_vs_cd274.png")
    plot_ifn_extra(cores, corr, figs / "fig_cldn4_vs_ifn_proteins.png")

    summary = {
        "question": (
            "CLDN4 protein vs CD274 protein and locked IFN proteins "
            "in CPTAC LUAD and LSCC public TMT. Additive; protein–protein only."
        ),
        "predictor": "CLDN4 protein (ENSG00000189143)",
        "primary_endpoint": "CD274 protein (ENSG00000120217)",
        "test": "Spearman, two-sided, pairwise-complete; bootstrap 95% CI 2000 resamples seed 20260817",
        "min_n": MIN_N,
        "cohorts": infos,
        "n_table": n_rows,
        "correlations": corr.to_dict(orient="records"),
        "coverage": cov.to_dict(orient="records"),
    }
    (tables / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    primary = corr[corr.symbol.isin(PRIMARY + ["IFN_core_protein", "MHC1_protein"])]
    print(ntab.to_string(index=False))
    print(primary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
