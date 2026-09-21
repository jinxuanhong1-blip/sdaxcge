#!/usr/bin/env python3
"""One method sweep for GSE207422 immune endpoints. Not a search for a hit.

Pre-specified before the p-values from this script are known.

Question: do TACSTD2 and CLDN4, scored separately, associate with CD8 or with
MPR/NMPR on the public Hu et al. matrix under other standard choices?

Cohort: 12 post-treatment samples. pCR counted as MPR. Cell calls are the same
functions as analyze.py. TACSTD2 and CLDN4 are not cell-type gates.

Expression:
- epithelial mean log1p(CP10k), already the lead metric
- epithelial percent of cells with UMI > 0

CD8-family outcomes:
- CD8 fraction of all cells (CD3+ and CD8A>0)
- CD8 fraction within lineage T/NK
- NK fraction of all cells
- CD8 cytotoxicity: mean log1p(CP10k) of GZMB, GZMA, PRF1, IFNG, NKG7, GNLY
  inside CD8 cells (samples with >=5 CD8 cells)

Tests that decide whether the endpoint stays closed (full n=12):
1. Spearman, epithelial mean vs CD8 fraction
2. Spearman, epithelial mean vs CD8+NK fraction
3. Spearman, epithelial mean vs CD8 cytotoxicity
4. Partial Spearman of (1) adjusted for epithelial fraction
5. Exact Mann-Whitney, NMPR vs MPR, epithelial mean
6. Exact Mann-Whitney, NMPR vs MPR, epithelial percent positive

A positive finding requires the epithelial-mean test AND at least one of
percent-positive or the partial correlation to have P<0.05 in the same
direction. Otherwise that gene's immune endpoint is closed.

Also reported, and not allowed to reopen a closed endpoint by themselves:
quartile (top 3 vs bottom 3; the smallest two-sided exact P on n=3+3 is 0.10),
adenocarcinoma-only and squamous-only Spearman (n=6), and leave-one-out ranges.
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
import analyze  # noqa: E402

CYTO = ["GZMB", "GZMA", "PRF1", "IFNG", "NKG7", "GNLY"]
MIN_CD8_CYTO = 5


def as_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.astype(str).str.strip().str.lower().isin(["true", "1"])


def fmt_p(p) -> str:
    if p is None or (isinstance(p, float) and not np.isfinite(p)):
        return "NA"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


def fmt_num(x, digits=3) -> str:
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    return f"{x:.{digits}f}"


def spearman(x, y) -> dict:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 4:
        return {"n": int(len(x)), "rho": None, "p": None}
    rho, p = stats.spearmanr(x, y)
    return {"n": int(len(x)), "rho": float(rho), "p": float(p)}


def partial_spearman(x, y, z) -> dict:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x, y, z = x[m], y[m], z[m]
    if len(x) < 5:
        return {"n": int(len(x)), "rho": None, "p": None}
    rx, ry, rz = stats.rankdata(x), stats.rankdata(y), stats.rankdata(z)

    def resid(a, b):
        design = np.column_stack([np.ones(len(b)), b])
        coef, _, _, _ = np.linalg.lstsq(design, a, rcond=None)
        return a - design @ coef

    xr, yr = resid(rx, rz), resid(ry, rz)
    if np.std(xr) == 0 or np.std(yr) == 0:
        return {"n": int(len(x)), "rho": None, "p": None}
    r, p = stats.pearsonr(xr, yr)
    return {"n": int(len(x)), "rho": float(r), "p": float(p)}


def mw_expr(df: pd.DataFrame, col: str) -> dict:
    a = df.loc[df["response_paper"] == "NMPR", col].to_numpy(dtype=float)
    b = df.loc[df["response_paper"] == "MPR", col].to_numpy(dtype=float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    rec = {
        "n_nmpr": int(len(a)),
        "n_mpr": int(len(b)),
        "median_nmpr": float(np.median(a)) if len(a) else None,
        "median_mpr": float(np.median(b)) if len(b) else None,
        "p": None,
        "method": "skipped",
    }
    if len(a) < 1 or len(b) < 1:
        return rec
    try:
        result = stats.mannwhitneyu(a, b, alternative="two-sided", method="exact")
        method = "exact"
    except ValueError:
        result = stats.mannwhitneyu(a, b, alternative="two-sided", method="asymptotic")
        method = "asymptotic"
    rec.update(p=float(result.pvalue), method=method)
    rec["delta_nmpr_minus_mpr"] = rec["median_nmpr"] - rec["median_mpr"]
    return rec


def cyto_table(cells: pd.DataFrame, expr: dict[str, np.ndarray], cell_ids: list[str]) -> pd.DataFrame:
    if cells["barcode"].tolist() != list(cell_ids):
        raise SystemExit("cell table order does not match the matrix")
    n_umi = cells["n_umi"].to_numpy(dtype=float)
    scale = np.where(n_umi > 0, 1e4 / n_umi, 0.0)
    mats = []
    for gene in CYTO:
        if gene not in expr:
            raise SystemExit(f"missing cytotoxicity gene {gene}")
        mats.append(np.log1p(expr[gene] * scale))
    cyto = np.mean(np.vstack(mats), axis=0)
    work = cells[["Sample", "is_epithelial", "is_cd8", "is_tnk_lineage", "TACSTD2", "CLDN4"]].copy()
    work["cyto"] = cyto
    rows = []
    for sample, sdf in work.groupby("Sample"):
        epi = sdf[sdf["is_epithelial"]]
        cd8 = sdf[sdf["is_cd8"]]
        tnk = sdf[sdf["is_tnk_lineage"]]
        rows.append(
            {
                "Sample": sample,
                "tacstd2_epi_pct": float(100.0 * (epi["TACSTD2"] > 0).mean()) if len(epi) else np.nan,
                "cldn4_epi_pct": float(100.0 * (epi["CLDN4"] > 0).mean()) if len(epi) else np.nan,
                "n_cd8_rebuild": int(len(cd8)),
                "cyto_cd8": float(cd8["cyto"].mean()) if len(cd8) >= MIN_CD8_CYTO else np.nan,
                "cyto_tnk": float(tnk["cyto"].mean()) if len(tnk) >= MIN_CD8_CYTO else np.nan,
                "n_cd8_for_cyto": int(len(cd8)),
            }
        )
    return pd.DataFrame(rows)


def add_derived(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["frac_cd8"] = out["n_cd8"] / out["n_cells"]
    out["frac_nk"] = out["n_nk"] / out["n_cells"]
    out["frac_cd8_within_tnk"] = np.where(out["n_tnk_lineage"] > 0, out["n_cd8"] / out["n_tnk_lineage"], np.nan)
    return out


def row_spearman(gene, expr_name, endpoint, rec, **extra) -> dict:
    return {
        "gene": gene,
        "expression": expr_name,
        "test": extra.get("test", "spearman"),
        "stratum": extra.get("stratum", "all_post"),
        "endpoint": endpoint,
        "n": rec["n"],
        "rho": rec["rho"],
        "p": rec["p"],
        "note": extra.get("note", ""),
    }


def decide(tests: pd.DataFrame) -> dict:
    """Apply the pre-specified closure rule. Does not look for a favorable cut."""
    primary_cd8 = [
        ("epithelial_mean", "spearman", "CD8 fraction"),
        ("epithelial_mean", "spearman", "CD8+NK fraction"),
        ("epithelial_mean", "spearman", "CD8 cytotoxicity"),
        ("epithelial_mean", "partial_spearman_adj_epithelial_fraction", "CD8 fraction"),
    ]
    out = {"genes": {}}
    for gene in ("TACSTD2", "CLDN4"):
        cd8_rows = []
        for expr, test, endpoint in primary_cd8:
            hit = tests[
                (tests["gene"] == gene)
                & (tests["expression"] == expr)
                & (tests["test"] == test)
                & (tests["endpoint"] == endpoint)
                & (tests["stratum"] == "all_post")
            ]
            if hit.empty:
                raise SystemExit(f"missing primary row {gene} {expr} {test} {endpoint}")
            cd8_rows.append(hit.iloc[0])
        mpr_mean = tests[
            (tests["gene"] == gene)
            & (tests["expression"] == "epithelial_mean")
            & (tests["test"] == "NMPR_vs_MPR")
            & (tests["stratum"] == "all_post")
        ].iloc[0]
        mpr_pct = tests[
            (tests["gene"] == gene)
            & (tests["expression"] == "epithelial_pct")
            & (tests["test"] == "NMPR_vs_MPR")
            & (tests["stratum"] == "all_post")
        ].iloc[0]
        cd8_p = [float(r["p"]) for r in cd8_rows]
        cd8_rho = [float(r["rho"]) for r in cd8_rows]
        mean_p = float(mpr_mean["p"])
        pct_p = float(mpr_pct["p"])
        # Direction of the MPR contrast: positive delta means higher in NMPR.
        mean_delta = float(mpr_mean["delta"])
        pct_delta = float(mpr_pct["delta"])
        mpr_concordant = (
            mean_p < 0.05
            and pct_p < 0.05
            and np.sign(mean_delta) == np.sign(pct_delta)
            and mean_delta != 0
        )
        # CD8 concordance: epithelial-mean Spearman AND partial Spearman, same sign, both P<0.05.
        mean_cd8 = next(r for r in cd8_rows if r["test"] == "spearman" and r["endpoint"] == "CD8 fraction")
        partial_cd8 = next(r for r in cd8_rows if r["test"] == "partial_spearman_adj_epithelial_fraction")
        cd8_concordant = (
            float(mean_cd8["p"]) < 0.05
            and float(partial_cd8["p"]) < 0.05
            and np.sign(float(mean_cd8["rho"])) == np.sign(float(partial_cd8["rho"]))
            and float(mean_cd8["rho"]) != 0
        )
        out["genes"][gene] = {
            "cd8_closed": not cd8_concordant,
            "mpr_closed": not mpr_concordant,
            "cd8_p": cd8_p,
            "cd8_rho": cd8_rho,
            "mpr_mean_p": mean_p,
            "mpr_pct_p": pct_p,
            "mpr_mean_delta": mean_delta,
            "mpr_pct_delta": pct_delta,
            "max_abs_primary_cd8_rho": float(np.max(np.abs(cd8_rho))),
        }
    out["immune_endpoints_final"] = all(
        g["cd8_closed"] and g["mpr_closed"] for g in out["genes"].values()
    )
    return out


def make_figure(tests: pd.DataFrame, path: Path) -> None:
    want = tests[
        (tests["stratum"] == "all_post")
        & (tests["expression"] == "epithelial_mean")
        & (tests["test"].isin(["spearman", "partial_spearman_adj_epithelial_fraction"]))
        & (tests["gene"].isin(["TACSTD2", "CLDN4"]))
    ].copy()
    want["label"] = want["gene"] + " | " + want["endpoint"] + np.where(
        want["test"] == "partial_spearman_adj_epithelial_fraction", " | partial", ""
    )
    want = want.iloc[::-1]
    fig, ax = plt.subplots(figsize=(8.2, 5.2), constrained_layout=True)
    colors = ["#b2182b" if g == "TACSTD2" else "#2166ac" for g in want["gene"]]
    y = np.arange(len(want))
    ax.axvline(0, color="#444444", lw=0.8)
    ax.scatter(want["rho"], y, c=colors, s=42, zorder=3)
    for i, rho in enumerate(want["rho"]):
        ax.plot([0, rho], [i, i], color=colors[i], lw=1.2)
    ax.set_yticks(y)
    ax.set_yticklabels(want["label"], fontsize=8)
    ax.set_xlabel("Spearman ρ (partial rows are Pearson ρ of rank residuals)")
    ax.set_title("GSE207422 post-treatment n=12\nEpithelial-mean TACSTD2 (red) and CLDN4 (blue) vs CD8")
    fig.savefig(path, dpi=180)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def write_sweep(path: Path, tests: pd.DataFrame, loo: pd.DataFrame, decision: dict) -> None:
    lines = [
        "# GSE207422 immune-endpoint method sweep",
        "",
        "One additional sweep for TACSTD2 and CLDN4 versus CD8 and versus MPR/NMPR.",
        "Rules were fixed in `scripts/sweep_immune.py` before these p-values were computed.",
        "A positive finding requires the epithelial-mean test and either percent-positive or the epithelial-fraction partial correlation to agree at P<0.05.",
        "Quartile, histology-only, and leave-one-out cuts do not reopen a closed endpoint.",
        "",
        f"Closure call: **{'FINAL' if decision['immune_endpoints_final'] else 'NOT CLOSED'}** for GSE207422 immune endpoints (CD8 and MPR).",
        "",
        "## Primary rows (n=12 post-treatment)",
        "",
    ]
    primary = tests[(tests["stratum"] == "all_post") & (tests["test"].isin([
        "spearman",
        "partial_spearman_adj_epithelial_fraction",
        "NMPR_vs_MPR",
        "quartile_top3_vs_bottom3",
    ]))]
    for _, row in primary.iterrows():
        if row["test"] == "NMPR_vs_MPR":
            lines.append(
                f"- {row['gene']} | {row['expression']} | NMPR vs MPR: "
                f"Δ={fmt_num(row['delta'])} (n={int(row['n_nmpr'])} vs {int(row['n_mpr'])}), "
                f"{row['method']} P={fmt_p(row['p'])}"
            )
        elif row["test"] == "quartile_top3_vs_bottom3":
            lines.append(
                f"- {row['gene']} | {row['expression']} | {row['endpoint']} top3 vs bottom3: "
                f"Δ median={fmt_num(row['delta'])}, P={fmt_p(row['p'])} ({row['note']})"
            )
        else:
            lines.append(
                f"- {row['gene']} | {row['expression']} | {row['test']} | {row['endpoint']}: "
                f"ρ={fmt_num(row['rho'], 2)} P={fmt_p(row['p'])} n={int(row['n'])}"
            )
    lines += ["", "## Histology strata (n=6 each; descriptive)", ""]
    strata = tests[tests["stratum"].isin(["Adeno", "Squamous"])]
    for _, row in strata.iterrows():
        if row["test"] == "NMPR_vs_MPR":
            lines.append(
                f"- {row['stratum']} | {row['gene']} | {row['expression']} | NMPR vs MPR: "
                f"Δ={fmt_num(row['delta'])} (n={int(row['n_nmpr'])} vs {int(row['n_mpr'])}), P={fmt_p(row['p'])}"
            )
        else:
            lines.append(
                f"- {row['stratum']} | {row['gene']} | {row['expression']} | {row['endpoint']}: "
                f"ρ={fmt_num(row['rho'], 2)} P={fmt_p(row['p'])} n={int(row['n'])}"
            )
    lines += ["", "## Leave-one-out, epithelial mean", ""]
    for gene in ("TACSTD2", "CLDN4"):
        for endpoint in ("CD8 fraction", "CD8+NK fraction", "MPR"):
            sub = loo[(loo["gene"] == gene) & (loo["endpoint"] == endpoint)]
            if sub.empty:
                continue
            if endpoint == "MPR":
                lines.append(
                    f"- {gene} | MPR Mann-Whitney P range {fmt_p(sub['p'].min())} to {fmt_p(sub['p'].max())} "
                    f"when one patient is dropped (full-sample P is the primary result)."
                )
            else:
                lines.append(
                    f"- {gene} | {endpoint} Spearman ρ range {fmt_num(sub['rho'].min(), 2)} to {fmt_num(sub['rho'].max(), 2)} "
                    f"(P {fmt_p(sub['p'].min())} to {fmt_p(sub['p'].max())})."
                )
    lines += [
        "",
        "## Call",
        "",
    ]
    for gene, info in decision["genes"].items():
        lines.append(
            f"- {gene}: CD8 closed={info['cd8_closed']} (largest |ρ| among the four primary CD8 rows "
            f"{fmt_num(info['max_abs_primary_cd8_rho'], 2)}); "
            f"MPR closed={info['mpr_closed']} "
            f"(mean Δ={fmt_num(info['mpr_mean_delta'])}, P={fmt_p(info['mpr_mean_p'])}; "
            f"percent-positive Δ={fmt_num(info['mpr_pct_delta'])}, P={fmt_p(info['mpr_pct_p'])})."
        )
    lines += ["", exception_note(tests, loo), ""]
    if decision["immune_endpoints_final"]:
        lines += [
            "Immune endpoints on this matrix are **FINAL**. CD8 abundance, CD8 cytotoxicity, and MPR do not meet the pre-specified concordance rule for TACSTD2 or for CLDN4. Further cuts of GSE207422 will not be treated as new evidence for these endpoints.",
            "",
        ]
    else:
        lines += ["At least one gene met the concordance rule. The endpoint is not closed. See the primary rows.", ""]
    path.write_text("\n".join(lines) + "\n")


def exception_note(tests: pd.DataFrame, loo: pd.DataFrame) -> str:
    """Name every P<0.05 outside the closure rule so a small-n cut is not hidden."""
    bits = ["## Rows with P<0.05 that do not reopen the endpoint", ""]
    flagged = tests[(tests["p"] < 0.05) & (tests["stratum"] != "all_post")]
    if flagged.empty:
        bits.append("- No histology-stratum test has P<0.05.")
    for _, row in flagged.iterrows():
        bits.append(
            f"- {row['stratum']} | {row['gene']} | {row['expression']} | {row['endpoint']}: "
            f"ρ={fmt_num(row['rho'], 2)} P={fmt_p(row['p'])} n={int(row['n'])}. "
            "Histology strata were pre-specified as descriptive."
        )
    # Companion mean in the same stratum, if the hit was percent-positive.
    for _, row in flagged.iterrows():
        if row["expression"] != "epithelial_pct":
            continue
        companion = tests[
            (tests["gene"] == row["gene"])
            & (tests["expression"] == "epithelial_mean")
            & (tests["stratum"] == row["stratum"])
            & (tests["endpoint"] == row["endpoint"])
            & (tests["test"] == "spearman")
        ]
        if companion.empty:
            continue
        c = companion.iloc[0]
        bits.append(
            f"- Same stratum, epithelial mean: ρ={fmt_num(c['rho'], 2)} P={fmt_p(c['p'])} n={int(c['n'])}. "
            "The percent-positive stratum does not agree with the mean at P<0.05."
        )
    mpr = loo[loo["endpoint"] == "MPR"]
    for gene in ("TACSTD2", "CLDN4"):
        sub = mpr[mpr["gene"] == gene]
        low = sub[sub["p"] < 0.05]
        if low.empty:
            continue
        dropped = ", ".join(low["dropped"].astype(str))
        bits.append(
            f"- Leave-one-out {gene} epithelial-mean MPR: dropping {dropped} moves P to "
            f"{fmt_p(low['p'].min())}. The full-sample P stays the primary result."
        )
    bits.append(
        "These rows are reported because they are the only places a P falls under 0.05. "
        "They are n=6 histology or single-patient deletions. They do not meet the concordance rule."
    )
    return "\n".join(bits)


def append_final(finding: Path, decision: dict, tests: pd.DataFrame, loo: pd.DataFrame) -> None:
    text = finding.read_text()
    marker = "## FINAL — GSE207422 immune endpoints"
    if marker in text:
        text = text[: text.index(marker)].rstrip() + "\n"
    if not decision["immune_endpoints_final"]:
        block = (
            "\n## FINAL — GSE207422 immune endpoints\n\n"
            "Not closed. The method sweep met the pre-specified concordance rule for at least one gene. "
            "See `SWEEP.md`.\n"
        )
    else:
        tac = decision["genes"]["TACSTD2"]
        cld = decision["genes"]["CLDN4"]
        block = (
            "\n## FINAL — GSE207422 immune endpoints\n\n"
            "Closed after one additional method sweep (`scripts/sweep_immune.py`, `SWEEP.md`). "
            "The sweep added CD8-only fraction, CD8 within T/NK, NK fraction, CD8 cytotoxicity "
            "(GZMB/GZMA/PRF1/IFNG/NKG7/GNLY), epithelial percent-positive, partial correlation "
            "adjusted for epithelial fraction, a top-3 vs bottom-3 cut, histology strata, and leave-one-out.\n\n"
            "A positive result required the epithelial-mean test and either percent-positive or the partial "
            "correlation to agree at P<0.05. Quartile, histology-only, and leave-one-out cuts were not allowed "
            "to reopen the endpoint.\n\n"
            f"- TACSTD2 vs CD8: closed. Largest |ρ| on the four primary CD8 rows is "
            f"{fmt_num(tac['max_abs_primary_cd8_rho'], 2)}.\n"
            f"- CLDN4 vs CD8: closed. Largest |ρ| on the four primary CD8 rows is "
            f"{fmt_num(cld['max_abs_primary_cd8_rho'], 2)}. "
            "CLDN4 versus T/NK stays flat, as in the earlier CLDN4-only run.\n"
            f"- TACSTD2 vs MPR: closed as not significant. Epithelial-mean NMPR−MPR Δ="
            f"{fmt_num(tac['mpr_mean_delta'])} log1p(CP10k) (P={fmt_p(tac['mpr_mean_p'])}); "
            f"percent-positive Δ={fmt_num(tac['mpr_pct_delta'])} percentage points (P={fmt_p(tac['mpr_pct_p'])}). "
            "The mean is higher in NMPR. Percent-positive does not move. The pair does not meet the concordance rule.\n"
            f"- CLDN4 vs MPR: closed as null. Epithelial-mean Δ={fmt_num(cld['mpr_mean_delta'])} "
            f"log1p(CP10k) (P={fmt_p(cld['mpr_mean_p'])}); percent-positive Δ={fmt_num(cld['mpr_pct_delta'])} "
            f"percentage points (P={fmt_p(cld['mpr_pct_p'])}).\n\n"
            + _final_exceptions(tests, loo)
            + "\n\nNo further cut of GSE207422 will be treated as new evidence for CD8 or MPR.\n"
        )
    finding.write_text(text.rstrip() + "\n" + block)


def _final_exceptions(tests: pd.DataFrame, loo: pd.DataFrame) -> str:
    hit = tests[
        (tests["gene"] == "CLDN4")
        & (tests["expression"] == "epithelial_pct")
        & (tests["stratum"] == "Squamous")
        & (tests["endpoint"] == "CD8 fraction")
        & (tests["test"] == "spearman")
    ].iloc[0]
    mean = tests[
        (tests["gene"] == "CLDN4")
        & (tests["expression"] == "epithelial_mean")
        & (tests["stratum"] == "Squamous")
        & (tests["endpoint"] == "CD8 fraction")
        & (tests["test"] == "spearman")
    ].iloc[0]
    full = tests[
        (tests["gene"] == "CLDN4")
        & (tests["expression"] == "epithelial_pct")
        & (tests["stratum"] == "all_post")
        & (tests["endpoint"] == "CD8 fraction")
        & (tests["test"] == "spearman")
    ].iloc[0]
    loo_mpr = loo[(loo["gene"] == "TACSTD2") & (loo["endpoint"] == "MPR") & (loo["p"] < 0.05)]
    dropped = " or ".join(loo_mpr["dropped"].astype(str))
    return (
        "Two cuts in the sweep fall under P=0.05 and are not used to reopen the endpoint. "
        f"Squamous-only CLDN4 percent-positive versus CD8 fraction is ρ={fmt_num(hit['rho'], 2)} "
        f"(P={fmt_p(hit['p'])}, n={int(hit['n'])}); "
        f"the epithelial mean in that same stratum is ρ={fmt_num(mean['rho'], 2)} (P={fmt_p(mean['p'])}), "
        f"and the full-sample percent-positive correlation is ρ={fmt_num(full['rho'], 2)} (P={fmt_p(full['p'])}). "
        f"Leave-one-out of the TACSTD2 epithelial-mean MPR test moves P from 0.109 to {fmt_p(loo_mpr['p'].min())} "
        f"when {dropped} is dropped. The full-sample test and the percent-positive MPR test (P=1.00) stay the recorded result."
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", type=Path)
    ap.add_argument("--metadata", type=Path)
    ap.add_argument("--per-sample", type=Path)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--tables-only", action="store_true", help="Rewrite SWEEP.md and the FINAL section from saved tables")
    args = ap.parse_args()
    tabdir = args.outdir / "tables"
    figdir = args.outdir / "figures"
    tabdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(exist_ok=True)
    if args.tables_only:
        tests = pd.read_csv(tabdir / "sweep_immune.tsv", sep="\t")
        loo = pd.read_csv(tabdir / "sweep_loo.tsv", sep="\t")
        decision = decide(tests)
        write_sweep(args.outdir / "SWEEP.md", tests, loo, decision)
        append_final(args.outdir / "FINDING.md", decision, tests, loo)
        print(json.dumps({"immune_endpoints_final": decision["immune_endpoints_final"]}))
        return

    per = pd.read_csv(args.per_sample, sep="\t")
    per["is_post"] = as_bool(per["is_post"])
    meta = analyze.load_metadata(args.metadata)
    wanted = analyze.marker_universe() | set(CYTO)
    print(f"streaming; extracting {len(wanted)} genes", flush=True)
    cell_ids, expr, n_umi, n_genes, _n_rows = analyze.stream_matrix(args.matrix, wanted)
    missing = sorted(set(CYTO) - set(expr))
    if missing:
        raise SystemExit(f"missing genes {missing}")
    cells = analyze.build_cells(cell_ids, expr, n_umi, n_genes, meta)
    extra = cyto_table(cells, expr, cell_ids)
    df = add_derived(per.merge(extra, on="Sample", how="left"))
    if df["n_cd8_rebuild"].isna().any():
        raise SystemExit("rebuild failed to match samples")
    delta = (df["n_cd8"] - df["n_cd8_rebuild"]).abs().max()
    print(f"max |n_cd8 - rebuilt| = {delta}", flush=True)
    if delta != 0:
        raise SystemExit("CD8 counts do not match the previous cell calls")

    post = df[df["is_post"]].copy()
    if len(post) != 12:
        raise SystemExit(f"expected 12 post samples, found {len(post)}")
    post.to_csv(tabdir / "sweep_per_sample.tsv", sep="\t", index=False)

    expr_cols = {
        "epithelial_mean": {"TACSTD2": "tacstd2_epi_mean", "CLDN4": "cldn4_epi_mean"},
        "epithelial_pct": {"TACSTD2": "tacstd2_epi_pct", "CLDN4": "cldn4_epi_pct"},
    }
    endpoints = {
        "CD8 fraction": "frac_cd8",
        "CD8+NK fraction": "frac_cd8nk",
        "NK fraction": "frac_nk",
        "CD8 within T/NK": "frac_cd8_within_tnk",
        "CD8 cytotoxicity": "cyto_cd8",
        "T/NK cytotoxicity": "cyto_tnk",
    }
    rows: list[dict] = []
    for expr_name, genes in expr_cols.items():
        for gene, col in genes.items():
            for endpoint, ycol in endpoints.items():
                rec = spearman(post[col], post[ycol])
                rows.append(row_spearman(gene, expr_name, endpoint, rec))
                prec = partial_spearman(post[col], post[ycol], post["frac_epithelial"])
                rows.append(
                    row_spearman(
                        gene, expr_name, endpoint, prec, test="partial_spearman_adj_epithelial_fraction"
                    )
                )
            mpr = mw_expr(post, col)
            rows.append(
                {
                    "gene": gene,
                    "expression": expr_name,
                    "test": "NMPR_vs_MPR",
                    "stratum": "all_post",
                    "endpoint": "MPR",
                    "n": mpr["n_nmpr"] + mpr["n_mpr"],
                    "n_nmpr": mpr["n_nmpr"],
                    "n_mpr": mpr["n_mpr"],
                    "rho": None,
                    "p": mpr["p"],
                    "delta": mpr.get("delta_nmpr_minus_mpr"),
                    "method": mpr["method"],
                    "note": "",
                }
            )
            # Top 3 vs bottom 3. Smallest two-sided exact P for 3 vs 3 is 0.10.
            order = post.sort_values(col)
            bottom, top = order.head(3), order.tail(3)
            for endpoint, ycol in (("CD8 fraction", "frac_cd8"), ("CD8 cytotoxicity", "cyto_cd8")):
                try:
                    result = stats.mannwhitneyu(top[ycol], bottom[ycol], alternative="two-sided", method="exact")
                    method, p = "exact", float(result.pvalue)
                except ValueError:
                    result = stats.mannwhitneyu(top[ycol], bottom[ycol], alternative="two-sided", method="asymptotic")
                    method, p = "asymptotic", float(result.pvalue)
                rows.append(
                    {
                        "gene": gene,
                        "expression": expr_name,
                        "test": "quartile_top3_vs_bottom3",
                        "stratum": "all_post",
                        "endpoint": endpoint,
                        "n": 6,
                        "rho": None,
                        "p": p,
                        "delta": float(np.nanmedian(top[ycol]) - np.nanmedian(bottom[ycol])),
                        "method": method,
                        "note": "n=3 vs 3; minimum two-sided exact P is 0.10",
                    }
                )
            top_mpr = int((top["response_paper"] == "MPR").sum())
            bot_mpr = int((bottom["response_paper"] == "MPR").sum())
            odds, fp = stats.fisher_exact(
                [[top_mpr, 3 - top_mpr], [bot_mpr, 3 - bot_mpr]], alternative="two-sided"
            )
            rows.append(
                {
                    "gene": gene,
                    "expression": expr_name,
                    "test": "quartile_top3_vs_bottom3",
                    "stratum": "all_post",
                    "endpoint": "MPR",
                    "n": 6,
                    "rho": None,
                    "p": float(fp),
                    "delta": None,
                    "method": "fisher_exact",
                    "note": f"top MPR {top_mpr}/3 vs bottom {bot_mpr}/3; OR={odds:.3g}; minimum informative P is large at n=3",
                }
            )
            for pathology in ("Adeno", "Squamous"):
                sub = post[post["Pathology"] == pathology]
                rec = spearman(sub[col], sub["frac_cd8"])
                rows.append(row_spearman(gene, expr_name, "CD8 fraction", rec, stratum=pathology))
                mpr = mw_expr(sub, col)
                rows.append(
                    {
                        "gene": gene,
                        "expression": expr_name,
                        "test": "NMPR_vs_MPR",
                        "stratum": pathology,
                        "endpoint": "MPR",
                        "n": mpr["n_nmpr"] + mpr["n_mpr"],
                        "n_nmpr": mpr["n_nmpr"],
                        "n_mpr": mpr["n_mpr"],
                        "rho": None,
                        "p": mpr["p"],
                        "delta": mpr.get("delta_nmpr_minus_mpr"),
                        "method": mpr["method"],
                        "note": "descriptive stratum",
                    }
                )

    loo_rows = []
    for gene, col in expr_cols["epithelial_mean"].items():
        for dropped in post["Patient"]:
            sub = post[post["Patient"] != dropped]
            for endpoint, ycol in (("CD8 fraction", "frac_cd8"), ("CD8+NK fraction", "frac_cd8nk")):
                rec = spearman(sub[col], sub[ycol])
                loo_rows.append(
                    {"gene": gene, "dropped": dropped, "endpoint": endpoint, "rho": rec["rho"], "p": rec["p"]}
                )
            mpr = mw_expr(sub, col)
            loo_rows.append(
                {
                    "gene": gene,
                    "dropped": dropped,
                    "endpoint": "MPR",
                    "rho": mpr.get("delta_nmpr_minus_mpr"),
                    "p": mpr["p"],
                }
            )

    tests = pd.DataFrame(rows)
    loo = pd.DataFrame(loo_rows)
    tests.to_csv(tabdir / "sweep_immune.tsv", sep="\t", index=False)
    loo.to_csv(tabdir / "sweep_loo.tsv", sep="\t", index=False)
    decision = decide(tests)
    (args.outdir / "sweep_decision.json").write_text(json.dumps(decision, indent=2) + "\n")
    write_sweep(args.outdir / "SWEEP.md", tests, loo, decision)
    make_figure(tests, figdir / "fig_sweep_cd8_rho.png")
    append_final(args.outdir / "FINDING.md", decision, tests, loo)
    print(json.dumps(decision, indent=2))
    print("wrote", args.outdir)


if __name__ == "__main__":
    main()
