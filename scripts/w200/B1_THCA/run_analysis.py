#!/usr/bin/env python3
"""B1 analog in TCGA-THCA: TACSTD2 (TROP2) × CLDN4 surface co-expression rank.

Claim B1 (user PPT): among cell-surface genes, CLDN4 is the *top* co-expression
partner of TACSTD2. This script asks the same question in TCGA-THCA only, and
also reports each gene's rank by surface-expression level.

Primary readout (matches B1_BRCA / claim B1):
  Rank every SURFY surfaceome gene by Spearman correlation with TACSTD2
  across primary tumors. Report where CLDN4 actually lands. Do not pre-select
  a favourable neighbourhood.

Secondary readouts (PAAD-style honesty checks + original "surface rank"):
  * TACSTD2–CLDN4 Pearson / Spearman with bootstrap 95% CI
  * Rank-based partial Spearman given ABSOLUTE purity
  * Classic papillary (8260/3) vs follicular-variant (8340/3) sensitivity
  * Surfaceome expression rank (median log2(TPM+1) among surface genes)

Inputs: scripts/w200/B1_THCA/download_data.py → data/raw/ (or B1_THCA_DATA).
Outputs: results/w200/B1_THCA/
"""
from __future__ import annotations

import hashlib
import json
import os
import warnings
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore", message="Unknown extension is not supported")

ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = Path(os.environ.get("B1_THCA_DATA", ROOT / "data" / "raw"))
OUT_DIR = ROOT / "results" / "w200" / "B1_THCA"

ANCHOR = "TACSTD2"
FOCUS = "CLDN4"
WINDOW = 200
RNG = np.random.default_rng(20260816)

# Morphology codes used in the histology sensitivity split.
CLASSIC_PAPILLARY = "8260/3"
FOLLICULAR_VARIANT = "8340/3"


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def one_per_patient(barcodes: list[str]) -> list[str]:
    """Keep the lexicographically first -01 vial per patient."""
    first: dict[str, str] = {}
    for c in sorted(barcodes):
        first.setdefault(c[:12], c)
    return sorted(first.values())


def spearman_ci(x: np.ndarray, y: np.ndarray, n_boot: int = 2000) -> tuple[float, float]:
    n = len(x)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = RNG.integers(0, n, n)
        boots[i] = stats.spearmanr(x[idx], y[idx]).statistic
    return float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))


def partial_spearman(x, y, z) -> tuple[float, float]:
    """Spearman partial correlation of x,y given z (residualised ranks)."""
    rx, ry, rz = (stats.rankdata(v) for v in (x, y, z))

    def resid(a, b):
        beta = np.polyfit(b, a, 1)
        return a - np.polyval(beta, b)

    r, p = stats.pearsonr(resid(rx, rz), resid(ry, rz))
    return float(r), float(p)


def corr_vs_vector(mat: np.ndarray, vec: np.ndarray) -> np.ndarray:
    vc = vec - vec.mean()
    mc = mat - mat.mean(axis=1, keepdims=True)
    num = mc @ vc
    den = np.sqrt((mc**2).sum(axis=1) * (vc**2).sum())
    with np.errstate(invalid="ignore", divide="ignore"):
        return num / den


def pvalues(r: np.ndarray, n: int) -> np.ndarray:
    r = np.clip(r, -0.999999999, 0.999999999)
    with np.errstate(invalid="ignore", divide="ignore"):
        t = r * np.sqrt((n - 2) / (1 - r**2))
    return 2 * stats.t.sf(np.abs(t), df=n - 2)


def bh_fdr(p: np.ndarray) -> np.ndarray:
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


def load_expression() -> tuple[pd.DataFrame, list[str]]:
    path = DATA_DIR / "TCGA-THCA.star_tpm.tsv.gz"
    expr = pd.read_csv(path, sep="\t", index_col=0)
    observed_min = float(np.nanmin(expr.to_numpy()))
    # log2(TPM+1) has floor 0; log2(TPM+0.001) would sit near -9.97.
    if observed_min < -1.0:
        raise SystemExit(
            f"Unexpected value floor {observed_min:.4f}; "
            "this script expects Xena GDC star_tpm as log2(TPM+1)."
        )

    pm = pd.read_csv(DATA_DIR / "gencode.v36.probemap", sep="\t")
    id2sym = pm.set_index("id")["gene"]
    expr.index = id2sym.reindex(expr.index).values
    expr = expr[~expr.index.isna()]
    # Collapse duplicate symbols by keeping the highest-mean row.
    means = expr.mean(axis=1)
    expr = expr.iloc[np.argsort(-means.to_numpy(), kind="mergesort")]
    expr = expr.loc[~expr.index.duplicated()].copy()

    tumor_cols = [c for c in expr.columns if len(c) >= 15 and c[13:15] == "01"]
    tumor_cols = one_per_patient(tumor_cols)
    if not tumor_cols:
        raise SystemExit("No primary-tumor samples found.")
    return expr, tumor_cols


def load_surfaceome() -> pd.DataFrame:
    surf = pd.read_excel(
        DATA_DIR / "table_S3_surfaceome.xlsx",
        sheet_name="in silico surfaceome only",
        engine="openpyxl",
        header=1,
    )
    surf = surf.rename(columns={"UniProt gene": "gene", "Ensembl gene": "ens"})
    surf = surf[surf["gene"].notna()].copy()
    surf["gene"] = surf["gene"].astype(str).str.strip()
    surf = surf[~surf["gene"].isin(["", "nan", "None"])]
    return surf


def load_morphology(tumor_cols: list[str]) -> pd.Series:
    cl = pd.read_csv(DATA_DIR / "TCGA-THCA.clinical.tsv.gz", sep="\t", low_memory=False)
    morph = cl.set_index(cl["sample"].str[:12])["morphology.diagnoses"]
    morph = morph[~morph.index.duplicated()]
    return morph.reindex([c[:12] for c in tumor_cols])


def load_purity(tumor_cols: list[str]) -> pd.Series:
    pur = pd.read_csv(DATA_DIR / "tcga_absolute_purity.txt", sep="\t")
    pur["sample15"] = pur["array"].str[:15]
    purity = pur.set_index("sample15")["purity"]
    purity = purity[~purity.index.duplicated()]
    return purity.reindex([c[:15] for c in tumor_cols])


def corr_block(x: np.ndarray, y: np.ndarray, label: str) -> dict:
    pr, pp = stats.pearsonr(x, y)
    sr, sp = stats.spearmanr(x, y)
    lo, hi = spearman_ci(x, y)
    return {
        "subset": label,
        "n": int(len(x)),
        "pearson_r": round(float(pr), 4),
        "pearson_p": float(pp),
        "spearman_rho": round(float(sr), 4),
        "spearman_p": float(sp),
        "spearman_rho_ci95_lo": round(lo, 4),
        "spearman_rho_ci95_hi": round(hi, 4),
    }


def rank_surface_expression(tum: pd.DataFrame, universe: list[str]) -> pd.DataFrame:
    """Rank surface genes by median log2(TPM+1) in primary tumors (desc)."""
    sub = tum.loc[universe]
    df = pd.DataFrame(
        {
            "gene": universe,
            "median_log2tpm1": sub.median(axis=1).to_numpy(),
            "mean_log2tpm1": sub.mean(axis=1).to_numpy(),
            "pct_log2tpm1_gt1": (sub > 1.0).mean(axis=1).to_numpy() * 100.0,
        }
    )
    df = df.sort_values("median_log2tpm1", ascending=False, kind="mergesort").reset_index(drop=True)
    df.insert(0, "expression_rank", np.arange(1, len(df) + 1))
    df["expression_percentile"] = 100.0 * (1 - (df["expression_rank"] - 1) / len(df))
    return df


def rank_coexpression(tum: pd.DataFrame, universe: list[str], anchor_vec: np.ndarray) -> pd.DataFrame:
    mat = tum.loc[universe].to_numpy(dtype=float)
    n = mat.shape[1]
    pear_r = corr_vs_vector(mat, anchor_vec)
    mat_rank = np.apply_along_axis(stats.rankdata, 1, mat)
    spear_r = corr_vs_vector(mat_rank, stats.rankdata(anchor_vec))
    df = pd.DataFrame(
        {
            "gene": universe,
            "spearman_r": spear_r,
            "spearman_p": pvalues(spear_r, n),
            "pearson_r": pear_r,
            "pearson_p": pvalues(pear_r, n),
            "mean_log2tpm1": mat.mean(axis=1),
        }
    )
    df["spearman_q"] = bh_fdr(df["spearman_p"].to_numpy())
    df["pearson_q"] = bh_fdr(df["pearson_p"].to_numpy())
    # Drop genes with undefined correlation (zero variance in the cohort).
    df = df[df["spearman_r"].notna()].copy()
    df = df.sort_values("spearman_r", ascending=False, kind="mergesort").reset_index(drop=True)
    df.insert(1, "spearman_rank", np.arange(1, len(df) + 1))
    df["pearson_rank"] = (
        df["pearson_r"].rank(ascending=False, method="min").astype("Int64")
    )
    return df


def focus_row(df: pd.DataFrame, gene: str, n: int) -> dict:
    hit = df[df["gene"] == gene]
    if hit.empty:
        return {"gene": gene, "in_universe": False}
    r = hit.iloc[0]
    s_rank = int(r["spearman_rank"])
    p_rank = int(r["pearson_rank"])
    return {
        "gene": gene,
        "in_universe": True,
        "is_top_spearman": s_rank == 1,
        "is_top_pearson": p_rank == 1,
        "spearman_rank": s_rank,
        "spearman_r": float(r["spearman_r"]),
        "spearman_q": float(r["spearman_q"]),
        "spearman_percentile": round(100 * (1 - (s_rank - 1) / n), 3),
        "pearson_rank": p_rank,
        "pearson_r": float(r["pearson_r"]),
        "pearson_q": float(r["pearson_q"]),
        "pearson_percentile": round(100 * (1 - (p_rank - 1) / n), 3),
    }


def write_writeup(
    n_samples: int,
    n_universe: int,
    n_surf_total: int,
    report: dict,
    coexp: pd.DataFrame,
    expr_rank: pd.DataFrame,
    corr_rows: list[dict],
    purity_rows: list[dict],
    med_anchor: float,
    med_focus: float,
) -> str:
    q = report
    top_gene = coexp.iloc[0]["gene"]
    if q.get("in_universe") and q["is_top_spearman"]:
        verdict = (
            f"YES — {FOCUS} is the top surface-gene co-expression partner of {ANCHOR}."
        )
    elif q.get("in_universe"):
        verdict = (
            f"NO — {FOCUS} is NOT the single top surface-gene co-expression "
            f"partner of {ANCHOR}."
        )
    else:
        verdict = f"{FOCUS} is not in the surface-gene universe."

    a_er = expr_rank[expr_rank["gene"] == ANCHOR]
    f_er = expr_rank[expr_rank["gene"] == FOCUS]
    a_txt = (
        f"#{int(a_er.iloc[0]['expression_rank'])} of {len(expr_rank)} "
        f"(median log2(TPM+1) = {a_er.iloc[0]['median_log2tpm1']:.2f})"
        if not a_er.empty
        else "not ranked"
    )
    f_txt = (
        f"#{int(f_er.iloc[0]['expression_rank'])} of {len(expr_rank)} "
        f"(median log2(TPM+1) = {f_er.iloc[0]['median_log2tpm1']:.2f})"
        if not f_er.empty
        else "not ranked"
    )

    main = corr_rows[0]
    lines = [
        f"# B1 analog — TCGA-THCA: {ANCHOR} × {FOCUS} surface co-expression rank",
        "",
        f"**Question:** Is `{FOCUS}` the top co-expression partner of `{ANCHOR}` "
        "among surface genes in TCGA-THCA?",
        "",
        f"**Answer:** {verdict}",
        "",
        f"- Cohort: TCGA-THCA primary tumors, n={n_samples} "
        "(sample-type 01, one vial per patient).",
        f"- Expression: UCSC Xena GDC hub `TCGA-THCA.star_tpm` = log2(TPM+1), "
        "GENCODE v36.",
        f"- Surface-gene universe: {n_universe} ranked genes "
        f"({n_surf_total} unique SURFY symbols in the workbook; "
        f"`{ANCHOR}` and zero-variance genes removed from the co-expression ranking).",
        "- Primary metric: Spearman correlation (honest full-universe rank; "
        "no pre-selected neighbourhood).",
        "",
        "## Co-expression rank (claim B1 analog)",
        "",
        "| metric | value |",
        "| --- | --- |",
        f"| {FOCUS} Spearman rank | **#{q['spearman_rank']} of {n_universe}** "
        f"({q['spearman_percentile']}th percentile) |",
        f"| {FOCUS} Spearman rho | {q['spearman_r']:.3f} "
        f"(FDR q={q['spearman_q']:.2e}) |",
        f"| {FOCUS} Pearson rank | #{q['pearson_rank']} of {n_universe} |",
        f"| {FOCUS} Pearson rho | {q['pearson_r']:.3f} |",
        f"| Actual #1 (Spearman) | {top_gene} "
        f"(rho={coexp.iloc[0]['spearman_r']:.3f}) |",
        f"| {ANCHOR}–{FOCUS} Spearman (bootstrap 95% CI) | "
        f"{main['spearman_rho']:.3f} "
        f"[{main['spearman_rho_ci95_lo']:.2f}, {main['spearman_rho_ci95_hi']:.2f}] |",
        "",
        "## Surface-expression rank (median log2(TPM+1))",
        "",
        "This is the original “surface rank” readout: how highly each gene is "
        "expressed among the surfaceome, independent of co-expression.",
        "",
        "| gene | expression rank |",
        "| --- | --- |",
        f"| {ANCHOR} | {a_txt} |",
        f"| {FOCUS} | {f_txt} |",
        "",
        f"Median log2(TPM+1) in the analysis cohort: {ANCHOR} = {med_anchor:.2f}, "
        f"{FOCUS} = {med_focus:.2f}. Both are highly expressed in THCA.",
        "",
        "## Sensitivity",
        "",
        "| analysis | n | Spearman rho | p |",
        "| --- | --- | --- | --- |",
    ]
    for r in corr_rows:
        lines.append(
            f"| {r['subset']} | {r['n']} | {r['spearman_rho']:.3f} "
            f"[{r['spearman_rho_ci95_lo']:.2f}, {r['spearman_rho_ci95_hi']:.2f}] | "
            f"{r['spearman_p']:.1e} |"
        )
    lines += [
        "",
        "Purity (ABSOLUTE, DNA-based):",
        "",
        "| quantity | n | rho | p |",
        "| --- | --- | --- | --- |",
    ]
    for r in purity_rows:
        lines.append(f"| {r['quantity']} | {r['n']} | {r['rho']:.3f} | {r['p']:.1e} |")

    lines += [
        "",
        f"## Top 15 surface-gene partners of {ANCHOR} (Spearman)",
        "",
        "| rank | gene | spearman_rho | pearson_rho | FDR q |",
        "| --- | --- | --- | --- | --- |",
    ]
    for _, r in coexp.head(15).iterrows():
        star = "  ← focus" if r["gene"] == FOCUS else ""
        lines.append(
            f"| {int(r['spearman_rank'])} | {r['gene']}{star} | "
            f"{r['spearman_r']:.3f} | {r['pearson_r']:.3f} | {r['spearman_q']:.2e} |"
        )

    lines += [
        "",
        "## Honest caveats",
        "",
        "1. **This is not a pan-cancer test of claim B1.** Claim B1 is stated "
        "for a TCGA pan-cancer surfaceome ranking. THCA is one analog cohort. "
        "A top (or near-top) rank here does not prove the pan-cancer claim; "
        "a mid-pack rank here does not disprove it outside thyroid.",
        "2. **Bulk co-expression is not proof of a protein complex or of "
        "co-occurrence on the same cell.** THCA is epithelium-rich, so two "
        "epithelial surface genes can correlate through shared epithelial "
        "content. ABSOLUTE purity adjustment is reported; DNA purity is an "
        "imperfect proxy for the epithelial mRNA fraction.",
        "3. **mRNA ≠ protein.** ADC-target relevance needs IHC / protein. "
        "This analysis cannot replace that.",
        "4. **Xena GDC `star_tpm` is treated as log2(TPM+1).** The matrix "
        "floor is 0.0 (not −9.97), which is the log2(TPM+1) signature. "
        "Ranks are invariant to any strictly increasing transform, so a "
        "mistaken invert-to-linear-TPM step would not change ranks, only "
        "the reported TPM numbers.",
        "5. **TCGA-THCA is almost all papillary thyroid carcinoma** "
        "(classic 8260/3 + follicular variant 8340/3). Follicular, poorly "
        "differentiated, and anaplastic carcinomas are essentially absent. "
        "Nothing here speaks to those histologies, or to ICI / ADC response. "
        "Pooled rho is higher than either histology subset "
        "(see Sensitivity): between-subtype mean shift inflates the "
        "headline number. Quote the within-histology figures if the claim "
        "is about a within-tumor-type relationship.",
        "6. Surfaceome membership is the 2018 SURFY in-silico set (2,886 "
        "proteins). Gene-symbol matching misses retired / renamed symbols; "
        "unmatched genes are dropped, not imputed.",
        "",
        "## Reproduce",
        "",
        "```bash",
        "pip install -r requirements.txt",
        "python3 scripts/w200/B1_THCA/download_data.py",
        "python3 scripts/w200/B1_THCA/run_analysis.py",
        "```",
        "",
        "Outputs: `coexpression_TACSTD2_surfaceome.csv`, `top200.csv`, "
        "`surface_expression_rank.csv`, `coexpression_main.csv`, "
        "`purity_adjusted.csv`, `summary.json`, `summary.md`, "
        "`scatter_TACSTD2_vs_CLDN4.png`, `top_partners_TACSTD2.png`.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    expr, tumor_cols = load_expression()
    tum = expr[tumor_cols]
    if ANCHOR not in tum.index or FOCUS not in tum.index:
        raise SystemExit(f"{ANCHOR} or {FOCUS} missing from the mapped matrix.")

    x = tum.loc[ANCHOR].to_numpy(dtype=float)
    y = tum.loc[FOCUS].to_numpy(dtype=float)
    n_samples = len(tumor_cols)

    surf = load_surfaceome()
    n_surf_proteins = int(len(surf))
    n_surf_total = int(surf["gene"].nunique())
    surf_genes = sorted(set(surf["gene"]))
    present = [g for g in surf_genes if g in tum.index]
    universe = [g for g in present if g != ANCHOR]
    n_universe = len(universe)

    # --- pairwise TACSTD2 × CLDN4 ------------------------------------------
    morph = load_morphology(tumor_cols)
    morph_vals = morph.to_numpy()
    corr_rows = [corr_block(x, y, "all_primary_tumors")]
    for code, label in (
        (CLASSIC_PAPILLARY, "classic_papillary_8260_3"),
        (FOLLICULAR_VARIANT, "follicular_variant_8340_3"),
    ):
        mask = morph_vals == code
        if mask.sum() >= 20:
            corr_rows.append(corr_block(x[mask], y[mask], label))

    purity = load_purity(tumor_cols).to_numpy(dtype=float)
    pmask = ~np.isnan(purity)
    r_part, p_part = partial_spearman(x[pmask], y[pmask], purity[pmask])
    sr_sub, sp_sub = stats.spearmanr(x[pmask], y[pmask])
    purity_rows = [
        {
            "quantity": "spearman_unadjusted_purity_subset",
            "n": int(pmask.sum()),
            "rho": round(float(sr_sub), 4),
            "p": float(sp_sub),
        },
        {
            "quantity": "spearman_partial_given_ABSOLUTE_purity",
            "n": int(pmask.sum()),
            "rho": round(r_part, 4),
            "p": p_part,
        },
        {
            "quantity": "spearman_TACSTD2_vs_purity",
            "n": int(pmask.sum()),
            "rho": round(float(stats.spearmanr(x[pmask], purity[pmask]).statistic), 4),
            "p": float(stats.spearmanr(x[pmask], purity[pmask]).pvalue),
        },
        {
            "quantity": "spearman_CLDN4_vs_purity",
            "n": int(pmask.sum()),
            "rho": round(float(stats.spearmanr(y[pmask], purity[pmask]).statistic), 4),
            "p": float(stats.spearmanr(y[pmask], purity[pmask]).pvalue),
        },
    ]
    pd.DataFrame(corr_rows).to_csv(OUT_DIR / "coexpression_main.csv", index=False)
    pd.DataFrame(purity_rows).to_csv(OUT_DIR / "purity_adjusted.csv", index=False)

    # --- surfaceome co-expression ranking (claim B1) -----------------------
    coexp = rank_coexpression(tum, universe, x)
    n_universe = len(coexp)
    coexp.to_csv(OUT_DIR / f"coexpression_{ANCHOR}_surfaceome.csv", index=False)
    coexp.head(WINDOW).to_csv(OUT_DIR / f"top{WINDOW}.csv", index=False)
    report = focus_row(coexp, FOCUS, n_universe)

    # --- surface-expression ranking ----------------------------------------
    expr_rank = rank_surface_expression(tum, present)
    expr_rank.to_csv(OUT_DIR / "surface_expression_rank.csv", index=False)

    med_anchor = float(np.median(x))
    med_focus = float(np.median(y))

    checksums = {
        p.name: md5(DATA_DIR / p)
        for p in (
            Path("TCGA-THCA.star_tpm.tsv.gz"),
            Path("gencode.v36.probemap"),
            Path("table_S3_surfaceome.xlsx"),
            Path("TCGA-THCA.clinical.tsv.gz"),
            Path("tcga_absolute_purity.txt"),
        )
        if (DATA_DIR / p).exists()
    }

    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "analysis": "B1-analog TACSTD2–CLDN4 surfaceome co-expression + expression rank",
        "cohort": "TCGA-THCA",
        "expression_dataset": "UCSC Xena GDC hub TCGA-THCA.star_tpm (log2(TPM+1), GENCODE v36)",
        "surfaceome_source": "Bausch-Fluck et al. 2018 in-silico surfaceome (table S3)",
        "anchor_gene": ANCHOR,
        "focus_gene": FOCUS,
        "sample_type_codes": ["01"],
        "n_samples": n_samples,
        "n_surfaceome_protein_rows": n_surf_proteins,
        "n_surfaceome_unique_symbols": n_surf_total,
        "surface_gene_universe": n_universe,
        "correlation_primary": "spearman",
        "window": WINDOW,
        "median_log2tpm1_TACSTD2": round(med_anchor, 3),
        "median_log2tpm1_CLDN4": round(med_focus, 3),
        "focus_result": report,
        "pairwise": corr_rows,
        "purity": purity_rows,
        "expression_rank": {
            g: (
                {
                    "rank": int(er.iloc[0]["expression_rank"]),
                    "out_of": int(len(expr_rank)),
                    "median_log2tpm1": float(er.iloc[0]["median_log2tpm1"]),
                }
                if not er.empty
                else None
            )
            for g, er in (
                (ANCHOR, expr_rank[expr_rank["gene"] == ANCHOR]),
                (FOCUS, expr_rank[expr_rank["gene"] == FOCUS]),
            )
        },
        "top10_spearman": coexp.head(10)[
            ["gene", "spearman_rank", "spearman_r", "spearman_q"]
        ].to_dict("records"),
        "input_checksums_md5": checksums,
    }
    (OUT_DIR / "summary.json").write_text(json.dumps(summary, indent=2))

    writeup = write_writeup(
        n_samples,
        n_universe,
        n_surf_total,
        report,
        coexp,
        expr_rank,
        corr_rows,
        purity_rows,
        med_anchor,
        med_focus,
    )
    (OUT_DIR / "summary.md").write_text(writeup)
    (OUT_DIR / "README.md").write_text(writeup)

    # --- figures -----------------------------------------------------------
    fig, ax = plt.subplots(figsize=(5.2, 5.0), dpi=150)
    classic = morph_vals == CLASSIC_PAPILLARY
    fv = morph_vals == FOLLICULAR_VARIANT
    other = ~(classic | fv)
    ax.scatter(
        x[classic], y[classic], s=14, alpha=0.55, color="#2b6a99",
        label=f"classic papillary 8260/3 (n={int(classic.sum())})",
    )
    ax.scatter(
        x[fv], y[fv], s=16, alpha=0.7, color="#d1495b", marker="^",
        label=f"follicular variant 8340/3 (n={int(fv.sum())})",
    )
    if other.any():
        ax.scatter(
            x[other], y[other], s=16, alpha=0.8, color="#7a7a7a", marker="s",
            label=f"other histology (n={int(other.sum())})",
        )
    r0 = corr_rows[0]
    ax.set_xlabel(f"{ANCHOR}  log2(TPM+1)")
    ax.set_ylabel(f"{FOCUS}  log2(TPM+1)")
    ax.set_title(
        f"TCGA-THCA primary tumors (n={n_samples})\n"
        f"Spearman rho = {r0['spearman_rho']:.2f} "
        f"[{r0['spearman_rho_ci95_lo']:.2f}, {r0['spearman_rho_ci95_hi']:.2f}]",
        fontsize=10,
    )
    ax.legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT_DIR / f"scatter_{ANCHOR}_vs_{FOCUS}.png")
    plt.close(fig)

    topn = coexp.head(20).iloc[::-1]
    colors = ["#d62728" if g == FOCUS else "#4c78a8" for g in topn["gene"]]
    fig, ax = plt.subplots(figsize=(6.2, 7.0), dpi=150)
    ax.barh(topn["gene"], topn["spearman_r"], color=colors)
    ax.set_xlabel(f"Spearman rho with {ANCHOR}")
    ax.set_title(
        f"Top surface-gene co-expression with {ANCHOR}\n"
        f"TCGA-THCA (focus: {FOCUS})"
    )
    fig.tight_layout()
    fig.savefig(OUT_DIR / f"top_partners_{ANCHOR}.png")
    plt.close(fig)

    print(writeup)
    print(f"\n[done] outputs written to {OUT_DIR}")


if __name__ == "__main__":
    main()
