#!/usr/bin/env python3
"""B1 analog in TCGA-UCS: TACSTD2 surfaceome rank of CLDN4.

The B1 claim is that CLDN4 is the *top* cell-surface gene co-expressed
with TACSTD2 (TROP2). This script scores every in-silico surfaceome gene
against TACSTD2 in TCGA uterine carcinosarcoma and reports where CLDN4
actually lands. It does not restrict to a hand-picked neighbourhood.

Primary matrix: Xena GDC STAR TPM, log2(TPM+1), GENCODE v36.
Surface universe: Bausch-Fluck et al. 2018 table S3 ("in silico
surfaceome only", 2,886 proteins).

Also reported, because n is small (n=57) and UCS is biphasic:
  * TACSTD2 × CLDN4 Spearman/Pearson with bootstrap 95% CI
  * bootstrap distribution of CLDN4's surface rank
  * rank-based partial correlation given ABSOLUTE purity
  * the same surface ranking on the legacy HiSeqV2 matrix (if present)

Outputs land in results/w200/B1_UCS/.
"""
from __future__ import annotations

import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = Path(os.environ.get("B1_UCS_DATA", ROOT / "data"))
OUT_DIR = ROOT / "results" / "w200" / "B1_UCS"
RNG = np.random.default_rng(20260816)
N_BOOT = 2000
N_BOOT_RANK = 1000
WINDOW = 200
ANCHOR = "TACSTD2"
FOCUS = "CLDN4"


def _strip_ensembl(x: str) -> str:
    x = str(x)
    return x.split(".")[0] if x.startswith("ENS") else x


def _corr_rows_vs_vec(mat: np.ndarray, vec: np.ndarray) -> np.ndarray:
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


def spearman_ci(x: np.ndarray, y: np.ndarray, n_boot: int = N_BOOT) -> tuple[float, float]:
    n = len(x)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = RNG.integers(0, n, n)
        boots[i] = stats.spearmanr(x[idx], y[idx]).statistic
    return float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))


def partial_spearman(x, y, z) -> tuple[float, float]:
    rx, ry, rz = (stats.rankdata(v) for v in (x, y, z))

    def resid(a, b):
        beta = np.polyfit(b, a, 1)
        return a - np.polyval(beta, b)

    r, p = stats.pearsonr(resid(rx, rz), resid(ry, rz))
    return float(r), float(p)


def primary_tumor_columns(columns) -> list[str]:
    # GDC barcodes are 16+ chars (vial letter); legacy HiSeqV2 uses 15-char
    # sample barcodes. Sample-type code is always characters 13-14.
    tumor = sorted(c for c in columns if len(str(c)) >= 15 and str(c)[13:15] == "01")
    patients: dict[str, str] = {}
    for c in tumor:
        patients.setdefault(str(c)[:12], c)
    return sorted(patients.values())


def load_star_tpm() -> tuple[pd.DataFrame, list[str]]:
    expr_path = DATA_DIR / "TCGA-UCS.star_tpm.tsv.gz"
    map_path = DATA_DIR / "gencode.v36.annotation.gtf.gene.probemap"
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    pmap = pd.read_csv(map_path, sep="\t")
    # Xena probemap columns are typically id, gene, chrom, chromStart, chromEnd, strand
    id_col = "id" if "id" in pmap.columns else pmap.columns[0]
    gene_col = "gene" if "gene" in pmap.columns else pmap.columns[1]
    id2sym = pmap.set_index(id_col)[gene_col]
    # keep Ensembl (version-stripped) as a column; set a symbol index after collapse
    expr = expr.copy()
    expr["ensembl"] = [_strip_ensembl(i) for i in expr.index]
    expr["symbol"] = id2sym.reindex(expr.index).values
    # some IDs may miss the versioned key; try stripped
    miss = expr["symbol"].isna()
    if miss.any():
        id2sym_stripped = pd.Series(
            id2sym.values, index=[_strip_ensembl(i) for i in id2sym.index]
        )
        id2sym_stripped = id2sym_stripped[~id2sym_stripped.index.duplicated()]
        expr.loc[miss, "symbol"] = id2sym_stripped.reindex(expr.loc[miss, "ensembl"]).values
    samples = primary_tumor_columns(expr.columns.drop(["ensembl", "symbol"]))
    return expr, samples


def load_surfaceome() -> pd.DataFrame:
    path = DATA_DIR / "table_S3_surfaceome.xlsx"
    surf = pd.read_excel(
        path, sheet_name="in silico surfaceome only", engine="openpyxl", header=1
    )
    surf = surf.rename(
        columns={
            "UniProt gene": "gene",
            "Ensembl gene": "ensembl",
            "Surfaceome Label": "surfaceome_label",
            "Surfaceome Label Source": "surfaceome_source",
        }
    )
    surf = surf[surf["gene"].notna()].copy()
    surf["gene"] = surf["gene"].astype(str)
    surf["ensembl"] = surf["ensembl"].apply(
        lambda x: _strip_ensembl(x) if pd.notna(x) else ""
    )
    return surf


def build_surface_matrix(
    expr: pd.DataFrame, surf: pd.DataFrame, samples: list[str]
) -> tuple[pd.DataFrame, dict]:
    """One row per surfaceome gene; expression taken by Ensembl then symbol."""
    expr_by_ens = expr.drop_duplicates("ensembl").set_index("ensembl")
    # symbol lookup: prefer the highest-mean protein-coding-like row
    tmp = expr.copy()
    tmp["_mean"] = tmp[samples].mean(axis=1)
    tmp = tmp[tmp["symbol"].notna()].sort_values("_mean", ascending=False)
    expr_by_sym = tmp.drop_duplicates("symbol").set_index("symbol")

    rows = []
    match_ens = match_sym = miss = 0
    for rec in surf.itertuples(index=False):
        gene = rec.gene
        ens = rec.ensembl
        src = None
        if ens and ens in expr_by_ens.index:
            vec = expr_by_ens.loc[ens, samples]
            src = "ensembl"
            match_ens += 1
        elif gene in expr_by_sym.index:
            vec = expr_by_sym.loc[gene, samples]
            src = "symbol"
            match_sym += 1
        else:
            miss += 1
            continue
        rows.append(
            {
                "gene": gene,
                "ensembl": ens,
                "match": src,
                **{s: float(vec[s]) for s in samples},
            }
        )
    mat = pd.DataFrame(rows)
    # collapse rare duplicate UniProt gene names by mean
    if mat["gene"].duplicated().any():
        num = mat.groupby("gene", sort=False)[samples].mean()
        meta = mat.drop_duplicates("gene").set_index("gene")[["ensembl", "match"]]
        mat = meta.join(num).reset_index()
    info = {
        "surfy_n": int(len(surf)),
        "matched_ensembl": match_ens,
        "matched_symbol_only": match_sym,
        "unmatched": miss,
        "universe_incl_anchor": int(len(mat)),
    }
    return mat.set_index("gene"), info


def rank_partners(mat: pd.DataFrame, anchor_vec: np.ndarray, samples: list[str]) -> pd.DataFrame:
    X = mat[samples].to_numpy(dtype=float)
    n = len(samples)
    pear_r = _corr_rows_vs_vec(X, anchor_vec)
    anchor_rank = stats.rankdata(anchor_vec)
    mat_rank = np.apply_along_axis(stats.rankdata, 1, X)
    spear_r = _corr_rows_vs_vec(mat_rank, anchor_rank)
    df = pd.DataFrame(
        {
            "gene": mat.index.to_numpy(),
            "spearman_r": spear_r,
            "spearman_p": _pvalues(spear_r, n),
            "pearson_r": pear_r,
            "pearson_p": _pvalues(pear_r, n),
            "mean_log2_expr": X.mean(axis=1),
            "pct_expressed": (X > 0).mean(axis=1) * 100.0,
            "match": mat["match"].to_numpy() if "match" in mat.columns else "na",
            "ensembl": mat["ensembl"].to_numpy() if "ensembl" in mat.columns else "",
        }
    )
    df["spearman_q"] = _bh_fdr(df["spearman_p"].to_numpy())
    df["pearson_q"] = _bh_fdr(df["pearson_p"].to_numpy())
    # Zero-variance genes have undefined correlation; they are not ranked.
    n_undefined = int(df["spearman_r"].isna().sum())
    df = df.dropna(subset=["spearman_r"]).copy()
    df.attrs["n_undefined"] = n_undefined
    df = df.sort_values("spearman_r", ascending=False, kind="mergesort").reset_index(drop=True)
    df.insert(1, "spearman_rank", np.arange(1, len(df) + 1))
    df["pearson_rank"] = (
        df["pearson_r"].rank(ascending=False, method="min").astype("Int64")
    )
    return df


def bootstrap_focus_rank(
    mat: pd.DataFrame, anchor_vec: np.ndarray, samples: list[str], focus: str
) -> np.ndarray:
    X = mat[samples].to_numpy(dtype=float)
    focus_i = list(mat.index).index(focus)
    n = len(samples)
    ranks = np.empty(N_BOOT_RANK, dtype=int)
    for b in range(N_BOOT_RANK):
        idx = RNG.integers(0, n, n)
        xb = X[:, idx]
        ab = anchor_vec[idx]
        ar = stats.rankdata(ab)
        mr = np.apply_along_axis(stats.rankdata, 1, xb)
        rho = _corr_rows_vs_vec(mr, ar)
        # rank of focus (1 = highest rho); ties → min rank
        ranks[b] = int((rho > rho[focus_i]).sum() + 1)
    return ranks


def load_purity(samples: list[str]) -> pd.Series:
    path = DATA_DIR / "tcga_absolute_purity.txt"
    pur = pd.read_csv(path, sep="\t")
    key = "array" if "array" in pur.columns else pur.columns[0]
    pur["sample15"] = pur[key].astype(str).str[:15]
    s = pur.set_index("sample15")["purity"]
    s = s[~s.index.duplicated()]
    return s.reindex([c[:15] for c in samples])


def load_hiseqv2() -> tuple[pd.DataFrame, list[str]] | None:
    path = DATA_DIR / "TCGA-UCS.HiSeqV2.gz"
    if not path.exists():
        return None
    expr = pd.read_csv(path, sep="\t", index_col=0)
    expr.index = expr.index.astype(str)
    if expr.index.duplicated().any():
        expr = expr.groupby(level=0).mean()
    samples = primary_tumor_columns(expr.columns)
    return expr, samples


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("[load] STAR TPM + GENCODE v36 probemap")
    expr, samples = load_star_tpm()
    print(f"       genes x all-cols = {expr.shape}; primary tumors (1/patient) = {len(samples)}")

    print("[load] SURFY table S3")
    surf = load_surfaceome()
    surf_mat, match_info = build_surface_matrix(expr, surf, samples)
    print(f"       match info: {match_info}")

    if ANCHOR not in surf_mat.index:
        # still need the anchor vector from the full matrix
        pass
    # anchor expression from Ensembl/symbol of the full matrix, not the
    # surface table (anchor is removed from the ranking universe)
    expr_by_sym = (
        expr[expr["symbol"].notna()]
        .assign(_mean=lambda d: d[samples].mean(axis=1))
        .sort_values("_mean", ascending=False)
        .drop_duplicates("symbol")
        .set_index("symbol")
    )
    if ANCHOR not in expr_by_sym.index:
        raise SystemExit(f"anchor {ANCHOR} not found after symbol mapping")
    if FOCUS not in expr_by_sym.index:
        raise SystemExit(f"focus {FOCUS} not found after symbol mapping")
    if FOCUS not in surf_mat.index:
        raise SystemExit(f"focus {FOCUS} not in matched surfaceome")

    anchor_vec = expr_by_sym.loc[ANCHOR, samples].to_numpy(dtype=float)
    focus_vec = expr_by_sym.loc[FOCUS, samples].to_numpy(dtype=float)

    universe = surf_mat.drop(index=[ANCHOR], errors="ignore")
    sd = universe[samples].std(axis=1)
    n_zero_var = int((sd == 0).sum())
    universe = universe.loc[sd > 0]
    n_universe = len(universe)
    print(
        f"[data] surface universe (anchor removed, zero-variance dropped "
        f"n={n_zero_var}): {n_universe}"
    )

    df = rank_partners(universe, anchor_vec, samples)
    df.to_csv(OUT_DIR / f"coexpression_{ANCHOR}_surfaceome.csv", index=False)
    df.head(WINDOW).to_csv(OUT_DIR / f"top{WINDOW}.csv", index=False)

    row = df[df["gene"] == FOCUS].iloc[0]
    s_rank = int(row["spearman_rank"])
    p_rank = int(row["pearson_rank"])

    # direct TACSTD2 × CLDN4
    pr, pp = stats.pearsonr(anchor_vec, focus_vec)
    sr, sp = stats.spearmanr(anchor_vec, focus_vec)
    lo, hi = spearman_ci(anchor_vec, focus_vec)
    direct = {
        "subset": "all_primary_tumors",
        "n": int(len(samples)),
        "pearson_r": float(pr),
        "pearson_p": float(pp),
        "spearman_rho": float(sr),
        "spearman_p": float(sp),
        "spearman_rho_ci95_lo": lo,
        "spearman_rho_ci95_hi": hi,
        "median_log2tpm1_TACSTD2": float(np.median(anchor_vec)),
        "median_log2tpm1_CLDN4": float(np.median(focus_vec)),
    }
    pd.DataFrame([direct]).to_csv(OUT_DIR / "coexpression_main.csv", index=False)

    # purity
    purity = load_purity(samples)
    pvec = purity.to_numpy(dtype=float)
    mask = ~np.isnan(pvec)
    r_part, p_part = partial_spearman(anchor_vec[mask], focus_vec[mask], pvec[mask])
    sr_sub, sp_sub = stats.spearmanr(anchor_vec[mask], focus_vec[mask])
    rho_t_p = stats.spearmanr(anchor_vec[mask], pvec[mask])
    rho_c_p = stats.spearmanr(focus_vec[mask], pvec[mask])
    purity_rows = [
        {
            "quantity": "spearman_unadjusted_purity_subset",
            "n": int(mask.sum()),
            "rho": float(sr_sub),
            "p": float(sp_sub),
        },
        {
            "quantity": "spearman_partial_given_ABSOLUTE_purity",
            "n": int(mask.sum()),
            "rho": r_part,
            "p": p_part,
        },
        {
            "quantity": "spearman_TACSTD2_vs_purity",
            "n": int(mask.sum()),
            "rho": float(rho_t_p.statistic),
            "p": float(rho_t_p.pvalue),
        },
        {
            "quantity": "spearman_CLDN4_vs_purity",
            "n": int(mask.sum()),
            "rho": float(rho_c_p.statistic),
            "p": float(rho_c_p.pvalue),
        },
    ]
    pd.DataFrame(purity_rows).to_csv(OUT_DIR / "purity_adjusted.csv", index=False)

    # bootstrap rank of CLDN4 (small-n honesty)
    print(f"[boot] {N_BOOT_RANK} resamples of CLDN4 surface rank")
    boot_ranks = bootstrap_focus_rank(universe, anchor_vec, samples, FOCUS)
    boot_rank_summary = {
        "n_boot": N_BOOT_RANK,
        "median_rank": int(np.median(boot_ranks)),
        "mean_rank": float(np.mean(boot_ranks)),
        "ci95_lo": int(np.percentile(boot_ranks, 2.5)),
        "ci95_hi": int(np.percentile(boot_ranks, 97.5)),
        "frac_rank1": float((boot_ranks == 1).mean()),
        "frac_top5": float((boot_ranks <= 5).mean()),
        "frac_top10": float((boot_ranks <= 10).mean()),
        "frac_top20": float((boot_ranks <= 20).mean()),
    }
    pd.Series(boot_ranks, name="cldn4_spearman_rank").to_csv(
        OUT_DIR / "cldn4_rank_bootstrap.csv", index=False
    )

    # transcriptome-wide rank among expressed genes (context, not the B1 claim)
    num = expr[samples].to_numpy(dtype=float)
    symbols = expr["symbol"].to_numpy()
    keep = (num > 1.0).mean(axis=1) >= 0.20
    keep = keep & pd.notna(pd.Series(symbols)).to_numpy()
    tw = pd.DataFrame(num[keep], index=pd.Index(symbols[keep], name="gene"), columns=samples)
    tw = tw.assign(_mean=tw.mean(axis=1)).sort_values("_mean", ascending=False)
    tw = tw[~tw.index.duplicated()].drop(columns="_mean")
    tw = tw[tw.index != ANCHOR]
    tw_rho = pd.Series(
        _corr_rows_vs_vec(
            np.apply_along_axis(stats.rankdata, 1, tw.to_numpy(dtype=float)),
            stats.rankdata(anchor_vec),
        ),
        index=tw.index,
    ).sort_values(ascending=False)
    tw_rank = int(tw_rho.index.get_loc(FOCUS)) + 1 if FOCUS in tw_rho.index else None
    tw_top = tw_rho.head(25).rename("spearman_rho_vs_TACSTD2").reset_index()
    tw_top.insert(0, "rank", range(1, len(tw_top) + 1))
    tw_top.to_csv(OUT_DIR / "top25_tacstd2_partners_transcriptome.csv", index=False)

    # HiSeqV2 robustness
    hiseq_result = None
    hv = load_hiseqv2()
    if hv is not None:
        hexpr, hsamples = hv
        print(f"[hiseq] {hexpr.shape[0]} genes x {len(hsamples)} primary tumors")
        surf_genes = sorted({g for g in surf["gene"] if g not in ("nan", "None", "")})
        huniv = [g for g in surf_genes if g in hexpr.index and g != ANCHOR]
        if hsamples and ANCHOR in hexpr.index and FOCUS in hexpr.index:
            h_anchor = hexpr.loc[ANCHOR, hsamples].to_numpy(dtype=float)
            hmat = hexpr.loc[huniv, hsamples]
            hmat = hmat.assign(match="symbol", ensembl="")
            hdf = rank_partners(hmat, h_anchor, hsamples)
            hrow = hdf[hdf["gene"] == FOCUS]
            if not hrow.empty:
                hiseq_result = {
                    "n_samples": len(hsamples),
                    "surface_gene_universe": len(huniv),
                    "spearman_rank": int(hrow.iloc[0]["spearman_rank"]),
                    "spearman_r": float(hrow.iloc[0]["spearman_r"]),
                    "spearman_q": float(hrow.iloc[0]["spearman_q"]),
                    "pearson_rank": int(hrow.iloc[0]["pearson_rank"]),
                    "pearson_r": float(hrow.iloc[0]["pearson_r"]),
                    "top1": hdf.iloc[0]["gene"],
                    "top1_rho": float(hdf.iloc[0]["spearman_r"]),
                }
                hdf.head(20).to_csv(OUT_DIR / "hiseqv2_top20.csv", index=False)

    n_universe = len(df)
    report = {
        "gene": FOCUS,
        "in_universe": True,
        "is_top_spearman": s_rank == 1,
        "is_top_pearson": p_rank == 1,
        "spearman_rank": s_rank,
        "spearman_r": float(row["spearman_r"]),
        "spearman_p": float(row["spearman_p"]),
        "spearman_q": float(row["spearman_q"]),
        "spearman_percentile": round(100 * (1 - (s_rank - 1) / n_universe), 3),
        "pearson_rank": p_rank,
        "pearson_r": float(row["pearson_r"]),
        "pearson_q": float(row["pearson_q"]),
        "pearson_percentile": round(100 * (1 - (p_rank - 1) / n_universe), 3),
    }

    # clinical histology (all mixed Müllerian / carcinosarcoma)
    clin = pd.read_csv(DATA_DIR / "TCGA-UCS.clinical.tsv.gz", sep="\t")
    hist = clin["primary_diagnosis.diagnoses"].value_counts().to_dict()

    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "cohort": "TCGA-UCS",
        "expression_dataset": (
            "UCSC Xena GDC hub TCGA-UCS.star_tpm.tsv.gz (log2(TPM+1), GENCODE v36)"
        ),
        "surfaceome_source": "Bausch-Fluck et al. 2018 in-silico surfaceome (table S3)",
        "anchor_gene": ANCHOR,
        "focus_gene": FOCUS,
        "sample_type_codes": ["01"],
        "n_samples": len(samples),
        "n_patients": len({c[:12] for c in samples}),
        "histology": hist,
        "surface_gene_universe": n_universe,
        "surfaceome_match": {**match_info, "zero_variance_dropped": n_zero_var},
        "correlation_primary": "spearman",
        "window": WINDOW,
        "focus_result": report,
        "direct_coexpression": direct,
        "purity": purity_rows,
        "cldn4_rank_bootstrap": boot_rank_summary,
        "transcriptome_wide": {
            "n_expressed_genes": int(len(tw_rho)),
            "CLDN4_rank": tw_rank,
            "CLDN4_percentile": (
                round(100 * (1 - tw_rank / len(tw_rho)), 2) if tw_rank else None
            ),
        },
        "hiseqv2_robustness": hiseq_result,
        "top10_spearman": df.head(10)[
            ["gene", "spearman_rank", "spearman_r", "spearman_q"]
        ].to_dict("records"),
    }
    (OUT_DIR / "summary.json").write_text(json.dumps(summary, indent=2))

    verdict = (
        f"YES -- {FOCUS} is the top surface-gene co-expression partner of {ANCHOR}."
        if report["is_top_spearman"]
        else (
            f"NO -- {FOCUS} is NOT the single top surface-gene co-expression "
            f"partner of {ANCHOR}."
        )
    )
    top_gene = df.iloc[0]["gene"]
    lines = [
        f"# {FOCUS} vs {ANCHOR} co-expression in TCGA-UCS (honest ranking)",
        "",
        f"**Question:** Is `{FOCUS}` the top co-expression partner of "
        f"`{ANCHOR}` among surface genes in TCGA-UCS?",
        "",
        f"**Answer:** {verdict}",
        "",
        f"- Cohort: TCGA-UCS, {len(samples)} primary-tumour samples "
        f"(sample-type 01; one per patient). All cases are mixed Müllerian / "
        f"carcinosarcoma.",
        f"- Surface-gene universe: {n_universe} surfaceome genes present in the "
        f"STAR TPM matrix (anchor removed).",
        f"- Primary metric: Spearman correlation on log2(TPM+1).",
        f"- Small-n note: n=57. CLDN4 surface-rank bootstrap median "
        f"{boot_rank_summary['median_rank']} "
        f"(95% {boot_rank_summary['ci95_lo']}–{boot_rank_summary['ci95_hi']}); "
        f"fraction of resamples with rank 1 = {boot_rank_summary['frac_rank1']:.2f}.",
        "",
        "| metric | value |",
        "| --- | --- |",
        f"| {FOCUS} Spearman rank | **#{s_rank} of {n_universe}** "
        f"({report['spearman_percentile']}th percentile) |",
        f"| {FOCUS} Spearman rho | {report['spearman_r']:.3f} "
        f"(p={report['spearman_p']:.2e}, FDR q={report['spearman_q']:.2e}) |",
        f"| {FOCUS} Spearman rho 95% CI | {lo:.3f} – {hi:.3f} |",
        f"| {FOCUS} Pearson rank | #{p_rank} of {n_universe} |",
        f"| {FOCUS} Pearson rho | {report['pearson_r']:.3f} |",
        f"| Actual #1 (Spearman) | {top_gene} (rho={df.iloc[0]['spearman_r']:.3f}) |",
        f"| Partial rho given ABSOLUTE purity | {r_part:.3f} (n={int(mask.sum())}) |",
        "",
        f"## Top 15 surface-gene partners of {ANCHOR} (Spearman)",
        "",
        "| rank | gene | spearman_rho | pearson_rho | FDR q |",
        "| --- | --- | --- | --- | --- |",
    ]
    for _, r in df.head(15).iterrows():
        star = "  <-- focus" if r["gene"] == FOCUS else ""
        lines.append(
            f"| {int(r['spearman_rank'])} | {r['gene']}{star} | "
            f"{r['spearman_r']:.3f} | {r['pearson_r']:.3f} | {r['spearman_q']:.2e} |"
        )
    (OUT_DIR / "summary.md").write_text("\n".join(lines) + "\n")

    # plots
    fig, ax = plt.subplots(figsize=(5.2, 5.0), dpi=150)
    ax.scatter(anchor_vec, focus_vec, s=22, alpha=0.75, color="#2b6a99", edgecolors="none")
    ax.set_xlabel("TACSTD2  log2(TPM+1)")
    ax.set_ylabel("CLDN4  log2(TPM+1)")
    ax.set_title(
        f"TCGA-UCS primary tumors (n={len(samples)})\n"
        f"Spearman rho = {sr:.2f} [{lo:.2f}, {hi:.2f}], "
        f"surface rank #{s_rank}",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(OUT_DIR / "scatter_TACSTD2_vs_CLDN4.png")
    plt.close(fig)

    topn = df.head(20).iloc[::-1]
    colors = ["#d62728" if g == FOCUS else "#4c78a8" for g in topn["gene"]]
    fig, ax = plt.subplots(figsize=(6.2, 7.0), dpi=150)
    ax.barh(topn["gene"], topn["spearman_r"], color=colors)
    ax.set_xlabel("Spearman rho with TACSTD2")
    ax.set_title(
        f"Top surface-gene co-expression with TACSTD2\n"
        f"TCGA-UCS (focus: CLDN4, rank #{s_rank})"
    )
    fig.tight_layout()
    fig.savefig(OUT_DIR / "top_partners_TACSTD2.png")
    plt.close(fig)

    print("\n".join(lines))
    print(json.dumps({k: summary[k] for k in (
        "n_samples", "surface_gene_universe", "focus_result",
        "cldn4_rank_bootstrap", "hiseqv2_robustness", "transcriptome_wide",
    )}, indent=2))
    print(f"\n[done] outputs written to {OUT_DIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
