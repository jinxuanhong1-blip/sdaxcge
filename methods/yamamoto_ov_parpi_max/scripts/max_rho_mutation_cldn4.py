#!/usr/bin/env python3
"""Maximize |Spearman rho| for TCGA-OV mutation burden vs CLDN4.

This is the public cohort behind Yamamoto et al., Mol Cancer Ther 2022,
Fig. 4F (lower claudin-4, higher mutational burden). Raw matrices are not
committed. Set OV_MAX_DATA (default /tmp/ov_max).

Spearman rho is invariant to monotonic transforms, so log counts and log
expression are not separate tests. The grid below is the eligible search:

Predictors
  mut_all      unique MAF rows (any Variant_Classification)
  mut_nonsyn   nonsynonymous classes listed in NONSYNONYMOUS
  tmb          PanCan TMB_NONSYNONYMOUS (mutations per Mb)
  mut_silent   Silent only; negative control, never eligible to win

Expression
  rsem         Firehose/PanCan RNA-seq V2 RSEM (Yamamoto mRNA)
  star         GDC STAR TPM, Xena TCGA-OV.star_tpm (already log2-scale)
  protein      CPTAC CLDN4 in data_protein_quantification (not RPPA;
               the RPPA matrix has CLDN7 and no CLDN4 antibody)

Sample filters, each pre-specified
  all_paired, purity>=0.7, purity>=0.8, msi_sensor<3.5,
  count<=500 on the predictor in use, somatic BRCA1/2 nonsilent present,
  somatic BRCA1/2 nonsilent absent

Association
  spearman
  partial Spearman residualizing ranks on ABSOLUTE purity,
  purity+ploidy, log2(KRT18) and log2(KRT19), or purity+those keratins

A cell can win only if n>=40, the predictor is not silent, and both
variables have non-zero variance. Ties break toward the Yamamoto anchor
(mut_all x RSEM, all paired samples, unadjusted Spearman), then larger n,
then unadjusted, then RSEM, then mut_all.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(os.environ.get("OV_MAX_DATA", "/tmp/ov_max"))
OUT_T = ROOT / "results" / "tables"
OUT_F = ROOT / "results" / "figures"
OUT_T.mkdir(parents=True, exist_ok=True)
OUT_F.mkdir(parents=True, exist_ok=True)

NONSYNONYMOUS = {
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Nonstop_Mutation",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "In_Frame_Del",
    "In_Frame_Ins",
    "Splice_Site",
    "Translation_Start_Site",
}
STAR_GENES = {
    "ENSG00000189143": "CLDN4",
    "ENSG00000170421": "KRT8",
    "ENSG00000111057": "KRT18",
    "ENSG00000171345": "KRT19",
}
MIN_N = 40


def sid15(series: pd.Series) -> pd.Series:
    return series.astype(str).str.slice(0, 15)


def load_rsem() -> pd.DataFrame:
    path = DATA / "datahub" / "data_mrna_seq_v2_rsem.txt"
    df = pd.read_csv(path, sep="\t", usecols=lambda c: True)
    # columns: Hugo_Symbol, Entrez_Gene_Id, samples...
    gene_col = df.columns[0]
    entrez_col = df.columns[1]
    df = df.set_index(gene_col).drop(columns=[entrez_col])
    symbols = ["CLDN4", "KRT8", "KRT18", "KRT19"]
    missing = [g for g in symbols if g not in df.index]
    if missing:
        raise SystemExit(f"RSEM missing {missing}")
    out = df.loc[symbols].apply(pd.to_numeric, errors="coerce")
    out.columns = sid15(pd.Series(out.columns))
    out = out.T.groupby(level=0).mean().T
    return out


def load_star() -> pd.DataFrame:
    path = DATA / "star" / "TCGA-OV.star_tpm.tsv.gz"
    wanted = dict(STAR_GENES)
    found: dict[str, list[str]] = {}
    with __import__("gzip").open(path, "rt") as fh:
        samples = fh.readline().rstrip("\n").split("\t")[1:]
        for line in fh:
            ensg = line.split("\t", 1)[0].split(".")[0]
            if ensg in wanted and ensg not in found:
                found[ensg] = line.rstrip("\n").split("\t")[1:]
            if len(found) == len(wanted):
                break
    missing = [STAR_GENES[e] for e in wanted if e not in found]
    if missing:
        raise SystemExit(f"STAR missing {missing}")
    frame = {}
    for ensg, values in found.items():
        frame[STAR_GENES[ensg]] = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy()
    out = pd.DataFrame(frame, index=sid15(pd.Series(samples)))
    out = out.groupby(level=0).mean().T
    return out


def load_protein() -> pd.Series:
    path = DATA / "datahub" / "data_protein_quantification.txt"
    df = pd.read_csv(path, sep="\t")
    gene_col = df.columns[0]
    hit = df[df[gene_col].astype(str).str.startswith("CLDN4|")]
    if hit.empty:
        raise SystemExit("CPTAC protein matrix has no CLDN4 row")
    row = hit.iloc[0].drop(labels=[gene_col])
    ser = pd.to_numeric(row, errors="coerce")
    ser.index = sid15(pd.Series(ser.index))
    ser = ser.groupby(level=0).mean()
    ser.name = "CLDN4_protein"
    return ser


def load_purity() -> pd.DataFrame:
    path = DATA / "datahub" / "tcga_absolute_purity.txt"
    df = pd.read_csv(path, sep="\t")
    df["sample_id"] = sid15(df["sample"])
    df = df[df["call status"].astype(str).eq("called")].copy()
    df["purity"] = pd.to_numeric(df["purity"], errors="coerce")
    df["ploidy"] = pd.to_numeric(df["ploidy"], errors="coerce")
    df["prefer"] = df["solution"].astype(str).eq("new").astype(int)
    df = df.sort_values(["sample_id", "prefer"], ascending=[True, False])
    df = df.drop_duplicates("sample_id", keep="first")
    return df.set_index("sample_id")[["purity", "ploidy"]]


def load_clinical() -> pd.DataFrame:
    path = DATA / "datahub_clinical" / "data_clinical_sample.txt"
    if not path.exists():
        path = Path("/tmp/ov_clin_sample.txt")
    sample = pd.read_csv(path, sep="\t", comment="#")
    sample["sample_id"] = sid15(sample["SAMPLE_ID"])
    keep = sample.set_index("sample_id")[["TMB_NONSYNONYMOUS", "MSI_SENSOR_SCORE", "ONCOTREE_CODE", "SAMPLE_TYPE"]].copy()
    keep["TMB_NONSYNONYMOUS"] = pd.to_numeric(keep["TMB_NONSYNONYMOUS"], errors="coerce")
    keep["MSI_SENSOR_SCORE"] = pd.to_numeric(keep["MSI_SENSOR_SCORE"], errors="coerce")
    return keep[~keep.index.duplicated(keep="first")]


def load_mutation_counts() -> pd.DataFrame:
    path = DATA / "datahub" / "data_mutations.txt"
    usecols = [
        "Hugo_Symbol",
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Tumor_Sample_Barcode",
        "Variant_Classification",
        "Tumor_Seq_Allele2",
    ]
    df = pd.read_csv(path, sep="\t", usecols=usecols, low_memory=False)
    df["sample_id"] = sid15(df["Tumor_Sample_Barcode"])
    df = df[df["sample_id"].str.match(r"TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-01$")]
    key = ["sample_id", "Chromosome", "Start_Position", "End_Position", "Tumor_Seq_Allele2", "Variant_Classification", "Hugo_Symbol"]
    df = df.drop_duplicates(key)
    brca = df[df["Hugo_Symbol"].isin(["BRCA1", "BRCA2"]) & df["Variant_Classification"].isin(NONSYNONYMOUS)]
    brca_samples = set(brca["sample_id"])
    rows = []
    for sample_id, sub in df.groupby("sample_id"):
        vc = sub["Variant_Classification"]
        rows.append(
            {
                "sample_id": sample_id,
                "mut_all": int(len(sub)),
                "mut_nonsyn": int(vc.isin(NONSYNONYMOUS).sum()),
                "mut_silent": int(vc.eq("Silent").sum()),
                "brca_nonsyn": int(sample_id in brca_samples),
            }
        )
    return pd.DataFrame(rows).set_index("sample_id")


def spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    res = stats.spearmanr(x, y)
    return float(res.statistic), float(res.pvalue)


def partial_spearman(y: np.ndarray, x: np.ndarray, z: np.ndarray) -> tuple[float, float]:
    """Rank-residual partial Spearman. z is n x k. p uses df = n-2-k."""
    rx = stats.rankdata(x).astype(float)
    ry = stats.rankdata(y).astype(float)
    if z.ndim == 1:
        z = z[:, None]
    rz = np.column_stack([stats.rankdata(z[:, j]).astype(float) for j in range(z.shape[1])])
    design = np.column_stack([np.ones(len(rx)), rz])

    def _resid(v: np.ndarray) -> np.ndarray:
        beta, *_ = np.linalg.lstsq(design, v, rcond=None)
        return v - design @ beta

    xr, yr = _resid(rx), _resid(ry)
    if np.std(xr) == 0 or np.std(yr) == 0:
        return np.nan, np.nan
    r = float(np.corrcoef(xr, yr)[0, 1])
    k = rz.shape[1]
    df = len(x) - 2 - k
    if df <= 0 or not np.isfinite(r) or abs(r) >= 1:
        return r, np.nan
    tstat = r * np.sqrt(df / (1.0 - r * r))
    p = float(2 * stats.t.sf(abs(tstat), df))
    return r, p


def fisher_ci(r: float, n: int) -> tuple[float, float]:
    if n <= 3 or not np.isfinite(r) or abs(r) >= 1:
        return np.nan, np.nan
    z = np.arctanh(np.clip(r, -0.999999, 0.999999))
    se = 1.0 / np.sqrt(n - 3)
    lo, hi = np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)
    return float(lo), float(hi)


def bh(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    q = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return q.tolist()
    idx = np.where(ok)[0]
    order = idx[np.argsort(p[idx])]
    m = len(order)
    ranked = p[order]
    qq = ranked * m / (np.arange(1, m + 1))
    qq = np.minimum.accumulate(qq[::-1])[::-1]
    qq = np.clip(qq, 0, 1)
    q[order] = qq
    return q.tolist()


def build_table() -> pd.DataFrame:
    print("[load] mutations, expression, purity, clinical")
    mut = load_mutation_counts()
    rsem = load_rsem()
    star = load_star()
    protein = load_protein()
    purity = load_purity()
    clinical = load_clinical()

    expr = pd.DataFrame(index=sorted(set(rsem.columns) | set(star.columns) | set(protein.index)))
    expr["rsem"] = rsem.loc["CLDN4"].reindex(expr.index)
    expr["star"] = star.loc["CLDN4"].reindex(expr.index)
    expr["protein"] = protein.reindex(expr.index)
    expr["krt8_rsem"] = rsem.loc["KRT8"].reindex(expr.index)
    expr["krt18_rsem"] = rsem.loc["KRT18"].reindex(expr.index)
    expr["krt19_rsem"] = rsem.loc["KRT19"].reindex(expr.index)
    expr["krt8_star"] = star.loc["KRT8"].reindex(expr.index)
    expr["krt18_star"] = star.loc["KRT18"].reindex(expr.index)
    expr["krt19_star"] = star.loc["KRT19"].reindex(expr.index)

    tab = mut.join(expr, how="left").join(purity, how="left").join(clinical, how="left")
    tab.index.name = "sample_id"
    tab.to_csv(OUT_T / "sample_level.tsv", sep="\t", float_format="%.6g")
    print(f"[load] samples with a mutation count: {len(tab)}")
    return tab


def apply_filter(tab: pd.DataFrame, predictor: str, filt: str) -> pd.DataFrame:
    xcol = "TMB_NONSYNONYMOUS" if predictor == "tmb" else predictor
    mask = tab[xcol].notna()
    if filt == "purity_ge_0.7":
        mask &= tab["purity"].ge(0.7)
    elif filt == "purity_ge_0.8":
        mask &= tab["purity"].ge(0.8)
    elif filt == "msi_lt_3.5":
        mask &= tab["MSI_SENSOR_SCORE"].notna() & tab["MSI_SENSOR_SCORE"].lt(3.5)
    elif filt == "count_le_500":
        mask &= tab[xcol].le(500)
    elif filt == "tmb_le_10":
        mask &= tab[xcol].le(10)
    elif filt == "brca_nonsyn_present":
        mask &= tab["brca_nonsyn"].eq(1)
    elif filt == "brca_nonsyn_absent":
        mask &= tab["brca_nonsyn"].eq(0)
    elif filt != "all_paired":
        raise SystemExit(f"unknown filter {filt}")
    return tab.loc[mask]


def evaluate(tab: pd.DataFrame) -> pd.DataFrame:
    predictors = ["mut_all", "mut_nonsyn", "tmb", "mut_silent"]
    expressions = ["rsem", "star", "protein"]
    filters = [
        "all_paired",
        "purity_ge_0.7",
        "purity_ge_0.8",
        "msi_lt_3.5",
        "count_le_500",
        "brca_nonsyn_present",
        "brca_nonsyn_absent",
    ]
    adjustments = ["none", "purity", "purity_ploidy", "krt18_krt19", "purity_krt18_krt19"]
    rows = []
    for pred in predictors:
        xcol = "TMB_NONSYNONYMOUS" if pred == "tmb" else pred
        for exp in expressions:
            ycol = exp
            for filt in filters:
                mask = tab[xcol].notna() & tab[ycol].notna()
                if filt == "purity_ge_0.7":
                    mask &= tab["purity"].ge(0.7)
                elif filt == "purity_ge_0.8":
                    mask &= tab["purity"].ge(0.8)
                elif filt == "msi_lt_3.5":
                    mask &= tab["MSI_SENSOR_SCORE"].notna() & tab["MSI_SENSOR_SCORE"].lt(3.5)
                elif filt == "count_le_500":
                    if pred == "tmb":
                        mask &= tab[xcol].le(10)  # ~500 muts / ~50 Mb exome, labeled below
                    else:
                        mask &= tab[xcol].le(500)
                elif filt == "brca_nonsyn_present":
                    mask &= tab["brca_nonsyn"].eq(1)
                elif filt == "brca_nonsyn_absent":
                    mask &= tab["brca_nonsyn"].eq(0)
                sub = tab.loc[mask]
                for adj in adjustments:
                    cov_cols: list[str] = []
                    if "purity" in adj:
                        cov_cols.append("purity")
                    if "ploidy" in adj:
                        cov_cols.append("ploidy")
                    if "krt18_krt19" in adj:
                        if exp == "star":
                            cov_cols += ["krt18_star", "krt19_star"]
                        elif exp == "protein":
                            cov_cols += ["krt18_rsem", "krt19_rsem"]
                        else:
                            cov_cols += ["krt18_rsem", "krt19_rsem"]
                    use = sub
                    if cov_cols:
                        use = sub.dropna(subset=cov_cols)
                    n = int(len(use))
                    rho, p = np.nan, np.nan
                    if n >= 5 and use[xcol].nunique() > 1 and use[ycol].nunique() > 1:
                        x = use[xcol].to_numpy(dtype=float)
                        y = use[ycol].to_numpy(dtype=float)
                        if adj == "none":
                            rho, p = spearman(x, y)
                        else:
                            z = use[cov_cols].to_numpy(dtype=float)
                            # Keratin covariates enter as log2(RSEM+1) or as-is for STAR,
                            # which is already on a log scale. Rank-residual makes the
                            # transform irrelevant for the partial Spearman.
                            if exp != "star":
                                for j, name in enumerate(cov_cols):
                                    if name.startswith("krt"):
                                        z[:, j] = np.log2(np.clip(z[:, j], 0, None) + 1.0)
                            rho, p = partial_spearman(y, x, z)
                    lo, hi = fisher_ci(rho, n)
                    count_label = "tmb_le_10" if (filt == "count_le_500" and pred == "tmb") else filt
                    eligible = (
                        pred != "mut_silent"
                        and n >= MIN_N
                        and np.isfinite(rho)
                    )
                    rows.append(
                        {
                            "predictor": pred,
                            "expression": exp,
                            "filter": count_label,
                            "adjustment": adj,
                            "n": n,
                            "rho": rho,
                            "abs_rho": abs(rho) if np.isfinite(rho) else np.nan,
                            "p": p,
                            "ci95_lo": lo,
                            "ci95_hi": hi,
                            "eligible": bool(eligible),
                        }
                    )
    grid = pd.DataFrame(rows)
    grid["q_grid"] = np.nan
    elig = grid["eligible"].to_numpy()
    qvals = bh(grid.loc[elig, "p"].tolist())
    grid.loc[elig, "q_grid"] = qvals
    return grid.sort_values(["eligible", "abs_rho"], ascending=[False, False])


def anchor_summary(tab: pd.DataFrame) -> dict:
    """Yamamoto Fig. 4F reproduction.

    On the RSEM overlap the nonsynonymous count has median 69 and maximum
    1844, which is the published 2–69 vs 70–1,899 split (their maximum was
    1,899). Counting every MAF row, including silent and UTR, moves the
    median to ~102 and is not the figure.
    """
    sub = tab.dropna(subset=["mut_nonsyn", "rsem"])
    x = sub["mut_nonsyn"].to_numpy(dtype=float)
    y = sub["rsem"].to_numpy(dtype=float)
    rho, p = spearman(x, y)
    med = float(np.median(x))
    low = y[x <= med]
    high = y[x > med]
    mw = stats.mannwhitneyu(low, high, alternative="two-sided")
    return {
        "definition": "nonsynonymous MAF count vs PanCan RSEM CLDN4, all paired primary samples",
        "n": int(len(sub)),
        "mut_min": float(np.min(x)),
        "mut_median": med,
        "mut_max": float(np.max(x)),
        "n_low_le_median": int((x <= med).sum()),
        "n_high_gt_median": int((x > med).sum()),
        "median_rsem_low": float(np.median(low)),
        "median_rsem_high": float(np.median(high)),
        "mannwhitney_p": float(mw.pvalue),
        "spearman_rho": rho,
        "spearman_p": p,
        "direction": "high mutation count has lower CLDN4" if np.median(high) < np.median(low) else "high mutation count has higher CLDN4",
    }


def pick_winner(grid: pd.DataFrame) -> pd.Series:
    elig = grid[grid["eligible"]].copy()
    if elig.empty:
        raise SystemExit("no eligible correlation cell")
    elig["is_anchor"] = (
        elig["predictor"].eq("mut_all")
        & elig["expression"].eq("rsem")
        & elig["filter"].eq("all_paired")
        & elig["adjustment"].eq("none")
    )
    elig["unadjusted"] = elig["adjustment"].eq("none")
    expr_rank = {"rsem": 0, "star": 1, "protein": 2}
    pred_rank = {"mut_all": 0, "mut_nonsyn": 1, "tmb": 2}
    elig["expr_rank"] = elig["expression"].map(expr_rank)
    elig["pred_rank"] = elig["predictor"].map(pred_rank)
    elig = elig.sort_values(
        ["abs_rho", "is_anchor", "n", "unadjusted", "expr_rank", "pred_rank"],
        ascending=[False, False, False, False, True, True],
    )
    return elig.iloc[0]


AXIS = {
    "mut_nonsyn": "nonsynonymous mutation count",
    "mut_all": "MAF-row mutation count",
    "mut_silent": "silent mutation count",
    "TMB_NONSYNONYMOUS": "nonsynonymous TMB per Mb",
    "rsem": "CLDN4 RSEM",
    "star": "CLDN4 STAR log2(TPM+1)",
    "protein": "CLDN4 CPTAC protein",
}


def scatter(tab: pd.DataFrame, xcol: str, ycol: str, rho: float, p: float, n: int, title: str, path: Path, logx: bool) -> None:
    sub = tab.dropna(subset=[xcol, ycol])
    fig, ax = plt.subplots(figsize=(5.4, 4.4))
    x = sub[xcol].to_numpy(dtype=float)
    y = sub[ycol].to_numpy(dtype=float)
    purity = sub["purity"].to_numpy(dtype=float) if "purity" in sub else np.full(len(sub), np.nan)
    sc = ax.scatter(x, y, c=purity, cmap="viridis", s=18, alpha=0.85, linewidths=0)
    if logx and np.nanmin(x) >= 0:
        ax.set_xscale("log")
    ax.set_xlabel(AXIS.get(xcol, xcol))
    ax.set_ylabel(AXIS.get(ycol, ycol))
    ax.set_title(title)
    ax.text(
        0.02,
        0.98,
        f"Spearman ρ = {rho:+.3f}\nn = {n}\np = {p:.2e}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85, "edgecolor": "0.8"},
    )
    cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("ABSOLUTE purity")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def forest(grid: pd.DataFrame, path: Path) -> None:
    show = grid[grid["eligible"]].nlargest(12, "abs_rho").iloc[::-1]
    labels = [
        f"{r.predictor} × {r.expression}\n{r.filter}; {r.adjustment} (n={r.n})"
        for r in show.itertuples()
    ]
    fig, ax = plt.subplots(figsize=(8.2, 6.2))
    ypos = np.arange(len(show))
    ax.errorbar(
        show["rho"],
        ypos,
        xerr=[show["rho"] - show["ci95_lo"], show["ci95_hi"] - show["rho"]],
        fmt="o",
        color="#1f4e79",
        ecolor="#7f97ad",
        capsize=3,
    )
    ax.axvline(0, color="0.5", lw=0.8)
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Spearman ρ (Fisher 95% CI)")
    ax.set_title("Largest |ρ| cells in the pre-specified grid")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main() -> None:
    tab = build_table()
    anchor = anchor_summary(tab)
    grid = evaluate(tab)
    grid.to_csv(OUT_T / "rho_grid.tsv", sep="\t", index=False, float_format="%.6g")
    winner = pick_winner(grid)
    # Counts that describe the published 2–69 vs 70–1899 bins.
    desc = {}
    paired = tab.dropna(subset=["mut_all", "rsem"])
    for col in ["mut_all", "mut_nonsyn", "mut_silent"]:
        desc[col] = {
            "n_with_rsem": int(paired[col].notna().sum()) if col != "mut_silent" else int(tab.dropna(subset=[col, "rsem"]).shape[0]),
            "min": float(tab[col].min()),
            "median": float(tab[col].median()),
            "max": float(tab[col].max()),
            "median_on_rsem_overlap": float(paired[col].median()),
            "max_on_rsem_overlap": float(paired[col].max()),
        }
    payload = {
        "anchor": anchor,
        "count_summaries": desc,
        "winner": winner.to_dict(),
        "n_grid": int(len(grid)),
        "n_eligible": int(grid["eligible"].sum()),
        "min_n": MIN_N,
        "note": (
            "RPPA has Claudin-7 and no CLDN4 antibody, so protein is the CPTAC "
            "quantification. Spearman rho does not change if either margin is "
            "transformed monotonically, so log10(count) is not a separate test."
        ),
    }
    (OUT_T / "rho_summary.json").write_text(json.dumps(payload, indent=2, default=float))
    xcol = "TMB_NONSYNONYMOUS" if winner["predictor"] == "tmb" else winner["predictor"]
    scatter(
        tab,
        "mut_nonsyn",
        "rsem",
        anchor["spearman_rho"],
        anchor["spearman_p"],
        anchor["n"],
        "Yamamoto anchor: nonsynonymous count vs RSEM",
        OUT_F / "rho_anchor_scatter.png",
        logx=True,
    )
    winner_tab = apply_filter(tab, str(winner["predictor"]), str(winner["filter"]))
    scatter(
        winner_tab.dropna(subset=[xcol, winner["expression"]]),
        xcol,
        winner["expression"],
        float(winner["rho"]),
        float(winner["p"]),
        int(winner["n"]),
        f"Max |ρ|: {winner['predictor']} × {winner['expression']}, {winner['filter']}",
        OUT_F / "rho_winner_scatter.png",
        logx=(winner["predictor"] != "tmb"),
    )
    forest(grid, OUT_F / "rho_grid_forest.png")
    print(json.dumps({"anchor_rho": anchor["spearman_rho"], "anchor_p": anchor["spearman_p"], "anchor_n": anchor["n"],
                      "winner": winner.to_dict()}, indent=2, default=float))


if __name__ == "__main__":
    main()
