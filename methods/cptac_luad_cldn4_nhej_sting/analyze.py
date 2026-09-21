#!/usr/bin/env python3
"""CPTAC LUAD protein: CLDN4 vs DNA-PKcs / Ku / STING / TBK1 / IRF3 / HLA.

Public freeze v1.2 TMT only. Treatment-naive surgical LUAD (Gillette et al.
Cell 2020). No ICI labels. Honest pairwise-complete n (CLDN4 protein is missing
in a large fraction of tumors). No TACSTD2 gate.

Primary tests are two-sided Spearman. Sensitivity is rank-partial Spearman
after residualizing ranks on WES purity (DNA). WGS purity is a second residual,
not the pre-specified one. BH-FDR is within the 9-gene primary family, computed
separately for unadjusted and for WES-partial p-values.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# Current Ensembl (GRCh37 and GRCh38): PRKDC = ENSG00000253729.
# ENSG00000101868 is POLA1, not DNA-PKcs, and is not analyzed.
GENES = {
    "CLDN4": {
        "ensembl": "ENSG00000189143",
        "label": "CLDN4",
        "family": "predictor",
    },
    "PRKDC": {
        "ensembl": "ENSG00000253729",
        "label": "DNA-PKcs (PRKDC)",
        "family": "NHEJ",
    },
    "XRCC5": {
        "ensembl": "ENSG00000079246",
        "label": "Ku80 (XRCC5)",
        "family": "NHEJ",
    },
    "XRCC6": {
        "ensembl": "ENSG00000196419",
        "label": "Ku70 (XRCC6)",
        "family": "NHEJ",
    },
    "STING1": {
        "ensembl": "ENSG00000184584",
        "label": "STING (STING1)",
        "family": "STING",
    },
    "TBK1": {
        "ensembl": "ENSG00000183735",
        "label": "TBK1",
        "family": "STING",
    },
    "IRF3": {
        "ensembl": "ENSG00000126456",
        "label": "IRF3",
        "family": "STING",
    },
    "HLA-A": {
        "ensembl": "ENSG00000206503",
        "label": "HLA-A",
        "family": "HLA",
    },
    "HLA-B": {
        "ensembl": "ENSG00000234745",
        "label": "HLA-B",
        "family": "HLA",
    },
    "HLA-C": {
        "ensembl": "ENSG00000204525",
        "label": "HLA-C",
        "family": "HLA",
    },
}

# Pre-specified primary endpoints. FDR family = these 9 only.
PRIMARY = ["PRKDC", "XRCC5", "XRCC6", "STING1", "TBK1", "IRF3", "HLA-A", "HLA-B", "HLA-C"]

# Exploratory HLA rows. Own FDR. Not used to call the primary result.
EXPLORATORY = {
    "HLA-E": {"ensembl": "ENSG00000204592", "label": "HLA-E", "family": "HLA-exploratory"},
    "HLA-F": {"ensembl": "ENSG00000204642", "label": "HLA-F", "family": "HLA-exploratory"},
    "HLA-G": {"ensembl": "ENSG00000204632", "label": "HLA-G", "family": "HLA-exploratory"},
    "HLA-DRA": {"ensembl": "ENSG00000204287", "label": "HLA-DRA", "family": "HLA-exploratory"},
    "B2M": {"ensembl": "ENSG00000166710", "label": "B2M", "family": "HLA-exploratory"},
}

COMPOSITES = {
    "NHEJ3": ["PRKDC", "XRCC5", "XRCC6"],
    "STING3": ["STING1", "TBK1", "IRF3"],
    "HLA_ABC": ["HLA-A", "HLA-B", "HLA-C"],
}

PROT_NAME = "LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"
PHENO_NAME = "LUAD_phenotype.txt"


def find_row(index: pd.Index, ensembl: str) -> str | None:
    hits = [str(i) for i in index if str(i) == ensembl or str(i).startswith(ensembl + ".")]
    return hits[0] if hits else None


def spearman_pair(x: pd.Series, y: pd.Series) -> dict:
    d = pd.concat([x.rename("x"), y.rename("y")], axis=1).dropna()
    n = int(len(d))
    rec = {"n": n, "rho": np.nan, "p": np.nan, "ci_low": np.nan, "ci_high": np.nan, "tested": False}
    if n < 4 or d["x"].nunique() < 2 or d["y"].nunique() < 2:
        rec["note"] = "not tested"
        return rec
    rho, p = stats.spearmanr(d["x"], d["y"])
    lo, hi = fisher_ci(float(rho), n)
    rec.update({"rho": float(rho), "p": float(p), "ci_low": lo, "ci_high": hi, "tested": True, "note": "spearman"})
    return rec


def fisher_ci(rho: float, n: int, alpha: float = 0.05) -> tuple[float, float]:
    if n < 4 or not np.isfinite(rho) or abs(rho) >= 1:
        return (np.nan, np.nan)
    z = np.arctanh(np.clip(rho, -0.999999, 0.999999))
    se = 1.0 / np.sqrt(n - 3)
    zcrit = stats.norm.ppf(1 - alpha / 2)
    return (float(np.tanh(z - zcrit * se)), float(np.tanh(z + zcrit * se)))


def rank_partial(x: pd.Series, y: pd.Series, z: pd.Series) -> dict:
    """Pearson correlation of rank residuals. Rank-partial Spearman."""
    d = pd.concat([x.rename("x"), y.rename("y"), z.rename("z")], axis=1).dropna()
    n = int(len(d))
    rec = {"n": n, "rho": np.nan, "p": np.nan, "ci_low": np.nan, "ci_high": np.nan, "tested": False}
    if n < 6 or d["x"].nunique() < 2 or d["y"].nunique() < 2 or d["z"].nunique() < 2:
        rec["note"] = "not tested"
        return rec
    rx = d["x"].rank().to_numpy(float)
    ry = d["y"].rank().to_numpy(float)
    rz = d["z"].rank().to_numpy(float)

    def resid(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        design = np.column_stack([np.ones(len(b)), b])
        beta, *_ = np.linalg.lstsq(design, a, rcond=None)
        return a - design @ beta

    xr = resid(rx, rz)
    yr = resid(ry, rz)
    if np.std(xr) == 0 or np.std(yr) == 0:
        rec["note"] = "no residual variance"
        return rec
    rho, p = stats.pearsonr(xr, yr)
    lo, hi = fisher_ci(float(rho), n)
    rec.update(
        {
            "rho": float(rho),
            "p": float(p),
            "ci_low": lo,
            "ci_high": hi,
            "tested": True,
            "note": "rank-partial spearman",
        }
    )
    return rec


def bh(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    q = np.full(p.shape, np.nan)
    finite = np.where(np.isfinite(p))[0]
    m = len(finite)
    if m == 0:
        return q.tolist()
    order = finite[np.argsort(p[finite])]
    raw = np.empty(m)
    for rank, idx in enumerate(order, start=1):
        raw[rank - 1] = p[idx] * m / rank
    cmin = 1.0
    adj = np.empty(m)
    for k in range(m - 1, -1, -1):
        cmin = min(cmin, raw[k])
        adj[k] = min(cmin, 1.0)
    for k, idx in enumerate(order):
        q[idx] = adj[k]
    return q.tolist()


def call_direction(rho: float, q: float) -> str:
    if not np.isfinite(rho) or not np.isfinite(q):
        return "NOT_TESTED"
    if q < 0.05 and rho > 0:
        return "POSITIVE"
    if q < 0.05 and rho < 0:
        return "INVERSE"
    return "NULL"


def zscore(s: pd.Series) -> pd.Series:
    v = pd.to_numeric(s, errors="coerce")
    mu = v.mean(skipna=True)
    sd = v.std(skipna=True, ddof=1)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.nan, index=s.index)
    return (v - mu) / sd


def load_series(mat: pd.DataFrame, spec: dict) -> tuple[pd.Series, dict]:
    row = find_row(mat.index, spec["ensembl"])
    rec = {
        "ensembl": spec["ensembl"],
        "row_id": row or "",
        "row_present": row is not None,
        "n_samples": int(mat.shape[1]),
        "n_observed": 0,
        "n_na": int(mat.shape[1]),
    }
    if row is None:
        rec["note"] = "row absent"
        return pd.Series(np.nan, index=mat.columns), rec
    s = pd.to_numeric(mat.loc[row], errors="coerce")
    rec["n_observed"] = int(s.notna().sum())
    rec["n_na"] = int(s.isna().sum())
    rec["na_frac"] = float(rec["n_na"] / rec["n_samples"]) if rec["n_samples"] else 1.0
    rec["median"] = float(s.median()) if rec["n_observed"] else np.nan
    rec["sd"] = float(s.std(ddof=1)) if rec["n_observed"] > 1 else np.nan
    rec["note"] = "row present" if rec["n_observed"] else "row present but all NA"
    return s, rec


def add_fdr(rows: list[dict], pkey: str, qkey: str, callkey: str) -> None:
    qs = bh([r.get(pkey, np.nan) if r.get("tested") else np.nan for r in rows])
    for r, q in zip(rows, qs):
        r[qkey] = q
        r[callkey] = call_direction(r.get("rho", np.nan), q)


def plot_forest(primary: pd.DataFrame, path: Path) -> None:
    colors = {"NHEJ": "#1f4e79", "STING": "#b85c38", "HLA": "#2f6f4e"}
    fig, ax = plt.subplots(figsize=(8.2, 5.6))
    y = np.arange(len(primary))[::-1]
    for i, (_, r) in enumerate(primary.iterrows()):
        yi = y[i]
        c = colors.get(r["family"], "#333333")
        ax.plot([r["ci_low"], r["ci_high"]], [yi, yi], color=c, lw=1.6, solid_capstyle="round")
        ax.scatter([r["rho"]], [yi], s=42, color=c, zorder=3, label=None)
        if np.isfinite(r.get("partial_wes_rho", np.nan)):
            ax.scatter([r["partial_wes_rho"]], [yi], s=36, facecolors="none", edgecolors=c, linewidths=1.3, zorder=4)
    ax.axvline(0, color="#888888", lw=0.8)
    labels = [f"{r['label']}   n={int(r['n'])}" for _, r in primary.iterrows()]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Spearman ρ  (filled = unadjusted; open = partial | WES purity)")
    ax.set_xlim(-0.75, 0.75)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title("CPTAC LUAD protein  ·  CLDN4 vs DNA-PKcs / Ku / STING / HLA", fontsize=11)
    # legend proxies
    ax.scatter([], [], s=42, color="#1f4e79", label="NHEJ")
    ax.scatter([], [], s=42, color="#b85c38", label="STING axis")
    ax.scatter([], [], s=42, color="#2f6f4e", label="HLA-A/B/C")
    ax.scatter([], [], s=36, facecolors="none", edgecolors="#333333", linewidths=1.3, label="partial | WES")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_scatters(wide: pd.DataFrame, primary_stats: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(3, 3, figsize=(9.2, 8.4))
    for ax, gene in zip(axes.ravel(), PRIMARY):
        sub = wide[["CLDN4", gene]].dropna()
        ax.scatter(sub["CLDN4"], sub[gene], s=16, c="#243447", alpha=0.75, edgecolors="none")
        row = primary_stats[primary_stats.gene == gene].iloc[0]
        ax.set_title(f"{GENES[gene]['label']}\nρ={row['rho']:.2f}  p={row['p']:.3g}  n={int(row['n'])}", fontsize=8)
        ax.set_xlabel("CLDN4 protein", fontsize=7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=7)
    fig.suptitle("CPTAC LUAD TMT  ·  pairwise-complete CLDN4 protein", fontsize=11)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def one_endpoint(name: str, spec: dict, x: pd.Series, y: pd.Series, wes: pd.Series, wgs: pd.Series, family: str) -> dict:
    base = spearman_pair(x, y)
    wes_r = rank_partial(x, y, wes)
    wgs_r = rank_partial(x, y, wgs)
    pur = spearman_pair(y, wes)
    rec = {
        "gene": name,
        "label": spec["label"],
        "family": family,
        "ensembl": spec["ensembl"],
        "n": base["n"],
        "rho": base["rho"],
        "p": base["p"],
        "ci_low": base["ci_low"],
        "ci_high": base["ci_high"],
        "tested": base["tested"],
        "partial_wes_n": wes_r["n"],
        "partial_wes_rho": wes_r["rho"],
        "partial_wes_p": wes_r["p"],
        "partial_wes_ci_low": wes_r["ci_low"],
        "partial_wes_ci_high": wes_r["ci_high"],
        "partial_wes_tested": wes_r["tested"],
        "partial_wgs_n": wgs_r["n"],
        "partial_wgs_rho": wgs_r["rho"],
        "partial_wgs_p": wgs_r["p"],
        "partial_wgs_tested": wgs_r["tested"],
        "endpoint_vs_wes_n": pur["n"],
        "endpoint_vs_wes_rho": pur["rho"],
        "endpoint_vs_wes_p": pur["p"],
    }
    return rec


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/cptac_luad_cldn4_nhej_sting")
    ap.add_argument("--outdir", default="methods/cptac_luad_cldn4_nhej_sting")
    args = ap.parse_args()
    data = Path(args.data)
    res = Path(args.outdir) / "results"
    res.mkdir(parents=True, exist_ok=True)

    mat = pd.read_csv(data / PROT_NAME, sep="\t", index_col=0)
    mat.index = mat.index.astype(str)
    ph = pd.read_csv(data / PHENO_NAME, sep="\t", index_col=0)
    ph.index = ph.index.astype(str)
    wes = pd.to_numeric(ph["WES_purity"], errors="coerce")
    wgs = pd.to_numeric(ph["WGS_purity"], errors="coerce")
    wes = wes.reindex(mat.columns)
    wgs = wgs.reindex(mat.columns)

    presence = []
    series: dict[str, pd.Series] = {}
    for name, spec in {**GENES, **EXPLORATORY}.items():
        s, rec = load_series(mat, spec)
        rec["gene"] = name
        rec["label"] = spec["label"]
        rec["family"] = spec["family"]
        presence.append(rec)
        series[name] = s

    cldn4 = series["CLDN4"]
    primary_rows = []
    for name in PRIMARY:
        primary_rows.append(
            one_endpoint(name, GENES[name], cldn4, series[name], wes, wgs, GENES[name]["family"])
        )
    add_fdr(primary_rows, "p", "q", "call_q")
    # WES-partial FDR on the same 9, only where tested
    wes_ps = []
    for r in primary_rows:
        wes_ps.append(r["partial_wes_p"] if r["partial_wes_tested"] else np.nan)
    for r, q in zip(primary_rows, bh(wes_ps)):
        r["partial_wes_q"] = q
        r["partial_wes_call"] = call_direction(r["partial_wes_rho"], q)

    explor_rows = []
    for name, spec in EXPLORATORY.items():
        explor_rows.append(one_endpoint(name, spec, cldn4, series[name], wes, wgs, spec["family"]))
    add_fdr(explor_rows, "p", "q", "call_q")

    # Composites: z-score each member on all observed tumors, mean where all members observed.
    comp_rows = []
    comp_series = {}
    for cname, members in COMPOSITES.items():
        z = pd.concat({m: zscore(series[m]) for m in members}, axis=1)
        comp = z.mean(axis=1, skipna=False)
        comp.name = cname
        comp_series[cname] = comp
        spec = {"ensembl": "+".join(GENES[m]["ensembl"] for m in members), "label": cname}
        rec = one_endpoint(cname, spec, cldn4, comp, wes, wgs, "composite")
        rec["members"] = ",".join(members)
        rec["n_member_complete"] = int(comp.notna().sum())
        comp_rows.append(rec)
    add_fdr(comp_rows, "p", "q", "call_q")
    for r, q in zip(comp_rows, bh([r["partial_wes_p"] if r["partial_wes_tested"] else np.nan for r in comp_rows])):
        r["partial_wes_q"] = q
        r["partial_wes_call"] = call_direction(r["partial_wes_rho"], q)

    # Missingness check: CLDN4 protein is NA in 31/110. Complete-case Spearman
    # does not use those tumors. MWU asks whether the endpoint differs in the
    # NA set. This is not a primary endpoint.
    miss_rows = []
    observed = cldn4.notna()
    for name in [*PRIMARY, *EXPLORATORY]:
        y = series[name]
        a = y[observed].dropna()
        b = y[~observed].dropna()
        rec = {
            "gene": name,
            "n_cldn4_observed": int(len(a)),
            "n_cldn4_na": int(len(b)),
            "median_cldn4_observed": float(a.median()) if len(a) else np.nan,
            "median_cldn4_na": float(b.median()) if len(b) else np.nan,
            "delta_median_na_minus_obs": np.nan,
            "p": np.nan,
            "tested": False,
        }
        if len(a) and len(b):
            rec["delta_median_na_minus_obs"] = float(b.median() - a.median())
        if len(a) >= 4 and len(b) >= 4 and a.nunique() > 1 and b.nunique() > 1:
            u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            rec["U"] = float(u)
            rec["p"] = float(p)
            rec["tested"] = True
        miss_rows.append(rec)
    for rec, qv in zip(miss_rows, bh([r["p"] if r["tested"] else np.nan for r in miss_rows])):
        rec["q"] = qv

    # QC: do the NHEJ subunits track each other?
    qc = []
    for a, b in [("PRKDC", "XRCC5"), ("PRKDC", "XRCC6"), ("XRCC5", "XRCC6"), ("HLA-A", "HLA-B"), ("HLA-A", "HLA-C"), ("HLA-B", "HLA-C"), ("STING1", "TBK1"), ("STING1", "IRF3"), ("TBK1", "IRF3")]:
        sp = spearman_pair(series[a], series[b])
        qc.append({"a": a, "b": b, **sp})

    cldn4_purity = {
        "CLDN4_vs_WES": spearman_pair(cldn4, wes),
        "CLDN4_vs_WGS": spearman_pair(cldn4, wgs),
    }

    presence_df = pd.DataFrame(presence)
    primary_df = pd.DataFrame(primary_rows)
    explor_df = pd.DataFrame(explor_rows)
    comp_df = pd.DataFrame(comp_rows)
    qc_df = pd.DataFrame(qc)
    miss_df = pd.DataFrame(miss_rows)

    presence_df.to_csv(res / "presence.tsv", sep="\t", index=False)
    primary_df.to_csv(res / "spearman_primary.tsv", sep="\t", index=False)
    explor_df.to_csv(res / "spearman_exploratory_hla.tsv", sep="\t", index=False)
    comp_df.to_csv(res / "spearman_composites.tsv", sep="\t", index=False)
    qc_df.to_csv(res / "qc_within_family.tsv", sep="\t", index=False)
    miss_df.to_csv(res / "cldn4_missingness_mwu.tsv", sep="\t", index=False)

    wide = pd.DataFrame({k: series[k] for k in ["CLDN4", *PRIMARY, *EXPLORATORY]})
    for cname, s in comp_series.items():
        wide[cname] = s
    wide["WES_purity"] = wes
    wide["WGS_purity"] = wgs
    wide.index.name = "case_id"
    wide.to_csv(res / "sample_scores.tsv", sep="\t")

    plot_forest(primary_df, res / "fig_forest_cldn4_protein.png")
    plot_scatters(wide, primary_df, res / "fig_scatter_primary.png")

    n_tumors = int(mat.shape[1])
    n_cldn4 = int(cldn4.notna().sum())
    summary = {
        "cohort": "CPTAC-LUAD",
        "source": "data_freeze_v1.2_reorganized LUAD tumor TMT gene abundance",
        "citation": "Gillette et al. Cell 2020",
        "treatment": "treatment-naive surgical resection; no ICI labels",
        "predictor": "CLDN4 protein ENSG00000189143",
        "dna_pkcs_id": "ENSG00000253729 (PRKDC). ENSG00000101868 is POLA1 and was not used.",
        "n_tumors": n_tumors,
        "n_cldn4_protein": n_cldn4,
        "n_cldn4_na": n_tumors - n_cldn4,
        "n_wes_purity": int(wes.notna().sum()),
        "n_wgs_purity": int(wgs.notna().sum()),
        "primary_fdr_family": PRIMARY,
        "cldn4_vs_purity": cldn4_purity,
        "primary": primary_rows,
        "composites": comp_rows,
        "exploratory_hla": explor_rows,
        "b2m": next(r for r in presence if r["gene"] == "B2M"),
    }
    (res / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    show = primary_df[
        ["gene", "n", "rho", "p", "q", "call_q", "partial_wes_n", "partial_wes_rho", "partial_wes_p", "partial_wes_q", "partial_wes_call"]
    ]
    print(show.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print("--- composites ---")
    print(comp_df[["gene", "n", "rho", "p", "q", "call_q", "partial_wes_n", "partial_wes_rho", "partial_wes_p", "partial_wes_q"]].to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print("--- exploratory ---")
    print(explor_df[["gene", "n", "rho", "p", "q", "call_q"]].to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print("--- cldn4 purity ---")
    print(json.dumps(cldn4_purity, indent=2))
    print("--- missingness MWU ---")
    print(miss_df[["gene", "n_cldn4_observed", "n_cldn4_na", "delta_median_na_minus_obs", "p", "q"]].to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print("--- presence cldn4 / b2m ---")
    print(presence_df[presence_df.gene.isin(["CLDN4", "B2M", "HLA-G"])][["gene", "row_present", "n_observed", "n_na"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
